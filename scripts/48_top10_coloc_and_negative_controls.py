from __future__ import annotations

import importlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from remote_tabix import query_region


benchmark_module = importlib.import_module("46_external_eqtl_variant_benchmark")
DATASETS = benchmark_module.DATASETS

EQL_COLUMNS = [
    "molecular_trait_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant",
    "ma_samples",
    "maf",
    "pvalue",
    "beta",
    "se",
    "type",
    "ac",
    "an",
    "r2",
    "molecular_trait_object_id",
    "gene_id",
    "median_tpm",
    "rsid",
]
GWAS_COLUMNS = [
    "chromosome",
    "position",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "effect_allele_frequency",
    "pvalue",
    "variant_id",
    "rsid",
    "FreqSE",
    "MinFreq",
    "MaxFreq",
    "Direction",
    "HetISq",
    "HetChiSq",
    "HetDf",
    "HetPVal",
    "hm_coordinate_conversion",
    "hm_code",
]
GWAS_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST90292001-GCST90293000/GCST90292538/harmonised/GCST90292538.h.tsv.gz"
)
REGION_HALF_WIDTH = 500_000


def urls(study_id: str, dataset_id: str) -> tuple[str, str]:
    base = (
        "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
        f"{study_id}/{dataset_id}/{dataset_id}"
    )
    return base + ".all.tsv.gz", base + ".all.tsv.gz.tbi"


def logsum(values: np.ndarray) -> float:
    if len(values) == 0:
        return -np.inf
    maximum = float(np.max(values))
    return maximum + float(np.log(np.exp(values - maximum).sum()))


def coloc_abf(
    frame: pd.DataFrame,
    p1: float,
    p2: float,
    p12: float,
    sd_prior: float,
) -> tuple[dict, pd.DataFrame]:
    work = frame.copy()
    for prefix in ["gwas", "eql"]:
        variance = work[f"{prefix}_se"].to_numpy(dtype=float) ** 2
        z = work[f"{prefix}_beta"].to_numpy(dtype=float) / np.sqrt(variance)
        r = sd_prior**2 / (sd_prior**2 + variance)
        work[f"{prefix}_labf"] = 0.5 * (
            np.log(1 - r) + r * z**2
        )
    evidence = {
        "H0": 1.0,
        "H1": p1 * np.exp(logsum(work["eql_labf"].to_numpy())),
        "H2": p2 * np.exp(logsum(work["gwas_labf"].to_numpy())),
        "H3": p1
        * p2
        * np.exp(logsum((work["eql_labf"] + work["gwas_labf"]).to_numpy())),
        "H4": p12
        * np.exp(logsum((work["eql_labf"] + work["gwas_labf"]).to_numpy())),
    }
    total = sum(evidence.values())
    posterior = {key: value / total for key, value in evidence.items()}
    return posterior, work


def select_hits() -> pd.DataFrame:
    table_dir = ROOT / "results" / "tables"
    top = pd.read_csv(table_dir / "external_eqtl_top10_triangulation.csv")
    matched = top[
        top["query_group"].eq("top10_matched")
        & top["any_nominal"].astype(bool)
    ].copy()
    lookups = pd.read_csv(table_dir / "external_eqtl_all_variant_lookups.csv")
    lookups = lookups[lookups["query_group"].eq("top10_matched")].copy()
    hits = (
        lookups.merge(
            matched[
                [
                    "variant_id",
                    "canonical_variant",
                    "target",
                    "gene_id",
                ]
            ],
            on=[
                "variant_id",
                "canonical_variant",
                "target",
                "gene_id",
            ],
            how="inner",
        )
        .sort_values("pvalue")
        .drop_duplicates(["target", "gene_id"])
        .copy()
    )
    dataset_map = {
        (study, condition): (study_id, dataset_id)
        for study_id, dataset_id, study, condition, _ in DATASETS
    }
    hits["study_id"] = [
        dataset_map[(study, condition)][0]
        for study, condition in zip(hits["study"], hits["condition"])
    ]
    hits["dataset_id"] = [
        dataset_map[(study, condition)][1]
        for study, condition in zip(hits["study"], hits["condition"])
    ]
    candidates = pd.read_csv(
        table_dir / "ranked_variant_gene_candidates_final.csv"
    )
    coordinate_lookup = candidates.drop_duplicates("variant_id").set_index(
        "variant_id"
    )[["chromosome", "position"]]
    hits["chromosome"] = hits["variant_id"].map(coordinate_lookup["chromosome"])
    hits["position"] = hits["variant_id"].map(coordinate_lookup["position"])
    return hits


def analyze_hit(hit: pd.Series, cache_dir: Path) -> tuple[dict, pd.DataFrame]:
    chromosome = str(int(hit["chromosome"]))
    position = int(hit["position"])
    gene_id = str(hit["gene_id"])
    start = max(1, position - REGION_HALF_WIDTH)
    end = position + REGION_HALF_WIDTH
    data_url, index_url = urls(hit["study_id"], hit["dataset_id"])
    eql_lines = query_region(
        data_url,
        index_url,
        chromosome,
        start,
        end,
        cache_dir,
    )
    gwas_lines = query_region(
        GWAS_BASE,
        GWAS_BASE + ".tbi",
        chromosome,
        start,
        end,
        cache_dir,
    )
    eqtl = pd.DataFrame(
        [line.split("\t") for line in eql_lines],
        columns=EQL_COLUMNS,
    )
    eqtl = eqtl[eqtl["molecular_trait_id"].eq(gene_id)].copy()
    for column in ["position", "maf", "pvalue", "beta", "se", "median_tpm"]:
        eqtl[column] = pd.to_numeric(eqtl[column], errors="coerce")
    gwas = pd.DataFrame(
        [line.split("\t") for line in gwas_lines],
        columns=GWAS_COLUMNS,
    )
    for column in [
        "position",
        "beta",
        "se",
        "effect_allele_frequency",
        "pvalue",
    ]:
        gwas[column] = pd.to_numeric(gwas[column], errors="coerce")
    eqtl = eqtl.dropna(subset=["rsid", "beta", "se", "maf", "position"])
    gwas = gwas.dropna(
        subset=["rsid", "beta", "se", "effect_allele_frequency", "position"]
    )
    eqtl = eqtl[
        eqtl["maf"].between(0.01, 0.99)
        & eqtl["se"].gt(0)
        & eqtl["beta"].notna()
    ]
    gwas = gwas[
        gwas["effect_allele_frequency"].between(0.01, 0.99)
        & gwas["se"].gt(0)
    ]
    merged = eqtl.merge(
        gwas,
        on="rsid",
        suffixes=("_eql", "_gwas"),
    )
    merged["eql_beta_aligned"] = np.where(
        merged["alt"].str.upper().eq(merged["effect_allele"].str.upper()),
        merged["beta_eql"],
        np.where(
            merged["alt"].str.upper().eq(merged["other_allele"].str.upper()),
            -merged["beta_eql"],
            np.nan,
        ),
    )
    merged = merged.dropna(subset=["eql_beta_aligned"])
    frame = pd.DataFrame(
        {
            "rsid": merged["rsid"],
            "position": merged["position_eql"],
            "gwas_beta": merged["beta_gwas"],
            "gwas_se": merged["se_gwas"],
            "eql_beta": merged["eql_beta_aligned"],
            "eql_se": merged["se_eql"],
            "gwas_p": merged["pvalue_gwas"],
            "eql_p": merged["pvalue_eql"],
        }
    ).drop_duplicates("rsid")
    summary = {
        "target": hit["target"],
        "gene_id": gene_id,
        "variant": hit["canonical_variant"],
        "rsid": hit["rsid"],
        "eqtl_study": hit["study"],
        "eqtl_condition": hit["condition"],
        "region": f"chr{chromosome}:{start}-{end}",
        "variants_harmonized": len(frame),
        "eqtl_pvalue": float(hit["pvalue"]),
        "coloc_status": "evaluable" if len(frame) >= 50 else "too_few_variants",
    }
    if len(frame) < 50:
        summary.update(
            {
                "H0": np.nan,
                "H1": np.nan,
                "H2": np.nan,
                "H3": np.nan,
                "H4": np.nan,
                "coloc_supported": False,
            }
        )
        return summary, frame.head(0)
    posterior, work = coloc_abf(
        frame,
        p1=1e-4,
        p2=1e-4,
        p12=1e-5,
        sd_prior=0.15,
    )
    summary.update(posterior)
    summary["coloc_supported"] = bool(posterior["H4"] >= 0.80)
    return summary, work


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    cache_dir = ROOT / "work" / "tabix_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    hits = select_hits()
    rows = []
    work_frames = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(analyze_hit, row, cache_dir): row["target"]
            for _, row in hits.iterrows()
        }
        for future in as_completed(futures):
            target = futures[future]
            try:
                summary, work = future.result()
                rows.append(summary)
                if len(work):
                    work = work.copy()
                    work["target"] = target
                    work_frames.append(work)
                print(
                    f"Completed coloc: {target}",
                    flush=True,
                )
            except Exception as exc:
                rows.append(
                    {
                        "target": target,
                        "coloc_status": "error",
                        "error": str(exc),
                    }
                )
    results = pd.DataFrame(rows)
    results.to_csv(
        table_dir / "top10_eqtl_coloc_results.csv",
        index=False,
    )
    if work_frames:
        pd.concat(work_frames, ignore_index=True).to_csv(
            table_dir / "top10_eqtl_coloc_variants.csv",
            index=False,
        )
    triangulation = pd.read_csv(
        table_dir / "external_eqtl_top10_triangulation.csv"
    )
    matched_hits = int(
        triangulation.loc[
            triangulation["query_group"].eq("top10_matched"),
            "any_nominal",
        ].sum()
    )
    control_hits = int(
        triangulation.loc[
            triangulation[
                "query_group"
            ].eq("mismatched_target_control"),
            "any_nominal",
        ].sum()
    )
    fisher_table = [
        [matched_hits, 10 - matched_hits],
        [control_hits, 10 - control_hits],
    ]
    odds_ratio, fisher_p = fisher_exact(fisher_table)
    summary = {
        "top10_pairs": 10,
        "matched_eqtl_hits": matched_hits,
        "mismatched_target_eqtl_hits": control_hits,
        "fisher_odds_ratio": float(odds_ratio),
        "fisher_p": float(fisher_p),
        "coloc_evaluable": int(
            results["coloc_status"].eq("evaluable").sum()
        ),
        "coloc_supported_h4_ge_0_8": int(
            results["coloc_supported"].fillna(False).astype(bool).sum()
        ),
        "note": (
            "Formal coloc.abf used the best exact-variant eQTL dataset for each "
            "matched top-10 pair. Mismatched target assignment provided the "
            "negative control."
        ),
    }
    (table_dir / "top10_eqtl_coloc_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(results.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
