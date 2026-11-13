from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from remote_tabix import query_region


CHROMOSOME = "5"
START = 149_500_000
END = 151_000_000
CAMK2A_GENE = "ENSG00000070808"
CAMK2A_RSID = "rs1549913"
EQL_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "QTS000001/QTD000006/QTD000006"
)
GWAS_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST90292001-GCST90293000/GCST90292538/harmonised/GCST90292538.h.tsv.gz"
)
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


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    cache_dir = ROOT / "work" / "tabix_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    eql_lines = query_region(
        EQL_BASE + ".all.tsv.gz",
        EQL_BASE + ".all.tsv.gz.tbi",
        CHROMOSOME,
        START,
        END,
        cache_dir,
    )
    gwas_lines = query_region(
        GWAS_BASE,
        GWAS_BASE + ".tbi",
        CHROMOSOME,
        START,
        END,
        cache_dir,
    )
    eqtl = pd.DataFrame(
        [line.split("\t") for line in eql_lines],
        columns=EQL_COLUMNS,
    )
    eqtl = eqtl[eqtl["molecular_trait_id"].eq(CAMK2A_GENE)].copy()
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
    if len(frame) < 10:
        raise RuntimeError("Too few harmonized CAMK2A variants for colocalization.")
    sensitivity_rows = []
    main_posterior = None
    main_work = None
    for sd_prior in [0.05, 0.10, 0.15, 0.20]:
        for p12 in [1e-5, 5e-5]:
            posterior, work = coloc_abf(
                frame,
                p1=1e-4,
                p2=1e-4,
                p12=p12,
                sd_prior=sd_prior,
            )
            sensitivity_rows.append(
                {
                    "sd_prior": sd_prior,
                    "p12": p12,
                    "variants": len(frame),
                    **posterior,
                }
            )
            if main_posterior is None:
                main_posterior = posterior
                main_work = work
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(
        table_dir / "camk2a_coloc_sensitivity.csv",
        index=False,
    )
    main_work.to_csv(
        table_dir / "camk2a_coloc_variants.csv",
        index=False,
    )
    exact = main_work[main_work["rsid"].eq(CAMK2A_RSID)]
    summary = {
        "method": "coloc.abf single-causal-variant approximation",
        "eqtl_dataset": "Alasoo_2018 macrophage IFN-gamma",
        "gwas_dataset": "GCST90292538",
        "region": f"chr{CHROMOSOME}:{START}-{END}",
        "eqtl_gene": "CAMK2A",
        "variants_harmonized": len(frame),
        "primary_sd_prior": 0.15,
        "primary_p12": 1e-5,
        "primary_posterior": main_posterior,
        "exact_variant": exact.iloc[0].to_dict() if len(exact) else None,
    }
    (table_dir / "camk2a_coloc_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(sensitivity.to_string(index=False))


if __name__ == "__main__":
    main()
