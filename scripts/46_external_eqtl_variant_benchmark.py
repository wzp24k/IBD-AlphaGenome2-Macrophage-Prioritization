from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy
from remote_tabix import query_point


DATASETS = [
    ("QTS000001", "QTD000001", "Alasoo_2018", "macrophage_naive", 84),
    ("QTS000001", "QTD000006", "Alasoo_2018", "macrophage_IFNg", 84),
    ("QTS000002", "QTD000021", "BLUEPRINT", "monocyte_naive", 191),
    ("QTS000021", "QTD000384", "Nedelec_2016", "macrophage_naive", 163),
    ("QTS000024", "QTD000409", "Quach_2016", "monocyte_naive", 200),
    ("QTS000024", "QTD000414", "Quach_2016", "monocyte_LPS", 184),
    ("QTS000026", "QTD000499", "Schmiedel_2018", "monocyte_CD16_naive", 90),
    ("QTS000026", "QTD000504", "Schmiedel_2018", "monocyte_naive", 91),
]
CORE_DATASETS = [DATASETS[0], DATASETS[1], DATASETS[3], DATASETS[5]]
COLUMNS = [
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


def urls(study_id: str, dataset_id: str) -> tuple[str, str]:
    base = (
        "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
        f"{study_id}/{dataset_id}/{dataset_id}"
    )
    return base + ".all.tsv.gz", base + ".all.tsv.gz.tbi"


def prepare_index(dataset: tuple) -> None:
    study_id, dataset_id, *_ = dataset
    _, index_url = urls(study_id, dataset_id)
    cache_dir = ROOT / "work" / "tabix_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    index_path = cache_dir / f"ftp.ebi.ac.uk_{dataset_id}.all.tsv.gz.tbi"
    if index_path.exists() and index_path.stat().st_size > 0:
        return
    response = requests.get(index_url, timeout=180)
    response.raise_for_status()
    index_path.write_bytes(response.content)


def candidate_variants() -> pd.DataFrame:
    frame = pd.read_csv(
        ROOT / "results" / "tables" / "ranked_variant_gene_candidates_final.csv"
    )
    frame = frame.dropna(
        subset=[
            "variant_id",
            "canonical_variant",
            "l2g_target_id",
            "l2g_target_symbol",
            "chromosome",
            "position",
            "referenceAllele",
            "alternateAllele",
        ]
    ).copy()
    frame["rsid"] = (
        frame["rsIds"]
        .fillna("")
        .astype(str)
        .str.split(r"[,;]")
        .str[0]
        .str.strip()
    )
    frame["query_group"] = "primary"
    frame["query_key"] = frame["variant_id"] + "|" + frame["l2g_target_id"]
    return frame


def top10_negative_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    top = (
        frame.sort_values(["priority_rank", "variant_id"])
        .drop_duplicates("variant_id")
        .head(10)
        .reset_index(drop=True)
    )
    negatives = top.copy()
    negatives["correct_target_id"] = top["l2g_target_id"]
    negatives["correct_target_symbol"] = top["l2g_target_symbol"]
    negatives["l2g_target_id"] = top["l2g_target_id"].shift(-1)
    negatives["l2g_target_symbol"] = top["l2g_target_symbol"].shift(-1)
    negatives.loc[negatives.index[-1], "l2g_target_id"] = top.iloc[0][
        "l2g_target_id"
    ]
    negatives.loc[negatives.index[-1], "l2g_target_symbol"] = top.iloc[0][
        "l2g_target_symbol"
    ]
    negatives["query_group"] = "mismatched_target_control"
    negatives["query_key"] = (
        negatives["variant_id"]
        + "|"
        + negatives["l2g_target_id"]
        + "|negative"
    )
    return negatives


def top10_matched_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    top = (
        frame.sort_values(["priority_rank", "variant_id"])
        .drop_duplicates("variant_id")
        .head(10)
        .copy()
    )
    top["query_group"] = "top10_matched"
    top["query_key"] = (
        top["variant_id"] + "|" + top["l2g_target_id"] + "|top10"
    )
    return top


def query_one(
    dataset: tuple,
    variant: pd.Series,
    cache_dir: Path,
) -> dict | None:
    study_id, dataset_id, study_label, condition, sample_size = dataset
    chromosome = str(int(variant["chromosome"]))
    position = int(variant["position"])
    gene_id = str(variant["l2g_target_id"])
    gene = str(variant["l2g_target_symbol"])
    ref = str(variant["referenceAllele"])
    alt = str(variant["alternateAllele"])
    rsid = str(variant["rsid"])
    data_url, index_url = urls(study_id, dataset_id)
    lines = query_point(
        data_url,
        index_url,
        chromosome,
        position,
        cache_dir,
    )
    expected_variants = {
        f"{chromosome}_{position}_{ref}_{alt}",
        f"chr{chromosome}_{position}_{ref}_{alt}",
    }
    for line in lines:
        fields = line.split("\t")
        if len(fields) != len(COLUMNS):
            continue
        record = dict(zip(COLUMNS, fields))
        if record["molecular_trait_id"] != gene_id:
            continue
        identifier_match = (
            (rsid and record["rsid"] == rsid)
            or record["variant"] in expected_variants
        )
        if not identifier_match:
            continue
        beta = float(record["beta"])
        if record["alt"].upper() != alt.upper():
            if record["alt"].upper() == ref.upper():
                beta = -beta
            else:
                continue
        return {
            "query_key": variant["query_key"],
            "query_group": variant["query_group"],
            "variant_id": variant["variant_id"],
            "canonical_variant": variant["canonical_variant"],
            "rsid": rsid,
            "target": gene,
            "gene_id": gene_id,
            "study": study_label,
            "condition": condition,
            "sample_size": sample_size,
            "pvalue": float(record["pvalue"]),
            "beta": beta,
            "se": float(record["se"]),
            "maf": float(record["maf"]) if record["maf"] != "NA" else None,
            "median_tpm": float(record["median_tpm"])
            if record["median_tpm"] != "NA"
            else None,
        }
    return None


def summarize_eqtl(
    results: pd.DataFrame,
    primary: pd.DataFrame,
    top10_matched: pd.DataFrame,
    negatives: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    if results.empty:
        results = pd.DataFrame(
            columns=[
                "query_key",
                "query_group",
                "variant_id",
                "canonical_variant",
                "target",
                "gene_id",
                "study",
                "condition",
                "pvalue",
            ]
        )
    best = (
        results.sort_values("pvalue")
        .drop_duplicates(["query_key", "study", "condition"])
        .copy()
    )
    best["nominal"] = best["pvalue"] < 0.05
    observed = (
        best.groupby(
            [
                "query_key",
                "query_group",
                "variant_id",
                "canonical_variant",
                "target",
                "gene_id",
            ],
            as_index=False,
        )
        .agg(
            datasets_queried=("study", "nunique"),
            nominal_datasets=("nominal", "sum"),
            minimum_pvalue=("pvalue", "min"),
            minimum_qvalue=("pvalue", "min"),
        )
    )
    template_parts = []
    for frame in [primary, top10_matched, negatives]:
        template_parts.append(
            frame[
                [
                    "query_key",
                    "query_group",
                    "variant_id",
                    "canonical_variant",
                    "l2g_target_symbol",
                    "l2g_target_id",
                ]
            ].rename(
                columns={
                    "l2g_target_symbol": "target",
                    "l2g_target_id": "gene_id",
                }
            )
        )
    query_template = pd.concat(
        template_parts,
        ignore_index=True,
    ).drop_duplicates("query_key")
    query_template["expected_datasets"] = np.where(
        query_template["query_group"].eq("primary"),
        len(CORE_DATASETS),
        len(DATASETS),
    )
    per_query = query_template.merge(
        observed,
        on=[
            "query_key",
            "query_group",
            "variant_id",
            "canonical_variant",
            "target",
            "gene_id",
        ],
        how="left",
    )
    per_query["datasets_queried"] = per_query["expected_datasets"]
    per_query["nominal_datasets"] = (
        per_query["nominal_datasets"].fillna(0).astype(int)
    )
    per_query["replicated_nominal"] = per_query["nominal_datasets"] >= 2
    per_query["any_nominal"] = per_query["nominal_datasets"] >= 1

    top_positive = per_query[
        per_query["query_group"].eq("top10_matched")
    ].copy()
    top_negative = per_query[
        per_query["query_group"].eq("mismatched_target_control")
    ].copy()
    top_positive["pair_type"] = "matched_target"
    top_negative["pair_type"] = "mismatched_target_control"
    top_triangulation = pd.concat(
        [top_positive, top_negative],
        ignore_index=True,
    )
    primary_queries = per_query[per_query["query_group"].eq("primary")]
    summary = {
        "variants_queried": int(primary["query_key"].nunique()),
        "core_datasets_queried": len(CORE_DATASETS),
        "top10_datasets_queried": len(DATASETS),
        "query_errors": int(
            results["query_error"].notna().sum()
            if "query_error" in results.columns
            else 0
        ),
        "variants_with_any_nominal_eqtl": int(
            primary_queries["any_nominal"].sum()
        ),
        "variants_with_replicated_nominal_eqtl": int(
            primary_queries["replicated_nominal"].sum()
        ),
        "top10_matched_nominal_pairs": int(
            top_positive["any_nominal"].sum()
        ),
        "top10_mismatched_nominal_pairs": int(
            top_negative["any_nominal"].sum()
        ),
        "note": (
            "Exact-variant lookups used the nominated target from the primary "
            "analysis. Mismatched-target top-10 pairs provide a negative control."
        ),
    }
    return per_query, top_triangulation, summary


def main() -> None:
    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    cache_dir = ROOT / "work" / "tabix_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    primary = candidate_variants()
    top_matched = top10_matched_pairs(primary)
    negatives = top10_negative_pairs(primary)

    lookup_path = table_dir / "external_eqtl_all_variant_lookups.csv"
    if "--reuse" in sys.argv and lookup_path.exists():
        results = pd.read_csv(lookup_path)
        print(f"Reused {len(results)} existing eQTL lookup records")
    else:
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = [
                executor.submit(prepare_index, dataset) for dataset in DATASETS
            ]
            for future in as_completed(futures):
                future.result()

        jobs = []
        for _, row in primary.iterrows():
            for dataset in CORE_DATASETS:
                jobs.append((dataset, row))
        for frame in [top_matched, negatives]:
            for _, row in frame.iterrows():
                for dataset in DATASETS:
                    jobs.append((dataset, row))

        records = []
        failures = []
        job_count = len(jobs)
        with ThreadPoolExecutor(max_workers=32) as executor:
            futures = {
                executor.submit(executor_query, dataset, row, cache_dir): (
                    dataset[2],
                    dataset[3],
                    row["query_key"],
                )
                for dataset, row in jobs
            }
            for completed, future in enumerate(as_completed(futures), start=1):
                study, condition, query_key = futures[future]
                try:
                    result = future.result()
                    if result is not None:
                        records.append(result)
                except Exception as exc:
                    failures.append(
                        {
                            "query_key": query_key,
                            "study": study,
                            "condition": condition,
                            "query_error": str(exc),
                        }
                    )
                if completed % 250 == 0:
                    print(
                        f"Completed {completed}/{job_count} eQTL lookups",
                        flush=True,
                    )

        results = pd.DataFrame(records + failures)
        results.to_csv(lookup_path, index=False)

    per_variant, top_triangulation, summary = summarize_eqtl(
        results,
        primary,
        top_matched,
        negatives,
    )
    per_variant.to_csv(
        table_dir / "external_eqtl_variant_benchmark.csv",
        index=False,
    )
    top_triangulation.to_csv(
        table_dir / "external_eqtl_top10_triangulation.csv",
        index=False,
    )
    (table_dir / "external_eqtl_benchmark_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(
        top_triangulation.sort_values(
            ["query_group", "minimum_pvalue"],
            ascending=[True, True],
        ).to_string(index=False)
    )


def executor_query(dataset: tuple, row: pd.Series, cache_dir: Path):
    return query_one(dataset, row, cache_dir)


if __name__ == "__main__":
    main()
