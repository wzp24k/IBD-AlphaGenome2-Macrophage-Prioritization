from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

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
VARIANTS = [
    ("ERGIC1", "ENSG00000113719", "rs564109", "5", 172897894, "A", "T"),
    ("GCKR", "ENSG00000084734", "rs780094", "2", 27518370, "T", "C"),
    ("ABCB5", "ENSG00000004846", "rs76380568", "7", 20542073, "T", "C"),
    ("IRF2BP2", "ENSG00000168264", "rs10910476", "1", 234599210, "C", "T"),
    ("CAMK2A", "ENSG00000070808", "rs1549913", "5", 150230409, "G", "A"),
]
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


def query_one(
    dataset: tuple,
    variant: tuple,
    cache_dir: Path,
) -> dict | None:
    study_id, dataset_id, study_label, condition, sample_size = dataset
    gene, gene_id, rsid, chromosome, position, ref, alt = variant
    data_url, index_url = urls(study_id, dataset_id)
    lines = query_point(
        data_url,
        index_url,
        chromosome,
        position,
        cache_dir,
    )
    for line in lines:
        fields = line.split("\t")
        if len(fields) != len(COLUMNS):
            continue
        record = dict(zip(COLUMNS, fields))
        if record["molecular_trait_id"] != gene_id:
            continue
        if record["rsid"] != rsid and record["variant"] != f"chr{chromosome}_{position}_{ref}_{alt}":
            continue
        beta = float(record["beta"])
        if record["alt"] != alt:
            beta = -beta
        return {
            "target": gene,
            "gene_id": gene_id,
            "rsid": rsid,
            "chromosome": chromosome,
            "position": position,
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
            "effect_allele_matched": record["alt"] == alt,
        }
    return None


def main() -> None:
    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    cache_dir = ROOT / "work" / "tabix_cache"
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(prepare_index, dataset) for dataset in DATASETS]
        for future in as_completed(futures):
            future.result()

    jobs = [
        (dataset, variant)
        for dataset in DATASETS
        for variant in VARIANTS
    ]
    records = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(query_one, dataset, variant, cache_dir): (
                dataset[2],
                dataset[3],
                variant[0],
            )
            for dataset, variant in jobs
        }
        for future in as_completed(futures):
            study, condition, target = futures[future]
            try:
                result = future.result()
                if result is not None:
                    records.append(result)
            except Exception as exc:
                records.append(
                    {
                        "target": target,
                        "study": study,
                        "condition": condition,
                        "query_error": str(exc),
                    }
                )
    result = pd.DataFrame(records)
    if result.empty:
        result = pd.DataFrame(
            columns=[
                "target",
                "study",
                "condition",
                "pvalue",
                "beta",
                "se",
            ]
        )
    result.to_csv(
        table_dir / "experimental_eqtl_catalogue_remaining_variants.csv",
        index=False,
    )
    coverage = []
    for variant in VARIANTS:
        gene, gene_id, rsid, chromosome, position, ref, alt = variant
        subset = result[result["target"].eq(gene)] if "target" in result else result
        nominal = subset[subset["pvalue"] < 0.05] if "pvalue" in subset else subset
        coverage.append(
            {
                "target": gene,
                "rsid": rsid,
                "datasets_queried": len(DATASETS),
                "eQTL_hits_nominal_p05": int(len(nominal)),
                "minimum_pvalue": float(subset["pvalue"].min())
                if len(subset) and "pvalue" in subset
                else None,
                "classification": "nominal_colon_or_myeloid_eQTL"
                if len(nominal)
                else "public_myeloid_eqtl_not_covered",
            }
        )
    coverage_frame = pd.DataFrame(coverage)
    coverage_frame.to_csv(
        table_dir / "experimental_eqtl_catalogue_coverage.csv",
        index=False,
    )
    summary = {
        "datasets_queried": len(DATASETS),
        "variants_queried": len(VARIANTS),
        "dataset_names": [dataset[2] + ":" + dataset[3] for dataset in DATASETS],
        "coverage": coverage_frame.to_dict("records"),
        "note": "Exact-variant nominal eQTL lookup; formal colocalization is only appropriate for pairs with nominal signal.",
    }
    (table_dir / "experimental_eqtl_catalogue_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    if not result.empty:
        print(result.sort_values(["target", "pvalue"]).to_string(index=False))


if __name__ == "__main__":
    main()
