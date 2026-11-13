from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_api_key, load_config


CHUNK_SIZE = 40
POOL_DIR = ROOT / "data" / "processed" / "ibd_conditioning"
POOL_FILES = {
    "active_colitis_macrophage_oliver": (
        POOL_DIR / "active_colitis_macrophage_oliver_256cells.h5ad",
        "Oliver2024",
        "active",
        "colon_macrophage",
    ),
    "healthy_colorectal_macrophage_oliver": (
        POOL_DIR / "healthy_colorectal_macrophage_oliver_256cells.h5ad",
        "Oliver2024",
        "healthy",
        "colon_macrophage",
    ),
    "active_crohn_macrophage_independent": (
        POOL_DIR / "active_crohn_macrophage_independent_256cells.h5ad",
        "Crohn2023",
        "active",
        "colon_macrophage",
    ),
    "healthy_colon_macrophage_crohn2023": (
        POOL_DIR / "healthy_colon_macrophage_crohn2023_256cells.h5ad",
        "Crohn2023",
        "healthy",
        "colon_macrophage",
    ),
    "active_ibd_myeloid_garrido": (
        POOL_DIR / "active_ibd_myeloid_garrido_256cells.h5ad",
        "Garrido2023",
        "active",
        "ibd_myeloid",
    ),
    "healthy_myeloid_garrido": (
        POOL_DIR / "healthy_myeloid_garrido_256cells.h5ad",
        "Garrido2023",
        "healthy",
        "ibd_myeloid",
    ),
}


def main() -> None:
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers

    configure_proxy()
    config = load_config()
    model = dna_client.create(load_api_key(), timeout=30)
    candidates = pd.read_csv(
        ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv"
    )
    primary = candidates[candidates["study_role"] == "primary_ibd"].drop_duplicates(
        "canonical_variant"
    )
    positive = pd.DataFrame(config["positive_controls"])[["variant"]].rename(
        columns={"variant": "canonical_variant"}
    )
    all_variants = (
        pd.concat(
            [primary[["canonical_variant"]], positive],
            ignore_index=True,
        )["canonical_variant"]
        .drop_duplicates()
        .astype(str)
        .tolist()
    )
    lead = pd.read_csv(ROOT / "results" / "tables" / "ranked_lead_variants.csv")
    top_variants = (
        pd.concat(
            [lead.head(28)[["canonical_variant"]], positive],
            ignore_index=True,
        )["canonical_variant"]
        .drop_duplicates()
        .astype(str)
        .tolist()
    )
    condition_variants = {
        "active_colitis_macrophage_oliver": all_variants,
        "healthy_colorectal_macrophage_oliver": all_variants,
        "active_crohn_macrophage_independent": top_variants,
        "healthy_colon_macrophage_crohn2023": top_variants,
        "active_ibd_myeloid_garrido": top_variants,
        "healthy_myeloid_garrido": top_variants,
    }
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[
        "SC_GENE_EXPRESSION_DENOISED"
    ]
    scratch_root = ROOT / "work" / "active_ibd_score_chunks"
    scratch_root.mkdir(parents=True, exist_ok=True)
    chunks = []
    for condition_name, variant_strings in condition_variants.items():
        cells_path, cohort, disease_state, cell_type = POOL_FILES[condition_name]
        cells = ad.read_h5ad(cells_path)
        condition_dir = scratch_root / condition_name
        condition_dir.mkdir(parents=True, exist_ok=True)
        for offset in range(0, len(variant_strings), CHUNK_SIZE):
            chunk_strings = variant_strings[offset : offset + CHUNK_SIZE]
            chunk_path = condition_dir / (
                f"chunk_{offset:05d}_{offset + len(chunk_strings) - 1:05d}.parquet"
            )
            if chunk_path.exists():
                chunks.append(pd.read_parquet(chunk_path))
                continue
            variants = [genome.Variant.from_str(value) for value in chunk_strings]
            intervals = [
                variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
                for variant in variants
            ]
            scores = None
            last_error = None
            for attempt in range(4):
                try:
                    scores = model.score_variants(
                        intervals=intervals,
                        variants=variants,
                        cell_conditionings=[cells] * len(variants),
                        cell_condition_names=[condition_name] * len(variants),
                        variant_scorers=[scorer],
                        progress_bar=True,
                        max_workers=4,
                    )
                    break
                except Exception as exc:
                    last_error = exc
                    if attempt == 3:
                        raise
                    wait_seconds = 10 * (2**attempt)
                    print(
                        f"{condition_name}: transient scoring error; "
                        f"retrying in {wait_seconds}s"
                    )
                    time.sleep(wait_seconds)
                    model = dna_client.create(load_api_key(), timeout=30)
            if scores is None:
                raise RuntimeError(f"Scoring failed: {last_error}")
            tidy = variant_scorers.tidy_scores(scores, filter_protein_coding=True)
            tidy["cohort"] = cohort
            tidy["disease_state"] = disease_state
            tidy["cell_pool_type"] = cell_type
            tidy["is_positive_control"] = tidy["variant_id"].astype(str).eq(
                "chr21:39094644:G>A"
            )
            for column in tidy.select_dtypes(include=["object", "str"]).columns:
                tidy[column] = tidy[column].map(
                    lambda value: value if pd.isna(value) else str(value)
                )
            tidy.to_parquet(chunk_path, index=False)
            chunks.append(tidy)
            print(
                f"{condition_name}: "
                f"{min(offset + len(chunk_strings), len(variant_strings))}/"
                f"{len(variant_strings)}"
            )
    output = pd.concat(chunks, ignore_index=True)
    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    output.to_parquet(
        table_dir / "active_ibd_custom_condition_scores.parquet", index=False
    )
    summary = (
        output.groupby(
            ["cohort", "disease_state", "cell_pool_type", "biosample_name"],
            as_index=False,
        )
        .agg(
            rows=("raw_score", "size"),
            variants=("variant_id", "nunique"),
            genes=("gene_name", "nunique"),
        )
    )
    summary.to_csv(table_dir / "active_ibd_condition_score_summary.csv", index=False)
    print(f"Saved {len(output)} disease-state score rows.")


if __name__ == "__main__":
    main()
