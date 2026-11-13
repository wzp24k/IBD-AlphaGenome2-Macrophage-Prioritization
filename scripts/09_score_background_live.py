from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_api_key, load_config


CHUNK_SIZE = 40
CONTROL_ROLES = [
    "primary_target",
    "immune_control",
    "epithelial_control",
    "epithelial_sensitivity",
]


def main() -> None:
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers

    configure_proxy()
    config = load_config()
    key = load_api_key()
    model = dna_client.create(key, timeout=30)
    background = pd.read_csv(
        ROOT / "data" / "processed" / "frequency_matched_background.csv"
    )
    variant_strings = background["canonical_variant"].astype(str).tolist()
    conditions = {
        name: value
        for name, value in config["conditions"].items()
        if value["role"] in CONTROL_ROLES
    }
    condition_names = list(conditions)
    condition_curies = [conditions[name]["curie"] for name in condition_names]
    condition_role = {
        conditions[name]["curie"]: conditions[name]["role"] for name in condition_names
    }
    scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS[key]
        for key in config["alphascores"]["live_outputs"]
    ]
    scratch_dir = ROOT / "work" / "background_score_chunks_v2"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    chunks = []
    for offset in range(0, len(variant_strings), CHUNK_SIZE):
        chunk_strings = variant_strings[offset : offset + CHUNK_SIZE]
        chunk_path = scratch_dir / (
            f"background_{offset:05d}_{offset + len(chunk_strings) - 1:05d}.parquet"
        )
        if chunk_path.exists():
            chunks.append(pd.read_parquet(chunk_path))
            continue
        variants = [genome.Variant.from_str(value) for value in chunk_strings]
        intervals = [
            variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
            for variant in variants
        ]
        batch_intervals = []
        batch_variants = []
        batch_conditions = []
        for variant, interval in zip(variants, intervals):
            for curie in condition_curies:
                batch_variants.append(variant)
                batch_intervals.append(interval)
                batch_conditions.append(curie)
        scores = model.score_variants(
            intervals=batch_intervals,
            variants=batch_variants,
            conditions=batch_conditions,
            variant_scorers=scorers,
            progress_bar=True,
            max_workers=4,
        )
        tidy = variant_scorers.tidy_scores(scores, filter_protein_coding=True)
        tidy["condition_role"] = tidy["biosample_curie"].map(condition_role)
        for column in tidy.select_dtypes(include=["object", "str"]).columns:
            tidy[column] = tidy[column].map(
                lambda value: value if pd.isna(value) else str(value)
            )
        tidy.to_parquet(chunk_path, index=False)
        chunks.append(tidy)
        print(
            f"Background scoring: completed "
            f"{min(offset + len(chunk_strings), len(variant_strings))}/"
            f"{len(variant_strings)} variants"
        )
    all_scores = pd.concat(chunks, ignore_index=True)
    output_dir = ROOT / "results" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    all_scores.to_parquet(
        output_dir / "background_live_condition_gene_scores.parquet", index=False
    )
    print(f"Saved {len(all_scores)} background score rows.")


if __name__ == "__main__":
    main()
