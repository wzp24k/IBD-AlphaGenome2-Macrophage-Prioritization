from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_api_key, load_config


CHUNK_SIZE = 40


def main() -> None:
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers

    configure_proxy()
    config = load_config()
    key = load_api_key()
    model = dna_client.create(key, timeout=30)
    input_path = ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv"
    candidates = pd.read_csv(input_path)
    primary = candidates[candidates["study_role"] == "primary_ibd"].drop_duplicates(
        "canonical_variant"
    )
    positive_controls = pd.DataFrame(config["positive_controls"])
    variants_df = pd.concat(
        [
            primary[["canonical_variant", "variant_id"]].rename(
                columns={"canonical_variant": "variant"}
            ),
            positive_controls[["variant"]],
        ],
        ignore_index=True,
    ).drop_duplicates("variant")
    variants_df["is_positive_control"] = variants_df["variant"].isin(
        positive_controls["variant"]
    )
    variant_strings = variants_df["variant"].astype(str).tolist()
    conditions = config["conditions"]
    condition_names = list(conditions)
    condition_curies = [conditions[name]["curie"] for name in condition_names]
    condition_role = {
        conditions[name]["curie"]: conditions[name]["role"] for name in condition_names
    }
    scorer_keys = config["alphascores"]["live_outputs"]
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[key] for key in scorer_keys]

    scratch_dir = ROOT / "work" / "live_score_chunks"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    chunks: list[pd.DataFrame] = []
    for offset in range(0, len(variant_strings), CHUNK_SIZE):
        chunk_strings = variant_strings[offset : offset + CHUNK_SIZE]
        chunk_path = scratch_dir / (
            f"live_{offset:05d}_{offset + len(chunk_strings) - 1:05d}.parquet"
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
        tidy["is_positive_control"] = tidy["variant_id"].isin(
            positive_controls["variant"]
        )
        for column in tidy.select_dtypes(include=["object"]).columns:
            tidy[column] = tidy[column].map(
                lambda value: value if pd.isna(value) else str(value)
            )
        tidy.to_parquet(chunk_path, index=False)
        chunks.append(tidy)
        print(
            f"Live scoring: completed "
            f"{min(offset + len(chunk_strings), len(variant_strings))}/"
            f"{len(variant_strings)} variants"
        )

    all_scores = pd.concat(chunks, ignore_index=True)
    output_dir = ROOT / "results" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    all_scores.to_parquet(output_dir / "live_condition_gene_scores.parquet", index=False)
    summary = (
        all_scores.groupby(
            ["output_type", "biosample_name", "biosample_curie", "condition_role"],
            dropna=False,
        )
        .agg(
            rows=("raw_score", "size"),
            variants=("variant_id", "nunique"),
            genes=("gene_name", "nunique"),
        )
        .reset_index()
    )
    summary.to_csv(output_dir / "live_condition_score_summary.csv", index=False)
    metadata = {
        "variants": len(variant_strings),
        "conditions": condition_names,
        "scorers": scorer_keys,
        "rows": len(all_scores),
    }
    (output_dir / "live_score_run_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )
    print(f"Saved {len(all_scores)} gene-level score rows.")


if __name__ == "__main__":
    main()
