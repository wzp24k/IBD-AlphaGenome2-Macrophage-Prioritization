from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_api_key


CHUNK_SIZE = 40


def main() -> None:
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers

    configure_proxy()
    key = load_api_key()
    model = dna_client.create(key, timeout=30)
    candidates = pd.read_csv(
        ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv"
    )
    primary = candidates[candidates["study_role"] == "primary_ibd"].drop_duplicates(
        "canonical_variant"
    )
    variant_strings = primary["canonical_variant"].astype(str).tolist()
    cells = ad.read_h5ad(
        ROOT
        / "data"
        / "processed"
        / "healthy_colorectal_macrophage_256cells.h5ad"
    )
    condition_name = "healthy_colorectal_macrophage"
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[
        "SC_GENE_EXPRESSION_DENOISED"
    ]

    scratch_dir = ROOT / "work" / "custom_macrophage_chunks"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    chunks = []
    for offset in range(0, len(variant_strings), CHUNK_SIZE):
        chunk_strings = variant_strings[offset : offset + CHUNK_SIZE]
        chunk_path = scratch_dir / (
            f"custom_{offset:05d}_{offset + len(chunk_strings) - 1:05d}.parquet"
        )
        if chunk_path.exists():
            chunks.append(pd.read_parquet(chunk_path))
            continue
        variants = [genome.Variant.from_str(value) for value in chunk_strings]
        intervals = [
            variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
            for variant in variants
        ]
        scores = model.score_variants(
            intervals=intervals,
            variants=variants,
            cell_conditionings=[cells] * len(variants),
            cell_condition_names=[condition_name] * len(variants),
            variant_scorers=[scorer],
            progress_bar=True,
            max_workers=4,
        )
        tidy = variant_scorers.tidy_scores(scores, filter_protein_coding=True)
        for column in tidy.select_dtypes(include=["object", "str"]).columns:
            tidy[column] = tidy[column].map(
                lambda value: value if pd.isna(value) else str(value)
            )
        tidy.to_parquet(chunk_path, index=False)
        chunks.append(tidy)
        print(
            f"Custom macrophage scoring: completed "
            f"{min(offset + len(chunk_strings), len(variant_strings))}/"
            f"{len(variant_strings)} variants"
        )
    all_scores = pd.concat(chunks, ignore_index=True)
    output_dir = ROOT / "results" / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    all_scores.to_parquet(
        output_dir / "custom_colorectal_macrophage_gene_scores.parquet",
        index=False,
    )
    print(f"Saved {len(all_scores)} custom-condition score rows.")


if __name__ == "__main__":
    main()
