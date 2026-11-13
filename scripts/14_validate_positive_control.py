from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_api_key


def main() -> None:
    from alphagenome.atlas import atlas
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers

    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    live = pd.read_parquet(table_dir / "live_condition_gene_scores.parquet")
    ets2 = live[
        (live["variant_id"].astype(str) == "chr21:39094644:G>A")
        & (live["gene_name"] == "ETS2")
        & live["imputation_mode"].isin(["observed", "imputed_supported"])
    ]
    summary = (
        ets2.groupby(["output_type", "biosample_name", "condition_role"], as_index=False)
        .agg(raw_score=("raw_score", "mean"), tracks=("raw_score", "size"))
    )
    summary.to_csv(table_dir / "positive_control_ets2_live.csv", index=False)

    key = load_api_key()
    atlas_client = atlas.create(key, timeout=30)
    variant = genome.Variant.from_str("chr21:39094644:G>A")
    atlas_scores = atlas_client.query_variant(
        variant, requested_scorers=["RNA_SEQ"]
    )["RNA_SEQ"]
    ets2_rows = np.where(atlas_scores.obs["gene_name"].to_numpy() == "ETS2")[0]
    row = int(ets2_rows[0])
    raw = np.asarray(atlas_scores.X[row])
    quantile = np.asarray(atlas_scores.layers["quantiles"][row])
    max_quantile_index = int(np.argmax(np.abs(quantile)))
    atlas_result = pd.DataFrame(
        [
            {
                "variant": "chr21:39094644:G>A",
                "rsid": "rs2836882",
                "target_gene": "ETS2",
                "atlas_max_abs_quantile": abs(quantile[max_quantile_index]),
                "atlas_signed_quantile": quantile[max_quantile_index],
                "atlas_raw_score": raw[max_quantile_index],
                "atlas_biosample": atlas_scores.var.iloc[max_quantile_index][
                    "biosample_name"
                ],
            }
        ]
    )
    atlas_result.to_csv(
        table_dir / "positive_control_ets2_atlas.csv", index=False
    )

    custom_path = table_dir / "positive_control_custom_ets2.csv"
    if not custom_path.exists():
        cells = ad.read_h5ad(
            ROOT
            / "data"
            / "processed"
            / "healthy_colorectal_macrophage_256cells.h5ad"
        )
        model = dna_client.create(key, timeout=30)
        interval = variant.reference_interval.resize(
            dna_client.SEQUENCE_LENGTH_1MB
        )
        score = model.score_variant(
            interval=interval,
            variant=variant,
            cell_conditioning=cells,
            cell_condition_name="healthy_colorectal_macrophage",
            variant_scorers=[
                variant_scorers.RECOMMENDED_VARIANT_SCORERS[
                    "SC_GENE_EXPRESSION_DENOISED"
                ]
            ],
        )
        custom = variant_scorers.tidy_scores([score], filter_protein_coding=True)
        custom = custom[custom["gene_name"] == "ETS2"]
        custom.to_csv(custom_path, index=False)
    else:
        custom = pd.read_csv(custom_path)

    def value(output_type: str, biosample: str) -> float:
        selected = summary[
            (summary["output_type"] == output_type)
            & (summary["biosample_name"] == biosample)
        ]
        return float(selected["raw_score"].iloc[0])

    checks = {
        "atlas_quantile_above_0_99": bool(
            atlas_result["atlas_max_abs_quantile"].iloc[0] >= 0.99
        ),
        "colon_macrophage_rna_down": value("RNA_SEQ", "colon macrophage") <= -0.05,
        "colon_macrophage_sc_down": value(
            "SC_GENE_EXPRESSION", "colon macrophage"
        ) <= -0.05,
        "colon_t_cell_effect_near_zero": abs(
            value("RNA_SEQ", "colon T cell")
        ) <= 0.02,
        "colon_epithelial_effect_near_zero": abs(
            value("RNA_SEQ", "colon epithelial cell")
        ) <= 0.02,
        "custom_macrophage_down": float(custom["raw_score"].mean()) <= -0.05,
    }
    result = {
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "colon_macrophage_effects": {
            "rna_seq": value("RNA_SEQ", "colon macrophage"),
            "single_cell": value("SC_GENE_EXPRESSION", "colon macrophage"),
        },
        "custom_colorectal_macrophage_effect": float(custom["raw_score"].mean()),
        "atlas_max_abs_quantile": float(
            atlas_result["atlas_max_abs_quantile"].iloc[0]
        ),
    }
    (table_dir / "positive_control_validation.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    if not result["all_checks_pass"]:
        raise SystemExit("Positive control validation failed.")


if __name__ == "__main__":
    main()
