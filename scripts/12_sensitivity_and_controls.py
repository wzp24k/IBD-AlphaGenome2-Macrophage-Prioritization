from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    lead = pd.read_csv(table_dir / "ranked_lead_variants.csv")
    ranked = pd.read_csv(table_dir / "ranked_variant_gene_candidates.csv")

    grid = []
    for atlas_threshold in [0.99, 0.975, 0.95, 0.90]:
        for specificity_threshold in [0.10, 0.25, 0.40]:
            selected = lead[
                (lead["atlas_gi_max_abs_q"] >= atlas_threshold)
                & (
                    (lead["rna_specificity"] >= specificity_threshold)
                    | (lead["sc_specificity"] >= specificity_threshold)
                )
                & (lead["custom_macrophage_abs_percentile"] >= 0.80)
            ]
            grid.append(
                {
                    "atlas_quantile_threshold": atlas_threshold,
                    "specificity_delta_threshold": specificity_threshold,
                    "candidate_loci": len(selected),
                    "candidate_genes": selected["l2g_target_symbol"].nunique(),
                }
            )
    pd.DataFrame(grid).to_csv(
        table_dir / "threshold_sensitivity.csv", index=False
    )

    topk = []
    for k in [3, 14, 28]:
        sets = {}
        for score in [
            "priority_score",
            "genetics_weighted_score",
            "macrophage_weighted_score",
        ]:
            sets[score] = set(
                lead.sort_values(score, ascending=False).head(k)["canonical_variant"]
            )
        topk.append(
            {
                "k": k,
                "fraction_of_lead_loci": k / len(lead),
                "overlap_priority_genetics": len(
                    sets["priority_score"] & sets["genetics_weighted_score"]
                ),
                "overlap_priority_macrophage": len(
                    sets["priority_score"] & sets["macrophage_weighted_score"]
                ),
                "overlap_genetics_macrophage": len(
                    sets["genetics_weighted_score"] & sets["macrophage_weighted_score"]
                ),
            }
        )
    pd.DataFrame(topk).to_csv(
        table_dir / "score_topk_overlap.csv", index=False
    )

    ets2_variant = "chr21:39094644:G>A"
    il2ra_variant = "chr10:6060049:C>T"
    live = pd.read_parquet(table_dir / "live_condition_gene_scores.parquet")
    live = live[
        live["imputation_mode"].isin(["observed", "imputed_supported"])
    ]
    controls = []
    for label, variant, gene, expected in [
        (
            "positive_ETS2",
            ets2_variant,
            "ETS2",
            "macrophage-selective down-regulation",
        ),
        (
            "non_macrophage_IL2RA",
            il2ra_variant,
            "IL2RA",
            "strong Atlas score without macrophage-selective live effect",
        ),
    ]:
        subset = live[
            (live["variant_id"].astype(str) == variant)
            & (live["gene_name"] == gene)
        ]
        summary = (
            subset.groupby(
                ["output_type", "biosample_name", "condition_role"],
                as_index=False,
            )
            .agg(raw_score=("raw_score", "mean"), tracks=("raw_score", "size"))
        )
        for row in summary.itertuples():
            controls.append(
                {
                    "control": label,
                    "expected_pattern": expected,
                    "variant": variant,
                    "target_gene": gene,
                    "output_type": row.output_type,
                    "biosample_name": row.biosample_name,
                    "condition_role": row.condition_role,
                    "raw_score": row.raw_score,
                    "tracks": row.tracks,
                }
            )
    custom = pd.read_csv(table_dir / "positive_control_custom_ets2.csv")
    controls.append(
        {
            "control": "positive_ETS2",
            "expected_pattern": "macrophage-selective down-regulation",
            "variant": ets2_variant,
            "target_gene": "ETS2",
            "output_type": "SC_GENE_EXPRESSION_DENOISED",
            "biosample_name": "healthy_colorectal_macrophage",
            "condition_role": "custom_primary_macrophage",
            "raw_score": custom["raw_score"].mean(),
            "tracks": len(custom),
        }
    )
    pd.DataFrame(controls).to_csv(
        table_dir / "control_variant_effects.csv", index=False
    )

    il2ra = ranked[ranked["canonical_variant"] == il2ra_variant].iloc[0]
    pd.DataFrame(
        [
            {
                "variant": il2ra_variant,
                "rsid": il2ra["rsIds"],
                "target_gene": il2ra["l2g_target_symbol"],
                "pip": il2ra["pip"],
                "l2g_score": il2ra["l2g_score"],
                "atlas_gi_max_abs_q": il2ra["atlas_gi_max_abs_q"],
                "rna_macrophage_effect": il2ra["rna_macrophage_effect"],
                "rna_control_abs_max": il2ra["rna_control_abs_max"],
                "sc_macrophage_effect": il2ra["sc_macrophage_effect"],
                "sc_control_abs_max": il2ra["sc_control_abs_max"],
                "custom_macrophage_effect": il2ra["custom_macrophage_effect"],
                "custom_macrophage_abs_percentile": il2ra[
                    "custom_macrophage_abs_percentile"
                ],
                "evidence_tier": il2ra["evidence_tier"],
            }
        ]
    ).to_csv(table_dir / "negative_control_il2ra.csv", index=False)
    print(pd.DataFrame(grid).to_string(index=False))
    print()
    print(pd.DataFrame(topk).to_string(index=False))


if __name__ == "__main__":
    main()
