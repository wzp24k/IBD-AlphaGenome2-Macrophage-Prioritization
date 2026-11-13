from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


def sign_with_tolerance(values: pd.Series, tolerance: float = 0.001) -> pd.Series:
    result = np.sign(values.fillna(0))
    result[values.abs() < tolerance] = 0
    return result


def condition_effects(scores: pd.DataFrame) -> pd.DataFrame:
    # Cell-conditioning-only outputs have an empty imputation label.
    usable = scores.copy()
    grouped = (
        usable.groupby(
            [
                "biosample_name",
                "cohort",
                "disease_state",
                "cell_pool_type",
                "variant_id",
                "gene_name",
            ],
            as_index=False,
        )
        .agg(raw_score=("raw_score", "mean"), track_count=("raw_score", "size"))
    )
    grouped["abs_score"] = grouped["raw_score"].abs()
    grouped["abs_effect_percentile"] = grouped.groupby("biosample_name")[
        "abs_score"
    ].rank(pct=True)
    return grouped


def merge_target_effects(
    effects: pd.DataFrame,
    candidates: pd.DataFrame,
    condition: str,
    prefix: str,
) -> pd.DataFrame:
    condition = effects[effects["biosample_name"] == condition]
    target = condition.merge(
        candidates[["canonical_variant", "l2g_target_symbol"]],
        left_on=["variant_id", "gene_name"],
        right_on=["canonical_variant", "l2g_target_symbol"],
        how="inner",
    )
    return target[
        [
            "canonical_variant",
            "raw_score",
            "abs_effect_percentile",
        ]
    ].rename(
        columns={
            "raw_score": f"{prefix}_effect",
            "abs_effect_percentile": f"{prefix}_abs_effect_percentile",
        }
    )


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    score_path = table_dir / "active_ibd_custom_condition_scores.parquet"
    if not score_path.exists():
        raise FileNotFoundError("Run scripts/17_score_active_ibd_conditions.py first.")
    scores = pd.read_parquet(score_path)
    effects = condition_effects(scores)
    ranked = pd.read_csv(table_dir / "ranked_variant_gene_candidates.csv")
    lead = pd.read_csv(table_dir / "ranked_lead_variants.csv")
    top28 = lead.head(28).copy()

    active_oliver = merge_target_effects(
        effects,
        ranked,
        "active_colitis_macrophage_oliver",
        "oliver_active",
    )
    healthy_oliver = merge_target_effects(
        effects,
        ranked,
        "healthy_colorectal_macrophage_oliver",
        "oliver_healthy",
    )
    updated = ranked.merge(active_oliver, on="canonical_variant", how="left")
    updated = updated.merge(healthy_oliver, on="canonical_variant", how="left")
    updated["oliver_disease_state_delta"] = (
        updated["oliver_active_effect"] - updated["oliver_healthy_effect"]
    )
    updated["oliver_active_same_direction"] = (
        sign_with_tolerance(updated["oliver_active_effect"])
        == sign_with_tolerance(updated["oliver_healthy_effect"])
    )

    cohort_conditions = {
        "crohn": {
            "active": "active_crohn_macrophage_independent",
            "healthy": "healthy_colon_macrophage_crohn2023",
        },
        "garrido": {
            "active": "active_ibd_myeloid_garrido",
            "healthy": "healthy_myeloid_garrido",
        },
    }
    top_effects = top28[["canonical_variant", "l2g_target_symbol"]].copy()
    for short_name, conditions in cohort_conditions.items():
        active = merge_target_effects(
            effects,
            top28,
            conditions["active"],
            f"{short_name}_active",
        )
        healthy = merge_target_effects(
            effects,
            top28,
            conditions["healthy"],
            f"{short_name}_healthy",
        )
        top_effects = top_effects.merge(active, on="canonical_variant", how="left")
        top_effects = top_effects.merge(healthy, on="canonical_variant", how="left")
        top_effects[f"{short_name}_disease_state_delta"] = (
            top_effects[f"{short_name}_active_effect"]
            - top_effects[f"{short_name}_healthy_effect"]
        )
    updated = updated.merge(top_effects.drop(columns=["l2g_target_symbol"]), on="canonical_variant", how="left")

    active_sign = sign_with_tolerance(updated["oliver_active_effect"])
    chronic_sign = sign_with_tolerance(updated["crohn_active_effect"])
    garrido_sign = sign_with_tolerance(updated["garrido_active_effect"])
    baseline_sign = sign_with_tolerance(updated["custom_macrophage_effect"])
    updated["active_oliver_baseline_concordant"] = (
        (active_sign == baseline_sign) & (active_sign != 0)
    )
    updated["crohn_active_oliver_concordant"] = (
        (chronic_sign == active_sign) & (chronic_sign != 0)
    )
    updated["garrido_active_oliver_concordant"] = (
        (garrido_sign == active_sign) & (garrido_sign != 0)
    )
    updated["independent_same_direction_count"] = (
        updated["crohn_active_oliver_concordant"].astype(int)
        + updated["garrido_active_oliver_concordant"].astype(int)
    )
    updated["active_disease_multicohort_supported"] = (
        updated["active_oliver_baseline_concordant"]
        & (updated["independent_same_direction_count"] >= 1)
    )
    updated["active_disease_strictly_replicated"] = (
        updated["active_oliver_baseline_concordant"]
        & (updated["independent_same_direction_count"] == 2)
    )
    updated["max_active_effect_magnitude"] = updated[
        ["oliver_active_effect", "crohn_active_effect", "garrido_active_effect"]
    ].abs().max(axis=1, skipna=True)
    updated["activity_state_class"] = np.select(
        [
            updated["active_disease_multicohort_supported"].fillna(False),
            updated["max_active_effect_magnitude"] < 0.005,
        ],
        ["persistent", "attenuated"],
        default="discordant",
    )
    updated["final_evidence_tier"] = np.where(
        (updated["evidence_tier"] == "Tier 1")
        & updated["active_disease_multicohort_supported"].fillna(False),
        "Tier 1",
        np.where(
            (updated["evidence_tier"] == "Tier 1")
            & (updated["max_active_effect_magnitude"] < 0.005),
            "Tier 1 attenuated",
            np.where(
                (updated["evidence_tier"] == "Tier 1"),
                "Tier 2",
                updated["evidence_tier"],
            ),
        ),
    )

    updated.to_csv(
        table_dir / "ranked_variant_gene_candidates_active_ibd.csv", index=False
    )
    lead_updated = lead.merge(
        updated[
            [
                "canonical_variant",
                "oliver_active_effect",
                "oliver_healthy_effect",
                "oliver_active_abs_effect_percentile",
                "oliver_disease_state_delta",
                "oliver_active_same_direction",
                "crohn_active_effect",
                "crohn_healthy_effect",
                "crohn_active_abs_effect_percentile",
                "crohn_disease_state_delta",
                "garrido_active_effect",
                "garrido_healthy_effect",
                "garrido_active_abs_effect_percentile",
                "garrido_disease_state_delta",
                "active_oliver_baseline_concordant",
                "crohn_active_oliver_concordant",
                "garrido_active_oliver_concordant",
                "independent_same_direction_count",
                "active_disease_multicohort_supported",
                "active_disease_strictly_replicated",
                "max_active_effect_magnitude",
                "activity_state_class",
                "final_evidence_tier",
            ]
        ],
        on="canonical_variant",
        how="left",
    )
    lead_updated.to_csv(
        table_dir / "ranked_lead_variants_active_ibd.csv", index=False
    )

    common = lead_updated.dropna(
        subset=["oliver_active_effect", "oliver_healthy_effect"]
    )
    rho_all, p_all = spearmanr(
        common["oliver_active_effect"],
        common["oliver_healthy_effect"],
        nan_policy="omit",
    )
    top = lead_updated.head(28)
    rho_top, p_top = spearmanr(
        top["oliver_active_effect"],
        top["oliver_healthy_effect"],
        nan_policy="omit",
    )
    summary = {
        "scored_variants": int(scores["variant_id"].nunique()),
        "lead_targets_with_oliver_active": int(
            common["oliver_active_effect"].notna().sum()
        ),
        "lead_active_healthy_direction_concordance": float(
            common["oliver_active_same_direction"].mean()
        ),
        "spearman_active_vs_healthy_all_leads": {
            "rho": float(rho_all),
            "p": float(p_all),
            "n": int(len(common)),
        },
        "spearman_active_vs_healthy_top28": {
            "rho": float(rho_top),
            "p": float(p_top),
            "n": int(len(top)),
        },
        "top28_multicohort_supported": int(
            top["active_disease_multicohort_supported"].fillna(False).sum()
        ),
        "top28_strictly_replicated": int(
            top["active_disease_strictly_replicated"].fillna(False).sum()
        ),
        "tier1_active_multicohort_supported": int(
            (
                (top["evidence_tier"] == "Tier 1")
                & top["active_disease_multicohort_supported"].fillna(False)
            ).sum()
        ),
        "final_validated_tier1": int(
            top["final_evidence_tier"].isin(
                ["Tier 1", "Tier 1 attenuated"]
            ).sum()
        ),
        "final_persistent_tier1": int(
            (top["final_evidence_tier"] == "Tier 1").sum()
        ),
        "final_attenuated_tier1": int(
            (top["final_evidence_tier"] == "Tier 1 attenuated").sum()
        ),
    }
    (table_dir / "active_ibd_validation_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print("\nTop 20 disease-state validation:")
    columns = [
        "lead_rank",
        "l2g_target_symbol",
        "evidence_tier",
        "priority_score",
        "oliver_active_effect",
        "oliver_healthy_effect",
        "crohn_active_effect",
        "crohn_healthy_effect",
        "garrido_active_effect",
        "garrido_healthy_effect",
        "independent_same_direction_count",
        "active_disease_multicohort_supported",
    ]
    print(lead_updated[columns].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
