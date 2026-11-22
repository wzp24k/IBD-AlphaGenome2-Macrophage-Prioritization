from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


TABLE_DIR = ROOT / "results" / "tables"


def as_bool(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "yes", "1"}


def main() -> None:
    clinical = pd.read_csv(TABLE_DIR / "top10_clinical_translation_annotation.csv")
    support = pd.read_csv(TABLE_DIR / "top10_cd_uc_target_support.csv")
    pqtl_mr = pd.read_csv(TABLE_DIR / "pqtl_cis_mr_ibd.csv")
    pqtl_coloc = pd.read_csv(TABLE_DIR / "pqtl_ibd_colocalization.csv")
    treatment = pd.read_csv(TABLE_DIR / "treatment_response_meta_analysis.csv")

    merged = clinical.merge(
        support.drop(columns=["priority_rank", "variant"]),
        on="target",
        how="left",
    )
    pqtl_mr["pqtl_positive_nominal"] = (
        pqtl_mr["beta"].gt(0) & pqtl_mr["p"].lt(0.05)
    )
    pqtl_target = (
        pqtl_mr.groupby(pqtl_mr["pqtl_dataset"].map(lambda _: "PRKCB"))
        .agg(
            pqtl_tests=("p", "size"),
            pqtl_positive_nominal=("pqtl_positive_nominal", "sum"),
            pqtl_min_p=("p", "min"),
        )
        .reset_index()
        .rename(columns={"pqtl_dataset": "target"})
    )
    coloc_target = (
        pqtl_coloc.assign(target="PRKCB")
        .groupby("target")
        .agg(
            pqtl_coloc_tests=("H4", "size"),
            pqtl_coloc_max_h4=("H4", "max"),
        )
        .reset_index()
    )
    treatment_target = treatment[treatment["feature"].eq("PRKCB_network_program")].copy()
    treatment_target["target"] = "PRKCB"
    treatment_target["treatment_response_support"] = (
        treatment_target["ci_low"].gt(0) | treatment_target["ci_high"].lt(0)
    )
    treatment_target = treatment_target[
        ["target", "smd", "ci_low", "ci_high", "p", "treatment_response_support"]
    ]
    merged = (
        merged.merge(pqtl_target, on="target", how="left")
        .merge(coloc_target, on="target", how="left")
        .merge(treatment_target, on="target", how="left")
    )

    for column in [
        "pqtl_tests",
        "pqtl_positive_nominal",
        "pqtl_coloc_tests",
        "pqtl_coloc_max_h4",
        "smd",
        "ci_low",
        "ci_high",
        "p",
        "treatment_response_support",
    ]:
        merged[column] = merged[column].fillna(0 if column != "pqtl_coloc_max_h4" else np.nan)

    rows = []
    for _, row in merged.iterrows():
        genetic_transport = min(float(row.get("external_clusters", 0)), 1.0)
        disease_specific = 1.0 if bool(row.get("cd_uc_both", False)) else 0.5 if bool(row.get("cd_uc_any", False)) else 0.0
        exact_eqtl = 1.0 if as_bool(row.get("exact_variant_eqtl", False)) else 0.0
        coloc = 1.0 if as_bool(row.get("formal_colocalization", False)) else 0.0
        perturbation_support = str(row.get("experimental_support", "")).strip().lower()
        perturbation = (
            1.0
            if perturbation_support in {
                "perturbseq_transcriptome",
                "target_engagement",
            }
            else 0.0
        )
        pqtl_positive = float(row.get("pqtl_positive_nominal", 0) or 0) > 0
        pqtl_min_p = float(row.get("pqtl_min_p", np.nan) or np.nan)
        if pqtl_positive and np.isfinite(pqtl_min_p) and pqtl_min_p <= 5e-8:
            pqtl_mr_score = 1.0
            pqtl_evidence_class = "Strong"
        elif pqtl_positive:
            pqtl_mr_score = 0.25
            pqtl_evidence_class = "Exploratory"
        else:
            pqtl_mr_score = 0.0
            pqtl_evidence_class = "None"
        pqtl_coloc_score = 1.0 if float(row.get("pqtl_coloc_max_h4", 0) or 0) >= 0.8 else 0.0
        treatment_score = 1.0 if bool(row.get("treatment_response_support", False)) else 0.0
        ibd_stage_value = row.get("ibd_indication_clinical_stage")
        ibd_stage_text = (
            "" if pd.isna(ibd_stage_value) else str(ibd_stage_value).strip()
        )
        clinical_stage = 1.0 if ibd_stage_text not in {"", "nan", "none"} else 0.0
        direction_compatible = 1.0 if as_bool(
            row.get("directionally_compatible_ibd_drug", False)
        ) else 0.0
        score = sum(
            [
                genetic_transport,
                disease_specific,
                exact_eqtl,
                coloc,
                perturbation,
                pqtl_mr_score,
                pqtl_coloc_score,
                treatment_score,
                clinical_stage,
                direction_compatible,
            ]
        )
        evidence_components = [
            genetic_transport,
            disease_specific,
            exact_eqtl,
            coloc,
            perturbation,
            pqtl_mr_score,
            pqtl_coloc_score,
            treatment_score,
        ]
        evidence_only = sum(evidence_components) / len(evidence_components) * 10
        clinical_weighted = (
            sum(evidence_components) + 2 * clinical_stage + 2 * direction_compatible
        ) / 12 * 10
        no_treatment = (
            sum(evidence_components) - treatment_score
        ) / (len(evidence_components) - 1) * 10
        rows.append(
            {
                "priority_rank": int(row["priority_rank"]),
                "target": row["target"],
                "genetic_transport_score": genetic_transport,
                "disease_specific_score": disease_specific,
                "exact_eqtl_score": exact_eqtl,
                "colocalization_score": coloc,
                "public_perturbation_score": perturbation,
                "pqtl_mr_score": pqtl_mr_score,
                "pqtl_evidence_class": pqtl_evidence_class,
                "pqtl_colocalization_score": pqtl_coloc_score,
                "treatment_response_score": treatment_score,
                "clinical_stage_score": clinical_stage,
                "direction_compatible_score": direction_compatible,
                "readiness_score_0_10": score,
                "evidence_only_score": evidence_only,
                "clinical_weighted_score": clinical_weighted,
                "no_treatment_score": no_treatment,
                "readiness_class": (
                    "Highest"
                    if score >= 7
                    else "High"
                    if score >= 5
                    else "Moderate"
                    if score >= 3
                    else "Exploratory"
                ),
                "ibd_indication_clinical_stage": ibd_stage_text,
                "any_indication_clinical_stage": row.get("clinical_stage", ""),
                "ibd_indication_drugs": row.get("ibd_indication_drugs", ""),
                "risk_effect_on_target": row.get("risk_effect_on_target", ""),
                "translational_priority": row.get("translational_priority", ""),
            }
        )

    output = pd.DataFrame(rows).sort_values(
        ["readiness_score_0_10", "priority_rank"],
        ascending=[False, True],
    )
    output.to_csv(
        TABLE_DIR / "clinical_development_readiness_matrix.csv",
        index=False,
    )
    for column in [
        "readiness_score_0_10",
        "evidence_only_score",
        "clinical_weighted_score",
        "no_treatment_score",
    ]:
        output[f"{column}_rank"] = output[column].rank(
            ascending=False,
            method="min",
        ).astype(int)
    rank_columns = [
        "readiness_score_0_10_rank",
        "evidence_only_score_rank",
        "clinical_weighted_score_rank",
        "no_treatment_score_rank",
    ]
    rank_sensitivity = output[
        ["priority_rank", "target", *rank_columns]
    ].copy()
    rank_sensitivity["rank_min"] = rank_sensitivity[rank_columns].min(axis=1)
    rank_sensitivity["rank_max"] = rank_sensitivity[rank_columns].max(axis=1)
    rank_sensitivity.to_csv(
        TABLE_DIR / "clinical_readiness_rank_sensitivity.csv",
        index=False,
    )
    prkcb_ranks = rank_sensitivity.loc[
        rank_sensitivity["target"].eq("PRKCB"),
        rank_columns,
    ].iloc[0]
    pdgfb_ranks = rank_sensitivity.loc[
        rank_sensitivity["target"].eq("PDGFB"),
        rank_columns,
    ].iloc[0]
    summary = {
        "targets": len(output),
        "highest_readiness": output.loc[
            output["readiness_class"].eq("Highest"), "target"
        ].tolist(),
        "top_score": float(output["readiness_score_0_10"].max()),
        "top_target": output.iloc[0]["target"],
        "pqtl_mr_support": output.loc[
            output["pqtl_mr_score"].gt(0), "target"
        ].tolist(),
        "pqtl_coloc_support": output.loc[
            output["pqtl_colocalization_score"].gt(0), "target"
        ].tolist(),
        "treatment_response_support": output.loc[
            output["treatment_response_score"].gt(0), "target"
        ].tolist(),
        "clinical_stage_targets": output.loc[
            output["clinical_stage_score"].gt(0), "target"
        ].tolist(),
        "prkcb_rank_min": int(prkcb_ranks.min()),
        "prkcb_rank_max": int(prkcb_ranks.max()),
        "pdgfb_rank_min": int(pdgfb_ranks.min()),
        "pdgfb_rank_max": int(pdgfb_ranks.max()),
        "ranking_schemes": len(rank_columns),
    }
    (TABLE_DIR / "clinical_development_readiness_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
