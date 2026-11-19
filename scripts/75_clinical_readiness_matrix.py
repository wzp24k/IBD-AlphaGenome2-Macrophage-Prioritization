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
        perturbation = 1.0 if pd.notna(row.get("experimental_support")) else 0.0
        pqtl_mr_score = 1.0 if float(row.get("pqtl_positive_nominal", 0) or 0) > 0 else 0.0
        pqtl_coloc_score = 1.0 if float(row.get("pqtl_coloc_max_h4", 0) or 0) >= 0.8 else 0.0
        treatment_score = 1.0 if bool(row.get("treatment_response_support", False)) else 0.0
        clinical_stage = 1.0 if pd.notna(row.get("clinical_stage")) else 0.0
        direction_compatible = 1.0 if pd.notna(row.get("directionally_compatible_clinical_actions")) and str(row.get("directionally_compatible_clinical_actions")).strip() not in {"", "nan"} else 0.0
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
                "pqtl_colocalization_score": pqtl_coloc_score,
                "treatment_response_score": treatment_score,
                "clinical_stage_score": clinical_stage,
                "direction_compatible_score": direction_compatible,
                "readiness_score_0_10": score,
                "readiness_class": (
                    "Highest"
                    if score >= 7
                    else "High"
                    if score >= 5
                    else "Moderate"
                    if score >= 3
                    else "Exploratory"
                ),
                "clinical_stage": row.get("clinical_stage", ""),
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
    }
    (TABLE_DIR / "clinical_development_readiness_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
