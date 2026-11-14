from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


URL = "https://api.platform.opentargets.org/api/v4/graphql"
QUERY = """
query TargetTranslation($id: String!) {
  target(ensemblId: $id) {
    id
    approvedSymbol
    approvedName
    tractability {
      label
      modality
      value
    }
    pathways {
      pathway
      pathwayId
      topLevelTerm
    }
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          name
          drugType
          maximumClinicalStage
          indications {
            count
            rows {
              maxClinicalStage
              disease {
                id
                name
              }
            }
          }
          mechanismsOfAction {
            rows {
              mechanismOfAction
              actionType
              targetName
              targets {
                id
                approvedSymbol
              }
            }
          }
        }
      }
    }
    associatedDiseases(page: {index: 0, size: 20}) {
      count
      rows {
        disease {
          id
          name
        }
        score
      }
    }
  }
}
"""
IBD_NAMES = {
    "inflammatory bowel disease",
    "crohn disease",
    "crohn's disease",
    "ulcerative colitis",
}
PATHWAY_TERMS = re.compile(
    r"TNF|chemokine|PDGF|MAPK|NF-kappa|PI3K|RAB|interleukin|"
    r"immune|cytokine|integrin|JAK|STAT",
    re.IGNORECASE,
)
STAGE_ORDER = {
    "APPROVAL": 6,
    "PHASE_4": 5,
    "PHASE_3": 4,
    "PHASE_2": 3,
    "PHASE_1": 2,
    "PRECLINICAL": 1,
    "": 0,
}


def stage_number(stage: str) -> int:
    return STAGE_ORDER.get(str(stage).upper(), 0)


def query_target(session: requests.Session, target_id: str) -> dict:
    response = session.post(
        URL,
        json={"query": QUERY, "variables": {"id": target_id}},
        timeout=60,
    )
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError("; ".join(error["message"] for error in payload["errors"]))
    return payload["data"]["target"]


def clinical_translation_status(
    ibd_indication: bool,
    clinical_stage_number: int,
    tractable: bool,
    validated: bool,
) -> tuple[str, str]:
    if ibd_indication:
        return (
            "Clinical-stage IBD-indication signal",
            "A clinical-stage drug lists an IBD, Crohn or ulcerative colitis indication; treatment use still requires directionality and trial-level evidence.",
        )
    if clinical_stage_number >= 5:
        return (
            "Approved or phase-4 drug in non-IBD indication",
            "Clinical-stage pharmacology exists, but the listed indication is not IBD.",
        )
    if clinical_stage_number >= 3:
        return (
            "Clinical-stage drug in non-IBD indication",
            "Clinical-stage pharmacology exists, but the listed indication is not IBD.",
        )
    if tractable and validated:
        return (
            "Tractable target without IBD clinical candidate",
            "Tractability and independent IBD evidence support preclinical prioritization only.",
        )
    if tractable:
        return (
            "Tractable target without clinical candidate",
            "Tractability supports assay or compound development, not treatment use.",
        )
    if validated:
        return (
            "Genetics-led target with limited tractability",
            "Independent IBD evidence exists, but no clinical-stage pharmacology was identified.",
        )
    return (
        "Limited translational annotation",
        "Neither tractability nor independent IBD evidence was sufficient for treatment positioning.",
    )


def main() -> None:
    configure_proxy()
    table_dir = ROOT / "results" / "tables"
    lead = (
        pd.read_csv(table_dir / "ranked_variant_gene_candidates_final.csv")
        .sort_values("priority_rank")
        .drop_duplicates("variant_id")
        .head(10)
    )
    eqtl = pd.read_csv(table_dir / "external_eqtl_top10_triangulation.csv")
    eqtl_hits = set(
        eqtl.loc[
            eqtl["query_group"].eq("top10_matched")
            & eqtl["any_nominal"].astype(bool),
            "target",
        ]
    )
    coloc = pd.read_csv(table_dir / "top10_eqtl_coloc_results.csv")
    coloc_hits = set(
        coloc.loc[
            coloc["coloc_supported"].fillna(False).astype(bool),
            "target",
        ]
    )

    session = requests.Session()
    rows = []
    errors = []
    for row in lead.itertuples():
        target_symbol = str(row.l2g_target_symbol)
        target_id = str(row.l2g_target_id)
        try:
            target = query_target(session, target_id)
        except Exception as exc:
            errors.append({"target": target_symbol, "error": str(exc)})
            continue

        tractability = [
            item
            for item in target.get("tractability") or []
            if item.get("value")
        ]
        modalities = sorted(
            {
                str(item.get("modality") or item.get("label") or "")
                for item in tractability
                if item.get("modality") or item.get("label")
            }
        )
        drugs = (
            (target.get("drugAndClinicalCandidates") or {}).get("rows") or []
        )
        drug_names = sorted(
            {
                drug["drug"]["name"]
                for drug in drugs
                if drug.get("drug") and drug["drug"].get("name")
            }
        )
        stages = [
            str(drug.get("maxClinicalStage") or "")
            for drug in drugs
        ]
        max_stage = max(stages, key=stage_number) if stages else ""
        ibd_indication_drugs = set()
        non_ibd_indication_drugs = set()
        ibd_indication_stages = []
        target_specific_actions = set()
        minimum_mechanism_target_count = None
        for candidate in drugs:
            drug = candidate.get("drug") or {}
            drug_name = str(drug.get("name") or "")
            indications = (drug.get("indications") or {}).get("rows") or []
            for indication in indications:
                disease_name = str(
                    (indication.get("disease") or {}).get("name") or ""
                ).lower()
                stage = str(indication.get("maxClinicalStage") or "")
                if disease_name in IBD_NAMES:
                    if drug_name:
                        ibd_indication_drugs.add(drug_name)
                    ibd_indication_stages.append(stage)
                elif drug_name:
                    non_ibd_indication_drugs.add(drug_name)
            mechanisms = (
                (drug.get("mechanismsOfAction") or {}).get("rows") or []
            )
            mechanism_targets = set()
            for mechanism in mechanisms:
                current_targets = {
                    str(
                        mechanism_target.get("id")
                        or mechanism_target.get("approvedSymbol")
                        or ""
                    )
                    for mechanism_target in mechanism.get("targets") or []
                }
                current_targets = {
                    value for value in current_targets if value
                }
                mechanism_targets.update(current_targets)
                target_name = str(mechanism.get("targetName") or "")
                if target_name == target_symbol or target_id in current_targets:
                    action = str(mechanism.get("actionType") or "")
                    if action:
                        target_specific_actions.add(action)
            if mechanism_targets:
                count = len({value for value in mechanism_targets if value})
                if count:
                    minimum_mechanism_target_count = (
                        count
                        if minimum_mechanism_target_count is None
                        else min(minimum_mechanism_target_count, count)
                    )
        max_ibd_stage = (
            max(ibd_indication_stages, key=stage_number)
            if ibd_indication_stages
            else ""
        )
        ibd_indication = bool(ibd_indication_drugs)
        pathways = [
            str(pathway.get("pathway", ""))
            for pathway in target.get("pathways") or []
        ]
        relevant_pathways = [
            pathway for pathway in pathways if PATHWAY_TERMS.search(pathway)
        ]
        associated = (
            (target.get("associatedDiseases") or {}).get("rows") or []
        )
        ibd_scores = [
            f"{item['disease']['name']}:{float(item['score']):.3g}"
            for item in associated
            if str(item.get("disease", {}).get("name", "")).lower()
            in IBD_NAMES
        ]
        top_diseases = [
            f"{item['disease']['name']}:{float(item['score']):.3g}"
            for item in associated[:3]
        ]
        independent_validation = bool(
            float(row.external_cluster_count or 0) > 0
            or pd.notna(row.experimental_support)
            or target_symbol in eqtl_hits
        )
        priority, rationale = clinical_translation_status(
            ibd_indication,
            stage_number(max_stage),
            bool(modalities),
            independent_validation,
        )
        rows.append(
            {
                "priority_rank": int(row.priority_rank),
                "target": target_symbol,
                "target_id": target_id,
                "variant": row.canonical_variant,
                "external_cluster_count": int(
                    float(row.external_cluster_count or 0)
                ),
                "experimental_support": row.experimental_support,
                "exact_variant_eqtl": target_symbol in eqtl_hits,
                "formal_colocalization": target_symbol in coloc_hits,
                "tractable_modalities": "; ".join(modalities),
                "clinical_candidate_count": int(
                    (target.get("drugAndClinicalCandidates") or {}).get(
                        "count", 0
                    )
                ),
                "clinical_stage": max_stage,
                "drug_names": "; ".join(drug_names[:8]),
                "ibd_indication_level_candidate": ibd_indication,
                "ibd_indication_clinical_stage": max_ibd_stage,
                "ibd_indication_drugs": "; ".join(
                    sorted(ibd_indication_drugs)
                ),
                "non_ibd_indication_drug_count": len(
                    non_ibd_indication_drugs
                ),
                "target_specific_clinical_actions": "; ".join(
                    sorted(target_specific_actions)
                ),
                "minimum_mechanism_target_count": minimum_mechanism_target_count,
                "clinical_treatment_recommendation_supported": False,
                "clinical_recommendation_boundary": (
                    "No IBD-indication clinical-stage drug and directionality "
                    "evidence are jointly available."
                ),
                "ibd_associated_diseases": "; ".join(ibd_scores),
                "top_associated_diseases": "; ".join(top_diseases),
                "ibd_relevant_pathways": " | ".join(relevant_pathways[:12]),
                "translational_priority": priority,
                "translation_rationale": rationale,
            }
        )

    frame = pd.DataFrame(rows).sort_values("priority_rank")
    frame.to_csv(
        table_dir / "top10_clinical_translation_annotation.csv",
        index=False,
    )
    summary = {
        "targets_annotated": len(frame),
        "targets_with_tractability": int(
            frame["tractable_modalities"].ne("").sum()
        ),
        "targets_with_clinical_candidates": int(
            frame["clinical_candidate_count"].gt(0).sum()
        ),
        "targets_with_ibd_indication_candidates": int(
            frame["ibd_indication_level_candidate"].sum()
        ),
        "ibd_indication_targets": sorted(
            frame.loc[
                frame["ibd_indication_level_candidate"].astype(bool),
                "target",
            ].tolist()
        ),
        "targets_with_non_ibd_clinical_candidates": int(
            frame["clinical_candidate_count"].gt(0).sum()
            - frame["ibd_indication_level_candidate"].sum()
        ),
        "targets_with_target_specific_clinical_actions": int(
            frame["target_specific_clinical_actions"].ne("").sum()
        ),
        "targets_with_eqtl": int(frame["exact_variant_eqtl"].sum()),
        "targets_with_formal_coloc": int(
            frame["formal_colocalization"].sum()
        ),
        "translational_status_counts": frame[
            "translational_priority"
        ].value_counts().to_dict(),
        "clinical_treatment_recommendations_supported": int(
            frame["clinical_treatment_recommendation_supported"].sum()
        ),
        "errors": errors,
        "note": (
            "Clinical-stage drugs were checked for IBD, Crohn and ulcerative "
            "colitis indications and target-specific mechanisms. No treatment "
            "recommendation is issued because an IBD indication and established "
            "therapeutic direction are not jointly available."
        ),
    }
    (table_dir / "top10_clinical_translation_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(frame.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
