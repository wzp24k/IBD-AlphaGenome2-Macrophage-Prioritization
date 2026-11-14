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


def translation_class(
    clinical_stage_number: int,
    tractable: bool,
    validated: bool,
) -> tuple[str, str]:
    if clinical_stage_number >= 3 and validated:
        return "High", "Clinical-stage target with independent IBD evidence"
    if clinical_stage_number >= 3:
        return "Medium", "Clinical-stage target without direct evidence in this study"
    if tractable and validated:
        return "Medium", "Tractable target with independent IBD evidence"
    if tractable:
        return "Exploratory", "Tractable target with early translational annotation"
    if validated:
        return "Exploratory", "Genetics-led target with limited tractability annotation"
    return "Low", "Limited tractability and independent validation"


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
        priority, rationale = translation_class(
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
        "targets_with_eqtl": int(frame["exact_variant_eqtl"].sum()),
        "targets_with_formal_coloc": int(
            frame["formal_colocalization"].sum()
        ),
        "translational_priority_counts": frame[
            "translational_priority"
        ].value_counts().to_dict(),
        "errors": errors,
        "note": (
            "Translational priority is a descriptive triage category based on "
            "tractability, clinical-stage annotations and available independent "
            "evidence. It is not a diagnostic, prognostic or treatment-selection "
            "score, and clinical-stage drugs are not necessarily IBD-specific."
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
