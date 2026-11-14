from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import requests
from openpyxl import load_workbook

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"
PRKCB_ID = "ENSG00000166501"
MAP3K8_ID = "ENSG00000107968"
PKC_GENES = {
    "PRKCA",
    "PRKCB",
    "PRKCD",
    "PRKCE",
    "PRKCG",
    "PRKCH",
    "PRKCI",
    "PRKCQ",
    "PRKCZ",
    "PRKD1",
    "PRKD2",
    "PRKD3",
}
TARGET_QUERY = """
query TargetDrugs($id: String!) {
  target(ensemblId: $id) {
    drugAndClinicalCandidates {
      rows {
        maxClinicalStage
        drug {
          name
          drugType
          maximumClinicalStage
          crossReferences {
            source
            ids
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
  }
}
"""


def query_target_drugs(session: requests.Session, target_id: str) -> dict:
    response = session.post(
        OPENTARGETS_URL,
        json={"query": TARGET_QUERY, "variables": {"id": target_id}},
        timeout=60,
    )
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError(
            "; ".join(error["message"] for error in payload["errors"])
        )
    return payload["data"]["target"]


def drug_selectivity_table(
    session: requests.Session,
    target_id: str,
) -> pd.DataFrame:
    target = query_target_drugs(session, target_id)
    rows = []
    for candidate in (
        target.get("drugAndClinicalCandidates") or {}
    ).get("rows") or []:
        drug = candidate.get("drug") or {}
        for mechanism in (
            drug.get("mechanismsOfAction") or {}
        ).get("rows") or []:
            targets = mechanism.get("targets") or []
            target_symbols = sorted(
                {
                    target_info.get("approvedSymbol")
                    for target_info in targets
                    if target_info.get("approvedSymbol")
                }
            )
            if not any(
                target_info.get("id") == target_id
                for target_info in targets
            ):
                continue
            pkc_targets = sorted(set(target_symbols) & PKC_GENES)
            cross_refs = drug.get("crossReferences") or []
            drugbank = next(
                (
                    reference.get("ids") or []
                    for reference in cross_refs
                    if reference.get("source") == "drugbank"
                ),
                [],
            )
            rows.append(
                {
                    "drug": drug.get("name"),
                    "drug_type": drug.get("drugType"),
                    "target_stage": candidate.get("maxClinicalStage"),
                    "mechanism": mechanism.get("mechanismOfAction"),
                    "action_type": mechanism.get("actionType"),
                    "target_symbols": "; ".join(target_symbols),
                    "pkc_gene_count": len(pkc_targets),
                    "pkc_targets": "; ".join(pkc_targets),
                    "prkcb_specific": pkc_targets == ["PRKCB"],
                    "drugbank_ids": "; ".join(drugbank),
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def cedar2_prkcb_table() -> tuple[pd.DataFrame, dict]:
    workbook_path = (
        ROOT
        / "data"
        / "external"
        / "clinical"
        / "CEDAR2_supplementary_data_1_29.xlsx"
    )
    if not workbook_path.exists():
        return pd.DataFrame(), {"status": "supplementary_workbook_not_found"}
    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    annotation = {}
    for row in workbook["ST11"].iter_rows(min_row=2, values_only=True):
        if row[3]:
            annotation[str(row[3])] = {
                "cell_type": row[1],
                "node_or_leaf": row[3],
                "nr_cells": row[4],
                "nr_eqtl": row[5],
                "abbreviation": row[6],
            }
    gut_rows = []
    for row in workbook["ST12"].iter_rows(min_row=3, values_only=True):
        if row[1] != PRKCB_ID:
            continue
        metadata = annotation.get(str(row[0]), {})
        gut_rows.append(
            {
                "dataset": "CEDAR2_gut",
                "node_or_leaf": row[0],
                "annotated_cell_type": metadata.get("cell_type"),
                "cell_abbreviation": metadata.get("abbreviation"),
                "gene": row[2],
                "gene_id": row[1],
                "top_variant": row[7],
                "beta": row[9],
                "nominal_pvalue": row[10],
                "window_adjusted_pvalue": row[11],
                "qvalue": row[12],
            }
        )
    blood_rows = []
    for row in workbook["ST5"].iter_rows(min_row=3, values_only=True):
        if row[1] != PRKCB_ID:
            continue
        blood_rows.append(
            {
                "dataset": "CEDAR2_blood",
                "node_or_leaf": row[0],
                "annotated_cell_type": row[0],
                "cell_abbreviation": row[0],
                "gene": row[2],
                "gene_id": row[1],
                "top_variant": row[7],
                "beta": row[9],
                "nominal_pvalue": row[10],
                "window_adjusted_pvalue": row[11],
                "qvalue": row[12],
            }
        )
    result = pd.DataFrame(gut_rows + blood_rows)
    if not result.empty:
        result["macrophage_like"] = result[
            "annotated_cell_type"
        ].fillna("").str.contains("macrophage", case=False)
        result["exact_rs7404003"] = result["top_variant"].eq(
            "chr16:23853326:C:A"
        )
    macrophage_nodes = sorted(
        {
            str(item["node_or_leaf"])
            for item in annotation.values()
            if str(item["cell_type"]).lower().find("macrophage") >= 0
        }
    )
    summary = {
        "status": "available",
        "macrophage_like_nodes": macrophage_nodes,
        "prkcb_significant_gut_rows": int(
            result["dataset"].eq("CEDAR2_gut").sum()
        )
        if not result.empty
        else 0,
        "prkcb_significant_macrophage_rows": int(
            (
                result["dataset"].eq("CEDAR2_gut")
                & result["macrophage_like"].astype(bool)
            ).sum()
        )
        if not result.empty
        else 0,
        "prkcb_exact_variant_blood_rows": int(
            (
                result["dataset"].eq("CEDAR2_blood")
                & result["exact_rs7404003"].astype(bool)
            ).sum()
        )
        if not result.empty
        else 0,
        "note": (
            "CEDAR2 significant-eQTL tables were queried. The full 1.5 GB "
            "macrophage nominal file was not exhaustively parsed because the "
            "significant macrophage node contained no PRKCB eQTL row."
        ),
    }
    return result, summary


def clinical_trial_tables() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    session = requests.Session()
    sot = session.get(
        CTGOV_URL,
        params={
            "query.intr": "sotrastaurin",
            "query.cond": "Inflammatory Bowel Diseases",
            "pageSize": 50,
            "format": "json",
        },
        timeout=60,
    )
    enz = session.get(
        CTGOV_URL,
        params={
            "query.intr": "enzastaurin",
            "query.cond": "Inflammatory Bowel Diseases",
            "pageSize": 50,
            "format": "json",
        },
        timeout=60,
    )
    sot.raise_for_status()
    enz.raise_for_status()
    sot_studies = sot.json().get("studies", [])
    enz_studies = enz.json().get("studies", [])
    trial_rows = []
    for study in sot_studies:
        protocol = study.get("protocolSection") or {}
        identification = protocol.get("identificationModule") or {}
        status = protocol.get("statusModule") or {}
        design = protocol.get("designModule") or {}
        conditions = protocol.get("conditionsModule") or {}
        arms = protocol.get("armsInterventionsModule") or {}
        trial_rows.append(
            {
                "drug": "sotrastaurin",
                "nct_id": identification.get("nctId"),
                "title": identification.get("briefTitle"),
                "status": status.get("overallStatus"),
                "phase": "; ".join(design.get("phases") or []),
                "conditions": "; ".join(conditions.get("conditions") or []),
                "interventions": "; ".join(
                    intervention.get("name", "")
                    for intervention in arms.get("interventions") or []
                ),
            }
        )
    for study in enz_studies:
        protocol = study.get("protocolSection") or {}
        identification = protocol.get("identificationModule") or {}
        status = protocol.get("statusModule") or {}
        design = protocol.get("designModule") or {}
        conditions = protocol.get("conditionsModule") or {}
        trial_rows.append(
            {
                "drug": "enzastaurin",
                "nct_id": identification.get("nctId"),
                "title": identification.get("briefTitle"),
                "status": status.get("overallStatus"),
                "phase": "; ".join(design.get("phases") or []),
                "conditions": "; ".join(conditions.get("conditions") or []),
                "interventions": "",
            }
        )
    results = pd.DataFrame(
        [
            {
                "drug": "sotrastaurin",
                "nct_id": "NCT00572585",
                "indication": "active moderate-to-severe ulcerative colitis",
                "phase": "Phase 2b",
                "dose": "300 mg bid for 2 weeks, then 200 mg bid for 2 weeks",
                "randomized": "42 sotrastaurin; 18 placebo",
                "clinical_remission": "6/36 (16.7%) vs 0/17 (0%)",
                "posterior_difference": "12.2% (95% CrI -2.1 to 27.4)",
                "partial_mayo_change_day29": "-3.7 vs -1.9",
                "any_adverse_event": "31/42 (73.8%) vs 10/18 (55.6%)",
                "serious_adverse_event": "3/42 (7.1%) vs 2/18 (11.1%)",
                "status": "completed but terminated early for strategic reasons",
                "publication": "No peer-reviewed primary publication listed by trial report",
                "source": "Novartis CAEB071A2210 trial report and ClinicalTrials.gov",
            },
            {
                "drug": "enzastaurin",
                "nct_id": "",
                "indication": "No IBD trial identified",
                "phase": "",
                "dose": "",
                "randomized": "",
                "clinical_remission": "",
                "posterior_difference": "",
                "partial_mayo_change_day29": "",
                "any_adverse_event": "",
                "serious_adverse_event": "",
                "status": "not tested in IBD in ClinicalTrials.gov search",
                "publication": "",
                "source": "ClinicalTrials.gov API",
            },
        ]
    )
    summary = {
        "sotrastaurin_ibd_trials": len(sot_studies),
        "enzastaurin_ibd_trials": len(enz_studies),
        "sotrastaurin_efficacy_signal": (
            "Posterior probability of any clinical-remission improvement was "
            "95.5%, but the 95% credibility interval crossed zero and the "
            "probability of a >=20% difference was only 14.2%."
        ),
        "sotrastaurin_safety": (
            "Higher any-AE and gastrointestinal/renal AE frequencies than "
            "placebo; serious AE frequency was not higher in the small trial."
        ),
        "treatment_recommendation": (
            "Not supported: small early-terminated phase 2b study, no primary "
            "publication, no confirmatory trial and no target-selective drug."
        ),
    }
    return pd.DataFrame(trial_rows), results, summary


def risk_direction_table() -> pd.DataFrame:
    table_dir = ROOT / "results" / "tables"
    lead = pd.read_csv(
        table_dir / "ranked_variant_gene_candidates_final.csv"
    ).sort_values("priority_rank")
    lookups = pd.read_csv(
        table_dir / "external_eqtl_all_variant_lookups.csv"
    )
    lookups = lookups[lookups["query_group"].eq("top10_matched")]
    selected = lead[lead["l2g_target_symbol"].isin(["PRKCB", "MAP3K8"])]
    rows = []
    for row in selected.itertuples():
        exact = lookups[
            lookups["variant_id"].eq(row.variant_id)
            & lookups["target"].eq(row.l2g_target_symbol)
        ].sort_values("pvalue")
        if exact.empty:
            continue
        eqtl_beta = float(exact.iloc[0]["beta"])
        risk_effect = (
            "risk_increases_target"
            if float(row.study_beta) * eqtl_beta > 0
            else "risk_decreases_target"
        )
        rows.append(
            {
                "target": row.l2g_target_symbol,
                "variant": row.canonical_variant,
                "study_beta_alt": row.study_beta,
                "eqtl_beta_alt": eqtl_beta,
                "risk_effect_on_target": risk_effect,
                "expected_therapeutic_action": (
                    "inhibition"
                    if risk_effect == "risk_increases_target"
                    else "activation"
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    session = requests.Session()
    selectivity = pd.concat(
        [
            drug_selectivity_table(session, PRKCB_ID),
            drug_selectivity_table(session, MAP3K8_ID),
        ],
        ignore_index=True,
    )
    selectivity.to_csv(
        table_dir / "prkcb_map3k8_drug_selectivity.csv",
        index=False,
    )
    cedar2, cedar2_summary = cedar2_prkcb_table()
    cedar2.to_csv(
        table_dir / "prkcb_cedar2_independent_eqtl.csv",
        index=False,
    )
    trial_search, trial_results, trial_summary = clinical_trial_tables()
    trial_search.to_csv(
        table_dir / "prkcb_ibd_trial_search.csv",
        index=False,
    )
    trial_results.to_csv(
        table_dir / "prkcb_ibd_trial_results.csv",
        index=False,
    )
    risk_direction = risk_direction_table()
    risk_direction.to_csv(
        table_dir / "prkcb_map3k8_risk_direction.csv",
        index=False,
    )
    summary = {
        "drug_selectivity_rows": len(selectivity),
        "prkcb_specific_drugs": sorted(
            selectivity.loc[
                selectivity["prkcb_specific"].astype(bool),
                "drug",
            ].unique()
        ),
        "pan_pkc_drugs": sorted(
            selectivity.loc[
                selectivity["pkc_gene_count"].gt(1),
                "drug",
            ].unique()
        ),
        "cedar2": cedar2_summary,
        "clinical_trials": trial_summary,
        "sotrastaurin_clinical_remission_difference_percent": 12.2,
        "sotrastaurin_remission_ci_low_percent": -2.1,
        "sotrastaurin_remission_ci_high_percent": 27.4,
        "sotrastaurin_probability_any_improvement_percent": 95.5,
        "sotrastaurin_probability_ge20_percent": 14.2,
        "clinical_treatment_recommendation_supported": False,
        "risk_direction": risk_direction.to_dict("records"),
        "mpra_design_status": (
            "rs7404003 is already included in the prespecified top-three "
            "lentiMPRA package with three 251-bp windows, both alleles, six "
            "donors and resting/LPS/IFN-gamma/IL-23 conditions. No result exists."
        ),
        "crispri_status": (
            "Gene-level PRKCB CRISPRi remains complementary, not allele-specific. "
            "Four independent guides plus safe-harbor and non-targeting controls "
            "should measure PRKCB, TNF, IL23A, CXCL10 and NF-kB readouts."
        ),
        "treatment_recommendation": (
            "No clinical treatment recommendation is supported. Sotrastaurin "
            "has a small early-terminated phase 2b UC study with an uncertain "
            "efficacy interval and no target-selective pharmacology."
        ),
    }
    (table_dir / "prkcb_map3k8_followup_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(selectivity.to_string(index=False))
    print(cedar2.to_string(index=False))
    print(trial_search.to_string(index=False))
    print(trial_results.to_string(index=False))
    print(risk_direction.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
