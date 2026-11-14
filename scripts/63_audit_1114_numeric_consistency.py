from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


TABLE_DIR = ROOT / "results" / "tables"
MANUSCRIPT = ROOT / "manuscript" / "manuscript.md"
OUTPUT = ROOT / "results" / "1114_numeric_consistency_audit.json"


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    summary = json.loads(
        (TABLE_DIR / "analysis_summary.json").read_text(encoding="utf-8")
    )
    active = json.loads(
        (TABLE_DIR / "active_ibd_validation_summary.json").read_text(
            encoding="utf-8"
        )
    )
    external = json.loads(
        (
            TABLE_DIR / "external_gwas_block_bootstrap_summary.json"
        ).read_text(encoding="utf-8")
    )
    mechanism = json.loads(
        (TABLE_DIR / "prkcb_causal_mechanism_summary.json").read_text(
            encoding="utf-8"
        )
    )
    susie = json.loads(
        (TABLE_DIR / "prkcb_susie_coloc_summary.json").read_text(
            encoding="utf-8"
        )
    )
    trial = json.loads(
        (TABLE_DIR / "prkcb_map3k8_followup_summary.json").read_text(
            encoding="utf-8"
        )
    )
    variant_sets = pd.read_csv(
        ROOT / "data" / "processed" / "variant_set_summary.csv"
    )
    lead = pd.read_csv(TABLE_DIR / "ranked_lead_variants_final.csv")
    pathways = pd.read_csv(
        TABLE_DIR / "external_gwas_block_bootstrap_pathways.csv"
    )
    perturb = pd.read_csv(
        TABLE_DIR / "perturbseq_gse327644_target_knockdown.csv"
    )
    differential = pd.read_csv(
        TABLE_DIR / "perturbseq_gse327644_differential_summary.csv"
    )
    gtex = pd.read_csv(
        TABLE_DIR / "experimental_gtex_prkcb_eqtl.csv"
    ).iloc[0]
    external_l2g = pd.read_csv(
        TABLE_DIR / "external_ibd_l2g_predictions.csv"
    )
    myeloid = pd.read_csv(
        TABLE_DIR / "top10_eqtl_coloc_variants.csv"
    )
    myeloid = myeloid[
        myeloid["target"].eq("PRKCB")
        & myeloid["rsid"].eq("rs7404003")
    ].iloc[0]
    claims = pd.read_csv(TABLE_DIR / "claim_calibration.csv")

    pdgfb_self = perturb[
        perturb["target"].eq("PDGFB")
        & perturb["condition"].eq("Combined")
    ].iloc[0]
    prkcq_self = perturb[
        perturb["target"].eq("PRKCQ")
        & perturb["condition"].eq("Combined")
    ].iloc[0]
    il23r_self = perturb[
        perturb["target"].eq("IL23R")
        & perturb["condition"].eq("Combined")
    ].iloc[0]
    pdgfb_de = differential[
        differential["target"].eq("PDGFB")
        & differential["condition"].eq("Combined")
    ].iloc[0]

    checks = {
        "primary_variants_395": summary["primary_variants"] == 395,
        "independent_loci_281": summary["independent_loci"] == 281,
        "atlas_scored_373": summary["variants_with_atlas_target"] == 373
        and "373 of 395" in text,
        "tier_counts_sum_281": (
            summary["tier_1"]
            + summary["tier_2"]
            + int(lead["evidence_tier"].eq("Tier 3").sum())
            + int(lead["evidence_tier"].eq("Tier 4").sum())
            == 281
        ),
        "study_rows_sum_1490": int(
            variant_sets["unique_noncoding_variants"].sum()
        )
        == 1490,
        "external_targets_510_73_19": (
            external["genes"] == 510
            and external["externally_supported_genes"] == 73
            and int(
                external_l2g.groupby("target_symbol")["cluster"]
                .nunique()
                .ge(2)
                .sum()
            )
            == 19
        ),
        "pathway_rows_4": len(pathways) == 4,
        "prkcb_mr_ivw_value": math.isclose(
            mechanism["ivw_beta"],
            0.051739900586056174,
            rel_tol=1e-8,
        )
        and "0.0517" in text,
        "prkcb_smr_and_susie": (
            mechanism["smr_p"] < 0.001
            and susie["primary_max_H4"] >= 0.8
            and "H4 = 1.0" in text
        ),
        "prkcb_direction_discrepancy_reported": bool(
            abs(float(gtex["nes"])) > 0
            and gtex["nes"] < 0
            and myeloid["eql_beta"] > 0
            and "tissue-context discrepancy" in text
            and "colon and myeloid eQTL directions were discordant" in text
        ),
        "pdgfb_perturbseq_values": (
            abs(pdgfb_self["mean_log2fc_vs_ntc"] + 0.636) < 0.001
            and abs(pdgfb_self["fdr"] - 0.00879) < 0.001
            and pdgfb_de["fdr05_genes"] == 460
        ),
        "prkcq_il23r_perturbseq_values": (
            abs(prkcq_self["mean_log2fc_vs_ntc"] + 0.104) < 0.001
            and abs(prkcq_self["fdr"] - 0.826) < 0.002
            and abs(il23r_self["mean_log2fc_vs_ntc"] - 0.107) < 0.001
            and abs(il23r_self["fdr"] - 0.845) < 0.002
        ),
        "clinical_trial_values": (
            trial["sotrastaurin_remission_ci_low_percent"] < 0
            < trial["sotrastaurin_remission_ci_high_percent"]
            and "12.2%" in text
            and "The publicly reported trial cited in this study is NCT00572585."
            in text
        ),
        "claim_calibration_eight_rows": len(claims) == 8,
        "top28_active_support_21_and_17": (
            active["top28_multicohort_supported"] == 21
            and active["top28_strictly_replicated"] == 17
            and "21 were concordant" in text
            and "Seventeen were concordant" in text
        ),
        "top28_selection_bias_wording": (
            "Among the preselected top 28 loci" in text
            and "not driven only by absolute-effect outliers" not in text
        ),
        "cedar2_exact_variant_precision": (
            "supports the locus-level direction" in text
            and "but not the exact variant in a colon macrophage node" in text
        ),
        "model_novelty_wording_removed": (
            "model novelty" not in text
            and text.count("model prediction") >= 2
        ),
    }
    checks = {key: bool(value) for key, value in checks.items()}
    report = {
        "source_manuscript": str(MANUSCRIPT),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "source_values": {
            "primary_variants": summary["primary_variants"],
            "independent_loci": summary["independent_loci"],
            "tier_1": summary["tier_1"],
            "tier_2": summary["tier_2"],
            "tier_3": int(lead["evidence_tier"].eq("Tier 3").sum()),
            "tier_4": int(lead["evidence_tier"].eq("Tier 4").sum()),
            "external_genes": external["genes"],
            "externally_supported_genes": external[
                "externally_supported_genes"
            ],
            "top28_multicohort_supported": active[
                "top28_multicohort_supported"
            ],
            "top28_strictly_replicated": active[
                "top28_strictly_replicated"
            ],
            "gtex_prkcb_nes": float(gtex["nes"]),
            "monocyte_prkcb_beta": float(myeloid["eql_beta"]),
            "claim_rows": len(claims),
        },
    }
    OUTPUT.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))
    if not report["all_checks_pass"]:
        raise SystemExit("Numeric consistency audit failed.")


if __name__ == "__main__":
    main()
