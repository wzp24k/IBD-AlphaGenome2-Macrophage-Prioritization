from __future__ import annotations

import hashlib
import json
import re
import sys
import zipfile
from datetime import datetime
from pathlib import Path

from docx import Document
from openpyxl import load_workbook

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1119manuscript.docx"
WORKBOOK = ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx"
FIGURE = ROOT / "results" / "figures" / "figure13_clinical_translation.png"
OUTPUT = ROOT / "results" / "clinical_translation_extension_audit.json"
MARKDOWN = ROOT / "results" / "CLINICAL_TRANSLATION_EXTENSION_AUDIT.md"


def write_markdown_report() -> None:
    support = json.loads(
        (ROOT / "results" / "tables" / "cd_uc_cross_ancestry_summary.json").read_text(
            encoding="utf-8"
        )
    )
    pqtl = json.loads(
        (ROOT / "results" / "tables" / "pqtl_translation_summary.json").read_text(
            encoding="utf-8"
        )
    )
    sun_pqtl = json.loads(
        (
            ROOT
            / "results"
            / "tables"
            / "prkcb_sun2018_pqtl_summary.json"
        ).read_text(encoding="utf-8")
    )
    treatment = json.loads(
        (ROOT / "results" / "tables" / "treatment_response_summary.json").read_text(
            encoding="utf-8"
        )
    )
    readiness = json.loads(
        (
            ROOT
            / "results"
            / "tables"
            / "clinical_development_readiness_summary.json"
        ).read_text(encoding="utf-8")
    )
    robustness = json.loads(
        (
            ROOT
            / "results"
            / "tables"
            / "clinical_translation_robustness_summary.json"
        ).read_text(encoding="utf-8")
    )
    incremental = json.loads(
        (
            ROOT
            / "results"
            / "tables"
            / "benchmark_incremental_model_summary.json"
        ).read_text(encoding="utf-8")
    )
    incremental_any = incremental["prespecified_by_task"][
        "External genetic support, any cohort"
    ]
    lines = [
        "# Clinical translation extension",
        "",
        f"Generated: {datetime.now().astimezone():%Y-%m-%d}.",
        "",
        "## 1. CD/UC and cross-ancestry validation",
        "",
        f"- Crohn targets with L2G >= 0.5: {support['cd_supported_targets_l2g_ge_0_5']}.",
        f"- UC targets with L2G >= 0.5: {support['uc_supported_targets_l2g_ge_0_5']}.",
        f"- Targets supported in both diseases: {support['both_disease_targets']}.",
        "- PRKCB cis-MR was positive in all four outcome cohorts.",
        "",
        "## 2. pQTL, colocalization and PheW-MR",
        "",
        f"- Sun 2018 plasma pQTL variants: {sun_pqtl['pqtl_variants']:,}.",
        f"- Independent pQTL instruments: {sun_pqtl['independent_instruments']}.",
        f"- Positive pQTL-MR outcomes: {sun_pqtl['positive_direction']}/"
        f"{sun_pqtl['mr_tests']}.",
        f"- Colocalization H4 >= 0.8: {sun_pqtl['coloc_h4_ge_0_8']}/"
        f"{sun_pqtl['coloc_tests']} outcomes.",
        f"- Maximum H4: {sun_pqtl['max_h4']:.3f}.",
        "- Legacy SomaScan and curated pQTL results remain sensitivity analyses.",
        f"- PheW-MR queried {pqtl['phewas_phenotypes_queried']} traits; "
        f"{pqtl['phewas_phenotypes_evaluable']} were evaluable.",
        "",
        "## 3. Public treatment-response programs",
        "",
        f"- Datasets: {', '.join(treatment['datasets'])}.",
        "- The PRKCB network program was lower in responders across both cohorts.",
        "- The macrophage marker program showed the same direction.",
        "- Direct PRKCB expression had a wider interval that included zero.",
        "",
        "## 4. Clinical development readiness",
        "",
        f"- PRKCB: {readiness['top_score']}/10 and the highest-ranked target.",
        "- IBD-indication clinical stage, not maximum stage across all indications, contributes to scoring.",
        f"- PRKCB ranked {readiness['prkcb_rank_min']} across four prespecified "
        "weighting schemes.",
        f"- PRKCB ranked first in {robustness['prkcb_probability_rank_1']:.1%} and top two in "
        f"{robustness['prkcb_probability_top_2']:.1%} of 100,000 randomized weight draws.",
        "",
        "## 5. Held-out incremental benchmark",
        "",
        f"- The prespecified PIP+context model improved AUROC by "
        f"{incremental_any['delta_auroc']:.3f} for external support in any cohort "
        f"(95% CI {incremental_any['delta_auroc_ci_low']:.3f} to "
        f"{incremental_any['delta_auroc_ci_high']:.3f}).",
        "- The advantage was not stable for support in both cohorts or "
        "replicated exact-variant eQTL.",
        "",
        "## Parallel robustness checks",
        "",
        f"- The robustness run used {robustness['cpu_count']} logical CPUs and "
        f"{robustness['permutation_chunks']} permutation chunks.",
        f"- Treatment-response permutations: {robustness['treatment_permutations_per_feature']} per feature.",
        f"- pQTL sensitivity retained a positive direction in "
        f"{robustness['pqtl_positive_sensitivity_fraction']:.0%} of tests.",
        f"- Maximum pQTL-colocalization H4 across prior combinations: "
        f"{robustness['pqtl_coloc_sensitivity_max_h4']:.3f}.",
        "",
        "## Interpretation boundary",
        "",
        "The extension is an evidence-integration analysis. The primary Sun 2018 "
        "pQTL is a one-instrument protein layer. Legacy pQTL and PheW-MR are "
        "exploratory. Treatment-response results are observational program-level evidence.",
    ]
    MARKDOWN.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    document = Document(DOCX)
    paragraphs = [paragraph.text.strip() for paragraph in document.paragraphs]
    nonempty = [text for text in paragraphs if text]
    discussion_index = nonempty.index("4. Discussion")
    body = "\n".join(nonempty[:discussion_index])
    table_mentions = sorted(set(map(int, re.findall(r"\bTable\s+(\d+)", body))))
    figure_mentions = sorted(set(map(int, re.findall(r"\bFigure\s+(\d+)", body))))
    table_captions = sorted(
        {
            int(match.group(1))
            for text in nonempty
            if (match := re.match(r"^Table\s+(\d+)\.", text))
        }
    )
    figure_captions = sorted(
        {
            int(match.group(1))
            for text in nonempty
            if (match := re.match(r"^Figure\s+(\d+)\.", text))
        }
    )
    with zipfile.ZipFile(DOCX) as archive:
        media_hashes = {
            hashlib.sha256(archive.read(name)).hexdigest()
            for name in archive.namelist()
            if name.startswith("word/media/")
        }
    figure_hash = hashlib.sha256(FIGURE.read_bytes()).hexdigest()
    workbook = load_workbook(WORKBOOK, read_only=True)
    narrative_audit = json.loads(
        (ROOT / "results" / "translational_narrative_audit.json").read_text(
            encoding="utf-8"
        )
    )
    readiness_rows = {}
    with (
        ROOT
        / "results"
        / "tables"
        / "clinical_development_readiness_matrix.csv"
    ).open(encoding="utf-8") as handle:
        import csv

        for row in csv.DictReader(handle):
            readiness_rows[row["target"]] = row
    expected_sheets = {
        "CD_UC_target_support",
        "PRKCB_cross_ancestry_MR",
        "pQTL_cis_MR",
        "pQTL_colocalization",
        "PheW_MR",
        "Treatment_response_programs",
        "Treatment_response_meta",
        "Clinical_readiness",
        "Clinical_readiness_rank_sensiti",
        "Benchmark_incremental_metrics",
        "Benchmark_incremental_deltas",
        "Sun2018_pQTL_cis",
        "Sun2018_pQTL_MR",
        "Sun2018_pQTL_coloc",
    }
    checks = {
        "twelve_tables": len(document.tables) == 12,
        "thirteen_figures": len(document.inline_shapes) == 13,
        "table_mentions_1_to_12": table_mentions == list(range(1, 13)),
        "figure_mentions_1_to_13": figure_mentions == list(range(1, 14)),
        "table_captions_1_to_12": table_captions == list(range(1, 13)),
        "figure_captions_1_to_13": figure_captions == list(range(1, 14)),
        "figure13_embedded": figure_hash in media_hashes,
        "methods_extension_present": "2.10 Clinical translation extension" in nonempty,
        "results_extension_present": "3.15 Clinical translation extension" in nonempty,
        "table12_positive": "Clinical development readiness matrix" in body,
        "workbook_extensions_present": expected_sheets.issubset(set(workbook.sheetnames)),
        "prkcb_public_perturbation_not_counted": float(
            readiness_rows["PRKCB"]["public_perturbation_score"]
        )
        == 0.0,
        "prkcb_uses_ibd_specific_stage": readiness_rows["PRKCB"][
            "ibd_indication_clinical_stage"
        ]
        == "PHASE_2",
        "pqtl_strong_support_in_score": readiness_rows["PRKCB"][
            "pqtl_evidence_class"
        ]
        == "Strong"
        and float(readiness_rows["PRKCB"]["pqtl_mr_score"]) == 1.0,
        "incremental_benchmark_present": "Repeated grouped out-of-fold models"
        in body,
        "public_causal_coverage_recheck_present": "Public MPRA and macrophage "
        "CRISPRi resources currently lack" in body,
        "translational_narrative_audit_pass": narrative_audit[
            "all_checks_pass"
        ],
    }
    report = {
        "document": str(DOCX),
        "tables": len(document.tables),
        "figures": len(document.inline_shapes),
        "table_mentions": table_mentions,
        "figure_mentions": figure_mentions,
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    OUTPUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown_report()
    print(json.dumps(report, indent=2))
    if not report["all_checks_pass"]:
        raise SystemExit("Clinical translation extension audit failed.")


if __name__ == "__main__":
    main()
