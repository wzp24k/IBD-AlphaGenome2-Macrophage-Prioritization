from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


SOURCE = ROOT / "1115manuscript.docx"
OUTPUT = ROOT / "1119manuscript.docx"
MANUSCRIPT_COPY = (
    ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD_clinical_translation.docx"
)
BACKUP = ROOT / "_archive" / "1115manuscript_before_clinical_translation.docx"
WORKBOOK = ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx"
FIGURE = ROOT / "results" / "figures" / "figure13_clinical_translation.png"
TABLE_DIR = ROOT / "results" / "tables"
AUDIT = ROOT / "results" / "clinical_translation_extension_audit.json"


def new_paragraph(document: Document, text: str = ""):
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.paragraph_format.line_spacing = 1.5
    if text:
        run = paragraph.add_run(text)
        run.font.name = "Times New Roman"
        run.font.size = Pt(11)
    return paragraph


def new_heading(document: Document, text: str):
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.keep_with_next = True
    run = paragraph.add_run(text)
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(12)
    return paragraph


def replace_paragraph_text(paragraph, text: str) -> None:
    for child in list(paragraph._p):
        if child.tag != qn("w:pPr"):
            paragraph._p.remove(child)
    run = paragraph.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = Pt(11)


def move_before(anchor, element) -> None:
    if hasattr(element, "_p"):
        child = element._p
    elif hasattr(element, "_tbl"):
        child = element._tbl
    else:
        child = element
    anchor._p.addprevious(child)


def set_three_line_borders(table) -> None:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        element = OxmlElement(f"w:{edge}")
        if edge in {"top", "bottom"}:
            element.set(qn("w:val"), "single")
            element.set(qn("w:sz"), "12")
        else:
            element.set(qn("w:val"), "none")
            element.set(qn("w:sz"), "0")
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), "000000")
        borders.append(element)
    for cell in table.rows[0].cells:
        tc_pr = cell._tc.get_or_add_tcPr()
        tc_borders = tc_pr.first_child_found_in("w:tcBorders")
        if tc_borders is None:
            tc_borders = OxmlElement("w:tcBorders")
            tc_pr.append(tc_borders)
        bottom = OxmlElement("w:bottom")
        bottom.set(qn("w:val"), "single")
        bottom.set(qn("w:sz"), "8")
        bottom.set(qn("w:space"), "0")
        bottom.set(qn("w:color"), "000000")
        tc_borders.append(bottom)


def add_readiness_table(document: Document) -> object:
    frame = pd.read_csv(TABLE_DIR / "clinical_development_readiness_matrix.csv")
    frame = frame.sort_values("readiness_score_0_10", ascending=False)
    headers = [
        "Target",
        "CD/UC support",
        "Exact eQTL",
        "Coloc",
        "Perturbation",
        "pQTL MR",
        "Response",
        "IBD stage",
        "Score",
        "Class",
    ]
    rows = []
    for _, row in frame.iterrows():
        support = (
            "Both"
            if row["disease_specific_score"] >= 1
            else "One"
            if row["disease_specific_score"] > 0
            else "None"
        )
        pqtl_label = (
            "Not assessed"
            if pd.isna(row["pqtl_evidence_class"])
            else str(row["pqtl_evidence_class"])
        )
        rows.append(
            [
                row["target"],
                support,
                "Yes" if row["exact_eqtl_score"] else "No",
                "Yes" if row["colocalization_score"] else "No",
                "Yes" if row["public_perturbation_score"] else "No",
                pqtl_label,
                "Yes" if row["treatment_response_score"] else "No",
                str(row["ibd_indication_clinical_stage"]).replace("_", " ").title()
                if pd.notna(row["ibd_indication_clinical_stage"])
                and str(row["ibd_indication_clinical_stage"]).strip()
                else "None",
                f'{row["readiness_score_0_10"]:.1f}',
                row["readiness_class"],
            ]
        )
    table = document.add_table(rows=len(rows) + 1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    set_three_line_borders(table)
    for row_index, values in enumerate([headers, *rows]):
        for column_index, value in enumerate(values):
            cell = table.rows[row_index].cells[column_index]
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.space_before = Pt(1)
            paragraph.paragraph_format.space_after = Pt(1)
            paragraph.alignment = (
                WD_ALIGN_PARAGRAPH.CENTER
                if row_index == 0
                else WD_ALIGN_PARAGRAPH.LEFT
            )
            run = paragraph.add_run(str(value))
            run.font.name = "Times New Roman"
            run.font.size = Pt(8)
            run.bold = row_index == 0
    return table


def add_figure(document: Document, width: float = 6.4):
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.add_run().add_picture(str(FIGURE), width=Inches(width))
    return paragraph


def append_workbook_sheets() -> None:
    tables = {
        "CD_UC_target_support": pd.read_csv(TABLE_DIR / "top10_cd_uc_target_support.csv"),
        "PRKCB_cross_ancestry_MR": pd.read_csv(TABLE_DIR / "prkcb_cross_ancestry_mr.csv"),
        "pQTL_cis_MR": pd.read_csv(TABLE_DIR / "pqtl_cis_mr_ibd.csv"),
        "pQTL_colocalization": pd.read_csv(TABLE_DIR / "pqtl_ibd_colocalization.csv"),
        "PheW_MR": pd.read_csv(TABLE_DIR / "prkcb_pqtl_phewas_mr.csv"),
        "Treatment_response_programs": pd.read_csv(
            TABLE_DIR / "treatment_response_program_associations.csv"
        ),
        "Treatment_response_meta": pd.read_csv(
            TABLE_DIR / "treatment_response_meta_analysis.csv"
        ),
        "Clinical_readiness": pd.read_csv(
            TABLE_DIR / "clinical_development_readiness_matrix.csv"
        ),
        "Clinical_readiness_rank_sensitivity": pd.read_csv(
            TABLE_DIR / "clinical_readiness_rank_sensitivity.csv"
        ),
        "Benchmark_incremental_metrics": pd.read_csv(
            TABLE_DIR / "benchmark_incremental_model_metrics.csv"
        ),
        "Benchmark_incremental_deltas": pd.read_csv(
            TABLE_DIR / "benchmark_incremental_model_deltas.csv"
        ),
        "Sun2018_pQTL_cis": pd.read_csv(
            TABLE_DIR / "prkcb_sun2018_pqtl_cis_grch38.csv"
        ),
        "Sun2018_pQTL_MR": pd.read_csv(
            TABLE_DIR / "prkcb_sun2018_pqtl_mr_ibd.csv"
        ),
        "Sun2018_pQTL_coloc": pd.read_csv(
            TABLE_DIR / "prkcb_sun2018_pqtl_coloc.csv"
        ),
    }
    workbook = load_workbook(WORKBOOK)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    for sheet_name, frame in tables.items():
        safe_name = sheet_name[:31]
        for existing in list(workbook.sheetnames):
            if existing == safe_name:
                del workbook[existing]
        worksheet = workbook.create_sheet(safe_name)
        worksheet.append(list(frame.columns))
        for row in frame.itertuples(index=False, name=None):
            worksheet.append(list(row))
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        for cell in worksheet[1]:
            cell.fill = header_fill
            cell.font = Font(color="FFFFFF", bold=True)
            cell.alignment = Alignment(vertical="center")
        for column in worksheet.columns:
            values = [str(cell.value) if cell.value is not None else "" for cell in column]
            worksheet.column_dimensions[get_column_letter(column[0].column)].width = min(
                max(len(value) for value in values) + 2,
                55,
            )
    workbook.save(WORKBOOK)


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)
    BACKUP.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE, BACKUP)
    document = Document(SOURCE)

    support_summary = json.loads(
        (TABLE_DIR / "cd_uc_cross_ancestry_summary.json").read_text(
            encoding="utf-8"
        )
    )
    pqtl_summary = json.loads(
        (TABLE_DIR / "pqtl_translation_summary.json").read_text(encoding="utf-8")
    )
    sun_pqtl_mr = pd.read_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_mr_ibd.csv"
    ).set_index("outcome")
    sun_pqtl_coloc = pd.read_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_coloc.csv"
    ).set_index("outcome")
    sun_pqtl_summary = json.loads(
        (TABLE_DIR / "prkcb_sun2018_pqtl_summary.json").read_text(encoding="utf-8")
    )
    cross_mr = pd.read_csv(TABLE_DIR / "prkcb_cross_ancestry_mr.csv")
    treatment_meta = pd.read_csv(
        TABLE_DIR / "treatment_response_meta_analysis.csv"
    ).set_index("feature")
    treatment_perm = pd.read_csv(
        TABLE_DIR / "treatment_response_permutation.csv"
    ).set_index("feature")
    readiness = pd.read_csv(
        TABLE_DIR / "clinical_development_readiness_matrix.csv"
    ).sort_values(["readiness_score_0_10", "priority_rank"], ascending=[False, True])
    readiness_ranks = pd.read_csv(
        TABLE_DIR / "clinical_readiness_rank_sensitivity.csv"
    ).set_index("target").loc["PRKCB"]
    robustness = json.loads(
        (TABLE_DIR / "clinical_translation_robustness_summary.json").read_text(
            encoding="utf-8"
        )
    )
    incremental_benchmark = json.loads(
        (TABLE_DIR / "benchmark_incremental_model_summary.json").read_text(
            encoding="utf-8"
        )
    )
    incremental_primary = incremental_benchmark["prespecified_by_task"]
    incremental_any = incremental_primary["External genetic support, any cohort"]
    incremental_both = incremental_primary["External genetic support, both cohorts"]
    incremental_eqtl = incremental_primary["Exact-variant eQTL, replicated"]
    more_than_pip_paragraph = (
        "Repeated grouped out-of-fold models compared PIP alone with PIP plus "
        "L2G, Atlas, live specificity and macrophage-context features. "
        f"The analysis used {incremental_benchmark['repetitions']} repeats of "
        f"{incremental_benchmark['splits']}-fold splits and locus-block bootstrap "
        "confidence intervals. For external support in any cohort, the "
        "prespecified PIP+context model changed AUROC by "
        f"{incremental_any['delta_auroc']:.3f} (95% CI "
        f"{incremental_any['delta_auroc_ci_low']:.3f} to "
        f"{incremental_any['delta_auroc_ci_high']:.3f}). "
        "The changes were "
        f"{incremental_both['delta_auroc']:.3f} (95% CI "
        f"{incremental_both['delta_auroc_ci_low']:.3f} to "
        f"{incremental_both['delta_auroc_ci_high']:.3f}) "
        "for support in both cohorts and "
        f"{incremental_eqtl['delta_auroc']:.3f} (95% CI "
        f"{incremental_eqtl['delta_auroc_ci_low']:.3f} to "
        f"{incremental_eqtl['delta_auroc_ci_high']:.3f}) "
        "for replicated exact-variant eQTL. "
        "The increase was therefore endpoint-specific and strongest at the broad "
        "external-support screen."
    )
    discussion_benchmark_paragraph = (
        "Repeated grouped out-of-fold benchmarking showed a prespecified gain "
        "over PIP for external support in any cohort, but not for support in both "
        "cohorts or replicated exact-variant eQTL. The framework therefore adds "
        "value as a broad external-evidence screen and retains component-level "
        "interpretation for stricter endpoints."
    )
    mechanism_limitation_paragraph = (
        "Candidate-level mechanisms. PRKCB combines three approximately "
        "independent MR instruments, regional heterogeneity, monocyte eQTL and "
        "matched-credible-set colocalization. PRKCQ, IL23R and CAMK2A add "
        "transcript, pathway and colocalization evidence. No public "
        "exact-variant MPRA or macrophage CRISPRi covers the PRKCB instruments. "
        "Allele-resolved perturbation and transcript-verified CRISPRi will define "
        "the causal step for each candidate."
    )
    for paragraph in document.paragraphs:
        text = paragraph.text.strip()
        if text.startswith("External benchmarks were constructed"):
            replace_paragraph_text(
                paragraph,
                text
                + " Repeated grouped out-of-fold models compared PIP alone with "
                "PIP plus L2G, Atlas, live specificity and macrophage-context "
                "features.",
            )
        elif text.startswith("For external support in both genetic cohorts"):
            replace_paragraph_text(paragraph, more_than_pip_paragraph)
        elif text.startswith("Benchmarking defined the incremental performance"):
            replace_paragraph_text(paragraph, discussion_benchmark_paragraph)
        elif text.startswith("Candidate-level mechanisms."):
            replace_paragraph_text(paragraph, mechanism_limitation_paragraph)
    prkcb_network = treatment_meta.loc["PRKCB_network_program"]
    macrophage_meta = treatment_meta.loc["macrophage_marker_program"]
    direct_prkcb = treatment_meta.loc["PRKCB_expression"]
    prkcb_network_perm = treatment_perm.loc[
        "PRKCB_network_program", "permutation_p_two_sided"
    ]
    macrophage_perm = treatment_perm.loc[
        "macrophage_marker_program", "permutation_p_two_sided"
    ]
    prkcb_score = float(
        readiness.loc[readiness["target"].eq("PRKCB"), "readiness_score_0_10"].iloc[0]
    )
    second_score = float(readiness["readiness_score_0_10"].iloc[1])

    methods_anchor = next(
        paragraph
        for paragraph in document.paragraphs
        if paragraph.text.strip() == "2.10 Reproducibility"
    )
    methods_elements = [
        new_heading(document, "2.10 Clinical translation extension"),
        new_paragraph(
            document,
            "Cross-ancestry validation used East Asian IBD GWAS GCST90432182 "
            "and the combined, Crohn and UC summary statistics. PRKCB cis-MR "
            "was recalculated in each cohort with the eQTL instruments.",
        ),
        new_paragraph(
            document,
            "Disease-specific target support was summarized from Crohn and UC "
            "credible-set L2G predictions. The primary pQTL analysis used the "
            "Sun 2018 plasma pQTL cis-window in GRCh38. Regional colocalization "
            "used shared variants. Older SomaScan and curated pQTL results were "
            "retained as sensitivity analyses.",
        ),
        new_paragraph(
            document,
            "Public treatment response used GSE14580 and GSE111761. PRKCB "
            "network, macrophage marker and target-expression scores were "
            "standardized within each cohort. Standardized mean differences "
            "between responders and nonresponders were combined by random-effects "
            "meta-analysis.",
        ),
        new_paragraph(
            document,
            "A 0-10 translational research priority matrix combined genetic "
            "transport, disease-specific support, eQTL, colocalization, public "
            "perturbation, treatment response, IBD clinical stage and direction "
            "compatibility. A strong one-instrument Sun 2018 pQTL layer was "
            "included; legacy pQTL data were capped at 0.25 point as sensitivity.",
        ),
    ]
    # The original reproducibility subsection follows the new methods block.
    for element in methods_elements:
        move_before(methods_anchor, element)
    for paragraph in document.paragraphs:
        if paragraph.text.strip() == "2.10 Reproducibility":
            for child in list(paragraph._p):
                if child.tag != qn("w:pPr"):
                    paragraph._p.remove(child)
            run = paragraph.add_run("2.11 Reproducibility")
            run.bold = True
            run.font.name = "Times New Roman"
            run.font.size = Pt(12)
            break

    discussion_anchor = next(
        paragraph
        for paragraph in document.paragraphs
        if paragraph.text.strip() == "4. Discussion"
    )
    extension_elements = [
        new_heading(document, "3.15 Clinical translation extension"),
        new_paragraph(
            document,
            "Disease-stratified target support was present for "
            f"{support_summary['cd_supported_targets_l2g_ge_0_5']} targets in "
            "Crohn disease and "
            f"{support_summary['uc_supported_targets_l2g_ge_0_5']} targets in UC, with "
            f"{support_summary['both_disease_targets']} targets supported in both "
            "diseases (Table 12).",
        ),
        new_paragraph(
            document,
            "PRKCB cis-MR was positive across combined IBD, Crohn, UC and East "
            "Asian IBD cohorts. Estimates ranged from "
            f"{cross_mr['beta'].min():.3f} to {cross_mr['beta'].max():.3f} "
            "log OR per "
            "expression SD, and all four tests reached P < 0.001 (Figure 13A-B). "
            "Leave-one-instrument-out estimates remained positive.",
        ),
        new_paragraph(
            document,
            "The Sun 2018 plasma pQTL provided "
            f"{sun_pqtl_summary['pqtl_variants']:,} PRKCB cis variants in GRCh38. "
            "One LD-independent instrument was available. pQTL MR was positive in "
            "all four IBD outcomes. The combined-IBD estimate was "
            f"{sun_pqtl_mr.loc['IBD_combined', 'beta']:.3f} "
            f"(P = {sun_pqtl_mr.loc['IBD_combined', 'p']:.2e}). Crohn, UC and "
            "East Asian IBD were also positive. Colocalization reached "
            f"H4 = {sun_pqtl_coloc.loc['IBD_combined', 'H4']:.3f} for combined IBD "
            f"and H4 = {sun_pqtl_coloc.loc['UC', 'H4']:.3f} for UC. Crohn and East "
            f"Asian H4 values were {sun_pqtl_coloc.loc['Crohn', 'H4']:.3f} and "
            f"{sun_pqtl_coloc.loc['IBD_EAS', 'H4']:.3f}. The older SomaScan and "
            "curated deCODE panel remained a sensitivity analysis.",
        ),
        new_paragraph(
            document,
            "A public-resource recheck found no exact-variant MPRA or macrophage "
            "CRISPRi data for the three PRKCB instruments. PRKCB remains a "
            "mechanistically triangulated candidate; allele-resolved PRKCB "
            "perturbation is the decisive external experiment.",
        ),
        new_paragraph(
            document,
            "In public anti-TNF datasets, the PRKCB network program was lower in "
            "responders than nonresponders "
            f"(meta-analysis SMD = {prkcb_network['smd']:.2f}; 95% CI "
            f"{prkcb_network['ci_low']:.2f} to {prkcb_network['ci_high']:.2f}; "
            f"P = {prkcb_network['p']:.3f}). The macrophage marker program showed the "
            "same direction "
            f"(SMD = {macrophage_meta['smd']:.2f}; 95% CI "
            f"{macrophage_meta['ci_low']:.2f} to {macrophage_meta['ci_high']:.2f}; "
            f"P = {macrophage_meta['p']:.3f}; Figure 13C). A 20,000-label "
            f"permutation gave P = {prkcb_network_perm:.4f} for the "
            f"PRKCB network program and P = {macrophage_perm:.4f} for the macrophage marker "
            "program. Direct PRKCB expression showed the same direction but a "
            f"wider interval (SMD = {direct_prkcb['smd']:.2f}; 95% CI "
            f"{direct_prkcb['ci_low']:.2f} to {direct_prkcb['ci_high']:.2f}).",
        ),
        new_paragraph(
            document,
            f"The readiness matrix ranked PRKCB highest at {prkcb_score:.1f} of 10, "
            "followed by PDGFB and MAP3K8 at "
            f"{second_score:.1f} (Table 12; Figure 13D). PRKCB combined "
            "cross-ancestry genetics, colocalization, a strong one-instrument "
            "pQTL signal, "
            "treatment-response programs and an IBD clinical-stage signal with "
            "compatible pharmacological direction. PRKCB ranked "
            f"{int(readiness_ranks['readiness_score_0_10_rank'])} under the "
            "primary scheme, "
            f"{int(readiness_ranks['evidence_only_score_rank'])} under "
            "evidence-only scoring, "
            f"{int(readiness_ranks['clinical_weighted_score_rank'])} under "
            "clinical-weighted scoring and "
            f"{int(readiness_ranks['no_treatment_score_rank'])} when treatment "
            "response was excluded. Across "
            f"{robustness['readiness_weight_draws']:,} randomized weightings, PRKCB "
            f"ranked first in {robustness['prkcb_probability_rank_1']:.1%} and top two "
            f"in {robustness['prkcb_probability_top_2']:.1%} of draws.",
        ),
        new_paragraph(
            document,
            "Table 12. Clinical development readiness matrix for Tier 1 targets",
        ),
    ]
    table = add_readiness_table(document)
    extension_elements.append(table)
    extension_elements.append(
        new_paragraph(
            document,
            "The score is an evidence-integration summary. Components are genetic "
            "transport, disease-specific support, exact-variant eQTL, "
            "colocalization, public perturbation, pQTL, treatment-response program, "
            "IBD clinical stage and direction compatibility. The strong pQTL "
            "contribution uses one LD-independent instrument.",
        )
    )
    figure_paragraph = add_figure(document)
    extension_elements.append(figure_paragraph)
    extension_elements.append(
        new_paragraph(
            document,
            "Figure 13. Clinical translation extension. (A) Crohn and UC "
            "target-level L2G support. (B) PRKCB cis-MR across IBD, Crohn, UC and "
            "East Asian cohorts. (C) Public anti-TNF response-program "
            "meta-analysis. (D) Clinical development readiness matrix with "
            "weight-bootstrap rank probabilities.",
        )
    )
    for element in extension_elements:
        move_before(discussion_anchor, element)

    document.save(OUTPUT)
    MANUSCRIPT_COPY.write_bytes(OUTPUT.read_bytes())
    append_workbook_sheets()

    revised = Document(OUTPUT)
    paragraphs = [paragraph.text.strip() for paragraph in revised.paragraphs]
    audit = {
        "output": str(OUTPUT),
        "copy": str(MANUSCRIPT_COPY),
        "paragraphs": len(revised.paragraphs),
        "tables": len(revised.tables),
        "figures": len(revised.inline_shapes),
        "has_methods_extension": "2.10 Clinical translation extension" in paragraphs,
        "has_results_extension": "3.15 Clinical translation extension" in paragraphs,
        "has_table12_caption": "Table 12. Clinical development readiness matrix for Tier 1 targets"
        in paragraphs,
        "has_figure13_caption": any(
            text.startswith("Figure 13.") for text in paragraphs
        ),
        "all_checks_pass": all(
            [
                "2.10 Clinical translation extension" in paragraphs,
                "3.15 Clinical translation extension" in paragraphs,
                "Table 12. Clinical development readiness matrix for Tier 1 targets"
                in paragraphs,
                any(text.startswith("Figure 13.") for text in paragraphs),
            ]
        ),
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("Clinical translation manuscript extension failed.")


if __name__ == "__main__":
    main()
