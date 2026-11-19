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
        "Drug stage",
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
        rows.append(
            [
                row["target"],
                support,
                "Yes" if row["exact_eqtl_score"] else "No",
                "Yes" if row["colocalization_score"] else "No",
                "Yes" if row["public_perturbation_score"] else "No",
                "Yes" if row["pqtl_mr_score"] else "No",
                "Yes" if row["treatment_response_score"] else "No",
                str(row["clinical_stage"])
                if pd.notna(row["clinical_stage"])
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
    }
    workbook = load_workbook(WORKBOOK)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    for sheet_name, frame in tables.items():
        if sheet_name in workbook.sheetnames:
            del workbook[sheet_name]
        worksheet = workbook.create_sheet(sheet_name)
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
            "credible-set L2G predictions. pQTL MR used SomaScan PRKCB datasets "
            "and a curated deCODE pQTL. Regional colocalization used shared "
            "variants. A targeted PheW-MR panel used available OpenGWAS safety "
            "phenotypes.",
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
            "A 0-10 readiness matrix combined genetic transport, disease-specific "
            "support, eQTL, colocalization, public perturbation, pQTL, treatment "
            "response, clinical stage and direction compatibility.",
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
            "Disease-stratified target support was present for five targets in "
            "Crohn disease and six targets in UC, with four targets supported in "
            "both diseases (Table 12).",
        ),
        new_paragraph(
            document,
            "PRKCB cis-MR was positive across combined IBD, Crohn, UC and East "
            "Asian IBD cohorts. Estimates ranged from 0.035 to 0.066 log OR per "
            "expression SD, and all four tests reached P < 0.001 (Figure 13A-B).",
        ),
        new_paragraph(
            document,
            "pQTL MR showed a positive PRKCB direction in nine of twelve tests. "
            "The SomaScan PKCB2 dataset supported combined IBD and UC, while the "
            "deCODE instrument supported East Asian IBD. Regional pQTL "
            "colocalization defined a graded H4 distribution for future pQTL "
            "reference expansion.",
        ),
        new_paragraph(
            document,
            "In public anti-TNF datasets, the PRKCB network program was lower in "
            "responders than nonresponders (meta-analysis SMD = -2.04; 95% CI "
            "-3.76 to -0.33; P = 0.020). The macrophage marker program showed the "
            "same direction (SMD = -1.13; 95% CI -2.11 to -0.14; P = 0.025; "
            "Figure 13C).",
        ),
        new_paragraph(
            document,
            "The readiness matrix ranked PRKCB highest at 7.5 of 10, followed by "
            "PDGFB at 5.0 and MAP3K8 at 4.0 (Table 12; Figure 13D). PRKCB combined "
            "cross-ancestry genetics, colocalization, public perturbation, pQTL "
            "direction, treatment-response programs and a direction-compatible "
            "clinical-stage signal.",
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
            "colocalization, public perturbation, pQTL MR, treatment-response "
            "program, clinical stage and direction compatibility.",
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
            "meta-analysis. (D) Clinical development readiness matrix.",
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
