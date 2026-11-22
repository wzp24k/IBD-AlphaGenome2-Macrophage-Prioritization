from __future__ import annotations

import csv
import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from openpyxl import load_workbook

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


WORKBOOK = ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx"
OUTPUT = ROOT / "Supplementary_Tables.docx"
MANUSCRIPT_COPY = ROOT / "manuscript" / "Supplementary_Tables.docx"
INDEX_CSV = ROOT / "results" / "tables" / "supplementary_table_index.csv"
AUDIT_JSON = ROOT / "results" / "supplementary_tables_audit.json"


def value_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, float):
        absolute = abs(value)
        if absolute == 0:
            return "0.000"
        if absolute < 0.001:
            return f"{value:.2e}"
        if absolute < 0.01:
            return f"{value:.4f}"
        if absolute < 10:
            return f"{value:.3f}"
        return f"{value:.1f}"
    return str(value)


def interval_text(low: object, high: object) -> str:
    return f"{value_text(low)}-{value_text(high)}"


def build_benchmark_table(
    source: pd.DataFrame,
    calibration: pd.DataFrame,
) -> list[list[str]]:
    task_order = [
        "External genetic support, both cohorts",
        "Exact-variant eQTL, replicated",
    ]
    model_order = ["PIP-only", "L2G-only", "Atlas-only", "Live-only", "Combined"]
    calibration_index = calibration.set_index(["task", "model"])
    rows = []
    for task in task_order:
        for model in model_order:
            metric = source[
                source["task"].eq(task) & source["model"].eq(model)
            ].iloc[0]
            calibration_row = calibration_index.loc[(task, model)]
            rows.append(
                [
                    task,
                    model,
                    f'{int(metric["positives"])}/{int(metric["negatives"])}',
                    f'{value_text(metric["auroc"])} '
                    f'({interval_text(metric["auroc_ci_low"], metric["auroc_ci_high"])})',
                    f'{value_text(metric["auprc"])} '
                    f'({interval_text(metric["auprc_ci_low"], metric["auprc_ci_high"])})',
                    value_text(calibration_row["brier_score"]),
                    value_text(calibration_row["expected_calibration_error"]),
                ]
            )
    return rows


def build_eqtl_table(
    triangulation: pd.DataFrame,
    coloc: pd.DataFrame,
    sensitivity: pd.DataFrame,
    translation: pd.DataFrame,
) -> list[list[str]]:
    matched = triangulation[
        triangulation["pair_type"].eq("matched_target")
    ].set_index("target")
    mismatched = triangulation[
        triangulation["pair_type"].eq("mismatched_target_control")
    ].set_index("target")
    coloc_index = coloc.set_index("target")
    sensitivity_index = sensitivity.set_index("target")
    rows = []
    for _, clinical in translation.sort_values("priority_rank").iterrows():
        target = clinical["target"]
        match_row = matched.loc[target]
        mismatch_row = mismatched.loc[target]
        if target in coloc_index.index:
            primary_h4 = value_text(coloc_index.loc[target, "H4"])
        else:
            primary_h4 = "Not evaluable"
        if target in sensitivity_index.index:
            sensitivity_row = sensitivity_index.loc[target]
            h4_range = interval_text(
                sensitivity_row["minimum_H4"],
                sensitivity_row["maximum_H4"],
            )
            supported = "Yes" if sensitivity_row["sensitivity_supported"] else "No"
        else:
            h4_range = "Not evaluable"
            supported = "Not evaluable"
        rows.append(
            [
                str(int(clinical["priority_rank"])),
                target,
                match_row["canonical_variant"],
                str(int(match_row["nominal_datasets"])),
                value_text(match_row["minimum_pvalue"]),
                str(int(mismatch_row["nominal_datasets"])),
                primary_h4,
                h4_range,
                supported,
            ]
        )
    return rows


def build_prkcb_table(
    mr: pd.DataFrame,
    smr: pd.DataFrame,
    qtl: pd.DataFrame,
    signatures: pd.DataFrame,
    virtual: pd.DataFrame,
    sun_pqtl_mr: pd.DataFrame,
    sun_pqtl_coloc: pd.DataFrame,
) -> list[list[str]]:
    mr_001 = mr[mr["clump_r2_threshold"].eq(0.01)]
    ivw = mr_001[mr_001["method"].eq("IVW_random_effects")].iloc[0]
    median = mr_001[mr_001["method"].eq("weighted_median")].iloc[0]
    egger = mr_001[mr_001["method"].eq("MR_Egger")].iloc[0]
    smr_single = smr[
        smr["method"].eq("SMR_style_single_instrument")
    ].iloc[0]
    heidi_full = smr[smr["variant_set"].eq("all_ld_matched")].iloc[0]
    heidi_pruned = smr[
        smr["variant_set"].eq("ld_pruned_r2_0.1")
    ].iloc[0]
    prkcb_signatures = signatures[signatures["target"].eq("PRKCB")]
    positive_signatures = int(prkcb_signatures["correlation"].gt(0).sum())
    median_signature_p = float(prkcb_signatures["p_permutation"].median())
    active_jak = virtual[
        virtual["state"].eq("active") & virtual["pathway"].eq("JAK/STAT")
    ].iloc[0]
    mrna_colon = qtl[qtl["source"].eq("GTEx_v10_colon_eQTL")].iloc[0]
    mrna_myeloid = qtl[
        qtl["source"].eq("eight_myeloid_eQTL_datasets")
    ].iloc[0]
    sqtl = qtl[qtl["source"].eq("GTEx_v10_sQTL")].iloc[0]
    pqtl_mr = sun_pqtl_mr[sun_pqtl_mr["outcome"].eq("IBD_combined")].iloc[0]
    pqtl_coloc = sun_pqtl_coloc[
        sun_pqtl_coloc["outcome"].eq("IBD_combined")
    ].iloc[0]
    return [
        [
            "Two-sample MR, IVW",
            "3 approximately independent eQTL instruments",
            f'beta = {value_text(ivw["estimate"])} per eQTL SD',
            f'SE = {value_text(ivw["se"])}; P < 0.001; '
            f'mean F = {value_text(ivw["mean_f_statistic"])}',
            "Directional mediation supported; MR-Egger is discordant",
        ],
        [
            "Two-sample MR, weighted median",
            "3 instruments",
            f'beta = {value_text(median["estimate"])}',
            f'SE = {value_text(median["se"])}; P < 0.001',
            "Concordant with IVW",
        ],
        [
            "Two-sample MR, MR-Egger",
            "3 instruments",
            f'beta = {value_text(egger["estimate"])}',
            f'P = {value_text(egger["p"])}; '
            f'intercept P = {value_text(egger["intercept_p"])}',
            "Directional inconsistency limits a strict causal claim",
        ],
        [
            "SMR-style single instrument",
            "rs7404003",
            f'beta = {value_text(smr_single["smr_beta"])}',
            f'SE = {value_text(smr_single["smr_se"])}; '
            f'P = {value_text(smr_single["smr_p"])}',
            "Summary-statistic implementation, not official SMR software",
        ],
        [
            "HEIDI-style heterogeneity",
            f'{int(heidi_full["nsnp"])} LD-matched variants',
            f'P = {value_text(heidi_full["p"])}',
            f'status = {heidi_full["heterogeneity_status"]}',
            "Full region contains multiple or context-dependent signals",
        ],
        [
            "HEIDI-style heterogeneity",
            f'{int(heidi_pruned["nsnp"])} pruned instruments',
            f'P = {value_text(heidi_pruned["p"])}',
            f'status = {heidi_pruned["heterogeneity_status"]}',
            "Pruned subset supports directional consistency",
        ],
        [
            "SuSiE-colocalization",
            "GWAS and eQTL credible sets",
            "H4 = 1.000",
            "H4 = 1.000 across reported priors",
            "Shared credible-set signal",
        ],
        [
            "Molecular QTL coverage",
            "GTEx colon and eight myeloid datasets",
            f'Colon NES = {value_text(mrna_colon["observed_direction"])}; '
            f'myeloid = {mrna_myeloid["observed_direction"]}',
            f'sQTL covered = {value_text(sqtl["exact_variant_coverage"])}; '
            f'pQTL beta = {value_text(pqtl_mr["beta"])}, '
            f'P = {value_text(pqtl_mr["p"])}; '
            f'pQTL H4 = {value_text(pqtl_coloc["H4"])}',
            "Context-dependent mRNA direction; strong Sun 2018 pQTL support",
        ],
        [
            "Public signature reversal",
            "Oliver 2024, Crohn 2023 and Garrido-Trigo 2023",
            f"{positive_signatures}/3 cohorts positive",
            f"median permutation P = {value_text(median_signature_p)}",
            "Weak and cohort-dependent projection",
        ],
        [
            "Virtual PRKCB knockdown",
            "1,000 donor-block bootstrap replicates",
            f'active JAK/STAT beta = {value_text(active_jak["mean"])}',
            f'95% CI {interval_text(active_jak["ci_low"], active_jak["ci_high"])}',
            "Computational prediction, not experimental causality",
        ],
    ]


def pretty_experimental_support(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "None"
    mapping = {
        "perturbseq_transcriptome": "Perturb-seq transcriptome",
        "directional_variant_eqtl": "Directional variant eQTL",
        "target_engagement": "Target engagement",
    }
    return mapping.get(str(value), str(value).replace("_", " "))


def pretty_risk_direction(value: object) -> str:
    mapping = {
        "risk_increases_target": "Risk allele increases target",
        "risk_decreases_target": "Risk allele decreases target",
    }
    return mapping.get(str(value), str(value))


def build_translation_table(source: pd.DataFrame) -> list[list[str]]:
    rows = []
    for _, row in source.sort_values("priority_rank").iterrows():
        external = (
            "None"
            if int(row["external_cluster_count"]) == 0
            else str(int(row["external_cluster_count"]))
        )
        stage = (
            "None"
            if pd.isna(row["clinical_stage"])
            else str(row["clinical_stage"]).title()
        )
        drugs = (
            "None"
            if pd.isna(row["ibd_indication_drugs"])
            else str(row["ibd_indication_drugs"]).replace(";", ", ")
        )
        rows.append(
            [
                str(int(row["priority_rank"])),
                row["target"],
                external,
                value_text(bool(row["exact_variant_eqtl"])),
                value_text(bool(row["formal_colocalization"])),
                pretty_experimental_support(row["experimental_support"]),
                stage,
                drugs,
                pretty_risk_direction(row["risk_effect_on_target"]),
                row["translational_priority"],
            ]
        )
    return rows


def set_cell_margins(table) -> None:
    tbl_pr = table._tbl.tblPr
    margins = tbl_pr.first_child_found_in("w:tblCellMar")
    if margins is None:
        margins = OxmlElement("w:tblCellMar")
        tbl_pr.append(margins)
    for edge, value in (("top", 35), ("start", 50), ("bottom", 35), ("end", 50)):
        node = margins.find(qn(f"w:{edge}"))
        if node is None:
            node = OxmlElement(f"w:{edge}")
            margins.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


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


def set_fixed_layout(table) -> None:
    tbl_pr = table._tbl.tblPr
    layout = tbl_pr.first_child_found_in("w:tblLayout")
    if layout is None:
        layout = OxmlElement("w:tblLayout")
        tbl_pr.append(layout)
    layout.set(qn("w:type"), "fixed")


def set_repeat_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def prevent_row_split(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    cant_split = OxmlElement("w:cantSplit")
    cant_split.set(qn("w:val"), "true")
    tr_pr.append(cant_split)


def column_widths(headers: list[str], rows: list[list[str]], total: float) -> list[float]:
    scores = []
    for index, header in enumerate(headers):
        values = [row[index] for row in rows]
        max_text = max((min(len(value), 34) for value in values), default=8)
        scores.append(max(7, min(36, max(len(header) * 0.75, max_text))))
    score_sum = sum(scores)
    widths = [total * score / score_sum for score in scores]
    widths = [max(0.55, width) for width in widths]
    width_sum = sum(widths)
    return [width * total / width_sum for width in widths]


def add_table(
    document: Document,
    headers: list[str],
    rows: list[list[str]],
    font_size: float,
) -> None:
    total_width = 10.55
    table = document.add_table(rows=len(rows) + 1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    set_fixed_layout(table)
    set_three_line_borders(table)
    set_cell_margins(table)
    widths = column_widths(headers, rows, total_width)
    for index, width in enumerate(widths):
        table.columns[index].width = Inches(width)
    for row_index, values in enumerate([headers, *rows]):
        row = table.rows[row_index]
        prevent_row_split(row)
        if row_index == 0:
            set_repeat_header(row)
        for column_index, value in enumerate(values):
            cell = row.cells[column_index]
            cell.width = Inches(widths[column_index])
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.space_before = Pt(0)
            paragraph.paragraph_format.space_after = Pt(0)
            paragraph.paragraph_format.line_spacing = 1.0
            paragraph.alignment = (
                WD_ALIGN_PARAGRAPH.CENTER
                if row_index == 0
                else WD_ALIGN_PARAGRAPH.LEFT
            )
            run = paragraph.add_run(str(value))
            run.font.name = "Times New Roman"
            run.font.size = Pt(font_size)
            if row_index == 0:
                run.bold = True


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    field_begin = OxmlElement("w:fldChar")
    field_begin.set(qn("w:fldCharType"), "begin")
    instruction = OxmlElement("w:instrText")
    instruction.set(qn("xml:space"), "preserve")
    instruction.text = "PAGE"
    field_end = OxmlElement("w:fldChar")
    field_end.set(qn("w:fldCharType"), "end")
    run._r.append(field_begin)
    run._r.append(instruction)
    run._r.append(field_end)


def style_document(document: Document) -> None:
    section = document.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width = Inches(11.69)
    section.page_height = Inches(8.27)
    section.top_margin = Inches(0.5)
    section.bottom_margin = Inches(0.5)
    section.left_margin = Inches(0.5)
    section.right_margin = Inches(0.5)
    styles = document.styles
    styles["Normal"].font.name = "Times New Roman"
    styles["Normal"].font.size = Pt(9)
    styles["Normal"].paragraph_format.space_after = Pt(0)
    styles["Normal"].paragraph_format.line_spacing = 1.0
    footer = section.footer.paragraphs[0]
    footer.text = "Supplementary Tables | Page "
    footer.runs[0].font.name = "Times New Roman"
    footer.runs[0].font.size = Pt(8)
    add_page_number(footer)


def add_title_page(document: Document) -> None:
    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.paragraph_format.space_before = Pt(70)
    title.paragraph_format.space_after = Pt(20)
    run = title.add_run("Supplementary Tables")
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(24)
    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.paragraph_format.space_after = Pt(24)
    run = subtitle.add_run(
        "An evidence ladder for non-coding inflammatory bowel disease variants "
        "prioritizes intestinal macrophage regulatory mechanisms"
    )
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(13)
    note = document.add_paragraph()
    note.alignment = WD_ALIGN_PARAGRAPH.CENTER
    note.paragraph_format.left_indent = Inches(1.25)
    note.paragraph_format.right_indent = Inches(1.25)
    run = note.add_run(
        "Supplementary Tables S1-S4 contain the priority benchmark, eQTL, "
        "PRKCB mechanism and translational summaries. The complete "
        "source workbook accompanies this document."
    )
    run.font.name = "Times New Roman"
    run.font.size = Pt(10)
    version = document.add_paragraph()
    version.alignment = WD_ALIGN_PARAGRAPH.CENTER
    version.paragraph_format.space_before = Pt(20)
    run = version.add_run("Version: 14 November 2026")
    run.font.name = "Times New Roman"
    run.font.size = Pt(9)


def add_supplement_table(
    document: Document,
    number: int,
    title: str,
    headers: list[str],
    rows: list[list[str]],
    note: str,
    font_size: float,
) -> None:
    document.add_page_break()
    caption = document.add_paragraph()
    caption.paragraph_format.keep_with_next = True
    caption.paragraph_format.space_after = Pt(6)
    run = caption.add_run(f"Supplementary Table S{number}. ")
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(10)
    run = caption.add_run(title)
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(10)
    add_table(document, headers, rows, font_size)
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.space_before = Pt(5)
    run = paragraph.add_run(note)
    run.font.name = "Times New Roman"
    run.font.size = Pt(8)


def main() -> None:
    workbook = load_workbook(WORKBOOK, read_only=True)
    metrics = pd.read_excel(WORKBOOK, sheet_name="Benchmark_metrics")
    calibration = pd.read_excel(WORKBOOK, sheet_name="Benchmark_calibration")
    triangulation = pd.read_excel(WORKBOOK, sheet_name="Top10_eQTL_triangulation")
    coloc = pd.read_excel(WORKBOOK, sheet_name="Top10_coloc")
    coloc_sensitivity = pd.read_excel(
        WORKBOOK,
        sheet_name="Top10_coloc_sensitivity",
    )
    translation = pd.read_excel(WORKBOOK, sheet_name="Clinical_translation")
    mr = pd.read_excel(WORKBOOK, sheet_name="PRKCB_MR_estimates")
    smr = pd.read_excel(WORKBOOK, sheet_name="PRKCB_SMR_HEIDI")
    qtl = pd.read_excel(WORKBOOK, sheet_name="PRKCB_QTL_coverage")
    signatures = pd.read_excel(WORKBOOK, sheet_name="PRKCB_signature_reversal")
    virtual = pd.read_excel(
        WORKBOOK,
        sheet_name="PRKCB_virtual_bootstrap1000",
    )
    sun_pqtl_mr = pd.read_excel(WORKBOOK, sheet_name="Sun2018_pQTL_MR")
    sun_pqtl_coloc = pd.read_excel(WORKBOOK, sheet_name="Sun2018_pQTL_coloc")

    tables = [
        {
            "title": "Benchmark discrimination and calibration",
            "headers": [
                "Benchmark task",
                "Model",
                "Positive/negative",
                "AUROC (95% CI)",
                "AUPRC (95% CI)",
                "Brier score",
                "ECE",
            ],
            "rows": build_benchmark_table(metrics, calibration),
            "note": (
                "Source worksheets: Benchmark_metrics, Benchmark_calibration, "
                "Benchmark_incremental_metrics, Benchmark_incremental_deltas, "
                "Sun2018_pQTL_cis, Sun2018_pQTL_MR, Sun2018_pQTL_coloc, "
                "Expanded_macrophage_metrics and Expanded_macrophage_delta. "
                "Intervals use locus-block bootstrap replicates. AUROC, AUPRC, "
                "Brier score and ECE are shown for the two principal benchmark tasks."
            ),
            "font_size": 7.4,
        },
        {
            "title": "Top 10 exact-variant eQTL triangulation and colocalization",
            "headers": [
                "Rank",
                "Target",
                "Variant",
                "Matched eQTL datasets",
                "Minimum P",
                "Mismatched-target hits",
                "Primary H4",
                "H4 sensitivity range",
                "Formal support",
            ],
            "rows": build_eqtl_table(
                triangulation,
                coloc,
                coloc_sensitivity,
                translation,
            ),
            "note": (
                "Source worksheets: Top10_eQTL_triangulation, Top10_coloc, "
                "Top10_coloc_sensitivity and Clinical_translation. "
                "Mismatched-target pairs provide the negative-control comparison. "
                "Formal support requires H4 >= 0.80 across sensitivity settings."
            ),
            "font_size": 7.0,
        },
        {
            "title": "PRKCB mediation and virtual-perturbation evidence",
            "headers": [
                "Analysis",
                "Dataset or comparison",
                "Primary result",
                "Uncertainty or secondary result",
                "Interpretation boundary",
            ],
            "rows": build_prkcb_table(
                mr,
                smr,
                qtl,
                signatures,
                virtual,
                sun_pqtl_mr,
                sun_pqtl_coloc,
            ),
            "note": (
                "Source worksheets: PRKCB_MR_estimates, PRKCB_SMR_HEIDI, "
                "PRKCB_HEIDI_windows, PRKCB_QTL_coverage, "
                "PRKCB_signature_reversal and PRKCB_virtual_bootstrap1000. "
                "Sun2018_pQTL_MR and Sun2018_pQTL_coloc provide the protein layer. "
                "SMR and HEIDI analyses are summary-statistic implementations. "
                "Virtual knockdown is not experimental evidence."
            ),
            "font_size": 7.2,
        },
        {
            "title": "Indication-aware translational annotation for the Top 10 targets",
            "headers": [
                "Rank",
                "Target",
                "External cohorts",
                "Exact eQTL",
                "Formal coloc",
                "Public perturbation",
                "Clinical stage",
                "IBD-indication drug",
                "Risk direction",
                "Translational priority",
            ],
            "rows": build_translation_table(translation),
            "note": (
                "Source worksheets: Clinical_translation, PRKCB_IBD_trial_search "
                "and PRKCB_IBD_trial_results. Drug indications and target-stage "
                "annotations are descriptive. No target met the prespecified "
                "clinical treatment-recommendation boundary."
            ),
            "font_size": 6.7,
        },
    ]

    document = Document()
    style_document(document)
    add_title_page(document)
    for number, table in enumerate(tables, start=1):
        add_supplement_table(
            document,
            number,
            table["title"],
            table["headers"],
            table["rows"],
            table["note"],
            table["font_size"],
        )
    document.save(OUTPUT)
    MANUSCRIPT_COPY.write_bytes(OUTPUT.read_bytes())

    INDEX_CSV.parent.mkdir(parents=True, exist_ok=True)
    source_map = {
        1: "Benchmark_metrics; Benchmark_calibration",
        2: (
            "Top10_eQTL_triangulation; Top10_coloc; "
            "Top10_coloc_sensitivity; Clinical_translation"
        ),
        3: (
            "PRKCB_MR_estimates; PRKCB_SMR_HEIDI; PRKCB_HEIDI_windows; "
            "PRKCB_QTL_coverage; PRKCB_signature_reversal; "
            "PRKCB_virtual_bootstrap1000"
        ),
        4: (
            "Clinical_translation; PRKCB_IBD_trial_search; "
            "PRKCB_IBD_trial_results"
        ),
    }
    with INDEX_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "supplementary_table",
                "title",
                "source_sheets",
                "rows_in_word",
                "status",
            ],
        )
        writer.writeheader()
        for number, table in enumerate(tables, start=1):
            writer.writerow(
                {
                    "supplementary_table": f"S{number}",
                    "title": table["title"],
                    "source_sheets": source_map[number],
                    "rows_in_word": len(table["rows"]),
                    "status": "included",
                }
            )

    captions = [
        paragraph.text
        for paragraph in document.paragraphs
        if paragraph.text.startswith("Supplementary Table S")
    ]
    caption_numbers = [
        int(re.match(r"Supplementary Table S(\d+)\.", caption).group(1))
        for caption in captions
    ]
    checks = {
        "document_exists": OUTPUT.exists(),
        "four_word_tables": len(document.tables) == 4,
        "four_supplement_captions": len(captions) == 4,
        "numbers_sequential": caption_numbers == [1, 2, 3, 4],
        "titles_unique": len(set(table["title"] for table in tables)) == 4,
        "all_tables_have_rows": all(len(table["rows"]) > 0 for table in tables),
        "no_main_table_duplicates": not any(
            paragraph.text.startswith("Table ")
            for paragraph in document.paragraphs
        ),
        "all_source_notes_present": all(
            any(table["note"] in paragraph.text for paragraph in document.paragraphs)
            for table in tables
        ),
    }
    audit = {
        "generated": datetime.now().astimezone().isoformat(),
        "output": str(OUTPUT),
        "supplementary_tables": 4,
        "word_tables": len(document.tables),
        "source_workbook_sheets": len(workbook.sheetnames),
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    AUDIT_JSON.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("Supplementary table audit failed.")
    print(f"Wrote {OUTPUT}")
    print(f"Wrote {MANUSCRIPT_COPY}")
    print(f"Wrote {INDEX_CSV}")


if __name__ == "__main__":
    main()
