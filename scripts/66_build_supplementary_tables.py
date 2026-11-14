from __future__ import annotations

import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path

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


@dataclass(frozen=True)
class TableSpec:
    source: str
    title: str
    note: str
    columns: tuple[str, ...] | None = None
    sort_by: str | None = None
    descending: bool = False
    max_rows: int | None = None


SPECS = [
    TableSpec(
        "Overview",
        "Analysis summary metrics",
        "Values are the primary counts and analysis settings retained in the source workbook.",
    ),
    TableSpec(
        "Live_summary",
        "Live AlphaGenome 2 prediction inventory",
        "Counts summarize returned predictions by output type, biosample and control role.",
    ),
    TableSpec(
        "Controls",
        "Live-model positive and negative control predictions",
        "Control effects are prediction outputs and are not count-based hypothesis tests.",
    ),
    TableSpec(
        "GTEx_PRKCB_eQTL",
        "GTEx v10 PRKCB colon eQTL",
        "NES is the normalized effect size for the ALT allele in sigmoid colon.",
        columns=(
            "snpId",
            "variant",
            "target",
            "tissueSiteDetailId",
            "pValue",
            "nes",
            "datasetId",
            "experiment",
        ),
    ),
    TableSpec(
        "ETS2_macrophage_MPRA",
        "ETS2 positive-control macrophage MPRA",
        "Effects are mean log2 risk-versus-nonrisk reporter output across RNA pools.",
    ),
    TableSpec(
        "THP1_OE_engagement",
        "THP-1 overexpression target engagement",
        "Self-expression deltas compare each overexpressed gene with empty-vector controls.",
    ),
    TableSpec(
        "eQTL_Catalogue_hits",
        "Nominal eQTL Catalogue hits",
        "Rows report nominal macrophage or monocyte associations queried for Tier 1 targets.",
    ),
    TableSpec(
        "eQTL_Catalogue_coverage",
        "eQTL Catalogue dataset coverage",
        "Missing coverage is classified as absence of a usable public record, not a negative association.",
    ),
    TableSpec(
        "PerturbSeq_targets",
        "Macrophage Perturb-seq target-transcript effects",
        "Values compare guide pseudobulk profiles with matched non-targeting controls.",
    ),
    TableSpec(
        "PerturbSeq_DE_summary",
        "Macrophage Perturb-seq differential-expression summary",
        "Responsive-gene counts require FDR < 0.05; absolute-effect counts use an additional 0.5 threshold.",
    ),
    TableSpec(
        "PerturbSeq_pathways",
        "Macrophage Perturb-seq pathway effects",
        "Pathway tests used equal-size empirical resampling from expressed genes.",
    ),
    TableSpec(
        "PerturbSeq_mixed_models",
        "Donor-aware Perturb-seq models",
        "The fallback label indicates a donor fixed-effect model when the random-intercept variance was singular.",
    ),
    TableSpec(
        "CAMK2A_coloc_sensitivity",
        "CAMK2A colocalization sensitivity",
        "Posterior probabilities vary the prior standard deviation and shared-causal-variant prior.",
    ),
    TableSpec(
        "Gut_scATAC_variants",
        "Gut single-cell ATAC annotation of Tier 1 variants",
        "Chromatin overlap is regulatory annotation, not allele-specific perturbation.",
        columns=(
            "target",
            "rsid",
            "chromosome",
            "position",
            "finemap_susie_pip",
            "overlapped_cCRE",
            "cCRE_celltype",
            "ABC_max_gene",
            "nearest_gene",
            "interpretation",
        ),
    ),
    TableSpec(
        "scATAC_topic_percentiles",
        "scATAC topic percentiles for the PDGFB regulatory element",
        "Percentiles rank the element score among the available gut chromatin topics.",
    ),
    TableSpec(
        "GSEA_growth_sensitivity",
        "GSEA sensitivity to proliferation residualization",
        "Residualized enrichment removes a sample-level E2F, G2-M and MYC principal component.",
    ),
    TableSpec(
        "PRKCQ_repeated_summary",
        "PRKCQ repeated-pool stability summary",
        "Repeated pools quantify sampling stability within available donors and are not biological replicates.",
    ),
    TableSpec(
        "External_GWAS_studies",
        "External IBD GWAS cohort and credible-set summary",
        "The UK Biobank and Million Veteran Program clusters were counted as non-overlapping external evidence.",
    ),
    TableSpec(
        "External_GWAS_bootstrap",
        "External IBD pathway locus-block bootstrap",
        "Intervals and FDR values use locus-block resampling of credible-set loci.",
    ),
    TableSpec(
        "Pathway_weighted_tests",
        "Pathway member-score comparisons",
        "Tests compare member and non-member target scores in the external-weighted gene ranking.",
        columns=(
            "pathway",
            "library",
            "members_in_background",
            "median_member_score",
            "median_nonmember_score",
            "p_value",
            "fdr_bh",
        ),
    ),
    TableSpec(
        "All_variant_pathway",
        "Variant-set pathway sensitivity",
        "Rows compare pathway scores under alternative supported-variant definitions.",
        columns=(
            "method",
            "pathway",
            "library",
            "members_in_background",
            "median_member_score",
            "median_nonmember_score",
            "p_value",
            "fdr_bh",
        ),
    ),
    TableSpec(
        "Benchmark_metrics",
        "Benchmark discrimination metrics",
        "AUROC and AUPRC intervals use locus-block bootstrap replicates.",
    ),
    TableSpec(
        "Benchmark_deltas",
        "Benchmark paired differences",
        "Deltas compare the combined score with each prespecified baseline on the same benchmark rows.",
    ),
    TableSpec(
        "Benchmark_calibration",
        "Benchmark calibration metrics",
        "Brier score, expected calibration error and grouped cross-validation metrics assess probability calibration.",
    ),
    TableSpec(
        "Benchmark_transport",
        "Cross-cohort benchmark transport",
        "Models were fitted to one cohort's feature labels and tested on the other cohort's labels.",
    ),
    TableSpec(
        "Benchmark_heldout",
        "Held-out credible-set and target-gene validation",
        "Grouped folds prevent the same credible set or target from contributing to both training and testing.",
    ),
    TableSpec(
        "Benchmark_coverage",
        "Benchmark label and feature coverage",
        "This two-column table summarizes the rows available for each benchmark task.",
    ),
    TableSpec(
        "Benchmark_priority_sensitivity",
        "Original-priority versus benchmark-score sensitivity",
        "Both scores were evaluated on the same rows to test construction dependence.",
    ),
    TableSpec(
        "Strict_eqtl_sensitivity",
        "Strict exact-variant eQTL sensitivity",
        "The strict label requires an exact-variant nominal association in the prespecified myeloid datasets.",
    ),
    TableSpec(
        "Top10_eQTL_triangulation",
        "Top 10 exact-variant eQTL triangulation",
        "Matched and mismatched target pairs were queried across eight myeloid datasets.",
        columns=(
            "target",
            "variant_id",
            "pair_type",
            "query_group",
            "datasets_queried",
            "nominal_datasets",
            "minimum_pvalue",
            "minimum_qvalue",
            "replicated_nominal",
            "any_nominal",
        ),
    ),
    TableSpec(
        "Top10_coloc",
        "Top 10 myeloid eQTL colocalization",
        "H4 is the posterior probability of a shared causal variant.",
        columns=(
            "target",
            "variant",
            "rsid",
            "eqtl_study",
            "eqtl_condition",
            "variants_harmonized",
            "eqtl_pvalue",
            "coloc_status",
            "H4",
            "coloc_supported",
        ),
    ),
    TableSpec(
        "Top10_coloc_sensitivity",
        "Top 10 colocalization prior sensitivity",
        "Support requires H4 >= 0.80 across the reported prior settings.",
    ),
    TableSpec(
        "Clinical_translation",
        "Indication-aware translational annotation",
        "Clinical-stage annotations describe hypotheses and do not establish treatment recommendations.",
        columns=(
            "priority_rank",
            "target",
            "variant",
            "external_cluster_count",
            "experimental_support",
            "exact_variant_eqtl",
            "formal_colocalization",
            "clinical_stage",
            "ibd_indication_drugs",
            "risk_effect_on_target",
            "translational_priority",
            "clinical_treatment_recommendation_supported",
            "clinical_recommendation_boundary",
        ),
    ),
    TableSpec(
        "PRKCB_drug_selectivity",
        "PRKCB drug-family selectivity annotation",
        "PRKCB-specific means the annotated drug targets PRKCB without the broader PKC family.",
    ),
    TableSpec(
        "PRKCB_IBD_trial_search",
        "PRKCB and PKC inhibitor IBD trial search",
        "Records were identified from ClinicalTrials.gov and the public trial report.",
    ),
    TableSpec(
        "PRKCB_IBD_trial_results",
        "Sotrastaurin ulcerative colitis trial results",
        "The credibility interval for clinical-remission difference crosses zero.",
    ),
    TableSpec(
        "PRKCB_CEDAR2_eQTL",
        "CEDAR2 PRKCB gut eQTL nodes",
        "Macrophage-like nodes are annotated explicitly; exact rs7404003 status is retained.",
        columns=(
            "dataset",
            "node_or_leaf",
            "annotated_cell_type",
            "cell_abbreviation",
            "gene",
            "top_variant",
            "beta",
            "nominal_pvalue",
            "window_adjusted_pvalue",
            "qvalue",
            "macrophage_like",
            "exact_rs7404003",
        ),
    ),
    TableSpec(
        "PRKCB_MAP3K8_direction",
        "PRKCB and MAP3K8 risk-allele direction",
        "The expected action follows the sign of the risk-allele effect on target expression.",
    ),
    TableSpec(
        "PRKCB_MR_estimates",
        "PRKCB two-sample Mendelian randomization estimates",
        "Estimates are reported per one-standard-deviation change in PRKCB expression.",
    ),
    TableSpec(
        "PRKCB_MR_instruments",
        "PRKCB Mendelian randomization instruments",
        "The F statistic and ratio summarize instrument strength and the Wald effect estimate.",
        columns=(
            "rsid",
            "gwas_beta",
            "gwas_se",
            "eql_beta",
            "eql_se",
            "gwas_p",
            "eql_p",
            "f_statistic",
            "ratio",
            "ratio_se",
            "r2_with_lead",
        ),
    ),
    TableSpec(
        "PRKCB_MR_Steiger",
        "PRKCB Steiger direction checks",
        "Correct direction indicates that the exposure explains more variance than the outcome.",
    ),
    TableSpec(
        "PRKCB_MR_leave_one_out",
        "PRKCB leave-one-instrument-out sensitivity",
        "Each row omits one approximately independent instrument and refits the MR summary.",
    ),
    TableSpec(
        "PRKCB_SMR_HEIDI",
        "PRKCB SMR and HEIDI-style sensitivity",
        "These are summary-statistic implementations and were not run with the official SMR executable.",
        columns=(
            "method",
            "variant",
            "eqtl_beta",
            "eqtl_se",
            "gwas_beta",
            "gwas_se",
            "smr_beta",
            "smr_se",
            "smr_p",
            "q",
            "p",
            "heterogeneity_status",
        ),
    ),
    TableSpec(
        "PRKCB_HEIDI_windows",
        "PRKCB HEIDI-style local-window sensitivity",
        "Each window limits the heterogeneity analysis to variants near the lead eQTL.",
    ),
    TableSpec(
        "PRKCB_SuSiE_credible_sets",
        "PRKCB SuSiE credible sets",
        "Fine-mapping used the external 1000 Genomes European LD reference.",
        columns=(
            "L",
            "trait",
            "credible_set",
            "n_snps",
            "top_snp",
            "top_pip",
            "credible_set_coverage",
        ),
    ),
    TableSpec(
        "PRKCB_SuSiE_coloc",
        "PRKCB SuSiE colocalization pair summary",
        "Only the eight credible-set pairs with the largest H4 values are shown; the full matrix remains in the workbook.",
        columns=(
            "nsnps",
            "hit1",
            "hit2",
            "PP.H0.abf",
            "PP.H1.abf",
            "PP.H2.abf",
            "PP.H3.abf",
            "PP.H4.abf",
        ),
        sort_by="PP.H4.abf",
        descending=True,
        max_rows=8,
    ),
    TableSpec(
        "PRKCB_QTL_coverage",
        "PRKCB molecular QTL coverage",
        "Absent sQTL, pQTL or caQTL records are classified as resource non-coverage.",
    ),
    TableSpec(
        "PRKCB_signature_reversal",
        "PRKCB public signature-reversal analysis",
        "Correlations and permutation tests compare PRKCB knockdown and disease signatures.",
    ),
    TableSpec(
        "PRKCB_virtual_pathways",
        "PRKCB virtual pathway effects",
        "Values are computational predictions from a shrunken co-expression network.",
    ),
    TableSpec(
        "PRKCB_virtual_bootstrap1000",
        "PRKCB virtual-knockdown donor bootstrap",
        "Intervals and positive fractions use 1,000 donor-block bootstrap replicates.",
    ),
    TableSpec(
        "PRKCB_virtual_shrinkage",
        "PRKCB virtual-network shrinkage sensitivity",
        "Each row varies the network shrinkage parameter and recomputes pathway effects.",
    ),
    TableSpec(
        "PRKCB_followup_summary",
        "PRKCB clinical and direction follow-up metrics",
        "The summary records drug selectivity, trial availability and risk-direction boundaries.",
    ),
]


MAIN_SHEET_DISPOSITION = {
    "Variant_sets": "Main Table 2",
    "Positive_control_atlas": "Main Table 4 support",
    "Active_IBD_summary": "Main text summary",
    "Candidate_evidence": "Main Table 11 support",
    "Experimental_evidence": "Main Table 8 support",
    "Experimental_summary": "Main text summary",
    "Claim_calibration": "Main Table 11",
    "Negative_control": "Main Table 5 support",
    "Threshold_sensitivity": "Main Table 6",
    "Score_overlap": "Main text summary",
}


FULL_DATA_DISPOSITION = {
    "Top20": "Full data only",
    "All_candidates": "Full data only",
    "Active_IBD_leads": "Full data only",
    "Pathway_corrected_ORA": "Full data only",
    "Pathway_all_L2G_ranks": "Full data only",
    "Pathway_prerank": "Full data only",
    "Cross_study_L2G_counts": "Full data only",
    "External_gene_rankings": "Full data only",
    "PRKCQ_repeated_scores": "Full data only",
    "Macrophage_CRISPR_genes": "Full data only",
    "Macrophage_CRISPR_guides": "Full data only",
    "PerturbSeq_UMI_thresholds": "Full data only",
    "PerturbSeq_guides": "Full data only",
    "PerturbSeq_leave_guide": "Full data only",
    "PerturbSeq_leave_donor": "Full data only",
    "GSE175685_GSEA": "Full data only",
    "Atlas_top_genes": "Full data only",
    "eQTL_variant_benchmark": "Full data only",
    "PRKCB_virtual_bootstrap": "Full data only",
}


HEADER_LABELS = {
    "p": "P",
    "p_value": "P",
    "pvalue": "P",
    "pValue": "P",
    "fdr": "FDR",
    "fdr_bh": "FDR",
    "qvalue": "q value",
    "auroc": "AUROC",
    "auprc": "AUPRC",
    "ci_low": "Lower 95% CI",
    "ci_high": "Upper 95% CI",
    "delta_auroc": "Delta AUROC",
    "delta_auprc": "Delta AUPRC",
    "prkcb": "PRKCB",
    "rsid": "rsID",
    "p_value": "P",
    "log2fc": "log2FC",
    "eqtl": "eQTL",
}


def clean_header(value: object, index: int) -> str:
    if value is None or str(value).strip() == "":
        return f"Column {index}"
    text = str(value).strip()
    if text in HEADER_LABELS:
        return HEADER_LABELS[text]
    replacements = [
        ("mean_log2fc", "Mean log2FC"),
        ("median_log2fc", "Median log2FC"),
        ("bootstrap_ci_low", "Bootstrap lower 95% CI"),
        ("bootstrap_ci_high", "Bootstrap upper 95% CI"),
        ("ci_low", "Lower 95% CI"),
        ("ci_high", "Upper 95% CI"),
        ("p_value", "P"),
        ("pvalue", "P"),
        ("p_val", "P"),
        ("fdr", "FDR"),
        ("auroc", "AUROC"),
        ("auprc", "AUPRC"),
        ("eqtl", "eQTL"),
        ("pqtl", "pQTL"),
        ("sqtl", "sQTL"),
        ("caqtl", "caQTL"),
        ("prkcb", "PRKCB"),
        ("prkcq", "PRKCQ"),
        ("map3k8", "MAP3K8"),
        ("rsid", "rsID"),
        ("nct_id", "NCT ID"),
        ("log2fc", "log2FC"),
        ("log2", "log2"),
        ("celltype", "Cell type"),
        ("nsnps", "N SNPs"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if text and text[0].islower():
        text = text[0].upper() + text[1:]
    return text


def format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, (datetime, date)):
        return value.isoformat()[:10]
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        absolute = abs(value)
        if absolute == 0:
            return "0.000"
        if absolute < 0.001 or absolute >= 10000:
            return f"{value:.3e}"
        if absolute < 0.01:
            return f"{value:.5f}"
        if absolute < 1:
            return f"{value:.4f}"
        if absolute < 100:
            return f"{value:.3f}"
        return f"{value:.1f}"
    text = str(value).strip()
    if text.lower() == "true":
        return "Yes"
    if text.lower() == "false":
        return "No"
    return re.sub(r"\s+", " ", text)


def read_sheet(workbook, source: str) -> tuple[list[str], list[list[object]]]:
    worksheet = workbook[source]
    raw_rows = list(worksheet.iter_rows(values_only=True))
    if not raw_rows:
        return [], []
    header_count = len(raw_rows[0])
    headers = [
        str(value).strip() if value is not None and str(value).strip() else f"Column {index + 1}"
        for index, value in enumerate(raw_rows[0])
    ]
    rows = [list(row[:header_count]) for row in raw_rows[1:] if any(value is not None for value in row)]
    while headers and headers[-1].startswith("Column "):
        column_index = len(headers) - 1
        if any(
            row[column_index] is not None
            for row in rows
            if column_index < len(row)
        ):
            break
        headers.pop()
    rows = [row[: len(headers)] for row in rows]
    return headers, rows


def select_rows(
    headers: list[str],
    rows: list[list[object]],
    spec: TableSpec,
) -> tuple[list[str], list[list[object]], int]:
    selected_headers = headers
    if spec.columns:
        missing = [column for column in spec.columns if column not in headers]
        if missing:
            raise KeyError(f"{spec.source}: missing columns {missing}")
        indices = [headers.index(column) for column in spec.columns]
        selected_headers = [clean_header(column, index + 1) for index, column in enumerate(spec.columns)]
        selected_rows = [
            [row[index] if index < len(row) else None for index in indices]
            for row in rows
        ]
    else:
        selected_headers = [
            clean_header(header, index + 1)
            for index, header in enumerate(selected_headers)
        ]
        selected_rows = [row[:] for row in rows]

    if spec.sort_by:
        if spec.sort_by not in headers:
            raise KeyError(f"{spec.source}: missing sort column {spec.sort_by}")
        source_index = headers.index(spec.sort_by)

        def sort_value(row: list[object]) -> float:
            value = row[source_index]
            return float(value) if isinstance(value, (int, float)) else float("-inf")

        selected_rows.sort(key=sort_value, reverse=spec.descending)

    full_count = len(selected_rows)
    if spec.max_rows is not None:
        selected_rows = selected_rows[: spec.max_rows]
    return selected_headers, selected_rows, full_count


def set_cell_margins(table, top: int = 30, start: int = 45, bottom: int = 30, end: int = 45) -> None:
    tbl_pr = table._tbl.tblPr
    margins = tbl_pr.first_child_found_in("w:tblCellMar")
    if margins is None:
        margins = OxmlElement("w:tblCellMar")
        tbl_pr.append(margins)
    for edge, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
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


def set_fixed_layout(table) -> None:
    tbl_pr = table._tbl.tblPr
    layout = tbl_pr.first_child_found_in("w:tblLayout")
    if layout is None:
        layout = OxmlElement("w:tblLayout")
        tbl_pr.append(layout)
    layout.set(qn("w:type"), "fixed")


def calculate_widths(headers: list[str], rows: list[list[object]], total_width: float) -> list[float]:
    scores = []
    for index, header in enumerate(headers):
        values = [format_value(row[index]) for row in rows[:60] if index < len(row)]
        max_value = max((min(len(value), 36) for value in values), default=4)
        score = max(5.5, min(36, max(len(header) * 0.72, max_value)))
        scores.append(score)
    score_sum = sum(scores)
    widths = [total_width * score / score_sum for score in scores]
    minimum = min(0.43, total_width / max(len(headers), 1))
    widths = [max(minimum, width) for width in widths]
    width_sum = sum(widths)
    return [width * total_width / width_sum for width in widths]


def add_table(
    document: Document,
    headers: list[str],
    rows: list[list[object]],
    *,
    font_size: float,
    total_width: float,
) -> None:
    table = document.add_table(rows=len(rows) + 1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    set_fixed_layout(table)
    set_three_line_borders(table)
    set_cell_margins(table)
    widths = calculate_widths(headers, rows, total_width)

    for column_index, width in enumerate(widths):
        table.columns[column_index].width = Inches(width)

    for row_index, values in enumerate([headers, *rows]):
        row = table.rows[row_index]
        row.height = None
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
            run = paragraph.add_run(format_value(value))
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
    section.top_margin = Inches(0.45)
    section.bottom_margin = Inches(0.45)
    section.left_margin = Inches(0.48)
    section.right_margin = Inches(0.48)

    styles = document.styles
    styles["Normal"].font.name = "Times New Roman"
    styles["Normal"].font.size = Pt(9)
    styles["Normal"].paragraph_format.space_after = Pt(0)
    styles["Normal"].paragraph_format.line_spacing = 1.0
    for style_name in ["Title", "Heading 1", "Heading 2"]:
        styles[style_name].font.name = "Times New Roman"

    footer = section.footer.paragraphs[0]
    footer.text = "Supplementary Tables | Page "
    footer.runs[0].font.name = "Times New Roman"
    footer.runs[0].font.size = Pt(8)
    add_page_number(footer)


def add_title_page(document: Document) -> None:
    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.paragraph_format.space_before = Pt(50)
    title.paragraph_format.space_after = Pt(18)
    run = title.add_run("Supplementary Tables")
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(24)

    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.paragraph_format.space_after = Pt(20)
    run = subtitle.add_run(
        "An evidence ladder for non-coding inflammatory bowel disease variants "
        "prioritizes intestinal macrophage regulatory mechanisms"
    )
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(13)

    note = document.add_paragraph()
    note.alignment = WD_ALIGN_PARAGRAPH.CENTER
    note.paragraph_format.left_indent = Inches(1.2)
    note.paragraph_format.right_indent = Inches(1.2)
    run = note.add_run(
        "This document contains review-ready supplementary summaries. "
        "The full row-level source tables and matrices are provided in "
        "IBD_AlphaGenome2_supplementary_tables.xlsx."
    )
    run.font.name = "Times New Roman"
    run.font.size = Pt(10)

    version = document.add_paragraph()
    version.alignment = WD_ALIGN_PARAGRAPH.CENTER
    version.paragraph_format.space_before = Pt(18)
    run = version.add_run("Version: 14 November 2026")
    run.font.name = "Times New Roman"
    run.font.size = Pt(9)


def add_contents(document: Document, index_rows: list[dict[str, object]]) -> None:
    document.add_page_break()
    heading = document.add_paragraph()
    heading.paragraph_format.keep_with_next = True
    run = heading.add_run("Contents of Supplementary Tables")
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(14)
    document.add_paragraph(
        "Main Tables 1-11 remain in the manuscript. Supplementary Tables S1-S52 "
        "contain additional analyses and are not duplicates of the main display tables."
    )
    headers = ["Table", "Title", "Source worksheet", "Rows"]
    rows = [
        [
            row["supplementary_table"],
            row["title"],
            row["source_sheet"],
            row["rows_in_word"],
        ]
        for row in index_rows
    ]
    add_table(document, headers, rows, font_size=7.8, total_width=10.55)


def add_supplement_tables(document: Document, tables: list[dict[str, object]]) -> None:
    for table_data in tables:
        document.add_page_break()
        caption = document.add_paragraph()
        caption.paragraph_format.keep_with_next = True
        caption.paragraph_format.space_after = Pt(5)
        number = table_data["supplementary_table"]
        run = caption.add_run(f"Supplementary Table {number}. ")
        run.bold = True
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)
        run = caption.add_run(str(table_data["title"]))
        run.bold = True
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)

        headers = table_data["headers"]
        rows = table_data["rows"]
        font_size = 7.4 if len(headers) >= 12 else 8.0
        add_table(
            document,
            headers,
            rows,
            font_size=font_size,
            total_width=10.55,
        )

        note = document.add_paragraph()
        note.paragraph_format.space_before = Pt(4)
        note.paragraph_format.space_after = Pt(0)
        note.paragraph_format.line_spacing = 1.0
        source_note = (
            f'Source worksheet: {table_data["source_sheet"]}. '
            f'{table_data["note"]}'
        )
        if table_data["status"] == "summarized":
            source_note += (
                f' This Word table shows {table_data["rows_in_word"]} of '
                f'{table_data["full_rows_in_workbook"]} source rows.'
            )
        run = note.add_run(source_note)
        run.font.name = "Times New Roman"
        run.font.size = Pt(8)


def add_workbook_index(document: Document, rows: list[dict[str, object]]) -> None:
    document.add_page_break()
    heading = document.add_paragraph()
    heading.paragraph_format.keep_with_next = True
    run = heading.add_run("Supplementary Workbook Index")
    run.bold = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(14)
    document.add_paragraph(
        "Every worksheet in the full supplementary workbook is listed below. "
        "Large row-level tables remain in the workbook to preserve complete source data."
    )
    headers = ["Worksheet", "Source rows", "Source columns", "Disposition"]
    table_rows = [
        [
            row["source_sheet"],
            row["full_rows_in_workbook"],
            row["source_columns"],
            row["disposition"],
        ]
        for row in rows
    ]
    add_table(document, headers, table_rows, font_size=7.4, total_width=10.55)


def infer_disposition(source: str, included: dict[str, dict[str, object]]) -> str:
    if source in included:
        return f'{included[source]["supplementary_table"]} ({included[source]["status"]})'
    if source in MAIN_SHEET_DISPOSITION:
        return MAIN_SHEET_DISPOSITION[source]
    if source in FULL_DATA_DISPOSITION:
        return FULL_DATA_DISPOSITION[source]
    raise KeyError(f"No disposition assigned for worksheet {source}")


def build_index(
    workbook,
    included: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    index_rows = []
    for source in workbook.sheetnames:
        worksheet = workbook[source]
        full_rows = max(worksheet.max_row - 1, 0)
        index_rows.append(
            {
                "source_sheet": source,
                "full_rows_in_workbook": full_rows,
                "source_columns": worksheet.max_column,
                "disposition": infer_disposition(source, included),
            }
        )
    return index_rows


def write_index_csv(
    table_rows: list[dict[str, object]],
    workbook_rows: list[dict[str, object]],
) -> None:
    INDEX_CSV.parent.mkdir(parents=True, exist_ok=True)
    workbook_by_source = {row["source_sheet"]: row for row in workbook_rows}
    with INDEX_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "supplementary_table",
                "title",
                "source_sheet",
                "status",
                "rows_in_word",
                "full_rows_in_workbook",
                "source_columns",
                "disposition",
            ],
        )
        writer.writeheader()
        for row in table_rows:
            workbook_row = workbook_by_source[row["source_sheet"]]
            writer.writerow(
                {
                    "supplementary_table": row["supplementary_table"],
                    "title": row["title"],
                    "source_sheet": row["source_sheet"],
                    "status": row["status"],
                    "rows_in_word": row["rows_in_word"],
                    "full_rows_in_workbook": row["full_rows_in_workbook"],
                    "source_columns": row["source_columns"],
                    "disposition": workbook_row["disposition"],
                }
            )
        included_sources = {row["source_sheet"] for row in table_rows}
        for source in workbook_by_source:
            if source in included_sources:
                continue
            row = workbook_by_source[source]
            writer.writerow(
                {
                    "supplementary_table": "",
                    "title": "",
                    "source_sheet": source,
                    "status": "not_expanded_in_word",
                    "rows_in_word": 0,
                    "full_rows_in_workbook": row["full_rows_in_workbook"],
                    "source_columns": row["source_columns"],
                    "disposition": row["disposition"],
                }
            )


def audit_document(
    document: Document,
    table_rows: list[dict[str, object]],
    workbook_rows: list[dict[str, object]],
) -> dict[str, object]:
    paragraphs = [paragraph.text.strip() for paragraph in document.paragraphs]
    captions = [
        text
        for text in paragraphs
        if text.startswith("Supplementary Table S")
    ]
    caption_numbers = []
    caption_titles = []
    for caption in captions:
        match = re.match(r"Supplementary Table S(\d+)\.\s+(.+)", caption)
        if match:
            caption_numbers.append(int(match.group(1)))
            caption_titles.append(match.group(2))
    expected_numbers = list(range(1, len(table_rows) + 1))
    workbook_sources = {row["source_sheet"] for row in workbook_rows}
    indexed_dispositions = [
        row["disposition"] for row in workbook_rows if row["disposition"]
    ]
    checks = {
        "document_exists": OUTPUT.exists(),
        "table_count_matches": len(document.tables) == len(table_rows) + 2,
        "supplement_numbers_sequential": caption_numbers == expected_numbers,
        "all_source_sheets_indexed": len(workbook_rows) == 81 and bool(workbook_sources),
        "all_dispositions_assigned": all(indexed_dispositions),
        "titles_unique": len(set(caption_titles)) == len(caption_titles),
        "no_empty_tables": all(
            len(table.rows) > 1 and len(table.columns) > 0
            for table in document.tables
        ),
        "no_main_table_captions_duplicated": not any(
            re.match(r"^Table \d+\.", text)
            for text in paragraphs
        ),
        "all_tables_have_notes": all(
            f'Source worksheet: {row["source_sheet"]}.'
            in "\n".join(paragraphs)
            for row in table_rows
        ),
        "all_tables_have_rows": all(int(row["rows_in_word"]) > 0 for row in table_rows),
    }
    return {
        "generated": datetime.now().astimezone().isoformat(),
        "output": str(OUTPUT),
        "supplementary_tables": len(table_rows),
        "word_tables_including_contents_and_index": len(document.tables),
        "workbook_sheets_indexed": len(workbook_rows),
        "included_sources": sum(row["status"] == "included" for row in table_rows),
        "summarized_sources": sum(row["status"] == "summarized" for row in table_rows),
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }


def main() -> None:
    if not WORKBOOK.exists():
        raise FileNotFoundError(WORKBOOK)
    workbook = load_workbook(WORKBOOK, read_only=True, data_only=True)
    included: dict[str, dict[str, object]] = {}
    table_rows: list[dict[str, object]] = []

    for number, spec in enumerate(SPECS, start=1):
        if spec.source not in workbook.sheetnames:
            raise KeyError(f"Missing worksheet: {spec.source}")
        headers, selected_rows, full_count = select_rows(
            *read_sheet(workbook, spec.source),
            spec,
        )
        status = "included" if len(selected_rows) == full_count else "summarized"
        row = {
            "supplementary_table": f"S{number}",
            "title": spec.title,
            "note": spec.note,
            "source_sheet": spec.source,
            "status": status,
            "headers": headers,
            "rows": selected_rows,
            "rows_in_word": len(selected_rows),
            "full_rows_in_workbook": full_count,
            "source_columns": workbook[spec.source].max_column,
        }
        table_rows.append(row)
        included[spec.source] = row

    workbook_rows = build_index(workbook, included)
    document = Document()
    style_document(document)
    add_title_page(document)
    add_contents(document, table_rows)
    add_supplement_tables(document, table_rows)
    add_workbook_index(document, workbook_rows)
    document.save(OUTPUT)
    MANUSCRIPT_COPY.parent.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_COPY.write_bytes(OUTPUT.read_bytes())

    write_index_csv(table_rows, workbook_rows)
    audit = audit_document(document, table_rows, workbook_rows)
    AUDIT_JSON.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("Supplementary table audit failed.")
    print(f"Wrote {OUTPUT}")
    print(f"Wrote {MANUSCRIPT_COPY}")
    print(f"Wrote {INDEX_CSV}")


if __name__ == "__main__":
    main()
