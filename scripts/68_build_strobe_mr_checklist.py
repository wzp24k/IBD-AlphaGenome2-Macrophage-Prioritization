from __future__ import annotations

import json
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


TEMPLATE = ROOT / "work" / "STROBE-MR-checklist-fillable-official.docx"
TEMPLATE_URL = "https://www.strobe-mr.org/download/strobe-mr-checklist-fillable-doc/"
OUTPUT = ROOT / "STROBE-MR_Checklist_Filled.docx"
MANUSCRIPT_COPY = ROOT / "manuscript" / "STROBE-MR_Checklist_Filled.docx"
AUDIT = ROOT / "results" / "strobe_mr_checklist_audit.json"


ROWS = {
    1: (
        "p. 2",
        "Addressed. MR is identified in the Abstract as a secondary PRKCB "
        "mechanism-triangulation analysis. The title does not use MR as the "
        "main study design because the primary contribution is a non-coding "
        "variant evidence ladder.",
    ),
    3: (
        "p. 3",
        "Addressed. Exposure: genetically proxied monocyte/macrophage PRKCB "
        "mRNA expression. Outcome: inflammatory bowel disease risk. PRKCB was "
        "prioritized from an IBD credible set and tested because the myeloid "
        "eQTL direction, PKC signaling and available pharmacological data "
        "made a directional expression-to-disease relation plausible. MR was "
        "used as one triangulation layer because observational evidence alone "
        "cannot distinguish direction.",
    ),
    4: (
        "p. 3",
        "Addressed. Objective: estimate the directional relation between "
        "genetically proxied PRKCB expression and IBD risk. Prespecified "
        "working hypothesis: higher myeloid PRKCB expression is associated "
        "with higher IBD risk. The manuscript states that MR estimates causal "
        "effects only under relevance, independence and exclusion-restriction "
        "assumptions.",
    ),
    6: (
        "p. 8",
        "Addressed. This was a two-sample summary-statistic MR analysis "
        "embedded as a secondary PRKCB mechanism analysis within a larger "
        "computational study.",
    ),
    7: (
        "p. 8",
        "Addressed. Exposure data were monocyte eQTL summary statistics "
        "(n = 191). Outcome data were the combined IBD GWAS summary statistics "
        "(n = 398,668). LD was estimated in 633 European 1000 Genomes "
        "high-coverage samples. All data were public and de-identified.",
    ),
    8: (
        "p. 8",
        "Addressed. The PRKCB region contained 877 variants with matched eQTL, "
        "GWAS and LD estimates. LD clumping at r2 < 0.01 retained three "
        "approximately independent instruments. No separate power calculation "
        "was performed; instrument number was limited by the available "
        "independent regional signals.",
    ),
    9: (
        "p. 8",
        "Addressed. Variants were ALT-harmonized across eQTL and GWAS data. "
        "LD clumping used r2 thresholds of 0.001, 0.01 and 0.10. Instrument "
        "strength was summarized by F statistics. The primary threshold was "
        "r2 < 0.01.",
    ),
    10: (
        "pp. 8 and 29",
        "Addressed. Exposure was PRKCB mRNA expression in monocytes. Outcome "
        "was IBD case-control status in the GWAS. Effects were analyzed on the "
        "log-odds scale per one-standard-deviation difference in genetically "
        "proxied expression.",
    ),
    11: (
        "p. 35",
        "Addressed. This study used publicly available de-identified summary "
        "statistics. No new participants were recruited and no new ethics "
        "approval or informed consent was required.",
    ),
    12: (
        "pp. 8 and 29",
        "Addressed. Relevance was supported by strong instrument F statistics. "
        "Independence was addressed by LD clumping and European ancestry "
        "matching. Exclusion restriction cannot be proven directly and was "
        "examined using MR-Egger intercept, HEIDI-style heterogeneity, "
        "leave-one-instrument-out analyses and SuSiE-colocalization.",
    ),
    13: (
        "p. 8",
        "Addressed. The main analysis used inverse-variance weighting, "
        "weighted median and MR-Egger. A Steiger-style check assessed "
        "direction.",
    ),
    14: (
        "p. 8",
        "Addressed. The exposure unit was one standard deviation of "
        "genetically proxied PRKCB expression. The outcome was log odds of "
        "IBD. Estimates are reported per expression standard deviation.",
    ),
    15: (
        "p. 8",
        "Addressed. Variant weights were the eQTL effect estimates. Variants "
        "were pruned using external European LD and retained when r2 was below "
        "the prespecified clumping threshold.",
    ),
    16: (
        "p. 8",
        "Addressed. Two-sample summary-data MR used inverse-variance "
        "weighting, weighted median and MR-Egger. No additional covariates "
        "were available or adjusted for in the summary-statistic analysis.",
    ),
    17: (
        "p. 8",
        "Addressed. Variants lacking complete eQTL, GWAS or LD information "
        "were excluded during harmonization. No outcome or exposure values "
        "were imputed.",
    ),
    18: (
        "p. 29",
        "No multiplicity adjustment was applied across the exploratory MR "
        "estimators. MR was a secondary triangulation component and the "
        "manuscript does not treat these estimates as a primary multiple-testing "
        "claim.",
    ),
    19: (
        "pp. 8, 29 and 30",
        "Addressed. Instrument F statistics and Steiger direction checks were "
        "reported. Regional heterogeneity was assessed by HEIDI-style tests, "
        "MR-Egger intercept and leave-one-out analyses. SuSiE-colocalization "
        "was used as an orthogonal assumption check.",
    ),
    20: (
        "pp. 8, 29 and 30",
        "Addressed. Sensitivity analyses varied LD clumping thresholds, "
        "reported weighted median and MR-Egger estimates, performed "
        "leave-one-instrument-out analyses, examined HEIDI-style windows and "
        "compared MR with SMR and SuSiE-colocalization.",
    ),
    21: (
        "p. 8",
        "Addressed. See items 9a and 9b.",
    ),
    22: (
        "p. 8",
        "Addressed. Analyses used Python 3.12.10 with NumPy 2.5.4, pandas "
        "3.1.1, SciPy 1.18.1 and statsmodels 0.15.0. SuSiE-colocalization was "
        "implemented in R. The analysis scripts are publicly available.",
    ),
    23: (
        "p. 8",
        "No prospective protocol registration was performed. The MR component "
        "was defined as a secondary analysis in the analysis code before "
        "manuscript preparation.",
    ),
    25: (
        "p. 29",
        "Addressed. See items 10a-10d.",
    ),
    26: (
        "p. 29",
        "Addressed. The region contributed 877 harmonized variants and three "
        "approximately independent instruments. The eQTL sample contained 191 "
        "participants and the IBD GWAS contained 398,668 participants. Reasons "
        "for variant exclusion were incomplete harmonization or LD pruning.",
    ),
    27: (
        "p. 29",
        "Addressed. Mean instrument F was 419. The three instrument-level eQTL "
        "effects were 1.138, 1.138 and 0.799. Full variant effects and allele "
        "frequencies are provided in Supplementary Table S3.",
    ),
    28: (
        "p. 29",
        "The outcome GWAS is a meta-analysis. Between-study heterogeneity was "
        "not re-estimated within this secondary summary-statistic MR "
        "component.",
    ),
    29: (
        "pp. 8 and 29",
        "i. Exposure and outcome data were restricted to European-ancestry "
        "analyses to improve variant-exposure comparability. ii. The monocyte "
        "eQTL and IBD GWAS are different public datasets. No known participant "
        "overlap exists, although individual-level overlap cannot be verified "
        "from summary statistics.",
    ),
    30: (
        "p. 29",
        "Addressed. Main MR results and uncertainty are reported in items "
        "11a-11d and Figure 12.",
    ),
    31: (
        "p. 29",
        "Addressed. Instrument-level eQTL and GWAS effect estimates are "
        "reported in Supplementary Table S3. The lead instrument rs7404003 had "
        "an eQTL beta of 1.138 and an IBD GWAS beta of 0.0585.",
    ),
    32: (
        "p. 29",
        "Addressed. IVW beta was 0.0517 log odds per SD (SE 0.00595; derived "
        "95% CI 0.0401-0.0634; OR 1.053). Weighted-median beta was 0.0514 "
        "(SE 0.00728; derived 95% CI 0.0371-0.0657). MR-Egger beta was -0.178 "
        "(SE 0.00891; derived 95% CI -0.1951 to -0.1602), with an intercept of "
        "0.260 (P = 0.025).",
    ),
    33: (
        "Not applicable",
        "Absolute-risk translation was not performed because no population "
        "absolute-risk model or time horizon was prespecified.",
    ),
    34: (
        "p. 30",
        "Addressed. Figure 12A shows the MR forest plot and Figure 12B shows "
        "regional eQTL versus IBD GWAS effects with the SMR slope.",
    ),
    35: (
        "p. 29",
        "Addressed. Assumption diagnostics are reported in items 12a-12b.",
    ),
    36: (
        "p. 29",
        "Addressed. Instrument relevance was strong. The MR-Egger intercept "
        "and full-region HEIDI-style test indicated directional pleiotropy or "
        "heterogeneity, so the exclusion restriction was not considered fully "
        "satisfied.",
    ),
    37: (
        "p. 29",
        "Addressed. Full-region HEIDI-style heterogeneity P was <0.001. The "
        "three-instrument pruned subset had Q = 1.750 and P = 0.417. "
        "MR-Egger intercept was 0.260 with P = 0.025.",
    ),
    38: (
        "pp. 29-30",
        "Addressed. Sensitivity and additional analyses are reported in items "
        "13a-13e.",
    ),
    39: (
        "p. 29",
        "Addressed. IVW and weighted-median estimates remained positive after "
        "omitting each instrument. Results were stable across r2 thresholds "
        "0.001, 0.01 and 0.10.",
    ),
    40: (
        "pp. 29 and 30",
        "Addressed. MR-Egger changed direction. The SMR-style beta was 0.0514 "
        "(P < 0.001). HEIDI-style heterogeneity persisted in local windows. "
        "SuSiE-colocalization reached H4 = 1.0 for matched credible sets.",
    ),
    41: (
        "p. 29",
        "Addressed. All three instruments passed the Steiger-style direction "
        "check. Bidirectional MR was not performed.",
    ),
    42: (
        "pp. 29 and 30",
        "Addressed. GTEx sigmoid-colon eQTL effects were negative (NES -0.171), "
        "whereas myeloid eQTL effects were positive. CEDAR2 did not identify a "
        "macrophage-like PRKCB eQTL. Public overexpression-signature reversal "
        "was weak and cohort-dependent.",
    ),
    43: (
        "p. 29",
        "Addressed. Leave-one-instrument-out estimates are reported in "
        "Supplementary Table S3. Figure 12 also shows forest and regional "
        "effect plots.",
    ),
    45: (
        "p. 29",
        "Addressed. The risk allele increased monocyte PRKCB expression. IVW "
        "and weighted median were positive, SMR was significant and "
        "SuSiE-colocalization reached H4 = 1.0. MR-Egger was discordant and "
        "regional heterogeneity was present.",
    ),
    46: (
        "pp. 33-34",
        "Addressed. Limitations include only three instruments, directional "
        "inconsistency in MR-Egger, regional heterogeneity, context-dependent "
        "colon versus myeloid direction, unverifiable sample overlap, "
        "gene-environment equivalence and the absence of experimental "
        "causality.",
    ),
    47: (
        "p. 33",
        "Addressed. See items 16a-16c.",
    ),
    48: (
        "pp. 32-33",
        "Addressed. PRKCB is interpreted as a mechanistic-triangulation "
        "candidate rather than a proven causal IBD gene. The manuscript "
        "explicitly limits causal language.",
    ),
    49: (
        "p. 33",
        "Addressed. PRKCB encodes a PKC family kinase with plausible myeloid "
        "inflammatory functions. The risk allele increases monocyte expression, "
        "consistent with an inhibitory hypothesis in myeloid cells. The "
        "opposite colon eQTL direction prevents a tissue-universal mechanism "
        "claim, and gene-environment equivalence is not assumed.",
    ),
    50: (
        "p. 33",
        "Addressed. PRKCB is a translational hypothesis because pan-PKC "
        "inhibitors are not PRKCB-selective and the sotrastaurin UC trial was "
        "small, early-terminated and inconclusive. No treatment recommendation "
        "is made.",
    ),
    51: (
        "p. 33",
        "Addressed. Findings apply primarily to European-ancestry monocyte or "
        "myeloid expression and a single locus. Generalizability to other "
        "ancestries, tissues, disease states, exposure timings and exposure "
        "levels is not established.",
    ),
    53: (
        "p. 35",
        "Addressed. The authors received no specific funding for this study.",
    ),
    54: (
        "p. 34",
        "Addressed. Processed variant sets, LD matrices, MR/SMR/SuSiE summaries "
        "and analysis code are included in the submission package or public "
        "repository. External source datasets are identified in the manuscript "
        "and are available under their original access terms.",
    ),
    55: (
        "p. 35",
        "Addressed. The authors declare no competing interests.",
    ),
}


def set_cell_text(cell, text: str, bold: bool = False) -> None:
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.paragraph_format.line_spacing = 1.0
    run = paragraph.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = Pt(7)
    run.bold = bold
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP


def main() -> None:
    if not TEMPLATE.exists():
        TEMPLATE.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(TEMPLATE_URL, TEMPLATE)
    document = Document(TEMPLATE)
    table = document.tables[0]

    title_paragraph = document.paragraphs[0]
    metadata = [
        (
            "Filled STROBE-MR checklist",
            True,
            12,
        ),
        (
            "Manuscript: An evidence ladder for non-coding inflammatory bowel "
            "disease variants prioritizes intestinal macrophage regulatory mechanisms",
            False,
            9,
        ),
        (
            "MR component covered: two-sample summary-data MR of genetically "
            "proxied monocyte PRKCB expression and IBD risk, analyzed as a "
            "secondary mechanism-triangulation module.",
            False,
            9,
        ),
        (
            "Page references are based on the current 39-page 1114manuscript.docx "
            "and use exact Word page numbers.",
            False,
            9,
        ),
    ]
    for text, bold, size in reversed(metadata):
        paragraph = title_paragraph.insert_paragraph_before()
        paragraph.paragraph_format.space_after = Pt(3)
        run = paragraph.add_run(text)
        run.font.name = "Times New Roman"
        run.font.size = Pt(size)
        run.bold = bold

    page_column = table.columns[3]
    text_column = table.columns[4]
    page_column.width = Inches(0.58)
    text_column.width = Inches(3.15)

    for row_index, (page, text) in ROWS.items():
        row = table.rows[row_index]
        set_cell_text(row.cells[3], page)
        set_cell_text(row.cells[4], text)

    for row_index, row in enumerate(table.rows):
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.line_spacing = 1.0
                for run in paragraph.runs:
                    run.font.name = "Times New Roman"
                    if row_index not in ROWS:
                        run.font.size = Pt(7.5)

    document.save(OUTPUT)
    MANUSCRIPT_COPY.parent.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_COPY.write_bytes(OUTPUT.read_bytes())

    item_numbers = []
    blank_items = []
    for row_index in range(1, len(table.rows)):
        item = table.rows[row_index].cells[0].text.strip()
        if item.isdigit():
            item_numbers.append(int(item))
            if not table.rows[row_index].cells[3].text.strip() or not table.rows[
                row_index
            ].cells[4].text.strip():
                blank_items.append(int(item))

    checks = {
        "official_table_present": len(document.tables) == 1,
        "all_20_main_items_present": item_numbers == list(range(1, 21)),
        "all_main_items_completed": not blank_items,
        "supplementary_subitems_completed": all(
            table.rows[index].cells[3].text.strip()
            and table.rows[index].cells[4].text.strip()
            for index in ROWS
        ),
        "mr_role_identified_as_secondary": any(
            "secondary mechanism-triangulation" in paragraph.text
            for paragraph in document.paragraphs
        ),
        "page_references_present": all(
            table.rows[index].cells[3].text.strip() for index in ROWS
        ),
    }
    audit = {
        "generated": datetime.now().astimezone().isoformat(),
        "template": str(TEMPLATE),
        "output": str(OUTPUT),
        "main_items": len(item_numbers),
        "completed_rows": len(ROWS),
        "blank_main_items": blank_items,
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("STROBE-MR checklist audit failed.")
    print(f"Wrote {OUTPUT}")
    print(f"Wrote {MANUSCRIPT_COPY}")


if __name__ == "__main__":
    main()
