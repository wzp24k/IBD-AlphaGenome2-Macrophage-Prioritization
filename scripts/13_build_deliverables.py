from __future__ import annotations

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


def clean_markdown(text: str) -> str:
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"\1", text)
    text = re.sub(r"\*([^*]+)\*", r"\1", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    return text.strip()


def flatten_json(path: Path) -> pd.DataFrame:
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for key, value in data.items():
        if isinstance(value, dict):
            for nested_key, nested_value in value.items():
                rows.append((f"{key}.{nested_key}", nested_value))
        else:
            rows.append((key, value))
    return pd.DataFrame(rows, columns=["metric", "value"])


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


def build_docx(markdown_path: Path, output_path: Path, figure_dir: Path) -> None:
    text = markdown_path.read_text(encoding="utf-8")
    lines = text.splitlines()
    document = Document()
    section = document.sections[0]
    section.top_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.right_margin = Inches(1)
    styles = document.styles
    styles["Normal"].font.name = "Times New Roman"
    styles["Normal"].font.size = Pt(12)
    styles["Normal"].paragraph_format.line_spacing = 1.5
    for style_name in ["Title", "Heading 1", "Heading 2", "Heading 3"]:
        styles[style_name].font.name = "Times New Roman"
    add_page_number(section.footer.paragraphs[0])

    title_seen = False
    index = 0
    pending_tables = []
    while index < len(lines):
        line = lines[index].rstrip()
        if line.startswith("|"):
            table_lines = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                table_lines.append(lines[index].strip())
                index += 1
            rows = []
            for row_index, table_line in enumerate(table_lines):
                cells = [clean_markdown(cell.strip()) for cell in table_line.strip("|").split("|")]
                if row_index == 1 and all(re.fullmatch(r":?-{2,}:?", cell) for cell in cells):
                    continue
                rows.append(cells)
            if rows:
                table = document.add_table(rows=len(rows), cols=len(rows[0]))
                table.style = "Light Grid Accent 1"
                for row_index, row in enumerate(rows):
                    for column_index, value in enumerate(row):
                        cell = table.cell(row_index, column_index)
                        cell.text = value
                        for paragraph in cell.paragraphs:
                            for run in paragraph.runs:
                                run.font.name = "Times New Roman"
                            run.font.size = Pt(9)
                document.add_paragraph()
            continue

        if line.startswith("# "):
            title = clean_markdown(line[2:])
            if not title_seen:
                paragraph = document.add_paragraph(style="Title")
                paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
                running = paragraph.add_run(title)
                running.bold = True
                running.font.size = Pt(18)
                document.add_paragraph()
                title_seen = True
        elif line.startswith("## "):
            heading = clean_markdown(line[3:])
            if heading.lower() in {"introduction", "results", "discussion", "methods"}:
                document.add_page_break()
            document.add_heading(heading, level=1)
        elif line.startswith("### "):
            document.add_heading(clean_markdown(line[4:]), level=2)
        elif line.startswith("- "):
            document.add_paragraph(clean_markdown(line[2:]), style="List Bullet")
        elif re.match(r"^\d+\. ", line):
            document.add_paragraph(
                clean_markdown(re.sub(r"^\d+\. ", "", line)),
                style="List Number",
            )
        elif line.strip() and line.strip() != "---":
            paragraph = document.add_paragraph(clean_markdown(line))
            paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        index += 1

    document.add_page_break()
    document.add_heading("Figures", level=1)
    captions = {
        "figure1_analysis_flow.png": "Figure 1. Four-stage technical roadmap from credible-set filtering to calibrated macrophage prioritization.",
        "figure2_quantile_vs_specificity.png": "Figure 2. Genome-calibrated Atlas RNA-seq effect versus live macrophage specificity.",
        "figure3_top_candidate_heatmap.png": "Figure 3. Top prioritized non-coding loci across genetic, Atlas, live and donor-derived evidence.",
        "figure4_ets2_positive_control.png": "Figure 4. Positive control at rs2836882 showing ETS2 down-regulation in intestinal macrophages.",
        "figure5_ets2_background_null.png": "Figure 5. Positive-control effect in the empirical common-variant background.",
        "figure6_active_ibd_validation.png": "Figure 6. Activity-state evaluation and repeated-pool stability of prioritized loci.",
        "figure7_all_l2g_pathways.png": "Figure 7. External GWAS locus-block bootstrap after retaining all L2G target genes.",
        "figure8_external_experimental_evidence.png": "Figure 8. External experimental evidence for Tier 1 candidates.",
        "figure9_perturbseq_validation.png": "Figure 9. GSE327644 primary human macrophage Perturb-seq.",
        "figure10_evidence_robustness.png": "Figure 10. Robustness and orthogonal secondary analyses.",
        "figure11_benchmark_translation.png": "Figure 11. External benchmarking, multi-cohort triangulation and translational annotation.",
    }
    for filename, caption in captions.items():
        path = figure_dir / filename
        if path.exists():
            document.add_picture(str(path), width=Inches(6.4))
            cap = document.add_paragraph(caption)
            cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            document.add_page_break()
    document.save(output_path)


def build_workbook(table_dir: Path, output_path: Path) -> None:
    evidence = json.loads(
        (table_dir / "candidate_evidence_summary.json").read_text(encoding="utf-8")
    )
    experimental = json.loads(
        (table_dir / "experimental_evidence_summary.json").read_text(
            encoding="utf-8"
        )
    )
    eqtl_catalogue = json.loads(
        (table_dir / "experimental_eqtl_catalogue_summary.json").read_text(
            encoding="utf-8"
        )
    )
    eqtl_benchmark = json.loads(
        (table_dir / "external_eqtl_benchmark_summary.json").read_text(
            encoding="utf-8"
        )
    )
    coloc_summary = json.loads(
        (table_dir / "top10_eqtl_coloc_summary.json").read_text(
            encoding="utf-8"
        )
    )
    clinical_summary = json.loads(
        (table_dir / "top10_clinical_translation_summary.json").read_text(
            encoding="utf-8"
        )
    )
    transport = pd.read_csv(
        table_dir / "benchmark_cohort_transportability.csv"
    )
    prkcq = pd.read_csv(
        table_dir / "prkcq_repeated_pool_classification.csv"
    ).iloc[0]
    sheets = {
        "Overview": pd.DataFrame(
            [
                ["Primary study", "GCST90292538"],
                ["Primary non-coding variants", 395],
                ["Independent loci", 281],
                ["Computational Tier 1", 10],
                ["Computational Tier 2", 57],
                ["Final Tier 1 total", evidence["final_tier1_total"]],
                ["Final Tier 1 persistent", 9],
                ["Final Tier 1 attenuated", 1],
                ["Tier 1 target-uncertain", evidence["tier1_target_uncertain"]],
                ["Tier 1 externally supported", evidence["tier1_external_supported"]],
                ["Tier 1 exact-variant experimental", experimental["tier1_variant_level_experimental"]],
                ["Tier 1 gene-level experimental", experimental["tier1_gene_level_experimental"]],
                ["Tier 1 tested-neutral", experimental["tier1_tested_neutral"]],
                ["Tier 1 Perturb-seq effective", experimental["tier1_perturbseq_effective"]],
                ["Tier 1 Perturb-seq inconclusive", experimental["tier1_perturbseq_inconclusive"]],
                ["Tier 1 variants in gut cCRE", experimental["tier1_scatac_ccre"]],
                ["Public myeloid eQTL datasets", eqtl_catalogue["datasets_queried"]],
                ["Frequency-matched background", 71],
                ["Positive control", "rs2836882 / ETS2"],
                ["PRKCQ repeated-pool class", prkcq["classification"]],
                ["External cohorts", "UKB_WGS_2025; MVP_2024"],
                ["Locus-block bootstrap replicates", 2000],
                [
                    "Exact-variant eQTL benchmark pairs",
                    eqtl_benchmark["variants_queried"],
                ],
                [
                    "Top10 matched eQTL hits",
                    f"{eqtl_benchmark['top10_matched_nominal_pairs']}/10",
                ],
                [
                    "Top10 mismatched-target eQTL controls",
                    f"{eqtl_benchmark['top10_mismatched_nominal_pairs']}/10",
                ],
                [
                    "Formal Top10 colocalizations",
                    f"{coloc_summary['coloc_supported_h4_ge_0_8']}/{coloc_summary['coloc_evaluable']}",
                ],
                [
                    "UKB-to-MVP transport AUROC",
                    float(transport.iloc[0]["auroc"]),
                ],
                [
                    "MVP-to-UKB transport AUROC",
                    float(transport.iloc[1]["auroc"]),
                ],
                [
                    "Top10 tractable targets",
                    clinical_summary["targets_with_tractability"],
                ],
                [
                    "Top10 targets with clinical candidates",
                    clinical_summary["targets_with_clinical_candidates"],
                ],
            ],
            columns=["Metric", "Value"],
        ),
        "Top20": pd.read_csv(
            table_dir / "ranked_lead_variants_final.csv"
        ).head(20),
        "All_candidates": pd.read_csv(
            table_dir / "ranked_variant_gene_candidates_final.csv"
        ),
        "Variant_sets": pd.read_csv(
            ROOT / "data" / "processed" / "variant_set_summary.csv"
        ),
        "Controls": pd.read_csv(table_dir / "control_variant_effects.csv"),
        "Positive_control_atlas": pd.read_csv(
            table_dir / "positive_control_ets2_atlas.csv"
        ),
        "Active_IBD_summary": flatten_json(
            table_dir / "active_ibd_validation_summary.json"
        ),
        "Active_IBD_leads": pd.read_csv(
            table_dir / "ranked_lead_variants_active_ibd.csv"
        ),
        "Pathway_corrected_ORA": pd.read_csv(
            table_dir / "all_l2g_pathway_ora_targeted.csv"
        ),
        "Pathway_all_L2G_ranks": pd.read_csv(
            table_dir / "all_l2g_gene_rankings.csv"
        ),
        "Pathway_prerank": pd.read_csv(
            table_dir / "all_l2g_pathway_prerank_targeted.csv"
        ),
        "Pathway_weighted_tests": pd.read_csv(
            table_dir / "pathway_weighted_rank_tests.csv"
        ),
        "All_variant_pathway": pd.read_csv(
            table_dir / "all_variant_pathway_rank_tests.csv"
        ),
        "Cross_study_L2G_counts": pd.read_csv(
            table_dir / "cross_study_l2g_target_counts.csv"
        ),
        "External_GWAS_studies": pd.read_csv(
            table_dir / "external_ibd_study_summary.csv"
        ),
        "External_gene_rankings": pd.read_csv(
            table_dir / "external_gwas_gene_rankings.csv"
        ),
        "External_GWAS_bootstrap": pd.read_csv(
            table_dir / "external_gwas_block_bootstrap_pathways.csv"
        ),
        "PRKCQ_repeated_scores": pd.read_csv(
            table_dir / "prkcq_repeated_pool_scores.csv"
        ),
        "PRKCQ_repeated_summary": pd.read_csv(
            table_dir / "prkcq_repeated_pool_summary.csv"
        ),
        "Candidate_evidence": flatten_json(
            table_dir / "candidate_evidence_summary.json"
        ),
        "Experimental_evidence": pd.read_csv(
            table_dir / "external_experimental_evidence.csv"
        ),
        "GTEx_PRKCB_eQTL": pd.read_csv(
            table_dir / "experimental_gtex_prkcb_eqtl.csv"
        ),
        "ETS2_macrophage_MPRA": pd.read_csv(
            table_dir / "experimental_ets2_mpra_context_effects.csv"
        ),
        "Macrophage_CRISPR_genes": pd.read_csv(
            table_dir / "experimental_primary_macrophage_crispr_genes.csv"
        ),
        "Macrophage_CRISPR_guides": pd.read_csv(
            table_dir / "experimental_primary_macrophage_crispr_guides.csv"
        ),
        "THP1_OE_engagement": pd.read_csv(
            table_dir / "experimental_thp1_overexpression_engagement.csv"
        ),
        "Experimental_summary": flatten_json(
            table_dir / "experimental_evidence_summary.json"
        ),
        "eQTL_Catalogue_hits": pd.read_csv(
            table_dir / "experimental_eqtl_catalogue_remaining_variants.csv"
        ),
        "eQTL_Catalogue_coverage": pd.read_csv(
            table_dir / "experimental_eqtl_catalogue_coverage.csv"
        ),
        "PerturbSeq_targets": pd.read_csv(
            table_dir / "perturbseq_gse327644_target_knockdown.csv"
        ),
        "PerturbSeq_DE_summary": pd.read_csv(
            table_dir / "perturbseq_gse327644_differential_summary.csv"
        ),
        "PerturbSeq_guides": pd.read_csv(
            table_dir / "perturbseq_gse327644_guide_effects.csv"
        ),
        "PerturbSeq_pathways": pd.read_csv(
            table_dir / "perturbseq_gse327644_pathway_effects.csv"
        ),
        "PerturbSeq_UMI_thresholds": pd.read_csv(
            table_dir / "perturbseq_umi_threshold_sensitivity.csv"
        ),
        "PerturbSeq_leave_guide": pd.read_csv(
            table_dir / "perturbseq_leave_one_guide_out.csv"
        ),
        "PerturbSeq_leave_donor": pd.read_csv(
            table_dir / "perturbseq_leave_one_donor_out.csv"
        ),
        "PerturbSeq_mixed_models": pd.read_csv(
            table_dir / "perturbseq_mixed_model_effects.csv"
        ),
        "CAMK2A_coloc_sensitivity": pd.read_csv(
            table_dir / "camk2a_coloc_sensitivity.csv"
        ),
        "Gut_scATAC_variants": pd.read_csv(
            table_dir / "external_gut_scatac_variant_annotation.csv"
        ),
        "scATAC_topic_percentiles": pd.read_csv(
            table_dir / "scatac_pdgfb_variant_topic_percentiles.csv"
        ),
        "GSE175685_GSEA": pd.read_csv(
            table_dir / "gse175685_prkcb_chsy3_gsea.csv"
        ),
        "GSEA_growth_sensitivity": pd.read_csv(
            table_dir / "gse175685_pathway_confounding_sensitivity.csv"
        ),
        "Claim_calibration": pd.read_csv(
            table_dir / "claim_calibration.csv"
        ),
        "Negative_control": pd.read_csv(table_dir / "negative_control_il2ra.csv"),
        "Threshold_sensitivity": pd.read_csv(table_dir / "threshold_sensitivity.csv"),
        "Score_overlap": pd.read_csv(table_dir / "score_topk_overlap.csv"),
        "Atlas_top_genes": pd.read_csv(table_dir / "atlas_variant_top_genes.csv"),
        "Live_summary": pd.read_csv(table_dir / "live_condition_score_summary.csv"),
        "Benchmark_metrics": pd.read_csv(
            table_dir / "benchmark_model_metrics.csv"
        ),
        "Benchmark_deltas": pd.read_csv(
            table_dir / "benchmark_model_deltas.csv"
        ),
        "Benchmark_calibration": pd.read_csv(
            table_dir / "benchmark_calibration_metrics.csv"
        ),
        "Benchmark_transport": pd.read_csv(
            table_dir / "benchmark_cohort_transportability.csv"
        ),
        "Benchmark_heldout": pd.read_csv(
            table_dir / "benchmark_heldout_validation.csv"
        ),
        "eQTL_variant_benchmark": pd.read_csv(
            table_dir / "external_eqtl_variant_benchmark.csv"
        ),
        "Top10_eQTL_triangulation": pd.read_csv(
            table_dir / "external_eqtl_top10_triangulation.csv"
        ),
        "Top10_coloc": pd.read_csv(
            table_dir / "top10_eqtl_coloc_results.csv"
        ),
        "Top10_coloc_sensitivity": pd.read_csv(
            table_dir / "top10_eqtl_coloc_sensitivity_summary.csv"
        ),
        "Clinical_translation": pd.read_csv(
            table_dir / "top10_clinical_translation_annotation.csv"
        ),
        "Benchmark_coverage": flatten_json(
            table_dir / "benchmark_coverage_summary.json"
        ),
        "Benchmark_priority_sensitivity": pd.read_csv(
            table_dir / "benchmark_priority_sensitivity.csv"
        ),
    }
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    workbook = load_workbook(output_path)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    for worksheet in workbook.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        for cell in worksheet[1]:
            cell.fill = header_fill
            cell.font = Font(color="FFFFFF", bold=True)
            cell.alignment = Alignment(vertical="center")
        for column in worksheet.columns:
            values = [str(cell.value) if cell.value is not None else "" for cell in column]
            width = min(max(len(value) for value in values) + 2, 55)
            worksheet.column_dimensions[get_column_letter(column[0].column)].width = width
    workbook.save(output_path)


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    figure_dir = ROOT / "results" / "figures"
    manuscript_dir = ROOT / "manuscript"
    manuscript_dir.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=2) as executor:
        docx_future = executor.submit(
            build_docx,
            manuscript_dir / "manuscript.md",
            manuscript_dir / "Manuscript_AlphaGenome2_IBD_general.docx",
            figure_dir,
        )
        workbook_future = executor.submit(
            build_workbook,
            table_dir,
            ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx",
        )
        docx_future.result()
        workbook_future.result()
    print("Final DOCX and XLSX written.")


if __name__ == "__main__":
    main()
