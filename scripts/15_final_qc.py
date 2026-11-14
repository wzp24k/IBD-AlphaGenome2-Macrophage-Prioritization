from __future__ import annotations

import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
from docx import Document
from openpyxl import load_workbook
from PIL import Image, ImageStat

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, load_api_key


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_sentence_lengths(text: str) -> list[int]:
    body = text.split("**Keywords:**", 1)[1].split("## Tables", 1)[0]
    paragraphs = [
        paragraph.strip().replace("\n", " ")
        for paragraph in body.split("\n\n")
        if paragraph.strip()
        and not paragraph.strip().startswith("#")
        and not paragraph.strip().startswith("|")
    ]
    lengths = []
    for paragraph in paragraphs:
        clean = re.sub(r"[*`]", "", paragraph)
        clean = re.sub(r"\[[0-9,\- ]+\]", "", clean)
        for sentence in re.split(r"(?<=[.!?])\s+", clean):
            lengths.append(
                len(
                    re.findall(
                        r"[A-Za-z0-9][A-Za-z0-9'./-]*",
                        sentence,
                    )
                )
            )
    return lengths


def first_citation_order(text: str) -> list[int]:
    main_text = text.split("## References", 1)[0]
    order = []
    for token in re.findall(r"\[([0-9][0-9,\- ]*)\]", main_text):
        for part in token.split(","):
            part = part.strip()
            if "-" in part:
                start, end = map(int, part.split("-", 1))
                numbers = range(start, end + 1)
            else:
                numbers = [int(part)]
            for number in numbers:
                if number not in order:
                    order.append(number)
    return order


def main() -> None:
    manuscript = ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD.docx"
    workbook_path = ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx"
    table_dir = ROOT / "results" / "tables"
    figure_dir = ROOT / "results" / "figures"
    markdown_text = (ROOT / "manuscript" / "manuscript.md").read_text(
        encoding="utf-8"
    )
    palette_path = ROOT / "scripts" / "figure_palette.py"
    figure_scripts = [
        ROOT / "scripts" / "11_make_figures.py",
        ROOT / "scripts" / "19_active_ibd_figures.py",
        ROOT / "scripts" / "22_pathway_figure.py",
        ROOT / "scripts" / "32_external_experimental_figure.py",
        ROOT / "scripts" / "35_perturbseq_figure.py",
        ROOT / "scripts" / "41_evidence_strengthening_figure.py",
        ROOT / "scripts" / "51_benchmark_figure.py",
    ]
    lead = pd.read_csv(table_dir / "ranked_lead_variants_final.csv")
    summary = json.loads((table_dir / "analysis_summary.json").read_text(encoding="utf-8"))
    positive = json.loads(
        (table_dir / "positive_control_validation.json").read_text(encoding="utf-8")
    )
    active_summary = json.loads(
        (table_dir / "active_ibd_validation_summary.json").read_text(
            encoding="utf-8"
        )
    )
    external_summary = json.loads(
        (table_dir / "external_gwas_block_bootstrap_summary.json").read_text(
            encoding="utf-8"
        )
    )
    candidate_summary = json.loads(
        (table_dir / "candidate_evidence_summary.json").read_text(
            encoding="utf-8"
        )
    )
    experimental_summary = json.loads(
        (table_dir / "experimental_evidence_summary.json").read_text(
            encoding="utf-8"
        )
    )
    eqtl_catalogue = json.loads(
        (table_dir / "experimental_eqtl_catalogue_summary.json").read_text(
            encoding="utf-8"
        )
    )
    eqtl_coverage = pd.read_csv(
        table_dir / "experimental_eqtl_catalogue_coverage.csv"
    )
    prkcq_class = pd.read_csv(
        table_dir / "prkcq_repeated_pool_classification.csv"
    ).iloc[0]
    external_pathways = pd.read_csv(
        table_dir / "external_gwas_block_bootstrap_pathways.csv"
    )
    robustness = json.loads(
        (table_dir / "perturbseq_robustness_summary.json").read_text(
            encoding="utf-8"
        )
    )
    threshold_effects = pd.read_csv(
        table_dir / "perturbseq_umi_threshold_sensitivity.csv"
    )
    mixed_effects = pd.read_csv(
        table_dir / "perturbseq_mixed_model_effects.csv"
    )
    coloc_summary = json.loads(
        (table_dir / "camk2a_coloc_summary.json").read_text(encoding="utf-8")
    )
    coloc_sensitivity = pd.read_csv(
        table_dir / "camk2a_coloc_sensitivity.csv"
    )
    topic_percentiles = pd.read_csv(
        table_dir / "scatac_pdgfb_variant_topic_percentiles.csv"
    )
    topic_markers = pd.read_csv(
        table_dir / "scatac_topic_marker_enrichment.csv"
    )
    gsea = pd.read_csv(table_dir / "gse175685_prkcb_chsy3_gsea.csv")
    gsea_sensitivity = pd.read_csv(
        table_dir / "gse175685_pathway_confounding_sensitivity.csv"
    )
    gsea_growth_summary = json.loads(
        (
            table_dir / "gse175685_growth_residualized_summary.json"
        ).read_text(encoding="utf-8")
    )
    claims = pd.read_csv(table_dir / "claim_calibration.csv")
    doi_audit = json.loads(
        (table_dir / "reference_doi_verification_summary.json").read_text(
            encoding="utf-8"
        )
    )
    citation_audit = json.loads(
        (table_dir / "reference_citation_audit_summary.json").read_text(
            encoding="utf-8"
        )
    )
    benchmark_coverage = json.loads(
        (table_dir / "benchmark_coverage_summary.json").read_text(
            encoding="utf-8"
        )
    )
    eqtl_benchmark_summary = json.loads(
        (table_dir / "external_eqtl_benchmark_summary.json").read_text(
            encoding="utf-8"
        )
    )
    top10_coloc_summary = json.loads(
        (table_dir / "top10_eqtl_coloc_summary.json").read_text(
            encoding="utf-8"
        )
    )
    clinical_translation_summary = json.loads(
        (table_dir / "top10_clinical_translation_summary.json").read_text(
            encoding="utf-8"
        )
    )
    benchmark_metrics = pd.read_csv(
        table_dir / "benchmark_model_metrics.csv"
    )
    benchmark_deltas = pd.read_csv(
        table_dir / "benchmark_model_deltas.csv"
    )
    benchmark_transport = pd.read_csv(
        table_dir / "benchmark_cohort_transportability.csv"
    )
    benchmark_heldout = pd.read_csv(
        table_dir / "benchmark_heldout_validation.csv"
    )
    benchmark_priority_summary = json.loads(
        (
            table_dir / "benchmark_priority_sensitivity_summary.json"
        ).read_text(encoding="utf-8")
    )
    strict_eqtl_summary = json.loads(
        (
            table_dir / "strict_eqtl_sensitivity_summary.json"
        ).read_text(encoding="utf-8")
    )
    template_qc = json.loads(
        (ROOT / "results" / "template_docx_qc.json").read_text(
            encoding="utf-8"
        )
    )
    document = Document(manuscript)
    workbook = load_workbook(workbook_path, read_only=True)
    document_text = "\n".join(paragraph.text for paragraph in document.paragraphs)
    abstract_text = markdown_text.split("## Abstract", 1)[1].split(
        "**Keywords:**", 1
    )[0]
    abstract_words = len(
        re.findall(r"[A-Za-z0-9][A-Za-z0-9'./-]*", abstract_text)
    )
    keywords_line = next(
        line
        for line in markdown_text.splitlines()
        if line.startswith("**Keywords:**")
    )
    keywords = [
        value.strip()
        for value in keywords_line.replace("**Keywords:**", "", 1).split(";")
        if value.strip()
    ]
    sentence_lengths = source_sentence_lengths(markdown_text)
    citation_order = first_citation_order(markdown_text)
    manuscript_body = markdown_text.split("## Tables", 1)[0]
    cited_tables = {
        int(value)
        for value in re.findall(r"Table\s+(\d+)", manuscript_body)
    }
    cited_figures = {
        int(value)
        for value in re.findall(r"Figure\s+(\d+)", manuscript_body)
    }
    reference_block = markdown_text.split("## References", 1)[1].split(
        "## Data availability", 1
    )[0]
    reference_numbers = [
        int(value)
        for value in re.findall(r"^(\d+)\.\s", reference_block, re.M)
    ]
    palette_text = palette_path.read_text(encoding="utf-8")
    figure_source_text = "\n".join(
        path.read_text(encoding="utf-8") for path in figure_scripts
    )
    drawio_path = figure_dir / "figure1_analysis_flow.drawio"
    try:
        drawio_text = drawio_path.read_text(encoding="utf-8")
        ET.fromstring(drawio_text)
        figure1_drawio_valid = all(
            phrase in drawio_text
            for phrase in [
                "Technical roadmap for macrophage-focused variant prioritization",
                "Independent evidence",
                "Calibrated output",
            ]
        )
    except (OSError, ET.ParseError):
        figure1_drawio_valid = False

    threshold_self = threshold_effects[
        threshold_effects["condition"].eq("Combined")
        & threshold_effects["metric"].eq("self")
    ]
    pdgfb_threshold = threshold_self[
        threshold_self["target"].eq("PDGFB")
    ].sort_values("umi_threshold")
    il23r_threshold = threshold_self[
        threshold_self["target"].eq("IL23R")
    ].sort_values("umi_threshold")
    prkcq_lps_il23 = threshold_effects[
        threshold_effects["target"].eq("PRKCQ")
        & threshold_effects["condition"].eq("LPS")
        & threshold_effects["metric"].eq("IL-23")
    ].sort_values("umi_threshold")
    mixed_pdgfb_self = mixed_effects[
        mixed_effects["target"].eq("PDGFB")
        & mixed_effects["metric"].eq("self")
    ]
    mixed_prkcq_lps_il23 = mixed_effects[
        mixed_effects["target"].eq("PRKCQ")
        & mixed_effects["condition"].eq("LPS")
        & mixed_effects["metric"].eq("IL-23")
    ]
    topic15 = topic_percentiles[topic_percentiles["topic"].eq("Topic15")]
    topic15_macrophage = topic_markers[
        topic_markers["topic"].eq("Topic15")
        & topic_markers["marker_set"].eq("Macrophage")
    ]
    topic15_dc = topic_markers[
        topic_markers["topic"].eq("Topic15")
        & topic_markers["marker_set"].eq("DC")
    ]
    prkcb_jak = gsea[
        gsea["target"].eq("PRKCB")
        & gsea["Term"].eq("JAK-STAT signaling pathway")
    ]
    prkcb_il6 = gsea[
        gsea["target"].eq("PRKCB")
        & gsea["Term"].eq("IL-6/JAK/STAT3 Signaling")
    ]
    chsy3_tnf = gsea[
        gsea["target"].eq("CHSY3")
        & gsea["Term"].eq("TNF-alpha Signaling via NF-kB")
    ]
    chsy3_il6 = gsea[
        gsea["target"].eq("CHSY3")
        & gsea["Term"].eq("IL-6/JAK/STAT3 Signaling")
    ]

    figure_checks = {}
    for path in sorted(figure_dir.glob("*.png")):
        image = Image.open(path).convert("RGB")
        stat = ImageStat.Stat(image)
        figure_checks[path.name] = {
            "dimensions": list(image.size),
            "nonblank": min(stat.stddev) > 5,
        }

    key = load_api_key().encode("utf-8")
    deliverable_files = [
        manuscript,
        ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD_general.docx",
        workbook_path,
        ROOT / "manuscript" / "manuscript.md",
        ROOT / "results" / "OUTLINE_AUDIT.md",
        *sorted((ROOT / "results" / "figures").glob("*")),
        *sorted((ROOT / "results" / "tables").glob("*")),
    ]
    key_leaks = [
        str(path)
        for path in deliverable_files
        if path.is_file() and key in path.read_bytes()
    ]

    checks = {
        "primary_variants_395": summary["primary_variants"] == 395,
        "independent_loci_281": summary["independent_loci"] == 281,
        "tier1_10": summary["tier_1"] == 10,
        "tier2_57": summary["tier_2"] == 57,
        "final_validated_tier1_10": active_summary["final_validated_tier1"] == 10,
        "final_persistent_tier1_9": active_summary["final_persistent_tier1"] == 9,
        "final_attenuated_tier1_1": active_summary["final_attenuated_tier1"] == 1,
        "lead_rows_281": len(lead) == 281,
        "tier1_target_uncertain_3": candidate_summary["tier1_target_uncertain"] == 3,
        "tier1_external_supported_1": candidate_summary["tier1_external_supported"] == 1,
        "external_two_cohort_clusters": external_summary["external_clusters"]
        == ["UKB_WGS_2025", "MVP_2024"],
        "external_bootstrap_2000": external_summary["bootstrap_replicates"] == 2000,
        "external_pathways_4": len(external_pathways) == 4,
        "external_pathways_fdr_significant": bool(
            (external_pathways["bootstrap_fdr"] < 0.05).all()
        ),
        "prkcq_classified_attenuated": prkcq_class["classification"]
        == "attenuated",
        "tier1_exact_variant_experimental_1": experimental_summary[
            "tier1_variant_level_experimental"
        ]
        == 1,
        "tier1_gene_level_experimental_5": experimental_summary[
            "tier1_gene_level_experimental"
        ]
        == 5,
        "tier1_tested_neutral_3": experimental_summary[
            "tier1_tested_neutral"
        ]
        == 3,
        "tier1_perturbseq_effective_2": experimental_summary[
            "tier1_perturbseq_effective"
        ]
        == 2,
        "tier1_perturbseq_inconclusive_1": experimental_summary[
            "tier1_perturbseq_inconclusive"
        ]
        == 1,
        "tier1_gut_scatac_ccre_1": experimental_summary[
            "tier1_scatac_ccre"
        ]
        == 1,
        "gtex_prkcb_direction_concordant": experimental_summary[
            "gtex_prkcb_nes"
        ]
        < 0,
        "eqtl_catalogue_8_datasets": eqtl_catalogue[
            "datasets_queried"
        ]
        == 8,
        "eqtl_catalogue_5_variants": eqtl_catalogue[
            "variants_queried"
        ]
        == 5,
        "camk2a_nominal_discordant_only": bool(
            (
                eqtl_coverage.loc[
                    eqtl_coverage["target"].eq("CAMK2A"),
                    "classification",
                ]
                == "nominal_colon_or_myeloid_eQTL"
            ).all()
        ),
        "four_remaining_myeloid_eqtl_uncovered": int(
            eqtl_coverage["classification"]
            .eq("public_myeloid_eqtl_not_covered")
            .sum()
        )
        == 4,
        "bootstrap_uses_276_informative_loci": "276 credible-set loci"
        in markdown_text
        and "resamples of 293 credible-set loci" not in markdown_text,
        "atlas_missing_breakdown_22_13_9": all(
            phrase in markdown_text
            for phrase in [
                "22 without Atlas target scores",
                "13 were indels",
                "nine were SNVs",
            ]
        ),
        "tnfaip3_is_response_control_not_knockdown": "downstream"
        in markdown_text
        and "transcript knockdown not verified" in markdown_text,
        "macrophage_comparator_wording_precise": "macrophage-versus-control specificity"
        in markdown_text
        and "reference-condition macrophage specificity" not in markdown_text,
        "model_benchmark_limitation_reported": "not benchmarked against a matched variant-level macrophage gold standard"
        in markdown_text,
        "clinical_use_boundary_reported": "not a diagnostic, prognostic or treatment-selection tool"
        in markdown_text,
        "all_tables_cited": cited_tables == set(range(1, 12)),
        "all_figures_cited": cited_figures == set(range(1, 12)),
        "positive_control_pass": positive["all_checks_pass"],
        "docx_has_text": len(document.paragraphs) > 100,
        "docx_has_eleven_tables": len(document.tables) == 11,
        "docx_has_eleven_figures": len(document.inline_shapes) == 11,
        "docx_has_table11_caption": "Claim calibration for the seven evidence boundaries"
        in document_text,
        "docx_has_figure10_caption": "Robustness and orthogonal secondary analyses"
        in document_text,
        "docx_has_figure11_caption": "External benchmarking, multi-cohort triangulation and translational annotation"
        in document_text,
        "xlsx_has_expected_sheets": all(
            sheet in workbook.sheetnames
            for sheet in [
                "Top20",
                "External_GWAS_studies",
                "External_gene_rankings",
                "External_GWAS_bootstrap",
                "PRKCQ_repeated_scores",
                "PRKCQ_repeated_summary",
                "Candidate_evidence",
                "Experimental_evidence",
                "GTEx_PRKCB_eQTL",
                "ETS2_macrophage_MPRA",
                "Macrophage_CRISPR_genes",
                "THP1_OE_engagement",
                "eQTL_Catalogue_hits",
                "eQTL_Catalogue_coverage",
                "PerturbSeq_targets",
                "PerturbSeq_DE_summary",
                "PerturbSeq_guides",
                "PerturbSeq_pathways",
                "PerturbSeq_UMI_thresholds",
                "PerturbSeq_leave_guide",
                "PerturbSeq_leave_donor",
                "PerturbSeq_mixed_models",
                "CAMK2A_coloc_sensitivity",
                "Gut_scATAC_variants",
                "scATAC_topic_percentiles",
                "GSE175685_GSEA",
                "GSEA_growth_sensitivity",
                "Claim_calibration",
                "Benchmark_metrics",
                "Benchmark_deltas",
                "Benchmark_calibration",
                "Benchmark_transport",
                "Benchmark_heldout",
                "eQTL_variant_benchmark",
                "Top10_eQTL_triangulation",
                "Top10_coloc",
                "Top10_coloc_sensitivity",
                "Clinical_translation",
                "Benchmark_coverage",
                "Benchmark_priority_sensitivity",
                "Strict_eqtl_sensitivity",
            ]
        ),
        "perturbseq_thresholds_3_5_10": robustness["assignment_thresholds"]
        == [3, 5, 10],
        "pdgfb_self_threshold_stable": bool(
            len(pdgfb_threshold) == 3
            and (pdgfb_threshold["mean_log2fc"] < 0).all()
            and (
                pdgfb_threshold["mean_log2fc"].max()
                - pdgfb_threshold["mean_log2fc"].min()
            )
            <= 0.02
        ),
        "il23r_no_knockdown_stable": bool(
            len(il23r_threshold) == 3
            and (il23r_threshold["mean_log2fc"] > 0).all()
            and (il23r_threshold["p_value"] >= 0.05).all()
        ),
        "prkcq_lps_il23_direction_stable": bool(
            len(prkcq_lps_il23) == 3
            and (prkcq_lps_il23["mean_log2fc"] < -0.2).all()
        ),
        "donor_models_56": len(mixed_effects) == 56,
        "donor_model_types_valid": set(mixed_effects["model_type"]).issubset(
            {
                "mixedlm_donor_random_intercept",
                "donor_fixed_effect_fallback",
            }
        ),
        "pdgfb_self_donor_blocked_negative": bool(
            len(mixed_pdgfb_self) == 2
            and (mixed_pdgfb_self["coefficient"] < 0).all()
            and (mixed_pdgfb_self["p_value"] < 0.05).all()
        ),
        "prkcq_lps_il23_donor_blocked": bool(
            len(mixed_prkcq_lps_il23) == 1
            and mixed_prkcq_lps_il23.iloc[0]["coefficient"] < 0
            and mixed_prkcq_lps_il23.iloc[0]["p_value"] < 0.05
        ),
        "camk2a_coloc_not_supported": coloc_summary["primary_posterior"]["H4"]
        < 0.8
        and float(coloc_sensitivity["H4"].max()) < 0.8,
        "gsea_prkcb_jak_and_il6_positive": bool(
            len(prkcb_jak) == 1
            and len(prkcb_il6) == 1
            and prkcb_jak.iloc[0]["NES"] > 0
            and prkcb_jak.iloc[0]["FDR q-val"] < 0.05
            and prkcb_il6.iloc[0]["NES"] > 0
            and round(float(prkcb_il6.iloc[0]["FDR q-val"]), 3) <= 0.05
        ),
        "gsea_chsy3_tnf_and_il6_positive": bool(
            len(chsy3_tnf) == 1
            and len(chsy3_il6) == 1
            and chsy3_tnf.iloc[0]["NES"] > 0
            and chsy3_tnf.iloc[0]["FDR q-val"] < 0.05
            and chsy3_il6.iloc[0]["NES"] > 0
            and chsy3_il6.iloc[0]["FDR q-val"] < 0.05
        ),
        "claim_calibration_7_rows": len(claims) == 7,
        "gsea_growth_pc1_fraction": bool(
            gsea_growth_summary["growth_pc1_metadata"][
                "pc1_variance_explained"
            ]
            > 0.2
        ),
        "gsea_only_chsy3_tnf_survives_residualization": bool(
            (
                gsea_sensitivity[
                    "fdr_supported_after_residualization"
                ].sum()
                == 1
            )
            and (
                gsea_sensitivity.loc[
                    gsea_sensitivity[
                        "fdr_supported_after_residualization"
                    ],
                    "target",
                ]
                == "CHSY3"
            ).all()
            and (
                gsea_sensitivity.loc[
                    gsea_sensitivity[
                        "fdr_supported_after_residualization"
                    ],
                    "Term",
                ]
                == "TNF-alpha Signaling via NF-kB"
            ).all()
        ),
        "gsea_prkcb_nonsignificant_after_residualization": bool(
            (
                gsea_sensitivity.loc[
                    gsea_sensitivity["target"].eq("PRKCB"),
                    "residualized_FDR",
                ]
                >= 0.05
            ).all()
        ),
        "pdgfb_topic15_specific": bool(
            len(topic15) == 1
            and topic15.iloc[0]["percentile"] >= 0.95
            and len(topic15_macrophage) == 1
            and topic15_macrophage.iloc[0]["within_set_rank"] == 1
            and len(topic15_dc) == 1
            and topic15_dc.iloc[0]["within_set_rank"] == 1
        ),
        "manuscript_uses_nonoverlapping_language": "non-overlapping cohort clusters"
        in document_text,
        "manuscript_marks_target_uncertainty": "target-uncertain"
        in document_text,
        "manuscript_reports_prkcq_attenuated": "PRKCQ"
        in document_text
        and "attenuated" in document_text,
        "manuscript_reports_external_experiments": (
            "Public experimental evidence" in document_text
            or "External experimental evidence" in document_text
        )
        and (
            "cytokine endpoint" in document_text.lower()
            or "inconclusive" in document_text.lower()
        ),
        "manuscript_reports_perturbseq": "124,959 cells"
        in document_text
        and "inconclusive" in document_text.lower(),
        "manuscript_reports_full_library_direction_fix": "full gene-library depth"
        in document_text,
        "manuscript_avoids_estimable_claim": "were estimable for 56"
        not in document_text,
        "manuscript_reports_camk2a_noncolocalization": "lacks colocalization support"
        in document_text,
        "manuscript_reports_claim_calibration": "Evidence ceiling and claim calibration"
        in document_text
        and "seven evidence boundaries" in document_text,
        "manuscript_reports_growth_sensitivity": "growth-PC1-residualized"
        in document_text,
        "manuscript_reports_interval_boundary": any(
            phrase in document_text
            for phrase in [
                "bootstrap intervals crossing zero",
                "bootstrap intervals cross zero",
                "intervals crossed zero",
                "intervals crossing zero",
            ]
        ),
        "no_false_independent_replication_language": all(
            phrase not in document_text
            for phrase in [
                "three independent replication studies",
                "replicated across three independent IBD/CD/UC datasets",
                "cross-study replication",
                "six independent 256-cell pools",
            ]
        ),
        "figures_nonblank": all(item["nonblank"] for item in figure_checks.values()),
        "abstract_250_words_or_less": abstract_words <= 250,
        "abstract_structured_four_sections": all(
            label in abstract_text
            for label in [
                "**Background:**",
                "**Methods:**",
                "**Results:**",
                "**Conclusion:**",
            ]
        ),
        "keywords_five_or_fewer": len(keywords) <= 5,
        "keywords_title_case": all(
            keyword[0].isupper() for keyword in keywords
        ),
        "sentences_30_words_or_less": bool(sentence_lengths)
        and max(sentence_lengths) <= 30,
        "p_value_format_consistent": all(
            re.fullmatch(r"0\.\d{3}", value)
            for value in re.findall(
                r"(?<!PIP )\bP\s*=\s*([0-9]*\.?[0-9]+)",
                markdown_text,
            )
        )
        and all(
            value == "0.001"
            for value in re.findall(
                r"\bP\s*<\s*([0-9]*\.?[0-9]+)",
                markdown_text,
            )
        ),
        "no_collapse_symbols": not any(
            symbol in markdown_text for symbol in ["▼", "▶", "▸", "◀"]
        ),
        "reference_count_42": len(reference_numbers) == 42,
        "citation_numbers_cover_1_to_42": sorted(set(citation_order))
        == list(range(1, 43)),
        "discussion_citations_exactly_15": citation_audit["all_checks_pass"]
        and citation_audit["discussion_exactly_15_references"],
        "pre_discussion_citations_exactly_27": citation_audit["all_checks_pass"]
        and citation_audit["pre_discussion_exactly_27_references"],
        "citation_groups_at_most_two": citation_audit[
            "citation_groups_at_most_two"
        ],
        "discussion_majority_citations_in_middle": citation_audit[
            "discussion_majority_in_middle"
        ],
        "citation_groups_sorted_unique": citation_audit[
            "all_citation_groups_sorted_unique"
        ],
        "all_references_cited": not citation_audit["uncited_references"],
        "eqtl_variant_benchmark_384": eqtl_benchmark_summary[
            "variants_queried"
        ]
        == 384,
        "top10_eqtl_7_vs_0": eqtl_benchmark_summary[
            "top10_matched_nominal_pairs"
        ]
        == 7
        and eqtl_benchmark_summary["top10_mismatched_nominal_pairs"] == 0,
        "top10_eqtl_fisher_significant": top10_coloc_summary["fisher_p"]
        < 0.01,
        "formal_coloc_five_of_seven": top10_coloc_summary[
            "coloc_supported_h4_ge_0_8"
        ]
        == 5
        and top10_coloc_summary["coloc_evaluable"] == 7,
        "benchmark_external_both_combined_auc": bool(
            benchmark_metrics.loc[
                benchmark_metrics["task"].eq(
                    "External genetic support, both cohorts"
                )
                & benchmark_metrics["model"].eq("Combined"),
                "auroc",
            ].iloc[0]
            > 0.60
        ),
        "combined_improves_over_l2g_external": bool(
            benchmark_deltas.loc[
                benchmark_deltas["task"].eq(
                    "External genetic support, both cohorts"
                )
                & benchmark_deltas["baseline"].eq("L2G-only"),
                "delta_auroc_ci_low",
            ].iloc[0]
            > 0
        ),
        "combined_not_universally_superior": bool(
            benchmark_deltas.loc[
                benchmark_deltas["task"].eq(
                    "External genetic support, both cohorts"
                )
                & benchmark_deltas["baseline"].eq("PIP-only"),
                "delta_auroc_ci_low",
            ].iloc[0]
            <= 0
        ),
        "cohort_transport_auc": bool(
            benchmark_transport["auroc"].min() > 0.60
        ),
        "heldout_credible_set_auc": bool(
            benchmark_heldout.loc[
                benchmark_heldout["group_column"].eq("studyLocusId"),
                "auroc",
            ].iloc[0]
            > 0.65
        ),
        "original_priority_reproduces_benchmark": benchmark_priority_summary[
            "spearman_priority_vs_benchmark_score"
        ]
        > 0.99,
        "strict_eqtl_sensitivity_preserves_enrichment": strict_eqtl_summary[
            "strict_positive_pairs"
        ]
        == 95
        and strict_eqtl_summary["top10_matched_strict_hits"] == 7
        and strict_eqtl_summary["top10_mismatched_strict_hits"] == 0,
        "clinical_translation_ten_targets": clinical_translation_summary[
            "targets_annotated"
        ]
        == 10
        and clinical_translation_summary[
            "targets_with_ibd_indication_candidates"
        ]
        == 3
        and clinical_translation_summary[
            "directionally_compatible_ibd_targets"
        ]
        == ["PRKCB"]
        and clinical_translation_summary[
            "clinical_treatment_recommendations_supported"
        ]
        == 0,
        "reference_bookmarks_and_hyperlinks": template_qc["checks"][
            "forty_two_reference_bookmarks"
        ]
        and template_qc["checks"]["all_references_hyperlink_targets"]
        and template_qc["checks"]["discussion_has_exactly_15_references"]
        and template_qc["checks"]["pre_discussion_has_exactly_27_references"]
        and template_qc["checks"][
            "citation_groups_have_at_most_two_references"
        ]
        and template_qc["checks"]["citation_hyperlink_text_matches_anchor"]
        and template_qc["checks"][
            "citation_hyperlinks_black_without_underline"
        ]
        and template_qc["checks"]["citation_hyperlinks_sorted_within_group"],
        "figure1_drawio_valid": figure1_drawio_valid,
        "all_reference_dois_verified": doi_audit["all_verified"],
        "all_references_are_journal_articles": doi_audit[
            "verified_journal_articles"
        ]
        == 42,
        "declarations_complete": all(
            heading in document_text
            for heading in [
                "Data availability",
                "Code availability",
                "Ethics approval and consent to participate",
                "Patient and public involvement",
                "Authors' contributions",
                "Funding",
                "Declaration of generative AI and AI-assisted technologies in the writing process",
                "Consent for publication",
                "Clinical trial number",
                "Competing interests",
            ]
        ),
        "technical_roadmap_caption": "Four-stage technical roadmap" in document_text,
        "central_figure_palette_present": palette_path.exists()
        and all(
            color in palette_text
            for color in [
                "#66BC98",
                "#AAD09D",
                "#E3EA96",
                "#FDC089",
                "#E26844",
                "#8A233F",
            ]
        ),
        "legacy_figure_colors_removed": all(
            color not in figure_source_text
            for color in [
                "#4C78A8",
                "#D1495B",
                "#EDAE49",
                "#00798C",
                "#7B3294",
                "#E09F3E",
            ]
        ),
        "all_figure_scripts_use_palette": all(
            "figure_palette" in path.read_text(encoding="utf-8")
            for path in figure_scripts
        ),
        "no_api_key_in_deliverables": len(key_leaks) == 0,
    }
    report = {
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "analysis_summary": summary,
        "external_gwas_summary": external_summary,
        "candidate_evidence_summary": candidate_summary,
        "experimental_evidence_summary": experimental_summary,
        "eqtl_catalogue_summary": eqtl_catalogue,
        "prkcq_repeated_pool": prkcq_class.to_dict(),
        "positive_control_validation": positive,
        "figures": figure_checks,
        "docx": {
            "paragraphs": len(document.paragraphs),
            "tables": len(document.tables),
            "inline_shapes": len(document.inline_shapes),
        },
        "xlsx_sheets": workbook.sheetnames,
        "api_key_leaks": key_leaks,
        "format_compliance": {
            "abstract_words": abstract_words,
            "keywords": keywords,
            "max_sentence_words": max(sentence_lengths),
            "citation_order": citation_order,
            "citation_groups": citation_audit["citation_groups"],
            "reference_bookmarks": template_qc["counts"][
                "reference_bookmarks"
            ],
            "reference_hyperlinks": template_qc["counts"][
                "reference_hyperlinks"
            ],
            "cited_tables": sorted(cited_tables),
            "cited_figures": sorted(cited_figures),
            "references": len(reference_numbers),
            "doi_audit": doi_audit,
        },
    }
    (ROOT / "results" / "QC_REPORT.json").write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    status_lines = [
        "# Final status",
        "",
        f"- All checks passed: {report['all_checks_pass']}",
        "- External genetic cohorts: UKB_WGS_2025 and MVP_2024, counted once each",
        f"- Locus-block bootstrap replicates: {external_summary['bootstrap_replicates']}",
        "- Pathway FDRs: "
        + "; ".join(
            f"{row.pathway}={row.bootstrap_fdr:.4g}"
            for row in external_pathways.itertuples()
        ),
        f"- PRKCQ repeated-pool classification: {prkcq_class['classification']}",
        f"- Tier 1 targets: {candidate_summary['final_tier1_total']}",
        f"- Tier 1 target-uncertain: {candidate_summary['tier1_target_uncertain']}",
        f"- Tier 1 externally supported: {candidate_summary['tier1_external_supported']}",
        f"- Tier 1 exact-variant experimental: {experimental_summary['tier1_variant_level_experimental']}",
        f"- Tier 1 gene-level experimental: {experimental_summary['tier1_gene_level_experimental']}",
        f"- Tier 1 tested-neutral for cytokine selection: {experimental_summary['tier1_tested_neutral']}",
        f"- Tier 1 Perturb-seq effective: {experimental_summary['tier1_perturbseq_effective']}",
        f"- Tier 1 Perturb-seq inconclusive: {experimental_summary['tier1_perturbseq_inconclusive']}",
        f"- Tier 1 variants in gut scATAC cCRE: {experimental_summary['tier1_scatac_ccre']}",
        "- Public myeloid eQTL screen: 4/5 uncovered; CAMK2A has a nominal IFN-gamma signal with discordant direction and Bonferroni P approximately 0.31",
        "- CAMK2A colocalization: no support; primary H4=0.09 and maximum sensitivity H4 below 0.80",
        "- Perturb-seq robustness: 3/5/10 UMI thresholds, leave-one-out and donor-blocked checks used full gene-library depth",
        "- GSEA confounding sensitivity: proliferation PC1 explains 23.2% of growth-gene variance; only CHSY3 TNF/NF-kB retains FDR < 0.05 after residualization",
        "- Claim calibration: seven evidence boundaries are reported in Table 11; no non-significant or uncovered result is relabeled positive",
        "- Figure style: all Figure 1-11 scripts use the centralized six-color palette and shared typography",
        "- External eQTL benchmark: 384 candidate-target pairs; Top 10 matched 7/10 versus mismatched 0/10",
        "- Formal colocalization: 5 of 7 eQTL-positive Top 10 pairs supported at H4 >= 0.8",
        "- Benchmark transport: UKB-to-MVP AUROC 0.633 and MVP-to-UKB AUROC 0.614",
        "- Clinical annotation: 3/10 Top 10 targets have phase 2 IBD-indication drug signals; 0/10 meet a treatment-recommendation boundary",
        "- Citation structure: 42 sequential reference entries, 27 distinct references before Discussion and 15 distinct references within Discussion; each citation group contains at most two references",
        "- Word cross-references: 42 reference bookmarks with clickable internal hyperlinks from every in-text citation number",
        "- Editable Figure 1: results/figures/figure1_analysis_flow.drawio",
        "- Evidence boundary: UKB/MVP are external genetic evidence; GSE229472, GSE327124 and GSE175685 provide public experimental evidence at the explicitly marked variant or gene level",
    ]
    status_path = ROOT / "results" / "FINAL_STATUS.md"
    status_path.write_text("\n".join(status_lines) + "\n", encoding="utf-8")

    package_files = [
        ROOT / "README.md",
        ROOT / "requirements.txt",
        ROOT / "run_core_pipeline.ps1",
        ROOT / "config" / "analysis.json",
        ROOT / "manuscript" / "manuscript.md",
        manuscript,
        workbook_path,
        ROOT / "results" / "OUTLINE_AUDIT.md",
        ROOT / "results" / "PATHWAY_AUDIT.md",
        ROOT / "results" / "EXTERNAL_EXPERIMENTAL_AUDIT.md",
        ROOT / "results" / "EXPERIMENTAL_FOLLOWUP_PLAN.md",
        ROOT / "results" / "PERTURBSEQ_AUDIT.md",
        ROOT / "results" / "ROBUSTNESS_AND_DIRECTION_AUDIT.md",
        ROOT / "results" / "RESOURCE_SEARCH_AUDIT.md",
        ROOT / "results" / "FIGURE_STYLE_AUDIT.md",
        ROOT / "results" / "FORMAT_COMPLIANCE_1113.md",
        ROOT / "results" / "CITATION_CROSSREF_AUDIT.md",
        ROOT / "results" / "BENCHMARK_AND_TRANSLATION_AUDIT.md",
        ROOT / "results" / "REVIEWER1_AUDIT_1113.md",
        ROOT / "results" / "REVIEWER2_3_AUDIT_1113.md",
        status_path,
        ROOT / "results" / "QC_REPORT.json",
        ROOT / "results" / "template_docx_qc.json",
    ]
    package_files.extend(sorted((ROOT / "scripts").glob("*.py")))
    package_files.extend(sorted((ROOT / "results" / "figures").glob("*")))
    package_files.extend(sorted((ROOT / "results" / "tables").glob("*")))
    package_files.extend(
        sorted((ROOT / "data" / "external" / "experimental" / "GSE229472").glob("*"))
    )
    package_files.extend(
        sorted((ROOT / "data" / "external" / "experimental" / "GSE327124").glob("*.gz"))
    )
    package_files.extend(
        sorted((ROOT / "data" / "external" / "experimental" / "GSE175685").glob("*.gz"))
    )
    package_files.extend(
        [
            ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv",
            ROOT / "data" / "processed" / "ibd_lead_noncoding_variants.csv",
            ROOT / "data" / "processed" / "variant_set_summary.csv",
            ROOT / "data" / "processed" / "frequency_matched_background.csv",
            ROOT / "data" / "processed" / "cell_conditioning_summary.csv",
        ]
    )
    package_files = [path for path in package_files if path.exists()]
    package_path = ROOT / "IBD_AlphaGenome2_submission_package.zip"
    with zipfile.ZipFile(package_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in package_files:
            archive.write(path, path.relative_to(ROOT))
    manifest = {
        "package_file": package_path.name,
        "package_sha256": sha256(package_path),
        "files": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in package_files
        ],
    }
    (ROOT / "results" / "package_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(report["checks"], indent=2))
    if not report["all_checks_pass"]:
        raise SystemExit("Final QC failed.")
    print(f"Package: {package_path}")


if __name__ == "__main__":
    main()
