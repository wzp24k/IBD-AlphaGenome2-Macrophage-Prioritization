$ErrorActionPreference = "Stop"
$python = Join-Path $PSScriptRoot "work\python312\python.exe"

$scripts = @(
    "03_query_atlas.py",
    "04_score_live_conditions.py",
    "05_annotate_gnomad.py",
    "06b_annotate_ensembl_background.py",
    "07_build_cell_conditioning.py",
    "08_score_custom_macrophage.py",
    "09_score_background_live.py",
    "10_rank_candidates.py",
    "11_make_figures.py",
    "12_sensitivity_and_controls.py",
    "16_build_active_ibd_conditioning.py",
    "17_score_active_ibd_conditions.py",
    "18_analyze_active_ibd.py",
    "27_fetch_external_ibd_gwas.py",
    "28_external_l2g_block_bootstrap.py",
    "29_prkcq_repeated_pool_validation.py",
    "30_finalize_candidate_evidence.py",
    "34_analyze_gse327644_perturbseq.py",
    "38_perturbseq_robustness.py",
    "36_integrate_gut_scatac_variants.py",
    "39_scatac_topic_analysis.py",
    "37_camk2a_coloc.py",
    "40_gse175685_gsea.py",
    "43_claim_calibration_and_gsea_sensitivity.py",
    "44_verify_references.py",
    "45_audit_and_normalize_citations.py",
    "46_external_eqtl_variant_benchmark.py",
    "47_benchmark_models_and_transportability.py",
    "48_top10_coloc_and_negative_controls.py",
    "49_clinical_translation_annotation.py",
    "50_top10_coloc_sensitivity.py",
    "51_benchmark_figure.py",
    "31_integrate_external_experimental_data.py",
    "33_eqtl_catalogue_remaining_variants.py",
    "32_external_experimental_figure.py",
    "35_perturbseq_figure.py",
    "41_evidence_strengthening_figure.py",
    "19_active_ibd_figures.py",
    "20_pathway_reenrichment.py",
    "21_all_l2g_pathway.py",
    "23_weighted_pathway_rank_tests.py",
    "24_all_variant_l2g_pathway.py",
    "22_pathway_figure.py",
    "14_validate_positive_control.py",
    "13_build_deliverables.py",
    "25_build_template_docx.py",
    "26_template_docx_qc.py",
    "15_final_qc.py"
)

foreach ($script in $scripts) {
    Write-Host "Running $script"
    & $python (Join-Path $PSScriptRoot "scripts\$script")
    if ($LASTEXITCODE -ne 0) {
        throw "Pipeline step failed: $script"
    }
}

Write-Host "Core pipeline completed."
