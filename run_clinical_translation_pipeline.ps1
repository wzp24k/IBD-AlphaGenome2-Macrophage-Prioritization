$ErrorActionPreference = "Stop"
$python = Join-Path $PSScriptRoot "work\python312\python.exe"
$env:OMP_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:NUMEXPR_NUM_THREADS = "1"

$scripts = @(
    "72_cd_uc_cross_ancestry_validation.py",
    "73_pqtl_mr_coloc_phewas.py",
    "74_treatment_response_programs.py",
    "75_clinical_readiness_matrix.py",
    "83_parallel_clinical_robustness.py",
    "76_clinical_translation_figure.py",
    "77_extend_manuscript_clinical_translation.py",
    "79_fix_1119_figure_numbering.py",
    "80_update_1119_abstract_and_key_messages.py",
    "81_update_strobe_mr_1119.py",
    "66_build_supplementary_tables.py",
    "71_build_title_page.py",
    "78_audit_clinical_translation_extension.py",
    "82_build_1119_submission_package.py"
)

foreach ($script in $scripts) {
    Write-Host "Running $script"
    & $python (Join-Path $PSScriptRoot "scripts\$script")
    if ($LASTEXITCODE -ne 0) {
        throw "Clinical translation pipeline step failed: $script"
    }
}

Write-Host "Clinical translation pipeline completed."
