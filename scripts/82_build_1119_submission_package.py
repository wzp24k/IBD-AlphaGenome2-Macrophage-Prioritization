from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, load_api_key


OUTPUT = ROOT / "IBD_1119_clinical_translation_submission_package.zip"
MANIFEST = ROOT / "results" / "package_manifest_1119.json"
DELIVERY = ROOT / "results" / "FINAL_DELIVERY_1119.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    files = [
        ROOT / "README.md",
        ROOT / "1119manuscript.docx",
        ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD_clinical_translation.docx",
        ROOT / "Title_Page.docx",
        ROOT / "Cover_Letter_Generic.docx",
        ROOT / "STROBE-MR_Checklist_1119.docx",
        ROOT / "Supplementary_Tables.docx",
        ROOT / "Graphical_Abstract.drawio",
        ROOT / "Graphical_Abstract.png",
        ROOT / "Figure_1.drawio",
        ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx",
        ROOT / "results" / "figures" / "figure13_clinical_translation.png",
        ROOT / "results" / "figures" / "figure13_clinical_translation.pdf",
        ROOT / "results" / "CLINICAL_TRANSLATION_EXTENSION_AUDIT.md",
        ROOT / "results" / "clinical_translation_extension_audit.json",
        ROOT / "results" / "1119_reference_audit.json",
        ROOT / "results" / "title_page_audit.json",
        ROOT / "results" / "strobe_mr_1119_audit.json",
        table_dir / "top10_cd_uc_target_support.csv",
        table_dir / "prkcb_cross_ancestry_harmonized.csv",
        table_dir / "prkcb_cross_ancestry_mr.csv",
        table_dir / "prkcb_cross_ancestry_leave_one_out.csv",
        table_dir / "cd_uc_cross_ancestry_summary.json",
        table_dir / "pqtl_cis_instruments.csv",
        table_dir / "pqtl_cis_mr_ibd.csv",
        table_dir / "pqtl_mr_sensitivity.csv",
        table_dir / "pqtl_ibd_harmonized.csv",
        table_dir / "pqtl_ibd_colocalization.csv",
        table_dir / "pqtl_coloc_prior_sensitivity.csv",
        table_dir / "prkcb_pqtl_phewas_mr.csv",
        table_dir / "pqtl_translation_summary.json",
        table_dir / "treatment_response_program_associations.csv",
        table_dir / "treatment_response_sample_scores.csv",
        table_dir / "treatment_response_meta_analysis.csv",
        table_dir / "treatment_response_permutation.csv",
        table_dir / "treatment_response_summary.json",
        table_dir / "clinical_development_readiness_matrix.csv",
        table_dir / "clinical_readiness_rank_sensitivity.csv",
        table_dir / "clinical_readiness_weight_bootstrap.csv",
        table_dir / "clinical_development_readiness_summary.json",
        table_dir / "clinical_translation_robustness_summary.json",
        table_dir / "benchmark_incremental_model_metrics.csv",
        table_dir / "benchmark_incremental_model_deltas.csv",
        table_dir / "benchmark_incremental_model_summary.json",
        table_dir / "public_prkcb_evidence_recheck.json",
        ROOT / "results" / "PUBLIC_PRKCB_EVIDENCE_RECHECK.md",
    ]
    script_names = [
        "72_cd_uc_cross_ancestry_validation.py",
        "73_pqtl_mr_coloc_phewas.py",
        "74_treatment_response_programs.py",
        "75_clinical_readiness_matrix.py",
        "76_clinical_translation_figure.py",
        "77_extend_manuscript_clinical_translation.py",
        "78_audit_clinical_translation_extension.py",
        "79_fix_1119_figure_numbering.py",
        "80_update_1119_abstract_and_key_messages.py",
        "81_update_strobe_mr_1119.py",
        "82_build_1119_submission_package.py",
        "83_parallel_clinical_robustness.py",
        "85_benchmark_incremental_model.py",
        "86_public_evidence_recheck.py",
    ]
    files.extend(ROOT / "scripts" / name for name in script_names)
    files = sorted({path.resolve(): path for path in files if path.exists()}.values())

    key = load_api_key().encode("utf-8")
    leaks = [
        str(path)
        for path in files
        if path.suffix.lower() not in {".parquet", ".h5ad"}
        and path.is_file()
        and key in path.read_bytes()
    ]
    if leaks:
        raise RuntimeError(f"API key detected in package files: {leaks}")

    with zipfile.ZipFile(OUTPUT, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, path.relative_to(ROOT))
    manifest = {
        "generated": datetime.now().astimezone().isoformat(),
        "package": str(OUTPUT),
        "sha256": sha256(OUTPUT),
        "files": [
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in files
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    reference = json.loads(
        (ROOT / "results" / "1119_reference_audit.json").read_text(encoding="utf-8")
    )
    clinical_audit = json.loads(
        (ROOT / "results" / "clinical_translation_extension_audit.json").read_text(
            encoding="utf-8"
        )
    )
    delivery = {
        "generated": datetime.now().astimezone().isoformat(),
        "manuscript": "1119manuscript.docx",
        "title_page": "Title_Page.docx",
        "supplement": "Supplementary_Tables.docx",
        "strobe_mr": "STROBE-MR_Checklist_1119.docx",
        "references": reference["references"],
        "each_reference_cited_once": not reference["repeated_reference_numbers"],
        "all_citations_hyperlinked": reference[
            "all_citation_occurrences_hyperlinked"
        ],
        "long_sentences_over_30": len(reference["long_sentences_over_30"]),
        "clinical_translation_audit_pass": clinical_audit["all_checks_pass"],
        "tables": clinical_audit["tables"],
        "figures": clinical_audit["figures"],
        "package": str(OUTPUT),
        "package_sha256": manifest["sha256"],
        "package_files": len(files),
    }
    DELIVERY.write_text(json.dumps(delivery, indent=2), encoding="utf-8")
    print(json.dumps(delivery, indent=2))


if __name__ == "__main__":
    main()
