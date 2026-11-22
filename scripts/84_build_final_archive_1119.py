from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "_archive" / "IBD_1119_FINAL_20261122.zip"
MANIFEST = ARCHIVE.with_suffix(".sha256.json")

FILES = [
    ROOT / "1119manuscript.docx",
    ROOT / "Title_Page.docx",
    ROOT / "Supplementary_Tables.docx",
    ROOT / "STROBE-MR_Checklist_1119.docx",
    ROOT / "IBD_1119_clinical_translation_submission_package.zip",
    ROOT / "PRKCB_Sun2018_cis_pQTL_20261122.zip",
    ROOT / "results" / "CLINICAL_TRANSLATION_EXTENSION_AUDIT.md",
    ROOT / "results" / "PUBLIC_PRKCB_EVIDENCE_RECHECK.md",
    ROOT / "results" / "clinical_translation_extension_audit.json",
    ROOT / "results" / "tables" / "public_prkcb_evidence_recheck.json",
    ROOT / "results" / "tables" / "benchmark_incremental_model_metrics.csv",
    ROOT / "results" / "tables" / "benchmark_incremental_model_deltas.csv",
    ROOT / "results" / "tables" / "benchmark_incremental_model_summary.json",
    ROOT / "results" / "tables" / "prkcb_sun2018_pqtl_summary.json",
    ROOT / "results" / "1119_reference_audit.json",
    ROOT / "results" / "translational_narrative_audit.json",
    ROOT / "results" / "title_page_audit.json",
    ROOT / "results" / "strobe_mr_1119_audit.json",
    ROOT / "results" / "package_manifest_1119.json",
    ROOT / "results" / "FINAL_DELIVERY_1119.json",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    missing = [path for path in FILES if not path.exists()]
    if missing:
        raise FileNotFoundError(missing)

    with zipfile.ZipFile(ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in FILES:
            archive.write(path, path.relative_to(ROOT))

    manifest = {
        "generated": datetime.now().astimezone().isoformat(),
        "archive": str(ARCHIVE),
        "sha256": sha256(ARCHIVE),
        "files": [
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in FILES
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
