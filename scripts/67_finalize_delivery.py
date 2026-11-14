from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    generated = datetime.now().astimezone()
    qc = json.loads(
        (ROOT / "results" / "QC_REPORT.json").read_text(encoding="utf-8")
    )
    supplement = json.loads(
        (ROOT / "results" / "supplementary_tables_audit.json").read_text(
            encoding="utf-8"
        )
    )
    reference_audit = json.loads(
        (ROOT / "results" / "1114_reference_audit.json").read_text(
            encoding="utf-8"
        )
    )
    figure_audit = json.loads(
        (ROOT / "results" / "1114_figure_table_audit.json").read_text(
            encoding="utf-8"
        )
    )

    files = [
        ROOT / "1114manuscript.docx",
        ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD.docx",
        ROOT / "Supplementary_Tables.docx",
        ROOT / "manuscript" / "Supplementary_Tables.docx",
        ROOT / "Figure_1.drawio",
        ROOT / "Cover_Letter_Generic.docx",
        ROOT / "Cover_Letter_Generic.txt",
        ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx",
        ROOT / "results" / "QC_REPORT.json",
        ROOT / "results" / "supplementary_tables_audit.json",
        ROOT / "results" / "tables" / "supplementary_table_index.csv",
        ROOT / "IBD_AlphaGenome2_submission_package.zip",
    ]
    files = [path for path in files if path.exists()]
    archive_dir = ROOT / "_archive"
    archive_dir.mkdir(parents=True, exist_ok=True)
    archive_path = archive_dir / (
        f"IBD_AlphaGenome2_FINAL_{generated:%Y%m%d_%H%M%S}.zip"
    )
    with zipfile.ZipFile(
        archive_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
    ) as archive:
        for path in files:
            archive.write(path, path.relative_to(ROOT))

    archive_manifest = {
        "generated": generated.isoformat(),
        "archive": str(archive_path),
        "archive_sha256": sha256(archive_path),
        "files": [
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in files
        ],
    }
    archive_manifest_path = archive_path.with_suffix(".sha256.json")
    archive_manifest_path.write_text(
        json.dumps(archive_manifest, indent=2),
        encoding="utf-8",
    )

    delivery = {
        "generated": generated.isoformat(),
        "target_manuscript": "1114manuscript.docx",
        "target_supplement": "Supplementary_Tables.docx",
        "qc_checks": len(qc["checks"]),
        "qc_all_pass": qc["all_checks_pass"],
        "supplementary_tables": supplement["supplementary_tables"],
        "workbook_sheets_indexed": supplement["workbook_sheets_indexed"],
        "supplementary_all_checks_pass": supplement["all_checks_pass"],
        "references": reference_audit["references"],
        "each_reference_cited_once": not reference_audit[
            "repeated_reference_numbers"
        ],
        "pre_discussion_references": reference_audit[
            "pre_discussion_count"
        ],
        "discussion_references": reference_audit["discussion_count"],
        "first_citation_order_sequential": reference_audit[
            "first_citation_order_1_to_42"
        ],
        "all_citation_occurrences_hyperlinked": reference_audit[
            "all_citation_occurrences_hyperlinked"
        ],
        "long_sentences_over_30": sum(
            row["words"] > 30
            for row in reference_audit.get("long_sentences_over_30", [])
        ),
        "table_first_mention_order": figure_audit[
            "table_first_mention_order"
        ],
        "figure_first_mention_order": figure_audit[
            "figure_first_mention_order"
        ],
        "numeric_consistency": qc["checks"]["1114_numeric_consistency"],
        "archive": str(archive_path),
        "archive_sha256": archive_manifest["archive_sha256"],
        "remote_status": "local_complete_remote_pending",
        "files": archive_manifest["files"],
    }
    (ROOT / "results" / "FINAL_DELIVERY.json").write_text(
        json.dumps(delivery, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(delivery, indent=2))


if __name__ == "__main__":
    main()
