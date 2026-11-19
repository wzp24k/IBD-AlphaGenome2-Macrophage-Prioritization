from __future__ import annotations

import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path

from docx import Document
from openpyxl import load_workbook

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1119manuscript.docx"
WORKBOOK = ROOT / "results" / "IBD_AlphaGenome2_supplementary_tables.xlsx"
FIGURE = ROOT / "results" / "figures" / "figure13_clinical_translation.png"
OUTPUT = ROOT / "results" / "clinical_translation_extension_audit.json"


def main() -> None:
    document = Document(DOCX)
    paragraphs = [paragraph.text.strip() for paragraph in document.paragraphs]
    nonempty = [text for text in paragraphs if text]
    discussion_index = nonempty.index("4. Discussion")
    body = "\n".join(nonempty[:discussion_index])
    table_mentions = sorted(set(map(int, re.findall(r"\bTable\s+(\d+)", body))))
    figure_mentions = sorted(set(map(int, re.findall(r"\bFigure\s+(\d+)", body))))
    table_captions = sorted(
        {
            int(match.group(1))
            for text in nonempty
            if (match := re.match(r"^Table\s+(\d+)\.", text))
        }
    )
    figure_captions = sorted(
        {
            int(match.group(1))
            for text in nonempty
            if (match := re.match(r"^Figure\s+(\d+)\.", text))
        }
    )
    with zipfile.ZipFile(DOCX) as archive:
        media_hashes = {
            hashlib.sha256(archive.read(name)).hexdigest()
            for name in archive.namelist()
            if name.startswith("word/media/")
        }
    figure_hash = hashlib.sha256(FIGURE.read_bytes()).hexdigest()
    workbook = load_workbook(WORKBOOK, read_only=True)
    expected_sheets = {
        "CD_UC_target_support",
        "PRKCB_cross_ancestry_MR",
        "pQTL_cis_MR",
        "pQTL_colocalization",
        "PheW_MR",
        "Treatment_response_programs",
        "Treatment_response_meta",
        "Clinical_readiness",
    }
    checks = {
        "twelve_tables": len(document.tables) == 12,
        "thirteen_figures": len(document.inline_shapes) == 13,
        "table_mentions_1_to_12": table_mentions == list(range(1, 13)),
        "figure_mentions_1_to_13": figure_mentions == list(range(1, 14)),
        "table_captions_1_to_12": table_captions == list(range(1, 13)),
        "figure_captions_1_to_13": figure_captions == list(range(1, 14)),
        "figure13_embedded": figure_hash in media_hashes,
        "methods_extension_present": "2.10 Clinical translation extension" in nonempty,
        "results_extension_present": "3.15 Clinical translation extension" in nonempty,
        "table12_positive": "Clinical development readiness matrix" in body,
        "workbook_extensions_present": expected_sheets.issubset(set(workbook.sheetnames)),
    }
    report = {
        "document": str(DOCX),
        "tables": len(document.tables),
        "figures": len(document.inline_shapes),
        "table_mentions": table_mentions,
        "figure_mentions": figure_mentions,
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    OUTPUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["all_checks_pass"]:
        raise SystemExit("Clinical translation extension audit failed.")


if __name__ == "__main__":
    main()
