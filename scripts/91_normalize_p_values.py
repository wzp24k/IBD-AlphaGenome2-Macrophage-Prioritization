from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from docx import Document

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1119manuscript.docx"
COPY = ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD_clinical_translation.docx"
AUDIT = ROOT / "results" / "p_value_format_audit.json"

PATTERN = re.compile(
    r"\bP\s*(=|≤|<)\s*"
    r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)"
)


def normalize(match: re.Match) -> str:
    operator = match.group(1)
    value = float(match.group(2))
    if value < 0.001 or (operator in {"<", "≤"} and value <= 0.001):
        return "P < 0.001"
    return f"P = {value:.3f}"


def replace_in_paragraph(paragraph) -> int:
    count = 0
    for run in paragraph.runs:
        if not run.text:
            continue
        new_text, replacements = PATTERN.subn(normalize, run.text)
        if replacements:
            run.text = new_text
            count += replacements
    return count


def main() -> None:
    document = Document(DOCX)
    count = 0
    for paragraph in document.paragraphs:
        count += replace_in_paragraph(paragraph)
    for table in document.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    count += replace_in_paragraph(paragraph)

    document.save(DOCX)
    COPY.write_bytes(DOCX.read_bytes())

    revised = Document(DOCX)
    remaining = []
    for paragraph in revised.paragraphs:
        for match in PATTERN.finditer(paragraph.text):
            canonical = normalize(match)
            if canonical != match.group(0):
                remaining.append(
                    {
                        "text": match.group(0),
                        "expected": canonical,
                        "paragraph": paragraph.text[:160],
                    }
                )
    audit = {
        "replacements": count,
        "remaining_noncanonical": remaining,
        "rule": "P >= 0.001 uses three decimals; P < 0.001 is written P < 0.001",
        "all_checks_pass": not remaining,
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("P-value formatting audit failed.")


if __name__ == "__main__":
    main()
