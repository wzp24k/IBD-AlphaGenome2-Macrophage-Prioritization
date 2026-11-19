from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.shared import Pt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


SOURCE = ROOT / "STROBE MR Checklist.docx"
OUTPUT = ROOT / "STROBE-MR_Checklist_1119.docx"
COPY = ROOT / "manuscript" / "STROBE-MR_Checklist_1119.docx"
AUDIT = ROOT / "results" / "strobe_mr_1119_audit.json"


def replace_paragraph(paragraph, text: str) -> None:
    first = paragraph.runs[0] if paragraph.runs else None
    size = first.font.size if first else Pt(9)
    bold = first.bold if first else None
    for child in list(paragraph._p):
        if child.tag != qn("w:pPr"):
            paragraph._p.remove(child)
    run = paragraph.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = size
    if bold is not None:
        run.bold = bold


def set_cell(cell, text: str) -> None:
    cell.text = text
    for paragraph in cell.paragraphs:
        for run in paragraph.runs:
            run.font.name = "Times New Roman"
            run.font.size = Pt(7)


def main() -> None:
    shutil.copy2(SOURCE, OUTPUT)
    document = Document(OUTPUT)
    replace_paragraph(
        document.paragraphs[0],
        "Page references follow the 41-page 1119manuscript.docx; the extended "
        "revision adds cross-ancestry, pQTL and targeted PheW-MR analyses.",
    )
    replace_paragraph(
        document.paragraphs[1],
        "MR component covered: two-sample summary-data MR of genetically proxied "
        "monocyte PRKCB expression and IBD risk, including cross-ancestry, pQTL "
        "and targeted PheW-MR extensions.",
    )
    table = document.tables[0]
    set_cell(
        table.rows[13].cells[4],
        "Addressed. The main analysis used inverse-variance weighting, weighted "
        "median and MR-Egger. Cross-ancestry PRKCB MR was repeated in East Asian "
        "IBD, and pQTL cis-MR used SomaScan and curated deCODE instruments. A "
        "targeted PheW-MR panel provided an initial safety screen. A Steiger-style "
        "check assessed direction.",
    )
    set_cell(
        table.rows[32].cells[4],
        "Addressed. Combined-IBD IVW beta was 0.0511 log odds per SD "
        "(P = 1.07e-17). Crohn beta was 0.0351 (P = 1.91e-5), UC beta was 0.0661 "
        "(P = 7.83e-19), and East Asian IBD beta was 0.0484 (P = 4.09e-4). "
        "SomaScan PKCB2 pQTL MR was positive for combined IBD and UC. Curated "
        "deCODE pQTL MR was positive for East Asian IBD.",
    )
    set_cell(
        table.rows[51].cells[4],
        "Addressed. Findings cover European/multi-ancestry IBD and an East Asian "
        "IBD replication for PRKCB. pQTL MR is exploratory because available "
        "protein instruments have limited regional colocalization support. "
        "Generalization to other myeloid tissues, disease states and exposure "
        "timings will require further cohorts.",
    )
    document.save(OUTPUT)
    COPY.parent.mkdir(parents=True, exist_ok=True)
    COPY.write_bytes(OUTPUT.read_bytes())
    revised = Document(OUTPUT)
    text = "\n".join(paragraph.text for paragraph in revised.paragraphs)
    table_text = "\n".join(
        cell.text
        for row in revised.tables[0].rows
        for cell in row.cells
    )
    checks = {
        "mentions_1119": "1119manuscript.docx" in text,
        "mentions_cross_ancestry": "Cross-ancestry" in table_text,
        "mentions_pqtl": "pQTL" in table_text,
        "mentions_phewas": "PheW-MR" in table_text,
        "same_template_table": len(revised.tables) == 1,
        "all_20_items_present": all(
            str(number) in [row.cells[0].text.strip() for row in revised.tables[0].rows]
            for number in range(1, 21)
        ),
    }
    audit = {
        "output": str(OUTPUT),
        "copy": str(COPY),
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("STROBE-MR 1119 update failed.")


if __name__ == "__main__":
    main()
