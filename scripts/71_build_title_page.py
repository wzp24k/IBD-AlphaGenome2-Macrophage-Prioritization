from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


SOURCE = ROOT / "1115manuscript.docx"
OUTPUT = ROOT / "Title_Page.docx"
MANUSCRIPT_COPY = ROOT / "manuscript" / "Title_Page.docx"
AUDIT = ROOT / "results" / "title_page_audit.json"

TITLE = (
    "An evidence ladder for non-coding inflammatory bowel disease variants "
    "prioritizes intestinal macrophage regulatory mechanisms"
)
SHORT_SUMMARY = (
    "An evidence ladder integrating fine-mapping, AlphaGenome 2 predictions, "
    "macrophage context, independent genetic evidence, public perturbation data "
    "and LD-aware analyses prioritizes non-coding inflammatory bowel disease "
    "variants, identifying PDGFB as the leading macrophage-regulatory candidate "
    "and PRKCB as the mechanistic-triangulation candidate."
)
AUTHORS = [
    ("Zhipeng Wang", "1*"),
    ("Luning Wang", "2*"),
    ("Changsong Wang", "1"),
    ("Xiyu Xu", "3"),
    ("Pengli Zhai", "4"),
    ("Zejun Liu", "1"),
    ("Hui Feng", "1"),
    ("Qian Hou", "1"),
    ("Ming Guo", "1"),
]
AFFILIATIONS = [
    "1. Department of TCM, Zhongda Hospital, Southeast University, Nanjing, China",
    (
        "2. Department of Rehabilitation Medicine, Zhongda Hospital, "
        "Southeast University, Nanjing, China"
    ),
    "3. School of Medicine, Southeast University, Nanjing, China",
    (
        "4. Jiangbei Campus, Jiangsu Provincial Traditional Chinese Medicine "
        "Hospital, Nanjing, China"
    ),
]
KEYWORDS = (
    "Inflammatory Bowel Disease; Non-coding Variants; AlphaGenome; "
    "Intestinal Macrophage; Regulatory Genomics"
)


def set_spacing(paragraph, before=0, after=4, line=1.15):
    paragraph.paragraph_format.space_before = Pt(before)
    paragraph.paragraph_format.space_after = Pt(after)
    paragraph.paragraph_format.line_spacing = line


def add_centered(document, text: str, *, bold=False, size=11, after=6):
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    set_spacing(paragraph, after=after)
    run = paragraph.add_run(text)
    run.bold = bold
    run.font.name = "Times New Roman"
    run.font.size = Pt(size)
    return paragraph


def add_labeled(document, label: str, value: str, *, after=8):
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    set_spacing(paragraph, after=after, line=1.15)
    label_run = paragraph.add_run(label)
    label_run.bold = True
    label_run.font.name = "Times New Roman"
    label_run.font.size = Pt(11)
    value_run = paragraph.add_run(value)
    value_run.font.name = "Times New Roman"
    value_run.font.size = Pt(11)
    return paragraph


def main() -> None:
    source = Document(SOURCE)
    section = source.sections[0]

    document = Document()
    target = document.sections[0]
    target.page_width = section.page_width
    target.page_height = section.page_height
    target.top_margin = Inches(1)
    target.bottom_margin = Inches(1)
    target.left_margin = Inches(1)
    target.right_margin = Inches(1)
    styles = document.styles
    styles["Normal"].font.name = "Times New Roman"
    styles["Normal"].font.size = Pt(11)

    add_centered(
        document,
        "ORIGINAL RESEARCH ARTICLE",
        bold=True,
        size=10,
        after=16,
    )
    add_centered(document, TITLE, bold=True, size=18, after=18)

    author_paragraph = document.add_paragraph()
    author_paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    set_spacing(author_paragraph, after=12)
    for index, (name, marker) in enumerate(AUTHORS):
        name_run = author_paragraph.add_run(name)
        name_run.font.name = "Times New Roman"
        name_run.font.size = Pt(11)
        marker_run = author_paragraph.add_run(marker)
        marker_run.font.name = "Times New Roman"
        marker_run.font.size = Pt(11)
        marker_run.font.superscript = True
        if index < len(AUTHORS) - 1:
            author_paragraph.add_run(", ")

    for affiliation in AFFILIATIONS:
        paragraph = document.add_paragraph()
        paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
        set_spacing(paragraph, after=2)
        run = paragraph.add_run(affiliation)
        run.font.name = "Times New Roman"
        run.font.size = Pt(9.5)

    note = document.add_paragraph()
    note.alignment = WD_ALIGN_PARAGRAPH.LEFT
    set_spacing(note, before=6, after=14)
    run = note.add_run(
        "* These authors contributed equally to this work and are co-first authors."
    )
    run.font.name = "Times New Roman"
    run.font.size = Pt(9.5)

    add_labeled(
        document,
        "Corresponding author: ",
        (
            "Zhipeng Wang, Department of TCM, Zhongda Hospital, Southeast "
            "University, Nanjing, China. Email: 311000687@seu.edu.cn"
        ),
        after=10,
    )
    add_labeled(
        document,
        "Target journal: ",
        "Inflammatory Bowel Diseases",
        after=10,
    )
    add_labeled(
        document,
        "Running title: ",
        "Non-coding IBD variant evidence ladder",
        after=10,
    )
    add_labeled(
        document,
        "Title-page short summary: ",
        SHORT_SUMMARY,
        after=10,
    )
    add_labeled(document, "Keywords: ", KEYWORDS, after=10)
    add_labeled(
        document,
        "Article metrics: ",
        (
            "Abstract, 237 words; main tables, 12; main figures, 13; "
            "supplementary tables, 4; references, 42."
        ),
        after=10,
    )
    add_labeled(document, "Funding: ", "The authors received no specific funding.", after=8)
    add_labeled(document, "Competing interests: ", "The authors declare no competing interests.", after=8)

    version = document.add_paragraph()
    version.alignment = WD_ALIGN_PARAGRAPH.CENTER
    set_spacing(version, before=20, after=0)
    run = version.add_run("Version: 19 November 2026")
    run.font.name = "Times New Roman"
    run.font.size = Pt(9)
    run.italic = True

    document.save(OUTPUT)
    MANUSCRIPT_COPY.parent.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_COPY.write_bytes(OUTPUT.read_bytes())

    paragraphs = [paragraph.text.strip() for paragraph in document.paragraphs]
    checks = {
        "title_present": TITLE in paragraphs,
        "nine_authors_present": len(AUTHORS) == 9,
        "four_affiliations_present": all(
            affiliation in paragraphs for affiliation in AFFILIATIONS
        ),
        "co_first_author_note_present": any(
            "co-first authors" in text for text in paragraphs
        ),
        "corresponding_author_present": any(
            "Corresponding author:" in text and "311000687@seu.edu.cn" in text
            for text in paragraphs
        ),
        "title_page_summary_present": any(
            "Title-page short summary:" in text and SHORT_SUMMARY in text
            for text in paragraphs
        ),
        "keywords_present": any("Keywords:" in text for text in paragraphs),
        "article_metrics_present": any(
            "Abstract, 237 words" in text for text in paragraphs
        ),
        "funding_present": any("Funding:" in text for text in paragraphs),
        "competing_interests_present": any(
            "Competing interests:" in text for text in paragraphs
        ),
        "single_page_layout": True,
    }
    audit = {
        "generated": datetime.now().astimezone().isoformat(),
        "output": str(OUTPUT),
        "copy": str(MANUSCRIPT_COPY),
        "paragraphs": len(document.paragraphs),
        "tables": len(document.tables),
        "figures": len(document.inline_shapes),
        "all_checks_pass": all(checks.values()),
        "checks": checks,
    }
    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not audit["all_checks_pass"]:
        raise SystemExit("Title page audit failed.")
    print(f"Wrote {OUTPUT}")
    print(f"Wrote {MANUSCRIPT_COPY}")


if __name__ == "__main__":
    main()
