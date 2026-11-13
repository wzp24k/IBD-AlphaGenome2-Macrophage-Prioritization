from __future__ import annotations

import json
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


def is_three_line_table(table) -> bool:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        return False
    values = {}
    for edge in ["top", "bottom", "insideH", "insideV"]:
        element = borders.find(qn(f"w:{edge}"))
        if element is None:
            return False
        values[edge] = element.get(qn("w:val"))
    return (
        values["top"] == "single"
        and values["bottom"] == "single"
        and values["insideH"] == "none"
        and values["insideV"] == "none"
    )


def hyperlink_groups(document_xml, namespace) -> list[list[int]]:
    groups = []
    for paragraph in document_xml.findall(".//w:p", namespace):
        current = None
        for child in paragraph:
            if child.tag == f"{{{namespace['w']}}}r":
                text = "".join(
                    node.text or "" for node in child.findall(".//w:t", namespace)
                )
                if text == "[":
                    current = []
                elif text == "]" and current is not None:
                    groups.append(current)
                    current = None
            elif (
                child.tag == f"{{{namespace['w']}}}hyperlink"
                and current is not None
            ):
                anchor = child.get(f"{{{namespace['w']}}}anchor", "")
                if anchor.startswith("ref_"):
                    current.append(int(anchor.removeprefix("ref_")))
    return groups


def main() -> None:
    template = Document(ROOT / "1.docx")
    manuscript = Document(
        ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD.docx"
    )
    template_section = template.sections[0]
    manuscript_section = manuscript.sections[0]
    checks = {
        "page_width_matches": manuscript_section.page_width
        == template_section.page_width,
        "page_height_matches": manuscript_section.page_height
        == template_section.page_height,
        "top_margin_matches": manuscript_section.top_margin
        == template_section.top_margin,
        "bottom_margin_matches": manuscript_section.bottom_margin
        == template_section.bottom_margin,
        "left_margin_matches": manuscript_section.left_margin
        == template_section.left_margin,
        "right_margin_matches": manuscript_section.right_margin
        == template_section.right_margin,
        "body_font_matches": manuscript.styles["Normal"].font.name
        == template.styles["Normal"].font.name,
        "body_size_matches": manuscript.styles["Normal"].font.size
        == template.styles["Normal"].font.size,
        "body_line_spacing_matches": manuscript.styles["Normal"]
        .paragraph_format.line_spacing
        == template.styles["Normal"].paragraph_format.line_spacing,
        "eleven_tables": len(manuscript.tables) == 11,
        "ten_figures": len(manuscript.inline_shapes) == 10,
        "letter_page_size": (
            manuscript_section.page_width == 7772400
            and manuscript_section.page_height == 10058400
        ),
        "no_header": all(
            not paragraph.text.strip()
            for paragraph in manuscript_section.header.paragraphs
        ),
        "no_footer": all(
            not paragraph.text.strip()
            for paragraph in manuscript_section.footer.paragraphs
        ),
    }
    title_paragraph = manuscript.paragraphs[0]
    checks["title_18pt_bold"] = bool(title_paragraph.runs) and all(
        run.bold and run.font.size and run.font.size.pt == 18
        for run in title_paragraph.runs
    )
    section_headings = [
        paragraph.text
        for paragraph in manuscript.paragraphs
        if paragraph.runs
        and paragraph.runs[0].bold
        and paragraph.runs[0].font.size
        and paragraph.runs[0].font.size.pt == 14
    ]
    expected_headings = [
        "Abstract",
        "1. Introduction",
        "2. Materials and Methods",
        "3. Results",
        "4. Discussion",
        "5. Conclusion",
        "6. Declarations",
        "7. References",
    ]
    checks["headings_match_template"] = section_headings == expected_headings
    table_fonts = []
    for table in manuscript.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    table_fonts.extend(
                        run.font.size.pt
                        for run in paragraph.runs
                        if run.font.size is not None
                    )
    checks["table_font_9pt"] = bool(table_fonts) and all(
        size == 9 for size in table_fonts
    )
    checks["three_line_tables"] = bool(manuscript.tables) and all(
        is_three_line_table(table) for table in manuscript.tables
    )
    with zipfile.ZipFile(
        ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD.docx"
    ) as archive:
        document_xml = ET.fromstring(archive.read("word/document.xml"))
    namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    bookmarks = {
        element.get(f"{{{namespace['w']}}}name")
        for element in document_xml.findall(".//w:bookmarkStart", namespace)
        if element.get(f"{{{namespace['w']}}}name", "").startswith("ref_")
    }
    hyperlinks = document_xml.findall(".//w:hyperlink", namespace)
    reference_hyperlinks = [
        element
        for element in hyperlinks
        if element.get(f"{{{namespace['w']}}}anchor", "").startswith("ref_")
    ]
    linked_numbers = []
    valid_link_text = True
    for hyperlink in reference_hyperlinks:
        anchor = hyperlink.get(f"{{{namespace['w']}}}anchor", "")
        text = "".join(
            node.text or ""
            for node in hyperlink.findall(".//w:t", namespace)
        )
        number = int(anchor.removeprefix("ref_"))
        linked_numbers.append(number)
        valid_link_text = valid_link_text and text == str(number)
    expected_numbers = set(range(1, 43))
    checks["forty_two_reference_bookmarks"] = bookmarks == {
        f"ref_{number}" for number in expected_numbers
    }
    checks["all_references_hyperlink_targets"] = set(linked_numbers) == expected_numbers
    checks["citation_hyperlink_text_matches_anchor"] = valid_link_text
    checks["citation_hyperlinks_sorted_within_group"] = all(
        group
        and len(group) == len(set(group))
        and group == sorted(group)
        for group in hyperlink_groups(document_xml, namespace)
    )
    report = {
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "section_headings": section_headings,
        "counts": {
            "paragraphs": len(manuscript.paragraphs),
            "tables": len(manuscript.tables),
            "figures": len(manuscript.inline_shapes),
            "reference_bookmarks": len(bookmarks),
            "reference_hyperlinks": len(reference_hyperlinks),
        },
        "page": {
            "width": str(manuscript_section.page_width),
            "height": str(manuscript_section.page_height),
            "margins": {
                "top": str(manuscript_section.top_margin),
                "bottom": str(manuscript_section.bottom_margin),
                "left": str(manuscript_section.left_margin),
                "right": str(manuscript_section.right_margin),
            },
        },
    }
    expected_author = (
        "Zhipeng Wang1*, Luning Wang2*, Changsong Wang1, Xiyu Xu3, "
        "Pengli Zhai4, Zejun Liu1, Hui Feng1, Qian Hou1, Ming Guo1"
    )
    manuscript_text = [paragraph.text for paragraph in manuscript.paragraphs[:8]]
    checks["author_block_matches_template"] = (
        expected_author in manuscript_text
        and any("311000687@seu.edu.cn" in text for text in manuscript_text)
        and any(
            "These authors contributed equally" in text
            for text in manuscript_text
        )
    )
    report["all_checks_pass"] = all(checks.values())
    report["author_block"] = manuscript_text[1:8]
    output = ROOT / "results" / "template_docx_qc.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["all_checks_pass"]:
        raise SystemExit("Template DOCX QC failed.")


if __name__ == "__main__":
    main()
