from __future__ import annotations

import copy
import re
import sys
from collections import OrderedDict
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


TEMPLATE = ROOT / "1.docx"
SOURCE = ROOT / "manuscript" / "manuscript.md"
OUTPUT = ROOT / "manuscript" / "Manuscript_AlphaGenome2_IBD.docx"
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_NUMBER_MAP = {
    2: 1,
    4: 2,
    5: 3,
    1: 4,
    7: 5,
    6: 6,
    3: 7,
    8: 8,
    9: 9,
    10: 10,
    11: 11,
}
FIGURE_NUMBER_MAP = {
    1: 1,
    4: 2,
    2: 3,
    3: 4,
    6: 5,
    5: 6,
    7: 7,
    8: 8,
    9: 9,
    10: 10,
}


def renumber_crossrefs(text: str) -> str:
    text = re.sub(
        r"Table\s+(\d+)",
        lambda match: f"Table {TABLE_NUMBER_MAP.get(int(match.group(1)), int(match.group(1)))}",
        text,
    )
    text = re.sub(
        r"Figure\s+(\d+)",
        lambda match: f"Figure {FIGURE_NUMBER_MAP.get(int(match.group(1)), int(match.group(1)))}",
        text,
    )
    return text


def clean_inline(text: str) -> str:
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    return text.strip()


def add_internal_hyperlink(paragraph, text: str, anchor: str, number: int) -> None:
    hyperlink = OxmlElement("w:hyperlink")
    hyperlink.set(qn("w:anchor"), anchor)
    hyperlink.set(qn("w:history"), "1")
    hyperlink.set(qn("w:tooltip"), f"Jump to reference {number}")
    run = paragraph.add_run(text)
    run.font.color.rgb = RGBColor(0x34, 0x34, 0x34)
    run.underline = False
    paragraph._p.remove(run._r)
    hyperlink.append(run._r)
    paragraph._p.append(hyperlink)


def add_citation_links(paragraph, token: str) -> None:
    numbers = [int(value) for value in token[1:-1].split(",")]
    paragraph.add_run("[")
    for index, number in enumerate(numbers):
        if index:
            paragraph.add_run(",")
        add_internal_hyperlink(paragraph, str(number), f"ref_{number}", number)
    paragraph.add_run("]")


def add_inline(paragraph, text: str) -> None:
    text = renumber_crossrefs(clean_inline(text))
    pattern = re.compile(r"(\[[0-9][0-9,]*\]|\*\*[^*]+\*\*|\*[^*]+\*)")
    cursor = 0
    for match in pattern.finditer(text):
        if match.start() > cursor:
            paragraph.add_run(text[cursor : match.start()])
        token = match.group(0)
        if token.startswith("["):
            add_citation_links(paragraph, token)
        elif token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        else:
            run = paragraph.add_run(token[1:-1])
            run.italic = True
        cursor = match.end()
    if cursor < len(text):
        paragraph.add_run(text[cursor:])


def set_table_borders(table) -> None:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        element = OxmlElement(f"w:{edge}")
        if edge in {"top", "bottom"}:
            element.set(qn("w:val"), "single")
            element.set(qn("w:sz"), "12")
        else:
            element.set(qn("w:val"), "none")
            element.set(qn("w:sz"), "0")
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), "000000")
        borders.append(element)
    for cell in table.rows[0].cells:
        tc_pr = cell._tc.get_or_add_tcPr()
        tc_borders = tc_pr.first_child_found_in("w:tcBorders")
        if tc_borders is None:
            tc_borders = OxmlElement("w:tcBorders")
            tc_pr.append(tc_borders)
        bottom = OxmlElement("w:bottom")
        bottom.set(qn("w:val"), "single")
        bottom.set(qn("w:sz"), "12")
        bottom.set(qn("w:space"), "0")
        bottom.set(qn("w:color"), "000000")
        tc_borders.append(bottom)


def add_table_caption(document, caption: str) -> None:
    caption = renumber_crossrefs(caption)
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.keep_with_next = True
    match = re.match(r"(Table\s+\d+\.)\s*(.*)", caption, re.IGNORECASE)
    if match:
        run = paragraph.add_run(match.group(1) + " ")
        run.bold = True
        add_inline(paragraph, match.group(2))
    else:
        add_inline(paragraph, caption)


def add_markdown_table(document, rows: list[list[str]], width: int = 1) -> None:
    if not rows:
        return
    n_cols = max(len(row) for row in rows)
    rows = [row + [""] * (n_cols - len(row)) for row in rows]
    table = document.add_table(rows=len(rows), cols=n_cols)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    table.style = "Normal Table"
    for row_index, row in enumerate(rows):
        for col_index, value in enumerate(row):
            cell = table.cell(row_index, col_index)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.space_before = Pt(2)
            paragraph.paragraph_format.space_after = Pt(2)
            paragraph.paragraph_format.line_spacing = 1.0
            paragraph.alignment = (
                WD_ALIGN_PARAGRAPH.CENTER if row_index == 0 else WD_ALIGN_PARAGRAPH.LEFT
            )
            run = paragraph.add_run(clean_inline(value))
            run.font.name = "Times New Roman"
            run.font.size = Pt(9)
            if row_index == 0:
                run.bold = True
    set_table_borders(table)
    document.add_paragraph()


def add_figure(document, path: Path, caption: str, width: float = 6.4) -> None:
    if not path.exists():
        return
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.add_run().add_picture(str(path), width=Inches(width))
    caption = renumber_crossrefs(caption)
    cap = document.add_paragraph()
    cap.paragraph_format.keep_with_next = False
    match = re.match(r"(Figure\s+\d+\.)\s*(.*)", caption, re.IGNORECASE)
    if match:
        run = cap.add_run(match.group(1) + " ")
        run.bold = True
        add_inline(cap, match.group(2))
    else:
        add_inline(cap, caption)


def parse_markdown(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    title = ""
    sections: OrderedDict[str, dict] = OrderedDict()
    current_section = None
    current_sub = None
    for line in lines:
        if line.startswith("# ") and not title:
            title = line[2:].strip()
            continue
        if line.startswith("## "):
            current_section = line[3:].strip()
            sections[current_section] = {"intro": [], "subs": OrderedDict()}
            current_sub = None
            continue
        if line.startswith("### "):
            current_sub = line[4:].strip()
            sections[current_section]["subs"][current_sub] = []
            continue
        if current_section is None:
            continue
        if current_sub is None:
            sections[current_section]["intro"].append(line)
        else:
            sections[current_section]["subs"][current_sub].append(line)
    return title, sections


def parse_table_block(lines: list[str]):
    rows = []
    notes = []
    in_table = False
    for line in lines:
        if line.strip().startswith("|"):
            cells = [clean_inline(cell.strip()) for cell in line.strip().strip("|").split("|")]
            if all(re.fullmatch(r":?-{2,}:?", cell) for cell in cells):
                continue
            rows.append(cells)
            in_table = True
        elif in_table and line.strip():
            notes.append(line.strip())
    return rows, notes


def add_body_paragraph(document, text: str, *, note: bool = False) -> None:
    if not text.strip():
        return
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.paragraph_format.line_spacing = 1.5
    if note:
        for run in paragraph.runs:
            run.font.size = Pt(10)
    add_inline(paragraph, text)
    if note:
        for run in paragraph.runs:
            run.font.size = Pt(10)


def add_section_heading(document, text: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.keep_with_next = True
    run = paragraph.add_run(text)
    run.bold = True
    run.font.size = Pt(14)


def add_subheading(document, text: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.keep_with_next = True
    run = paragraph.add_run(text)
    run.bold = True


def add_reference_bookmark(paragraph, number: int, bookmark_id: int) -> None:
    bookmark_start = OxmlElement("w:bookmarkStart")
    bookmark_start.set(qn("w:id"), str(bookmark_id))
    bookmark_start.set(qn("w:name"), f"ref_{number}")
    bookmark_end = OxmlElement("w:bookmarkEnd")
    bookmark_end.set(qn("w:id"), str(bookmark_id))
    paragraph._p.insert(0, bookmark_start)
    paragraph._p.append(bookmark_end)


def main() -> None:
    title, sections = parse_markdown(SOURCE)
    document = Document(TEMPLATE)
    author_block_xml = [
        copy.deepcopy(paragraph._p) for paragraph in document.paragraphs[1:8]
    ]
    body = document._element.body
    for child in list(body):
        if child.tag != qn("w:sectPr"):
            body.remove(child)

    title_p = document.add_paragraph()
    title_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title_run = title_p.add_run(title)
    title_run.bold = True
    title_run.font.size = Pt(18)
    for author_paragraph in author_block_xml:
        body.insert(len(body) - 1, author_paragraph)

    abstract = sections.get("Abstract", {"intro": [], "subs": {}})
    add_section_heading(document, "Abstract")
    keywords = ""
    for line in abstract["intro"]:
        if line.startswith("**Keywords:**"):
            keywords = line
        elif line.strip():
            add_body_paragraph(document, line)
    if keywords:
        add_body_paragraph(document, keywords)

    section_titles = {
        "Introduction": "1. Introduction",
        "Materials and Methods": "2. Materials and Methods",
        "Results": "3. Results",
        "Discussion": "4. Discussion",
    }
    section_lookup = {
        "Introduction": "Introduction",
        "Materials and Methods": "Methods",
        "Results": "Results",
        "Discussion": "Discussion",
    }

    table_source = sections.get("Tables", {"intro": [], "subs": OrderedDict()})["subs"]
    tables = {}
    for caption, lines in table_source.items():
        rows, notes = parse_table_block(lines)
        number = int(re.search(r"Table\s+(\d+)", caption).group(1))
        tables[number] = (caption, rows, notes)

    figure_source = sections.get("Figure legends", {"intro": [], "subs": OrderedDict()})["intro"]
    figures = {}
    for line in figure_source:
        match = re.match(r"\*\*Figure\s+(\d+)\.\*\*\s*(.*)", line)
        if match:
            figures[int(match.group(1))] = f"Figure {match.group(1)}. {match.group(2)}"

    table_map = {
        1: "Tiered macrophage-priority map",
        2: "Study design and credible-set landscape",
        3: "External IBD genetic validation",
        4: "Controls establish cell-context behavior",
        5: "Controls establish cell-context behavior",
        6: "Threshold, weighting and background sensitivity",
        7: "Donor-derived and activity-state evaluation",
        8: "Public experimental evidence",
        9: "Macrophage Perturb-seq and robustness",
        10: "Gene-expression programs and confounding",
        11: "Evidence ceiling and claim calibration",
    }
    figure_map = {
        1: "Study design and credible-set landscape",
        2: "Atlas calibration and live specificity are complementary",
        3: "Tiered macrophage-priority map",
        4: "Controls establish cell-context behavior",
        5: "Threshold, weighting and background sensitivity",
        6: "Donor-derived and activity-state evaluation",
        7: "External IBD genetic validation",
        8: "Public experimental evidence",
        9: "Macrophage Perturb-seq and robustness",
        10: "Gene-expression programs and confounding",
    }
    subsection_count = {key: 0 for key in section_lookup}
    rendered_tables = set()
    rendered_figures = set()
    conclusion_text = ""

    def render_assets_for_subsection(sub_name: str, section_key: str) -> None:
        for number, mapped in table_map.items():
            if mapped == sub_name and number not in rendered_tables and number in tables:
                caption, rows, notes = tables[number]
                add_table_caption(document, caption)
                add_markdown_table(document, rows)
                for note in notes:
                    add_body_paragraph(document, note, note=True)
                rendered_tables.add(number)
        for number, mapped in figure_map.items():
            if mapped != sub_name or number in rendered_figures or number not in figures:
                continue
            width = 5.4 if number == 3 else 6.4
            add_figure(
                document,
                FIGURE_DIR / {
                    1: "figure1_analysis_flow.png",
                    2: "figure2_quantile_vs_specificity.png",
                    3: "figure3_top_candidate_heatmap.png",
                    4: "figure4_ets2_positive_control.png",
                    5: "figure5_ets2_background_null.png",
                    6: "figure6_active_ibd_validation.png",
                    7: "figure7_all_l2g_pathways.png",
                    8: "figure8_external_experimental_evidence.png",
                    9: "figure9_perturbseq_validation.png",
                    10: "figure10_evidence_robustness.png",
                }[number],
                figures[number],
                width=width,
            )
            rendered_figures.add(number)

    for section_name, section_heading in section_titles.items():
        section = sections.get(section_lookup[section_name])
        if not section:
            continue
        add_section_heading(document, section_heading)
        for line in section["intro"]:
            if section_name == "Discussion" and line.strip().startswith(
                "In conclusion,"
            ):
                conclusion_text = line.strip()
            elif line.strip():
                add_body_paragraph(document, line)
        for sub_name, lines in section["subs"].items():
            subsection_count[section_name] += 1
            numbered_sub = f"{section_heading.split('.')[0]}.{subsection_count[section_name]} {sub_name}"
            add_subheading(document, numbered_sub)
            for line in lines:
                if section_name == "Discussion" and line.strip().startswith(
                    "In conclusion,"
                ):
                    conclusion_text = line.strip()
                elif line.strip():
                    add_body_paragraph(document, line)
            render_assets_for_subsection(sub_name, section_name)
    add_section_heading(document, "5. Conclusion")
    if conclusion_text:
        add_body_paragraph(document, conclusion_text)

    add_section_heading(document, "6. Declarations")
    declaration_sections = [
        ("Data availability", "Data availability"),
        ("Code availability", "Code availability"),
        (
            "Ethics approval and consent to participate",
            "Ethics approval and consent to participate",
        ),
        ("Patient and public involvement", "Patient and public involvement"),
        ("Authors' contributions", "Author contributions"),
        ("Funding", "Funding"),
        (
            "Declaration of generative AI and AI-assisted technologies in the writing process",
            "Declaration of generative AI and AI-assisted technologies in the writing process",
        ),
        ("Consent for publication", "Consent for publication"),
        ("Clinical trial number", "Clinical trial number"),
        ("Competing interests", "Competing interests"),
    ]
    for display_title, source_key in declaration_sections:
        add_subheading(document, display_title)
        for line in sections.get(source_key, {"intro": []})["intro"]:
            if line.strip():
                add_body_paragraph(document, line)

    add_section_heading(document, "7. References")
    references = sections.get("References", {"intro": []})["intro"]
    existing_bookmark_ids = []
    for bookmark in document.element.body.iter(qn("w:bookmarkStart")):
        value = bookmark.get(qn("w:id"))
        if value and value.isdigit():
            existing_bookmark_ids.append(int(value))
    next_bookmark_id = max(existing_bookmark_ids, default=0) + 1
    for line in references:
        if line.strip():
            paragraph = document.add_paragraph(style="RefEntry")
            add_inline(paragraph, line.strip())
            reference_match = re.match(r"^(\d+)\.\s+", line.strip())
            if reference_match:
                add_reference_bookmark(
                    paragraph,
                    int(reference_match.group(1)),
                    next_bookmark_id,
                )
                next_bookmark_id += 1

    # Ensure any tables/figures not anchored to a detected subsection are appended
    # before the declarations, preserving captions and notes.
    if len(rendered_tables) < len(tables):
        add_section_heading(document, "Supplementary display items")
        for number, (caption, rows, notes) in sorted(tables.items()):
            if number in rendered_tables:
                continue
            add_table_caption(document, caption)
            add_markdown_table(document, rows)
            for note in notes:
                add_body_paragraph(document, note, note=True)

    document.save(OUTPUT)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
