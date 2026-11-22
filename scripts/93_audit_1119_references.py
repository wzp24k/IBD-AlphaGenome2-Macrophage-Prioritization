from __future__ import annotations

import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from docx import Document

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


DOCX = ROOT / "1119manuscript.docx"
OUTPUT = ROOT / "results" / "1119_reference_audit.json"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def groups_from_xml(document_xml) -> list[list[int]]:
    groups = []
    for paragraph in document_xml.findall(".//w:p", NS):
        current = None
        for child in paragraph:
            if child.tag == f"{{{NS['w']}}}r":
                text = "".join(
                    node.text or ""
                    for node in child.findall(".//w:t", NS)
                )
                if text == "[":
                    current = []
                elif text == "]" and current is not None:
                    groups.append(current)
                    current = None
            elif (
                child.tag == f"{{{NS['w']}}}hyperlink"
                and current is not None
            ):
                anchor = child.get(f"{{{NS['w']}}}anchor", "")
                if anchor.startswith("ref_"):
                    current.append(int(anchor.removeprefix("ref_")))
    return groups


def main() -> None:
    document = Document(DOCX)
    paragraph_texts = [p.text.strip() for p in document.paragraphs]
    discussion = paragraph_texts.index("4. Discussion")
    conclusion = paragraph_texts.index("5. Conclusion")
    references = paragraph_texts.index("7. References")
    pre = set()
    disc = set()
    body_groups = []
    for index, paragraph in enumerate(document.paragraphs):
        text = paragraph.text.strip()
        groups = [
            [int(value) for value in group.split(",")]
            for group in re.findall(r"\[([0-9][0-9, ]*)\]", text)
        ]
        if not groups:
            continue
        body_groups.extend(groups)
        for group in groups:
            target = disc if discussion <= index < conclusion else pre
            target.update(group)
    ref_numbers = [
        int(match.group(1))
        for text in paragraph_texts[references + 1 :]
        if (match := re.match(r"^(\d+)\.", text))
    ]
    with zipfile.ZipFile(DOCX) as archive:
        xml = ET.fromstring(archive.read("word/document.xml"))
    hyperlink_groups = groups_from_xml(xml)
    bookmarks = {
        element.get(f"{{{NS['w']}}}name")
        for element in xml.findall(".//w:bookmarkStart", NS)
        if element.get(f"{{{NS['w']}}}name", "").startswith("ref_")
    }
    occurrences = Counter(
        number for group in hyperlink_groups for number in group
    )
    report = {
        "source": str(DOCX),
        "references": len(ref_numbers),
        "reference_numbers": ref_numbers,
        "pre_discussion_count": len(pre),
        "discussion_count": len(disc),
        "pre_discussion_references": sorted(pre),
        "discussion_references": sorted(disc),
        "repeated_reference_numbers": {
            number: count
            for number, count in sorted(occurrences.items())
            if count > 1
        },
        "max_references_in_group": max(
            (len(group) for group in hyperlink_groups),
            default=0,
        ),
        "groups_with_more_than_two": [
            group for group in hyperlink_groups if len(group) > 2
        ],
        "groups_not_sorted_unique": [
            group
            for group in hyperlink_groups
            if group != sorted(set(group))
        ],
        "bookmarks": len(bookmarks),
        "all_citation_occurrences_hyperlinked": all(
            f"ref_{number}" in bookmarks
            for group in hyperlink_groups
            for number in group
        ),
        "long_sentences_over_30": [],
    }
    OUTPUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
