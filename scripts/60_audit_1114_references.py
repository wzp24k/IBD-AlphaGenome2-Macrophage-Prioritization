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


SOURCE = ROOT / "1114manuscript.docx"
OUTPUT_JSON = ROOT / "results" / "1114_reference_audit.json"
OUTPUT_TEXT = ROOT / "work" / "1114_manuscript_extracted.txt"
NAMESPACE = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
}


def paragraph_citations(text: str) -> list[list[int]]:
    groups = []
    for token in re.findall(r"\[([0-9][0-9,\- ]*)\]", text):
        numbers = []
        for part in token.split(","):
            part = part.strip()
            if "-" in part:
                start, end = map(int, part.split("-", 1))
                numbers.extend(range(start, end + 1))
            elif part:
                numbers.append(int(part))
        groups.append(numbers)
    return groups


def xml_citation_groups(document_xml) -> list[list[int]]:
    groups = []
    for paragraph in document_xml.findall(".//w:p", NAMESPACE):
        current = None
        for child in paragraph:
            if child.tag == f"{{{NAMESPACE['w']}}}r":
                text = "".join(
                    node.text or ""
                    for node in child.findall(".//w:t", NAMESPACE)
                )
                if text == "[":
                    current = []
                elif text == "]" and current is not None:
                    groups.append(current)
                    current = None
            elif (
                child.tag == f"{{{NAMESPACE['w']}}}hyperlink"
                and current is not None
            ):
                anchor = child.get(f"{{{NAMESPACE['w']}}}anchor", "")
                if anchor.startswith("ref_"):
                    current.append(int(anchor.removeprefix("ref_")))
    return groups


def sentence_lengths(text: str) -> list[dict]:
    clean = re.sub(r"\[[0-9,\- ]+\]", "", text)
    clean = re.sub(r"[*`]", "", clean)
    rows = []
    for sentence in re.split(r"(?<=[.!?])\s+", clean):
        count = len(
            re.findall(r"[A-Za-z0-9][A-Za-z0-9'./-]*", sentence)
        )
        if count:
            rows.append({"words": count, "sentence": sentence.strip()})
    return rows


def main() -> None:
    document = Document(SOURCE)
    paragraphs = [paragraph.text.strip() for paragraph in document.paragraphs]
    nonempty = [
        (index, text)
        for index, text in enumerate(paragraphs)
        if text
    ]
    text = "\n".join(paragraphs)
    OUTPUT_TEXT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_TEXT.write_text(text, encoding="utf-8")

    discussion_index = next(
        index
        for index, text in nonempty
        if text == "4. Discussion"
    )
    conclusion_index = next(
        index
        for index, text in nonempty
        if text == "5. Conclusion"
    )
    references_index = next(
        index
        for index, text in nonempty
        if text == "7. References"
    )

    body_rows = []
    for index, text in nonempty:
        if index >= references_index:
            continue
        groups = paragraph_citations(text)
        if groups:
            body_rows.append(
                {
                    "paragraph_index": index,
                    "section": (
                        "discussion"
                        if discussion_index <= index < conclusion_index
                        else (
                            "conclusion"
                            if conclusion_index <= index < references_index
                            else "pre_discussion"
                        )
                    ),
                    "groups": groups,
                    "text": text,
                }
            )

    references = []
    for index, text in nonempty:
        if index <= references_index:
            continue
        match = re.match(r"^(\d+)\.\s+", text)
        if match:
            references.append(
                {
                    "number": int(match.group(1)),
                    "paragraph_index": index,
                    "text": text,
                }
            )

    occurrences = Counter()
    paragraph_occurrences = []
    for row in body_rows:
        local = Counter(
            number
            for group in row["groups"]
            for number in group
        )
        occurrences.update(local)
        paragraph_occurrences.append(
            {
                "paragraph_index": row["paragraph_index"],
                "citations": dict(sorted(local.items())),
            }
        )

    with zipfile.ZipFile(SOURCE) as archive:
        xml = ET.fromstring(archive.read("word/document.xml"))
    bookmarks = {
        element.get(f"{{{NAMESPACE['w']}}}name")
        for element in xml.findall(".//w:bookmarkStart", NAMESPACE)
        if element.get(f"{{{NAMESPACE['w']}}}name", "").startswith("ref_")
    }
    hyperlinks = xml.findall(".//w:hyperlink", NAMESPACE)
    ref_hyperlinks = [
        element
        for element in hyperlinks
        if element.get(f"{{{NAMESPACE['w']}}}anchor", "").startswith("ref_")
    ]
    hyperlink_numbers = [
        int(
            element.get(f"{{{NAMESPACE['w']}}}anchor", "").removeprefix(
                "ref_"
            )
        )
        for element in ref_hyperlinks
    ]
    hyperlink_groups = xml_citation_groups(xml)

    pre_numbers = sorted(
        {
            number
            for row in body_rows
            if row["section"] == "pre_discussion"
            for group in row["groups"]
            for number in group
        }
    )
    discussion_numbers = sorted(
        {
            number
            for row in body_rows
            if row["section"] == "discussion"
            for group in row["groups"]
            for number in group
        }
    )
    conclusion_numbers = sorted(
        {
            number
            for row in body_rows
            if row["section"] == "conclusion"
            for group in row["groups"]
            for number in group
        }
    )
    all_groups = [
        group for row in body_rows for group in row["groups"]
    ]
    first_seen = []
    for row in body_rows:
        for group in row["groups"]:
            for number in group:
                if number not in first_seen:
                    first_seen.append(number)
    repeated_numbers = {
        number: count
        for number, count in sorted(occurrences.items())
        if count > 1
    }
    abstract_index = next(
        index
        for index, text in nonempty
        if text == "Abstract"
    )
    long_sentences = []
    for index, paragraph in enumerate(document.paragraphs):
        text = paragraph.text.strip()
        if not text:
            continue
        if index < abstract_index or index >= references_index:
            continue
        if paragraph.style and paragraph.style.name.lower().startswith(
            ("title", "heading")
        ):
            continue
        if all(
            run.bold
            for run in paragraph.runs
            if run.text.strip()
        ) and len(
            re.findall(r"[A-Za-z0-9][A-Za-z0-9'./-]*", text)
        ) < 20:
            continue
        for row in sentence_lengths(text):
            if row["words"] > 30:
                long_sentences.append(
                    {
                        "paragraph_index": index,
                        **row,
                    }
                )
    report = {
        "source": str(SOURCE),
        "paragraphs": len(paragraphs),
        "tables": len(document.tables),
        "references": len(references),
        "reference_numbers": [row["number"] for row in references],
        "pre_discussion_unique": pre_numbers,
        "discussion_unique": discussion_numbers,
        "conclusion_unique": conclusion_numbers,
        "pre_discussion_count": len(pre_numbers),
        "discussion_count": len(discussion_numbers),
        "body_citation_occurrences": dict(sorted(occurrences.items())),
        "first_citation_order": first_seen,
        "first_citation_order_1_to_42": first_seen == list(range(1, 43)),
        "repeated_reference_numbers": repeated_numbers,
        "max_references_in_group": max(
            (len(group) for group in all_groups),
            default=0,
        ),
        "groups_with_more_than_two": [
            group for group in all_groups if len(group) > 2
        ],
        "groups_not_sorted_unique": [
            group
            for group in all_groups
            if group != sorted(set(group))
        ],
        "bookmarks": len(bookmarks),
        "hyperlink_targets": len(set(hyperlink_numbers)),
        "hyperlink_occurrences": len(hyperlink_numbers),
        "citation_occurrences": sum(occurrences.values()),
        "hyperlink_groups": hyperlink_groups,
        "hyperlink_groups_with_more_than_two": [
            group for group in hyperlink_groups if len(group) > 2
        ],
        "hyperlink_groups_not_sorted_unique": [
            group
            for group in hyperlink_groups
            if group != sorted(set(group))
        ],
        "all_citations_hyperlinked": set(hyperlink_numbers)
        == set(occurrences),
        "all_citation_occurrences_hyperlinked": Counter(hyperlink_numbers)
        == occurrences,
        "long_sentences_over_30": long_sentences,
        "paragraph_citations": paragraph_occurrences,
        "body_rows": body_rows,
        "reference_rows": references,
    }
    OUTPUT_JSON.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                key: report[key]
                for key in [
                    "references",
                    "pre_discussion_count",
                    "discussion_count",
                    "max_references_in_group",
                    "groups_with_more_than_two",
                    "groups_not_sorted_unique",
                    "bookmarks",
                    "hyperlink_targets",
                    "all_citations_hyperlinked",
                    "all_citation_occurrences_hyperlinked",
                    "first_citation_order_1_to_42",
                    "repeated_reference_numbers",
                ]
            },
            indent=2,
        )
    )
    print(
        "Long sentences over 30 words:",
        len(long_sentences),
    )
    for row in long_sentences[:30]:
        print(row["words"], row["sentence"][:220])


if __name__ == "__main__":
    main()
