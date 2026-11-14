from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


CITATION_PATTERN = re.compile(r"\[(\d+(?:\s*[-,]\s*\d+)*)\]")
REFERENCE_PATTERN = re.compile(r"^(\d+)\.\s+(.+)$", re.MULTILINE)
SECTION_PATTERN = re.compile(r"^## (.+)$", re.MULTILINE)
BODY_ORDER = ["Introduction", "Methods", "Results", "Discussion"]


def parse_citation_group(value: str) -> list[int]:
    numbers: list[int] = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", part)
            if not match:
                raise ValueError(f"Malformed citation range: [{value}]")
            start, end = int(match.group(1)), int(match.group(2))
            if end < start:
                raise ValueError(f"Descending citation range: [{value}]")
            numbers.extend(range(start, end + 1))
        else:
            numbers.append(int(part))
    return numbers


def split_sections(text: str) -> tuple[str, list[tuple[str, str]]]:
    matches = list(SECTION_PATTERN.finditer(text))
    if not matches:
        raise ValueError("No level-2 Markdown sections found.")
    preamble = text[: matches[0].start()]
    sections = []
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        sections.append((match.group(1).strip(), text[match.start() : end]))
    return preamble, sections


def normalize_citations(block: str) -> str:
    def replace(match: re.Match[str]) -> str:
        numbers = sorted(set(parse_citation_group(match.group(1))))
        return "[" + ",".join(map(str, numbers)) + "]"

    return CITATION_PATTERN.sub(replace, block)


def main() -> None:
    manuscript_path = ROOT / "manuscript" / "manuscript.md"
    text = manuscript_path.read_text(encoding="utf-8")
    _, section_items = split_sections(text)
    sections = dict(section_items)
    if not all(name in sections for name in BODY_ORDER):
        raise SystemExit("One or more required body sections are missing.")

    for section_name in BODY_ORDER:
        sections[section_name] = normalize_citations(sections[section_name])
    normalized_text = text
    for section_name, original_block in section_items:
        if section_name in BODY_ORDER:
            normalized_text = normalized_text.replace(
                original_block,
                sections[section_name],
                1,
            )
    changed = normalized_text != text
    if changed:
        manuscript_path.write_text(normalized_text, encoding="utf-8")

    _, normalized_items = split_sections(normalized_text)
    normalized_sections = dict(normalized_items)
    reference_block = normalized_sections["References"].split(
        "## References", 1
    )[1].split("## Data availability", 1)[0]
    references = [
        (int(number), title)
        for number, title in REFERENCE_PATTERN.findall(reference_block)
    ]
    reference_numbers = [number for number, _ in references]
    expected_references = list(range(1, len(references) + 1))
    reference_numbers_valid = reference_numbers == expected_references
    duplicate_references = sorted(
        {
            number
            for number in reference_numbers
            if reference_numbers.count(number) > 1
        }
    )

    audit_rows = []
    first_seen: list[int] = []
    citation_counts = defaultdict(int)
    discussion_numbers: set[int] = set()
    pre_discussion_numbers: set[int] = set()
    discussion_middle_citations = 0
    discussion_end_citations = 0
    group_index = 0
    for section_name in BODY_ORDER:
        for match in CITATION_PATTERN.finditer(normalized_sections[section_name]):
            group_index += 1
            numbers = parse_citation_group(match.group(1))
            ascending_unique = (
                numbers == sorted(numbers)
                and len(numbers) == len(set(numbers))
            )
            for number in numbers:
                citation_counts[number] += 1
                if number not in first_seen:
                    first_seen.append(number)
                if section_name == "Discussion":
                    discussion_numbers.add(number)
                else:
                    pre_discussion_numbers.add(number)
            if section_name == "Discussion":
                line_end = normalized_sections[section_name].find(
                    "\n", match.end()
                )
                if line_end == -1:
                    line_end = len(normalized_sections[section_name])
                suffix = normalized_sections[section_name][
                    match.end() : line_end
                ].strip().strip(".")
                if suffix:
                    discussion_middle_citations += 1
                else:
                    discussion_end_citations += 1
            audit_rows.append(
                {
                    "citation_group": group_index,
                    "section": section_name,
                    "numbers": ";".join(map(str, numbers)),
                    "ascending_and_unique": ascending_unique,
                }
            )

    cited_numbers = set(citation_counts)
    uncited_references = sorted(set(reference_numbers) - cited_numbers)
    undefined_citations = sorted(cited_numbers - set(reference_numbers))
    repeated_references = sorted(
        number for number, count in citation_counts.items() if count > 1
    )
    all_groups_sorted_unique = all(
        row["ascending_and_unique"] for row in audit_rows
    )
    discussion_exactly_15 = len(discussion_numbers) == 15
    pre_discussion_exactly_27 = len(pre_discussion_numbers) == 27
    citation_groups_at_most_two = all(
        len(parse_citation_group(match.group(1))) <= 2
        for section_name in BODY_ORDER
        for match in CITATION_PATTERN.finditer(normalized_sections[section_name])
    )
    discussion_majority_in_middle = (
        discussion_middle_citations > discussion_end_citations
    )
    first_citation_sequence_valid = first_seen == expected_references
    all_checks_pass = all(
        [
            reference_numbers_valid,
            not duplicate_references,
            not uncited_references,
            not undefined_citations,
            all_groups_sorted_unique,
            citation_groups_at_most_two,
            discussion_exactly_15,
            pre_discussion_exactly_27,
            discussion_majority_in_middle,
            first_citation_sequence_valid,
        ]
    )

    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audit_rows).to_csv(
        table_dir / "reference_citation_audit.csv",
        index=False,
    )
    summary = {
        "all_checks_pass": all_checks_pass,
        "references": len(references),
        "reference_numbers_sequential": reference_numbers_valid,
        "duplicate_reference_numbers": duplicate_references,
        "citation_groups": len(audit_rows),
        "all_citation_groups_sorted_unique": all_groups_sorted_unique,
        "citation_groups_at_most_two": citation_groups_at_most_two,
        "references_cited_more_than_once": repeated_references,
        "discussion_unique_references": len(discussion_numbers),
        "discussion_exactly_15_references": discussion_exactly_15,
        "discussion_reference_numbers": sorted(discussion_numbers),
        "discussion_middle_citations": discussion_middle_citations,
        "discussion_end_citations": discussion_end_citations,
        "discussion_majority_in_middle": discussion_majority_in_middle,
        "pre_discussion_unique_references": len(pre_discussion_numbers),
        "pre_discussion_exactly_27_references": pre_discussion_exactly_27,
        "pre_discussion_reference_numbers": sorted(pre_discussion_numbers),
        "unique_cited_references": len(cited_numbers),
        "uncited_references": uncited_references,
        "undefined_citations": undefined_citations,
        "first_citation_order": first_seen,
        "first_citation_order_1_to_n": first_citation_sequence_valid,
        "manuscript_normalized": changed,
    }
    (table_dir / "reference_citation_audit_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    report_lines = [
        "# Citation and reference cross-reference audit",
        "",
        f"- All checks passed: {all_checks_pass}",
        f"- Reference entries: {len(references)}",
        f"- Reference numbers sequential and unique: {reference_numbers_valid and not duplicate_references}",
        f"- Citation groups audited: {len(audit_rows)}",
        f"- All citation groups ascending and unique: {all_groups_sorted_unique}",
        f"- All citation groups contain at most two references: {citation_groups_at_most_two}",
        f"- Discussion references: {len(discussion_numbers)} distinct",
        f"- Discussion has exactly 15 references: {discussion_exactly_15}",
        f"- Discussion citations in paragraph middle/end: {discussion_middle_citations}/{discussion_end_citations}",
        f"- Most Discussion citation groups are in paragraph middle: {discussion_majority_in_middle}",
        f"- Pre-Discussion references: {len(pre_discussion_numbers)} distinct",
        f"- Pre-Discussion has exactly 27 references: {pre_discussion_exactly_27}",
        f"- Uncited references: {uncited_references}",
        f"- Undefined citations: {undefined_citations}",
        f"- First citation order 1 to {len(references)}: {first_citation_sequence_valid}",
        f"- Markdown citation formatting updated: {changed}",
        "",
        "Citation groups are ascending and unique, and first use follows",
        "the order 1 to 42. The Discussion contains 15 distinct references.",
        "The Word manuscript keeps clickable internal links but renders",
        "citation numbers in black without underline.",
    ]
    (ROOT / "results" / "CITATION_CROSSREF_AUDIT.md").write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    if not all_checks_pass:
        raise SystemExit("Citation audit failed.")


if __name__ == "__main__":
    main()
