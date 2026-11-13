from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


REFERENCE_HEADER = "## References"
REFERENCE_END = "## Data availability"
CITATION_PATTERN = re.compile(r"\[(\d+(?:\s*[-,]\s*\d+)*)\]")
REFERENCE_PATTERN = re.compile(r"^(\d+)\.\s+(.+)$", re.MULTILINE)
SECTION_PATTERN = re.compile(r"^## (.+)$", re.MULTILINE)
BODY_ORDER = ["Introduction", "Methods", "Results", "Discussion"]
RETENTION_PRIORITY = {
    "Introduction": 4,
    "Results": 3,
    "Methods": 2,
    "Discussion": 1,
}


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
            start, end = (int(match.group(1)), int(match.group(2)))
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


def citation_occurrences(
    blocks: dict[str, str],
    source_order: list[str],
) -> list[dict]:
    occurrences = []
    for section_index, section_name in enumerate(source_order):
        block = blocks[section_name]
        for group_index, match in enumerate(CITATION_PATTERN.finditer(block)):
            numbers = sorted(set(parse_citation_group(match.group(1))))
            occurrences.append(
                {
                    "id": len(occurrences),
                    "section": section_name,
                    "section_index": section_index,
                    "group_index": group_index,
                    "start": match.start(),
                    "end": match.end(),
                    "numbers": numbers,
                    "priority": RETENTION_PRIORITY[section_name],
                }
            )
    return occurrences


def filter_selected_citations(
    block: str,
    occurrences: list[dict],
    selected_ids: dict[int, int],
) -> str:
    replacements = []
    for occurrence in occurrences:
        kept_numbers = [
            number
            for number in occurrence["numbers"]
            if selected_ids[number] == occurrence["id"]
        ]
        start = occurrence["start"]
        end = occurrence["end"]
        if kept_numbers:
            replacement = "[" + ",".join(map(str, kept_numbers)) + "]"
        else:
            replacement = ""
            if start > 0 and block[start - 1].isspace():
                start -= 1
        replacements.append((start, end, replacement))
    for start, end, replacement in sorted(replacements, reverse=True):
        block = block[:start] + replacement + block[end:]
    return block


def remap_citations(block: str, mapping: dict[int, int]) -> str:
    def replace(match: re.Match[str]) -> str:
        numbers = sorted(
            mapping[int(number)] for number in match.group(1).split(",")
        )
        return "[" + ",".join(map(str, numbers)) + "]"

    return CITATION_PATTERN.sub(replace, block)


def main() -> None:
    manuscript_path = ROOT / "manuscript" / "manuscript.md"
    text = manuscript_path.read_text(encoding="utf-8")
    preamble, section_items = split_sections(text)
    sections = dict(section_items)
    original_titles = [title for title, _ in section_items]
    if not all(name in sections for name in BODY_ORDER):
        raise SystemExit("One or more required body sections are missing.")
    body_start = min(original_titles.index(name) for name in BODY_ORDER)
    prefix_titles = original_titles[:body_start]

    reference_block = sections["References"].split(
        REFERENCE_HEADER, 1
    )[1].split(REFERENCE_END, 1)[0]
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

    source_order = [name for name in original_titles if name in BODY_ORDER]
    occurrences = citation_occurrences(sections, source_order)
    occurrences_by_section = defaultdict(list)
    occurrences_per_reference = defaultdict(list)
    for occurrence in occurrences:
        occurrences_by_section[occurrence["section"]].append(occurrence)
        for number in occurrence["numbers"]:
            occurrences_per_reference[number].append(occurrence)

    duplicate_references_before = sorted(
        number
        for number, items in occurrences_per_reference.items()
        if len(items) > 1
    )
    selected_ids = {}
    for number, items in occurrences_per_reference.items():
        selected_ids[number] = max(
            items,
            key=lambda item: (
                item["priority"],
                -item["section_index"],
                -item["group_index"],
            ),
        )["id"]

    filtered_blocks = {
        section_name: filter_selected_citations(
            sections[section_name],
            occurrences_by_section[section_name],
            selected_ids,
        )
        for section_name in BODY_ORDER
    }

    mapping = {}
    for section_name in BODY_ORDER:
        for match in CITATION_PATTERN.finditer(filtered_blocks[section_name]):
            for number in parse_citation_group(match.group(1)):
                mapping.setdefault(number, len(mapping) + 1)

    cited_numbers = set(mapping)
    uncited_references = sorted(set(reference_numbers) - cited_numbers)
    undefined_citations = sorted(cited_numbers - set(reference_numbers))
    first_citation_sequence_valid = (
        not uncited_references
        and not undefined_citations
        and list(mapping.values()) == expected_references
    )

    remapped_blocks = {
        section_name: remap_citations(block, mapping)
        for section_name, block in filtered_blocks.items()
    }
    inverse_mapping = {new: old for old, new in mapping.items()}
    reference_lookup = dict(references)
    remapped_reference_block = "\n\n" + "\n\n".join(
        f"{new_number}. {reference_lookup[inverse_mapping[new_number]]}"
        for new_number in expected_references
    ) + "\n\n"
    remapped_references_section = (
        REFERENCE_HEADER + "\n" + remapped_reference_block
    )

    post_sections = []
    for title in original_titles:
        if title in BODY_ORDER or title in {"Article type", "Abstract"}:
            continue
        if title == "References":
            post_sections.append(remapped_references_section)
        else:
            post_sections.append(sections[title])
    normalized_text = (
        preamble
        + "".join(sections[title] for title in prefix_titles)
        + "".join(remapped_blocks[name] for name in BODY_ORDER)
        + "".join(post_sections)
    )
    changed = normalized_text != text
    if changed:
        manuscript_path.write_text(normalized_text, encoding="utf-8")

    _, normalized_section_items = split_sections(normalized_text)
    normalized_sections = dict(normalized_section_items)
    final_occurrences = citation_occurrences(normalized_sections, BODY_ORDER)
    final_counts = defaultdict(int)
    for occurrence in final_occurrences:
        for number in occurrence["numbers"]:
            final_counts[number] += 1
    duplicate_references_after = sorted(
        number for number, count in final_counts.items() if count > 1
    )

    first_seen = []
    audit_rows = []
    for group_number, occurrence in enumerate(final_occurrences, start=1):
        for number in occurrence["numbers"]:
            if number not in first_seen:
                first_seen.append(number)
        audit_rows.append(
            {
                "citation_group": group_number,
                "section": occurrence["section"],
                "numbers": ";".join(map(str, occurrence["numbers"])),
                "ascending": occurrence["numbers"]
                == sorted(occurrence["numbers"]),
                "duplicate_within_group": len(occurrence["numbers"])
                != len(set(occurrence["numbers"])),
            }
        )
    first_citation_sequence_valid = first_seen == expected_references
    all_checks_pass = all(
        [
            reference_numbers_valid,
            not duplicate_references,
            not uncited_references,
            not undefined_citations,
            not duplicate_references_after,
            all(count == 1 for count in final_counts.values()),
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
        "citation_groups": len(final_occurrences),
        "references_cited_more_than_once_before": duplicate_references_before,
        "duplicate_reference_count_before": len(duplicate_references_before),
        "duplicate_references_after": duplicate_references_after,
        "every_reference_cited_once": all(
            count == 1 for count in final_counts.values()
        ),
        "unique_cited_references": len(cited_numbers),
        "uncited_references": uncited_references,
        "undefined_citations": undefined_citations,
        "first_citation_order": first_seen,
        "first_citation_order_1_to_n": first_citation_sequence_valid,
        "old_to_new_reference_map": {
            str(old): new for old, new in sorted(mapping.items())
        },
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
        f"- References cited more than once before: {len(duplicate_references_before)}",
        f"- References cited more than once after: {len(duplicate_references_after)}",
        f"- Every reference cited exactly once: {summary['every_reference_cited_once']}",
        f"- Uncited references: {uncited_references}",
        f"- Undefined citations: {undefined_citations}",
        f"- First citation order 1 to {len(references)}: {first_citation_sequence_valid}",
        f"- Markdown citation and reference numbering updated: {changed}",
        "",
        "Each reference is retained at the most directly supportive occurrence,",
        "renumbered in Word section order, and cited once. The Word manuscript adds",
        "a bookmark to each reference entry and an internal hyperlink to every",
        "in-text citation number. DOI and bibliographic text are unchanged.",
    ]
    (ROOT / "results" / "CITATION_CROSSREF_AUDIT.md").write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    if not all_checks_pass:
        raise SystemExit("Citation normalization or cross-reference audit failed.")


if __name__ == "__main__":
    main()
