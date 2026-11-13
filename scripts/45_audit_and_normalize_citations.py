from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


REFERENCE_HEADER = "## References"
REFERENCE_END = "## Data availability"
CITATION_PATTERN = re.compile(r"\[(\d+(?:\s*[-,]\s*\d+)*)\]")
REFERENCE_PATTERN = re.compile(r"^(\d+)\.\s+(.+)$", re.MULTILINE)


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


def main() -> None:
    manuscript_path = ROOT / "manuscript" / "manuscript.md"
    text = manuscript_path.read_text(encoding="utf-8")
    main_text, remainder = text.split(REFERENCE_HEADER, 1)
    reference_block, tail = remainder.split(REFERENCE_END, 1)

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
    cited_numbers: set[int] = set()
    duplicate_groups = 0
    unsorted_groups = 0
    expanded_groups = 0
    malformed_groups = 0

    for match_index, match in enumerate(CITATION_PATTERN.finditer(main_text), start=1):
        original = match.group(0)
        try:
            numbers = parse_citation_group(match.group(1))
        except ValueError:
            malformed_groups += 1
            audit_rows.append(
                {
                    "citation_group": match_index,
                    "original": original,
                    "normalized": "",
                    "numbers": "",
                    "duplicate_within_group": False,
                    "ascending": False,
                    "expanded": False,
                }
            )
            continue

        unique_numbers = sorted(set(numbers))
        has_duplicates = len(unique_numbers) != len(numbers)
        ascending = numbers == unique_numbers
        expanded = "-" in match.group(1)
        duplicate_groups += int(has_duplicates)
        unsorted_groups += int(not ascending)
        expanded_groups += int(expanded)
        for number in unique_numbers:
            cited_numbers.add(number)
            if number not in first_seen:
                first_seen.append(number)
        audit_rows.append(
            {
                "citation_group": match_index,
                "original": original,
                "normalized": "[" + ",".join(map(str, unique_numbers)) + "]",
                "numbers": ";".join(map(str, unique_numbers)),
                "duplicate_within_group": has_duplicates,
                "ascending": ascending,
                "expanded": expanded,
            }
        )

    def normalize(match: re.Match[str]) -> str:
        numbers = sorted(set(parse_citation_group(match.group(1))))
        return "[" + ",".join(map(str, numbers)) + "]"

    normalized_main = CITATION_PATTERN.sub(normalize, main_text)
    normalized_text = normalized_main + REFERENCE_HEADER + remainder
    changed = normalized_text != text
    if changed:
        manuscript_path.write_text(normalized_text, encoding="utf-8")

    uncited_references = sorted(set(reference_numbers) - cited_numbers)
    undefined_citations = sorted(cited_numbers - set(reference_numbers))
    first_citation_sequence_valid = first_seen == expected_references
    all_checks_pass = all(
        [
            reference_numbers_valid,
            not duplicate_references,
            not malformed_groups,
            not duplicate_groups,
            not unsorted_groups,
            not uncited_references,
            not undefined_citations,
            first_citation_sequence_valid,
        ]
    )

    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    audit_frame = pd.DataFrame(audit_rows)
    audit_frame.to_csv(table_dir / "reference_citation_audit.csv", index=False)

    summary = {
        "all_checks_pass": all_checks_pass,
        "references": len(references),
        "reference_numbers_sequential": reference_numbers_valid,
        "duplicate_reference_numbers": duplicate_references,
        "citation_groups": len(audit_rows),
        "expanded_range_groups": expanded_groups,
        "duplicate_groups_before_normalization": duplicate_groups,
        "unsorted_groups_before_normalization": unsorted_groups,
        "malformed_groups": malformed_groups,
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
        f"- Range groups expanded to individual ascending numbers: {expanded_groups}",
        f"- Duplicate groups before normalization: {duplicate_groups}",
        f"- Unsorted groups before normalization: {unsorted_groups}",
        f"- Uncited references: {uncited_references}",
        f"- Undefined citations: {undefined_citations}",
        f"- First citation order 1 to {len(references)}: {first_citation_sequence_valid}",
        f"- Markdown citation formatting updated: {changed}",
        "",
        "The Word manuscript adds a bookmark to each reference entry and an internal",
        "hyperlink to every in-text citation number. DOI and bibliographic text are",
        "unchanged.",
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
