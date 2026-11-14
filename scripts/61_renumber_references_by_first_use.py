from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


MANUSCRIPT = ROOT / "manuscript" / "manuscript.md"
BACKUP = (
    ROOT
    / "manuscript"
    / "archive"
    / "manuscript_pre_1114_reference_renumber.md"
)
TABLE_DIR = ROOT / "results" / "tables"
CITATION_PATTERN = re.compile(r"\[(\d+(?:\s*,\s*\d+)*)\]")
REFERENCE_PATTERN = re.compile(r"^(\d+)\.\s+(.+)$", re.MULTILINE)
SECTION_PATTERN = re.compile(r"^## (.+)$", re.MULTILINE)
BODY_ORDER = ["Introduction", "Methods", "Results", "Discussion"]


def split_sections(text: str) -> tuple[str, list[tuple[str, str]]]:
    matches = list(SECTION_PATTERN.finditer(text))
    preamble = text[: matches[0].start()]
    sections = []
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        sections.append((match.group(1).strip(), text[match.start() : end]))
    return preamble, sections


def parse_group(value: str) -> list[int]:
    return [int(part.strip()) for part in value.split(",") if part.strip()]


def replace_citations(
    block: str,
    mapping: dict[int, int],
) -> str:
    def replace(match: re.Match[str]) -> str:
        old_numbers = parse_group(match.group(1))
        mapped = sorted({mapping[number] for number in old_numbers})
        return "[" + ",".join(str(number) for number in mapped) + "]"

    return CITATION_PATTERN.sub(replace, block)


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    _, section_items = split_sections(text)
    sections = dict(section_items)
    for required in [*BODY_ORDER, "References"]:
        if required not in sections:
            raise SystemExit(f"Missing required section: {required}")

    first_seen: list[int] = []
    occurrence_count: dict[int, int] = {}
    for section_name in BODY_ORDER:
        for match in CITATION_PATTERN.finditer(sections[section_name]):
            for number in parse_group(match.group(1)):
                occurrence_count[number] = occurrence_count.get(number, 0) + 1
                if number not in first_seen:
                    first_seen.append(number)
    if len(first_seen) != 42:
        raise SystemExit(
            f"Expected 42 distinct cited references, found {len(first_seen)}."
        )
    repeated = {
        number: count
        for number, count in occurrence_count.items()
        if count > 1
    }
    if repeated:
        raise SystemExit(f"References cited more than once: {repeated}")
    mapping = {
        old_number: new_number
        for new_number, old_number in enumerate(first_seen, start=1)
    }

    reference_block = sections["References"].split(
        "## References",
        1,
    )[1].split("## Data availability", 1)[0]
    references = [
        (int(number), title)
        for number, title in REFERENCE_PATTERN.findall(reference_block)
    ]
    if sorted(number for number, _ in references) != list(range(1, 43)):
        raise SystemExit("Reference list is not sequential before renumbering.")
    reference_lookup = {
        old_number: title
        for old_number, title in references
    }

    BACKUP.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MANUSCRIPT, BACKUP)
    updated_text = text
    for section_name, original_block in section_items:
        if section_name not in BODY_ORDER:
            continue
        updated_block = replace_citations(original_block, mapping)
        updated_text = updated_text.replace(
            original_block,
            updated_block,
            1,
        )

    _, updated_items = split_sections(updated_text)
    updated_sections = dict(updated_items)
    old_reference_block = updated_sections["References"].split(
        "## References",
        1,
    )[1].split("## Data availability", 1)[0]
    new_reference_lines = []
    for new_number in range(1, 43):
        old_number = next(
            old
            for old, mapped in mapping.items()
            if mapped == new_number
        )
        new_reference_lines.append(
            f"{new_number}. {reference_lookup[old_number]}"
        )
    new_reference_block = "\n\n" + "\n\n".join(new_reference_lines) + "\n\n"
    updated_text = updated_text.replace(
        old_reference_block,
        new_reference_block,
        1,
    )
    MANUSCRIPT.write_text(updated_text, encoding="utf-8")

    mapping_rows = []
    for old_number in range(1, 43):
        mapping_rows.append(
            {
                "old_reference": old_number,
                "new_reference": mapping[old_number],
                "title": reference_lookup[old_number],
            }
        )
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(mapping_rows).sort_values("new_reference").to_csv(
        TABLE_DIR / "reference_renumbering_1114.csv",
        index=False,
    )
    summary = {
        "references": len(references),
        "repeated_before": repeated,
        "first_use_old_order": first_seen,
        "mapping_old_to_new": {
            str(old): new for old, new in sorted(mapping.items())
        },
        "discussion_new_numbers": sorted(
            {
                number
                for match in CITATION_PATTERN.finditer(
                    updated_sections["Discussion"]
                )
                for number in parse_group(match.group(1))
            }
        ),
        "backup": str(BACKUP.relative_to(ROOT)),
    }
    (
        TABLE_DIR / "reference_renumbering_1114_summary.json"
    ).write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
