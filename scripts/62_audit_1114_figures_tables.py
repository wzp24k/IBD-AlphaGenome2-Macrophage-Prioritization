from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from docx import Document

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


SOURCE = ROOT / "1114manuscript.docx"
OUTPUT = ROOT / "results" / "1114_figure_table_audit.json"


def main() -> None:
    document = Document(SOURCE)
    paragraphs = [paragraph.text.strip() for paragraph in document.paragraphs]
    references_index = next(
        index
        for index, text in enumerate(paragraphs)
        if text == "7. References"
    )
    body = "\n".join(paragraphs[:references_index])
    table_mentions = [
        int(value) for value in re.findall(r"\bTable\s+(\d+)\b", body)
    ]
    figure_mentions = [
        int(value) for value in re.findall(r"\bFigure\s+(\d+)\b", body)
    ]
    table_captions = sorted(
        {
            int(match.group(1))
            for text in paragraphs
            if (
                match := re.match(
                    r"^Table\s+(\d+)\.\s+",
                    text,
                    flags=re.IGNORECASE,
                )
            )
        }
    )
    figure_captions = sorted(
        {
            int(match.group(1))
            for text in paragraphs
            if (
                match := re.match(
                    r"^Figure\s+(\d+)\.\s+",
                    text,
                    flags=re.IGNORECASE,
                )
            )
        }
    )
    report = {
        "source": str(SOURCE),
        "tables": len(document.tables),
        "figures": len(document.inline_shapes),
        "table_mentions": table_mentions,
        "table_mention_counts": {
            str(number): table_mentions.count(number)
            for number in sorted(set(table_mentions))
        },
        "table_captions": table_captions,
        "all_tables_mentioned": sorted(set(table_mentions))
        == list(range(1, 12)),
        "all_table_captions_present": table_captions == list(range(1, 12)),
        "figure_mentions": figure_mentions,
        "figure_mention_counts": {
            str(number): figure_mentions.count(number)
            for number in sorted(set(figure_mentions))
        },
        "figure_captions": figure_captions,
        "all_figures_mentioned": sorted(set(figure_mentions))
        == list(range(1, 13)),
        "all_figure_captions_present": figure_captions == list(range(1, 13)),
        "table_first_mention_order": list(dict.fromkeys(table_mentions)),
        "figure_first_mention_order": list(dict.fromkeys(figure_mentions)),
    }
    OUTPUT.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))
    if not all(
        [
            report["all_tables_mentioned"],
            report["all_table_captions_present"],
            report["all_figures_mentioned"],
            report["all_figure_captions_present"],
        ]
    ):
        raise SystemExit("Figure or table audit failed.")


if __name__ == "__main__":
    main()
