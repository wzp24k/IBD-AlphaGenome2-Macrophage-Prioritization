from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


def main() -> None:
    configure_proxy()
    manuscript = ROOT / "manuscript" / "manuscript.md"
    text = manuscript.read_text(encoding="utf-8")
    reference_text = text.split("## References", 1)[1].split(
        "## Data availability", 1
    )[0]
    entries = []
    for line in reference_text.splitlines():
        match = re.match(r"^(\d+)\.\s+(.+)$", line.strip())
        if match:
            entries.append((int(match.group(1)), match.group(2)))
    session = requests.Session()
    rows = []
    for number, reference in entries:
        doi_match = re.search(r"doi:(10\.[^\s]+)", reference)
        doi = doi_match.group(1).rstrip(".") if doi_match else ""
        row = {
            "number": number,
            "doi": doi,
            "crossref_status": None,
            "crossref_type": None,
            "crossref_journal": None,
            "crossref_title": None,
            "verified_journal_article": False,
        }
        if doi:
            try:
                response = session.get(
                    f"https://api.crossref.org/works/{doi}",
                    timeout=30,
                )
                row["crossref_status"] = response.status_code
                if response.status_code == 200:
                    message = response.json()["message"]
                    row["crossref_type"] = message.get("type")
                    journals = message.get("container-title", [])
                    titles = message.get("title", [])
                    row["crossref_journal"] = journals[0] if journals else ""
                    row["crossref_title"] = titles[0] if titles else ""
                    row["verified_journal_article"] = (
                        message.get("type") == "journal-article"
                    )
            except (requests.RequestException, ValueError, KeyError) as exc:
                row["crossref_status"] = f"error:{type(exc).__name__}"
        rows.append(row)
        time.sleep(0.15)
    frame = pd.DataFrame(rows)
    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(table_dir / "reference_doi_verification.csv", index=False)
    summary = {
        "references": len(frame),
        "with_doi": int(frame["doi"].ne("").sum()),
        "verified_journal_articles": int(
            frame["verified_journal_article"].sum()
        ),
        "all_verified": bool(
            len(frame) == 42
            and frame["verified_journal_article"].all()
            and frame["doi"].is_unique
        ),
        "failures": frame.loc[
            ~frame["verified_journal_article"],
            ["number", "doi", "crossref_status", "crossref_type"],
        ].to_dict("records"),
    }
    (table_dir / "reference_doi_verification_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    if not summary["all_verified"]:
        raise SystemExit("Reference DOI verification failed.")


if __name__ == "__main__":
    main()
