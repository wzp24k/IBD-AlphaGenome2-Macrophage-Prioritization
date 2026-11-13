from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, load_config


EXCLUDED_CONSEQUENCE_PATTERNS = [
    r"^missense_variant$",
    r"^synonymous_variant$",
    r"^stop_gained$",
    r"^stop_lost$",
    r"^start_lost$",
    r"^frameshift_variant$",
    r"^inframe_(insertion|deletion)$",
    r"^protein_altering_variant$",
    r"^splice_(acceptor|donor)_variant$",
    r"^splice_region_variant$",
]


def is_noncoding(consequence: str | None) -> bool:
    if not isinstance(consequence, str) or not consequence:
        return False
    terms = {term.strip() for term in consequence.split(",")}
    return not any(
        re.match(pattern, term)
        for term in terms
        for pattern in EXCLUDED_CONSEQUENCE_PATTERNS
    )


def build_canonical_variant(row: pd.Series) -> str:
    chromosome = str(row["chromosome"])
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"
    return (
        f"{chromosome}:{int(row['position'])}:"
        f"{row['referenceAllele']}>{row['alternateAllele']}"
    )


def main() -> None:
    config = load_config()
    raw_path = ROOT / "data" / "raw" / "open_targets_credible_sets_v2.csv"
    if not raw_path.exists():
        raise FileNotFoundError("Run scripts/01_fetch_credible_sets.py first.")

    frame = pd.read_csv(raw_path)
    numeric_columns = ["pip", "study_beta", "study_standardError", "r2Overall"]
    for column in numeric_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["canonical_variant"] = frame.apply(build_canonical_variant, axis=1)
    frame["is_noncoding"] = frame["most_severe_consequence"].map(is_noncoding)

    threshold = float(config["thresholds"]["pip_min"])
    eligible = frame[
        frame["is_noncoding"]
        & (frame["pip"] >= threshold)
        & frame["is95CredibleSet"].fillna(False).astype(bool)
    ].copy()

    eligible = eligible.sort_values(
        ["studyId", "studyLocusId", "pip"],
        ascending=[True, True, False],
    )
    eligible["rank_in_locus"] = eligible.groupby(
        ["studyId", "studyLocusId"]
    ).cumcount() + 1
    lead = eligible[eligible["rank_in_locus"] == 1].copy()

    processed_dir = ROOT / "data" / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)
    eligible.to_csv(processed_dir / "ibd_noncoding_credible_variants.csv", index=False)
    lead.to_csv(processed_dir / "ibd_lead_noncoding_variants.csv", index=False)

    summary = (
        eligible.groupby(["studyId", "study_role"], as_index=False)
        .agg(
            noncoding_variant_rows=("canonical_variant", "size"),
            unique_noncoding_variants=("canonical_variant", "nunique"),
            loci=("studyLocusId", "nunique"),
        )
        .sort_values(["study_role", "studyId"])
    )
    summary.to_csv(processed_dir / "variant_set_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Eligible variant rows: {len(eligible)}")
    print(f"Lead variants: {len(lead)}")


if __name__ == "__main__":
    main()
