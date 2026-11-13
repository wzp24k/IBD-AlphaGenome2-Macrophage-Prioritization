from __future__ import annotations

import json
import os
import random
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


GENOME = "hg38"
WINDOW_LENGTH = 251
WINDOW_OFFSETS = [-75, 0, 75]
BARCODES_PER_CONSTRUCT = 8

VARIANTS = [
    {
        "panel": "candidate",
        "label": "PDGFB",
        "rsid": "rs2143178",
        "chromosome": "chr22",
        "position": 39264824,
        "reference": "T",
        "alternate": "C",
        "target_gene": "PDGFB",
        "atlas_target_quantile": 0.999853,
        "priority_rank": 3,
    },
    {
        "panel": "candidate",
        "label": "PRKCB",
        "rsid": "rs7404003",
        "chromosome": "chr16",
        "position": 23853326,
        "reference": "C",
        "alternate": "A",
        "target_gene": "PRKCB",
        "atlas_target_quantile": -0.992305,
        "priority_rank": 5,
    },
    {
        "panel": "candidate",
        "label": "CHSY3",
        "rsid": "rs971891",
        "chromosome": "chr5",
        "position": 130237973,
        "reference": "C",
        "alternate": "T",
        "target_gene": "CHSY3",
        "atlas_target_quantile": -0.999255,
        "priority_rank": 8,
    },
    {
        "panel": "positive_control",
        "label": "ETS2_POS",
        "rsid": "rs2836882",
        "chromosome": "chr21",
        "position": 39094644,
        "reference": "G",
        "alternate": "A",
        "target_gene": "ETS2",
        "atlas_target_quantile": None,
        "priority_rank": None,
    },
    {
        "panel": "negative_control",
        "label": "IL2RA_NEG",
        "rsid": "rs12722489",
        "chromosome": "chr10",
        "position": 6060049,
        "reference": "C",
        "alternate": "T",
        "target_gene": "IL2RA",
        "atlas_target_quantile": None,
        "priority_rank": None,
    },
]


def fetch_json(session: requests.Session, url: str) -> dict:
    last_error = None
    for attempt in range(3):
        try:
            response = session.get(url, timeout=60)
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(2**attempt)
    raise RuntimeError(f"UCSC sequence request failed: {last_error}")


def fetch_reference(session: requests.Session, variant: dict) -> tuple[str, int]:
    position = int(variant["position"])
    start = position - 1 - 200
    end = position + 200
    url = (
        "https://api.genome.ucsc.edu/getData/sequence"
        f"?genome={GENOME};chrom={variant['chromosome']};"
        f"start={start};end={end}"
    )
    data = fetch_json(session, url)
    sequence = data["dna"].upper()
    if len(sequence) != 401:
        raise ValueError(
            f"Unexpected sequence length for {variant['rsid']}: {len(sequence)}"
        )
    if sequence[200] != variant["reference"]:
        raise ValueError(
            f"Reference mismatch for {variant['rsid']}: "
            f"expected {variant['reference']}, observed {sequence[200]}"
        )
    return sequence, start


def window_records(
    variant: dict,
    reference_sequence: str,
    reference_start: int,
) -> list[dict]:
    rows = []
    for offset in WINDOW_OFFSETS:
        sequence_start = reference_start + 75 + offset
        sequence = reference_sequence[
            75 + offset : 75 + offset + WINDOW_LENGTH
        ]
        variant_index = int(variant["position"]) - 1 - sequence_start
        if sequence[variant_index] != variant["reference"]:
            raise ValueError(f"Window reference mismatch for {variant['rsid']}")
        for allele, allele_role in [
            (variant["reference"], "reference"),
            (variant["alternate"], "alternate"),
        ]:
            allele_sequence = list(sequence)
            allele_sequence[variant_index] = allele
            context_code = f"{offset:+04d}".replace("+", "p").replace("-", "m")
            construct_id = (
                f"{variant['label']}_{variant['rsid']}_{context_code}_{allele_role}"
            )
            rows.append(
                {
                    "construct_id": construct_id,
                    "panel": variant["panel"],
                    "label": variant["label"],
                    "rsid": variant["rsid"],
                    "target_gene": variant["target_gene"],
                    "genome": GENOME,
                    "chromosome": variant["chromosome"],
                    "position_1based": int(variant["position"]),
                    "reference_allele": variant["reference"],
                    "alternate_allele": variant["alternate"],
                    "allele": allele,
                    "allele_role": allele_role,
                    "window_offset_bp": offset,
                    "window_start_0based": sequence_start,
                    "window_end_0based_exclusive": sequence_start
                    + WINDOW_LENGTH,
                    "variant_index_in_window_0based": variant_index,
                    "sequence_length": WINDOW_LENGTH,
                    "sequence": "".join(allele_sequence),
                    "alpha_genome_atlas_target_quantile": variant[
                        "atlas_target_quantile"
                    ],
                    "priority_rank": variant["priority_rank"],
                }
            )
    return rows


def make_scrambles(design: pd.DataFrame) -> list[dict]:
    source = design[
        design["panel"].eq("candidate")
        & design["allele_role"].eq("reference")
        & design["window_offset_bp"].eq(0)
    ].copy()
    rows = []
    for seed, source_row in zip([2143178, 971891], source.itertuples()):
        bases = list(source_row.sequence)
        rng = random.Random(seed)
        for _ in range(1000):
            rng.shuffle(bases)
            shuffled = "".join(bases)
            if (
                shuffled != source_row.sequence
                and not any(
                    motif in shuffled
                    for motif in ["AAAAAA", "TTTTTT", "GGGGGG", "CCCCCC"]
                )
            ):
                break
        rows.append(
            {
                "construct_id": f"SCRAMBLE_{source_row.label}",
                "panel": "scramble",
                "label": f"SCRAMBLE_{source_row.label}",
                "rsid": "none",
                "target_gene": "none",
                "genome": GENOME,
                "chromosome": "none",
                "position_1based": pd.NA,
                "reference_allele": "none",
                "alternate_allele": "none",
                "allele": "none",
                "allele_role": "scramble",
                "window_offset_bp": 0,
                "window_start_0based": pd.NA,
                "window_end_0based_exclusive": pd.NA,
                "variant_index_in_window_0based": pd.NA,
                "sequence_length": WINDOW_LENGTH,
                "sequence": shuffled,
                "alpha_genome_atlas_target_quantile": pd.NA,
                "priority_rank": pd.NA,
            }
        )
    return rows


def valid_barcode(sequence: str) -> bool:
    gc = sum(base in "GC" for base in sequence) / len(sequence)
    forbidden = ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG", "GAAGAC", "GTCTTC"]
    return (
        0.4 <= gc <= 0.6
        and len(set(sequence)) == 4
        and not any(base * 5 in sequence for base in "ACGT")
        and not any(motif in sequence for motif in forbidden)
    )


def make_barcodes(construct_ids: list[str]) -> pd.DataFrame:
    rng = random.Random(20261113)
    seen = set()
    rows = []
    for construct_id in construct_ids:
        count = 0
        while count < BARCODES_PER_CONSTRUCT:
            barcode = "".join(rng.choice("ACGT") for _ in range(20))
            if barcode in seen or not valid_barcode(barcode):
                continue
            seen.add(barcode)
            count += 1
            rows.append(
                {
                    "construct_id": construct_id,
                    "barcode_id": f"{construct_id}_BC{count:02d}",
                    "barcode": barcode,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    configure_proxy()
    output_dir = ROOT / "experimental_design" / "top3_mdmpra"
    output_dir.mkdir(parents=True, exist_ok=True)
    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    rows = []
    checks = []
    for variant in VARIANTS:
        reference_sequence, reference_start = fetch_reference(session, variant)
        checks.append(
            {
                "rsid": variant["rsid"],
                "chromosome": variant["chromosome"],
                "position_1based": variant["position"],
                "expected_reference": variant["reference"],
                "observed_reference": reference_sequence[200],
                "reference_match": reference_sequence[200]
                == variant["reference"],
            }
        )
        rows.extend(
            window_records(variant, reference_sequence, reference_start)
        )
    design = pd.DataFrame(rows)
    scrambles = pd.DataFrame(make_scrambles(design))
    design = pd.concat([design, scrambles], ignore_index=True)
    barcodes = make_barcodes(design["construct_id"].tolist())

    design.to_csv(output_dir / "construct_manifest.csv", index=False)
    barcodes.to_csv(output_dir / "barcode_manifest.csv", index=False)
    pd.DataFrame(checks).to_csv(
        table_dir / "top3_mdmpra_reference_checks.csv",
        index=False,
    )
    design[
        [
            "construct_id",
            "label",
            "rsid",
            "allele_role",
            "window_offset_bp",
            "sequence",
        ]
    ].to_csv(output_dir / "synthesis_order.csv", index=False)
    with (output_dir / "construct_sequences.fasta").open(
        "w", encoding="ascii"
    ) as handle:
        for row in design.itertuples():
            handle.write(f">{row.construct_id}\n{row.sequence}\n")

    summary = {
        "status": "designed_not_performed",
        "genome": GENOME,
        "window_length_bp": WINDOW_LENGTH,
        "window_offsets_bp": WINDOW_OFFSETS,
        "candidate_variants": [
            {
                key: variant[key]
                for key in [
                    "label",
                    "rsid",
                    "chromosome",
                    "position",
                    "reference",
                    "alternate",
                    "target_gene",
                    "priority_rank",
                ]
            }
            for variant in VARIANTS
            if variant["panel"] == "candidate"
        ],
        "constructs": int(len(design)),
        "candidate_constructs": int(
            design["panel"].eq("candidate").sum()
        ),
        "control_constructs": int(
            design["panel"].isin(
                ["positive_control", "negative_control"]
            ).sum()
        ),
        "scramble_constructs": int(design["panel"].eq("scramble").sum()),
        "proposed_barcodes": int(len(barcodes)),
        "barcodes_per_construct": BARCODES_PER_CONSTRUCT,
        "reference_checks": checks,
        "outputs": [
            "experimental_design/top3_mdmpra/construct_manifest.csv",
            "experimental_design/top3_mdmpra/synthesis_order.csv",
            "experimental_design/top3_mdmpra/construct_sequences.fasta",
            "experimental_design/top3_mdmpra/barcode_manifest.csv",
        ],
    }
    (table_dir / "top3_mdmpra_design_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
