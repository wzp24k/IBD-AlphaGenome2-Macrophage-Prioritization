from __future__ import annotations

import random
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, post_graphql


ENSEMBL_URL = "https://rest.ensembl.org/overlap/region/human"
OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
AF_QUERY_TEMPLATE = """
query {
  %s
}
"""


def af_bin(af: float) -> str:
    if af < 0.05:
        return "0.01-0.05"
    if af < 0.1:
        return "0.05-0.10"
    if af < 0.2:
        return "0.10-0.20"
    if af < 0.3:
        return "0.20-0.30"
    return "0.30-0.50"


def random_region_pool(
    session: requests.Session,
    forbidden_regions: dict[str, list[tuple[int, int]]],
    target_pool_size: int = 5000,
    seed: int = 20261112,
) -> pd.DataFrame:
    rng = random.Random(seed)
    chrom_lengths = {
        str(chrom): length
        for chrom, length in {
            1: 248_956_422,
            2: 242_193_529,
            3: 198_295_559,
            4: 190_214_555,
            5: 181_538_259,
            6: 170_805_979,
            7: 159_345_973,
            8: 145_138_636,
            9: 138_394_717,
            10: 133_797_422,
            11: 135_086_622,
            12: 133_275_309,
            13: 114_364_328,
            14: 107_043_718,
            15: 101_991_189,
            16: 90_338_345,
            17: 83_257_441,
            18: 80_373_285,
            19: 58_617_616,
            20: 64_444_167,
            21: 46_709_983,
            22: 50_818_468,
        }.items()
    }
    records: list[dict] = []
    for window in range(20):
        chrom = rng.choice(list(chrom_lengths))
        start = rng.randint(10_000_000, chrom_lengths[chrom] - 30_000)
        stop = start + 20_000
        url = (
            f"{ENSEMBL_URL}/{chrom}:{start}..{stop}"
            "?feature=variation;content-type=application/json"
        )
        response = session.get(url, timeout=120)
        response.raise_for_status()
        variants = response.json()
        window_records = []
        for row in variants:
            alleles = row.get("alleles") or []
            if (
                row.get("feature_type") != "variation"
                or len(row.get("alleles") or []) < 2
                or len(alleles[0]) != 1
            ):
                continue
            snv_alts = [
                allele
                for allele in alleles[1:]
                if len(allele) == 1 and allele != alleles[0]
            ]
            if not snv_alts:
                continue
            pos = int(row["start"])
            if int(row["end"]) != pos:
                continue
            if any(
                left <= pos <= right
                for left, right in forbidden_regions.get(chrom, [])
            ):
                continue
            alt = rng.choice(snv_alts)
            ref = alleles[0]
            window_records.append(
                {
                    "canonical_variant": f"chr{chrom}:{pos}:{ref}>{alt}",
                    "open_targets_variant_id": f"{chrom}_{pos}_{ref}_{alt}",
                    "ensembl_rsid": row.get("id"),
                    "most_severe_consequence": row.get("consequence_type"),
                    "sampling_window": window,
                }
            )
        records.extend(window_records)
        print(
            f"Ensembl background pool: window {window + 1}, "
            f"cumulative variants {len(records)}"
        )
        if len(records) >= target_pool_size:
            break
        time.sleep(1)
    return pd.DataFrame(records).drop_duplicates("canonical_variant")


def annotate_frequency(
    session: requests.Session,
    pool: pd.DataFrame,
    batch_size: int = 40,
) -> pd.DataFrame:
    output = []
    records = pool.to_dict("records")
    for offset in range(0, len(records), batch_size):
        batch = records[offset : offset + batch_size]
        fields = []
        alias_map = {}
        for index, row in enumerate(batch):
            alias = f"v{index}"
            variant_id = row["open_targets_variant_id"]
            alias_map[alias] = row
            fields.append(
                f'{alias}: variant(variantId: "{variant_id}") '
                "{ id alleleFrequencies { populationName alleleFrequency } }"
            )
        data = post_graphql(
            session,
            OPEN_TARGETS_URL,
            AF_QUERY_TEMPLATE % "\n  ".join(fields),
            attempts=4,
        )["data"]
        for alias, row in alias_map.items():
            payload = data.get(alias) or {}
            frequencies = {
                item["populationName"]: item["alleleFrequency"]
                for item in payload.get("alleleFrequencies") or []
                if item.get("alleleFrequency") is not None
            }
            af_max = max(frequencies.values()) if frequencies else None
            enriched = dict(row)
            enriched.update(
                {
                    "af_max": af_max,
                    "af_matched": min(af_max, 1 - af_max) if af_max is not None else None,
                    "af_common": af_max is not None
                    and 0.01 <= min(af_max, 1 - af_max) <= 0.5,
                    "af_nfe": frequencies.get("nfe_adj"),
                    "af_eas": frequencies.get("eas_adj"),
                    "af_afr": frequencies.get("afr_adj"),
                }
            )
            output.append(enriched)
        print(f"Open Targets AF annotation: {min(offset + len(batch), len(records))}/{len(records)}")
        time.sleep(0.5)
    return pd.DataFrame(output)


def match_by_af(
    candidates: pd.DataFrame,
    pool: pd.DataFrame,
    sample_size: int,
    seed: int = 20261112,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    candidate_common = candidates[candidates["af_common"].fillna(False)].copy()
    pool_common = pool[pool["af_common"].fillna(False)].copy()
    candidate_common["af_bin"] = candidate_common["af_matched"].map(af_bin)
    pool_common["af_bin"] = pool_common["af_matched"].map(af_bin)
    proportions = candidate_common["af_bin"].value_counts(normalize=True)
    target = (proportions * sample_size).round().astype(int)
    selected = []
    for label, count in target.items():
        available = pool_common[pool_common["af_bin"] == label]
        if available.empty:
            continue
        take = min(int(count), len(available))
        indices = rng.choice(available.index.to_numpy(), size=take, replace=False)
        selected.append(available.loc[indices])
    if not selected:
        return pd.DataFrame()
    output = pd.concat(selected, ignore_index=True)
    if len(output) < sample_size:
        remaining = pool_common[~pool_common.index.isin(output.index)]
        fill = remaining.sample(
            n=min(sample_size - len(output), len(remaining)),
            random_state=seed,
        )
        output = pd.concat([output, fill], ignore_index=True)
    return output.sample(frac=1, random_state=seed).reset_index(drop=True)


def main() -> None:
    configure_proxy()
    candidates = pd.read_csv(
        ROOT / "results" / "tables" / "primary_variants_frequency.csv"
    )
    primary = pd.read_csv(
        ROOT / "data" / "processed" / "ibd_noncoding_credible_variants.csv"
    )
    primary = primary[primary["study_role"] == "primary_ibd"].drop_duplicates(
        "canonical_variant"
    )
    forbidden: dict[str, list[tuple[int, int]]] = {}
    for row in primary.itertuples():
        forbidden.setdefault(str(row.chromosome), []).append(
            (int(row.position) - 1_000_000, int(row.position) + 1_000_000)
        )

    session = requests.Session()
    pool = random_region_pool(session, forbidden_regions=forbidden)
    annotated_pool = annotate_frequency(session, pool)
    background = match_by_af(candidates, annotated_pool, sample_size=200)
    output_dir = ROOT / "data" / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)
    annotated_pool.to_csv(output_dir / "ensembl_background_pool_annotated.csv", index=False)
    background.to_csv(output_dir / "frequency_matched_background.csv", index=False)
    print(f"Background pool: {len(annotated_pool)}")
    print(f"Common background pool: {int(annotated_pool['af_common'].fillna(False).sum())}")
    print(f"Frequency-matched background: {len(background)}")
    print(background["af_bin"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
