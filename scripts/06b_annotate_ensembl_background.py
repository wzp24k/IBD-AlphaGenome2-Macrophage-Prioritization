from __future__ import annotations

import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


ENSEMBL_BATCH_URL = "https://rest.ensembl.org/variation/human?pops=1"
BATCH_SIZE = 100


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


def post_ids(ids: list[str]) -> dict:
    session = requests.Session()
    for attempt in range(5):
        response = session.post(
            ENSEMBL_BATCH_URL,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
            json={"ids": ids},
            timeout=240,
        )
        if response.status_code == 200:
            return response.json()
        if response.status_code == 429:
            time.sleep(5 * (attempt + 1))
            continue
        if response.status_code >= 500:
            time.sleep(5 * (attempt + 1))
            continue
        response.raise_for_status()
    raise RuntimeError("Ensembl batch variation query failed after retries.")


def alt_af_from_record(record: dict, alt: str) -> float | None:
    populations = record.get("populations") or []
    preferred = [
        "gnomADg:ALL",
        "gnomADe:ALL",
        "1000GENOMES:phase_3:ALL",
    ]
    for population in preferred:
        values = [
            row.get("frequency")
            for row in populations
            if row.get("population") == population
            and row.get("allele") == alt
            and row.get("frequency") is not None
        ]
        if values:
            return float(values[0])
    values = [
        float(row["frequency"])
        for row in populations
        if row.get("allele") == alt and row.get("frequency") is not None
    ]
    return max(values) if values else None


def annotate_batch(
    batch: list[dict],
    batch_index: int,
    cache_dir: Path,
) -> list[dict]:
    cache_path = cache_dir / f"batch_{batch_index:03d}.csv"
    if cache_path.exists():
        return pd.read_csv(cache_path).to_dict("records")
    ids = [row["ensembl_rsid"] for row in batch]
    payload = post_ids(ids)
    output = []
    for row in batch:
        record = payload.get(row["ensembl_rsid"]) or {}
        alt = row["canonical_variant"].split(">")[-1]
        alt_af = alt_af_from_record(record, alt)
        maf = min(alt_af, 1 - alt_af) if alt_af is not None else None
        enriched = dict(row)
        enriched.update(
            {
                "ensembl_alt_af": alt_af,
                "ensembl_maf": maf,
                "ensembl_common": maf is not None and 0.01 <= maf <= 0.5,
                "ensembl_var_class": record.get("var_class"),
                "ensembl_consequence": record.get("most_severe_consequence"),
            }
        )
        output.append(enriched)
    pd.DataFrame(output).to_csv(cache_path, index=False)
    return output


def match_by_af(
    candidates: pd.DataFrame,
    pool: pd.DataFrame,
    sample_size: int = 200,
    seed: int = 20261112,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    candidate_common = candidates[candidates["af_common"].fillna(False)].copy()
    pool_common = pool[pool["ensembl_common"].fillna(False)].copy()
    candidate_common["af_bin"] = candidate_common["af_matched"].map(af_bin)
    pool_common["af_bin"] = pool_common["ensembl_maf"].map(af_bin)
    target = (
        candidate_common["af_bin"].value_counts(normalize=True) * sample_size
    ).round().astype(int)
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
    return (
        output.drop_duplicates("canonical_variant")
        .sample(frac=1, random_state=seed)
        .reset_index(drop=True)
    )


def main() -> None:
    configure_proxy()
    input_path = (
        ROOT / "data" / "processed" / "ensembl_background_pool_annotated.csv"
    )
    pool = pd.read_csv(input_path)
    pool = pool[pool["ensembl_rsid"].notna()].drop_duplicates("ensembl_rsid")
    records = pool.to_dict("records")
    batches = [records[i : i + BATCH_SIZE] for i in range(0, len(records), BATCH_SIZE)]
    annotated: list[dict] = []
    cache_dir = ROOT / "work" / "ensembl_frequency_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            executor.submit(annotate_batch, batch, index, cache_dir): index
            for index, batch in enumerate(batches)
        }
        completed = 0
        for future in as_completed(futures):
            try:
                annotated.extend(future.result())
            except Exception as exc:
                print(f"Ensembl batch failed: {exc}")
            completed += 1
            print(
                f"Ensembl frequency annotation: "
                f"{completed}/{len(batches)} batches, {len(annotated)} variants"
            )
    annotated_pool = pd.DataFrame(annotated)
    valid_snv = annotated_pool["canonical_variant"].astype(str).str.fullmatch(
        r"chr[0-9]+:[0-9]+:[ACGT]>[ACGT]"
    )
    annotated_pool = annotated_pool[valid_snv].copy()
    candidates = pd.read_csv(
        ROOT / "results" / "tables" / "primary_variants_frequency.csv"
    )
    background = match_by_af(candidates, annotated_pool, sample_size=200)
    output_dir = ROOT / "data" / "processed"
    annotated_pool.to_csv(
        output_dir / "ensembl_background_frequency_annotated.csv", index=False
    )
    background.to_csv(output_dir / "frequency_matched_background.csv", index=False)
    print(f"Common background pool: {int(annotated_pool['ensembl_common'].fillna(False).sum())}")
    print(f"Frequency-matched background: {len(background)}")
    print(background["af_bin"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
