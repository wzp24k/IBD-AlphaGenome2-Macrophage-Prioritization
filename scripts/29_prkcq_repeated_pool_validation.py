from __future__ import annotations

import importlib
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy, load_api_key


REPLICATES = 6
SEED = 20261112
PRKCQ_VARIANT = "chr10:6607973:C>T"


def main() -> None:
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers

    configure_proxy()
    conditioning = importlib.import_module("16_build_active_ibd_conditioning")
    model = dna_client.create(load_api_key(), timeout=30)
    variant = genome.Variant.from_str(PRKCQ_VARIANT)
    interval = variant.reference_interval.resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[
        "SC_GENE_EXPRESSION_DENOISED"
    ]
    source_cache = {}
    jobs = []
    for spec in conditioning.pool_specs():
        path = spec["path"]
        if path not in source_cache:
            source_cache[path] = conditioning.ad.read_h5ad(path, backed="r")
        backed = source_cache[path]
        eligible = spec["selector"](backed.obs.copy())
        if "n_counts" in eligible:
            eligible = eligible[
                pd.to_numeric(eligible["n_counts"], errors="coerce").fillna(0)
                >= 300
            ]
        for replicate in range(REPLICATES):
            sampled = conditioning.balanced_sample(
                eligible,
                n_cells=256,
                donor_column=spec["donor_column"],
                seed=SEED + 1000 * replicate,
            )
            cells, count_mode = conditioning.extract_counts(backed, sampled.index)
            jobs.append(
                {
                    "pool_name": spec["pool_name"],
                    "cohort": spec["cohort"],
                    "disease_state": spec["disease_state"],
                    "replicate": replicate,
                    "donors": sampled[spec["donor_column"]].nunique(),
                    "median_counts": pd.to_numeric(
                        sampled.get("n_counts", pd.Series(dtype=float)),
                        errors="coerce",
                    ).median(),
                    "count_mode": count_mode,
                    "cells": cells,
                }
            )

    def score_job(job):
        nonlocal_model = model
        scores = None
        last_error = None
        for attempt in range(3):
            try:
                scores = nonlocal_model.score_variant(
                    interval=interval,
                    variant=variant,
                    cell_conditioning=job["cells"],
                    cell_condition_name=job["pool_name"],
                    variant_scorers=[scorer],
                )
                break
            except Exception as exc:
                last_error = exc
                time.sleep(10 * (attempt + 1))
                nonlocal_model = dna_client.create(load_api_key(), timeout=30)
        if scores is None:
            raise RuntimeError(f"PRKCQ scoring failed: {last_error}")
        tidy = variant_scorers.tidy_scores(
            [scores], filter_protein_coding=True
        )
        target = tidy[tidy["gene_name"] == "PRKCQ"]
        if target.empty:
            raise ValueError(f"PRKCQ missing for {job['pool_name']}")
        return {
            "pool_name": job["pool_name"],
            "cohort": job["cohort"],
            "disease_state": job["disease_state"],
            "replicate": job["replicate"],
            "donors": job["donors"],
            "median_counts": job["median_counts"],
            "count_mode": job["count_mode"],
            "prkcq_effect": float(target["raw_score"].mean()),
        }

    cache_dir = ROOT / "work" / "prkcq_repeat_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for completed, job in enumerate(jobs, start=1):
        cache_path = cache_dir / (
            f"{job['pool_name']}_replicate_{job['replicate']}.csv"
        )
        if cache_path.exists():
            record = pd.read_csv(cache_path).iloc[0].to_dict()
        else:
            record = score_job(job)
            pd.DataFrame([record]).to_csv(cache_path, index=False)
        records.append(record)
        print(f"PRKCQ repeated pools: {completed}/{len(jobs)}")
    result = pd.DataFrame(records).sort_values(
        ["cohort", "disease_state", "replicate"]
    )
    table_dir = ROOT / "results" / "tables"
    result.to_csv(
        table_dir / "prkcq_repeated_pool_scores.csv", index=False
    )
    summary = (
        result.groupby(["cohort", "disease_state"], as_index=False)
        .agg(
            replicates=("prkcq_effect", "size"),
            median_effect=("prkcq_effect", "median"),
            mean_effect=("prkcq_effect", "mean"),
            sd_effect=("prkcq_effect", "std"),
            q025=("prkcq_effect", lambda x: x.quantile(0.025)),
            q975=("prkcq_effect", lambda x: x.quantile(0.975)),
            sign_positive=("prkcq_effect", lambda x: int((x > 0).sum())),
            sign_negative=("prkcq_effect", lambda x: int((x < 0).sum())),
        )
    )
    summary.to_csv(
        table_dir / "prkcq_repeated_pool_summary.csv", index=False
    )
    active = summary[summary["disease_state"] == "active"]
    max_abs_median = float(active["median_effect"].abs().max())
    all_ci_cross_zero = bool(
        ((active["q025"] <= 0) & (active["q975"] >= 0)).all()
    )
    if max_abs_median < 0.005 and all_ci_cross_zero:
        classification = "indeterminate"
    elif max_abs_median < 0.005:
        classification = "attenuated"
    else:
        classification = "persistent"
    pd.DataFrame(
        [
            {
                "variant": PRKCQ_VARIANT,
                "target": "PRKCQ",
                "classification": classification,
                "max_active_abs_median": max_abs_median,
                "all_active_ci_cross_zero": all_ci_cross_zero,
            }
        ]
    ).to_csv(
        table_dir / "prkcq_repeated_pool_classification.csv", index=False
    )
    print(summary.to_string(index=False))
    print(f"PRKCQ classification: {classification}")


if __name__ == "__main__":
    main()
