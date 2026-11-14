from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


benchmark_module = importlib.import_module(
    "47_benchmark_models_and_transportability"
)


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    frame = pd.read_csv(table_dir / "benchmark_pair_level_dataset.csv")
    tasks = [
        ("External genetic support, both cohorts", "external_both", None),
        ("Replicated exact-variant eQTL", "eqtl_replicated", None),
        ("Macrophage TNF CRISPR", "crispr_tnf_supported", "crispr_evaluable"),
    ]
    rows = []
    for task, label, evaluable in tasks:
        subset = frame.copy()
        if evaluable:
            subset = subset[subset[evaluable].astype(bool)].copy()
        for score_name, score_column in [
            ("Original priority score", "priority_score"),
            ("Benchmark equal-weight score", "framework_pair_score"),
        ]:
            metrics = benchmark_module.bootstrap_metrics(
                subset,
                label,
                score_column,
            )
            rows.append(
                {
                    "task": task,
                    "score": score_name,
                    **metrics,
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(
        table_dir / "benchmark_priority_sensitivity.csv",
        index=False,
    )
    correlation = spearmanr(
        frame["priority_score"].fillna(0),
        frame["framework_pair_score"].fillna(0),
    ).statistic
    summary = {
        "spearman_priority_vs_benchmark_score": float(correlation),
        "note": (
            "The original variant priority score and the benchmark equal-weight "
            "score were compared on the same rows. The analysis checks that "
            "conclusions do not depend on post hoc score construction."
        ),
    }
    (table_dir / "benchmark_priority_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(result.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
