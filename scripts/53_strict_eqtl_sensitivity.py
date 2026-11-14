from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


benchmark_module = importlib.import_module(
    "47_benchmark_models_and_transportability"
)


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    frame = pd.read_csv(table_dir / "benchmark_pair_level_dataset.csv")
    frame["eqtl_strict_bonferroni8"] = (
        pd.to_numeric(frame["minimum_pvalue"], errors="coerce")
        < 0.05 / 8
    ).astype(int)
    scores = {
        "PIP-only": "pip",
        "L2G-only": "l2g_score",
        "Atlas-only": "atlas_global_max_abs_q",
        "Live-only": "live_specificity",
        "Combined": "framework_pair_score",
        "Original priority": "priority_score",
    }
    rows = []
    for model, column in scores.items():
        metrics = benchmark_module.bootstrap_metrics(
            frame,
            "eqtl_strict_bonferroni8",
            column,
        )
        rows.append({"model": model, **metrics})
    result = pd.DataFrame(rows)
    result.to_csv(
        table_dir / "strict_eqtl_sensitivity_metrics.csv",
        index=False,
    )

    top = pd.read_csv(
        table_dir / "external_eqtl_top10_triangulation.csv"
    )
    top["strict"] = (
        pd.to_numeric(top["minimum_pvalue"], errors="coerce") < 0.05 / 8
    )
    matched = int(
        top.loc[
            top["query_group"].eq("top10_matched"),
            "strict",
        ].sum()
    )
    control = int(
        top.loc[
            top["query_group"].eq("mismatched_target_control"),
            "strict",
        ].sum()
    )
    odds_ratio, p_value = fisher_exact(
        [[matched, 10 - matched], [control, 10 - control]]
    )
    summary = {
        "strict_threshold": 0.05 / 8,
        "strict_positive_pairs": int(
            frame["eqtl_strict_bonferroni8"].sum()
        ),
        "top10_matched_strict_hits": matched,
        "top10_mismatched_strict_hits": control,
        "fisher_odds_ratio": float(odds_ratio),
        "fisher_p": float(p_value),
        "note": (
            "A conservative Bonferroni threshold across eight myeloid datasets "
            "was applied to exact-variant eQTL labels. Formal colocalization "
            "was not redefined by this label."
        ),
    }
    (table_dir / "strict_eqtl_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(result.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
