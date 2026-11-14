from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


coloc_module = importlib.import_module("48_top10_coloc_and_negative_controls")


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    variants = pd.read_csv(table_dir / "top10_eqtl_coloc_variants.csv")
    rows = []
    for target, frame in variants.groupby("target"):
        for sd_prior in [0.05, 0.10, 0.15, 0.20]:
            for p12 in [1e-5, 5e-5]:
                posterior, _ = coloc_module.coloc_abf(
                    frame,
                    p1=1e-4,
                    p2=1e-4,
                    p12=p12,
                    sd_prior=sd_prior,
                )
                rows.append(
                    {
                        "target": target,
                        "sd_prior": sd_prior,
                        "p12": p12,
                        "variants_harmonized": len(frame),
                        **posterior,
                    }
                )
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(
        table_dir / "top10_eqtl_coloc_sensitivity.csv",
        index=False,
    )
    summary_rows = []
    for target, frame in sensitivity.groupby("target"):
        summary_rows.append(
            {
                "target": target,
                "minimum_H4": float(frame["H4"].min()),
                "maximum_H4": float(frame["H4"].max()),
                "primary_H4": float(
                    frame.loc[
                        frame["sd_prior"].eq(0.15)
                        & frame["p12"].eq(1e-5),
                        "H4",
                    ].iloc[0]
                ),
                "sensitivity_supported": bool(
                    (
                        frame.loc[
                            frame["sd_prior"].eq(0.15)
                            & frame["p12"].eq(1e-5),
                            "H4",
                        ].iloc[0]
                        >= 0.80
                    )
                    and frame["H4"].max() >= 0.80
                ),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values("target")
    summary.to_csv(
        table_dir / "top10_eqtl_coloc_sensitivity_summary.csv",
        index=False,
    )
    output = {
        "targets_evaluated": len(summary),
        "targets_supported": int(summary["sensitivity_supported"].sum()),
        "targets": summary.to_dict("records"),
        "note": (
            "Sensitivity varied sd_prior from 0.05 to 0.20 and p12 from 1e-5 "
            "to 5e-5. Support required primary H4 >= 0.8 and at least one "
            "sensitivity setting with H4 >= 0.8."
        ),
    }
    (table_dir / "top10_eqtl_coloc_sensitivity_summary.json").write_text(
        json.dumps(output, indent=2),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
