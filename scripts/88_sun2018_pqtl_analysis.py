from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


TABLE_DIR = ROOT / "results" / "tables"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PQTL = load_module("pqtl73", Path(__file__).parent / "73_pqtl_mr_coloc_phewas.py")


def independent_instruments(pqtl: pd.DataFrame) -> pd.DataFrame:
    ld = pd.read_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_ld_matrix.csv",
        index_col=0,
    )
    candidates = (
        pqtl[pqtl["p"].lt(5e-8)]
        .drop_duplicates("rsid")
        .sort_values("p")
        .reset_index(drop=True)
    )
    selected = []
    for _, row in candidates.iterrows():
        rsid = str(row["rsid"])
        if rsid not in ld.index:
            continue
        if any(
            abs(float(ld.loc[rsid, str(old["rsid"])])) > 0.01
            for old in selected
            if str(old["rsid"]) in ld.columns
        ):
            continue
        selected.append(row)
    return pd.DataFrame(selected)


def main() -> None:
    pqtl = pd.read_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_cis_grch38.csv"
    )
    pqtl["position"] = pqtl["position_grch38"]
    outcomes = PQTL.load_gwas_catalog_regions()
    instruments = independent_instruments(pqtl)
    instruments.to_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_instruments.csv",
        index=False,
    )

    mr_rows = []
    coloc_rows = []
    for outcome_name, outcome in outcomes.items():
        harmonized = PQTL.harmonize(
            instruments,
            outcome,
            merge_key="rsid",
        )
        mr = PQTL.mr_estimate(harmonized)
        mr.update(
            {
                "outcome": outcome_name,
                "pqtl_dataset": "Sun2018_QTD000584",
            }
        )
        mr_rows.append(mr)
        coloc = PQTL.run_coloc(
            pqtl,
            outcome,
            "Sun2018_QTD000584",
            outcome_name,
        )
        coloc_rows.append(coloc)

    mr_table = pd.DataFrame(mr_rows)
    coloc_table = pd.DataFrame(coloc_rows)
    mr_table.to_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_mr_ibd.csv",
        index=False,
    )
    coloc_table.to_csv(
        TABLE_DIR / "prkcb_sun2018_pqtl_coloc.csv",
        index=False,
    )
    summary = {
        "pqtl_dataset": "Sun2018_QTD000584",
        "assembly": "GRCh38",
        "region": "chr16:23,300,000-23,950,000",
        "samples": 3301,
        "pqtl_variants": len(pqtl),
        "independent_instruments": len(instruments),
        "instrument_rsids": instruments["rsid"].tolist(),
        "mr_tests": len(mr_table),
        "positive_direction": int((mr_table["beta"] > 0).sum()),
        "nominal_p_lt_0_05": int((mr_table["p"] < 0.05).sum()),
        "coloc_tests": len(coloc_table),
        "coloc_h4_ge_0_8": int((coloc_table["H4"] >= 0.8).sum()),
        "max_h4": float(coloc_table["H4"].max()),
    }
    (TABLE_DIR / "prkcb_sun2018_pqtl_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(mr_table.to_string(index=False))
    print(coloc_table.to_string(index=False))


if __name__ == "__main__":
    main()
