from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


INPUT = ROOT / "8_Healthy_Pan-GI_atlas_Myeloid_20241119.h5ad"
COLORECTAL_ORGANS = {
    "colon",
    "transverse colon",
    "sigmoid colon",
    "ascending colon",
    "descending colon",
    "caecum",
    "rectum",
    "appendix",
}
ILEAL_ORGANS = {"ileum", "small intestine", "jejunum", "duodenum"}


def stratified_sample(
    obs: pd.DataFrame,
    organs: set[str],
    n_cells: int = 256,
    per_donor_cap: int = 8,
    seed: int = 20261112,
) -> pd.DataFrame:
    eligible = obs[
        (obs["level_2_annot"].astype(str) == "Macrophage")
        & (obs["control_vs_disease"].astype(str) == "control")
        & (obs["organ_unified"].astype(str).isin(organs))
        & (pd.to_numeric(obs["n_counts"], errors="coerce") >= 1_000)
    ].copy()
    rng = np.random.default_rng(seed)
    eligible["random_order"] = rng.random(len(eligible))
    eligible = eligible.sort_values(["random_order"])
    selected = (
        eligible.groupby("donorID_unified", sort=False)
        .head(per_donor_cap)
        .copy()
    )
    if len(selected) < n_cells:
        missing = eligible[~eligible.index.isin(selected.index)]
        selected = pd.concat([selected, missing.head(n_cells - len(selected))])
    elif len(selected) > n_cells:
        donor_balance = selected.groupby("donorID_unified").size().sort_values()
        take_indices = []
        cursor = 0
        ordered = selected.sort_values("random_order")
        while len(take_indices) < n_cells:
            progressed = False
            for donor in donor_balance.index:
                rows = ordered[(ordered["donorID_unified"] == donor)]
                if cursor < len(rows) and len(take_indices) < n_cells:
                    take_indices.append(rows.index[cursor])
                    progressed = True
            cursor += 1
            if not progressed:
                break
        selected = selected.loc[take_indices]
    return selected.sample(frac=1, random_state=seed)


def build_pool(backed: ad.AnnData, selected: pd.DataFrame, path: Path) -> None:
    subset = backed[selected.index, :].to_memory()
    if subset.raw is not None:
        raw_obs = subset.obs.copy()
        subset = subset.raw.to_adata()
        subset.obs = raw_obs
    values = subset.X
    if not np.issubdtype(values.dtype, np.number):
        raise TypeError("Expected a numeric count matrix.")
    sample_values = values.data if hasattr(values, "data") else values.ravel()
    if sample_values.size and not np.allclose(
        sample_values[: min(sample_values.size, 100_000)],
        np.round(sample_values[: min(sample_values.size, 100_000)]),
    ):
        raise ValueError("Expression conditioning requires raw integer UMI counts.")
    subset.X = values.astype(np.float32)
    subset.raw = None
    subset.var_names = subset.var_names.astype(str)
    subset.obs_names_make_unique()
    subset.write_h5ad(path, compression="gzip")


def main() -> None:
    output_dir = ROOT / "data" / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)
    backed = ad.read_h5ad(INPUT, backed="r")
    obs = backed.obs.copy()
    summary = []
    for name, organs in (
        ("healthy_colorectal_macrophage", COLORECTAL_ORGANS),
        ("healthy_ileal_macrophage", ILEAL_ORGANS),
    ):
        selected = stratified_sample(obs, organs)
        path = output_dir / f"{name}_256cells.h5ad"
        build_pool(backed, selected, path)
        summary.append(
            {
                "pool": name,
                "cells": len(selected),
                "donors": selected["donorID_unified"].nunique(),
                "samples": selected["sampleID"].nunique(),
                "organs": ",".join(sorted(selected["organ_unified"].unique())),
                "median_counts": selected["n_counts"].median(),
                "min_counts": selected["n_counts"].min(),
            }
        )
        print(f"Wrote {name}: {len(selected)} cells")
    pd.DataFrame(summary).to_csv(
        output_dir / "cell_conditioning_summary.csv", index=False
    )


if __name__ == "__main__":
    main()
