from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


EXTERNAL = ROOT / "data" / "external" / "ibd_single_cell"
OLIVER = EXTERNAL / "48ecc234-7073-4de9-8be6-1c2f08d4c146 (1).h5ad"
CROHN = EXTERNAL / "1ddca82d-2948-4dbe-8bba-be9de7115bb7.h5ad"
GARRIDO = EXTERNAL / "5e75dd8f-5c5f-4056-ade2-117748429022.h5ad"
COLORECTAL = {
    "colon",
    "caecum",
    "sigmoid colon",
    "transverse colon",
    "ascending colon",
    "descending colon",
    "rectum",
    "appendix",
    "intestine",
}


def balanced_sample(
    obs: pd.DataFrame,
    n_cells: int = 256,
    donor_column: str = "donor_id",
    seed: int = 20261112,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    eligible = obs.copy()
    eligible["_random_order"] = rng.random(len(eligible))
    eligible = eligible.sort_values("_random_order")
    donor_rows = {
        donor: rows
        for donor, rows in eligible.groupby(donor_column, sort=False)
    }
    selected = []
    cursor = 0
    while len(selected) < n_cells and donor_rows:
        progressed = False
        for donor in list(donor_rows):
            rows = donor_rows[donor]
            if cursor < len(rows):
                selected.append(rows.index[cursor])
                progressed = True
                if len(selected) == n_cells:
                    break
            else:
                del donor_rows[donor]
        if not progressed:
            break
        cursor += 1
    return eligible.loc[selected].drop(columns=["_random_order"])


def extract_counts(backed: ad.AnnData, indices: pd.Index) -> tuple[ad.AnnData, str]:
    subset = backed[indices, :].to_memory()
    if subset.raw is not None:
        obs = subset.obs.copy()
        subset = subset.raw.to_adata()
        subset.obs = obs
        mode = "raw.X"
    else:
        subset = subset[
            :, ~subset.var["feature_is_filtered"].fillna(False).to_numpy()
        ].copy()
        mode = "X"
    values = subset.X
    data = values.data if hasattr(values, "data") else np.asarray(values).ravel()
    if data.size and not np.allclose(
        data[: min(data.size, 100_000)],
        np.round(data[: min(data.size, 100_000)]),
    ):
        raise ValueError("Selected matrix does not contain integer counts.")
    subset.X = values.astype(np.float32)
    subset.raw = None
    subset.var_names = subset.var_names.astype(str)
    subset.obs_names_make_unique()
    return subset, mode


def pool_specs():
    oliver_disease = {
        "pool_name": "active_colitis_macrophage_oliver",
        "source": "Oliver2024",
        "cohort": "Oliver",
        "disease_state": "active",
        "disease_labels": ["ulcerative colitis", "inflammatory bowel disease"],
        "path": OLIVER,
        "selector": lambda o: o[
            (o["level_2_annot"].astype(str) == "Macrophage")
            & (o["organ_unified"].astype(str).isin(COLORECTAL))
            & (o["disease"].astype(str).isin(["ulcerative colitis", "inflammatory bowel disease"]))
        ],
        "donor_column": "donorID_unified",
    }
    oliver_healthy = {
        "pool_name": "healthy_colorectal_macrophage_oliver",
        "source": "Oliver2024",
        "cohort": "Oliver",
        "disease_state": "healthy",
        "disease_labels": ["normal"],
        "path": OLIVER,
        "selector": lambda o: o[
            (o["level_2_annot"].astype(str) == "Macrophage")
            & (o["organ_unified"].astype(str).isin(COLORECTAL))
            & (o["disease"].astype(str) == "normal")
        ],
        "donor_column": "donorID_unified",
    }
    crohn_active = {
        "pool_name": "active_crohn_macrophage_independent",
        "source": "Crohn2023",
        "cohort": "Crohn2023",
        "disease_state": "active",
        "disease_labels": ["Crohn disease"],
        "path": CROHN,
        "selector": lambda o: o[
            (o["cell_type"].astype(str) == "macrophage")
            & (o["disease"].astype(str) == "Crohn disease")
            & (o["Type"].astype(str) == "Infl")
        ],
        "donor_column": "donor_id",
    }
    crohn_healthy = {
        "pool_name": "healthy_colon_macrophage_crohn2023",
        "source": "Crohn2023",
        "cohort": "Crohn2023",
        "disease_state": "healthy",
        "disease_labels": ["normal"],
        "path": CROHN,
        "selector": lambda o: o[
            (o["cell_type"].astype(str) == "macrophage")
            & (o["disease"].astype(str) == "normal")
            & (o["Type"].astype(str) == "Heal")
        ],
        "donor_column": "donor_id",
    }
    garrido_active = {
        "pool_name": "active_ibd_myeloid_garrido",
        "source": "Garrido2023",
        "cohort": "Garrido2023",
        "disease_state": "active",
        "disease_labels": ["Crohn disease", "ulcerative colitis"],
        "path": GARRIDO,
        "selector": lambda o: o[
            (o["cell_type"].astype(str) == "myeloid cell")
            & (o["disease"].astype(str).isin(["Crohn disease", "ulcerative colitis"]))
        ],
        "donor_column": "donor_id",
    }
    garrido_healthy = {
        "pool_name": "healthy_myeloid_garrido",
        "source": "Garrido2023",
        "cohort": "Garrido2023",
        "disease_state": "healthy",
        "disease_labels": ["normal"],
        "path": GARRIDO,
        "selector": lambda o: o[
            (o["cell_type"].astype(str) == "myeloid cell")
            & (o["disease"].astype(str) == "normal")
        ],
        "donor_column": "donor_id",
    }
    return [
        oliver_disease,
        oliver_healthy,
        crohn_active,
        crohn_healthy,
        garrido_active,
        garrido_healthy,
    ]


def main() -> None:
    output_dir = ROOT / "data" / "processed" / "ibd_conditioning"
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    cache: dict[Path, ad.AnnData] = {}
    for spec in pool_specs():
        path = spec["path"]
        if path not in cache:
            cache[path] = ad.read_h5ad(path, backed="r")
        backed = cache[path]
        selected = spec["selector"](backed.obs.copy())
        if "n_counts" in selected:
            selected = selected[
                pd.to_numeric(selected["n_counts"], errors="coerce").fillna(0)
                >= 300
            ]
        sampled = balanced_sample(
            selected,
            n_cells=256,
            donor_column=spec["donor_column"],
        )
        cells, count_mode = extract_counts(backed, sampled.index)
        output_path = output_dir / f"{spec['pool_name']}_256cells.h5ad"
        cells.write_h5ad(output_path, compression="gzip")
        summary.append(
            {
                "pool_name": spec["pool_name"],
                "source": spec["source"],
                "cohort": spec["cohort"],
                "disease_state": spec["disease_state"],
                "cells": len(sampled),
                "donors": sampled[spec["donor_column"]].nunique(),
                "diseases": ",".join(sorted(sampled["disease"].astype(str).unique())),
                "tissues": ",".join(sorted(sampled["tissue"].astype(str).unique())),
                "count_mode": count_mode,
                "genes": cells.shape[1],
                "median_counts": pd.to_numeric(
                    sampled.get("n_counts", pd.Series(dtype=float)),
                    errors="coerce",
                ).median(),
            }
        )
        print(f"Wrote {spec['pool_name']}: {len(sampled)} cells")
    pd.DataFrame(summary).to_csv(
        output_dir / "ibd_conditioning_pool_summary.csv", index=False
    )


if __name__ == "__main__":
    main()
