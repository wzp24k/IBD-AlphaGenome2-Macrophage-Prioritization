from __future__ import annotations

import importlib
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "32")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("MKL_NUM_THREADS", "32")

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


mechanism = importlib.import_module("55_prkcb_causal_mechanism")
virtual = importlib.import_module(
    "57_prkcb_perturbation_and_virtual_knockdown"
)
SEED = 20261114


def mr_leave_one_out() -> pd.DataFrame:
    table_dir = ROOT / "results" / "tables"
    instruments = pd.read_csv(table_dir / "prkcb_mr_instruments.csv")
    instruments = instruments[
        instruments["clump_r2_threshold"].eq(0.01)
    ].copy()
    rows = []
    rng = np.random.default_rng(SEED)
    for omitted in instruments["rsid"]:
        subset = instruments[instruments["rsid"].ne(omitted)].copy()
        bx = subset["eql_beta"].to_numpy(dtype=float)
        by = subset["gwas_beta"].to_numpy(dtype=float)
        sy = subset["gwas_se"].to_numpy(dtype=float)
        ivw = mechanism.mr_ivw(bx, by, sy)
        weighted = mechanism.mr_weighted_median(bx, by, sy, rng)
        rows.append(
            {
                "omitted_variant": omitted,
                "nsnp": len(subset),
                "ivw_beta": ivw["estimate"],
                "ivw_se": ivw["se"],
                "ivw_p": ivw["p"],
                "weighted_median_beta": weighted["estimate"],
                "weighted_median_p": weighted["p"],
                "mean_f_statistic": float(
                    subset["f_statistic"].mean()
                ),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(
        table_dir / "prkcb_mr_leave_one_out.csv",
        index=False,
    )
    return frame


def heidi_window_sensitivity() -> pd.DataFrame:
    table_dir = ROOT / "results" / "tables"
    variants = pd.read_csv(
        table_dir / "prkcb_ld_analysis_variants.csv"
    )
    ld = pd.read_csv(table_dir / "prkcb_ld_matrix.csv", index_col=0)
    lead_position = int(
        variants.loc[
            variants["rsid"].eq(mechanism.PRKCB_RSID),
            "position",
        ].iloc[0]
    )
    rows = []
    for half_width in [50_000, 100_000, 150_000, 250_000]:
        subset = variants[
            (variants["position"] - lead_position).abs() <= half_width
        ].copy()
        subset = subset.reset_index(drop=True)
        if len(subset) < 3:
            continue
        sub_ld = ld.loc[
            subset["rsid"].tolist(),
            subset["rsid"].tolist(),
        ].to_numpy(dtype=float)
        _, heidi = mechanism.smr_and_heidi(
            subset,
            sub_ld,
            mechanism.PRKCB_RSID,
        )
        primary = heidi[
            heidi["variant_set"].eq("all_ld_matched")
        ].iloc[0]
        rows.append(
            {
                "window_half_width_kb": half_width / 1_000,
                "variants": len(subset),
                "q": primary["q"],
                "effective_df": primary["effective_df"],
                "heterogeneity_p": primary["p"],
                "status": primary["heterogeneity_status"],
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(
        table_dir / "prkcb_heidi_window_sensitivity.csv",
        index=False,
    )
    return frame


def build_state_matrix(state: str) -> tuple[np.ndarray, list[str], np.ndarray]:
    matrices = []
    lookups = []
    donors_all = []
    for cohort in virtual.COHORTS:
        short = cohort.split("2024")[0].split("2023")[0]
        name = f"{short}_{state}"
        adata = ad.read_h5ad(virtual.POOLS[name][0])
        raw = adata.raw.X if adata.raw is not None else adata.X
        if adata.raw is not None:
            var = adata.raw.var
            if "gene_symbols" in var.columns:
                symbols = var["gene_symbols"].astype(str).str.upper().to_numpy()
            elif "feature_name" in var.columns:
                symbols = var["feature_name"].astype(str).str.upper().to_numpy()
            else:
                symbols = np.asarray([str(value).upper() for value in var.index])
        elif "gene_symbols" in adata.var.columns:
            symbols = adata.var["gene_symbols"].astype(str).str.upper().to_numpy()
        elif "feature_name" in adata.var.columns:
            symbols = adata.var["feature_name"].astype(str).str.upper().to_numpy()
        else:
            symbols = np.asarray([str(value).upper() for value in adata.var_names])
        lookup = {
            gene: index
            for index, gene in enumerate(symbols)
            if gene in set(virtual.FOCUS_GENES)
        }
        raw_dense = virtual.dense_counts(raw)
        library = raw_dense.sum(axis=1, keepdims=True)
        matrices.append(
            np.log1p(raw_dense / np.maximum(library, 1) * 1e4)
        )
        lookups.append(lookup)
        donor_column = "donor_id"
        if donor_column not in adata.obs:
            donor_column = next(
                candidate
                for candidate in [
                    "donorID_unified",
                    "sourceID",
                    "biosample_id",
                ]
                if candidate in adata.obs
            )
        donors_all.extend(
            f"{cohort}:{donor}"
            for donor in adata.obs[donor_column].astype(str)
        )
    genes = sorted(
        set.intersection(*(set(lookup) for lookup in lookups))
    )
    matrix = np.vstack(
        [
            values[:, [lookup[gene] for gene in genes]]
            for values, lookup in zip(matrices, lookups)
        ]
    )
    donors = np.asarray(donors_all)
    matrix = virtual.residualize_by_donor(matrix, donors)
    return matrix, genes, donors


def shrinkage_sensitivity(
    state_matrices: dict[str, tuple[np.ndarray, list[str], np.ndarray]],
) -> pd.DataFrame:
    rows = []
    for state, (matrix, genes, _) in state_matrices.items():
        for shrinkage in [0.02, 0.05, 0.10, 0.20, 0.30]:
            covariance = virtual.shrunk_covariance(
                matrix,
                shrinkage=shrinkage,
            )
            precision = np.linalg.pinv(covariance, hermitian=True)
            delta = virtual.conditional_knockdown(
                precision,
                genes,
                "PRKCB",
            )
            scores = virtual.pathway_scores(delta, "PRKCB")
            scores["state"] = state
            scores["shrinkage"] = shrinkage
            rows.append(scores)
    frame = pd.concat(rows, ignore_index=True)
    frame.to_csv(
        ROOT / "results" / "tables" / "prkcb_virtual_shrinkage_sensitivity.csv",
        index=False,
    )
    return frame


def bootstrap_one(
    matrix: np.ndarray,
    genes: list[str],
    donors: np.ndarray,
    reproduction: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + reproduction)
    unique_donors = np.unique(donors)
    sampled = rng.choice(
        unique_donors,
        size=len(unique_donors),
        replace=True,
    )
    boot = np.vstack([matrix[donors == donor] for donor in sampled])
    covariance = virtual.shrunk_covariance(boot, shrinkage=0.1)
    precision = np.linalg.pinv(covariance, hermitian=True)
    delta = virtual.conditional_knockdown(
        precision,
        genes,
        "PRKCB",
    )
    scores = virtual.pathway_scores(delta, "PRKCB")
    scores["bootstrap"] = reproduction
    return scores


def bootstrap_1000(
    state_matrices: dict[str, tuple[np.ndarray, list[str], np.ndarray]],
    replicates: int = 1_000,
) -> pd.DataFrame:
    jobs = []
    for state, (matrix, genes, donors) in state_matrices.items():
        for reproduction in range(replicates):
            jobs.append((state, matrix, genes, donors, reproduction))
    rows = []
    workers = min(32, max(2, os.cpu_count() or 8))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                bootstrap_one,
                matrix,
                genes,
                donors,
                reproduction,
            ): (state, reproduction)
            for state, matrix, genes, donors, reproduction in jobs
        }
        for completed, future in enumerate(
            as_completed(futures),
            start=1,
        ):
            state, reproduction = futures[future]
            frame = future.result()
            frame["state"] = state
            rows.append(frame)
            if completed % 250 == 0:
                print(
                    f"Completed {completed}/{len(futures)} bootstrap fits",
                    flush=True,
                )
    frame = pd.concat(rows, ignore_index=True)
    frame.to_csv(
        ROOT / "results" / "tables" / "prkcb_virtual_bootstrap_1000.csv",
        index=False,
    )
    return frame


def main() -> None:
    table_dir = ROOT / "results" / "tables"
    print("Running MR leave-one-out", flush=True)
    mr_loo = mr_leave_one_out()
    print("Running HEIDI window sensitivity", flush=True)
    heidi_windows = heidi_window_sensitivity()
    print("Building virtual-network matrices", flush=True)
    state_matrices = {
        state: build_state_matrix(state)
        for state in ["healthy", "active"]
    }
    print("Running shrinkage sensitivity", flush=True)
    shrinkage = shrinkage_sensitivity(state_matrices)
    print("Running 1000-replicate donor bootstrap", flush=True)
    bootstrap = bootstrap_1000(state_matrices, replicates=1_000)
    bootstrap_summary = (
        bootstrap.groupby(["state", "target", "pathway"], as_index=False)
        .agg(
            mean=("mean_delta_z", "mean"),
            sd=("mean_delta_z", "std"),
            ci_low=("mean_delta_z", lambda values: np.quantile(values, 0.025)),
            ci_high=("mean_delta_z", lambda values: np.quantile(values, 0.975)),
            positive_fraction=("mean_delta_z", lambda values: np.mean(values > 0)),
        )
    )
    bootstrap_summary.to_csv(
        table_dir / "prkcb_virtual_bootstrap_1000_summary.csv",
        index=False,
    )
    jak = bootstrap_summary[
        bootstrap_summary["state"].eq("active")
        & bootstrap_summary["pathway"].eq("JAK/STAT")
    ].iloc[0]
    shrinkage_jak = shrinkage[
        shrinkage["state"].eq("active")
        & shrinkage["pathway"].eq("JAK/STAT")
    ]
    summary = {
        "mr_leave_one_out_all_positive": bool(
            (mr_loo["ivw_beta"] > 0).all()
            and (mr_loo["weighted_median_beta"] > 0).all()
        ),
        "heidi_window_p": heidi_windows[
            ["window_half_width_kb", "heterogeneity_p", "status"]
        ].to_dict("records"),
        "virtual_shrinkage_jak_stat": shrinkage_jak[
            ["shrinkage", "mean_delta_z"]
        ].to_dict("records"),
        "virtual_bootstrap_1000_active_jak_stat": {
            "mean": float(jak["mean"]),
            "ci_low": float(jak["ci_low"]),
            "ci_high": float(jak["ci_high"]),
            "positive_fraction": float(jak["positive_fraction"]),
        },
        "bootstrap_replicates_per_state": 1_000,
        "workers": min(32, max(2, os.cpu_count() or 8)),
    }
    (table_dir / "prkcb_mechanism_robustness_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
