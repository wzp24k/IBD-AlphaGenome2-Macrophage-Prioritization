from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, spearmanr, ttest_ind

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


SEED = 20261114
POOL_DIR = ROOT / "data" / "processed" / "ibd_conditioning"
POOLS = {
    "Oliver_active": (
        POOL_DIR / "active_colitis_macrophage_oliver_256cells.h5ad",
        "Oliver2024",
        "active",
    ),
    "Oliver_healthy": (
        POOL_DIR / "healthy_colorectal_macrophage_oliver_256cells.h5ad",
        "Oliver2024",
        "healthy",
    ),
    "Crohn_active": (
        POOL_DIR / "active_crohn_macrophage_independent_256cells.h5ad",
        "Crohn2023",
        "active",
    ),
    "Crohn_healthy": (
        POOL_DIR / "healthy_colon_macrophage_crohn2023_256cells.h5ad",
        "Crohn2023",
        "healthy",
    ),
    "Garrido_active": (
        POOL_DIR / "active_ibd_myeloid_garrido_256cells.h5ad",
        "Garrido2023",
        "active",
    ),
    "Garrido_healthy": (
        POOL_DIR / "healthy_myeloid_garrido_256cells.h5ad",
        "Garrido2023",
        "healthy",
    ),
}
COHORTS = ["Oliver2024", "Crohn2023", "Garrido2023"]
PATHWAY_SETS = {
    "TNF/NF-kB": [
        "TNF",
        "IL1B",
        "IL6",
        "CXCL8",
        "CXCL10",
        "NFKB1",
        "NFKB2",
        "RELA",
        "RELB",
        "NFKBIA",
        "TNFAIP3",
        "ICAM1",
        "CCL2",
        "CCL3",
        "CCL4",
        "CXCL9",
        "CXCL11",
        "IRF1",
        "STAT1",
    ],
    "JAK/STAT": [
        "JAK1",
        "JAK2",
        "JAK3",
        "TYK2",
        "STAT1",
        "STAT2",
        "STAT3",
        "STAT4",
        "STAT5A",
        "STAT5B",
        "STAT6",
        "IFNGR1",
        "IFNGR2",
        "IL6R",
        "IL10RA",
        "IL10RB",
        "SOCS1",
        "SOCS3",
        "PTPN11",
        "IRF9",
    ],
    "IL-23": [
        "IL23A",
        "IL12B",
        "IL23R",
        "IL12RB1",
        "JAK2",
        "TYK2",
        "STAT3",
        "STAT4",
        "RORC",
        "IL17A",
        "IL17F",
        "IL22",
        "CCL20",
        "CXCL1",
        "CXCL2",
        "CXCL8",
        "SOCS3",
    ],
    "Chemokine": [
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "CCL2",
        "CCL3",
        "CCL4",
        "CCL5",
        "CCL20",
        "CXCL1",
        "CXCL2",
        "CXCL8",
        "CCR7",
    ],
}
FOCUS_GENES = sorted(
    set(
        [
            "PRKCB",
            "MAP3K8",
            *[gene for genes in PATHWAY_SETS.values() for gene in genes],
        ]
    )
)


def dense_counts(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray().astype(np.float32, copy=False)
    return np.asarray(matrix, dtype=np.float32)


def feature_symbols(adata) -> np.ndarray:
    for column in ["gene_symbols", "feature_name"]:
        if column in adata.var.columns:
            values = adata.var[column].astype(str).str.upper()
            if column == "feature_name" and values.str.startswith("ENSG").mean() > 0.8:
                continue
            return values.to_numpy()
    if "gene_symbols" in adata.raw.var.columns:
        return adata.raw.var["gene_symbols"].astype(str).str.upper().to_numpy()
    if "feature_name" in adata.raw.var.columns:
        return adata.raw.var["feature_name"].astype(str).str.upper().to_numpy()
    return pd.Index(adata.var_names.astype(str)).str.upper().to_numpy()


def load_pool(
    path: Path,
    selected_genes: set[str] | None = None,
) -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    adata = ad.read_h5ad(path)
    if adata.raw is not None:
        matrix = adata.raw.X
        var = adata.raw.var
        symbols = feature_symbols(
            type("VarView", (), {"var": var, "raw": None, "var_names": var.index})()
        )
    else:
        matrix = adata.X
        symbols = feature_symbols(adata)
    counts = dense_counts(matrix)
    obs = adata.obs.copy()
    donor_column = "donor_id"
    if donor_column not in obs:
        for candidate in ["donorID_unified", "sourceID", "biosample_id"]:
            if candidate in obs:
                donor_column = candidate
                break
    donors = obs[donor_column].astype(str).to_numpy()
    if selected_genes is not None:
        keep = np.array([symbol in selected_genes for symbol in symbols])
        counts = counts[:, keep]
        symbols = symbols[keep]
    cell_frame = pd.DataFrame(counts, columns=symbols)
    cell_frame["donor"] = donors
    counts_frame = (
        cell_frame.melt(
            id_vars="donor",
            var_name="gene",
            value_name="count",
        )
        .groupby(["donor", "gene"], as_index=False)["count"]
        .sum()
    )
    return counts_frame, donors, symbols.tolist(), np.array([])


def pseudobulk_log_cpm(counts_frame: pd.DataFrame) -> pd.DataFrame:
    matrix = counts_frame.pivot_table(
        index="donor",
        columns="gene",
        values="count",
        aggfunc="sum",
        fill_value=0,
    ).astype(float)
    library = matrix.sum(axis=1)
    normalized = np.log2(matrix.div(library, axis=0) * 1e6 + 1)
    return normalized


def load_donor_pseudobulk(path: Path) -> pd.DataFrame:
    adata = ad.read_h5ad(path)
    matrix = adata.raw.X if adata.raw is not None else adata.X
    counts = dense_counts(matrix)
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
    donors = adata.obs[donor_column].astype(str).to_numpy()
    rows = []
    index = []
    for donor in sorted(set(donors)):
        values = counts[donors == donor].sum(axis=0)
        rows.append(values)
        index.append(donor)
    frame = pd.DataFrame(np.vstack(rows), index=index, columns=symbols)
    frame = frame.T.groupby(level=0).sum().T
    library = frame.sum(axis=1)
    return np.log2(frame.div(library, axis=0) * 1e6 + 1)


def disease_signature(
    active: pd.DataFrame,
    healthy: pd.DataFrame,
    cohort: str,
) -> pd.DataFrame:
    shared = sorted(set(active.columns) & set(healthy.columns))
    active = active[shared]
    healthy = healthy[shared]
    mean_active = active.mean(axis=0)
    mean_healthy = healthy.mean(axis=0)
    statistic, p_value = ttest_ind(
        active.to_numpy(),
        healthy.to_numpy(),
        axis=0,
        equal_var=False,
        nan_policy="omit",
    )
    result = pd.DataFrame(
        {
            "cohort": cohort,
            "gene": shared,
            "disease_log2fc": (mean_active - mean_healthy).to_numpy(),
            "t_statistic": statistic,
            "p_value": p_value,
            "active_donors": active.shape[0],
            "healthy_donors": healthy.shape[0],
        }
    )
    result["gene"] = result["gene"].astype(str).str.upper()
    result = result.drop_duplicates("gene").set_index("gene")
    return result


def signature_projection(
    perturbation: pd.DataFrame,
    disease: pd.DataFrame,
    target: str,
    rng: np.random.Generator,
) -> dict:
    target_rows = perturbation[perturbation["overexpressed_gene"].eq(target)]
    target_delta = (
        target_rows.set_index("response_gene")["expression_delta"]
        .astype(float)
        .groupby(level=0)
        .mean()
    )
    target_delta.index = target_delta.index.astype(str).str.upper()
    shared = sorted(set(target_delta.index) & set(disease.index))
    x = target_delta.loc[shared].to_numpy()
    y = disease.loc[shared, "disease_log2fc"].to_numpy()
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if len(x) < 100:
        return {
            "cohort": disease["cohort"].iloc[0],
            "target": target,
            "genes": len(x),
            "correlation": np.nan,
            "p_permutation": np.nan,
            "direction": "insufficient_genes",
        }
    correlation, p_correlation = pearsonr(x, y)
    rank_correlation, _ = spearmanr(x, y)
    weights = np.sqrt(np.maximum(0, -np.log10(
        np.clip(disease.loc[shared, "p_value"].to_numpy()[valid], 1e-300, 1)
    )))
    weights = np.where(np.isfinite(weights), weights, 0)
    weights = weights / weights.sum() if weights.sum() else np.ones(len(x)) / len(x)
    weighted_score = float(np.sum(weights * np.sign(y) * x))
    null = np.empty(10000)
    n = len(x)
    for index in range(len(null)):
        sampled = rng.choice(len(target_delta), size=n, replace=False)
        null[index] = float(np.corrcoef(target_delta.iloc[sampled].to_numpy(), y)[0, 1])
    p_permutation = float((1 + np.sum(null >= correlation)) / (len(null) + 1))
    return {
        "cohort": disease["cohort"].iloc[0],
        "target": target,
        "genes": n,
        "correlation": float(correlation),
        "rank_correlation": float(rank_correlation),
        "correlation_p": float(p_correlation),
        "p_permutation": p_permutation,
        "weighted_knockdown_reversal": weighted_score,
        "direction": (
            "knockdown_reverses_disease"
            if correlation > 0
            else "knockdown_tracks_disease"
        ),
    }


def zscore_columns(matrix: np.ndarray) -> np.ndarray:
    mean = matrix.mean(axis=0, keepdims=True)
    standard = matrix.std(axis=0, ddof=0, keepdims=True)
    standard[standard < 1e-8] = 1
    return (matrix - mean) / standard


def residualize_by_donor(matrix: np.ndarray, donors: np.ndarray) -> np.ndarray:
    residual = matrix.copy()
    for donor in np.unique(donors):
        mask = donors == donor
        if mask.sum() <= 1:
            residual[mask] = 0
        else:
            residual[mask] -= residual[mask].mean(axis=0, keepdims=True)
    return zscore_columns(residual)


def shrunk_covariance(matrix: np.ndarray, shrinkage: float = 0.1) -> np.ndarray:
    covariance = np.cov(matrix, rowvar=False, ddof=1)
    diagonal = np.diag(np.diag(covariance))
    return (1 - shrinkage) * covariance + shrinkage * diagonal


def conditional_knockdown(
    precision: np.ndarray,
    genes: list[str],
    target: str,
    magnitude_sd: float = 1.0,
) -> pd.Series:
    target_index = genes.index(target)
    diagonal = np.diag(precision)[target_index]
    delta = -precision[:, target_index] / diagonal * (-magnitude_sd)
    result = pd.Series(delta, index=genes)
    result[target] = -magnitude_sd
    return result


def pathway_scores(
    delta: pd.Series,
    target: str,
) -> pd.DataFrame:
    rows = []
    for pathway, genes in PATHWAY_SETS.items():
        shared = [gene for gene in genes if gene in delta.index and gene != target]
        values = delta.loc[shared]
        rows.append(
            {
                "target": target,
                "pathway": pathway,
                "genes": len(shared),
                "mean_delta_z": float(values.mean()) if len(values) else np.nan,
                "down_fraction": float((values < 0).mean()) if len(values) else np.nan,
                "members": ";".join(shared),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_network(
    matrix: np.ndarray,
    genes: list[str],
    donors: np.ndarray,
    target: str,
    reps: int,
    rng: np.random.Generator,
    seed: int,
) -> pd.DataFrame:
    unique_donors = np.unique(donors)
    rows = []
    for reproduction in range(reps):
        local_rng = np.random.default_rng(seed + reproduction)
        sampled = local_rng.choice(unique_donors, size=len(unique_donors), replace=True)
        blocks = [matrix[donors == donor] for donor in sampled]
        boot = np.vstack(blocks)
        covariance = shrunk_covariance(boot)
        precision = np.linalg.pinv(covariance, hermitian=True)
        delta = conditional_knockdown(precision, genes, target)
        scores = pathway_scores(delta, target)
        scores["bootstrap"] = reproduction
        rows.append(scores)
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    rng = np.random.default_rng(SEED)
    table_dir = ROOT / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    selected = set(FOCUS_GENES)
    pool_data = {}
    for name, (path, cohort, state) in POOLS.items():
        counts, donors, symbols, _ = load_pool(path, selected)
        pool_data[name] = {
            "counts": counts,
            "donors": donors,
            "cohort": cohort,
            "state": state,
            "symbols": symbols,
        }
        print(f"Loaded {name}: {len(donors)} cells")

    signature_frames = []
    for cohort in COHORTS:
        short = cohort.split("2024")[0].split("2023")[0]
        active = load_donor_pseudobulk(POOLS[f"{short}_active"][0])
        healthy = load_donor_pseudobulk(POOLS[f"{short}_healthy"][0])
        signature = disease_signature(
            active,
            healthy,
            cohort,
        )
        signature_frames.append(signature)
    signatures = pd.concat(signature_frames)
    signatures.reset_index().to_csv(
        table_dir / "prkcb_public_disease_signatures.csv",
        index=False,
    )

    perturbation = pd.read_csv(
        table_dir / "experimental_thp1_overexpression_response.csv"
    )
    reversal_rows = []
    for cohort, disease in signatures.groupby("cohort"):
        for target in ["PRKCB", "CHSY3", "ETS2", "JAK2", "STAT3"]:
            reversal_rows.append(
                signature_projection(
                    perturbation,
                    disease,
                    target,
                    rng,
                )
            )
    reversal = pd.DataFrame(reversal_rows)
    reversal.to_csv(
        table_dir / "prkcb_signature_reversal.csv",
        index=False,
    )

    # Build state-specific cell matrices from targeted pathway genes.
    virtual_rows = []
    pathway_rows = []
    bootstrap_frames = []
    for state in ["healthy", "active"]:
        matrices = []
        lookups = []
        donors_state = []
        for cohort in COHORTS:
            short = cohort.split("2024")[0].split("2023")[0]
            name = f"{short}_{state}"
            # Cell-level matrices are reconstructed from the restricted h5ad.
            adata = ad.read_h5ad(POOLS[name][0])
            raw = adata.raw.X if adata.raw is not None else adata.X
            symbols = (
                adata.raw.var["gene_symbols"].astype(str).str.upper().to_numpy()
                if adata.raw is not None
                and "gene_symbols" in adata.raw.var.columns
                else (
                    adata.raw.var["feature_name"].astype(str).str.upper().to_numpy()
                    if adata.raw is not None
                    and "feature_name" in adata.raw.var.columns
                    else (
                        adata.var["gene_symbols"].astype(str).str.upper().to_numpy()
                        if "gene_symbols" in adata.var.columns
                        else adata.var["feature_name"].astype(str).str.upper().to_numpy()
                    )
                )
            )
            lookup = {
                gene: index
                for index, gene in enumerate(symbols)
                if gene in selected
            }
            raw_dense = dense_counts(raw)
            library = raw_dense.sum(axis=1, keepdims=True)
            matrix = np.log1p(
                raw_dense / np.maximum(library, 1) * 1e4
            )
            lookups.append(lookup)
            matrices.append(matrix)
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
            donors_state.extend(
                f"{cohort}:{donor}"
                for donor in adata.obs[donor_column].astype(str)
            )
        genes = sorted(
            set.intersection(*(set(lookup) for lookup in lookups))
        )
        if "PRKCB" not in genes:
            raise RuntimeError("PRKCB is absent from the virtual network.")
        aligned = [
            matrix[:, [lookup[gene] for gene in genes]]
            for matrix, lookup in zip(matrices, lookups)
        ]
        matrix = np.vstack(aligned)
        donors_array = np.asarray(donors_state)
        matrix = residualize_by_donor(matrix, donors_array)
        covariance = shrunk_covariance(matrix)
        precision = np.linalg.pinv(covariance, hermitian=True)
        for target in ["PRKCB", "MAP3K8"]:
            delta = conditional_knockdown(precision, genes, target)
            scores = pathway_scores(delta, target)
            scores["state"] = state
            scores["network_genes"] = len(genes)
            scores["method"] = "linear_gaussian_virtual_knockdown"
            pathway_rows.append(scores)
            top = delta.drop(labels=[target], errors="ignore").sort_values()
            for gene, value in top.head(25).items():
                virtual_rows.append(
                    {
                        "state": state,
                        "target": target,
                        "gene": gene,
                        "predicted_log_expression_delta": float(value),
                    }
                )
            # Random-gene null, matched by expression variance.
            candidates = [
                gene
                for gene in genes
                if gene not in PATHWAY_SETS.get("TNF/NF-kB", [])
            ]
            random_scores = []
            for pathway, pathway_genes in PATHWAY_SETS.items():
                indices = [
                    genes.index(gene)
                    for gene in pathway_genes
                    if gene in genes and gene != target
                ]
                random_indices = rng.choice(
                    len(genes),
                    size=min(2000, len(genes) * 20),
                    replace=True,
                )
                null_scores = []
                for random_index in random_indices:
                    if genes[random_index] == target:
                        continue
                    random_delta = (
                        precision[indices, random_index]
                        / precision[random_index, random_index]
                    )
                    null_scores.append(float(np.mean(random_delta)))
                observed = scores.loc[
                    scores["pathway"].eq(pathway),
                    "mean_delta_z",
                ].iloc[0]
                random_scores.append(
                    {
                        "state": state,
                        "target": target,
                        "pathway": pathway,
                        "observed_mean_delta_z": observed,
                        "random_null_mean": float(np.mean(null_scores)),
                        "random_null_sd": float(np.std(null_scores, ddof=1)),
                        "empirical_p_lower": float(
                            (1 + np.sum(np.asarray(null_scores) <= observed))
                            / (len(null_scores) + 1)
                        ),
                    }
                )
            scores = scores.merge(
                pd.DataFrame(random_scores)[
                    ["state", "target", "pathway", "empirical_p_lower"]
                ],
                on=["state", "target", "pathway"],
                how="left",
            )
            pathway_rows[-1] = scores
        bootstrap = bootstrap_network(
            matrix,
            genes,
            donors_array,
            "PRKCB",
            reps=200,
            rng=rng,
            seed=SEED,
        )
        bootstrap["state"] = state
        bootstrap_frames.append(bootstrap)
        print(f"Completed virtual network: {state}")

    virtual_genes = pd.DataFrame(virtual_rows)
    virtual_genes.to_csv(
        table_dir / "prkcb_virtual_knockdown_genes.csv",
        index=False,
    )
    virtual_pathways = pd.concat(pathway_rows, ignore_index=True)
    virtual_pathways.to_csv(
        table_dir / "prkcb_virtual_pathway_effects.csv",
        index=False,
    )
    bootstrap = pd.concat(bootstrap_frames, ignore_index=True)
    bootstrap.to_csv(
        table_dir / "prkcb_virtual_pathway_bootstrap.csv",
        index=False,
    )
    bootstrap_summary = (
        bootstrap.groupby(["state", "target", "pathway"], as_index=False)
        .agg(
            bootstrap_mean=("mean_delta_z", "mean"),
            bootstrap_sd=("mean_delta_z", "std"),
            ci_low=("mean_delta_z", lambda values: np.quantile(values, 0.025)),
            ci_high=("mean_delta_z", lambda values: np.quantile(values, 0.975)),
            positive_fraction=("mean_delta_z", lambda values: np.mean(values > 0)),
        )
    )
    bootstrap_summary.to_csv(
        table_dir / "prkcb_virtual_pathway_bootstrap_summary.csv",
        index=False,
    )

    primary = reversal[reversal["target"].eq("PRKCB")]
    cohort_consistency = float((primary["correlation"] > 0).mean())
    virtual_primary = virtual_pathways[virtual_pathways["target"].eq("PRKCB")]
    summary = {
        "method": "public perturbation signature projection and linear-Gaussian virtual knockdown",
        "virtual_perturbation_is_experimental": False,
        "cohorts": COHORTS,
        "prkcb_signature_reversal_cohorts_positive": int(
            (primary["correlation"] > 0).sum()
        ),
        "prkcb_signature_reversal_cohort_consistency": cohort_consistency,
        "median_prkcb_signature_correlation": float(
            primary["correlation"].median()
        ),
        "median_prkcb_permutation_p": float(
            primary["p_permutation"].median()
        ),
        "states_evaluated": ["healthy", "active"],
        "network_genes": int(
            virtual_primary["network_genes"].max()
        ),
        "pathways": {
            state: virtual_primary[virtual_primary["state"].eq(state)]
            .set_index("pathway")["mean_delta_z"]
            .to_dict()
            for state in ["healthy", "active"]
        },
        "boundary": (
            "The virtual knockdown is a conditional Gaussian network "
            "prediction. It is not CellOracle, Perturb-seq or primary "
            "macrophage CRISPR evidence. Public overexpression data are "
            "gene-level and not allele-specific."
        ),
    }
    (table_dir / "prkcb_virtual_network_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print(reversal.to_string(index=False))
    print(virtual_pathways.to_string(index=False))


if __name__ == "__main__":
    main()
