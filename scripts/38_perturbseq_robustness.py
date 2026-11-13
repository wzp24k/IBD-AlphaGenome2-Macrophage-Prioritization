from __future__ import annotations

import importlib
import json
import re
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.sparse import csc_matrix, csr_matrix
from scipy.stats import ttest_1samp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT


BASE = importlib.import_module("34_analyze_gse327644_perturbseq")
TARGET_GENES = BASE.TARGET_GENES
PATHWAYS = BASE.PATHWAYS
RESPONSE_GENES = sorted(
    set(TARGET_GENES).union(
        gene for genes in PATHWAYS.values() for gene in genes
    )
)
THRESHOLDS = [3, 10]


def decode(values) -> np.ndarray:
    return np.array([value.decode("utf-8") for value in values])


def sample_metadata(path: Path) -> dict:
    match = re.search(r"_(L\d+)_Mac_(Ctrl|LPS)_", path.name)
    if not match:
        raise ValueError(f"Cannot parse sample name: {path.name}")
    donor, condition = match.groups()
    return {
        "sample": f"{donor}_{condition}",
        "donor": donor,
        "condition": condition,
    }


def process_sample(
    path: Path,
    guide_map: dict[str, str],
    chunk_size: int,
) -> pd.DataFrame:
    metadata = sample_metadata(path)
    with h5py.File(path, "r") as handle:
        feature_type = decode(handle["matrix/features/feature_type"][:])
        feature_id = decode(handle["matrix/features/id"][:])
        feature_name = decode(handle["matrix/features/name"][:])
        gene_indices = np.where(feature_type == "Gene Expression")[0]
        guide_indices = np.where(feature_type == "CRISPR Guide Capture")[0]
        gene_names = feature_name[gene_indices]
        response_positions = [
            index
            for index, gene in enumerate(gene_names)
            if gene in RESPONSE_GENES
        ]
        response_gene_ids = feature_id[gene_indices[response_positions]]
        response_gene_names = gene_names[response_positions]
        guide_names = feature_name[guide_indices]
        group_names = []
        for guide in guide_names:
            target = guide_map.get(guide)
            if target == "Non-Targeting Control":
                group_names.append("NTC")
            elif target in TARGET_GENES:
                group_names.append(guide)
            else:
                group_names.append(None)
        group_keys = sorted({group for group in group_names if group is not None})
        group_index = {group: index for index, group in enumerate(group_keys)}
        guide_map_index = np.full(len(guide_names), -1, dtype=np.int16)
        for index, group in enumerate(group_names):
            if group is not None:
                guide_map_index[index] = group_index[group]
        n_groups = len(group_keys)
        counts = {
            threshold: np.zeros(
                (len(response_positions), n_groups), dtype=np.float64
            )
            for threshold in THRESHOLDS
        }
        assigned = {
            threshold: np.zeros(n_groups, dtype=np.int64)
            for threshold in THRESHOLDS
        }
        library_counts = {
            threshold: np.zeros(n_groups, dtype=np.float64)
            for threshold in THRESHOLDS
        }
        indptr = handle["matrix/indptr"][:]
        n_features, n_cells = handle["matrix/shape"][:]
        for start in range(0, n_cells, chunk_size):
            stop = min(start + chunk_size, n_cells)
            indices = handle["matrix/indices"][indptr[start] : indptr[stop]]
            values = handle["matrix/data"][indptr[start] : indptr[stop]]
            lengths = np.diff(indptr[start : stop + 1])
            columns = np.repeat(np.arange(stop - start), lengths)
            matrix = csc_matrix(
                (values, (indices, columns)),
                shape=(n_features, stop - start),
            )
            guides = matrix[guide_indices, :].toarray().astype(np.int16)
            top_two = np.partition(guides, -2, axis=0)[-2:, :]
            second, maximum = top_two[0], top_two[1]
            dominant = guides.argmax(axis=0)
            mapped = guide_map_index[dominant]
            selected_gene_matrix = matrix[
                gene_indices[response_positions], :
            ]
            all_gene_matrix = matrix[gene_indices, :]
            cell_library_sizes = np.asarray(
                all_gene_matrix.sum(axis=0)
            ).ravel()
            for threshold in THRESHOLDS:
                confident = (
                    (maximum >= threshold)
                    & (second <= 0.5 * np.maximum(maximum, 1))
                    & (mapped >= 0)
                )
                if not confident.any():
                    continue
                cells = np.where(confident)[0]
                groups = mapped[cells]
                assigned[threshold] += np.bincount(groups, minlength=n_groups)
                library_counts[threshold] += np.bincount(
                    groups,
                    weights=cell_library_sizes[cells],
                    minlength=n_groups,
                )
                indicator = csr_matrix(
                    (
                        np.ones(len(cells), dtype=np.float64),
                        (groups, cells),
                    ),
                    shape=(n_groups, stop - start),
                )
                counts[threshold] += (
                    selected_gene_matrix @ indicator.T
                ).toarray()
        records = []
        for threshold in THRESHOLDS:
            for group, group_id in group_index.items():
                for local_index, (gene_index, gene_id, gene_name) in enumerate(zip(
                    response_positions,
                    response_gene_ids,
                    response_gene_names,
                )):
                    records.append(
                        {
                            **metadata,
                            "umi_threshold": threshold,
                            "group": group,
                            "gene_id": gene_id,
                            "gene_name": gene_name,
                            "count": counts[threshold][local_index, group_id],
                            "cells": assigned[threshold][group_id],
                            "library_counts": library_counts[threshold][group_id],
                        }
                    )
    return pd.DataFrame(records)


def process_and_cache(
    path: Path,
    guide_map: dict[str, str],
    work_dir: Path,
) -> pd.DataFrame:
    metadata = sample_metadata(path)
    cache_path = work_dir / f"{metadata['sample']}_thresholds.csv.gz"
    if cache_path.exists():
        cached = pd.read_csv(cache_path)
        if "library_counts" in cached.columns:
            return cached
    frame = process_sample(path, guide_map, 8000)
    frame.to_csv(cache_path, index=False, compression="gzip")
    return frame


def threshold_effect(
    pseudobulk: pd.DataFrame,
    target: str,
    threshold: int,
) -> list[dict]:
    wide = pseudobulk.pivot_table(
        index=["sample", "donor", "condition", "group"],
        columns="gene_name",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    group_library = (
        pseudobulk.groupby(["sample", "donor", "condition", "group"])[
            "library_counts"
        ].max()
    )
    guides = [
        group
        for group in wide.index.get_level_values("group").unique()
        if group.startswith(target + "-")
    ]
    target_counts = (
        wide[wide.index.get_level_values("group").isin(guides)]
        .groupby(level=["sample", "donor", "condition"])
        .sum()
    )
    target_library = (
        group_library[group_library.index.get_level_values("group").isin(guides)]
        .groupby(level=["sample", "donor", "condition"])
        .sum()
    )
    target_log = np.log2(
        target_counts.div(target_library, axis=0) * 1_000_000 + 1
    )
    ntc_counts = wide[wide.index.get_level_values("group") == "NTC"].droplevel(
        "group"
    )
    ntc_library = group_library[
        group_library.index.get_level_values("group") == "NTC"
    ].droplevel("group")
    ntc_log = np.log2(
        ntc_counts.div(ntc_library, axis=0) * 1_000_000 + 1
    )
    difference = target_log.sub(ntc_log, axis=0)
    rows = []
    metrics = {"self": target}
    for pathway, genes in PATHWAYS.items():
        present = [gene for gene in genes if gene in difference.columns]
        if present:
            difference[pathway] = difference[present].mean(axis=1)
            metrics[pathway] = pathway
    for condition in ["Ctrl", "LPS", "Combined"]:
        state = (
            difference
            if condition == "Combined"
            else difference[
                difference.index.get_level_values("condition") == condition
            ]
        )
        for metric_name, column in metrics.items():
            values = state[column].to_numpy(dtype=float)
            rows.append(
                {
                    "umi_threshold": threshold,
                    "target": target,
                    "condition": condition,
                    "metric": metric_name,
                    "pairs": len(values),
                    "mean_log2fc": float(np.mean(values)),
                    "p_value": float(
                        ttest_1samp(values, 0, nan_policy="omit").pvalue
                    )
                    if len(values) > 1
                    else np.nan,
                }
            )
    return rows


def leave_out(
    pseudobulk: pd.DataFrame,
    target: str,
    mode: str,
) -> pd.DataFrame:
    wide = pseudobulk.pivot_table(
        index=["sample", "donor", "condition", "group"],
        columns="gene_name",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    if wide.shape[1] < 1000:
        raise ValueError(
            "Leave-one-out analysis requires full-gene pseudobulk, not a "
            "selected-gene panel."
        )
    log_cpm = np.log2(wide.div(wide.sum(axis=1), axis=0) * 1_000_000 + 1)
    ntc = log_cpm[log_cpm.index.get_level_values("group") == "NTC"]
    ntc = ntc.droplevel("group")
    guides = sorted(
        group
        for group in log_cpm.index.get_level_values("group").unique()
        if group.startswith(target + "-")
    )
    full_aggregated = log_cpm[
        log_cpm.index.get_level_values("group").isin(guides)
    ].groupby(level=["sample", "donor", "condition"]).apply(
        lambda frame: np.log2(
            frame.sum() / frame.sum().sum() * 1_000_000 + 1
        )
    )
    full_difference = full_aggregated.sub(ntc, axis=0)
    for pathway, genes in PATHWAYS.items():
        present = [gene for gene in genes if gene in full_difference.columns]
        if present:
            full_difference[pathway] = full_difference[present].mean(axis=1)
    rows = []
    if mode == "guide":
        subsets = [
            (f"without_{guide}", [item for item in guides if item != guide])
            for guide in guides
        ]
    else:
        donors = sorted(
            log_cpm.index.get_level_values("donor").unique()
        )
        subsets = [
            (f"without_{donor}", guides)
            for donor in donors
        ]
    for label, selected_guides in subsets:
        selected = log_cpm[
            log_cpm.index.get_level_values("group").isin(selected_guides)
        ]
        if mode == "donor":
            donor = label.replace("without_", "")
            selected = selected[
                selected.index.get_level_values("donor") != donor
            ]
            selected_ntc = ntc[ntc.index.get_level_values("donor") != donor]
        else:
            selected_ntc = ntc
        aggregated = selected.groupby(level=["sample", "donor", "condition"]).apply(
            lambda frame: np.log2(
                frame.sum() / frame.sum().sum() * 1_000_000 + 1
            )
        )
        difference = aggregated.sub(selected_ntc, axis=0)
        metrics = {"self": target}
        for pathway, genes in PATHWAYS.items():
            present = [gene for gene in genes if gene in difference.columns]
            if present:
                difference[pathway] = difference[present].mean(axis=1)
                metrics[pathway] = pathway
        for condition in ["Ctrl", "LPS", "Combined"]:
            state = (
                difference
                if condition == "Combined"
                else difference[
                    difference.index.get_level_values("condition") == condition
                ]
            )
            for metric, column in metrics.items():
                values = state[column].to_numpy(dtype=float)
                full_state = (
                    full_difference
                    if condition == "Combined"
                    else full_difference[
                        full_difference.index.get_level_values("condition")
                        == condition
                    ]
                )
                rows.append(
                    {
                        "target": target,
                        "omitted": label,
                        "condition": condition,
                        "metric": metric,
                        "pairs": len(values),
                        "mean_log2fc": float(np.mean(values)),
                        "same_direction_as_full": (
                            np.sign(np.mean(values))
                            == np.sign(full_state[column].mean())
                        ),
                    }
                )
    return pd.DataFrame(rows)


def mixed_models(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    wide = pseudobulk.pivot_table(
        index=["sample", "donor", "condition", "group"],
        columns="gene_name",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    if wide.shape[1] < 1000:
        raise ValueError(
            "Mixed-model analysis requires full-gene pseudobulk, not a "
            "selected-gene panel."
        )
    log_cpm = np.log2(wide.div(wide.sum(axis=1), axis=0) * 1_000_000 + 1)
    rows = []
    for target in TARGET_GENES:
        guides = [
            group
            for group in log_cpm.index.get_level_values("group").unique()
            if group.startswith(target + "-")
        ]
        for condition in ["Ctrl", "LPS"]:
            condition_frame = log_cpm[
                log_cpm.index.get_level_values("condition") == condition
            ]
            ko = condition_frame[
                condition_frame.index.get_level_values("group").isin(guides)
            ].copy()
            ko["is_ko"] = 1
            ntc = condition_frame[
                condition_frame.index.get_level_values("group") == "NTC"
            ].copy()
            ntc["is_ko"] = 0
            data = pd.concat([ko, ntc])
            for metric, genes in [
                ("self", [target]),
                *[
                    (
                        pathway,
                        [
                            gene
                            for gene in pathway_genes
                            if gene in data.columns
                        ],
                    )
                    for pathway, pathway_genes in PATHWAYS.items()
                ],
            ]:
                present = [gene for gene in genes if gene in data.columns]
                if not present:
                    continue
                data["metric_value"] = data[present].mean(axis=1)
                model_type = "mixedlm_donor_random_intercept"
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        model = smf.mixedlm(
                            "metric_value ~ is_ko",
                            data=data.reset_index(),
                            groups=data.reset_index()["donor"],
                        ).fit(method="lbfgs", maxiter=200)
                    coefficient = float(model.params["is_ko"])
                    p_value = float(model.pvalues["is_ko"])
                    random_variance = float(model.cov_re.iloc[0, 0])
                    if (
                        not np.isfinite(coefficient)
                        or not np.isfinite(p_value)
                        or not np.isfinite(random_variance)
                        or random_variance <= 1e-8
                        or not getattr(model, "converged", True)
                    ):
                        raise ValueError("Singular or non-finite mixed model")
                except Exception:
                    model_type = "donor_fixed_effect_fallback"
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        fallback = smf.ols(
                            "metric_value ~ is_ko + C(donor)",
                            data=data.reset_index(),
                        ).fit()
                    coefficient = float(fallback.params["is_ko"])
                    p_value = float(fallback.pvalues["is_ko"])
                    random_variance = 0.0
                rows.append(
                    {
                        "target": target,
                        "condition": condition,
                        "metric": metric,
                        "coefficient": coefficient,
                        "p_value": p_value,
                        "random_intercept_variance": random_variance,
                        "model_type": model_type,
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    raw_dir = (
        ROOT
        / "data"
        / "external"
        / "experimental"
        / "GSE327644"
        / "raw"
    )
    table_dir = ROOT / "results" / "tables"
    work_dir = ROOT / "work" / "gse327644_robustness"
    work_dir.mkdir(parents=True, exist_ok=True)
    guides = BASE.guide_information()
    guide_map = dict(
        zip(guides["guide_id"], guides["guide_target_gene_name"])
    )
    for guide_id in guides.loc[
        guides["guide_role"].eq("Control"), "guide_id"
    ]:
        guide_map[guide_id] = "Non-Targeting Control"
    paths = sorted(raw_dir.glob("*.h5"))
    parts = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(process_and_cache, path, guide_map, work_dir): path
            for path in paths
        }
        for future in as_completed(futures):
            sample = future.result()
            parts.append(sample)
            print(
                f"processed {sample['sample'].iloc[0]}",
                flush=True,
            )
    threshold_pseudobulk = pd.concat(parts, ignore_index=True)
    threshold_pseudobulk.to_csv(
        table_dir / "perturbseq_threshold_pseudobulk.csv.gz",
        index=False,
        compression="gzip",
    )
    threshold_rows = []
    for threshold, frame in threshold_pseudobulk.groupby("umi_threshold"):
        for target in TARGET_GENES:
            threshold_rows.extend(
                threshold_effect(frame, target, int(threshold))
            )
    threshold_effects = pd.DataFrame(threshold_rows)
    lead_targets = pd.read_csv(
        table_dir / "perturbseq_gse327644_target_knockdown.csv"
    )
    lead_pathways = pd.read_csv(
        table_dir / "perturbseq_gse327644_pathway_effects.csv"
    )
    threshold5_rows = []
    for row in lead_targets.itertuples():
        threshold5_rows.append(
            {
                "umi_threshold": 5,
                "target": row.target,
                "condition": row.condition,
                "metric": "self",
                "pairs": row.pairs,
                "mean_log2fc": row.mean_log2fc_vs_ntc,
                "p_value": row.p_value,
            }
        )
    for row in lead_pathways.itertuples():
        threshold5_rows.append(
            {
                "umi_threshold": 5,
                "target": row.target,
                "condition": row.condition,
                "metric": row.pathway,
                "pairs": 4 if row.condition != "Combined" else 8,
                "mean_log2fc": row.mean_log2fc,
                "p_value": row.pathway_p,
            }
        )
    threshold_effects = pd.concat(
        [threshold_effects, pd.DataFrame(threshold5_rows)],
        ignore_index=True,
    ).sort_values(["target", "condition", "metric", "umi_threshold"])
    threshold_effects.to_csv(
        table_dir / "perturbseq_umi_threshold_sensitivity.csv",
        index=False,
    )
    lead_pseudobulk = pd.read_csv(
        table_dir / "perturbseq_gse327644_pseudobulk.csv.gz"
    )
    guide_rows = []
    donor_rows = []
    for target in TARGET_GENES:
        guide_rows.append(leave_out(lead_pseudobulk, target, "guide"))
        donor_rows.append(leave_out(lead_pseudobulk, target, "donor"))
    guide_out = pd.concat(guide_rows, ignore_index=True)
    donor_out = pd.concat(donor_rows, ignore_index=True)
    guide_out.to_csv(
        table_dir / "perturbseq_leave_one_guide_out.csv",
        index=False,
    )
    donor_out.to_csv(
        table_dir / "perturbseq_leave_one_donor_out.csv",
        index=False,
    )
    mixed = mixed_models(lead_pseudobulk)
    mixed.to_csv(
        table_dir / "perturbseq_mixed_model_effects.csv",
        index=False,
    )
    primary_metrics = threshold_effects[
        threshold_effects["condition"].eq("Combined")
        & threshold_effects["metric"].isin(
            ["self", "IL-23", "IL-6/JAK/STAT3"]
        )
        & threshold_effects["target"].isin(
            ["PDGFB", "PRKCQ", "IL23R"]
        )
    ]
    summary = {
        "assignment_thresholds": sorted(
            threshold_effects["umi_threshold"].unique().tolist()
        ),
        "threshold_primary_metrics": primary_metrics.to_dict("records"),
        "leave_one_guide_same_direction_fraction": guide_out.groupby(
            ["target", "metric"], as_index=False
        )["same_direction_as_full"].mean().to_dict("records"),
        "leave_one_donor_same_direction_fraction": donor_out.groupby(
            ["target", "metric"], as_index=False
        )["same_direction_as_full"].mean().to_dict("records"),
        "mixed_model": mixed.to_dict("records"),
    }
    (table_dir / "perturbseq_robustness_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(
        primary_metrics.to_string(index=False)
    )
    print(mixed.to_string(index=False))


if __name__ == "__main__":
    main()
