from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT, configure_proxy


TARGET_GENES = [
    "PDGFB",
    "PRKCQ",
    "IL23R",
    "TNFSF15",
    "MAP3K8",
    "TNF",
    "TNFAIP3",
]
PATHWAYS = {
    "IL-6/JAK/STAT3": [
        "IL6",
        "IL6R",
        "JAK1",
        "JAK2",
        "STAT3",
        "IL10",
        "IL2RA",
        "MAP3K8",
        "PTPN11",
        "IFNGR2",
    ],
    "IL-23": ["IL23A", "IL12B", "IL23R", "JAK2", "TYK2", "STAT3", "STAT4"],
    "JAK-STAT": [
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
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional sample stems such as L01_Ctrl L01_LPS.",
    )
    return parser.parse_args()


def decode(values) -> np.ndarray:
    return np.array([value.decode("utf-8") for value in values])


def guide_information() -> pd.DataFrame:
    path = (
        ROOT
        / "data"
        / "external"
        / "experimental"
        / "GSE327644"
        / "GSE327644_PRJNA1450051_marson_mac_perturbseq_guide_template.txt.gz"
    )
    return pd.read_csv(path, sep="\t")


def sample_metadata(path: Path) -> dict:
    match = re.search(r"_(L\d+)_Mac_(Ctrl|LPS)_", path.name)
    if not match:
        raise ValueError(f"Cannot parse sample metadata from {path.name}")
    donor, condition = match.groups()
    return {
        "sample": f"{donor}_{condition}",
        "donor": donor,
        "condition": condition,
    }


def process_sample(
    path: Path,
    guide_map: dict[str, str],
    selected_genes: set[str],
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = sample_metadata(path)
    records = []
    qc_rows = []
    with h5py.File(path, "r") as handle:
        feature_type = decode(handle["matrix/features/feature_type"][:])
        feature_id = decode(handle["matrix/features/id"][:])
        feature_name = decode(handle["matrix/features/name"][:])
        gene_mask = feature_type == "Gene Expression"
        guide_mask = feature_type == "CRISPR Guide Capture"
        gene_indices = np.where(gene_mask)[0]
        guide_indices = np.where(guide_mask)[0]
        if not np.array_equal(gene_indices, np.arange(len(gene_indices))):
            raise ValueError("Gene features are not a contiguous leading block.")
        guide_rows = feature_name[guide_indices]
        group_names = []
        for guide_id in guide_rows:
            target = guide_map.get(guide_id)
            if target == "Non-Targeting Control":
                group_names.append("NTC")
            elif target in selected_genes:
                group_names.append(guide_id)
            else:
                group_names.append(None)
        guide_to_group = {
            name: group
            for name, group in zip(guide_rows, group_names)
            if group is not None
        }
        guides_per_gene = {
            gene: [
                guide
                for guide, target in guide_map.items()
                if target == gene and guide in guide_to_group
            ]
            for gene in selected_genes
        }
        group_keys = sorted(
            {
                group
                for groups in guides_per_gene.values()
                for group in groups
            }
        )
        group_index = {group: index for index, group in enumerate(group_keys)}
        group_index["NTC"] = len(group_keys)
        n_groups = len(group_keys) + 1
        guide_feature_to_group = np.full(len(guide_rows), -1, dtype=np.int16)
        for row_index, guide_id in enumerate(guide_rows):
            group = guide_to_group.get(guide_id)
            if group is not None:
                guide_feature_to_group[row_index] = group_index[group]

        indptr = handle["matrix/indptr"][:]
        n_features, n_cells = handle["matrix/shape"][:]
        n_gene_features = len(gene_indices)
        total_counts = np.zeros((n_gene_features, n_groups), dtype=np.float64)
        assigned_cells = np.zeros(n_groups, dtype=np.int64)
        guide_high_confidence = 0
        for start in range(0, n_cells, chunk_size):
            stop = min(start + chunk_size, n_cells)
            begin_ptr = indptr[start]
            end_ptr = indptr[stop]
            indices = handle["matrix/indices"][begin_ptr:end_ptr]
            data = handle["matrix/data"][begin_ptr:end_ptr]
            lengths = np.diff(indptr[start : stop + 1])
            columns = np.repeat(np.arange(stop - start), lengths)
            matrix = csc_matrix(
                (data, (indices, columns)),
                shape=(n_features, stop - start),
            )
            guide_matrix = (
                matrix[guide_indices, :].toarray().astype(np.int16, copy=False)
            )
            if guide_matrix.shape[0] > 1:
                top_two = np.partition(guide_matrix, -2, axis=0)[-2:, :]
                second = top_two[0]
                maximum = top_two[1]
            else:
                second = np.zeros(guide_matrix.shape[1], dtype=np.int16)
                maximum = guide_matrix[0]
            max_row = guide_matrix.argmax(axis=0)
            mapped_group = guide_feature_to_group[max_row]
            confidence = (
                (maximum >= 5)
                & (second <= 0.5 * np.maximum(maximum, 1))
                & (mapped_group >= 0)
            )
            guide_high_confidence += int(confidence.sum())
            if not confidence.any():
                continue
            selected_columns = np.where(confidence)[0]
            selected_groups = mapped_group[selected_columns]
            assigned_cells += np.bincount(selected_groups, minlength=n_groups)
            group_indicator = csr_matrix(
                (
                    np.ones(len(selected_columns), dtype=np.float64),
                    (selected_groups, selected_columns),
                ),
                shape=(n_groups, stop - start),
            )
            gene_matrix = matrix[:n_gene_features, :]
            total_counts += (gene_matrix @ group_indicator.T).toarray()
            print(
                f"{metadata['sample']}: {stop}/{n_cells} cells; "
                f"high-confidence guides={guide_high_confidence}",
                flush=True,
            )
        for group, index in group_index.items():
            if group == "NTC":
                guides = [
                    guide
                    for guide, target in guide_map.items()
                    if target == "Non-Targeting Control"
                    and guide in guide_to_group
                ]
                target_gene = "NTC"
            else:
                guides = [group]
                target_gene = guide_map[group]
            qc_rows.append(
                {
                    **metadata,
                    "group": group,
                    "target_gene": target_gene,
                    "guides": ",".join(guides),
                    "cells": int(assigned_cells[index]),
                }
            )
            for row_index, gene_id in enumerate(feature_id[gene_indices]):
                records.append(
                    {
                        **metadata,
                        "group": group,
                        "target_gene": target_gene,
                        "gene_id": gene_id,
                        "gene_name": feature_name[gene_indices[row_index]],
                        "count": float(total_counts[row_index, index]),
                        "cells": int(assigned_cells[index]),
                    }
                )
    return pd.DataFrame(records), pd.DataFrame(qc_rows)


def process_and_cache(
    path: Path,
    guide_map: dict[str, str],
    selected_genes: set[str],
    chunk_size: int,
    cache_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = sample_metadata(path)
    pseudobulk_path = cache_dir / f"{metadata['sample']}_pseudobulk.csv.gz"
    qc_path = cache_dir / f"{metadata['sample']}_qc.csv"
    if pseudobulk_path.exists() and qc_path.exists():
        cached_qc = pd.read_csv(qc_path)
        if cached_qc.loc[cached_qc["group"].eq("NTC"), "cells"].sum() > 0:
            return pd.read_csv(pseudobulk_path), cached_qc
    pseudobulk, qc = process_sample(
        path,
        guide_map,
        selected_genes,
        chunk_size,
    )
    pseudobulk.to_csv(pseudobulk_path, index=False, compression="gzip")
    qc.to_csv(qc_path, index=False)
    return pseudobulk, qc


def paired_effects(
    pseudobulk: pd.DataFrame,
    target_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = pseudobulk.pivot_table(
        index=["sample", "donor", "condition", "group"],
        columns="gene_name",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    library_size = wide.sum(axis=1)
    log_cpm = np.log2(wide.div(library_size, axis=0) * 1_000_000 + 1)
    target_rows = []
    differential_rows = []
    for target in target_genes:
        guide_groups = [
            group
            for group in wide.index.get_level_values("group").unique()
            if group.startswith(target + "-")
        ]
        if not guide_groups:
            continue
        target_counts = wide[
            wide.index.get_level_values("group").isin(guide_groups)
        ].groupby(level=["sample", "donor", "condition"]).sum()
        target_log = np.log2(
            target_counts.div(target_counts.sum(axis=1), axis=0) * 1_000_000 + 1
        )
        ntc_counts = wide[
            wide.index.get_level_values("group") == "NTC"
        ].groupby(level=["sample", "donor", "condition"]).sum()
        ntc_log = np.log2(
            ntc_counts.div(ntc_counts.sum(axis=1), axis=0) * 1_000_000 + 1
        )
        paired = target_log.sub(ntc_log, axis=0)
        for state in ["Ctrl", "LPS", "Combined"]:
            subset = (
                paired
                if state == "Combined"
                else paired[paired.index.get_level_values("condition") == state]
            )
            if len(subset) == 0:
                continue
            values_matrix = subset.to_numpy(dtype=float)
            if len(subset) >= 2:
                statistics, p_values = ttest_1samp(
                    values_matrix,
                    0,
                    axis=0,
                    nan_policy="omit",
                )
            else:
                statistics = np.full(values_matrix.shape[1], np.nan)
                p_values = np.full(values_matrix.shape[1], np.nan)
            if target in subset.columns:
                target_index = subset.columns.get_loc(target)
                target_rows.append(
                    {
                        "target": target,
                        "condition": state,
                        "pairs": len(subset),
                        "mean_log2fc_vs_ntc": float(
                            np.nanmean(values_matrix[:, target_index])
                        ),
                        "median_log2fc_vs_ntc": float(
                            np.nanmedian(values_matrix[:, target_index])
                        ),
                        "p_value": float(
                            p_values[target_index]
                        ),
                    }
                )
            differential_rows.append(
                pd.DataFrame(
                    {
                        "target": target,
                        "condition": state,
                        "response_gene": subset.columns,
                        "mean_log2fc_vs_ntc": np.nanmean(
                            values_matrix, axis=0
                        ),
                        "median_log2fc_vs_ntc": np.nanmedian(
                            values_matrix, axis=0
                        ),
                        "t_statistic": statistics,
                        "p_value": p_values,
                        "pairs": len(subset),
                    }
                )
            )
    differential = pd.concat(differential_rows, ignore_index=True)
    differential["fdr"] = np.nan
    for (target, state), indices in differential.groupby(
        ["target", "condition"]
    ).groups.items():
        p_values = np.nan_to_num(
            differential.loc[indices, "p_value"].to_numpy(dtype=float),
            nan=1.0,
        )
        differential.loc[indices, "fdr"] = multipletests(
            p_values, method="fdr_bh"
        )[1]
    target_frame = pd.DataFrame(target_rows)
    target_frame["fdr"] = np.nan
    for index, row in target_frame.iterrows():
        match = differential[
            differential["target"].eq(row["target"])
            & differential["condition"].eq(row["condition"])
            & differential["response_gene"].eq(row["target"])
        ]
        if not match.empty:
            target_frame.loc[index, "fdr"] = match.iloc[0]["fdr"]
    return target_frame, differential


def pathway_effects(
    differential: pd.DataFrame,
) -> pd.DataFrame:
    rng = np.random.default_rng(20261113)
    rows = []
    for (target, state), frame in differential.groupby(["target", "condition"]):
        if state == "Combined":
            continue
        frame = frame.dropna(subset=["mean_log2fc_vs_ntc"])
        frame = frame.sort_values("mean_log2fc_vs_ntc", ascending=False)
        background = set(frame["response_gene"])
        all_values = frame["mean_log2fc_vs_ntc"].to_numpy(dtype=float)
        for pathway, genes in PATHWAYS.items():
            members = [gene for gene in genes if gene in background]
            if len(members) < 3:
                continue
            values = frame[frame["response_gene"].isin(members)][
                "mean_log2fc_vs_ntc"
            ].to_numpy(dtype=float)
            null = np.array(
                [
                    rng.choice(all_values, size=len(members), replace=False).mean()
                    for _ in range(10_000)
                ]
            )
            observed = float(values.mean())
            rows.append(
                {
                    "target": target,
                    "condition": state,
                    "pathway": pathway,
                    "members": len(members),
                    "mean_log2fc": observed,
                    "median_log2fc": float(np.median(values)),
                    "pathway_z": float(
                        (observed - null.mean()) / null.std(ddof=1)
                    ),
                    "pathway_p": float(
                        (1 + np.sum(np.abs(null - null.mean()) >= abs(observed - null.mean())))
                        / (1 + len(null))
                    ),
                }
            )
    result = pd.DataFrame(rows)
    result["pathway_fdr"] = np.nan
    for (_, _), indices in result.groupby(["target", "condition"]).groups.items():
        result.loc[indices, "pathway_fdr"] = multipletests(
            result.loc[indices, "pathway_p"], method="fdr_bh"
        )[1]
    return result


def guide_level_effects(
    pseudobulk: pd.DataFrame,
    target_genes: list[str],
) -> pd.DataFrame:
    wide = pseudobulk.pivot_table(
        index=["sample", "donor", "condition", "group"],
        columns="gene_name",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    log_cpm = np.log2(wide.div(wide.sum(axis=1), axis=0) * 1_000_000 + 1)
    ntc = log_cpm[
        log_cpm.index.get_level_values("group") == "NTC"
    ].droplevel("group")
    rows = []
    for target in target_genes:
        for guide in sorted(
            group
            for group in log_cpm.index.get_level_values("group").unique()
            if group.startswith(target + "-")
        ):
            guide_log = log_cpm[
                log_cpm.index.get_level_values("group") == guide
            ].droplevel("group")
            difference = guide_log.sub(ntc, axis=0)
            for condition in ["Ctrl", "LPS"]:
                state = difference[
                    difference.index.get_level_values("condition") == condition
                ]
                if state.empty:
                    continue
                row = {
                    "target": target,
                    "guide": guide,
                    "condition": condition,
                    "donors": int(
                        state.index.get_level_values("donor").nunique()
                    ),
                    "target_self_log2fc": float(
                        state[target].mean()
                    )
                    if target in state.columns
                    else np.nan,
                }
                for pathway, genes in PATHWAYS.items():
                    present = [gene for gene in genes if gene in state.columns]
                    row[pathway] = float(
                        state[present].mean(axis=1).mean()
                    )
                rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    configure_proxy()
    guide_dir = ROOT / "data" / "external" / "experimental" / "GSE327644"
    raw_dir = guide_dir / "raw"
    work_dir = ROOT / "work" / "gse327644_perturbseq"
    table_dir = ROOT / "results" / "tables"
    work_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    guides = guide_information()
    guide_map = dict(
        zip(
            guides["guide_id"],
            guides["guide_target_gene_name"],
        )
    )
    for guide_id in guides.loc[guides["guide_role"].eq("Control"), "guide_id"]:
        guide_map[guide_id] = "Non-Targeting Control"
    selected_genes = set(TARGET_GENES)
    files = sorted(raw_dir.glob("*.h5"))
    if args.samples:
        requested = set(args.samples)
        files = [
            path
            for path in files
            if "_".join(
                [
                    re.search(r"_(L\d+)_Mac_", path.name).group(1),
                    re.search(r"_Mac_(Ctrl|LPS)_", path.name).group(1),
                ]
            )
            in requested
        ]
    if not files:
        raise FileNotFoundError(f"No Perturb-seq H5 samples in {raw_dir}")
    pseudobulk_parts = []
    qc_parts = []
    with ThreadPoolExecutor(max_workers=min(args.workers, len(files))) as executor:
        futures = {
            executor.submit(
                process_and_cache,
                path,
                guide_map,
                selected_genes,
                args.chunk_size,
                work_dir,
            ): path
            for path in files
        }
        for future in as_completed(futures):
            pseudobulk, qc = future.result()
            pseudobulk_parts.append(pseudobulk)
            qc_parts.append(qc)
    pseudobulk = pd.concat(pseudobulk_parts, ignore_index=True)
    qc = pd.concat(qc_parts, ignore_index=True)
    pseudobulk.to_csv(
        table_dir / "perturbseq_gse327644_pseudobulk.csv.gz",
        index=False,
        compression="gzip",
    )
    qc.to_csv(
        table_dir / "perturbseq_gse327644_cell_assignment.csv",
        index=False,
    )
    target_effects, differential = paired_effects(pseudobulk, TARGET_GENES)
    target_effects.to_csv(
        table_dir / "perturbseq_gse327644_target_knockdown.csv",
        index=False,
    )
    differential.to_csv(
        table_dir / "perturbseq_gse327644_differential_effects.csv.gz",
        index=False,
        compression="gzip",
    )
    pathway = pathway_effects(differential)
    pathway.to_csv(
        table_dir / "perturbseq_gse327644_pathway_effects.csv",
        index=False,
    )
    guide_effects = guide_level_effects(pseudobulk, TARGET_GENES)
    guide_effects.to_csv(
        table_dir / "perturbseq_gse327644_guide_effects.csv",
        index=False,
    )
    de_summary = (
        differential.groupby(["target", "condition"], as_index=False)
        .agg(
            pairs=("pairs", "first"),
            genes_tested=("response_gene", "nunique"),
            fdr05_genes=("fdr", lambda values: int((values < 0.05).sum())),
            fdr05_abs05_genes=(
                "fdr",
                lambda values: int(
                    (
                        (values < 0.05)
                        & (
                            differential.loc[values.index, "mean_log2fc_vs_ntc"]
                            .abs()
                            > 0.5
                        )
                    ).sum()
                ),
            ),
        )
    )
    de_summary.to_csv(
        table_dir / "perturbseq_gse327644_differential_summary.csv",
        index=False,
    )
    summary = {
        "samples": int(pseudobulk["sample"].nunique()),
        "selected_guides": int(
            qc[~qc["group"].eq("NTC")]["group"].nunique()
        ),
        "selected_targets": sorted(
            qc.loc[~qc["group"].eq("NTC"), "target_gene"].unique().tolist()
        ),
        "total_assigned_cells": int(qc["cells"].sum()),
        "target_knockdown": target_effects.to_dict("records"),
        "pathway_effects": pathway.to_dict("records"),
        "differential_summary": de_summary.to_dict("records"),
    }
    (table_dir / "perturbseq_gse327644_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
