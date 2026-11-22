from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from figure_palette import (
    BURGUNDY,
    GREEN,
    INK,
    LIGHT_GREEN,
    LIME,
    NEUTRAL,
    ORANGE,
    PEACH,
    WHITE,
    configure_style,
)


MODEL_COLORS = {
    "PIP-only": NEUTRAL,
    "L2G-only": GREEN,
    "Atlas-only": LIME,
    "Live-only": ORANGE,
    "Combined": BURGUNDY,
}
MODEL_OFFSETS = {
    "PIP-only": -0.16,
    "L2G-only": -0.08,
    "Atlas-only": 0.0,
    "Live-only": 0.08,
    "Combined": 0.16,
}


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def forest_panel(axis, metrics: pd.DataFrame, metric: str, title: str) -> None:
    tasks = [
        "External genetic support, both cohorts",
        "Exact-variant eQTL, replicated",
    ]
    y_positions = np.arange(len(tasks))[::-1]
    for y, task in zip(y_positions, tasks):
        subset = metrics[metrics["task"].eq(task)]
        for model in MODEL_COLORS:
            row = subset[subset["model"].eq(model)]
            if row.empty:
                continue
            row = row.iloc[0]
            x = row[metric]
            if metric == "auroc":
                low = row[f"{metric}_ci_low"]
                high = row[f"{metric}_ci_high"]
            else:
                low = row[f"{metric}_ci_low"]
                high = row[f"{metric}_ci_high"]
            axis.errorbar(
                x,
                y + MODEL_OFFSETS[model],
                xerr=[[x - low], [high - x]],
                fmt="o",
                color=MODEL_COLORS[model],
                markersize=5,
                capsize=2,
                linewidth=1.2,
                label=model if y == y_positions[0] else None,
            )
    axis.set_yticks(y_positions)
    axis.set_yticklabels(["External support\nin both cohorts", "Replicated\nexact-variant eQTL"])
    axis.set_xlabel(title)
    axis.set_xlim(0, 1)
    axis.grid(axis="y", visible=False)
    axis.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.42),
        ncol=5,
        frameon=False,
        fontsize=7,
        handletextpad=0.3,
        columnspacing=0.7,
    )


def pr_panel(
    axis,
    frame: pd.DataFrame,
    label_column: str,
    title: str,
    score_columns: dict[str, str],
) -> None:
    labels = frame[label_column].astype(int).to_numpy()
    for model, column in score_columns.items():
        precision, recall, _ = precision_recall_curve(
            labels,
            frame[column].fillna(0).to_numpy(),
        )
        axis.plot(
            recall,
            precision,
            color=MODEL_COLORS[model],
            linewidth=1.8,
            label=model,
        )
    prevalence = labels.mean()
    axis.axhline(
        prevalence,
        color=INK,
        linestyle=":",
        linewidth=1,
    )
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("Recall")
    axis.set_ylabel("Precision")
    axis.set_title(title, loc="left", fontweight="bold")
    axis.legend(frameon=False, fontsize=7, loc="upper right")


def top10_panel(axis, triangulation: pd.DataFrame, coloc: pd.DataFrame) -> None:
    matched = int(
        triangulation.loc[
            triangulation["query_group"].eq("top10_matched"),
            "any_nominal",
        ].sum()
    )
    mismatched = int(
        triangulation.loc[
            triangulation["query_group"].eq("mismatched_target_control"),
            "any_nominal",
        ].sum()
    )
    axis.bar(
        [0, 1],
        [matched, mismatched],
        color=[BURGUNDY, NEUTRAL],
        width=0.58,
    )
    axis.set_xticks([0, 1])
    axis.set_xticklabels(["Matched target", "Mismatched control"])
    axis.set_ylabel("Top-10 pairs with exact eQTL")
    axis.set_ylim(0, 10)
    axis.text(
        0,
        matched + 0.25,
        f"{matched}/10",
        ha="center",
        fontsize=13,
        fontweight="bold",
    )
    axis.text(
        1,
        mismatched + 0.25,
        f"{mismatched}/10",
        ha="center",
        fontsize=13,
        fontweight="bold",
    )
    supported = int(coloc["coloc_supported"].fillna(False).astype(bool).sum())
    axis.text(
        0.5,
        8.8,
        f"Formal colocalization H4 >= 0.8: {supported}/{len(coloc)}",
        ha="center",
        fontsize=9,
        color=INK,
    )
    axis.text(
        0.5,
        8.15,
        "Fisher P = 0.0031",
        ha="center",
        fontsize=9,
        color=INK,
    )
    axis.set_title(
        "Exact-variant eQTL and negative control",
        loc="left",
        fontweight="bold",
    )


def translation_panel(axis, frame: pd.DataFrame) -> None:
    ordered = frame.sort_values("priority_rank").copy()
    columns = [
        "external_cluster_count",
        "exact_variant_eqtl",
        "formal_colocalization",
        "experimental_support",
        "tractable_modalities",
        "ibd_indication_level_candidate",
        "directionally_compatible_ibd_drug",
        "clinical_stage",
    ]
    matrix = []
    for row in ordered.itertuples():
        stage = {
            "APPROVAL": 4,
            "PHASE_4": 4,
            "PHASE_3": 3,
            "PHASE_2": 2,
            "PHASE_1": 1,
            "WITHDRAWAL": 0,
            "": 0,
        }.get(str(row.clinical_stage), 0)
        matrix.append(
            [
                min(float(row.external_cluster_count), 2),
                1 if row.exact_variant_eqtl else 0,
                1 if row.formal_colocalization else 0,
                1 if pd.notna(row.experimental_support) else 0,
                1 if row.tractable_modalities else 0,
                1 if row.ibd_indication_level_candidate else 0,
                1 if row.directionally_compatible_ibd_drug else 0,
                stage,
            ]
        )
    matrix = np.asarray(matrix, dtype=float)
    cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
        "translation",
        [WHITE, LIME, PEACH, ORANGE, BURGUNDY],
    )
    axis.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=4)
    axis.set_xticks(np.arange(len(columns)))
    axis.set_xticklabels(
        [
            "External",
            "eQTL",
            "Coloc",
            "Experiment",
            "Tractable",
            "IBD drug",
            "Direction",
            "Clinical",
        ],
        rotation=35,
        ha="right",
        fontsize=7,
    )
    axis.set_yticks(np.arange(len(ordered)))
    axis.set_yticklabels(ordered["target"], fontsize=7)
    for row_index in range(matrix.shape[0]):
        for col_index in range(matrix.shape[1]):
            value = int(matrix[row_index, col_index])
            axis.text(
                col_index,
                row_index,
                str(value),
                ha="center",
                va="center",
                fontsize=6,
                color=WHITE if value >= 3 else INK,
            )
    axis.set_title(
        "Clinical-stage context (not a treatment recommendation)",
        loc="left",
        fontweight="bold",
        fontsize=9,
    )
    axis.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(ordered), 1), minor=True)
    axis.grid(which="minor", color=WHITE, linewidth=2)
    axis.tick_params(which="minor", length=0)


def main() -> None:
    configure_style()
    table_dir = ROOT / "results" / "tables"
    figure_dir = ROOT / "results" / "figures"
    benchmark = pd.read_csv(table_dir / "benchmark_pair_level_dataset.csv")
    metrics = pd.read_csv(table_dir / "benchmark_model_metrics.csv")
    triangulation = pd.read_csv(
        table_dir / "external_eqtl_top10_triangulation.csv"
    )
    coloc = pd.read_csv(table_dir / "top10_eqtl_coloc_results.csv")
    translation = pd.read_csv(
        table_dir / "top10_clinical_translation_annotation.csv"
    )

    fig = plt.figure(figsize=(14.0, 8.8))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0])
    forest_auc = fig.add_subplot(grid[0, 0])
    forest_ap = fig.add_subplot(grid[0, 1])
    top10_axis = fig.add_subplot(grid[0, 2])
    pr_external = fig.add_subplot(grid[1, 0])
    pr_eqtl = fig.add_subplot(grid[1, 1])
    translation_axis = fig.add_subplot(grid[1, 2])

    forest_panel(forest_auc, metrics, "auroc", "AUROC")
    forest_panel(forest_ap, metrics, "auprc", "AUPRC")
    top10_panel(top10_axis, triangulation, coloc)
    pr_panel(
        pr_external,
        benchmark,
        "external_both",
        "External support in both cohorts",
        {
            "PIP-only": "pip",
            "L2G-only": "l2g_score",
            "Atlas-only": "atlas_global_max_abs_q",
            "Live-only": "live_specificity",
            "Combined": "framework_pair_score",
        },
    )
    pr_panel(
        pr_eqtl,
        benchmark,
        "eqtl_replicated",
        "Replicated exact-variant eQTL",
        {
            "PIP-only": "pip",
            "L2G-only": "l2g_score",
            "Atlas-only": "atlas_global_max_abs_q",
            "Live-only": "live_specificity",
            "Combined": "framework_pair_score",
        },
    )
    translation_panel(translation_axis, translation)
    for axis, label in zip(
        [
            forest_auc,
            forest_ap,
            top10_axis,
            pr_external,
            pr_eqtl,
            translation_axis,
        ],
        list("ABCDEF"),
    ):
        axis.text(
            -0.12,
            1.05,
            label,
            transform=axis.transAxes,
            fontsize=11,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
    fig.tight_layout()
    save(fig, figure_dir / "figure11_benchmark_translation")
    print("Wrote figure11_benchmark_translation")


if __name__ == "__main__":
    main()
