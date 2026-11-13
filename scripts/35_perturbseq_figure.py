from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from figure_palette import (
    AXIS,
    BURGUNDY,
    DIVERGING,
    GREEN,
    INK,
    ORANGE,
    PEACH,
    WHITE,
    configure_style,
    label_panels,
)


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    configure_style()
    table_dir = ROOT / "results" / "tables"
    figure_dir = ROOT / "results" / "figures"
    guide = pd.read_csv(table_dir / "perturbseq_gse327644_guide_effects.csv")
    knockdown = pd.read_csv(
        table_dir / "perturbseq_gse327644_target_knockdown.csv"
    )
    differential_summary = pd.read_csv(
        table_dir / "perturbseq_gse327644_differential_summary.csv"
    )
    pathway = pd.read_csv(
        table_dir / "perturbseq_gse327644_pathway_effects.csv"
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.2, 9.2),
        gridspec_kw={"width_ratios": [1.15, 1.0], "height_ratios": [1.1, 1.0]},
    )
    target_order = [
        "TNFAIP3",
        "PDGFB",
        "PRKCQ",
        "TNFSF15",
        "IL23R",
        "MAP3K8",
        "TNF",
    ]
    sns.stripplot(
        data=guide,
        x="target_self_log2fc",
        y="target",
        hue="condition",
        order=target_order,
        hue_order=["Ctrl", "LPS"],
        dodge=True,
        jitter=0.10,
        size=4.5,
        alpha=0.72,
        palette={"Ctrl": GREEN, "LPS": ORANGE},
        edgecolor=AXIS,
        linewidth=0.35,
        ax=axes[0, 0],
    )
    axes[0, 0].axvline(0, color=AXIS, linewidth=0.8)
    axes[0, 0].set_xlabel("Guide-level target transcript log2FC versus NTC")
    axes[0, 0].set_ylabel("")
    axes[0, 0].set_title(
        "Perturb-seq target engagement by guide",
        loc="left",
        fontweight="bold",
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    axes[0, 0].legend(
        handles,
        labels,
        frameon=False,
        fontsize=8,
        loc="lower right",
        title=None,
    )

    combined = differential_summary[
        differential_summary["condition"].eq("Combined")
    ].set_index("target").loc[target_order].reset_index()
    sns.barplot(
        data=combined,
        y="target",
        x="fdr05_genes",
        color=GREEN,
        ax=axes[0, 1],
    )
    for index, row in combined.iterrows():
        axes[0, 1].text(
            row["fdr05_genes"] + 25,
            index,
            f"{int(row['fdr05_genes'])} / {int(row['fdr05_abs05_genes'])}",
            va="center",
            fontsize=7,
        )
    axes[0, 1].set_xlabel("FDR < 0.05 genes (absolute log2FC > 0.5)")
    axes[0, 1].set_ylabel("")
    axes[0, 1].set_title(
        "Transcriptome-wide perturbation response",
        loc="left",
        fontweight="bold",
    )

    prkcq = guide[
        guide["target"].eq("PRKCQ") & guide["condition"].eq("LPS")
    ].copy()
    prkcq = prkcq.melt(
        id_vars=["guide"],
        value_vars=["IL-6/JAK/STAT3", "IL-23"],
        var_name="pathway",
        value_name="log2fc",
    )
    sns.barplot(
        data=prkcq,
        x="guide",
        y="log2fc",
        hue="pathway",
        palette={"IL-6/JAK/STAT3": GREEN, "IL-23": PEACH},
        edgecolor=AXIS,
        linewidth=0.4,
        ax=axes[1, 0],
    )
    axes[1, 0].axhline(0, color=AXIS, linewidth=0.8)
    axes[1, 0].set_xlabel("")
    axes[1, 0].set_ylabel("Mean pathway-gene log2FC versus NTC")
    axes[1, 0].set_title(
        "PRKCQ knockout in LPS macrophages",
        loc="left",
        fontweight="bold",
    )
    axes[1, 0].legend(frameon=False, fontsize=8, title=None)

    lps = pathway[pathway["condition"].eq("LPS")]
    heat = lps.pivot(index="target", columns="pathway", values="mean_log2fc")
    heat = heat.loc[
        [target for target in target_order if target in heat.index]
    ]
    sns.heatmap(
        heat,
        cmap=DIVERGING,
        center=0,
        annot=True,
        fmt=".2f",
        linewidths=0.6,
        linecolor=WHITE,
        cbar_kws={"label": "Mean pathway-gene log2FC"},
        ax=axes[1, 1],
    )
    axes[1, 1].set_xlabel("")
    axes[1, 1].set_ylabel("")
    axes[1, 1].set_title(
        "LPS pathway response by perturbation",
        loc="left",
        fontweight="bold",
    )
    label_panels(axes.flat, x=-0.07, y=1.03)
    fig.tight_layout()
    save(fig, figure_dir / "figure9_perturbseq_validation")
    print("Wrote figure9_perturbseq_validation")


if __name__ == "__main__":
    main()
