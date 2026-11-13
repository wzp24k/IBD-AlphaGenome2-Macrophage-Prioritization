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
    INK,
    NEUTRAL,
    PATHWAY_COLORS,
    PEACH,
    configure_style,
    label_panels,
)


def main() -> None:
    configure_style()
    table_dir = ROOT / "results" / "tables"
    figure_dir = ROOT / "results" / "figures"
    summary = pd.read_csv(
        table_dir / "external_gwas_block_bootstrap_pathways.csv"
    )
    summary = summary.rename(
        columns={"bootstrap_fdr": "fdr"}
    )
    summary["neglog10_fdr"] = -np.log10(summary["fdr"])
    rankings = pd.read_csv(table_dir / "external_gwas_gene_rankings.csv")
    top = rankings.sort_values("external_score", ascending=False).head(25).copy()
    il2_genes = {"IL10", "SPRED2", "MUC1", "GPR65", "IL2RA", "IRF8", "MAP3K8"}
    il6_genes = {"ITGA4", "IFNGR2", "IL2RA", "MAP3K8", "PTPN11"}
    top["pathway"] = np.select(
        [
            top["target_symbol"].isin(il2_genes & il6_genes),
            top["target_symbol"].isin(il2_genes),
            top["target_symbol"].isin(il6_genes),
        ],
        ["Both", "IL-2/STAT5", "IL-6/JAK/STAT3"],
        default="Other",
    )
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.2, 4.9),
        gridspec_kw={"width_ratios": [1.05, 1.25]},
    )
    summary = summary.iloc[::-1].reset_index(drop=True)
    y_positions = np.arange(len(summary))
    colors = [
        BURGUNDY if value <= 0.05 else PEACH if value <= 0.10 else NEUTRAL
        for value in summary["fdr"]
    ]
    for y_position, row, color in zip(y_positions, summary.itertuples(), colors):
        axes[0].plot(
            [row.bootstrap_ci_low, row.bootstrap_ci_high],
            [y_position, y_position],
            color=color,
            linewidth=2.2,
        )
        axes[0].scatter(
            row.observed_auc_difference,
            y_position,
            color=color,
            edgecolor=INK,
            linewidth=0.45,
            s=42,
            zorder=3,
        )
        axes[0].text(
            0.585,
            y_position,
            f"FDR={row.fdr:.3g}",
            ha="right",
            va="center",
            fontsize=7,
            color=INK,
        )
    axes[0].axvline(0, color=AXIS, linewidth=0.9)
    axes[0].set_yticks(y_positions, summary["pathway"])
    axes[0].set_xlim(-0.03, 0.60)
    axes[0].set_xlabel("AUC difference (95% locus-block bootstrap CI)")
    axes[0].set_ylabel("")
    axes[0].set_title(
        "External-GWAS locus-block bootstrap",
        loc="left",
        fontweight="bold",
    )

    palette = PATHWAY_COLORS
    sns.barplot(
        data=top,
        y="target_symbol",
        x="external_score",
        hue="pathway",
        palette=palette,
        dodge=False,
        ax=axes[1],
    )
    axes[1].set_ylabel("")
    axes[1].set_xlabel("Primary priority x external-cluster support")
    axes[1].set_title("Top weighted target genes", loc="left", fontweight="bold")
    axes[1].legend(title="Set membership", frameon=False, fontsize=7, loc="lower right")
    label_panels(axes, x=-0.08, y=1.03)
    fig.tight_layout()
    fig.savefig(figure_dir / "figure7_all_l2g_pathways.png", dpi=300, bbox_inches="tight")
    fig.savefig(figure_dir / "figure7_all_l2g_pathways.pdf", bbox_inches="tight")
    plt.close(fig)
    print("Wrote figure7_all_l2g_pathways")


if __name__ == "__main__":
    main()
