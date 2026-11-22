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
    GREEN,
    LIME,
    NEUTRAL,
    ORANGE,
    PEACH,
    WHITE,
    configure_style,
    label_panels,
)


def save(fig, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    configure_style()
    table_dir = ROOT / "results" / "tables"
    figure_dir = ROOT / "results" / "figures"
    support = pd.read_csv(table_dir / "top10_cd_uc_target_support.csv")
    cross = pd.read_csv(table_dir / "prkcb_cross_ancestry_mr.csv")
    treatment = pd.read_csv(table_dir / "treatment_response_meta_analysis.csv")
    readiness = pd.read_csv(table_dir / "clinical_development_readiness_matrix.csv")

    support = support.sort_values("priority_rank")
    heat = support.set_index("target")[["max_l2g_crohn", "max_l2g_uc"]]
    heat.columns = ["Crohn", "UC"]
    heat = heat.fillna(0)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(10.4, 8.2),
        gridspec_kw={"width_ratios": [1.0, 1.15], "height_ratios": [1.0, 1.0]},
    )
    sns.heatmap(
        heat,
        cmap="YlGnBu",
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor=WHITE,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Maximum CD/UC L2G score", "shrink": 0.8},
        ax=axes[0, 0],
    )
    axes[0, 0].set_title("Disease-specific target support", loc="left", fontweight="bold")
    axes[0, 0].set_xlabel("")
    axes[0, 0].set_ylabel("")

    cross = cross.copy()
    cross["label"] = cross["study"].replace(
        {
            "IBD_combined": "IBD combined",
            "Crohn": "Crohn",
            "UC": "UC",
            "IBD_EAS": "IBD East Asian",
        }
    )
    colors = [BURGUNDY if study == "IBD_EAS" else NEUTRAL for study in cross["study"]]
    axes[0, 1].errorbar(
        cross["beta"],
        np.arange(len(cross)),
        xerr=1.96 * cross["se"],
        fmt="o",
        color=AXIS,
        ecolor=AXIS,
        capsize=3,
        markersize=4,
    )
    for index, (_, row) in enumerate(cross.iterrows()):
        axes[0, 1].scatter(
            row["beta"],
            index,
            s=55,
            color=colors[index],
            edgecolor=WHITE,
            zorder=3,
        )
    axes[0, 1].axvline(0, color=AXIS, lw=0.8)
    axes[0, 1].set_yticks(np.arange(len(cross)))
    axes[0, 1].set_yticklabels(cross["label"])
    axes[0, 1].set_xlabel("PRKCB cis-MR beta (log OR per SD)")
    axes[0, 1].set_title("Cross-ancestry PRKCB validation", loc="left", fontweight="bold")

    treatment = treatment.copy()
    treatment["label"] = treatment["feature"].replace(
        {
            "PRKCB_expression": "PRKCB expression",
            "PDGFB_expression": "PDGFB expression",
            "PRKCB_network_program": "PRKCB network program",
            "macrophage_marker_program": "Macrophage marker program",
        }
    )
    y = np.arange(len(treatment))
    axes[1, 0].errorbar(
        treatment["smd"],
        y,
        xerr=[
            treatment["smd"] - treatment["ci_low"],
            treatment["ci_high"] - treatment["smd"],
        ],
        fmt="o",
        color=ORANGE,
        ecolor=ORANGE,
        capsize=3,
        markersize=5,
    )
    axes[1, 0].axvline(0, color=AXIS, lw=0.8)
    axes[1, 0].set_yticks(y)
    axes[1, 0].set_yticklabels(treatment["label"])
    axes[1, 0].set_xlabel("Responder versus non-responder SMD")
    axes[1, 0].set_title("Public anti-TNF response programs", loc="left", fontweight="bold")

    readiness = readiness.sort_values("readiness_score_0_10", ascending=True)
    colors = [
        BURGUNDY
        if value >= 7
        else ORANGE
        if value >= 5
        else GREEN
        if value >= 3
        else NEUTRAL
        for value in readiness["readiness_score_0_10"]
    ]
    axes[1, 1].barh(
        readiness["target"],
        readiness["readiness_score_0_10"],
        color=colors,
        edgecolor=WHITE,
    )
    for index, value in enumerate(readiness["readiness_score_0_10"]):
        axes[1, 1].text(value + 0.12, index, f"{value:.1f}", va="center", fontsize=7)
    axes[1, 1].set_xlim(0, 8.6)
    axes[1, 1].set_xlabel("Clinical development readiness score (0-10)")
    axes[1, 1].set_title("Target readiness matrix", loc="left", fontweight="bold")
    axes[1, 1].text(
        0.98,
        0.03,
        "PRKCB rank 1: 70.1%; top 2: 90.4%",
        transform=axes[1, 1].transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color=AXIS,
    )

    for axis in axes.flat:
        axis.tick_params(labelsize=7.5)
        axis.xaxis.label.set_size(8)
        axis.yaxis.label.set_size(8)
        axis.title.set_size(9)
    label_panels(axes.flat, x=-0.06, y=1.03)
    fig.tight_layout()
    save(fig, figure_dir / "figure13_clinical_translation")


if __name__ == "__main__":
    main()
