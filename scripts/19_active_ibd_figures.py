from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from figure_palette import (
    AXIS,
    BURGUNDY,
    DIVERGING,
    GREEN,
    INK,
    NEUTRAL,
    ORANGE,
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
    figure_dir.mkdir(parents=True, exist_ok=True)
    lead = pd.read_csv(table_dir / "ranked_lead_variants_active_ibd.csv")
    repeated = pd.read_csv(table_dir / "prkcq_repeated_pool_scores.csv")
    repeated_summary = pd.read_csv(table_dir / "prkcq_repeated_pool_summary.csv")
    top = lead.head(20).copy()
    columns = [
        "custom_macrophage_effect",
        "oliver_healthy_effect",
        "oliver_active_effect",
        "crohn_healthy_effect",
        "crohn_active_effect",
        "garrido_healthy_effect",
        "garrido_active_effect",
    ]
    labels = [
        "Healthy\natlas",
        "Oliver\nhealthy",
        "Oliver\nactive",
        "Crohn\nhealthy",
        "Crohn\nactive",
        "IBD myeloid\nnormal",
        "IBD myeloid\nactive",
    ]
    heat = top.set_index("l2g_target_symbol")[columns].apply(
        pd.to_numeric, errors="coerce"
    )
    heat.index = [
        f"{gene} ({rank})"
        for gene, rank in zip(heat.index, top["lead_rank"])
    ]
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(15.6, 6.5),
        gridspec_kw={"width_ratios": [1.55, 1.05, 0.92]},
    )
    sns.heatmap(
        heat,
        cmap=DIVERGING,
        center=0,
        linewidths=0.5,
        linecolor=WHITE,
        cbar_kws={
            "label": "Predicted target-gene expression effect",
            "shrink": 0.82,
        },
        ax=axes[0],
    )
    axes[0].set_xticklabels(labels, rotation=35, ha="right")
    axes[0].set_ylabel("")
    axes[0].set_title("Top loci across healthy and active IBD macrophages", loc="left", fontweight="bold")

    scatter = lead.head(28).dropna(
        subset=["custom_macrophage_effect", "oliver_active_effect"]
    )
    palette = {True: BURGUNDY, False: NEUTRAL}
    sns.scatterplot(
        data=scatter,
        x="custom_macrophage_effect",
        y="oliver_active_effect",
        hue="active_disease_multicohort_supported",
        palette=palette,
        s=55,
        edgecolor=WHITE,
        linewidth=0.4,
        ax=axes[1],
    )
    axes[1].axhline(0, color=AXIS, lw=0.8)
    axes[1].axvline(0, color=AXIS, lw=0.8)
    texts = []
    for row in scatter.head(12).itertuples():
        texts.append(
            axes[1].text(
                row.custom_macrophage_effect,
                row.oliver_active_effect,
                row.l2g_target_symbol,
                fontsize=7,
            )
        )
    adjust_text(
        texts,
        arrowprops={"arrowstyle": "-", "color": AXIS, "lw": 0.5},
        ax=axes[1],
    )
    axes[1].set_xlabel("Healthy macrophage effect")
    axes[1].set_ylabel("Active colitis macrophage effect")
    axes[1].set_title("Direction stability in active disease", loc="left", fontweight="bold")
    axes[1].legend(
        title="Multi-cohort support",
        frameon=False,
        fontsize=8,
        loc="best",
    )

    cohort_labels = {
        "Oliver": "Oliver",
        "Crohn2023": "Kong",
        "Garrido2023": "Garrido",
    }
    state_labels = {"healthy": "healthy", "active": "active"}
    repeated["condition_label"] = repeated.apply(
        lambda row: f"{cohort_labels[row['cohort']]}\n{state_labels[row['disease_state']]}",
        axis=1,
    )
    repeated_summary["condition_label"] = repeated_summary.apply(
        lambda row: f"{cohort_labels[row['cohort']]}\n{state_labels[row['disease_state']]}",
        axis=1,
    )
    condition_order = [
        f"{cohort_labels[cohort]}\n{state_labels[state]}"
        for cohort in ["Oliver", "Crohn2023", "Garrido2023"]
        for state in ["healthy", "active"]
    ]
    state_palette = {"healthy": GREEN, "active": ORANGE}
    sns.stripplot(
        data=repeated,
        x="condition_label",
        y="prkcq_effect",
        hue="disease_state",
        order=condition_order,
        hue_order=["healthy", "active"],
        palette=state_palette,
        dodge=False,
        jitter=0.14,
        alpha=0.62,
        size=4,
        edgecolor=WHITE,
        linewidth=0.3,
        ax=axes[2],
    )
    summary_lookup = repeated_summary.set_index("condition_label")
    for index, condition in enumerate(condition_order):
        row = summary_lookup.loc[condition]
        color = state_palette[row["disease_state"]]
        axes[2].vlines(
            index,
            row["q025"],
            row["q975"],
            color=color,
            linewidth=1.8,
            zorder=3,
        )
        axes[2].scatter(
            index,
            row["median_effect"],
            color=color,
            edgecolor=INK,
            linewidth=0.5,
            s=36,
            zorder=4,
        )
    axes[2].axhline(0, color=AXIS, linewidth=0.8)
    axes[2].axhline(
        0.005,
        color=ORANGE,
        linewidth=0.8,
        linestyle=":",
    )
    axes[2].axhline(
        -0.005,
        color=ORANGE,
        linewidth=0.8,
        linestyle=":",
    )
    axes[2].set_xlabel("")
    axes[2].set_ylabel("PRKCQ single-cell expression effect")
    axes[2].set_title(
        "PRKCQ repeated-pool stability",
        loc="left",
        fontweight="bold",
    )
    axes[2].text(
        0.03,
        0.04,
        "Attenuated in active disease",
        transform=axes[2].transAxes,
        fontsize=8,
        color=INK,
    )
    handles, labels = axes[2].get_legend_handles_labels()
    axes[2].legend(
        handles,
        labels,
        frameon=False,
        fontsize=8,
        loc="upper right",
        title=None,
    )
    label_panels(axes, x=-0.06, y=1.03)
    fig.tight_layout()
    save(fig, figure_dir / "figure6_active_ibd_validation")
    print(f"Wrote {figure_dir / 'figure6_active_ibd_validation.png'}")


if __name__ == "__main__":
    main()
