from __future__ import annotations

import json
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
    CONCORDANCE,
    GREEN,
    INK,
    ORANGE,
    TARGET_COLORS,
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
    threshold = pd.read_csv(
        table_dir / "perturbseq_umi_threshold_sensitivity.csv"
    )
    robustness = json.loads(
        (table_dir / "perturbseq_robustness_summary.json").read_text(
            encoding="utf-8"
        )
    )
    mixed = pd.read_csv(table_dir / "perturbseq_mixed_model_effects.csv")
    coloc = pd.read_csv(table_dir / "camk2a_coloc_sensitivity.csv")
    topics = pd.read_csv(
        table_dir / "scatac_pdgfb_variant_topic_percentiles.csv"
    )
    gsea = pd.read_csv(table_dir / "gse175685_prkcb_chsy3_gsea.csv")
    gsea_sensitivity = pd.read_csv(
        table_dir / "gse175685_pathway_confounding_sensitivity.csv"
    )
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.0))

    self_threshold = threshold[
        threshold["condition"].eq("Combined")
        & threshold["metric"].eq("self")
        & threshold["target"].isin(["PDGFB", "PRKCQ", "IL23R"])
    ]
    sns.lineplot(
        data=self_threshold,
        x="umi_threshold",
        y="mean_log2fc",
        hue="target",
        palette=TARGET_COLORS,
        marker="o",
        ax=axes[0, 0],
    )
    axes[0, 0].axhline(0, color=AXIS, linewidth=0.8)
    axes[0, 0].set_xticks([3, 5, 10])
    axes[0, 0].set_xlabel("Dominant-guide UMI threshold")
    axes[0, 0].set_ylabel("Target transcript log2FC")
    axes[0, 0].set_title(
        "Perturb-seq assignment robustness",
        loc="left",
        fontweight="bold",
    )
    axes[0, 0].legend(frameon=False, fontsize=8, title=None)

    coloc_plot = coloc.melt(
        id_vars=["sd_prior", "p12", "variants"],
        value_vars=["H3", "H4"],
        var_name="hypothesis",
        value_name="posterior",
    )
    sns.lineplot(
        data=coloc_plot,
        x="p12",
        y="posterior",
        hue="hypothesis",
        palette={"H3": GREEN, "H4": BURGUNDY},
        style="sd_prior",
        markers=True,
        ax=axes[0, 1],
    )
    axes[0, 1].set_xscale("log")
    axes[0, 1].axhline(0.8, color=ORANGE, linestyle="--", linewidth=0.8)
    axes[0, 1].set_ylim(0, 1)
    axes[0, 1].set_xlabel("Colocalization prior p12")
    axes[0, 1].set_ylabel("Posterior probability")
    axes[0, 1].set_title(
        "CAMK2A formal colocalization",
        loc="left",
        fontweight="bold",
    )

    gsea_terms = [
        "JAK-STAT signaling pathway",
        "IL-2/STAT5 Signaling",
        "IL-6/JAK/STAT3 Signaling",
        "TNF-alpha Signaling via NF-kB",
    ]
    gsea_subset = gsea[gsea["Term"].isin(gsea_terms)].copy()
    sns.barplot(
        data=gsea_subset,
        y="Term",
        x="NES",
        hue="target",
        palette=TARGET_COLORS,
        ax=axes[0, 2],
    )
    sns.scatterplot(
        data=gsea_sensitivity,
        y="Term",
        x="residualized_NES",
        hue="target",
        palette=TARGET_COLORS,
        marker="X",
        s=90,
        edgecolor=WHITE,
        linewidth=1.0,
        legend=False,
        ax=axes[0, 2],
    )
    axes[0, 2].axvline(0, color=AXIS, linewidth=0.8)
    axes[0, 2].set_xlabel("Normalized enrichment score")
    axes[0, 2].set_ylabel("")
    axes[0, 2].set_title(
        "GSEA before and after growth residualization",
        loc="left",
        fontweight="bold",
    )
    axes[0, 2].legend(frameon=False, fontsize=8, title=None)

    topics_top = topics.head(12).sort_values("percentile")
    axes[1, 0].barh(
        topics_top["topic"],
        topics_top["percentile"],
        color=GREEN,
    )
    axes[1, 0].axvline(0.95, color=ORANGE, linestyle="--", linewidth=0.8)
    axes[1, 0].set_xlim(0, 1)
    axes[1, 0].set_xlabel("Percentile among scATAC variants")
    axes[1, 0].set_ylabel("")
    axes[1, 0].set_title(
        "PDGFB variant topic specificity",
        loc="left",
        fontweight="bold",
    )

    selected_mixed = mixed[
        (
            mixed["target"].eq("PDGFB")
            & mixed["metric"].eq("self")
        )
        | (
            mixed["target"].eq("PRKCQ")
            & mixed["condition"].eq("LPS")
            & mixed["metric"].isin(["IL-23", "IL-6/JAK/STAT3"])
        )
    ].copy()
    selected_mixed["label"] = (
        selected_mixed["target"]
        + " "
        + selected_mixed["condition"]
        + " "
        + selected_mixed["metric"]
    )
    axes[1, 1].barh(
        selected_mixed["label"],
        selected_mixed["coefficient"],
        color=np.where(
            selected_mixed["coefficient"] < 0,
            GREEN,
            ORANGE,
        ),
    )
    axes[1, 1].axvline(0, color=AXIS, linewidth=0.8)
    axes[1, 1].set_xlabel("Donor-blocked model coefficient")
    axes[1, 1].set_ylabel("")
    axes[1, 1].set_title(
        "Donor-blocked model recheck",
        loc="left",
        fontweight="bold",
    )

    guide_fraction = pd.DataFrame(
        robustness["leave_one_guide_same_direction_fraction"]
    )
    donor_fraction = pd.DataFrame(
        robustness["leave_one_donor_same_direction_fraction"]
    )
    guide_fraction["analysis"] = "leave-one-guide"
    donor_fraction["analysis"] = "leave-one-donor"
    fractions = pd.concat([guide_fraction, donor_fraction])
    fractions = fractions[
        fractions["target"].isin(["PDGFB", "PRKCQ", "IL23R"])
        & fractions["metric"].isin(["self", "IL-23", "IL-6/JAK/STAT3"])
    ].copy()
    fractions["label"] = (
        fractions["target"] + " " + fractions["metric"]
    )
    fraction_table = fractions.pivot_table(
        index="label",
        columns="analysis",
        values="same_direction_as_full",
        aggfunc="mean",
    )
    sns.heatmap(
        fraction_table,
        annot=True,
        fmt=".2f",
        cmap=CONCORDANCE,
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor=WHITE,
        cbar_kws={"label": "Direction concordance"},
        ax=axes[1, 2],
    )
    axes[1, 2].set_xlabel("")
    axes[1, 2].set_ylabel("")
    axes[1, 2].set_title(
        "Leave-one-out direction stability",
        loc="left",
        fontweight="bold",
    )
    label_panels(axes.flat, x=-0.06, y=1.03)
    fig.tight_layout()
    save(fig, figure_dir / "figure10_evidence_robustness")
    print("Wrote figure10_evidence_robustness")


if __name__ == "__main__":
    main()
