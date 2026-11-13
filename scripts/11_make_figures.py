from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.patches import FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ROOT
from figure_palette import (
    AXIS,
    BURGUNDY,
    DIVERGING,
    GREEN,
    LIME,
    LIGHT_GREEN,
    NEUTRAL,
    ORANGE,
    PEACH,
    TIER_COLORS,
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
    lead = pd.read_csv(table_dir / "ranked_lead_variants.csv")
    all_ranked = pd.read_csv(table_dir / "ranked_variant_gene_candidates.csv")
    summary = json.loads((table_dir / "analysis_summary.json").read_text(encoding="utf-8"))

    counts = pd.DataFrame(
        {
            "stage": [
                "Credible-set variant rows",
                "Non-coding, PIP >= 0.1",
                "Independent lead loci",
                "Atlas-scored targets",
                "Live-model targets",
                "Tier 1 candidates",
            ],
            "count": [
                5192,
                summary["primary_variants"],
                summary["independent_loci"],
                summary["variants_with_atlas_target"],
                summary["variants_with_live_target"],
                summary["tier_1"],
            ],
        }
    )
    fig, ax = plt.subplots(figsize=(11.0, 7.0))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def roadmap_box(
        x: float,
        y: float,
        width: float,
        height: float,
        title: str,
        body: str,
        facecolor: str,
        textcolor: str = "#343434",
    ) -> None:
        ax.add_patch(
            Rectangle(
                (x, y),
                width,
                height,
                facecolor=facecolor,
                edgecolor=WHITE,
                linewidth=2,
            )
        )
        ax.text(
            x + width / 2,
            y + height * 0.80,
            title,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=textcolor,
        )
        ax.text(
            x + width / 2,
            y + height * 0.56,
            body,
            ha="center",
            va="top",
            fontsize=8.2,
            color=textcolor,
            linespacing=1.35,
        )

    def flow(start, end) -> None:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.4,
                color=AXIS,
            )
        )

    roadmap_box(
        28,
        83,
        44,
        12,
        "1  Input and credible-set filtering",
        "5,192 credible-set rows -> 395 non-coding variants\n281 independent loci",
        LIGHT_GREEN,
    )
    roadmap_box(
        6,
        57,
        39,
        20,
        "2  Regulatory prediction",
        "Open Targets L2G\nAlphaGenome Atlas quantiles\nLive colon-macrophage scoring",
        LIME,
    )
    roadmap_box(
        55,
        57,
        39,
        20,
        "2  Macrophage context",
        "Macrophage versus T-cell and epithelial controls\nDonor-derived macrophages\nActive-disease pools",
        PEACH,
    )
    flow((50, 83), (25.5, 77))
    flow((50, 83), (74.5, 77))
    roadmap_box(
        18,
        32,
        64,
        15,
        "3  Independent evidence",
        "UKB 2025 + MVP 2024; locus-block bootstrap\nPublic eQTL, scATAC, CRISPR, Perturb-seq and overexpression",
        ORANGE,
    )
    flow((25.5, 57), (42, 47))
    flow((74.5, 57), (58, 47))
    roadmap_box(
        20,
        6,
        60,
        18,
        "4  Calibrated output",
        "10 Tier 1 and 57 Tier 2 loci\nFive targets lack supportive experimental evidence\nClaims bounded by evidence level",
        BURGUNDY,
        textcolor=WHITE,
    )
    flow((50, 32), (50, 24))
    ax.set_title(
        "Technical roadmap for macrophage-focused variant prioritization",
        loc="left",
        fontweight="bold",
        pad=10,
    )
    save(fig, figure_dir / "figure1_analysis_flow")

    scatter = lead.dropna(subset=["atlas_gi_max_abs_q", "live_specificity"]).copy()
    palette = TIER_COLORS
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    sns.scatterplot(
        data=scatter,
        x="atlas_gi_max_abs_q",
        y="live_specificity",
        hue="evidence_tier",
        palette=palette,
        size="priority_score",
        sizes=(25, 170),
        alpha=0.75,
        linewidth=0.3,
        edgecolor=WHITE,
        ax=ax,
    )
    ax.axhline(0, color=AXIS, lw=0.8, ls="--")
    ax.axhline(0.25, color=BURGUNDY, lw=0.8, ls=":")
    ax.axvline(0.95, color=BURGUNDY, lw=0.8, ls=":")
    texts = []
    for row in scatter.head(12).itertuples():
        texts.append(
            ax.text(
                row.atlas_gi_max_abs_q,
                row.live_specificity,
                row.l2g_target_symbol,
                fontsize=7,
            )
        )
    adjust_text(
        texts,
        arrowprops={"arrowstyle": "-", "color": AXIS, "lw": 0.5},
        ax=ax,
    )
    ax.set_xlabel("AlphaGenome Atlas maximum absolute RNA-seq quantile")
    ax.set_ylabel("Macrophage-versus-control specificity percentile delta")
    ax.set_title("Genome-calibrated effect versus macrophage selectivity", loc="left", fontweight="bold")
    ax.legend(
        title="Evidence",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=7,
    )
    save(fig, figure_dir / "figure2_quantile_vs_specificity")

    heat_columns = [
        "pip_percentile",
        "l2g_percentile",
        "atlas_gi_max_abs_q",
        "rna_specificity",
        "sc_specificity",
        "custom_macrophage_abs_percentile",
        "priority_score",
    ]
    labels = [
        "PIP",
        "L2G",
        "Atlas RNA-q",
        "RNA specificity",
        "scRNA specificity",
        "Patient macrophage",
        "Priority",
    ]
    heat = (
        lead.head(25)
        .set_index("l2g_target_symbol")[heat_columns]
        .apply(pd.to_numeric, errors="coerce")
    )
    heat.index = [f"{gene} ({variant.split(':')[0]})" for gene, variant in zip(heat.index, lead.head(25)["canonical_variant"])]
    fig, ax = plt.subplots(figsize=(7.6, 7.2))
    sns.heatmap(
        heat,
        cmap=DIVERGING,
        center=0,
        vmin=-0.5,
        vmax=1,
        linewidths=0.5,
        linecolor=WHITE,
        cbar_kws={"label": "Scaled evidence or percentile"},
        ax=ax,
    )
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("")
    ax.set_title("Top prioritized non-coding loci", loc="left", fontweight="bold")
    save(fig, figure_dir / "figure3_top_candidate_heatmap")

    ets2 = all_ranked  # Keep the plotting function independent of candidate status.
    live_path = table_dir / "live_condition_gene_scores.parquet"
    live = pd.read_parquet(live_path)
    ets2_live = live[
        (live["variant_id"].astype(str) == "chr21:39094644:G>A")
        & (live["gene_name"] == "ETS2")
        & live["imputation_mode"].isin(["observed", "imputed_supported"])
    ]
    ets2_summary = (
        ets2_live.groupby(["output_type", "biosample_name"], as_index=False)
        .agg(raw_score=("raw_score", "mean"))
    )
    condition_order = [
        "colon macrophage",
        "colon T cell",
        "colon epithelial cell",
        "colon colonic epithelial cell",
        "ileum macrophage",
        "ileum T cell",
        "ileum epithelial cell",
    ]
    ets2_summary["biosample_name"] = pd.Categorical(
        ets2_summary["biosample_name"], condition_order, ordered=True
    )
    ets2_summary = ets2_summary.sort_values(["output_type", "biosample_name"])
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=True)
    for ax, output_type in zip(axes, ["RNA_SEQ", "SC_GENE_EXPRESSION"]):
        subset = ets2_summary[ets2_summary["output_type"] == output_type]
        sns.barplot(
            data=subset,
            y="biosample_name",
            x="raw_score",
            color=GREEN,
            ax=ax,
        )
        ax.axvline(0, color=AXIS, lw=0.8)
        ax.set_title(output_type.replace("_", " "), fontsize=9)
        ax.set_xlabel("Predicted change in ETS2 expression")
        ax.set_ylabel("")
    custom = pd.read_csv(table_dir / "positive_control_custom_ets2.csv")
    custom_value = float(custom["raw_score"].mean())
    fig.suptitle(
        f"Positive control rs2836882 / ETS2; patient macrophage score = {custom_value:.2f}",
        x=0.02,
        ha="left",
        fontweight="bold",
        fontsize=10,
    )
    label_panels(axes, x=-0.08, y=1.03)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, figure_dir / "figure4_ets2_positive_control")

    background_path = table_dir / "background_live_condition_gene_scores.parquet"
    if background_path.exists():
        background = pd.read_parquet(background_path)
        # Cell-conditioning-only DENOISED outputs have an empty imputation label.
        bg_group = background.groupby(
            ["variant_id", "gene_name", "condition_role"], as_index=False
        ).agg(raw_score=("raw_score", "mean"))
        bg_pivot = bg_group.pivot_table(
            index=["variant_id", "gene_name"],
            columns="condition_role",
            values="raw_score",
        ).reset_index()
        controls = [
            column
            for column in ["immune_control", "epithelial_control", "epithelial_sensitivity"]
            if column in bg_pivot.columns
        ]
        bg_pivot["macrophage_abs"] = bg_pivot["primary_target"].abs()
        bg_pivot["control_abs_max"] = bg_pivot[controls].abs().max(axis=1)
        values = np.sort(bg_pivot["macrophage_abs"].dropna().to_numpy())
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(10.2, 4.2),
            gridspec_kw={"width_ratios": [1.25, 1]},
        )
        ecdf = np.arange(1, len(values) + 1) / len(values)
        axes[0].step(values, ecdf, where="post", color=GREEN, lw=1.8)
        median = float(np.median(values))
        p95 = float(np.quantile(values, 0.95))
        maximum = float(values.max())
        for value, label in [
            (median, f"median {median:.4f}"),
            (p95, f"95th {p95:.4f}"),
            (maximum, f"max {maximum:.4f}"),
        ]:
            axes[0].axvline(value, color=AXIS, ls=":", lw=0.8)
            axes[0].text(
                value,
                0.08,
                label,
                rotation=90,
                va="bottom",
                ha="right",
                fontsize=7,
            )
        axes[0].set_xlim(0, maximum * 1.08)
        axes[0].set_ylim(0, 1.02)
        axes[0].set_xlabel("Absolute predicted macrophage effect")
        axes[0].set_ylabel("Cumulative fraction of background variants")
        axes[0].set_title("Frequency-matched background", loc="left", fontweight="bold")

        comparison = pd.DataFrame(
            [
                ("Background median", median, LIGHT_GREEN),
                ("Background 95th", p95, GREEN),
                ("Background max", maximum, NEUTRAL),
                ("ETS2 / reference macrophage", 0.146828, ORANGE),
                ("ETS2 / donor macrophage", 0.306536, BURGUNDY),
            ],
            columns=["label", "value", "color"],
        )
        axes[1].barh(
            comparison["label"],
            comparison["value"],
            color=comparison["color"],
        )
        axes[1].set_xscale("log")
        axes[1].set_xlim(1e-4, 1.0)
        for row_index, value in enumerate(comparison["value"]):
            axes[1].text(
                value * 1.12,
                row_index,
                f"{value:.4g}",
                va="center",
                fontsize=7,
            )
        axes[1].set_xlabel("Absolute effect, log scale")
        axes[1].set_ylabel("")
        axes[1].set_title("Positive controls versus background", loc="left", fontweight="bold")
        label_panels(axes, x=-0.08, y=1.03)
        fig.tight_layout()
        save(fig, figure_dir / "figure5_ets2_background_null")

    print(f"Figures written to {figure_dir}")


if __name__ == "__main__":
    main()
