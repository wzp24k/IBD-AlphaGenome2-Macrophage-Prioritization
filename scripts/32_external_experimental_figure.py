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
    COVERAGE,
    GREEN,
    INK,
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

    evidence = pd.read_csv(table_dir / "external_experimental_evidence.csv")
    gtex = pd.read_csv(
        table_dir / "experimental_gtex_prkcb_eqtl.csv"
    ).iloc[0]
    mpra = pd.read_csv(
        table_dir / "experimental_ets2_mpra_context_effects.csv"
    )
    guides = pd.read_csv(
        table_dir / "experimental_primary_macrophage_crispr_guides.csv"
    )
    engagement = pd.read_csv(
        table_dir / "experimental_thp1_overexpression_engagement.csv"
    )
    summary = json.loads(
        (table_dir / "experimental_evidence_summary.json").read_text(
            encoding="utf-8"
        )
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.2, 9.2),
        gridspec_kw={"width_ratios": [1.0, 1.15], "height_ratios": [1.25, 1.0]},
    )

    coverage = evidence.set_index("target")[
        [
            "variant_level_experimental_sources",
            "gene_level_experimental_sources",
        ]
    ].copy()
    coverage_matrix = pd.DataFrame(
        {
            "GTEx colon\neQTL": coverage[
                "variant_level_experimental_sources"
            ]
            .fillna("")
            .str.contains("GTEx")
            .astype(int),
            "Gut scATAC\ncCRE": coverage[
                "variant_level_experimental_sources"
            ]
            .fillna("")
            .str.contains("GSE273194")
            .astype(int),
            "Primary macrophage\nCRISPR": coverage[
                "gene_level_experimental_sources"
            ]
            .fillna("")
            .str.contains("GSE327124")
            .astype(int),
            "THP-1\noverexpression": coverage[
                "gene_level_experimental_sources"
            ]
            .fillna("")
            .str.contains("GSE175685")
            .astype(int),
        }
    )
    coverage_matrix = coverage_matrix.loc[
        evidence.sort_values("lead_rank")["target"]
    ]
    sns.heatmap(
        coverage_matrix,
        cmap=COVERAGE,
        cbar=False,
        linewidths=0.8,
        linecolor=WHITE,
        vmin=0,
        vmax=1,
        ax=axes[0, 0],
    )
    axes[0, 0].set_xticklabels(
        coverage_matrix.columns, rotation=25, ha="right"
    )
    axes[0, 0].set_ylabel("")
    axes[0, 0].set_title(
        "Public evidence coverage of Tier 1 targets",
        loc="left",
        fontweight="bold",
    )
    axes[0, 0].text(
        0.02,
        -0.31,
        f"Variant-level public evidence: {summary['tier1_variant_level_experimental']}/10; "
        f"gene-level perturbation: {summary['tier1_gene_level_experimental']}/10; "
        f"gut cCRE: {summary['tier1_scatac_ccre']}/10",
        transform=axes[0, 0].transAxes,
        fontsize=8,
        color=INK,
    )

    variant_labels = ["PRKCB\nGTEx eQTL", "ETS2 control\nmacrophage MPRA"]
    variant_effects = [
        float(gtex["nes"]),
        float(summary["ets2_mpra"]["mean_risk_allele_log2_effect"]),
    ]
    mpra_range = np.array(
        [
            row.mean_log2_risk_vs_nonrisk
            for row in mpra.itertuples()
        ]
    )
    variant_errors = [
        0,
        [variant_effects[1] - mpra_range.min(), mpra_range.max() - variant_effects[1]],
    ]
    colors = [GREEN, ORANGE]
    axes[0, 1].bar(variant_labels, variant_effects, color=colors, width=0.58)
    axes[0, 1].errorbar(
        [1],
        [variant_effects[1]],
        yerr=np.array([[variant_errors[1][0]], [variant_errors[1][1]]]),
        fmt="none",
        color=AXIS,
        capsize=4,
        linewidth=1.0,
    )
    axes[0, 1].axhline(0, color=AXIS, linewidth=0.8)
    axes[0, 1].set_ylabel("Normalized effect")
    axes[0, 1].set_title(
        "Variant-level public evidence",
        loc="left",
        fontweight="bold",
    )
    axes[0, 1].text(
        0.02,
        0.94,
        f"GTEx P<0.001; NES={float(gtex['nes']):.2f}\n"
        "MPRA risk-vs-nonrisk in 2/3 contexts positive; mean +0.17 log2",
        transform=axes[0, 1].transAxes,
        va="top",
        fontsize=8,
        color=INK,
    )

    selected = ["PDGFB", "PRKCQ", "IL23R", "MAP3K8", "TNFAIP3", "TNF"]
    tnf_guides = guides[
        guides["screen"].eq("TNF") & guides["Gene"].isin(selected)
    ].copy()
    gene_order = (
        tnf_guides.groupby("Gene")["guide_score"]
        .median()
        .sort_values()
        .index.tolist()
    )
    sns.stripplot(
        data=tnf_guides,
        x="guide_score",
        y="Gene",
        order=gene_order,
        color=GREEN,
        alpha=0.65,
        size=5,
        jitter=0.12,
        edgecolor=AXIS,
        linewidth=0.35,
        ax=axes[1, 0],
    )
    medians = tnf_guides.groupby("Gene")["guide_score"].median()
    for index, gene in enumerate(gene_order):
        axes[1, 0].scatter(
            medians[gene],
            index,
            color=BURGUNDY,
            edgecolor=INK,
            linewidth=0.4,
            s=36,
            zorder=4,
        )
    axes[1, 0].axvline(0, color=AXIS, linewidth=0.8)
    axes[1, 0].set_xlabel("TNF screen log2(high/low)")
    axes[1, 0].set_ylabel("")
    axes[1, 0].set_title(
        "Primary human macrophage CRISPR screen",
        loc="left",
        fontweight="bold",
    )

    engagement = engagement.sort_values("self_expression_delta")
    axes[1, 1].barh(
        engagement["gene"],
        engagement["self_expression_delta"],
        color=BURGUNDY,
        alpha=0.85,
    )
    axes[1, 1].axvline(0, color=AXIS, linewidth=0.8)
    axes[1, 1].set_xlabel("Target-gene expression delta versus empty vector")
    axes[1, 1].set_ylabel("")
    axes[1, 1].set_title(
        "THP-1 overexpression target engagement",
        loc="left",
        fontweight="bold",
    )
    for index, row in enumerate(engagement.itertuples()):
        if row.gene in {"PRKCB", "CHSY3"}:
            axes[1, 1].text(
                row.self_expression_delta + 0.08,
                index,
                "Tier 1 target",
                va="center",
                fontsize=7,
                color=INK,
            )

    label_panels(axes.flat, x=-0.07, y=1.03)
    fig.tight_layout()
    save(fig, figure_dir / "figure8_external_experimental_evidence")
    print("Wrote figure8_external_experimental_evidence")


if __name__ == "__main__":
    main()
