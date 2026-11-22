from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    STATE_COLORS,
    WHITE,
    configure_style,
    label_panels,
)


METHOD_LABELS = {
    "IVW_random_effects": "IVW random effects",
    "weighted_median": "Weighted median",
    "MR_Egger": "MR-Egger",
}
METHOD_COLORS = {
    "IVW_random_effects": BURGUNDY,
    "weighted_median": GREEN,
    "MR_Egger": ORANGE,
}
TARGET_COLORS = {
    "PRKCB": BURGUNDY,
    "CHSY3": ORANGE,
    "ETS2": GREEN,
    "JAK2": PEACH,
    "STAT3": LIGHT_GREEN,
}


def mr_forest_panel(axis: plt.Axes, estimates: pd.DataFrame) -> None:
    primary = estimates[estimates["clump_r2_threshold"].eq(0.01)].copy()
    methods = ["IVW_random_effects", "weighted_median", "MR_Egger"]
    rows = primary.set_index("method").loc[methods]
    y = np.arange(len(rows))[::-1]
    for position, (method, row) in zip(y, rows.iterrows()):
        if not np.isfinite(row["estimate"]):
            axis.scatter([0], [position], marker="x", color=NEUTRAL, s=35)
            axis.text(
                0.004,
                position,
                "not estimable with 2 instruments",
                va="center",
                fontsize=7,
                color=NEUTRAL,
            )
            continue
        low = row["estimate"] - 1.96 * row["se"]
        high = row["estimate"] + 1.96 * row["se"]
        axis.errorbar(
            row["estimate"],
            position,
            xerr=[[row["estimate"] - low], [high - row["estimate"]]],
            fmt="o",
            color=METHOD_COLORS[method],
            markersize=6,
            capsize=3,
            linewidth=1.5,
        )
        axis.text(
            high + 0.008,
            position,
            f"P={row['p']:.2g}",
            va="center",
            fontsize=7,
        )
    axis.axvline(0, color=INK, linestyle="--", linewidth=1)
    axis.set_yticks(y)
    axis.set_yticklabels([METHOD_LABELS[method] for method in methods])
    axis.set_xlim(-0.22, 0.14)
    axis.set_xlabel("MR causal estimate (log OR per eQTL SD)")
    axis.set_title(
        "PRKCB two-sample MR\n3 approximately independent instruments",
        loc="left",
        fontweight="bold",
        fontsize=9,
    )
    axis.grid(axis="y", visible=False)


def regional_panel(
    axis: plt.Axes,
    variants: pd.DataFrame,
    ld: pd.DataFrame,
    smr: pd.DataFrame,
) -> None:
    frame = variants.merge(
        ld[["rsid", "r2_with_rs7404003"]],
        on="rsid",
        how="left",
    )
    scatter = axis.scatter(
        frame["eql_beta"],
        frame["gwas_beta"],
        c=frame["r2_with_rs7404003"],
        cmap="RdYlGn_r",
        vmin=0,
        vmax=1,
        s=8,
        alpha=0.62,
        linewidths=0,
    )
    x = np.linspace(0, frame["eql_beta"].max() * 1.02, 100)
    beta = float(smr.iloc[0]["smr_beta"])
    axis.plot(x, beta * x, color=BURGUNDY, linewidth=1.4)
    lead = frame[frame["rsid"].eq("rs7404003")]
    if len(lead):
        axis.scatter(
            lead["eql_beta"],
            lead["gwas_beta"],
            marker="*",
            s=110,
            color=BURGUNDY,
            edgecolor=WHITE,
            linewidth=0.5,
            zorder=5,
        )
        axis.annotate(
            "rs7404003",
            (float(lead.iloc[0]["eql_beta"]), float(lead.iloc[0]["gwas_beta"])),
            xytext=(-48, 10),
            textcoords="offset points",
            fontsize=7,
        )
    axis.set_xlabel("Monocyte eQTL beta")
    axis.set_ylabel("IBD GWAS beta")
    axis.set_title(
        "Regional effect coupling\nand single-variant SMR",
        loc="left",
        fontweight="bold",
        fontsize=9,
    )
    cbar = axis.figure.colorbar(scatter, ax=axis, fraction=0.046, pad=0.03)
    cbar.set_label("r$^2$ with PRKCB lead", fontsize=7)
    cbar.ax.tick_params(labelsize=6)


def coloc_panel(axis: plt.Axes, sensitivity: pd.DataFrame) -> None:
    primary = sensitivity[
        sensitivity["L"].eq(5) & sensitivity["p12"].eq(5e-6)
    ].copy()
    gwas_sets = sorted(primary["idx1"].dropna().astype(int).unique())
    eqtl_sets = sorted(primary["idx2"].dropna().astype(int).unique())
    matrix = np.zeros((len(gwas_sets), len(eqtl_sets)))
    for row_index, gwas_set in enumerate(gwas_sets):
        for col_index, eqtl_set in enumerate(eqtl_sets):
            value = primary[
                primary["idx1"].eq(gwas_set)
                & primary["idx2"].eq(eqtl_set)
            ]["PP.H4.abf"]
            if len(value) and np.isfinite(value.iloc[0]):
                matrix[row_index, col_index] = float(value.iloc[0])
    image = axis.imshow(
        matrix,
        vmin=0,
        vmax=1,
        cmap=plt.matplotlib.colors.LinearSegmentedColormap.from_list(
            "coloc",
            [WHITE, LIME, PEACH, ORANGE, BURGUNDY],
        ),
        aspect="auto",
    )
    axis.set_xticks(np.arange(len(eqtl_sets)))
    axis.set_xticklabels([f"eQTL CS{value}" for value in eqtl_sets], fontsize=7)
    axis.set_yticks(np.arange(len(gwas_sets)))
    axis.set_yticklabels([f"GWAS CS{value}" for value in gwas_sets])
    axis.set_xlabel("eQTL credible set")
    axis.set_ylabel("GWAS credible set")
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            axis.text(
                col,
                row,
                "1.00" if value >= 0.995 else "<1e-6",
                ha="center",
                va="center",
                fontsize=7,
                color=WHITE if value >= 0.6 else INK,
            )
    axis.set_title(
        "SuSiE-colocalization at L=5\nshared credible-set pairs",
        loc="left",
        fontweight="bold",
        fontsize=9,
    )
    sensitivity_lines = []
    for L in [3, 5]:
        frame = sensitivity[sensitivity["L"].eq(L) & sensitivity["p12"].eq(5e-6)]
        sensitivity_lines.append(
            f"L={L}: {int((frame['PP.H4.abf'] >= 0.8).sum())}/{len(frame)}"
        )
    axis.text(
        0.5,
        -0.28,
        "  |  ".join(sensitivity_lines)
        + " H4>=0.8; p12 1e-6 to 1e-5 unchanged",
        transform=axis.transAxes,
        ha="center",
        va="top",
        fontsize=6.5,
        color=INK,
    )
    axis.set_xticks(np.arange(-0.5, len(eqtl_sets), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(gwas_sets), 1), minor=True)
    axis.grid(which="minor", color=WHITE, linewidth=2)
    axis.tick_params(which="minor", length=0)


def qtl_panel(axis: plt.Axes, coverage: pd.DataFrame) -> None:
    ordered = coverage.copy()
    matrix = np.column_stack(
        [
            ordered["exact_variant_coverage"].astype(int),
            ordered["association_supported"].astype(float) * 1.0,
        ]
    )
    cmap = plt.matplotlib.colors.ListedColormap(
        [WHITE, PEACH, GREEN]
    )
    if matrix.max() <= 1:
        shown = matrix
    else:
        shown = matrix
    axis.imshow(shown, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    axis.set_xticks([0, 1])
    axis.set_xticklabels(["Exact-variant\ncoverage", "Association\nsupported"], fontsize=7)
    axis.set_yticks(np.arange(len(ordered)))
    axis.set_yticklabels(
        [
            f"{row.molecular_layer} ({row.source.split('_')[0]})"
            for row in ordered.itertuples()
        ],
        fontsize=7,
    )
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = int(matrix[row, col])
            axis.text(
                col,
                row,
                "yes" if value else "no",
                ha="center",
                va="center",
                fontsize=7,
                color=INK,
            )
    axis.set_title(
        "Molecular QTL coverage boundary",
        loc="left",
        fontweight="bold",
        fontsize=9,
    )
    axis.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(ordered), 1), minor=True)
    axis.grid(which="minor", color=WHITE, linewidth=2)
    axis.tick_params(which="minor", length=0)


def signature_panel(axis: plt.Axes, reversal: pd.DataFrame) -> None:
    order = ["Oliver2024", "Crohn2023", "Garrido2023"]
    targets = ["PRKCB", "CHSY3", "ETS2", "JAK2", "STAT3"]
    width = 0.15
    x = np.arange(len(order))
    for target_index, target in enumerate(targets):
        values = []
        for cohort in order:
            row = reversal[
                reversal["cohort"].eq(cohort) & reversal["target"].eq(target)
            ]
            values.append(float(row.iloc[0]["correlation"]) if len(row) else np.nan)
        axis.bar(
            x + (target_index - 2) * width,
            values,
            width=width,
            color=TARGET_COLORS[target],
            label=target,
        )
    axis.axhline(0, color=INK, linewidth=1)
    axis.set_xticks(x)
    axis.set_xticklabels(["Oliver", "Crohn", "Garrido"])
    axis.set_ylabel("Disease-signature correlation")
    axis.set_title(
        "Public overexpression signature\nprojected onto disease macrophages",
        loc="left",
        fontweight="bold",
        fontsize=9,
    )
    axis.legend(frameon=False, ncol=5, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.16))


def virtual_panel(axis: plt.Axes, summary: pd.DataFrame) -> None:
    frame = summary[summary["target"].eq("PRKCB")].copy()
    pathways = ["TNF/NF-kB", "JAK/STAT", "IL-23", "Chemokine"]
    x = np.arange(len(pathways))
    for state_index, state in enumerate(["healthy", "active"]):
        values = []
        lows = []
        highs = []
        for pathway in pathways:
            row = frame[
                frame["state"].eq(state) & frame["pathway"].eq(pathway)
            ]
            if len(row):
                row = row.iloc[0]
                values.append(row["bootstrap_mean"])
                lows.append(row["bootstrap_mean"] - row["ci_low"])
                highs.append(row["ci_high"] - row["bootstrap_mean"])
            else:
                values.append(np.nan)
                lows.append(np.nan)
                highs.append(np.nan)
        axis.bar(
            x + (state_index - 0.5) * 0.36,
            values,
            width=0.34,
            color=STATE_COLORS[state],
            label=state.capitalize(),
        )
        axis.errorbar(
            x + (state_index - 0.5) * 0.36,
            values,
            yerr=np.asarray([lows, highs]),
            fmt="none",
            ecolor=INK,
            capsize=2,
            linewidth=0.8,
        )
    axis.axhline(0, color=INK, linewidth=1)
    axis.set_xticks(x)
    axis.set_xticklabels(["TNF/NF-kB", "JAK/STAT", "IL-23", "Chemokine"], fontsize=7)
    axis.set_ylabel("Predicted log-expression change")
    axis.set_title(
        "PRKCB virtual knockdown\nconditional network prediction",
        loc="left",
        fontweight="bold",
        fontsize=9,
    )
    axis.legend(frameon=False, fontsize=7)


def main() -> None:
    configure_style()
    table_dir = ROOT / "results" / "tables"
    figure_dir = ROOT / "results" / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    estimates = pd.read_csv(table_dir / "prkcb_mr_estimates.csv")
    variants = pd.read_csv(table_dir / "prkcb_ld_analysis_variants.csv")
    ld = pd.read_csv(table_dir / "prkcb_ld_with_lead.csv")
    smr = pd.read_csv(table_dir / "prkcb_smr_estimates.csv")
    sensitivity = pd.read_csv(
        table_dir / "prkcb_susie_coloc_sensitivity.csv"
    )
    coverage = pd.read_csv(
        table_dir / "prkcb_molecular_qtl_coverage.csv"
    )
    reversal = pd.read_csv(table_dir / "prkcb_signature_reversal.csv")
    virtual = pd.read_csv(
        table_dir / "prkcb_virtual_pathway_bootstrap_summary.csv"
    )

    fig = plt.figure(figsize=(15.0, 9.2))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0])
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
        fig.add_subplot(grid[1, 2]),
    ]
    mr_forest_panel(axes[0], estimates)
    regional_panel(axes[1], variants, ld, smr)
    coloc_panel(axes[2], sensitivity)
    qtl_panel(axes[3], coverage)
    signature_panel(axes[4], reversal)
    virtual_panel(axes[5], virtual)
    label_panels(axes)
    fig.tight_layout()
    base = figure_dir / "figure12_prkcb_mechanism"
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {base.with_suffix('.png')}")


if __name__ == "__main__":
    main()
