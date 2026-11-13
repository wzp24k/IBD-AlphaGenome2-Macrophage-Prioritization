from __future__ import annotations

from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib as mpl
import seaborn as sns


GREEN = "#66BC98"
LIGHT_GREEN = "#AAD09D"
LIME = "#E3EA96"
PEACH = "#FDC089"
ORANGE = "#E26844"
BURGUNDY = "#8A233F"

WHITE = "#FFFFFF"
PAPER = "#FCFCF8"
INK = "#343434"
AXIS = "#5C5C5C"
GRID = "#E7E7E1"
NEUTRAL = "#B9B9B2"

PALETTE = [
    BURGUNDY,
    ORANGE,
    GREEN,
    PEACH,
    LIGHT_GREEN,
    LIME,
]

TIER_COLORS = {
    "Tier 1": BURGUNDY,
    "Tier 2": ORANGE,
    "Tier 3": GREEN,
    "Tier 4": NEUTRAL,
}

STATE_COLORS = {
    "healthy": GREEN,
    "active": ORANGE,
    "Ctrl": GREEN,
    "LPS": ORANGE,
    True: BURGUNDY,
    False: NEUTRAL,
}

TARGET_COLORS = {
    "PDGFB": BURGUNDY,
    "PRKCB": GREEN,
    "CHSY3": ORANGE,
    "PRKCQ": ORANGE,
    "IL23R": PEACH,
    "TNFSF15": GREEN,
    "MAP3K8": LIGHT_GREEN,
    "TNF": LIME,
    "TNFAIP3": AXIS,
}

PATHWAY_COLORS = {
    "Both": BURGUNDY,
    "IL-2/STAT5": ORANGE,
    "IL-6/JAK/STAT3": GREEN,
    "IL-23": PEACH,
    "Other": NEUTRAL,
}

DIVERGING = LinearSegmentedColormap.from_list(
    "ibd_diverging",
    [GREEN, LIGHT_GREEN, PAPER, PEACH, ORANGE, BURGUNDY],
)

SEQUENTIAL = LinearSegmentedColormap.from_list(
    "ibd_sequential",
    [PAPER, LIME, PEACH, ORANGE, BURGUNDY],
)

COVERAGE = ListedColormap([PAPER, GREEN])

CONCORDANCE = LinearSegmentedColormap.from_list(
    "ibd_concordance",
    [PAPER, PEACH, ORANGE, BURGUNDY],
)


def configure_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "figure.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "font.family": "DejaVu Sans",
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": AXIS,
            "axes.titlecolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "axes.grid.axis": "both",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def label_panels(axes, labels=None, x: float = -0.08, y: float = 1.04) -> None:
    axes = list(axes)
    if labels is None:
        labels = [chr(ord("A") + index) for index in range(len(axes))]
    for axis, label in zip(axes, labels):
        axis.text(
            x,
            y,
            label,
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color=INK,
            clip_on=False,
        )
