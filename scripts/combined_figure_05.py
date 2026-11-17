"""Combined Figure 5 -- the crossover, every network on one page.

One row per registered arm, three sowing dates across:

    row 1   hand-curated PSoup network      Oct 15 | Nov 15 | Dec 15
    row 2   Flash-P (30 nodes)              Oct 15 | Nov 15 | Dec 15
    row 3   Flash-P (69 nodes)              Oct 15 | Nov 15 | Dec 15

The contrast is the high-tillering representative (HT, max2) against the
low-tillering one (LT, nced3). The plotted quantity is

    (LT - HT) / HT * 100

so a POSITIVE value means the low-tillering genotype yields more. The blue
region is therefore the top half of every panel and the red the bottom.

Axes are shared across all panels, so a difference in where the points move
from the blue region into the red one is a real difference and not a rescaling
artefact. Palette, quadratic fit, shaded win regions and circled representative
years follow the paper's Figure 5, so the panels read as one figure.

The point of the figure: if the rows have the same shape, a network the agent
built from the literature carries the same agronomic conclusion as the
hand-curated one.

No figure title -- the figure is meant to be dropped into a document that
supplies its own caption. Font sizes are set for that: large enough to stay
legible after Word scales the image down to column width.

Writes 04_figures/figures_combined/Fig05_Crossover_AllNetworks.{png,svg}
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS as REGISTRY, FIG_ROOT  # noqa: E402

RUNS = ROOT / "03_runs" / "apsim"
OUT = FIG_ROOT / "figures_combined"

# (directory name, row label), one row per arm, registry order -- the reference
# arm is listed first in scripts/arms.json, so it is the top row. The row label
# carries the network's node and edge counts on a second line, because that size
# is the only thing distinguishing the two Flash-P rows from each other.
ARMS = [(a.name, a.fig_label_2line) for a in REGISTRY.values()
        if (RUNS / a.name / "DensityManagement_NoTiller.db").exists()]
SOWINGS = ["Oct15", "Nov15", "Dec15"]
SOWING_LABELS = {"Oct15": "Oct 15", "Nov15": "Nov 15", "Dec15": "Dec 15"}

HIGH, LOW = "max2", "nced3"
REP_DENSITY = 2.5
G_TO_T_PER_HA = 0.01

# Okabe-Ito, as in the paper's Figure 5
LOW_COLOR = "#009E73"
WIN_LOW_COLOR = "#2196F3"
WIN_HIGH_COLOR = "#E53935"
YEAR_COLORS = {"dry": "#B2182B", "crossover": "#542788", "wet": "#1F4E79"}
ANNOT_LABEL = {"dry": "Dry", "crossover": "Intermediate", "wet": "Wet"}
# Each callout goes in a different direction so the three never collide: the
# Dry point is the highest in its panel (label straight up), the Wet point the
# lowest and furthest right (label up-right, into the empty corner), and the
# Intermediate point sits on the zero line (label down-left, under the curve).
ANNOT_OFFSET = {"dry": (20, 56), "crossover": (-62, -58), "wet": (44, 40)}

# Sized for a figure that will be scaled down into a Word column.
FS_AXIS_LABEL = 26
FS_TICK = 22
FS_PANEL_TAG = 28
FS_ROW_LABEL = 30
FS_LEGEND = 23
FS_ANNOT = 17


def load(arm: str) -> pd.DataFrame:
    db = RUNS / arm / "DensityManagement_NoTiller.db"
    if not db.exists():
        raise SystemExit(f"missing {db} -- run the APSIM stage first")
    with sqlite3.connect(db) as c:
        # named columns rather than SELECT *: this now runs over one ~87 MB db per arm
        h = pd.read_sql(
            'SELECT "Clock.Today", Density, Genotype, SowingDate, '
            '"Sorghum.Grain.Wt" FROM HarvestReport', c)
    h["Density"] = pd.to_numeric(h["Density"], errors="coerce")
    h["Year"] = pd.to_datetime(h["Clock.Today"], errors="coerce").dt.year
    h["GrainYield_tha"] = h["Sorghum.Grain.Wt"] * G_TO_T_PER_HA
    return h[h["Density"] == REP_DENSITY]


def contrast(h: pd.DataFrame, sow: str) -> pd.DataFrame:
    """Per-year HT yield, and LT's percentage difference from it.

    Positive dLow means the low-tillering genotype yields more than the
    high-tillering one in that season.
    """
    piv = (h[h["SowingDate"] == sow]
           .pivot_table(index="Year", columns="Genotype",
                        values="GrainYield_tha", aggfunc="first")
           .dropna(subset=[HIGH, LOW]))
    piv["dLow"] = (piv[LOW] - piv[HIGH]) / piv[HIGH] * 100
    return piv


def rep_years(piv: pd.DataFrame) -> dict[str, int]:
    """The paper's per-sowing dry / crossover / wet ranking."""
    return {
        "dry":       int(piv["dLow"].idxmax()),
        "crossover": int(piv["dLow"].abs().idxmin()),
        "wet":       int((-piv["dLow"] * piv[HIGH]).idxmax()),
    }


def main() -> None:
    plt.rcParams.update({
        "font.size": FS_TICK,
        "axes.labelsize": FS_AXIS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
    })
    OUT.mkdir(parents=True, exist_ok=True)

    data = {arm: load(arm) for arm, _ in ARMS}

    # One y-range for all six panels, from the data itself. The shaded win
    # regions are drawn far beyond it and clipped, so they fill the panel
    # without dictating the scale. Extra headroom at the top: the "Dry" year is
    # the highest point in its panel and its callout is drawn above it.
    allv = np.concatenate([contrast(data[a], s)["dLow"].values
                           for a, _ in ARMS for s in SOWINGS])
    span = allv.max() - allv.min()
    ylim = (allv.min() - 0.16 * span, allv.max() + 0.34 * span)

    # Matching x headroom, so the Wet callout in the bottom-right of each panel
    # has somewhere to sit instead of running off the edge.
    allx = np.concatenate([contrast(data[a], s)[HIGH].values
                           for a, _ in ARMS for s in SOWINGS])
    xspan = allx.max() - allx.min()
    xlim = (allx.min() - 0.06 * xspan, allx.max() + 0.13 * xspan)

    nrow = len(ARMS)
    fig, axes = plt.subplots(nrow, 3, figsize=(24, 6.5 * nrow),
                             sharex=True, sharey=True, squeeze=False)

    for r, (arm, _) in enumerate(ARMS):
        for c, sow in enumerate(SOWINGS):
            ax = axes[r, c]
            piv = contrast(data[arm], sow)

            ax.axhline(0, color="black", linewidth=1.2, linestyle="--", alpha=0.5)
            ax.axhspan(0, 1e4,  alpha=0.05, color=WIN_LOW_COLOR)    # LT wins
            ax.axhspan(-1e4, 0, alpha=0.05, color=WIN_HIGH_COLOR)   # HT wins

            x = piv[HIGH].values
            y = piv["dLow"].values
            ax.scatter(x, y, color=LOW_COLOR, s=90, alpha=0.65, zorder=3)

            if len(x) > 3:
                coef = np.polyfit(x, y, 2)
                xs = np.linspace(x.min(), x.max(), 200)
                ax.plot(xs, np.polyval(coef, xs), color=LOW_COLOR,
                        linewidth=3.5, alpha=0.8)

            for key, yr in rep_years(piv).items():
                if yr not in piv.index:
                    continue
                row = piv.loc[yr]
                col = YEAR_COLORS[key]
                ax.scatter([row[HIGH]], [row["dLow"]], s=300, facecolors="none",
                           edgecolors=col, linewidths=3.0, zorder=5)

                dx, dy = ANNOT_OFFSET[key]
                # The Intermediate year lands wherever the curve crosses zero,
                # which is far left in some panels and mid-panel in others. A
                # fixed leftward offset runs off the axis in the former, so
                # flip it whenever the point is in the left third.
                if key == "crossover" and \
                        row[HIGH] < xlim[0] + (xlim[1] - xlim[0]) / 3:
                    dx = -dx

                ax.annotate(f"{yr}\n({ANNOT_LABEL[key]})",
                            xy=(row[HIGH], row["dLow"]),
                            xytext=(dx, dy), textcoords="offset points",
                            fontsize=FS_ANNOT, fontweight="bold", color=col,
                            ha="center",
                            arrowprops=dict(arrowstyle="->", color=col, lw=1.8),
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                      edgecolor=col, alpha=0.9))

            # Sowing date, top centre of every panel.
            ax.text(0.5, 0.965, SOWING_LABELS[sow], transform=ax.transAxes,
                    fontsize=FS_PANEL_TAG, fontweight="bold",
                    ha="center", va="top")

            ax.set_ylim(ylim)
            ax.set_xlim(xlim)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(axis="both", labelsize=FS_TICK)
            if r == nrow - 1:
                ax.set_xlabel("High-tillering (HT) yield (t/ha)")

        axes[r, 0].set_ylabel("Yield difference,\nLT − HT (%)")

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=LOW_COLOR,
                   markersize=13, alpha=0.7, label="One season"),
        plt.Line2D([0], [0], color=LOW_COLOR, linewidth=3.5, label="Quadratic fit"),
        Patch(facecolor=WIN_LOW_COLOR, alpha=0.18,
              label="Δ (LT − HT) > 0: Low tillering yields more"),
        Patch(facecolor=WIN_HIGH_COLOR, alpha=0.18,
              label="Δ (LT − HT) < 0: High tillering yields more"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=FS_LEGEND, bbox_to_anchor=(0.5, -0.012))

    plt.tight_layout(rect=[0.045, 0.045, 1, 1])

    # Row labels sit outside the axes entirely, so they cannot collide with the
    # y-axis labels. Placed after tight_layout, once the axes are in final
    # position, and read off those positions so each is centred on its row.
    for r, (_, row_label) in enumerate(ARMS):
        box = axes[r, 0].get_position()
        fig.text(0.006, box.y0 + box.height / 2, row_label, rotation=90,
                 fontsize=FS_ROW_LABEL, fontweight="bold",
                 ha="center", va="center", linespacing=1.35,
                 multialignment="center")

    png = OUT / "Fig05_Crossover_AllNetworks.png"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")


if __name__ == "__main__":
    main()
