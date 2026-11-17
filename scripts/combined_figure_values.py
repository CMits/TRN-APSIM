"""Combined value figures -- the two numbers APSIM actually receives.

The rank figure asks whether the two networks ORDER the genotypes the same way.
These two ask the harder question: are the VALUES the same? One figure per
APSIM parameter, both built from the same data and the same layout:

  FigPTN   mean potential tiller number, max(0, intercept + slope*SD) averaged
           over the realised supply/demand range -- what the crop model sees
  FigSlope tillerSdSlope -- how strongly a genotype converts supply/demand into
           tillers, i.e. its sucrose responsiveness carried through the network

  panel a  dumbbell -- one row per genotype, sorted by the hand-curated value,
           with both networks' values marked and joined. A short bar is a
           genotype the two networks agree on. Each row carries its rank in
           both networks (hand -> flashp) so ordering is readable off the same
           panel as magnitude.

  panel b  value against value, with the 1:1 line. Spearman rho (ordering) and
           Pearson r (magnitude) are printed beside it -- the pair matters,
           because rank can hold while the numbers drift.

Ranks for the tiller-number figure are the ones already in Arm_Comparison.csv,
so this figure and the rank figure cannot disagree. Slope ranks are computed
here, ties broken by tiller number.

The drawing lives in draw_pair(), which takes prepared axes instead of making
its own figure, so combined_figure_pairs.py can stack two comparisons into one
four-panel figure without a second copy of the code.

No figure title -- the figures are dropped into a document that supplies its own
caption. They are built at 6.5 x 5.77 in, the A4 text block at 2 cm margins, so
Word inserts them at 100% and the point sizes below are the printed ones.

Panel a's x-axis names the parameter without repeating that it is a mean over
the realised supply/demand range -- the axis has 2 in of width and the caption
is the place for the definition.

Reads  03_runs/Arm_Comparison_long.csv
Writes 04_figures/figures_combined/FigPTN_TillerNumber_{pair}.{png,svg}
       04_figures/figures_combined/FigSlope_SupplyDemandSlope_{pair}.{png,svg}
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS as REGISTRY, REFERENCE, FIG_ROOT  # noqa: E402

SRC = ROOT / "03_runs" / "Arm_Comparison_long.csv"
OUT = FIG_ROOT / "figures_combined"

# Only the high/low contrast is highlighted. `lbo` is also carried into APSIM as the
# Medium representative, but the results show only high vs low, so it is drawn as an
# ordinary genotype like the other 29.
REPS = {"max2": ("#E53935", "High tillering"),
        "nced3": ("#4CAF50", "Low tillering")}

EXTREME_N = 5

# ----------------------------------------------------------------- print size
# Built at the width it OCCUPIES on the page, so Word places it at 100% and the
# point sizes below are the point sizes the reader gets. These used to be drawn
# on a 19.5 in canvas carrying the same type, which Word shrank to a third: 13 pt
# genotype names arriving as 4.3 pt. 6.5 in is the A4 text block at 2 cm margins.
FIG_W = 6.5
DPI = 600

# ------------------------------------------------------------------- geometry
# An explicit inch budget, not tight_layout: panel b is aspect="equal", which
# tight_layout warns it cannot lay out, and placing the axes directly makes the
# saved file exactly FIG_W wide, so Word inserts it at 100%.
#
# NAME_W is set by TYPE -- the y tick labels are "brc12hb214053   28→27" at
# FS_NAME, 1.09 in before the dumbbells get any width at all. DUMB_W then has
# to hold the slope x-axis label, which is why that one is set on two lines.
PAD_L, PAD_R = 0.03, 0.06
PAD_T, PAD_B = 0.22, 0.08
NAME_W = 1.16                 # genotype + rank arrow, left of panel a
DUMB_W = 2.02                 # panel a's plotting area
GUTTER = 0.50                 # panel b's y-axis label and tick labels
XAXIS_H = 0.62                # panel a's tick labels and x-axis label
XLABEL_H = 0.50               # panel b's, above the legend block
DUMB_H = 4.85                 # 33 genotype rows at a pitch that clears FS_NAME
ROW_H = DUMB_H + XAXIS_H
FIG_H = PAD_T + ROW_H + PAD_B

# Ink -- strokes and marker areas -- stays proportional to the old canvas; only
# TYPE is set in true printed points, which is the whole change.
DESIGN_W = 19.5
INK = FIG_W / DESIGN_W


def lw(pt: float) -> float:
    """A design stroke width, rescaled to the printed canvas."""
    return pt * INK


def area(pt2: float) -> float:
    """A design marker size. scatter's `s` is an AREA, so it scales as INK²."""
    return pt2 * INK * INK


# True printed points, at the size they appear in the document.
FS_AXIS = 9.0
FS_TICK = 8.0
FS_TAG = 12.0
FS_NAME = 6.4
FS_LEGEND = 7.8
FS_STAT = 7.8

# The two parameters APSIM receives. Shared with combined_figure_pairs.py so the
# stacked four-panel figures cannot end up describing a different quantity from
# the single-pair ones. `stat_unit` is the short form: the statistics block sits
# in a 2.7 in column, which the full slope unit string does not fit into.
METRICS = (
    dict(hand="hand_TTN_mean", flashp="flashp_TTN_mean",
         hrank="hand_rank", frank="flashp_rank",
         axis_label="potential tiller number",
         unit="tillers", stat_unit="tillers",
         tag="PTN", stem="FigPTN_TillerNumber"),
    dict(hand="hand_slope", flashp="flashp_slope",
         hrank="hand_slope_rank", frank="flashp_slope_rank",
         axis_label="tillerSdSlope\n(tillers per unit supply/demand)",
         unit="tillerSdSlope", stat_unit="slope units",
         tag="SLOPE", stem="FigSlope_SupplyDemandSlope"),
)


def row_axes(fig, top_in: float = PAD_T):
    """The three axes of one comparison, placed `top_in` inches below the top.

    (dumbbell, value-value scatter, legend + statistics). Shared with
    combined_figure_pairs.py, which stacks one of these bands per non-reference
    arm, so a stacked row and its standalone twin have identical geometry.
    """
    w, h = fig.get_size_inches()

    def rect(x, y, dx, dy):          # y measured DOWN from the figure top
        return [x / w, 1.0 - (y + dy) / h, dx / w, dy / h]

    x_r = PAD_L + NAME_W + DUMB_W + GUTTER
    w_r = w - x_r - PAD_R            # panel b is square: aspect="equal"
    return (fig.add_axes(rect(PAD_L + NAME_W, top_in, DUMB_W, DUMB_H)),
            fig.add_axes(rect(x_r, top_in, w_r, w_r)),
            fig.add_axes(rect(x_r, top_in + w_r + XLABEL_H,
                              w_r, ROW_H - w_r - XLABEL_H)))


def set_style() -> None:
    plt.rcParams.update({"font.size": FS_TICK,
                         "axes.labelsize": FS_AXIS,
                         "xtick.labelsize": FS_TICK,
                         "ytick.labelsize": FS_TICK})


def rankdata(a: np.ndarray) -> np.ndarray:
    """Average ranks, ties shared -- both parameters have genuine ties at zero."""
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    r[order] = np.arange(1, len(a) + 1)
    for v in np.unique(a):
        m = a == v
        r[m] = r[m].mean()
    return r


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def draw_pair(ax_d, ax_v, ax_t, df: pd.DataFrame, metric: dict,
              labels: tuple[str, str], fig_labels: tuple[str, str],
              colors: tuple[str, str],
              tags: tuple[str, str] = ("a", "b")) -> dict:
    """One reference-vs-arm comparison, drawn into three prepared axes.

    labels is the bare name, used where space is tight (the panel b axes, the
    rank column head). fig_labels carries the node/edge counts as well and goes
    in the legend, which is where the reader identifies the network -- both
    FlashP arms are called "Flash-P", so the size is what tells them apart.
    tags names the two panels, so a stacked figure can label its lower half c/d.
    """
    hand, flashp = metric["hand"], metric["flashp"]
    hrank, frank = metric["hrank"], metric["frank"]
    unit = metric["unit"]
    ax_t.axis("off")

    # Sort by rank, not by value. Both parameters have tied values (five zero
    # slopes in the hand-curated arm), and sorting on the value leaves those
    # rows in arbitrary order, which makes the rank column look scrambled.
    # Rank is derived from the value, so this is still value-monotone.
    d = df.sort_values(hrank, ascending=False).reset_index(drop=True)
    n = len(d)
    diff = d[flashp] - d[hand]

    rho = spearman(df[hand].values, df[flashp].values)
    r_p = float(np.corrcoef(df[hand], df[flashp])[0, 1])
    top_shared = len(set(df.nlargest(EXTREME_N, hand)["apsim_name"]) &
                     set(df.nlargest(EXTREME_N, flashp)["apsim_name"]))
    bot_shared = len(set(df.nsmallest(EXTREME_N, hand)["apsim_name"]) &
                     set(df.nsmallest(EXTREME_N, flashp)["apsim_name"]))

    # ---------------------------------------------------------------- panel a
    ax_d.axvline(0, color="#BBBBBB", linewidth=lw(1.2), zorder=0)
    for i, row in d.iterrows():
        rep = REPS.get(row["apsim_name"])
        ax_d.plot([row[hand], row[flashp]], [i, i], color="#9E9E9E",
                  linewidth=lw(3.4 if rep else 2.2), alpha=0.75, zorder=2,
                  solid_capstyle="round")
        ax_d.scatter([row[hand]], [i], s=area(170 if rep else 105),
                     color=colors[0], zorder=4,
                     edgecolors="black" if rep else "none", linewidths=lw(1.3))
        ax_d.scatter([row[flashp]], [i], s=area(170 if rep else 105),
                     color=colors[1], zorder=4,
                     edgecolors="black" if rep else "none", linewidths=lw(1.3))

    ax_d.set_yticks(range(n))
    ax_d.set_yticklabels([f"{r['apsim_name']}   {int(r[hrank])}→{int(r[frank])}"
                          for _, r in d.iterrows()], fontsize=FS_NAME)
    for i, row in d.iterrows():
        rep = REPS.get(row["apsim_name"])
        if rep:
            lab = ax_d.get_yticklabels()[i]
            lab.set_color(rep[0])
            lab.set_fontweight("bold")

    ax_d.set_ylim(-1.0, n)
    ax_d.invert_yaxis()
    ax_d.set_xlabel(metric["axis_label"])
    ax_d.grid(axis="x", alpha=0.25, linewidth=lw(1.0))
    ax_d.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax_d.spines[s].set_visible(False)
    ax_d.tick_params(axis="y", length=0)
    ax_d.text(-NAME_W / DUMB_W, 1.005, tags[0], transform=ax_d.transAxes,
              fontsize=FS_TAG,
              fontweight="bold", ha="left", va="bottom")
    ax_d.text(1.0, 1.005,
              f"genotype   rank {labels[0]}→{labels[1]}",
              transform=ax_d.transAxes,
              fontsize=FS_NAME + 0.6, style="italic", color="#555555",
              ha="right", va="bottom")

    # ---------------------------------------------------------------- panel b
    lo = min(df[hand].min(), df[flashp].min())
    hi = max(df[hand].max(), df[flashp].max())
    pad = 0.12 * (hi - lo)
    lim = (lo - pad, hi + pad)
    ax_v.plot(lim, lim, linestyle="--", color="black", linewidth=lw(1.6),
              alpha=0.55, zorder=1)

    for _, row in df.iterrows():
        rep = REPS.get(row["apsim_name"])
        ax_v.scatter(row[hand], row[flashp], s=area(250 if rep else 110),
                     color=rep[0] if rep else "#9E9E9E",
                     alpha=1.0 if rep else 0.75,
                     edgecolors="black" if rep else "none",
                     linewidths=lw(1.4), zorder=5 if rep else 3)

    # Points pile up near the 1:1 line, so every callout is pushed well off it,
    # into the empty half below the diagonal, with the two low cultivars kept
    # far apart vertically so their boxes cannot touch. Offsets are fractions of
    # the axis span, not points: aspect="equal" rescales the axes to fit the
    # data, so a point offset that clears the axis on one figure runs off it on
    # the other.
    rep_offset = {"max2": (-0.32, 0.0), "nced3": (0.36, 0.02)}
    span = lim[1] - lim[0]
    for name, (color, role) in REPS.items():
        row = df[df["apsim_name"] == name].iloc[0]
        fx, fy = rep_offset[name]
        ax_v.annotate(f"{name}\n({role})", xy=(row[hand], row[flashp]),
                      xytext=(row[hand] + fx * span, row[flashp] + fy * span),
                      textcoords="data",
                      fontsize=FS_NAME + 1, fontweight="bold", color=color,
                      ha="center", va="center",
                      arrowprops=dict(arrowstyle="->", color=color, lw=lw(1.8)),
                      bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                                edgecolor=color, alpha=0.92))

    ax_v.set_xlim(lim)
    ax_v.set_ylim(lim)
    ax_v.set_aspect("equal", adjustable="box")
    ax_v.set_xlabel(f"{labels[0]}  ({unit})")
    ax_v.set_ylabel(f"{labels[1]}  ({unit})")
    ax_v.spines["top"].set_visible(False)
    ax_v.spines["right"].set_visible(False)
    ax_v.text(-0.175, 1.02, tags[1], transform=ax_v.transAxes, fontsize=FS_TAG,
              fontweight="bold", ha="left", va="bottom")

    # ------------------------------------------------------- legend + numbers
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=colors[0],
               markersize=15 * INK, label=fig_labels[0]),
        Line2D([0], [0], marker="o", linestyle="", color=colors[1],
               markersize=15 * INK, label=fig_labels[1]),
        Line2D([0], [0], color="#9E9E9E", lw=lw(3.4),
               label="difference between networks"),
        Line2D([0], [0], color="black", lw=lw(1.6), linestyle="--", alpha=0.6,
               label="1:1 (identical value)"),
    ]
    ax_t.legend(handles=handles, loc="upper left", frameon=False,
                fontsize=FS_LEGEND, bbox_to_anchor=(0.02, 1.02),
                handlelength=1.6, labelspacing=0.8)

    stats = "\n".join([
        f"Spearman ρ = {rho:.3f}   (ordering, n = {n})",
        f"Pearson r = {r_p:.3f}   (magnitude)",
        f"median |difference| = {diff.abs().median():.3f} {metric['stat_unit']}",
        f"highest {EXTREME_N} shared: {top_shared}/{EXTREME_N}",
        f"lowest {EXTREME_N} shared: {bot_shared}/{EXTREME_N}",
    ])
    ax_t.text(0.02, 0.60, stats, transform=ax_t.transAxes, ha="left", va="top",
              fontsize=FS_STAT, linespacing=1.85,
              bbox=dict(boxstyle="round,pad=0.6", facecolor="#F5F5F5",
                        edgecolor="#BBBBBB"))

    return dict(rho=rho, r_p=r_p, top=top_shared, bot=bot_shared)


def save(fig, png: Path) -> None:
    # No bbox_inches="tight": the axes are placed in inches, so the file is
    # already exactly FIG_W wide and cropping it would only make Word rescale.
    fig.savefig(png, dpi=DPI)
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True)
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")


def report(metric: dict, st: dict) -> None:
    print(f"  {metric['tag']}: rho {st['rho']:.3f}  r {st['r_p']:.3f}  "
          f"top{EXTREME_N} {st['top']}/{EXTREME_N}  "
          f"bottom{EXTREME_N} {st['bot']}/{EXTREME_N}")


def make_figure(df: pd.DataFrame, metric: dict, out_name: str,
                labels: tuple[str, str], fig_labels: tuple[str, str],
                colors: tuple[str, str]) -> None:
    set_style()
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    st = draw_pair(*row_axes(fig), df, metric, labels, fig_labels, colors)
    save(fig, OUT / out_name)
    report(metric, st)


def derive(df: pd.DataFrame) -> pd.DataFrame:
    """Add the two columns these figures rank on, for both arms."""
    # Mean clamped TTN over the realised SD range -- exactly the quantity the
    # ranks in the CSV were built on, so the two figures stay consistent.
    for arm in ("hand", "flashp"):
        df[f"{arm}_TTN_mean"] = (df[f"{arm}_TTN_lo"] + df[f"{arm}_TTN_hi"]) / 2
        # Slope ranks are new. Ties (five zero slopes in the hand-curated arm)
        # are broken by tiller number, so equal-slope genotypes are still
        # ordered by what they do in the crop.
        df[f"{arm}_slope_rank"] = (
            df.sort_values([f"{arm}_slope", f"{arm}_TTN_mean"])
              .assign(_r=range(1, len(df) + 1))
              .sort_index()["_r"])
    return df


def pair_frame(long: pd.DataFrame, other: str) -> pd.DataFrame:
    """Reference arm vs one other arm, in the wide shape these figures expect."""
    a = long[long["arm"] == REFERENCE].set_index("apsim_name")
    b = long[long["arm"] == other].set_index("apsim_name").loc[a.index]
    df = pd.DataFrame({"apsim_name": a.index, "genotype": a["genotype"].values})
    for src, pre in ((a, "hand"), (b, "flashp")):
        for col in ("intercept", "slope", "TTN_lo", "TTN_hi", "live_frac"):
            df[f"{pre}_{col}"] = src[col].values
        df[f"{pre}_rank"] = src["rank"].values
    df["rank_shift"] = df["flashp_rank"] - df["hand_rank"]
    return derive(df.sort_values("hand_rank").reset_index(drop=True))


def pair_style(other: str) -> tuple[tuple[str, str], tuple[str, str],
                                    tuple[str, str]]:
    """(labels, fig_labels, colors) for the reference arm versus `other`."""
    return ((REGISTRY[REFERENCE].label, REGISTRY[other].label),
            (REGISTRY[REFERENCE].fig_label, REGISTRY[other].fig_label),
            (REGISTRY[REFERENCE].color, REGISTRY[other].color))


def load_long() -> pd.DataFrame:
    if not SRC.exists():
        raise SystemExit(f"missing {SRC} -- run build_apsim_parameters.py first")
    return pd.read_csv(SRC)


def others_present(long: pd.DataFrame) -> list[str]:
    """Non-reference arms that actually have rows, in registry order."""
    present = [a for a in REGISTRY if a in set(long["arm"])]
    rest = [a for a in present if a != REFERENCE]
    if REFERENCE not in present or not rest:
        raise SystemExit(
            f"need {REFERENCE} and at least one other arm in {SRC.name}")
    return rest


def main() -> None:
    long = load_long()
    OUT.mkdir(parents=True, exist_ok=True)
    for other in others_present(long):
        labels, fig_labels, colors = pair_style(other)
        # the original pair keeps its filenames so existing references still resolve
        suffix = "BothNetworks" if other == "flashp_suc" else other
        print(f"\n=== {REFERENCE} vs {other} ===")
        df = pair_frame(long, other)
        for m in METRICS:
            make_figure(df, m, f"{m['stem']}_{suffix}.png",
                        labels, fig_labels, colors)


if __name__ == "__main__":
    main()
