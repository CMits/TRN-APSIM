"""Combined rank figure -- does the FlashP network order the genotypes the same way?

The whole coupling collapses to two APSIM numbers per genotype, so the only
thing that has to survive a change of network is the ORDER: which mutants
tiller most, which tiller least. This figure shows that order directly.

  panel a   slopegraph -- every genotype's rank in the hand-curated network
            joined to its rank in the FlashP-built one. A horizontal line is a
            genotype the two networks agree on; a steep line is one they do not.

  panel b   rank against rank, with the 1:1 line. Spearman rho and the overlap
            of the extreme groups are printed beside it.

Ranking is by clamped mean potential tiller number over the REALISED supply/
demand range (11.70..42.36, read off the published APSIM run), not by intercept:
the intercept is a threshold in TTN = max(0, intercept + slope*SD), not a tiller
count, so ranking on it alone is wrong.

The drawing lives in draw_pair(), which takes prepared axes instead of making
its own figure, so combined_figure_pairs.py can stack two comparisons into one
four-panel figure without a second copy of the code.

No figure title -- the figure is dropped into a document that supplies its own
caption. It is built at 6.5 x 5.85 in, the A4 text block at 2 cm margins, so
Word inserts it at 100% and the point sizes below are the printed ones.

Reads  03_runs/Arm_Comparison_long.csv
Writes 04_figures/figures_combined/FigRank_GenotypeRanks_{pair}.{png,svg}
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS as REGISTRY, REFERENCE, FIG_ROOT  # noqa: E402

SRC = ROOT / "03_runs" / "Arm_Comparison_long.csv"
OUT = FIG_ROOT / "figures_combined"

# The two cultivars the agronomic comparison contrasts, coloured by the role they
# play there. `lbo` is also carried into APSIM as the Medium representative, but it
# is deliberately NOT highlighted here: the results show only the high/low contrast,
# so lbo is drawn as an ordinary genotype like the other 29.
REPS = {"max2": ("#E53935", "High tillering"),
        "nced3": ("#4CAF50", "Low tillering")}

# Movement classes, kept so a future network that does reorder the genotypes is
# flagged rather than silently averaged away. On the networks built so far the
# largest shift is 3 ranks, so the >= 6 class is empty and the legend omits it.
MOVE_COLORS = {"stable": "#9E9E9E", "small": "#E8A33D", "large": "#7B3294"}
EXTREME_N = 5
BAND_HIGH, BAND_LOW = "#E53935", "#2196F3"

# ----------------------------------------------------------------- print size
# The figure is built at the width it OCCUPIES on the page, so Word places it at
# 100% and the point sizes below are the point sizes the reader gets. It used to
# be drawn on a 19.5 in canvas carrying the same type, which Word then shrank to
# a third: 13 pt genotype names arriving as 4.3 pt.
#
# 6.5 in is the A4 text block at 2 cm margins.
FIG_W = 6.5
DPI = 600

# ------------------------------------------------------------------- geometry
# An explicit inch budget, not tight_layout. Panel b is aspect="equal", which
# tight_layout warns it cannot lay out, and the slack it did distribute went to
# the wrong column: panel a's genotype names ran into their own rank numbers.
# Placing the axes directly also makes the saved file exactly FIG_W wide, so
# Word inserts it at 100% and the point sizes above are the printed ones.
#
# SLOPE_W is set by TYPE, not by the lines: panel a is a name column, a rank
# column, one unit of slopegraph, another rank and name column, then the band
# labels, and at FS_NAME the longest name ("brc12hb214053") is 0.62 in on its own.
PAD_L, PAD_R = 0.03, 0.06     # panel a keeps its names inside its own x-limits
PAD_T, PAD_B = 0.22, 0.08
SLOPE_W = 3.54
GUTTER = 0.46                 # panel b's y-axis label and tick labels
XLABEL_H = 0.55               # panel b's x-axis label, above the legend block
ROW_H = 5.55                  # 37 rank units at a pitch that clears FS_NAME
FIG_H = PAD_T + ROW_H + PAD_B

# Ink -- stroke widths, marker areas, callout offsets -- was proportioned to the
# old canvas and stays proportional to it; only TYPE is set in true printed
# points, which is the whole change. Everything measured in points goes through
# lw() or area() so the drawing keeps the weights it was tuned with.
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
FS_HEAD = 9.5
FS_LEGEND = 7.8
FS_STAT = 7.8


def row_axes(fig, top_in: float = PAD_T):
    """The three axes of one comparison, placed `top_in` inches below the top.

    (slopegraph, rank-rank scatter, legend + statistics). Shared with
    combined_figure_pairs.py, which stacks one of these bands per non-reference
    arm, so a stacked row and its standalone twin have identical geometry.
    """
    w, h = fig.get_size_inches()

    def rect(x, y, dx, dy):          # y measured DOWN from the figure top
        return [x / w, 1.0 - (y + dy) / h, dx / w, dy / h]

    x_r = PAD_L + SLOPE_W + GUTTER
    w_r = w - x_r - PAD_R            # panel b is square: aspect="equal"
    return (fig.add_axes(rect(PAD_L, top_in, SLOPE_W, ROW_H)),
            fig.add_axes(rect(x_r, top_in, w_r, w_r)),
            fig.add_axes(rect(x_r, top_in + w_r + XLABEL_H,
                              w_r, ROW_H - w_r - XLABEL_H)))


def set_style() -> None:
    plt.rcParams.update({"font.size": FS_TICK,
                         "axes.labelsize": FS_AXIS,
                         "xtick.labelsize": FS_TICK,
                         "ytick.labelsize": FS_TICK})


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Ranks are already 1..n with no ties, so the shortcut formula is exact."""
    d = a - b
    n = len(a)
    return 1.0 - 6.0 * float((d ** 2).sum()) / (n * (n * n - 1))


def column_head(a) -> tuple[str, str]:
    """(name, size) for one slopegraph column.

    The two headers sit a single x-unit apart, so a long name set on one line
    runs straight into its neighbour: the name is wrapped, and the node/edge
    count is returned separately so it can be set smaller still.
    """
    lines = textwrap.wrap(a.label, 16, break_on_hyphens=False) or [a.label]
    return "\n".join(lines), a.size_label


def move_class(shift: int) -> str:
    s = abs(shift)
    return "stable" if s <= 2 else ("small" if s <= 5 else "large")


def pair_frame(long: pd.DataFrame, other: str) -> pd.DataFrame:
    """Reference arm vs one other arm, in the wide shape this figure was built for.

    The slopegraph and the rank-rank scatter are both inherently pairwise, so a
    third arm gets its own panel pair rather than a third column.
    """
    a = long[long["arm"] == REFERENCE].set_index("apsim_name")
    b = long[long["arm"] == other].set_index("apsim_name")
    df = pd.DataFrame({
        "apsim_name": a.index,
        "genotype": a["genotype"].values,
        "hand_rank": a["rank"].values,
        "flashp_rank": b.loc[a.index, "rank"].values,
    })
    df["rank_shift"] = df["flashp_rank"] - df["hand_rank"]
    return df.sort_values("hand_rank").reset_index(drop=True)


def pair_style(other: str) -> tuple[tuple[str, str],
                                    tuple[tuple[str, str], tuple[str, str]]]:
    """(labels, fig_labels) for the reference arm versus `other`."""
    return ((REGISTRY[REFERENCE].label, REGISTRY[other].label),
            (column_head(REGISTRY[REFERENCE]), column_head(REGISTRY[other])))


def draw_pair(ax_s, ax_r, ax_t, df: pd.DataFrame, labels: tuple[str, str],
              fig_labels: tuple[tuple[str, str], tuple[str, str]],
              tags: tuple[str, str] = ("a", "b")) -> dict:
    """One reference-vs-arm comparison, drawn into three prepared axes.

    labels is the bare name, used for the panel b axes. fig_labels is one
    (wrapped name, size) pair per slopegraph column, from column_head().
    tags names the two panels, so a stacked figure can label its lower half c/d.
    """
    ax_t.axis("off")
    df = df.copy()
    df["move"] = df["rank_shift"].apply(move_class)
    n = len(df)

    rho = spearman(df["hand_rank"].values, df["flashp_rank"].values)
    top_hand = set(df.nlargest(EXTREME_N, "hand_rank")["apsim_name"])
    top_fp = set(df.nlargest(EXTREME_N, "flashp_rank")["apsim_name"])
    bot_hand = set(df.nsmallest(EXTREME_N, "hand_rank")["apsim_name"])
    bot_fp = set(df.nsmallest(EXTREME_N, "flashp_rank")["apsim_name"])
    counts = df["move"].value_counts()

    # ---------------------------------------------------------------- panel a
    # Extreme bands first, so every line is drawn over them.
    ax_s.axhspan(n - EXTREME_N + 0.5, n + 0.5, color=BAND_HIGH, alpha=0.06, zorder=0)
    ax_s.axhspan(0.5, EXTREME_N + 0.5, color=BAND_LOW, alpha=0.06, zorder=0)

    for _, r in df.iterrows():
        name = r["apsim_name"]
        rep = REPS.get(name)
        color = rep[0] if rep else MOVE_COLORS[r["move"]]
        width = 4.2 if rep else (2.9 if r["move"] == "large" else 2.0)
        alpha = 1.0 if rep else (0.85 if r["move"] == "large" else 0.55)
        ax_s.plot([0, 1], [r["hand_rank"], r["flashp_rank"]],
                  color=color, linewidth=lw(width), alpha=alpha,
                  solid_capstyle="round", zorder=4 if rep else 2)
        ax_s.scatter([0, 1], [r["hand_rank"], r["flashp_rank"]],
                     s=area(95 if rep else 45), color=color,
                     alpha=alpha, zorder=5 if rep else 3)

        weight = "bold" if rep else "normal"
        # Rank number in its own narrow column next to the line, name outside it.
        for x_name, x_num, ha, rank in ((-0.206, -0.049, "right", r["hand_rank"]),
                                        (1.206, 1.049, "left", r["flashp_rank"])):
            ax_s.text(x_name, rank, name, ha=ha, va="center",
                      fontsize=FS_NAME, fontweight=weight, color=color)
            ax_s.text(x_num, rank, str(rank), ha=ha, va="center",
                      fontsize=FS_NAME - 0.6, fontweight=weight,
                      color=color, alpha=0.75)

    for x, (head, size) in zip((0, 1), fig_labels):
        ax_s.text(x, n + 2.72, head, ha="center", va="bottom",
                  fontsize=FS_HEAD, fontweight="bold", linespacing=1.25,
                  multialignment="center")
        if size:
            ax_s.text(x, n + 1.95, size, ha="center", va="bottom",
                      fontsize=FS_NAME, color="#555555")

    # Band labels go beyond the longest genotype name, not on top of it.
    ax_s.text(2.00, n - (EXTREME_N - 1) / 2, f"highest-\ntillering {EXTREME_N}",
              ha="left", va="center", fontsize=FS_NAME + 1,
              fontweight="bold", color=BAND_HIGH)
    ax_s.text(2.00, (EXTREME_N + 1) / 2, f"lowest-\ntillering {EXTREME_N}",
              ha="left", va="center", fontsize=FS_NAME + 1,
              fontweight="bold", color=BAND_LOW)

    ax_s.set_xlim(-0.90, 2.57)
    ax_s.set_ylim(-1.6, n + 3.4)
    ax_s.axis("off")
    ax_s.text(-0.90, n + 2.72, tags[0], fontsize=FS_TAG, fontweight="bold",
              ha="left", va="bottom")
    ax_s.text(0.5, -1.2, "rank by potential tiller number  (1 = lowest)",
              ha="center", va="center",
              fontsize=FS_TICK, style="italic", color="#555555")

    # ---------------------------------------------------------------- panel b
    lim = (0.0, n + 1.0)
    ax_r.add_patch(Rectangle((n - EXTREME_N + 0.5, n - EXTREME_N + 0.5),
                             EXTREME_N + 0.5, EXTREME_N + 0.5,
                             facecolor=BAND_HIGH, alpha=0.07, zorder=0))
    ax_r.add_patch(Rectangle((0.5, 0.5), EXTREME_N, EXTREME_N,
                             facecolor=BAND_LOW, alpha=0.07, zorder=0))
    ax_r.plot(lim, lim, linestyle="--", color="black", linewidth=lw(1.6),
              alpha=0.55, zorder=1)

    for _, r in df.iterrows():
        rep = REPS.get(r["apsim_name"])
        color = rep[0] if rep else MOVE_COLORS[r["move"]]
        ax_r.scatter(r["hand_rank"], r["flashp_rank"],
                     s=area(250 if rep else 110), color=color,
                     alpha=1.0 if rep else 0.8,
                     edgecolors="black" if rep else "none",
                     linewidths=lw(1.4), zorder=5 if rep else 3)

    # Only the three APSIM cultivars are labelled -- naming all 32 twice would
    # repeat panel a and bury the 1:1 relationship. The points hug the diagonal,
    # so every callout is pushed off it, on a leader. Offsets are fractions of
    # the axis span rather than points, so they survive a change of canvas size
    # -- the same reason combined_figure_values.py places its callouts that way.
    rep_offset = {"max2": (-0.40, -0.075), "nced3": (0.40, 0.055)}
    for name, (color, role) in REPS.items():
        r = df[df["apsim_name"] == name].iloc[0]
        fx, fy = rep_offset[name]
        span = lim[1] - lim[0]
        ax_r.annotate(f"{name}\n({role})",
                      xy=(r["hand_rank"], r["flashp_rank"]),
                      xytext=(r["hand_rank"] + fx * span,
                              r["flashp_rank"] + fy * span),
                      textcoords="data",
                      fontsize=FS_NAME + 1, fontweight="bold", color=color,
                      ha="center", va="center",
                      arrowprops=dict(arrowstyle="->", color=color, lw=lw(1.8)),
                      bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                                edgecolor=color, alpha=0.92))

    ax_r.set_xlim(lim)
    ax_r.set_ylim(lim)
    ax_r.set_xlabel(f"{labels[0]} rank")
    ax_r.set_ylabel(f"{labels[1]} rank")
    ax_r.set_aspect("equal", adjustable="box")
    ax_r.spines["top"].set_visible(False)
    ax_r.spines["right"].set_visible(False)
    ax_r.tick_params(labelsize=FS_TICK)
    ax_r.text(-0.155, 1.02, tags[1], transform=ax_r.transAxes, fontsize=FS_TAG,
              fontweight="bold", ha="left", va="bottom")

    # ------------------------------------------------------- legend + numbers
    # Only the movement classes that actually occur are listed. An entry reading
    # "(0 genotypes)" is a key to a colour that is nowhere on the figure, and it
    # invites the reader to look for movers that are not there.
    handles = [
        Line2D([0], [0], color=MOVE_COLORS[cls], lw=lw(3.4),
               label=f"|Δrank| {txt}   ({counts.get(cls, 0)} genotype"
                     f"{'' if counts.get(cls, 0) == 1 else 's'})")
        for cls, txt in (("stable", "≤ 2"), ("small", "3–5"), ("large", "≥ 6"))
        if counts.get(cls, 0) > 0
    ] + [
        Patch(facecolor=BAND_HIGH, alpha=0.20, label="highest-tillering 5"),
        Patch(facecolor=BAND_LOW, alpha=0.20, label="lowest-tillering 5"),
    ]
    ax_t.legend(handles=handles, loc="upper left", frameon=False,
                fontsize=FS_LEGEND, bbox_to_anchor=(0.02, 1.02),
                handlelength=1.6, labelspacing=0.75)

    stats = "\n".join([
        f"Spearman ρ = {rho:.3f}   (n = {n})",
        f"highest-tillering {EXTREME_N} shared: "
        f"{len(top_hand & top_fp)}/{EXTREME_N}",
        f"lowest-tillering {EXTREME_N} shared: "
        f"{len(bot_hand & bot_fp)}/{EXTREME_N}",
        f"identical or ±1 rank: "
        f"{int((df['rank_shift'].abs() <= 1).sum())}/{n}",
    ])
    ax_t.text(0.02, 0.62, stats, transform=ax_t.transAxes,
              ha="left", va="top", fontsize=FS_STAT, linespacing=1.85,
              bbox=dict(boxstyle="round,pad=0.6", facecolor="#F5F5F5",
                        edgecolor="#BBBBBB"))

    return dict(rho=rho, top=len(top_hand & top_fp), bot=len(bot_hand & bot_fp),
                movers=sorted(df[df["move"] == "large"]["apsim_name"]))


def save(fig, png: Path) -> None:
    # No bbox_inches="tight": the axes are placed in inches, so the file is
    # already exactly FIG_W wide and cropping it would only make Word rescale.
    fig.savefig(png, dpi=DPI)
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True)
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")


def report(st: dict) -> None:
    print(f"  rho {st['rho']:.3f}   top{EXTREME_N} {st['top']}/{EXTREME_N}"
          f"   bottom{EXTREME_N} {st['bot']}/{EXTREME_N}")
    print(f"  movers >=6: {st['movers']}")


def make_figure(df: pd.DataFrame, labels: tuple[str, str],
                fig_labels: tuple[tuple[str, str], tuple[str, str]],
                png: Path) -> None:
    set_style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    st = draw_pair(*row_axes(fig), df, labels, fig_labels)
    save(fig, png)
    report(st)


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
    for other in others_present(long):
        labels, fig_labels = pair_style(other)
        # the original pair keeps its filename so existing references still resolve
        stem = ("FigRank_GenotypeRanks_BothNetworks" if other == "flashp_suc"
                else f"FigRank_GenotypeRanks_{other}")
        print(f"\n=== {REFERENCE} vs {other} ===")
        make_figure(pair_frame(long, other), labels, fig_labels,
                    OUT / f"{stem}.png")


if __name__ == "__main__":
    main()
