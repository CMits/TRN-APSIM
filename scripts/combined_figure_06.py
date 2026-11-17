"""How often does low tillering out-yield high tillering? One bar per network.

  04_figures/figures_combined/Fig06_LowTilleringWins_AllNetworks.{png,svg}

Three sowing dates across the x-axis, one bar per network inside each. The bar is
the percentage of simulated seasons in which the low-tillering cultivar (nced3)
produced a higher grain yield than the high-tillering one (max2), at 2.5 pl/m2.

WHY NOT A STACK OF THE PER-ARM Fig06_Win_Rate.png
Because that figure is BYTE-IDENTICAL for every arm -- nced3 clamps to zero
potential tillers and max2 saturates far above the ~3.3 fertile tillers the crop
can realise, so both contrast cultivars sit outside the band where the network
matters (README, "What this figure cannot show"). Stacking three copies of the
same panel would invite the reader to see agreement the figure cannot evidence.
One bar per network says the same thing in one row, and the agreement line
printed under the axes says it in words.

WHY THE NUMBERS DIFFER FROM Fig06_Win_Rate.png
That figure reports 67 / 59 / 69 %; this one reports 64 / 56 / 66 %. It keeps a
season in which only ONE of the two cultivars was harvested and hands that
season's "win" to whichever has data, unopposed -- 3 such seasons at Oct 15, 3 at
Nov 15, 5 at Dec 15, and 28 = 25 + 3. An axis that asserts "low tillering
out-yielded high tillering" cannot count a season where there is nothing to
out-yield, so this figure requires both (dropna(subset=[HIGH, LOW])).

That is also the definition compare_crossover.py already uses, so these bars are
the picture of the %win Low column of 03_runs/Crossover_Comparison.txt -- 64 / 56
/ 66 over n = 39 / 43 / 44, for every arm. Diff the printed report against that
file to check a run.

COLOUR MEANS NETWORK HERE, NOT CULTIVAR
The per-arm figure colours its two bars red (High) and green (Low). This one has
one bar per network, so it uses the registry arm colours -- the same ones that key
the networks on FigPTN and FigSlope. HT and LT appear only as text in the footer,
in grey, so no reader maps a bar colour onto a cultivar.

Reads  03_runs/apsim/<arm>/DensityManagement_NoTiller.db  (via combined_figure_05)
Writes 04_figures/figures_combined/Fig06_LowTilleringWins_AllNetworks.{png,svg}
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS as REGISTRY, FIG_ROOT  # noqa: E402

# The loader, the contrast cultivars and the sowing list come from the crossover
# figure rather than being restated here, so the two figures cannot end up
# disagreeing about which rows are in scope. load() issues a named-column query
# rather than SELECT * -- these are ~87 MB databases, one per arm.
from combined_figure_05 import (  # noqa: E402
    load, RUNS, HIGH, LOW, SOWINGS, SOWING_LABELS)

OUT = FIG_ROOT / "figures_combined"

# Registry order, restricted to arms whose APSIM run actually exists. Rebuilt here
# rather than imported because this figure needs .color and .fig_label_2line,
# which combined_figure_05's (name, label) tuples do not carry.
ARMS = [a for a in REGISTRY.values()
        if (RUNS / a.name / "DensityManagement_NoTiller.db").exists()]

# ----------------------------------------------------------------- print size
# Built at the width it OCCUPIES on the page, so Word places it at 100% and the
# point sizes below are the point sizes the reader gets. Same discipline as
# combined_figure_rank.py, including the rect() helper and the absence of both
# tight_layout and bbox_inches="tight".
#
# Every inch below was measured against DejaVu Sans, matplotlib's default and what
# this environment resolves to. A font.family override would invalidate the left
# margin and the value-label clearances.
FIG_W = 6.5                      # A4 text block at 2 cm margins
DPI = 600

LEFT = 0.75                      # two lines of rotated y-label, "100", tick pad
PAD_R, PAD_T, PAD_B = 0.10, 0.05, 0.06
LEG_H, LEG_GAP = 0.40, 0.06      # the legend gets its own band -- see below
AX_H = 2.55
XTICK_H = 0.26                   # "Oct 15" at FS_XTICK, plus tick length and pad
FOOT_H = 0.32                    # two footer lines at FS_FOOT
NOTE_H = 0.22                    # the agreement line, below the footer

AX_W = FIG_W - LEFT - PAD_R
FIG_H = (PAD_T + LEG_H + LEG_GAP + AX_H + XTICK_H + FOOT_H + NOTE_H + PAD_B)

# Bars span 76% of a group; the remaining 24% is the gap between groups. At three
# arms that is 0.06 in inside a group against 0.45 in between them -- the 8:1
# ratio is what makes it read as three groups rather than nine bars.
GROUP_FRAC = 0.76
BAR_FRAC = 0.88                  # of the slot, leaving a hairline between bars

# True printed points.
FS_AXIS = 10.5                   # y-axis label, set on two lines
FS_XTICK = 11.0                  # the sowing dates read as group headings
FS_TICK = 10.0
FS_VALUE = 12.0                  # bold, above each bar
FS_LEG = 8.5
FS_FOOT = 7.5
FS_NOTE = 8.0

# The value label is the percentage ALONE. "64%\n(25/39)" at FS_VALUE overlaps its
# neighbour by 0.20 in at three arms, and the denominator is a property of the
# sowing date rather than of the network, so it belongs in the footer where it is
# printed once. "64%" clears by 0.08 in at three arms; a value above ~90% would
# need FS_VALUE dropped to 11, which the data (max 66%) does not approach.
GRID = "#DDDDDD"
FOOT_COLOR = "#555555"


def rect(w: float, h: float, x: float, y: float, dx: float, dy: float) -> list:
    """Axes rectangle in figure coordinates. y is measured DOWN from the top."""
    return [x / w, 1.0 - (y + dy) / h, dx / w, dy / h]


def win_rate(h: pd.DataFrame, sow: str) -> tuple[int, int, float, float, float]:
    """(wins, n, percent, mean HT, mean LT) for one sowing date.

    dropna(subset=...) rather than how="all": a season in which only one cultivar
    was harvested cannot answer "did low out-yield high", so it is not a season.
    Strictly ">", matching compare_crossover.py -- there are no exact ties in the
    data, but the convention should be visible rather than inferred.
    """
    piv = (h[h["SowingDate"] == sow]
           .pivot_table(index="Year", columns="Genotype",
                        values="GrainYield_tha", aggfunc="first")
           .dropna(subset=[HIGH, LOW]))
    if piv.empty:
        raise SystemExit(
            f"{sow}: no season has BOTH {HIGH} and {LOW} -- the win rate would be "
            f"undefined. Check that every arm's APSIM run produced both cultivars.")
    wins = int((piv[LOW] > piv[HIGH]).sum())
    return (wins, len(piv), 100.0 * wins / len(piv),
            float(piv[HIGH].mean()), float(piv[LOW].mean()))


def gather() -> dict[str, dict[str, tuple]]:
    """arm -> sowing -> win_rate(), one database read per arm."""
    if not ARMS:
        raise SystemExit(
            f"no APSIM output under {RUNS} -- run build_apsim_runs.py and APSIM first")
    return {a.name: {sow: win_rate(load(a.name), sow) for sow in SOWINGS}
            for a in ARMS}


def seasons_per_sowing(stats: dict) -> dict[str, int]:
    """The shared denominator, asserted to be the same in every arm.

    The footer prints one "n seasons" per group rather than one per bar, which is
    only honest while the arms agree. They do today (39 / 43 / 44 everywhere) --
    the harvest failures are a property of the weather, not of the network.
    """
    out = {}
    for sow in SOWINGS:
        ns = {stats[a.name][sow][1] for a in ARMS}
        if len(ns) != 1:
            raise SystemExit(
                f"{sow}: arms disagree on the number of paired seasons ({sorted(ns)}). "
                f"The shared footer would be wrong -- move the count back onto each bar.")
        out[sow] = ns.pop()
    return out


def draw(ax, axl, stats: dict, seasons: dict[str, int]) -> None:
    slot = GROUP_FRAC / len(ARMS)

    ax.axhline(50, linestyle=(0, (4, 3)), color="#888888", linewidth=0.9, zorder=1)
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)

    for j, a in enumerate(ARMS):
        xs, ys = [], []
        for i, sow in enumerate(SOWINGS):
            xs.append(i - GROUP_FRAC / 2 + slot * (j + 0.5))
            ys.append(stats[a.name][sow][2])
        ax.bar(xs, ys, width=slot * BAR_FRAC, color=a.color, zorder=3,
               edgecolor="black", linewidth=0.5)
        for x, y in zip(xs, ys):
            ax.text(x, y + 1.8, f"{y:.0f}%", ha="center", va="bottom",
                    fontsize=FS_VALUE, fontweight="bold", zorder=4)

    ax.set_xlim(-0.5, len(SOWINGS) - 0.5)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(range(len(SOWINGS)))
    ax.set_xticklabels([SOWING_LABELS[s] for s in SOWINGS],
                       fontsize=FS_XTICK, fontweight="bold")
    ax.tick_params(axis="y", labelsize=FS_TICK)
    ax.tick_params(axis="x", length=0, pad=4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Two lines, not one: the label is 3.57 in set on a single line, against a
    # 2.55 in axes height to rotate into.
    ax.set_ylabel("Years low tillering out-yielded\nhigh tillering (%)",
                  fontsize=FS_AXIS, linespacing=1.4)

    # The 50% line, named in the gap BETWEEN the first two groups -- the only
    # place at that height where no bar can reach, at any arm count, since the
    # groups always span the same 76% of their slot.
    ax.text(0.5, 51.5, "even split", ha="center", va="bottom",
            fontsize=FS_FOOT, style="italic", color="#888888", zorder=5,
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                      edgecolor="none"))

    # Footer: the denominator once per group, then the means it is taken over. HT
    # and LT appear only here, and only in grey -- the bars are keyed to networks.
    ref = ARMS[0].name
    for i, sow in enumerate(SOWINGS):
        _, _, _, m_high, m_low = stats[ref][sow]
        ax.text(i, -0.10, f"{seasons[sow]} seasons\n"
                          f"mean yield   HT {m_high:.2f}   LT {m_low:.2f} t/ha",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=FS_FOOT, color=FOOT_COLOR, linespacing=1.55)

    axl.axis("off")
    # ncol = one column per arm. A fourth arm would overflow 6.5 in and needs two
    # rows; fig_label_2line already keeps the node counts off the first line,
    # which is what makes three fit at all.
    axl.legend(handles=[Patch(facecolor=a.color, edgecolor="black", linewidth=0.5,
                              label=a.fig_label_2line) for a in ARMS],
               loc="center", ncol=len(ARMS), frameon=False, fontsize=FS_LEG,
               handlelength=1.3, handletextpad=0.6, columnspacing=1.6)


def agreement_note(fig, stats: dict) -> str:
    """State the cross-arm spread on the figure, in words.

    Bars of identical height in different colours look like a plotting bug. They
    are in fact the result, so the figure has to say so itself rather than leaving
    it to a caption that may not travel with it.
    """
    spread = max(max(stats[a.name][sow][2] for a in ARMS)
                 - min(stats[a.name][sow][2] for a in ARMS) for sow in SOWINGS)
    note = (f"The {len(ARMS)} networks give identical values at every sowing date."
            if spread < 0.05 else
            f"Largest difference between networks: {spread:.1f} percentage points.")
    fig.text(LEFT / FIG_W, (PAD_B + NOTE_H / 2) / FIG_H, note,
             ha="left", va="center",
             fontsize=FS_NOTE, style="italic", color=FOOT_COLOR)
    return note


def save(fig, png: Path) -> None:
    # No bbox_inches="tight": the axes are placed in inches, so the file is
    # already exactly FIG_W wide and cropping it would only make Word rescale.
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=DPI)
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True)
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")


def report(stats: dict, seasons: dict[str, int], note: str) -> None:
    print(f"\n  {'sowing':<8}{'arm':<14}{'wins':>6}{'n':>5}{'%':>6}"
          f"{'meanHT':>9}{'meanLT':>9}")
    for sow in SOWINGS:
        for a in ARMS:
            w, n, p, mh, ml = stats[a.name][sow]
            print(f"  {SOWING_LABELS[sow]:<8}{a.name:<14}{w:>6}{n:>5}{p:>5.0f}%"
                  f"{mh:>9.2f}{ml:>9.2f}")
    print(f"\n  {note}")
    print("  cross-check: the %win Low column of 03_runs/Crossover_Comparison.txt")


def main() -> None:
    stats = gather()
    seasons = seasons_per_sowing(stats)

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    r = lambda x, y, dx, dy: rect(FIG_W, FIG_H, x, y, dx, dy)  # noqa: E731
    ax = fig.add_axes(r(LEFT, PAD_T + LEG_H + LEG_GAP, AX_W, AX_H))
    axl = fig.add_axes(r(LEFT, PAD_T, AX_W, LEG_H))

    draw(ax, axl, stats, seasons)
    note = agreement_note(fig, stats)
    save(fig, OUT / "Fig06_LowTilleringWins_AllNetworks.png")
    report(stats, seasons, note)


if __name__ == "__main__":
    main()
