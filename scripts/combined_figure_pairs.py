"""All-networks versions of the three pairwise figures -- one page instead of two.

The rank, tiller-number and slope figures are all inherently PAIRWISE: a
slopegraph, a dumbbell and a value-value scatter each compare the reference
network with exactly one other. With three arms that means two figures per
question, which makes the reader flip between pages to answer the only thing
they actually want to know -- do BOTH Flash-P networks land where the
hand-curated one does?

So this script stacks them. Each non-reference arm gets one row of the original
layout, panels tagged a/b then c/d down the page:

    a  b     hand-curated  vs  Flash-P (30 nodes, 44 edges)
    c  d     hand-curated  vs  Flash-P (69 nodes, 116 edges)

Nothing is redrawn here. The panels come from draw_pair() in
combined_figure_rank.py and combined_figure_values.py, the same functions the
single-pair figures use, and they are placed by those modules' own row_axes(),
so a stacked panel and its standalone twin are the same size to the inch. This
script only stacks the bands, assigns the panel tags and writes the heading that
says which Flash-P network the row is about -- both are called "Flash-P" and the
node/edge counts are what tell them apart.

Reads  03_runs/Arm_Comparison_long.csv
Writes 04_figures/figures_combined/FigRank_GenotypeRanks_AllNetworks.{png,svg}
       04_figures/figures_combined/FigPTN_TillerNumber_AllNetworks.{png,svg}
       04_figures/figures_combined/FigSlope_SupplyDemandSlope_AllNetworks.{png,svg}
"""

from __future__ import annotations

import string
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS as REGISTRY, REFERENCE, FIG_ROOT  # noqa: E402
import combined_figure_rank as rankfig        # noqa: E402
import combined_figure_values as valfig       # noqa: E402

OUT = FIG_ROOT / "figures_combined"

# One band of the single-pair layout per non-reference arm, stacked. The width
# and the band height come from those modules, so the stacked panels are the
# printed size the standalone ones were tuned at rather than a rescaling of it.
FIG_W = rankfig.FIG_W
DPI = rankfig.DPI

# Above each band, naming the network the row compares against. The single-pair
# figures do not need this -- there is only one row -- so it is the one piece of
# vertical budget this script adds.
HEAD_H = 0.50
GAP = 0.26          # one band's x-axis label to the next band's heading
FS_HEAD = 10.0


def _figure(row_h: float, k: int):
    """A k-band figure and the inch offset of each band's top."""
    pitch = HEAD_H + row_h + GAP
    h = rankfig.PAD_T + k * pitch - GAP + rankfig.PAD_B
    fig = plt.figure(figsize=(FIG_W, h))
    tops = [rankfig.PAD_T + i * pitch + HEAD_H for i in range(k)]
    return fig, tops


def _tags(row: int) -> tuple[str, str]:
    """('a','b') for the first arm, ('c','d') for the second, and so on."""
    return (string.ascii_lowercase[2 * row], string.ascii_lowercase[2 * row + 1])


def _heading(fig, top_in: float, other: str, keyed: bool) -> None:
    """Which comparison this band is, set in the gap above it.

    keyed says whether the arm colour is a key on this figure. On the value
    figures it is -- it is the colour of that row's Flash-P dots -- so the
    heading is set in it. On the rank figure nothing is drawn in the arm colour
    (lines are coloured by how far a genotype moved), so a coloured heading
    would be a key to nothing, and one of the two arm colours is close enough to
    the low-tillering cultivar's green to be read as it. There it is dark grey.
    """
    h = fig.get_size_inches()[1]
    fig.text(rankfig.PAD_L / FIG_W, 1.0 - (top_in - 0.26) / h,
             f"{REGISTRY[REFERENCE].label}  vs  {REGISTRY[other].fig_label}",
             fontsize=FS_HEAD, fontweight="bold",
             color=REGISTRY[other].color if keyed else "#333333",
             ha="left", va="bottom")


def build_rank(long, others: list[str]) -> None:
    rankfig.set_style()
    fig, tops = _figure(rankfig.ROW_H, len(others))
    for i, (other, top) in enumerate(zip(others, tops)):
        labels, fig_labels = rankfig.pair_style(other)
        st = rankfig.draw_pair(*rankfig.row_axes(fig, top),
                               rankfig.pair_frame(long, other),
                               labels, fig_labels, _tags(i))
        _heading(fig, top, other, keyed=False)
        print(f"  {_tags(i)[0]}/{_tags(i)[1]}  {REFERENCE} vs {other}")
        rankfig.report(st)
    _save(fig, OUT / "FigRank_GenotypeRanks_AllNetworks.png")


def build_values(long, others: list[str], metric: dict) -> None:
    valfig.set_style()
    fig, tops = _figure(valfig.ROW_H, len(others))
    for i, (other, top) in enumerate(zip(others, tops)):
        labels, fig_labels, colors = valfig.pair_style(other)
        st = valfig.draw_pair(*valfig.row_axes(fig, top),
                              valfig.pair_frame(long, other), metric,
                              labels, fig_labels, colors, _tags(i))
        _heading(fig, top, other, keyed=True)
        print(f"  {_tags(i)[0]}/{_tags(i)[1]}  {REFERENCE} vs {other}")
        valfig.report(metric, st)
    _save(fig, OUT / f"{metric['stem']}_AllNetworks.png")


def _save(fig, png: Path) -> None:
    # No bbox_inches="tight", for the same reason as the single-pair figures:
    # the bands are placed in inches, so the file is already exactly FIG_W wide.
    fig.savefig(png, dpi=DPI)
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True)
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")


def main() -> None:
    long = valfig.load_long()
    others = valfig.others_present(long)
    if len(others) < 2:
        # One non-reference arm means the stacked figure IS the single-pair
        # figure; writing it again under a second name would just be a duplicate
        # the reader has to reconcile.
        print(f"only one non-reference arm ({others[0]}) -- "
              f"nothing to stack, the pairwise figures already say it all")
        return

    OUT.mkdir(parents=True, exist_ok=True)
    print(f"\n=== ranks, {len(others)} arms vs {REFERENCE} ===")
    build_rank(long, others)
    for m in valfig.METRICS:
        print(f"\n=== {m['tag']}, {len(others)} arms vs {REFERENCE} ===")
        build_values(long, others, m)


if __name__ == "__main__":
    main()
