"""Methods figures 4 and 5 -- the supply/demand mapping, Arabidopsis to sorghum.

  04_figures/methods/FigM04_Arabidopsis_SupplyDemand.(png|svg)
  04_figures/methods/FigM05_Sorghum_Tillering_SupplyDemand.(png|svg)
  04_figures/methods/FigM04_05_SupplyDemand_Combined.(png|svg)
  04_figures/methods/FigM06_SupplyDemand_FromZero.(png|svg)

FigM01 shows Step 1 (the networks' sucrose response) and FigM03 the Step-2 anchor.
These two draw the back half of the chain: the anchored Arabidopsis response placed
on APSIM's continuous supply/demand axis, and the sorghum tillering line that comes
out of it. One panel per network in both.

The third output is those same two as one 3x2 with a single shared legend --
Arabidopsis in the left column, sorghum in the right, one network per row, so the
transformation reads across a row. It is drawn from the same derive() result and the
same per-panel helpers as the two singles, which stay available separately.

They are the reproduction's equivalents of the paper's Fig04_Weighted_Branches_SD
and Fig05_APSIM_Lines_Final, drawn for all three networks rather than one.

TWO THINGS THE CAPTIONS SHOULD CARRY
1. There is NO sorghum wild type anywhere in this chain. What sets the sorghum
   scale is the min-max rescale in step 5, into the bands INT_BAND (-1.363, 1.5)
   and SLOPE_BAND (0, 0.4) taken from the published parameter set. That is why the
   intercept range comes out identical for all three networks: it is imposed, not
   derived.
2. For the same reason the global WT reference used in step 2 CANCELS. It is a
   factor common to every genotype, and a common factor cancels from both the
   numerator and the denominator of
       rescale(v) = lo + (v - min v)/(max v - min v) * (hi - lo)
   Substituting 1.0 or 7.5 for 3.0 branches moves every fitted intercept and slope
   by less than 2e-15. Only the genotypes' RELATIVE pattern reaches APSIM.

X RANGE
Both figures are drawn over the REALISED supply/demand range, 11.70 .. 42.36 --
what APSIM actually reaches in the published run (build_apsim_parameters.SD_REALISED),
rather than the wider 10.46 .. 47.48 the line is fitted between. Nothing is
distorted by that: b_hat is linear in S/D by construction, since step 5 fits a
straight line through b_low at SD_min and b_high at SD_max, so the drawn range is
an exact sub-segment.

THE ZERO CLAMP (FigM05)
APSIM computes TTN = max(0, intercept + slope * S/D), so a negative intercept is a
THRESHOLD, not a tiller count. The figure plots the CLAMPED value -- what the crop
model actually receives -- so a cultivar that never tillers lies flat on zero. Lines
that cross zero inside the range get the kink at the crossing, not a clamped
endpoint. 01_audit/sd_mapping_lines.csv keeps the unclamped intercept and slope.

Reads  00_reference/paper_data/{Genotype_Linear_Equations,weighted_arabidopsis_branches}.csv
       the arm registry, and each network via build_apsim_parameters
Writes 04_figures/methods/FigM04_Arabidopsis_SupplyDemand.{png,svg}
       04_figures/methods/FigM05_Sorghum_Tillering_SupplyDemand.{png,svg}
       04_figures/methods/FigM04_05_SupplyDemand_Combined.{png,svg}
       04_figures/methods/FigM06_SupplyDemand_FromZero.{png,svg}
       01_audit/sd_mapping_lines.csv
"""

from __future__ import annotations

import csv
import string
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "00_reference" / "psoup_source"))

import build_apsim_parameters as B          # noqa: E402
from arms import ARMS, FIG_ROOT              # noqa: E402
from methods_figure_sucrose import (         # noqa: E402
    POLYCHROME, pretty_name, eligible_genotypes, published_genotypes,
    stack_geometry,
)

PD = ROOT / "00_reference" / "paper_data"
AUDIT = ROOT / "01_audit"
OUT = FIG_ROOT / "methods"

# True printed point sizes -- the figure is built at the width it occupies on the
# page (FIG_W), so Word does not rescale it. Same settings as FigM01.
# Height comes from stack_geometry(), shared with FigM01, so all three figures are
# the same shape and shorten together when the arm set is restricted.
FIG_W = 6.5
FS_AXIS = 10.0
FS_TICK = 9.0
FS_TITLE = 10.0
FS_TAG = 12.0
FS_LEGEND = 8.0

# The combined Nx2 has its own inch budget: shorter panel rows than the stacks,
# because only the bottom row needs an x-label.
BOTTOM_IN_COMBINED = 0.495
ROW_IN_COMBINED = 2.25
COLUMN_HEADS = ("Arabidopsis", "Sorghum")

# Five legend columns rather than the stacks' four. The genotype labels are short
# enough that five fit the 6.5 in width once the handles and column gaps are pulled
# in (see LEGEND_KW), and five columns put the 32 genotypes on seven rows instead of
# eight -- a row of legend bought back for the panels.
LEGEND_NCOL_COMBINED = 5
LEGEND_ROW_IN = 0.155         # one legend row at FS_LEGEND
HEAD_CLEAR_IN = 0.515         # legend to the column headers

# One y-scale down the Arabidopsis column, so the networks are directly
# comparable. Flip to False for free axes there, as in standalone FigM04.
COMBINED_SHARE_Y_ARABIDOPSIS = True


def combined_geometry(nrow: int, legend_rows: int) -> tuple[float, float, float]:
    """(figure height, gridspec top, gridspec bottom) for the Nx2 combined figure.

    The legend band is sized from the rows it actually needs rather than being a
    constant, so changing LEGEND_NCOL_COMBINED reflows the whole figure instead of
    leaving a gap or crowding the column headers.
    """
    legend_in = HEAD_CLEAR_IN + LEGEND_ROW_IN * legend_rows
    h = legend_in + BOTTOM_IN_COMBINED + ROW_IN_COMBINED * nrow
    return h, 1.0 - legend_in / h, BOTTOM_IN_COMBINED / h

SD_LO, SD_HI = B.SD_REALISED                 # 11.70 .. 42.36

# FigM04_05 is drawn twice, from the same derive() result and the same panel
# helpers, differing only in this flag:
#
#   from_zero=False   the realised range 11.70 .. 42.36, endpoints ticked, and the
#                     CLAMPED tiller number -- what APSIM actually receives.
#   from_zero=True    the axis extended back to S/D 0 and nothing clamped, so both
#                     parameters in Table 1 are visible on the page: the intercept
#                     where each line meets the y-axis, the slope as its gradient.
#
# The second exists because the first cannot answer "why does a genotype with
# intercept 1.088 start at 5.8 tillers?" -- the answer is that 11.70 is not zero
# and the slope term has already contributed slope * 11.70 by the left edge.


def set_sd_axis(ax, from_zero: bool) -> None:
    """Tick the endpoints, always.

    With matplotlib's automatic ticks the leftmost one is 20 and the axis simply
    begins at 11.70 with nothing to mark it, so it reads as an axis starting at
    zero -- which makes every line look as though it starts at its intercept.

    The from-zero figure starts at a true zero, so that problem does not arise
    there: it takes plain ticks every 5 out to 45. The lines still stop at SD_HI.
    """
    if from_zero:
        ax.set_xlim(0, 45)
        ax.set_xticks(range(0, 46, 5))
        return
    ticks = [SD_LO, 20, 30, SD_HI]
    ax.set_xlim(SD_LO, SD_HI)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.1f}" if t in (SD_LO, SD_HI) else f"{t:g}"
                        for t in ticks])


def derive() -> tuple[dict, list[str], float, float]:
    """{arm: {genotype: step-3-to-5 record}} plus the genotype list and SD bounds.

    Everything comes from build_apsim_parameters itself -- Step 1 per arm, then the
    same steps_3_to_5 the pipeline runs -- so these figures cannot drift from the
    parameters actually emitted.
    """
    lin = list(csv.DictReader((PD / "Genotype_Linear_Equations.csv").open(encoding="utf-8-sig")))
    sd_min, sd_max = float(lin[0]["SD_min"]), float(lin[0]["SD_max"])
    anchor = {r["genotype"]: float(r["medium"]) for r in
              csv.DictReader((PD / "weighted_arabidopsis_branches.csv").open(encoding="utf-8-sig"))}
    genotypes = eligible_genotypes()

    out = {}
    for name, a in ARMS.items():
        resp = (B.psoup_response(a.response_csv, a.phenotype) if a.kind == "r_psoup"
                else B.network_response(a.equations, a.phenotype, a.sucrose_node, a.gene_map))
        out[name] = B.steps_3_to_5(resp, anchor, genotypes, sd_min, sd_max)
    return out, genotypes, sd_min, sd_max


def colour_map() -> dict[str, str]:
    """Indexed within the PUBLISHED 34, so a genotype keeps one colour everywhere."""
    ordered = sorted(published_genotypes(), key=pretty_name)
    return {g: POLYCHROME[i % len(POLYCHROME)] for i, g in enumerate(ordered)}


def bhat_at(rec: dict, sd: float, sd_min: float, sd_max: float) -> float:
    """b_hat interpolated on the S/D axis.

    Exact, not approximate: step 5 defines b_hat as the straight line through
    b_low at SD_min and b_high at SD_max, so this is that same line evaluated.
    """
    t = (sd - sd_min) / (sd_max - sd_min)
    return rec["b_low"] + t * (rec["b_high"] - rec["b_low"])


def ttn_at(rec: dict, sd: float) -> float:
    """The APSIM line BEFORE the clamp -- intercept is a threshold, not a count."""
    return rec["intercept"] + rec["slope"] * sd


def set_style() -> None:
    plt.rcParams.update({"font.size": FS_TICK,
                         "axes.labelsize": FS_AXIS,
                         "xtick.labelsize": FS_TICK,
                         "ytick.labelsize": FS_TICK})


def save_figure(fig, png: Path) -> None:
    """600 dpi PNG plus a transparent SVG, both cropped to the drawn content."""
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")


def new_stack():
    """One stacked panel per arm, sharing the x-axis, legend across the top.

    Stacked rather than side by side because the figure is built at 6.5 in: three
    panels in a row would leave under 2 in each. Same layout as FigM01.
    """
    nrow = len(ARMS)
    fig_h, top, bottom = stack_geometry(nrow)
    fig = plt.figure(figsize=(FIG_W, fig_h))
    gs = fig.add_gridspec(nrow, 1, hspace=0.34, top=top, bottom=bottom,
                          left=0.135, right=0.985)
    return fig, [fig.add_subplot(gs[i, 0]) for i in range(nrow)]


def finish(fig, axes, handles, xlabel, png: Path) -> None:
    for i, ax in enumerate(axes):
        set_sd_axis(ax, from_zero=False)
        if i == len(axes) - 1:
            ax.set_xlabel(xlabel)
        else:
            ax.tick_params(labelbottom=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=FS_TICK)

    leg = fig.legend(handles=handles, loc="upper center", ncol=4,
                     fontsize=FS_LEGEND, frameon=False, handlelength=1.4,
                     columnspacing=1.0, labelspacing=0.32, handletextpad=0.4,
                     borderaxespad=0.0, bbox_to_anchor=(0.5, 1.0))
    for t in leg.get_texts():
        t.set_fontstyle("italic")

    save_figure(fig, png)


def draw_arabidopsis_panel(ax, rows, ordered, colours, sd_min, sd_max,
                           from_zero: bool = False) -> list:
    """One FigM04 panel -- b_hat for every genotype, straight over the S/D range.

    Shared with the combined figure so the two cannot draw the same quantity in
    two different ways.
    """
    x0 = 0.0 if from_zero else SD_LO
    hs = []
    for g in ordered:
        y = [bhat_at(rows[g], x0, sd_min, sd_max),
             bhat_at(rows[g], SD_HI, sd_min, sd_max)]
        (ln,) = ax.plot([x0, SD_HI], y, color=colours[g], lw=1.0,
                        label=pretty_name(g))
        hs.append(ln)
    return hs


def draw_sorghum_panel(ax, rows, ordered, colours, from_zero: bool = False) -> list:
    """One FigM05 panel -- the sorghum tillering line, clamped or raw.

    from_zero=False is exactly what APSIM uses, TTN = max(0, intercept + slope *
    S/D). A line that crosses zero gets the kink at the crossing rather than a
    clamped endpoint, so the flat and rising parts meet where they actually do. A
    cultivar negative throughout sits flat on zero -- which is the truth about it:
    it never tillers.

    from_zero=True draws the RAW line instead, nothing clamped, extended back to
    S/D 0. The clamp is what hides the intercept: it flattens every negative one
    onto the axis, and a reader who then compares the figure with Table 1 cannot
    find the parameter anywhere on the page. Here the line meets the y-axis at its
    intercept, its gradient is its slope, and where it runs below zero it is
    showing the threshold the crop has to clear before it tillers at all. Read it
    as the parameters; read the clamped figure for what the crop receives.
    """
    hs = []
    for g in ordered:
        rec = rows[g]
        y0, y1 = ttn_at(rec, SD_LO), ttn_at(rec, SD_HI)
        if from_zero:
            xs, ys = [0.0, SD_HI], [ttn_at(rec, 0.0), y1]
        elif y0 < 0 and y1 < 0:
            xs, ys = [SD_LO, SD_HI], [0.0, 0.0]
        elif y0 >= 0 and y1 >= 0:
            xs, ys = [SD_LO, SD_HI], [y0, y1]
        else:
            cross = SD_LO + (SD_HI - SD_LO) * (0 - y0) / (y1 - y0)
            xs = [SD_LO, cross, SD_HI]
            ys = [max(y0, 0.0), 0.0, max(y1, 0.0)]
        (ln,) = ax.plot(xs, ys, color=colours[g], lw=1.0, label=pretty_name(g))
        hs.append(ln)
    if not from_zero:
        # The zero rule marks the clamp floor, so it only means something where the
        # clamp is applied. In the from-zero figure the lines cross it freely and it
        # would read as a second, unexplained line.
        ax.axhline(0.0, color="#666666", lw=0.7, zorder=0)
    return hs


def _padded(values) -> tuple[float, float]:
    pad = 0.06 * (max(values) - min(values))
    return min(values) - pad, max(values) + pad


def sorghum_ylim(derived, genotypes, from_zero: bool = False) -> tuple[float, float]:
    """One scale for all three arms: after the step-5 rescale they land in the same band.

    from_zero has to reach below zero: nothing is clamped there, so the axis has
    to hold the lowest intercept (-1.363, the bottom of the published band).
    """
    if from_zero:
        vals = [ttn_at(derived[n][g], sd) for n in ARMS
                for g in genotypes for sd in (0.0, SD_LO, SD_HI)]
    else:
        vals = [max(ttn_at(derived[n][g], sd), 0.0) for n in ARMS
                for g in genotypes for sd in (SD_LO, SD_HI)]
    return _padded(vals)


def arabidopsis_ylim(derived, genotypes, sd_min, sd_max,
                     from_zero: bool = False) -> tuple[float, float]:
    """The same for b_hat -- used ONLY by the combined figure.

    Standalone FigM04 keeps free y-axes (see figure_arabidopsis). In the combined
    layout the point is comparing rows, so the column shares one scale; the
    narrowest arm then uses about 60% of its panel height, which still reads.
    Set COMBINED_SHARE_Y_ARABIDOPSIS to False to go back to free axes there.
    """
    sds = (0.0, SD_LO, SD_HI) if from_zero else (SD_LO, SD_HI)
    return _padded([bhat_at(derived[n][g], sd, sd_min, sd_max) for n in ARMS
                    for g in genotypes for sd in sds])


def figure_arabidopsis(derived, genotypes, colours, sd_min, sd_max, png: Path) -> None:
    """FigM04 -- the anchored Arabidopsis response, placed on the S/D axis."""
    set_style()
    ordered = sorted(genotypes, key=pretty_name)
    fig, axes = new_stack()
    handles = None
    for ax, (name, a), tag in zip(axes, ARMS.items(), string.ascii_lowercase):
        hs = draw_arabidopsis_panel(ax, derived[name], ordered, colours,
                                    sd_min, sd_max)
        handles = handles or hs
        # Free y-axes: the three arms span very different ranges here, and forcing
        # one scale would flatten two panels.
        ax.set_ylabel("Weighted shoot branches\n(Arabidopsis units)")
        ax.set_title(a.fig_label, fontsize=FS_TITLE, pad=5)
        ax.text(-0.135, 1.02, tag, transform=ax.transAxes, fontsize=FS_TAG,
                ha="left", va="bottom")
    finish(fig, axes, handles, "Supply/demand ratio (APSIM)", png)


def figure_sorghum(derived, genotypes, colours, png: Path) -> None:
    """FigM05 -- the sorghum tillering lines APSIM receives."""
    set_style()
    ordered = sorted(genotypes, key=pretty_name)
    fig, axes = new_stack()

    # Shared y here, unlike FigM04: after the step-5 rescale every arm lands in the
    # same band, so a common scale is meaningful and makes the panels comparable.
    ylim = sorghum_ylim(derived, genotypes)

    handles = None
    for ax, (name, a), tag in zip(axes, ARMS.items(), string.ascii_lowercase):
        hs = draw_sorghum_panel(ax, derived[name], ordered, colours)
        handles = handles or hs
        ax.set_ylim(*ylim)
        ax.set_ylabel("Potential tiller number")
        ax.set_title(a.fig_label, fontsize=FS_TITLE, pad=5)
        ax.text(-0.135, 1.02, tag, transform=ax.transAxes, fontsize=FS_TAG,
                ha="left", va="bottom")
    finish(fig, axes, handles, "Supply/demand ratio (APSIM)", png)


def figure_combined(derived, genotypes, colours, sd_min, sd_max, png: Path,
                    from_zero: bool = False) -> None:
    """FigM04 and FigM05 as one Nx2 -- the transformation read across each row.

    Left column the anchored Arabidopsis response on APSIM's S/D axis, right column
    the sorghum line the step-5 rescale turns it into, one row per network. Both
    columns colour genotypes through the same colour_map, which indexes within the
    published 34, so one legend serves every panel.

    Built at FIG_W, the width it occupies on the page, so Word places it at 100%
    and the FS_* constants above are the printed point sizes.
    """
    set_style()
    ordered = sorted(genotypes, key=pretty_name)

    nrow = len(ARMS)
    legend_rows = -(-len(ordered) // LEGEND_NCOL_COMBINED)   # ceil
    fig_h, top, bottom = combined_geometry(nrow, legend_rows)
    fig = plt.figure(figsize=(FIG_W, fig_h))
    # hspace is tighter than new_stack's because only the bottom row carries an
    # x-label; left is wider because that margin holds the row label, the column's
    # y-label and its tick labels side by side.
    # top leaves the legend clear of the column headers, which sit a title pad
    # above the first row rather than tucked under the legend.
    gs = fig.add_gridspec(nrow, 2, hspace=0.20, wspace=0.34, top=top,
                          bottom=bottom, left=0.155, right=0.985)
    axes = [[fig.add_subplot(gs[r, c]) for c in range(2)] for r in range(nrow)]

    # The y-label goes on the middle row where there is one. With two rows there is
    # no middle, and labelling both reads better than hanging it off one of them.
    label_rows = {nrow // 2} if nrow >= 3 else set(range(nrow))

    ylim_a = (arabidopsis_ylim(derived, genotypes, sd_min, sd_max, from_zero)
              if COMBINED_SHARE_Y_ARABIDOPSIS else None)
    ylim_s = sorghum_ylim(derived, genotypes, from_zero)

    handles = None
    for r, (name, a) in enumerate(ARMS.items()):
        rows = derived[name]
        ax_a, ax_s = axes[r]
        hs = draw_arabidopsis_panel(ax_a, rows, ordered, colours, sd_min,
                                    sd_max, from_zero)
        draw_sorghum_panel(ax_s, rows, ordered, colours, from_zero)
        handles = handles or hs
        if ylim_a is not None:
            ax_a.set_ylim(*ylim_a)
        ax_s.set_ylim(*ylim_s)

        for c, ax in enumerate((ax_a, ax_s)):
            set_sd_axis(ax, from_zero)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=FS_TICK)
            if r == 0:
                # pad clears the panel tag, which sits at 1.02 in axes fraction.
                ax.set_title(COLUMN_HEADS[c], fontsize=FS_TITLE, pad=15)
            if r == len(axes) - 1:
                ax.set_xlabel("Supply/demand ratio (APSIM)")
            else:
                ax.tick_params(labelbottom=False)
            # a, b, c ... COLUMN-major: down the Arabidopsis column first,
            # then down the sorghum column.
            ax.text(-0.17, 1.02, string.ascii_lowercase[c * nrow + r],
                    transform=ax.transAxes,
                    fontsize=FS_TAG, ha="left", va="bottom")

        # One y-label per column where the stack is tall enough: the column shares
        # its scale, so repeating the label on every row would only cost margin.
        if r in label_rows:
            # Short here, unlike standalone FigM04: the column header already says
            # Arabidopsis, so the units qualifier is saying it twice.
            ax_a.set_ylabel("Shoot branches")
            ax_s.set_ylabel("Potential tiller number")

        # The network is a property of the ROW, so it is labelled once at the far
        # left rather than twice as a panel title. The position is exact because the
        # gridspec is placed explicitly and never re-laid out.
        pos = ax_a.get_position()
        fig.text(0.030, (pos.y0 + pos.y1) / 2, a.fig_label_2line, rotation=90,
                 fontsize=FS_TITLE, ha="center", va="center", linespacing=1.4)

    # Tighter handles and column gaps than the stacks' 4-column legend, which is
    # what buys the fifth column inside the same 6.5 in.
    leg = fig.legend(handles=handles, loc="upper center",
                     ncol=LEGEND_NCOL_COMBINED,
                     fontsize=FS_LEGEND, frameon=False, handlelength=1.1,
                     columnspacing=0.7, labelspacing=0.32, handletextpad=0.3,
                     borderaxespad=0.0, bbox_to_anchor=(0.5, 1.0))
    for t in leg.get_texts():
        t.set_fontstyle("italic")

    save_figure(fig, png)

def write_audit(derived, genotypes, sd_min, sd_max) -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    path = AUDIT / "sd_mapping_lines.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["arm", "genotype", "apsim_name", "b_low", "b_medium", "b_high",
                    "slope_raw", "intercept_raw", "intercept", "slope",
                    "ttn_at_sd_lo", "ttn_at_sd_hi", "clamped_over_realised_range"])
        for name in ARMS:
            rows = derived[name]
            for g in sorted(genotypes, key=B.apsim_name):
                r = rows[g]
                lo, hi = ttn_at(r, SD_LO), ttn_at(r, SD_HI)
                w.writerow([
                    name, g, B.apsim_name(g),
                    f"{r['b_low']:.6f}", f"{r['b_medium']:.6f}", f"{r['b_high']:.6f}",
                    f"{r['slope_raw']:.6f}", f"{r['intercept_raw']:.6f}",
                    f"{r['intercept']:.6f}", f"{r['slope']:.6f}",
                    f"{lo:.6f}", f"{hi:.6f}", lo < 0 and hi < 0,
                ])
    print(f"Saved: {path}")


def check(derived, genotypes) -> None:
    """Assertions that tie these figures to the parameters the pipeline emitted."""
    # The step-5 rescale forces these bands, identically in every arm. If this fails,
    # steps_3_to_5 was not driven the way build_apsim_parameters drives it.
    for name in ARMS:
        ints = [derived[name][g]["intercept"] for g in genotypes]
        slopes = [derived[name][g]["slope"] for g in genotypes]
        assert abs(min(ints) - B.INT_BAND[0]) < 1e-9 and abs(max(ints) - B.INT_BAND[1]) < 1e-9, name
        assert abs(min(slopes) - B.SLOPE_BAND[0]) < 1e-9 and abs(max(slopes) - B.SLOPE_BAND[1]) < 1e-9, name

    # And against the CSV the pipeline actually wrote, which stores 3 dp.
    for name in ARMS:
        emitted = {}
        p = ROOT / "03_runs" / f"Genotype_APSIM_Parameters_{name}.csv"
        if not p.exists():
            print(f"  (skipped cross-check: {p.name} not present)")
            continue
        for r in csv.DictReader(p.open(encoding="utf-8-sig")):
            emitted[r["Genotype"]] = (float(r["Intercept_Line"].rsplit("=", 1)[1]),
                                      float(r["Slope_Line"].rsplit("=", 1)[1]))
        worst = max(max(abs(derived[name][g]["intercept"] - emitted[B.apsim_name(g)][0]),
                        abs(derived[name][g]["slope"] - emitted[B.apsim_name(g)][1]))
                    for g in genotypes)
        assert worst < 1e-3, f"{name}: figure disagrees with emitted parameters by {worst}"
        # FigM06 asks the reader to read the intercept off the y-axis, so the value
        # the line is drawn through at S/D 0 has to BE the emitted intercept.
        assert all(abs(ttn_at(derived[name][g], 0.0) - emitted[B.apsim_name(g)][0]) < 1e-3
                   for g in genotypes), f"{name}: FigM06 y-intercepts are not the parameters"
        print(f"  {name:<14}matches Genotype_APSIM_Parameters_{name}.csv "
              f"(max deviation {worst:.1e})")

    for name in ARMS:
        dead = [B.apsim_name(g) for g in genotypes
                if ttn_at(derived[name][g], SD_LO) < 0 and ttn_at(derived[name][g], SD_HI) < 0]
        print(f"  {name:<14}clamped at TTN=0 across the realised range: "
              f"{len(dead)}  {sorted(dead)}")


def main() -> None:
    derived, genotypes, sd_min, sd_max = derive()
    print(f"genotypes: {len(genotypes)}   fitted S/D {sd_min:.2f}..{sd_max:.2f}   "
          f"drawn over the realised {SD_LO:.2f}..{SD_HI:.2f}")
    check(derived, genotypes)
    write_audit(derived, genotypes, sd_min, sd_max)

    colours = colour_map()
    figure_arabidopsis(derived, genotypes, colours, sd_min, sd_max,
                       OUT / "FigM04_Arabidopsis_SupplyDemand.png")
    figure_sorghum(derived, genotypes, colours,
                   OUT / "FigM05_Sorghum_Tillering_SupplyDemand.png")
    # Same data, same colours, same helpers -- only the layout differs, so the
    # combined figure cannot disagree with the two it is built from.
    figure_combined(derived, genotypes, colours, sd_min, sd_max,
                    OUT / "FigM04_05_SupplyDemand_Combined.png")
    # The same figure with the axis carried back to S/D 0 and nothing clamped, so
    # the intercept and slope of Table 1 can be read straight off it. See the
    # from_zero note above the panel helpers for why both are kept.
    figure_combined(derived, genotypes, colours, sd_min, sd_max,
                    OUT / "FigM06_SupplyDemand_FromZero.png", from_zero=True)
    print(f"\nfonts are true size: axis {FS_AXIS:.0f} pt, ticks {FS_TICK:.0f} pt, "
          f"legend {FS_LEGEND:.0f} pt on a {FIG_W} in figure")


if __name__ == "__main__":
    main()
