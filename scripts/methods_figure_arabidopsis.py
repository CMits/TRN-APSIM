"""Supplementary data figure: the raw Arabidopsis shoot-branch observations.

  04_figures/methods/FigM03_Arabidopsis_Branches.(png|svg)

One violin + box per genotype for every published observation the paper uses, and
a diamond where a genotype has only one. Ordered by the derived anchor, descending,
matching the paper's own Fig 2.

WHAT THIS IS FOR
This is a DATA figure, not a methods figure. The Step-2 normalisation it feeds is
explained in prose (README, and the manuscript methods) because two sentences do
it better than a panel: wild type ranges from 0 to 9 branches across the 48
experiments, so each observation is expressed as a fold-change against the WT of
its own experiment, those are averaged per genotype, and the result is multiplied
by the median WT (3.0) purely to restore branch-number units.

What the figure adds that prose cannot is the shape of the evidence base: it is
thin and very uneven. 17 of the 32 genotypes rest on a SINGLE published
observation, and the anchor for those is that one number rescaled by its
experiment's WT. Given that the shared anchor alone reproduces every arm's final
TTN ranking (rho = 0.996 / 0.997 / 0.992), that is worth showing.

Nothing derived is drawn on the panel. The anchor is not the same quantity as a
branch count -- it is a fold-change rescaled by the global WT median -- so a mark
at that height on this y-axis would invite a comparison that does not hold. The
anchor, the fold-change it came from, the replication, and the WT that set it are
all in 01_audit/arabidopsis_anchor_audit.csv, which this script also writes.

CLEANING
Reproduces the paper's own steps: Arabidopsis rows only; WT auto-filled where the
genotype is named only in the free-text Mutant column; environmental treatments
(+SL, +Low_R_FR, +High_R_FR) dropped so that only untreated plants define the
baseline. 225 Arabidopsis observations -> 163 used, 62 dropped.

WHAT THIS DOES *NOT* CHANGE
The pipeline still reads the frozen `weighted_arabidopsis_branches.csv`, and it
should: the anchor has to be byte-identical across arms or a difference in the
final parameters could come from re-deriving the meta-analysis rather than from
the network (build_apsim_parameters.py:15-23). The recomputation here CHECKS that
file -- 34/34 genotypes exactly -- it does not replace it. Nothing downstream
reads this script's output.

Reads  00_reference/paper_data/DatasetS2_Arabidopsis_branches.csv   (raw, cp1252)
       00_reference/paper_data/weighted_arabidopsis_branches.csv    (the frozen anchor)
       00_reference/paper_data/Genotype_Linear_Equations.csv        (the published 34)
Writes 04_figures/methods/FigM03_Arabidopsis_Branches.{png,svg}
       01_audit/arabidopsis_anchor_audit.csv
"""

from __future__ import annotations

import csv
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "00_reference" / "psoup_source"))

import build_apsim_parameters as B          # noqa: E402
from arms import FIG_ROOT                   # noqa: E402
from methods_figure_sucrose import (        # noqa: E402
    POLYCHROME, pretty_name, eligible_genotypes, published_genotypes,
)

PD = ROOT / "00_reference" / "paper_data"
AUDIT = ROOT / "01_audit"
OUT = FIG_ROOT / "methods"

RAW = PD / "DatasetS2_Arabidopsis_branches.csv"
FROZEN = PD / "weighted_arabidopsis_branches.csv"

# Divide-by-zero floor for experiments whose WT is ~0, from the paper's own
# meta-analysis. Without it a WT of 0 makes every ratio in that experiment
# infinite and the genotype's anchor explodes.
WT_MIN_DENOM = 0.5

WT_LIKE = {"wt", "wild type", "wild-type", "col-0", "col 0"}

# Sized for a figure that gets dropped into a Word page and scaled down to the
# text width (~6.5 in) from the 13 in it is drawn at. At that ~50% reduction these
# land back on the 7 / 18 / 14 pt sizes the panel used to carry AT FULL SIZE, i.e.
# roughly double what a reader saw on the page before. Line and marker weights
# below are scaled with them, otherwise the strokes vanish under the larger text.
FS_TICK = 13.0
FS_LABEL = 18.0
FS_N = 10.0
FS_LEGEND = 14.0



def load_observations() -> tuple[list[tuple[str, str, float]], dict[str, int]]:
    """Every Arabidopsis observation the paper uses: (genotype, experiment, branches).

    Reproduces the cleaning in the paper's own notebook cell, step for step, so
    this figure and the published Fig 2 are drawn from the same 163 rows.
    """
    # cp1252, NOT utf-8-sig as everywhere else in scripts/: the Publication column
    # contains "Gonzalez-Grandio" with a Latin-1 accented byte, and utf-8 raises
    # UnicodeDecodeError on it. The file is frozen ground truth, so it is read as
    # it is rather than being re-encoded.
    rows = list(csv.DictReader(RAW.open(encoding="cp1252")))
    arab = [r for r in rows if (r.get("Species") or "").strip() == "Arabidopsis"]

    kept: list[tuple[str, str, float]] = []
    dropped = defaultdict(int)
    for r in arab:
        geno = (r.get("Psoup Treatment") or "").strip()
        mutant = (r.get("Mutant") or "").strip().lower()

        # Some WT rows name the genotype only in the free-text Mutant column.
        if not geno and mutant in WT_LIKE:
            geno = "WT"

        if not geno or geno.lower() == "nan":
            dropped["blank genotype"] += 1
            continue
        # Environmental treatments (+SL, +Low_R_FR, +High_R_FR) are excluded: the
        # anchor is the UNTREATED / control-sucrose plant. The source spells these
        # inconsistently ("WT, +High_R_FR" vs "BRC1_2_n,+ High_R_FR"), so the test
        # is a bare "+" rather than a regex on the separator.
        if "+" in geno:
            dropped["environmental treatment"] += 1
            continue

        value = (r.get("Actual tiller number") or "").strip()
        if not value:
            dropped["no branch count"] += 1
            continue

        # Same genotype-name normalisation the rest of the pipeline uses
        # (MAX_4_n -> MAX4_n, comma spacing), reused rather than re-written.
        geno = B.normalise(geno)
        # Merges the Luo 2019 publication-string typo, so its rows share one
        # experiment key and therefore one WT denominator.
        exp = re.sub(r",(?!\s)", ", ", (r.get("Publication") or "").strip())
        kept.append((geno, exp, float(value)))

    print(f"raw rows {len(rows)}  |  Arabidopsis {len(arab)}  |  used {len(kept)}")
    print("dropped: " + ", ".join(f"{n} {k}" for k, n in sorted(dropped.items()))
          + f"  (total {sum(dropped.values())})")
    return kept, dropped


def recompute_anchor(obs):
    """The Step-2 derivation, returned in full so the figure can draw every stage.

    Returns (anchor, ratios, wt_exp, wt_raw, wt_ref):
      ratios[g]   the fold-changes -- each observation over the WT of its OWN experiment
      wt_exp      the FLOORED WT denominator per experiment (what the ratios divide by)
      wt_raw      the unfloored per-experiment WT (what panel a's inset shows)
      wt_ref      the global WT reference, median over experiments
      anchor[g]   mean(ratios[g]) * wt_ref

    This is the recipe build_apsim_parameters.py:21-23 describes in prose but does
    not run. Reproducing it here turns that parenthetical into a measured number.
    """
    wt_by_exp = defaultdict(list)
    for geno, exp, v in obs:
        if geno == "WT":
            wt_by_exp[exp].append(v)
    # One WT per experiment in this dataset, but mean() rather than [0] so a second
    # WT row would be collapsed rather than silently ignored.
    wt_raw = {e: sum(v) / len(v) for e, v in wt_by_exp.items()}
    wt_exp = {e: max(v, WT_MIN_DENOM) for e, v in wt_raw.items()}

    # Absolute scale: the median WT across experiments. Median, not mean, because a
    # few experiments report atypically high or low WT and the mean chases them.
    wt_ref = statistics.median(wt_raw.values())

    ratios = defaultdict(list)
    for geno, exp, v in obs:
        if exp in wt_exp:                  # an experiment with no WT cannot contribute
            ratios[geno].append(v / wt_exp[exp])

    # Arithmetic mean, not geometric. A ratio is a multiplicative effect size, so
    # pooling on the log scale would normally be the right call -- but one
    # observation here is exactly zero (Umehara 2008 Exp 1, WT = 0), which leaves
    # the geometric mean undefined. The ordering is insensitive to the choice
    # anyway (Spearman rho 0.982 between the two, largest shift 4 places), so this
    # is kept for robustness to zero counts rather than because it is the
    # meta-analysis convention for ratios -- it is not.
    anchor = {g: sum(r) / len(r) * wt_ref for g, r in ratios.items()}
    return anchor, dict(ratios), wt_exp, wt_raw, wt_ref


def frozen_anchor() -> dict[str, float]:
    """The `medium` column the pipeline actually uses -- Step 2's input."""
    return {r["genotype"]: float(r["medium"]) for r in
            csv.DictReader(FROZEN.open(encoding="utf-8-sig"))}


def colour_map_34() -> dict[str, str]:
    """Genotype -> colour, indexed within the PUBLISHED 34.

    Same rule as methods_figure_sucrose.colour_map(): taking the index from the 34
    rather than the plotted 32 is what keeps a genotype on the colour it has in
    FigM01 and in the paper's own Fig01 / Fig02.
    """
    ordered = sorted(published_genotypes(), key=pretty_name)
    return {g: POLYCHROME[i % len(POLYCHROME)] for i, g in enumerate(ordered)}


def write_audit(obs, anchor, ratios, wt_exp, frozen, in_fig):
    by_geno = defaultdict(list)
    exps = defaultdict(set)
    for g, e, v in obs:
        by_geno[g].append(v)
        exps[g].add(e)

    AUDIT.mkdir(parents=True, exist_ok=True)
    path = AUDIT / "arabidopsis_anchor_audit.csv"
    exact = both = 0
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["genotype", "apsim_name", "in_figure", "n_observations",
                    "n_experiments", "raw_min", "raw_median", "raw_max",
                    "wt_of_experiment", "mean_fold_change",
                    "recomputed_anchor", "published_anchor_medium",
                    "abs_diff", "matches_exactly"])
        # All 34, so the two FHY3_FAR1 lines stay on the record even though the
        # figure shows the 32 the reproduction carries end to end.
        for g in sorted(frozen, key=B.apsim_name):
            vals = by_geno.get(g, [])
            rs = ratios.get(g, [])
            rec = anchor.get(g)
            pub = frozen[g]
            if rec is not None:
                both += 1
                d = abs(rec - pub)
                ok = d < 1e-6
                exact += ok
            else:
                d, ok = None, None
            # Only meaningful for a genotype measured once: which experiment's WT
            # set its anchor outright. With several observations there is no single one.
            wt_one = f"{wt_exp[next(iter(exps[g]))]:.3f}" if len(vals) == 1 else ""
            w.writerow([
                g, B.apsim_name(g), g in in_fig, len(vals), len(exps.get(g, ())),
                f"{min(vals):.3f}" if vals else "",
                f"{statistics.median(vals):.3f}" if vals else "",
                f"{max(vals):.3f}" if vals else "",
                wt_one,
                f"{sum(rs)/len(rs):.6f}" if rs else "",
                f"{rec:.6f}" if rec is not None else "",
                f"{pub:.6f}",
                f"{d:.6f}" if d is not None else "",
                ok if ok is not None else "",
            ])
    print(f"Saved: {path}")
    print(f"recomputed anchor reproduces the frozen file exactly: {exact}/{both}")
    return exact, both


def _distributions(ax, order, pos, data, colours):
    """Violin + box where a genotype has several observations, diamond where it has one."""
    n_violin = n_diamond = 0
    for g in order:
        vals = data[g]
        col = colours.get(g, "#666666")
        x = pos[g]
        if len(vals) > 1:
            n_violin += 1
            vp = ax.violinplot([vals], positions=[x], widths=0.8,
                               showmeans=False, showmedians=False, showextrema=False)
            for body in vp["bodies"]:
                body.set_facecolor(col)
                body.set_edgecolor(col)
                body.set_alpha(0.25)
                body.set_linewidth(1.0)
            ax.boxplot([vals], positions=[x], widths=0.35, patch_artist=True,
                       showfliers=True,
                       medianprops=dict(color="black", linewidth=1.6),
                       boxprops=dict(facecolor=col, edgecolor=col, linewidth=1.2, alpha=0.9),
                       whiskerprops=dict(color=col, linewidth=1.2),
                       capprops=dict(color=col, linewidth=1.2),
                       flierprops=dict(marker="o", markersize=4,
                                       markerfacecolor=col, markeredgecolor=col))
        else:
            n_diamond += 1
            ax.scatter([x], [vals[0]], marker="D", s=110, facecolor=col,
                       edgecolor="black", linewidth=0.8, zorder=3)
    return n_violin, n_diamond


def make_figure(obs, colours, order, png: Path):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": FS_TICK,
        "axes.labelsize": FS_LABEL,
        "axes.linewidth": 1.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.major.size": 5.0,
        "ytick.major.size": 5.0,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
    })

    raw = defaultdict(list)
    for g, _e, v in obs:
        raw[g].append(v)

    pos = {g: i for i, g in enumerate(order)}
    # Narrower and taller than before: the 32 rotated genotype labels are now
    # ~2x the height, so they need the extra vertical room, and a less extreme
    # aspect ratio survives being fitted to a Word text column.
    fig, ax = plt.subplots(figsize=(13.0, 7.0))
    n_violin, n_diamond = _distributions(ax, order, pos, raw, colours)

    ax.set_xticks(list(pos.values()))
    ax.set_xticklabels([pretty_name(g) for g in order], rotation=45, ha="right",
                       fontstyle="italic", fontsize=FS_TICK)
    ax.set_ylabel("Shoot branch number")
    ax.set_xlim(-0.7, len(order) - 0.3)
    ax.margins(x=0.01)

    ax.legend(handles=[
        mpatches.Patch(facecolor="0.7", edgecolor="0.4", alpha=0.6,
                       label="Multiple observations (n > 1)"),
        Line2D([0], [0], marker="D", linestyle="", markerfacecolor="0.7",
               markeredgecolor="black", markersize=11,
               label="Single observation (n = 1)"),
    ], loc="upper right", fontsize=FS_LEGEND, frameon=True, framealpha=0.9,
        edgecolor="0.6")

    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")
    print(f"genotypes drawn: {len(order)}  ({n_violin} violins, {n_diamond} diamonds)")


def main() -> None:
    obs, _dropped = load_observations()
    anchor, ratios, wt_exp, wt_raw, wt_ref = recompute_anchor(obs)

    floored = sorted(e for e, v in wt_raw.items() if v < WT_MIN_DENOM)
    print(f"global WT reference (median WT): {wt_ref:g} branches "
          f"over {len(wt_raw)} experiments with a WT")
    print(f"WT across experiments: {min(wt_raw.values()):g} .. {max(wt_raw.values()):g}"
          f"   below the {WT_MIN_DENOM:g} floor: {len(floored)}")
    for e in floored:
        print(f"    WT = {wt_raw[e]:g}  ->  {WT_MIN_DENOM:g}   ({e})")
    wt_mean_r = sum(ratios["WT"]) / len(ratios["WT"])
    print(f"WT's own mean fold-change: {wt_mean_r:.4f} -- exactly 1.0 but for the "
          f"floor, so WT anchors at {wt_mean_r * wt_ref:.3f}, not {wt_ref:g}")

    frozen = frozen_anchor()
    plotted = eligible_genotypes()          # the 32 carried end to end, as FigM01
    missing = [g for g in plotted if g not in anchor]
    if missing:
        raise SystemExit(f"no raw Arabidopsis observations for: {missing}")

    # The claim panel b rests on: the mean mark, times wt_ref, IS the anchor.
    for g in plotted:
        m = sum(ratios[g]) / len(ratios[g])
        assert abs(m * wt_ref - anchor[g]) < 1e-12, f"{g}: mean x wt_ref != anchor"

    write_audit(obs, anchor, ratios, wt_exp, frozen, set(plotted))

    # Ordered by the derived anchor, descending -- as the paper's own Fig 2. The
    # order therefore comes from panel b's statistic, which is why panel a's boxes
    # do not descend tidily: that mismatch is exactly why both panels are shown.
    order = sorted(plotted, key=lambda g: -anchor[g])
    make_figure(obs, colour_map_34(), order,
                OUT / "FigM03_Arabidopsis_Branches.png")


if __name__ == "__main__":
    main()


