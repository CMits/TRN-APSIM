"""Methods figure -- Step 1: the sucrose response of every genotype, in all three networks.

This is the figure of the FIRST step of the chain, before any of the downstream
arithmetic can hide a difference: each genotype propagated through each network at
sucrose 0.5 / 1.0 / 1.5, one line per genotype, one panel per network. It is the
three-network version of the paper's own Fig01_PSoup_Sucrose_Response.png.

WHICH NUMBERS
All three arms are simulated with PSoup's update rule via PSoupFormSimulator,
by calling build_apsim_parameters.network_response() / psoup_response() -- the very
functions that feed steps 3-5 and the APSIM parameters. That matters twice over:

  * the equation FORM is held constant across the arms, so a difference between
    panels is attributable to network topology and nothing else;
  * FlashP's own algebraic engine is avoided. Its clamped inhibition
    min(1/max(prod(inh), 0.1), 10) saturates, collapsing 22 of 37 genotypes onto one
    curve and leaving 25/111 combinations unconverged on the shoot network
    (see psoup_form_sim.py). The 01_audit/*_sucrose_response.csv sweeps ARE that
    engine; they are the solver diagnostic, not this figure's input.

So this figure will NOT match 01_audit/flashp*_sucrose_response.csv, and is not
meant to. It does match 01_audit/handcurated_sucrose_response.csv, because the
hand-curated arm is read from that R sweep either way.

GENOTYPES
The 32 the reproduction actually carries end to end -- the paper's 34 minus the two
FHY3_FAR1 lines, which have no confirmed grass orthologue. Same `eligible` set
build_apsim_parameters.py builds, derived here the same way rather than hardcoded.

Reads  02_networks/*/network/algebraic_equations.json (via the arm registry)
       01_audit/handcurated_sucrose_response.csv
Writes 04_figures/methods/FigM01_Sucrose_Response_AllNetworks.{png,svg}
       01_audit/step1_sucrose_response_psoupform.csv
       01_audit/sucrose_wiring.csv
"""

from __future__ import annotations

import csv
import json
import re
import string
import sys
from collections import defaultdict, deque
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "00_reference" / "psoup_source"))

import build_apsim_parameters as B  # noqa: E402
from arms import ARMS, FIG_ROOT     # noqa: E402

PD = ROOT / "00_reference" / "paper_data"
AUDIT = ROOT / "01_audit"
OUT = FIG_ROOT / "methods"

LEVELS = ["low", "medium", "high"]
LEVEL_LABELS = ["Low", "Medium", "High"]

# Prettified phenotype names for the y-axes. Each arm reports a differently-named
# node; the axis should read as English, not as a node id.
Y_LABEL = {
    "Sustained_growth": "Sustained growth",
    "Tiller_Number": "Tiller number",
    "Shoot_Branching": "Shoot branching",
}

# ---------------------------------------------------------------------------
# Colours and labels, carried over from the paper's own Fig 1 cell
# (NewScripts/Pipeline.ipynb) so a genotype keeps its identity across figures.
# Copied rather than imported: nothing in this folder reads from the parent project.
#
# paletteer Polychrome::palette36 minus the near-white entries, plus 4 Kelly colours.
# ---------------------------------------------------------------------------
POLYCHROME = [
    '#5A5156', '#F6222E', '#FE00FA', '#16FF32', '#3283FE', '#FEAF16', '#B00068',
    '#1CFFCE', '#90AD1C', '#2ED9FF', '#DEA0FD', '#AA0DFE', '#F8A19F', '#325A9B',
    '#C4451C', '#1C8356', '#85660D', '#B10DA1', '#FBE426', '#1CBE4F', '#FA0087',
    '#FC1CBF', '#F7E1A0', '#C075A6', '#782AB6', '#AAF400', '#BDCDFF', '#822E1C',
    '#B5EFB5', '#7ED7D1', '#1C7F93', '#D85FF7', '#683B79', '#66B0FF', '#3B00FB',
    '#FFB300', '#803E75', '#A6BDD7', '#CEA262',
]

# Sizes are TRUE POINT SIZES on the page. The figure is built at the width it will
# occupy in the thesis (FIG_W below), so Word does not rescale it and 10 pt here is
# 10 pt printed. Building wide and letting Word shrink it was what made the earlier
# version illegible: a 21 in figure squeezed into a 6.5 in column divides every
# font by 3.2, turning a 13 pt legend into 4 pt.
FIG_W = 6.5                   # inches -- portrait, full text width

FS_AXIS = 10.0
FS_TICK = 9.0
FS_TITLE = 10.0
FS_TAG = 12.0
FS_LEGEND = 8.0

# Height is built from three fixed INCH measurements rather than being a constant
# with fractional margins hung off it, so a figure over two arms comes out shorter
# instead of stretching two panels across three rows' worth of page. These three
# reproduce the original 9.0 in / top=0.845 / bottom=0.055 exactly at three arms.
LEGEND_IN = 1.395             # band above the panels -- the 8-row genotype legend
BOTTOM_IN = 0.495             # room for the x tick labels and the x-label
ROW_IN = 2.37                 # one panel row


def stack_geometry(nrow: int) -> tuple[float, float, float]:
    """(figure height, gridspec top, gridspec bottom) for an nrow stacked figure.

    Shared with methods_figure_sd.py so FigM01, FigM04 and FigM05 stay the same
    shape as each other whatever the arm set is.
    """
    h = LEGEND_IN + BOTTOM_IN + ROW_IN * nrow
    return h, 1.0 - LEGEND_IN / h, BOTTOM_IN / h


def pretty_name(name: str) -> str:
    """PSoup genotype -> paper label: 'BRC1_2_n' -> 'brc1,2'. (Pipeline.ipynb)"""
    s = re.sub(r'_n(?=$|[,_\s])', '', name)
    s = s.replace('_', ',')
    return s.lower()


def eligible_genotypes() -> list[str]:
    """The 32 genotypes carried through every arm, built the way build_apsim_parameters does."""
    lin = list(csv.DictReader(
        (PD / "Genotype_Linear_Equations.csv").open(encoding="utf-8-sig")))
    return [r["Genotype"] for r in lin if B.EXCLUDE_TOKEN not in r["Genotype"]]


def published_genotypes() -> list[str]:
    """All 34 published genotypes -- the colour index is taken from THIS list.

    Indexing colours by rank within the 34 rather than the plotted 32 is what keeps
    every genotype on the same colour as in Figures/Fig01_PSoup_Sucrose_Response.png.
    """
    lin = list(csv.DictReader(
        (PD / "Genotype_Linear_Equations.csv").open(encoding="utf-8-sig")))
    return [r["Genotype"] for r in lin]


def colour_map() -> dict[str, str]:
    ordered = sorted(published_genotypes(), key=pretty_name)
    return {g: POLYCHROME[i % len(POLYCHROME)] for i, g in enumerate(ordered)}


def step1_response(a) -> dict[str, dict[str, float]]:
    """{genotype: {level: value}} for one arm, from the pipeline's own Step-1 code."""
    if a.kind == "r_psoup":
        return B.psoup_response(a.response_csv, a.phenotype)
    return B.network_response(a.equations, a.phenotype, a.sucrose_node, a.gene_map)


# ---------------------------------------------------------------------------
# Sucrose wiring diagnostic
# ---------------------------------------------------------------------------
def sucrose_reach(eq_path: Path, source: str) -> set[str]:
    """Every node reachable from the sucrose node, over activator + inhibitor lists."""
    spec = json.loads(Path(eq_path).read_text(encoding="utf-8"))
    fwd: dict[str, set[str]] = defaultdict(set)
    for e in spec["equations"]:
        for u in (e.get("a") or []) + (e.get("inh") or []):
            fwd[u].add(e["n"])
    seen, q = {source}, deque([source])
    while q:
        for y in fwd[q.popleft()]:
            if y not in seen:
                seen.add(y)
                q.append(y)
    return seen


def sucrose_edges(a) -> list[tuple[str, int]]:
    """(target, sign) for every edge leaving the arm's sucrose node."""
    if a.equations is None:                       # hand-curated: two plain CSVs
        rows = csv.DictReader(
            (a.network_dir / "Edges.csv").open(encoding="utf-8-sig"))
        return [(r["To"].strip(), -1 if r["Influence"].strip() == "inhibition" else 1)
                for r in rows if r["From"].strip() == a.sucrose_node]
    net = json.loads((a.network_dir / "network" / "network.json").read_text(encoding="utf-8"))
    return [(e["t"], e["x"]) for e in net["edges"] if e["s"] == a.sucrose_node]


def write_audit(responses: dict, genotypes: list[str]) -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)

    tidy = AUDIT / "step1_sucrose_response_psoupform.csv"
    with tidy.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["arm", "genotype", "apsim_name", "sucrose_label",
                    "SUC_n", "value", "value_rel_medium"])
        for name, a in ARMS.items():
            resp = responses[name]
            for g in genotypes:
                v = resp[g]
                for lv, suc in zip(LEVELS, (0.5, 1.0, 1.5)):
                    w.writerow([name, g, B.apsim_name(g), lv, suc,
                                f"{v[lv]:.9f}", f"{v[lv]/v['medium']:.9f}"])
    print(f"Saved: {tidy}")

    wiring = AUDIT / "sucrose_wiring.csv"
    with wiring.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["arm", "sucrose_node", "phenotype", "out_edges", "targets",
                    "reaches_phenotype", "nodes_downstream",
                    "wt_low", "wt_medium", "wt_high", "amplitude", "distinct_shapes"])
        for name, a in ARMS.items():
            resp = responses[name]
            edges = sucrose_edges(a)
            targets = ";".join(f"{'-|' if x < 0 else '->'}{t}" for t, x in edges)
            if a.equations is None:
                # The curated model has no algebraic_equations.json to walk; its five
                # edges include SUC -> Sustained_growth, so reach is direct by inspection.
                reaches, ndown = True, ""
            else:
                r = sucrose_reach(a.equations, a.sucrose_node)
                reaches, ndown = a.phenotype in r, len(r) - 1
            shapes = {(round(resp[g]["low"] / resp[g]["medium"], 6),
                       round(resp[g]["high"] / resp[g]["medium"], 6)) for g in genotypes}
            wt = resp["WT"]
            w.writerow([name, a.sucrose_node, a.phenotype, len(edges), targets,
                        reaches, ndown,
                        f"{wt['low']:.6f}", f"{wt['medium']:.6f}", f"{wt['high']:.6f}",
                        f"{wt['high']/wt['low']:.3f}", len(shapes)])
    print(f"Saved: {wiring}")


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def make_figure(responses: dict, genotypes: list[str], png: Path) -> None:
    plt.rcParams.update({"font.size": FS_TICK,
                         "axes.labelsize": FS_AXIS,
                         "xtick.labelsize": FS_TICK,
                         "ytick.labelsize": FS_TICK})
    colours = colour_map()
    # Legend order is alphabetical by the printed label, as in the paper's Fig 1.
    ordered = sorted(genotypes, key=pretty_name)

    # Panels STACKED, not side by side. At 6.5 in three panels in a row leave under
    # 2 in each; stacked, every panel gets the full column width and they share the
    # Low/Medium/High axis, which only the bottom panel needs to label.
    # One row per registered arm, so a restricted ARMS_ONLY set gives a shorter
    # figure rather than a blank trailing panel.
    nrow = len(ARMS)
    fig_h, top, bottom = stack_geometry(nrow)
    fig = plt.figure(figsize=(FIG_W, fig_h))
    gs = fig.add_gridspec(nrow, 1, hspace=0.34, top=top, bottom=bottom,
                          left=0.135, right=0.985)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(nrow)]

    handles = None
    for i, (ax, (name, a), tag) in enumerate(
            zip(axes, ARMS.items(), string.ascii_lowercase)):
        resp = responses[name]
        hs = []
        for g in ordered:
            y = [resp[g][lv] for lv in LEVELS]
            (ln,) = ax.plot([0, 1, 2], y, color=colours[g], lw=1.0,
                            marker="o", ms=2.6, mew=0, label=pretty_name(g))
            hs.append(ln)
        handles = handles or hs

        ax.set_xticks([0, 1, 2])
        ax.set_xlim(-0.12, 2.12)
        last = i == len(axes) - 1
        if last:
            ax.set_xticklabels(LEVEL_LABELS)
            ax.set_xlabel("Sucrose level")
        else:
            ax.set_xticklabels([])
        ax.set_ylabel(Y_LABEL.get(a.phenotype, a.phenotype))
        # Name and topology, from the registry -- arms.py counts the nodes and edges
        # from the network itself, so this label cannot drift from what was simulated.
        ax.set_title(a.fig_label, fontsize=FS_TITLE, pad=5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=FS_TICK)
        ax.text(-0.135, 1.02, tag, transform=ax.transAxes, fontsize=FS_TAG,
                ha="left", va="bottom")

    # Four columns puts the 32 genotypes on eight rows. Long labels such as
    # "brc1,2, pin,1,3,4,7" need about 1.35 in each at 8 pt, so four is the most
    # that fits the 6.5 in width without truncation.
    leg = fig.legend(handles=handles, loc="upper center", ncol=4,
                     fontsize=FS_LEGEND, frameon=False, handlelength=1.4,
                     columnspacing=1.0, labelspacing=0.32, handletextpad=0.4,
                     borderaxespad=0.0, bbox_to_anchor=(0.5, 1.0))
    for t in leg.get_texts():
        t.set_fontstyle("italic")

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(png.with_suffix(".svg"), format="svg", transparent=True,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {png.with_suffix('.svg')}")


def report(responses: dict, genotypes: list[str]) -> None:
    print(f"\n{'arm':<15}{'SUC edges':>10}{'WT low':>9}{'WT high':>9}"
          f"{'amp':>7}{'shapes':>8}   spread at high")
    for name, a in ARMS.items():
        resp = responses[name]
        wt = resp["WT"]
        shapes = {(round(resp[g]["low"] / resp[g]["medium"], 6),
                   round(resp[g]["high"] / resp[g]["medium"], 6)) for g in genotypes}
        highs = [resp[g]["high"] for g in genotypes]
        print(f"{name:<15}{len(sucrose_edges(a)):>10}{wt['low']:>9.3f}{wt['high']:>9.3f}"
              f"{wt['high']/wt['low']:>6.2f}x{len(shapes):>8}"
              f"   [{min(highs):.2f}, {max(highs):.2f}]")


def main() -> None:
    genotypes = eligible_genotypes()
    responses = {name: step1_response(a) for name, a in ARMS.items()}

    missing = {name: [g for g in genotypes if g not in r]
               for name, r in responses.items()}
    for name, gs in missing.items():
        if gs:
            raise SystemExit(f"{name}: no Step-1 response for {gs}")

    print(f"genotypes: {len(genotypes)}   arms: {', '.join(ARMS)}")
    write_audit(responses, genotypes)
    make_figure(responses, genotypes,
                OUT / "FigM01_Sucrose_Response_AllNetworks.png")
    report(responses, genotypes)


if __name__ == "__main__":
    main()
