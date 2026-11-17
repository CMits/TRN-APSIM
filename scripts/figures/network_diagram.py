"""The three networks, side by side.

  04_figures/methods/FigM02_Networks_SideBySide.(png|svg)
      a  Hand-curated (32n/56e)   b  Flash-P (30n/44e)   c  Flash-P (69n/120e)

All three panels are drawn in the Flash-P house style and share one legend, one
node-size scale and one set of horizontal tier bands, so equivalent biology sits
at the same height in every panel and the eye can travel straight across.

WHY THIS FIGURE
FigM01 plots each genotype's sucrose response per network, and the things it
raises are all topological:

  * hand-curated resolves only 12 distinct lines out of 32 -- eleven genotypes
    coincide exactly because they all act through one SL -> MAX2 -> SMXL6,7,8 ->
    BRC1_2 chain. Distinct lines rise 12 -> 16 -> 22 with network size.
  * sucrose amplitude falls the other way, 4.31x -> 2.54x -> 1.26x. Not because
    of node count -- mean in-degree is ~2.2 in all three -- but because of fan-in
    at the phenotype node: Sustained_growth has 4 regulators, Tiller_Number 4,
    Shoot_Branching 10. Under PSoup's mean(activators) * 2/(1+mean(inhibitors))
    rule, sucrose is one voice in ten instead of one in four.

Both are visible here directly: count the sucrose routes (highlighted), and count
the arrows entering the phenotype node.

PROVENANCE
Adapted from FLASH-P_Paper/ComparisonShootBranching/scripts/12_network_diagram.py.
Vendored rather than imported: this folder reads nothing from a sibling project
(README.md:8-9), as with 00_reference/flashp_source and scripts/figures.

Changes from the original, and nothing else:
  * one panel per registered arm instead of two, each with its own xlim (the
    networks are 32 / 30 / 69 nodes; one shared x-scale would shrink the two
    small panels to illegibility). ylim IS shared, which keeps the bands aligned.
  * the arm registry (scripts/arms.py) replaces lib/registry + lib/netio;
    a local reader covers both on-disk formats.
  * the shared-concept ring means "present in every panel drawn" rather than
    "in both"; the legend says which in words, since over two arms it is a
    different and larger set than over three,
    resolved through 01_audit/ortholog_map.json plus a local alias table.
  * the sucrose node, its outgoing edges and the phenotype node are emphasised.

Everything worth changing lives in the EDIT ME block. Node positions are seeded
once and then frozen to scripts/figures/network_layout.json; edit two numbers
there and re-run to nudge a node. Pass --relayout to re-seed from scratch.
"""
from __future__ import annotations

import csv
import json
import math
import string
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# Keep SVG labels as real <text> so the figure stays editable in
# Illustrator / Inkscape, and use the house Arial/Helvetica stack.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import PathPatch, Polygon  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from matplotlib.transforms import Affine2D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS, FIG_ROOT, RESTRICTED  # noqa: E402

STYLE_PATH = HERE / "style" / "flashp_style.json"
LAYOUT_PATH = HERE / "network_layout.json"
ORTHOLOGS = ROOT / "01_audit" / "ortholog_map.json"
OUT_DIR = FIG_ROOT / "methods"
OVERLAP_CSV = ROOT / "01_audit" / "network_concept_overlap.csv"

# =========================================================================== #
# ===============================  EDIT ME  ================================= #
# =========================================================================== #

# Built at the size it will occupy on a LANDSCAPE page, so Word does not rescale
# it and every point size below is the printed point size. Building wide and letting
# Word shrink it is what made the earlier 32 in version illegible: squeezed into a
# 9.7 in page it divided every font by 3.3, turning 19 pt node labels into 5.8 pt.
#
# Keep FIG_H large enough that the figure stays WIDTH-limited. place_axes() solves
# scale = min(h_avail/yspan, w_avail/sum_spans); if height binds, the panels shrink
# horizontally and the labels lose the room the stagger just bought them.
FIG_W, FIG_H = 10.1, 7.4      # inches -- A4 landscape, ~0.8 in side margins
DPI = 300
WORD_WIDTH_IN = 10.1          # placed at its own width, i.e. no rescaling at all.

# True printed point sizes. What makes 9 pt possible at this width is the THREE-
# level label stagger (see LABEL_STAGGER): a label only has to clear the neighbour
# two columns away rather than the one next door, which triples its horizontal
# allowance. The widest tier across the three panels is 35 node columns, so at
# 9.7 in each column is ~0.22 in and a 3-level stagger gives each label ~0.66 in --
# enough for a ten-character name at 9 pt (~0.62 in). X_GAP cannot buy this room:
# widening it grows the data span and shrinks the scale by the same factor.
FONT_NODE = 8.0               # node labels
FONT_TIER = 9.0               # tier band names in the left margin
FONT_TITLE = 11.0             # panel titles
FONT_LEGEND = 8.5
FONT_TAG = 13.0               # a / b / c panel tags

X_GAP = 3.00                  # data units between adjacent nodes in a tier
Y_GAP = 11.60                 # data units between tiers
R_BASE = 0.64                 # node radius (GENE / default) in data units
BORDER_LW = 1.07              # node border, points
RING_LW = 0.95                # "shared concept" ring, points
RING_SCALE = 1.24             # ring size relative to the node

EDGE_LW = 0.66                # edge line, points
EDGE_ALPHA = 0.9
HEAD_SCALE = 1.6              # arrowhead / inhibition-bar size multiplier
EDGE_PAD = 0.08               # gap between node border and edge end, data units

# --- sucrose / phenotype emphasis ---------------------------------------- #
# The question this figure exists to answer, so it is drawn, not annotated.
KEY_BORDER_LW = 2.21          # border on the sucrose and phenotype nodes
KEY_RING_COLOR = "#111111"

# --- wiring -------------------------------------------------------------- #
ROUTING = "orthogonal"        # or "curved"
FAN_SPREAD = 90               # % of the node face edges fan across
LANE_ZONE = 7.40              # top of the routing lanes, below a tier's labels
LANE_STEP = 0.27              # horizontal lane pitch (auto-compressed)
LANE_PAD = 0.08               # x-clearance before two runs may share a lane
LANE_DEPTH = 1.80             # depth available for lanes below LANE_ZONE
SIDE_PAD = 0.60               # outermost node -> first side corridor
SIDE_STEP = 0.22              # side-corridor lane pitch
PEEL_SPAN = 5                 # bands spanned before an edge peels to the side
CORR_SUB = (0.0, 0.20, -0.20)  # sub-lanes within one free corridor
CORNER_R = 0.15               # bend fillet radius, data units
CLEARANCE = 0.10              # extra gap an edge must keep from a foreign node

EDGE_CURVE = 0.10             # ROUTING="curved" only
SAME_TIER_CURVE = 0.45        # ROUTING="curved" only

LABEL_DROP = 0.28             # label offset below the node, data units
LABEL_STAGGER = 3.40          # drop per stagger LEVEL (see LABEL_LEVELS)
LABEL_LEVELS = 3              # labels cycle over this many drops across a tier.
                              # 3 rather than 2 is what makes a legible font size
                              # possible at page width: a label then only has to
                              # clear the node two columns away, not its immediate
                              # neighbour, so it may be ~3 columns wide instead of 2.

PANEL_LEFT, PANEL_RIGHT = 0.085, 0.995
PANEL_BOTTOM, PANEL_TOP = 0.150, 0.930
PANEL_GAP = 0.020             # gutter between panels, figure fraction

# Nudge individual nodes after the auto-layout: ("flashp_shoot", "BRC1_2"): (x, y)
POS_OVERRIDE: dict = {}

# The shared tier grammar. This is what makes the panels comparable: the same
# eight bands, in the same order, in all three. Every node must appear exactly
# once in its network's column -- the script refuses to draw otherwise.
#
# hand_curated is lifted verbatim from the reference figure's `psoup` column;
# our 02_networks/Hand_Curated is byte-identical to that network, so its frozen
# coordinates transfer too and panel a reproduces the published one exactly.
TIERS = [
    ("Environmental input", {
        "hand_curated": ["Decap_signal", "Low_R_FR", "High_R_FR"],
        "flashp_suc": ["Decapitation_Signal", "Low_R_FR", "High_R_FR"],
        "flashp_shoot": ["Decap_signal", "Low_R_FR", "High_R_FR",
                         "Light", "Nitrate", "Phosphate"],
    }),
    ("Mobile signal", {
        "hand_curated": ["Aux_Shoot", "SUC", "Aux_Bud"],
        "flashp_suc": ["Auxin_Shoot", "Sucrose", "Auxin_Bud"],
        "flashp_shoot": ["Aux_Shoot", "SUC", "Aux_Bud", "miR156"],
    }),
    ("Hormone biosynthesis", {
        "hand_curated": ["D27", "MAX3", "MAX4", "MAX1", "CLAMT", "LBO", "CXE"],
        "flashp_suc": ["D27", "MAX3", "MAX4", "MAX1", "CLAMT", "LBO", "CXE"],
        # The shoot network carries both the Arabidopsis and the grass-orthologue
        # names for the same SL biosynthesis steps (MAX3=D17=CCD7, MAX4=D10=CCD8),
        # so both are placed here, adjacent, rather than one being hidden.
        "flashp_shoot": ["D27", "MAX3", "D17", "MAX4", "D10", "MAX1", "CLAMT",
                         "LBO", "CXE", "IPT", "LOG", "CKX", "SD1"],
    }),
    ("Hormone pool", {
        "hand_curated": ["SL", "CK", "GA"],
        "flashp_suc": ["Strigolactone", "Cytokinin", "Gibberellin"],
        "flashp_shoot": ["SL", "CK", "GA"],
    }),
    ("Perception & signalling", {
        "hand_curated": ["PIFs", "PhyB", "MAX2", "D14", "Perception_SL",
                         "SMXL6_7_8", "SMXL6_7_8_Transcription", "FHY3_FAR1"],
        # The empty slot opposite PSoup's Perception_SL is deliberate: this
        # Flash-P build has no perception node, it wires SL -> D14 -| SMXL678.
        "flashp_suc": ["PIFS", "PHYB", "MAX2", "D14", "SMXL678", "FHY3FAR1"],
        "flashp_shoot": ["PIFS", "PHYB", "MAX2", "D3", "SBMAX2", "D14",
                         "Perception_SL", "SMXL6_7_8", "SMXL6_7_8_Transcription",
                         "D53", "FHY3_FAR1", "GID1", "DELLA", "TIR1", "AUXIAA",
                         "ARF", "AHK3", "AHP", "ARR_A", "ARR_B"],
    }),
    ("Transcriptional integrator", {
        "hand_curated": ["SPL9_15", "BRC1_2", "HB21_40_53"],
        "flashp_suc": ["SPL9_15", "BRC1", "GT1"],
        "flashp_shoot": ["SPL9_15", "IPA1", "BRC1_2", "TB1", "GT1",
                         "HB21_40_53", "MOC1", "NGR5", "DEP1",
                         "TCP15", "TIE1", "TIN1"],
    }),
    ("Transport & downstream", {
        "hand_curated": ["PAT_Bud", "Bud_release", "NCED3", "ABA"],
        "flashp_suc": ["PIN1", "Bud_Release", "NCED3", "ABA"],
        "flashp_shoot": ["PIN1", "PIN3", "DW3", "PAT_Bud", "Bud_release",
                         "Axillary_Meristem_Initiation", "Axillary_Bud_Outgrowth",
                         "Sustained_growth", "NCED3", "ABA"],
    }),
    ("Phenotype", {
        "hand_curated": ["Sustained_growth"],
        "flashp_suc": ["Tiller_Number"],
        "flashp_shoot": ["Shoot_Branching"],
    }),
]

# Tiers whose left-to-right order above is authoritative: the auto-layout leaves
# them alone. Every tier is pinned as shipped -- drop a name to hand that tier
# back to the barycentre sweep.
PIN_ORDER = {
    "Environmental input",
    "Mobile signal",
    "Hormone biosynthesis",
    "Hormone pool",
    "Perception & signalling",
    "Transcriptional integrator",
    "Transport & downstream",
    "Phenotype",
}

# Long ids broken for legibility. Anything not listed is drawn as-is.
DISPLAY = {
    "Low_R_FR": "Low R:FR",
    "High_R_FR": "High R:FR",
    "Shoot_Branching": "Shoot\nbranching",
    "Sustained_growth": "Sustained\ngrowth",
    "Tiller_Number": "Tiller\nnumber",
    "Decap_signal": "Decap.",
    "Decapitation_Signal": "Decap.",
    "Aux_Shoot": "Aux\n(shoot)",
    "Aux_Bud": "Aux\n(bud)",
    "Auxin_Shoot": "Aux\n(shoot)",
    "Auxin_Bud": "Aux\n(bud)",
    "PAT_Bud": "PAT\n(bud)",
    "Bud_release": "Bud\nrelease",
    "Bud_Release": "Bud\nrelease",
    "Perception_SL": "Percep.\nSL",
    "SMXL6_7_8": "SMXL6/7/8",
    "SMXL6_7_8_Transcription": "SMXL\nmRNA",
    "SMXL678": "SMXL6/7/8",
    "FHY3_FAR1": "FHY3/FAR1",
    "FHY3FAR1": "FHY3/FAR1",
    "HB21_40_53": "HB21/40/53",
    "BRC1_2": "BRC1/2",
    "SPL9_15": "SPL9/15",
    "Strigolactone": "SL",
    "Cytokinin": "CK",
    "Gibberellin": "GA",
    "Axillary_Bud_Outgrowth": "Ax. bud\noutgrowth",
    "Axillary_Meristem_Initiation": "Ax.\nmeristem",
    "Sustained_growth": "Sustained\ngrowth",
}

# Panel order, straight from the registry -- so a restricted ARMS_ONLY set draws
# only those networks instead of crashing on a name that is no longer loaded.
# Titles come from the registry too, so the node/edge counts printed beside a
# label are counted from the network and cannot drift.
PANELS = list(ARMS)

# Regression gate: what each network is expected to contain. Not derived from
# the files being drawn, so a network changing underneath this figure is caught.
EXPECT = {"hand_curated": (32, 56), "flashp_suc": (30, 44), "flashp_shoot": (69, 120)}


def _all_phrase(n: int) -> str:
    """How to say "the shared-concept ring covers every panel", for n panels.

    The ring is an intersection over the networks actually drawn, so over two it
    is a genuinely different -- larger -- set than over three. The label has to
    say which, or the figure quietly claims more than it shows.
    """
    return {2: "in both", 3: "in all three"}.get(n, f"in all {n}")

# Concepts that ortholog_map.json does not cover: hormone / process / phenotype
# spellings that differ between the benchmark convention and the psoup one.
# Keyed by canonical name -> the id used in each network that has it.
EXTRA_ALIASES = {
    "SUCROSE": {"SUC", "Sucrose"},
    "STRIGOLACTONE": {"SL", "Strigolactone"},
    "CYTOKININ": {"CK", "Cytokinin"},
    "GIBBERELLIN": {"GA", "Gibberellin"},
    "ABA": {"ABA"},
    "AUXIN_BUD": {"Aux_Bud", "Auxin_Bud"},
    "AUXIN_SHOOT": {"Aux_Shoot", "Auxin_Shoot"},
    "BUD_RELEASE": {"Bud_release", "Bud_Release"},
    "DECAPITATION": {"Decap_signal", "Decapitation_Signal"},
    "LOW_R_FR": {"Low_R_FR"},
    "HIGH_R_FR": {"High_R_FR"},
    "PHYB": {"PHYB", "PhyB"},
    "PIFS": {"PIFS", "PIFs"},
    "PHENOTYPE": {"Sustained_growth", "Tiller_Number", "Shoot_Branching"},
}

# =========================================================================== #
# ============================  END OF EDIT ME  ============================= #
# =========================================================================== #

TYPE_LABELS = {
    "GENE": "Gene",
    "HORMONE": "Hormone",
    "METABOLITE": "Metabolite",
    "REGULATORY_RNA": "Regulatory RNA",
    "ENVIRONMENT": "Environment",
    "PHENOTYPE": "Phenotype",
    "PROCESS": "Process",
    "PROTEIN_COMPLEX": "Protein complex",
}
LEGEND_TYPE_ORDER = ["GENE", "HORMONE", "METABOLITE", "PROTEIN_COMPLEX",
                     "REGULATORY_RNA", "PROCESS", "ENVIRONMENT", "PHENOTYPE"]

# Flash-P stores node type as a single letter; the hand-curated CSV already uses
# the style file's own vocabulary.
TY_LETTER = {"G": "GENE", "H": "HORMONE", "M": "METABOLITE", "E": "ENVIRONMENT",
             "P": "PHENOTYPE", "PR": "PROCESS", "PC": "PROTEIN_COMPLEX",
             "R": "REGULATORY_RNA"}

BG = "#fbfdfc"
INK = "#13251c"
MUTED = "#5b6b63"
BAND = "#eef3f0"
RING_SHARED = "#13251c"
RING_UNIQUE = "#c2cec8"


def pick(m, key):
    """Resolve a {default, by} discrete mapping -- the flashp-graph.js rule."""
    by = m.get("by") or {}
    return by[key] if by.get(key) is not None else m["default"]


# --------------------------------------------------------------------------- #
# reading the networks
#
# Two on-disk formats, sniffed the way lib/netio.py does it in the source
# project: a FlashP export carries network/network.json, the hand-curated PSoup
# model carries Edges.csv beside PSoup_Network_nodes.csv.
# --------------------------------------------------------------------------- #

def load_network(arm):
    d = arm.network_dir
    nj = d / "network" / "network.json"
    if nj.is_file():
        net = json.loads(nj.read_text(encoding="utf-8"))
        nodes = {n["id"]: {"type": TY_LETTER[n["ty"]], "label": n.get("fn", n["id"])}
                 for n in net["nodes"]}
        edges = [(e["s"], e["t"], 1 if e["x"] > 0 else -1) for e in net["edges"]]
        return nodes, edges

    rows = list(csv.DictReader((d / "PSoup_Network_nodes.csv").open(encoding="utf-8-sig")))
    nodes = {r["node_id"].strip(): {"type": r["type"].strip(),
                                    "label": r.get("label", "").strip()}
             for r in rows}
    edges = []
    for r in csv.DictReader((d / "Edges.csv").open(encoding="utf-8-sig")):
        sign = -1 if r["Influence"].strip().lower().startswith("inhib") else 1
        edges.append((r["From"].strip(), r["To"].strip(), sign))
    return nodes, edges


def canonical_map(nets):
    """{network: {node: CANONICAL}} across all three, from ortholog_map + aliases.

    ortholog_map.json already carries the benchmark and psoup_style spelling of
    16 gene concepts (BRC1<->BRC1_2, SMXL678<->SMXL6_7_8, PIN1<->PIN_1_3_4_7,
    GT1<->HB21_40_53, ...). EXTRA_ALIASES covers the hormone / process / phenotype
    names that map does not reach. Anything still unmatched falls back to its own
    id, so a node is never silently merged with an unrelated one.
    """
    om = json.loads(ORTHOLOGS.read_text(encoding="utf-8"))
    by_name: dict[str, str] = {}
    for concept, v in om.items():
        for key in ("benchmark", "psoup_style", "sorghum", "arabidopsis"):
            name = v.get(key)
            if name:
                by_name.setdefault(name, concept)
    for concept, names in EXTRA_ALIASES.items():
        for name in names:
            by_name[name] = concept

    out = {}
    for net_id, (nodes, _e) in nets.items():
        out[net_id] = {n: by_name.get(n, n) for n in nodes}
    return out


# --------------------------------------------------------------------------- #
# unit shapes (centred on the origin, nominal radius 1)
# --------------------------------------------------------------------------- #

def _rounded_rect(hw, hh, rad):
    k = rad * 0.5523
    v = [(-hw + rad, -hh), (hw - rad, -hh),
         (hw - rad + k, -hh), (hw, -hh + rad - k), (hw, -hh + rad),
         (hw, hh - rad),
         (hw, hh - rad + k), (hw - rad + k, hh), (hw - rad, hh),
         (-hw + rad, hh),
         (-hw + rad - k, hh), (-hw, hh - rad + k), (-hw, hh - rad),
         (-hw, -hh + rad),
         (-hw, -hh + rad - k), (-hw + rad - k, -hh), (-hw + rad, -hh)]
    c = [MPath.MOVETO, MPath.LINETO,
         MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
         MPath.LINETO,
         MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
         MPath.LINETO,
         MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
         MPath.LINETO,
         MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    return MPath(v, c)


SHAPES = {
    "ellipse": MPath.circle(center=(0.0, 0.0), radius=1.0),
    "diamond": MPath([(0, 1.15), (1.15, 0), (0, -1.15), (-1.15, 0), (0, 1.15)],
                     [MPath.MOVETO] + [MPath.LINETO] * 3 + [MPath.CLOSEPOLY]),
    "triangle": MPath([(0, 1.2), (1.04, -0.6), (-1.04, -0.6), (0, 1.2)],
                      [MPath.MOVETO] + [MPath.LINETO] * 2 + [MPath.CLOSEPOLY]),
    # Cytoscape 'vee': a downward chevron -- reads as "something coming in".
    "vee": MPath([(-1.15, 1.05), (0, 0.15), (1.15, 1.05), (0, -1.15), (-1.15, 1.05)],
                 [MPath.MOVETO] + [MPath.LINETO] * 3 + [MPath.CLOSEPOLY]),
    "round-rectangle": _rounded_rect(1.75, 0.92, 0.34),
}


def shape_of(style, ntype):
    return SHAPES[pick(style["node"]["shape"], ntype)]


def radius_of(style, ntype):
    """Node radius in data units, keeping the house 65 / 80 / 30 size ratios."""
    base = float(style["node"]["size"]["default"])
    return R_BASE * float(pick(style["node"]["size"], ntype)) / base


def unit_patch(style, ntype, xy, scale, **kw):
    r = radius_of(style, ntype) * scale
    p = shape_of(style, ntype).transformed(Affine2D().scale(r).translate(*xy))
    return PathPatch(p, **kw)


# --------------------------------------------------------------------------- #
# layout
# --------------------------------------------------------------------------- #

def tier_of(net_id):
    out = {}
    for i, (_, members) in enumerate(TIERS):
        for n in members[net_id]:
            if n in out:
                raise SystemExit(f"{net_id}: '{n}' listed in more than one tier")
            out[n] = i
    return out


def seed_positions(net_id, nodes, edges):
    """Barycentre sweep: y from the tier, x ordered to reduce crossings."""
    tier = tier_of(net_id)
    missing = sorted(set(nodes) - set(tier))
    extra = sorted(set(tier) - set(nodes))
    if missing or extra:
        raise SystemExit(
            f"{net_id}: TIERS is out of step with the network.\n"
            f"  not placed in any tier: {missing}\n"
            f"  in TIERS but not in the network: {extra}")

    rows = [list(members[net_id]) for _, members in TIERS]
    nb = {n: [] for n in nodes}
    for s, t, _ in edges:
        if tier[s] != tier[t]:
            nb[s].append(t)
            nb[t].append(s)

    def slots():
        # Integer offset, so every node x is a multiple of X_GAP and the odd
        # multiples of X_GAP/2 are guaranteed-free vertical corridors.
        out = {}
        for row in rows:
            off = -((len(row) - 1) // 2)
            for i, n in enumerate(row):
                out[n] = (off + i) * X_GAP
        return out

    free = [i for i, (name, _) in enumerate(TIERS) if name not in PIN_ORDER]
    for it in range(24):
        pos = slots()
        order = free if it % 2 == 0 else list(reversed(free))
        for ti in order:
            row = rows[ti]
            key = {}
            for n in row:
                near = [pos[m] for m in nb[n]]
                key[n] = sum(near) / len(near) if near else pos[n]
            rows[ti] = sorted(row, key=lambda n: (key[n], n))

    pos = slots()
    return rows, {n: (pos[n], -tier[n] * Y_GAP) for n in nodes}


def load_layout(nets, relayout):
    """Seed once, then honour scripts/figures/network_layout.json on later runs."""
    saved = {}
    if LAYOUT_PATH.is_file() and not relayout:
        saved = json.loads(LAYOUT_PATH.read_text(encoding="utf-8"))

    out, rows, changed = {}, {}, False
    for net_id, (nodes, edges) in nets.items():
        seeded_rows, seeded = seed_positions(net_id, nodes, edges)
        rows[net_id] = seeded_rows
        have = saved.get(net_id, {})
        if set(have) == set(nodes):
            out[net_id] = {n: tuple(have[n]) for n in nodes}
        else:
            out[net_id] = seeded
            changed = True

    for (net_id, n), xy in POS_OVERRIDE.items():
        out[net_id][n] = tuple(xy)

    if changed:
        LAYOUT_PATH.write_text(
            json.dumps({k: {n: [round(x, 4), round(y, 4)] for n, (x, y) in v.items()}
                        for k, v in out.items()}, indent=2) + "\n", encoding="utf-8")
        print("  wrote", LAYOUT_PATH.name,
              "(edit it to nudge a node; --relayout to re-seed)")
    else:
        print("  using node positions from", LAYOUT_PATH.name)
    return rows, out


# --------------------------------------------------------------------------- #
# edge routing -- axis-aligned, reproducing ELK edgeRouting=ORTHOGONAL on a
# fixed layout (ELK itself would recompute the layers and destroy the bands).
# --------------------------------------------------------------------------- #

def half_extent(style, ntype):
    """Half-width / half-height of the node's bounding box, in data units."""
    r = radius_of(style, ntype)
    if pick(style["node"]["shape"], ntype) == "round-rectangle":
        return 1.75 * r, 0.92 * r
    return r, r


def _pack(items):
    """Greedy interval colouring. items = [(lo, hi, key)] -> {key: lane}, n_lanes."""
    lanes, out = [], {}
    for lo, hi, key in sorted(items, key=lambda t: (t[0], t[1])):
        for k, used in enumerate(lanes):
            if all(hi < a - LANE_PAD or lo > b + LANE_PAD for a, b in used):
                used.append((lo, hi))
                out[key] = k
                break
        else:
            lanes.append([(lo, hi)])
            out[key] = len(lanes) - 1
    return out, len(lanes)


def _fillet(pts, r):
    """Round the corners of an axis-aligned polyline."""
    if len(pts) < 3:
        return list(pts)
    out = [pts[0]]
    for i in range(1, len(pts) - 1):
        a, b, c = pts[i - 1], pts[i], pts[i + 1]
        la = math.hypot(b[0] - a[0], b[1] - a[1])
        lc = math.hypot(c[0] - b[0], c[1] - b[1])
        rr = min(r, la / 2.0, lc / 2.0)
        if rr < 1e-6:
            out.append(b)
            continue
        ua = ((a[0] - b[0]) / la, (a[1] - b[1]) / la)
        uc = ((c[0] - b[0]) / lc, (c[1] - b[1]) / lc)
        p, q = (b[0] + ua[0] * rr, b[1] + ua[1] * rr), (b[0] + uc[0] * rr, b[1] + uc[1] * rr)
        out.append(p)
        for j in range(1, 6):                    # quadratic bezier through the bend
            t = j / 6.0
            out.append(((1 - t) ** 2 * p[0] + 2 * (1 - t) * t * b[0] + t * t * q[0],
                        (1 - t) ** 2 * p[1] + 2 * (1 - t) * t * b[1] + t * t * q[1]))
        out.append(q)
    out.append(pts[-1])
    return out


def lane_y(gap, k, n):
    """y of horizontal lane k of n in the routing zone below tier `gap`."""
    step = min(LANE_STEP, LANE_DEPTH / max(1, n - 1)) if n > 1 else LANE_STEP
    return -gap * Y_GAP - LANE_ZONE - k * step


def route_edges(style, net_id, nodes, edges, pos):
    """Fixed node positions -> one axis-aligned polyline per edge."""
    tier = tier_of(net_id)
    ext = {n: half_extent(style, nodes[n]["type"]) for n in nodes}
    xs = [pos[n][0] for n in nodes]
    x_lo, x_hi, x_mid = min(xs), max(xs), (min(xs) + max(xs)) / 2.0

    recs = []
    for i, (s, t, sign) in enumerate(edges):
        r = {"i": i, "s": s, "t": t, "sign": sign, "span": tier[t] - tier[s]}
        r["kind"] = ("self" if s == t else
                     "flat" if r["span"] == 0 else
                     "direct" if abs(r["span"]) == 1 else
                     "peel" if abs(r["span"]) >= PEEL_SPAN else "corridor")
        recs.append(r)

    # --- 1. fan the endpoints across the node faces ------------------------ #
    groups = {}
    for r in recs:
        if r["kind"] == "self":
            continue
        groups.setdefault(("src", r["s"], r["span"] >= 0), []).append(r)
        groups.setdefault(("tgt", r["t"], r["span"] > 0), []).append(r)
    for (role, n, _), grp in groups.items():
        other = "t" if role == "src" else "s"
        grp.sort(key=lambda r: (pos[r[other]][0], r["i"]))
        hw = ext[n][0]
        for j, r in enumerate(grp):
            off = 0.0 if len(grp) == 1 else (j / (len(grp) - 1.0) - 0.5) * (FAN_SPREAD / 100.0) * 2 * hw
            r["off_src" if role == "src" else "off_tgt"] = off

    for r in recs:
        if r["kind"] == "self":
            continue
        (sx, sy), (tx, ty) = pos[r["s"]], pos[r["t"]]
        sh, th = ext[r["s"]][1], ext[r["t"]][1]
        down_out = r["span"] >= 0            # leaves the source's bottom face
        down_in = r["span"] > 0              # enters the target's top face
        r["p0"] = (sx + r.get("off_src", 0.0),
                   sy - sh - EDGE_PAD if down_out else sy + sh + EDGE_PAD)
        r["p1"] = (tx + r.get("off_tgt", 0.0),
                   ty + th + EDGE_PAD if down_in else ty - th - EDGE_PAD)
        ts, tt = tier[r["s"]], tier[r["t"]]
        if r["span"] >= 0:
            r["g1"] = ts
            r["g2"] = tt - 1 if r["span"] > 1 else ts
        else:
            r["g1"] = ts - 1
            r["g2"] = tt

    # --- 2. drop the multi-band edges through free in-field corridors ------ #
    slots_ = sorted(
        {round((2 * k + 1) * X_GAP / 2.0 + d, 4)
         for k in range(int((x_lo - X_GAP) / X_GAP) - 1,
                        int((x_hi + X_GAP) / X_GAP) + 1)
         for d in CORR_SUB})
    busy = {c: [] for c in slots_}
    for r in sorted((r for r in recs if r["kind"] == "corridor"),
                    key=lambda r: (abs(r["span"]), r["i"])):
        want = r["p0"][0]
        lo, hi = min(r["g1"], r["g2"]), max(r["g1"], r["g2"])
        free = [c for c in slots_ if all(hi < a or lo > b for a, b in busy[c])]
        if not free:                         # everything taken -> use the margin
            r["kind"] = "peel"
            continue
        c = min(free, key=lambda c: (abs(c - want), c))
        busy[c].append((lo, hi))
        r["mid_x"] = c

    # --- 3. side corridors for the long haul ------------------------------- #
    sides = {"L": [], "R": []}
    for r in recs:
        if r["kind"] != "peel":
            continue
        r["side"] = "L" if (r["p0"][0] + r["p1"][0]) / 2.0 < x_mid else "R"
        sides[r["side"]].append((min(r["g1"], r["g2"]), max(r["g1"], r["g2"]), r["i"]))
    side_lane, n_side = {}, {}
    for k, items in sides.items():
        side_lane[k], n_side[k] = _pack(items)
    for r in recs:
        if r["kind"] != "peel":
            continue
        k = side_lane[r["side"]][r["i"]]
        r["mid_x"] = (x_lo - SIDE_PAD - k * SIDE_STEP if r["side"] == "L"
                      else x_hi + SIDE_PAD + k * SIDE_STEP)

    # --- 4. horizontal lanes, packed per gap ------------------------------- #
    runs = {}
    for r in recs:
        if r["kind"] == "self":
            continue
        if r["kind"] in ("flat", "direct"):
            a, b = sorted((r["p0"][0], r["p1"][0]))
            runs.setdefault(r["g1"], []).append((a, b, (r["i"], 1)))
        else:
            a, b = sorted((r["p0"][0], r["mid_x"]))
            runs.setdefault(r["g1"], []).append((a, b, (r["i"], 1)))
            a, b = sorted((r["mid_x"], r["p1"][0]))
            runs.setdefault(r["g2"], []).append((a, b, (r["i"], 2)))
    hlane, n_lane = {}, {}
    for g, items in runs.items():
        got, n = _pack(items)
        hlane.update({(g,) + k: v for k, v in got.items()})
        n_lane[g] = n

    # --- 5. build the polylines -------------------------------------------- #
    out = []
    for r in recs:
        if r["kind"] == "self":                       # D14 -| D14
            (x, y) = pos[r["s"]]
            hw, hh = ext[r["s"]]
            w, h = hw + 0.42, hh + 0.34
            pts = [(x + hw + EDGE_PAD, y), (x + w, y), (x + w, y + h),
                   (x, y + h), (x, y + hh + EDGE_PAD)]
            out.append((_fillet(pts, CORNER_R), r["sign"]))
            continue
        L1 = lane_y(r["g1"], hlane[(r["g1"], r["i"], 1)], n_lane[r["g1"]])
        if r["kind"] in ("flat", "direct"):
            pts = [r["p0"], (r["p0"][0], L1), (r["p1"][0], L1), r["p1"]]
        else:
            L2 = lane_y(r["g2"], hlane[(r["g2"], r["i"], 2)], n_lane[r["g2"]])
            pts = [r["p0"], (r["p0"][0], L1), (r["mid_x"], L1),
                   (r["mid_x"], L2), (r["p1"][0], L2), r["p1"]]
        # Drop a zero-length jog so a straight drop stays straight.
        clean = [pts[0]]
        for p in pts[1:]:
            if math.hypot(p[0] - clean[-1][0], p[1] - clean[-1][1]) > 1e-9:
                clean.append(p)
        out.append((_fillet(clean, CORNER_R), r["sign"]))

    stats = {"lanes": n_lane, "side": n_side,
             "kinds": {k: sum(1 for r in recs if r["kind"] == k)
                       for k in ("self", "flat", "direct", "corridor", "peel")}}
    return out, stats


def route_curved(style, net_id, nodes, edges, pos):
    """ROUTING = "curved": centre-to-centre beziers."""
    tier = tier_of(net_id)
    out = []
    for s, t, sign in edges:
        p0, p2 = pos[s], pos[t]
        r0 = radius_of(style, nodes[s]["type"])
        if s == t:
            cx, cy = p0[0] + r0 * 1.25, p0[1] + r0 * 1.25
            pts = [(cx + r0 * 1.05 * math.cos(a), cy + r0 * 1.05 * math.sin(a))
                   for a in [(-70 + i * 2.2) * math.pi / 180.0 for i in range(120)]]
        else:
            dx, dy = p2[0] - p0[0], p2[1] - p0[1]
            rad = SAME_TIER_CURVE if tier[s] == tier[t] else EDGE_CURVE
            mx, my = (p0[0] + p2[0]) / 2.0, (p0[1] + p2[1]) / 2.0
            cx, cy = mx - dy * rad, my + dx * rad
            pts = [((1 - u) ** 2 * p0[0] + 2 * (1 - u) * u * cx + u * u * p2[0],
                    (1 - u) ** 2 * p0[1] + 2 * (1 - u) * u * cy + u * u * p2[1])
                   for u in (i / 149.0 for i in range(150))]
        r2 = radius_of(style, nodes[t]["type"])
        i = next((k for k, q in enumerate(pts)
                  if math.hypot(q[0] - p0[0], q[1] - p0[1]) >= r0 + EDGE_PAD), 0)
        j = next((k for k in range(len(pts) - 1, -1, -1)
                  if math.hypot(pts[k][0] - p2[0], pts[k][1] - p2[1]) >= r2 + EDGE_PAD),
                 len(pts) - 1)
        out.append((pts[i:j + 1] if j - i >= 8 else pts[len(pts) // 2 - 4:len(pts) // 2 + 5],
                    sign))
    return out, {"lanes": {}, "side": {}, "kinds": {}}


def foreign_hits(style, net_id, nodes, pos, edges, routes):
    """How many routed edges pass over a node that is not one of their endpoints."""
    ext = {n: half_extent(style, nodes[n]["type"]) for n in nodes}
    hits = []
    for (pts, _), (s, t, _sg) in zip(routes, edges):
        for n, (x, y) in pos.items():
            if n in (s, t):
                continue
            hw, hh = ext[n]
            hw, hh = hw + CLEARANCE, hh + CLEARANCE
            if any(abs(px - x) <= hw and abs(py - y) <= hh for px, py in pts):
                hits.append((s, t, n))
                break
    return hits


def draw_route(ax, style, pts, sign):
    key = "negative" if sign < 0 else "positive"
    alpha = EDGE_ALPHA
    ax.add_line(Line2D([q[0] for q in pts], [q[1] for q in pts],
                       color=pick(style["edge"]["line"], key), lw=EDGE_LW,
                       alpha=alpha, solid_capstyle="round",
                       solid_joinstyle="round", zorder=1))

    # Head geometry ported from flashp-prop-core.js::headD -- filled triangle for
    # activation, perpendicular bar for inhibition.
    head_c = pick(style["edge"]["arrowColor"], key)
    tip = pts[-1]
    back = pts[max(0, len(pts) - 4)]
    m = math.hypot(tip[0] - back[0], tip[1] - back[1]) or 1.0
    ca, sa = (tip[0] - back[0]) / m, (tip[1] - back[1]) / m
    u = R_BASE * HEAD_SCALE
    if sign > 0:
        length, half = 0.50 * u, 0.27 * u
        tri = [(tip[0], tip[1]),
               (tip[0] - length * ca + half * sa, tip[1] - length * sa - half * ca),
               (tip[0] - length * ca - half * sa, tip[1] - length * sa + half * ca)]
        ax.add_patch(Polygon(tri, closed=True, facecolor=head_c, edgecolor="none",
                             alpha=alpha, zorder=1.5))
    else:
        bar = 0.39 * u
        ax.add_line(Line2D([tip[0] + bar * sa, tip[0] - bar * sa],
                           [tip[1] - bar * ca, tip[1] + bar * ca],
                           color=head_c, lw=EDGE_LW + 1.2 * (10.1 / 32.0), alpha=alpha,
                           solid_capstyle="butt", zorder=1.5))


# --------------------------------------------------------------------------- #
# panel
# --------------------------------------------------------------------------- #

def draw_panel(ax, style, net_id, nodes, edges, rows, pos, shared, routes,
               xlim, ylim, key_nodes):
    labels = []
    ax.set_facecolor(BG)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")

    for i in range(len(TIERS)):
        if i % 2 == 0:
            y = -i * Y_GAP
            ax.axhspan(y - Y_GAP / 2, y + Y_GAP / 2, color=BAND, zorder=0, lw=0)

    for pts, sign in routes:
        draw_route(ax, style, pts, sign)

    # Stagger labels by left-to-right rank within the tier, so neighbouring
    # labels never sit on the same baseline.
    rank = {n: i for row in rows for i, n in enumerate(row)}
    for n, meta in nodes.items():
        ty = meta["type"]
        xy = pos[n]
        is_shared = n in shared
        is_key = n in key_nodes.values()
        ax.add_patch(unit_patch(style, ty, xy, RING_SCALE, facecolor="none",
                                edgecolor=RING_SHARED if is_shared else RING_UNIQUE,
                                lw=RING_LW if is_shared else 1.0,
                                linestyle="solid" if is_shared else (0, (2, 2)),
                                zorder=2))
        ax.add_patch(unit_patch(style, ty, xy, 1.0,
                                facecolor=pick(style["node"]["fill"], ty),
                                edgecolor=KEY_RING_COLOR if is_key
                                else pick(style["node"]["border"], ty),
                                lw=KEY_BORDER_LW if is_key else BORDER_LW,
                                zorder=3.5 if is_key else 3))
        drop = LABEL_DROP + LABEL_STAGGER * (rank[n] % LABEL_LEVELS)
        # The house renderer puts a 0.85-opacity #fbfdfc plate behind the label.
        # A real bbox reproduces that and leaves the SVG label editable <text>.
        t = ax.text(xy[0], xy[1] - radius_of(style, ty) * RING_SCALE - drop,
                    DISPLAY.get(n, n), ha="center", va="top", fontsize=FONT_NODE,
                    fontweight="bold" if is_key else "semibold",
                    color=INK, zorder=4, linespacing=0.95,
                    bbox=dict(boxstyle="round,pad=0.10", facecolor=BG,
                              edgecolor="none", alpha=0.85))
        labels.append((n, t))

    return labels


def fit_labels(fig, axes, panel_labels, xlims, rounds=4):
    """Widen each panel's x limits until no node label runs off its canvas.

    A label is centred on its node and is routinely wider than a column, so the
    outermost ones overhang the node field by an amount that depends on the very
    scale it is setting. Measuring and re-measuring converges in a few passes.
    Per panel, not shared: the three networks are 32 / 30 / 69 nodes wide.
    """
    for _ in range(rounds):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        done = True
        for ax in axes:
            net_id = ax._net_id
            lo, hi = xlims[net_id]
            inv = ax.transData.inverted()
            need = 0.0
            for _n, t in panel_labels[net_id]:
                bb = t.get_window_extent(renderer=r)
                x0, _ = inv.transform((bb.x0, bb.y0))
                x1, _ = inv.transform((bb.x1, bb.y0))
                need = max(need, lo - x0, x1 - hi)
            if need > 0.01:
                xlims[net_id] = (lo - need * 1.05, hi + need * 1.05)
                ax.set_xlim(*xlims[net_id])
                done = False
        if done:
            break
    return xlims


def place_axes(fig, axes, xlims, ylim):
    """Size every panel from ONE inches-per-data-unit, then place them.

    The three networks span different widths, so each panel needs its own xlim.
    But `aspect="equal"` then fights any width we invent: matplotlib shrinks the
    axes box until the data aspect fits, which silently pushes the tier bands to
    different heights in different panels -- exactly the alignment the figure
    exists to provide. So the scale is solved for once and every rectangle is
    derived from it: one scale means the bands line up, node glyphs are the same
    physical size everywhere, and equal aspect is satisfied without shrinking.
    """
    n = len(axes)
    w_avail = (PANEL_RIGHT - PANEL_LEFT) * FIG_W - PANEL_GAP * FIG_W * (n - 1)
    h_avail = (PANEL_TOP - PANEL_BOTTOM) * FIG_H
    yspan = ylim[1] - ylim[0]
    spans = [xlims[ax._net_id][1] - xlims[ax._net_id][0] for ax in axes]
    scale = min(h_avail / yspan, w_avail / sum(spans))       # inches per data unit

    h_frac = (yspan * scale) / FIG_H
    widths = [(s * scale) / FIG_W for s in spans]
    # Centre the row of panels, and hang them from PANEL_TOP so the titles align.
    total = sum(widths) + PANEL_GAP * (n - 1)
    x = PANEL_LEFT + ((PANEL_RIGHT - PANEL_LEFT) - total) / 2.0
    bottom = PANEL_TOP - h_frac
    for ax, w in zip(axes, widths):
        ax.set_position([x, bottom, w, h_frac])
        x += w + PANEL_GAP
    return scale


def label_collisions(fig, panel_labels):
    """Pairwise overlap of the drawn label plates, in device pixels."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    out = {}
    for net_id, labels in panel_labels.items():
        boxes = [(n, t.get_window_extent(renderer=r)) for n, t in labels]
        hits = []
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                (a, ba), (b, bb) = boxes[i], boxes[j]
                if ba.overlaps(bb):
                    hits.append((a, b))
        out[net_id] = hits
    return out


def write_overlap_csv(nets, canon, present_in_all):
    """Per-node canonical concept and whether every network carries it."""
    OVERLAP_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OVERLAP_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["arm", "node", "node_type", "canonical_concept", "in_all_arms"])
        for net_id in PANELS:
            nodes, _ = nets[net_id]
            for n in sorted(nodes):
                c = canon[net_id][n]
                w.writerow([net_id, n, nodes[n]["type"], c, c in present_in_all])
    print(f"  wrote {OVERLAP_CSV.relative_to(ROOT)}")


# --------------------------------------------------------------------------- #

def main():
    relayout = "--relayout" in sys.argv
    style = json.loads(STYLE_PATH.read_text(encoding="utf-8"))

    nets = {net_id: load_network(ARMS[net_id]) for net_id in PANELS}
    for net_id, (nodes, edges) in nets.items():
        got = (len(nodes), len(edges))
        if got != EXPECT[net_id]:
            raise SystemExit(f"{net_id}: expected {EXPECT[net_id]} nodes/edges, got {got}")

    # Shared concepts across every drawn network, resolved through the ortholog map.
    canon = canonical_map(nets)
    concept_sets = {net_id: set(canon[net_id].values()) for net_id in nets}
    present_in_all = set.intersection(*concept_sets.values())
    shared = {net_id: {n for n in nets[net_id][0]
                       if canon[net_id][n] in present_in_all} for net_id in nets}
    print(f"  concepts per network: "
          f"{', '.join(f'{k} {len(v)}' for k, v in concept_sets.items())}")
    print(f"  present {_all_phrase(len(PANELS))}: {len(present_in_all)}")
    write_overlap_csv(nets, canon, present_in_all)

    key_nodes = {net_id: {"sucrose": ARMS[net_id].sucrose_node,
                          "phenotype": ARMS[net_id].phenotype} for net_id in PANELS}

    rows, pos = load_layout(nets, relayout)

    router = route_edges if ROUTING == "orthogonal" else route_curved
    routes, stats = {}, {}
    for net_id, (nodes, edges) in nets.items():
        routes[net_id], stats[net_id] = router(style, net_id, nodes, edges, pos[net_id])
        hits = foreign_hits(style, net_id, nodes, pos[net_id], edges, routes[net_id])
        st = stats[net_id]
        n_suc = sum(1 for s, _t, _x in edges if s == key_nodes[net_id]["sucrose"])
        n_phen = sum(1 for _s, t, _x in edges if t == key_nodes[net_id]["phenotype"])
        print(f"  {net_id}: {len(edges)} edges routed "
              f"({', '.join(f'{k} {v}' for k, v in st['kinds'].items() if v)})"
              if st["kinds"] else f"  {net_id}: {len(edges)} edges")
        print(f"      sucrose '{key_nodes[net_id]['sucrose']}' out-edges: {n_suc}"
              f"   |   '{key_nodes[net_id]['phenotype']}' in-edges (fan-in): {n_phen}")
        if st["lanes"]:
            print(f"      lanes per band gap: "
                  f"{', '.join(str(st['lanes'][g]) for g in sorted(st['lanes']))}"
                  f" | side corridors L{st['side'].get('L', 0)} R{st['side'].get('R', 0)}")
        print(f"      edges crossing a foreign node: {len(hits)}"
              + ("" if not hits else "  -> " + ", ".join(
                  f"{s}->{t} over {n}" for s, t, n in hits[:6])
                 + (" ..." if len(hits) > 6 else "")))

    # ylim is SHARED so the tier bands line up across panels -- that alignment is
    # the whole point. xlim is PER PANEL: the networks are 32 / 30 / 69 nodes and
    # a common x-scale would shrink the two small panels to illegibility.
    ys = ([p[1] for v in pos.values() for p in v.values()]
          + [q[1] for v in routes.values() for pts, _ in v for q in pts])
    ylim = (min(min(ys), -(len(TIERS) - 1) * Y_GAP) - Y_GAP * 0.30,
            max(ys) + Y_GAP * 0.28)
    pad = X_GAP * 0.18
    xlims = {}
    for net_id in PANELS:
        xs = ([p[0] for p in pos[net_id].values()]
              + [q[0] for pts, _ in routes[net_id] for q in pts])
        xlims[net_id] = (min(xs) - pad, max(xs) + pad)

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI, facecolor=BG)
    h = PANEL_TOP - PANEL_BOTTOM
    axes, panel_labels = [], {}
    # Provisional rectangles -- place_axes() solves the real geometry once the
    # labels have been measured and the x limits have stopped moving.
    x = PANEL_LEFT
    prov = (PANEL_RIGHT - PANEL_LEFT - PANEL_GAP * (len(PANELS) - 1)) / len(PANELS)
    for tag, net_id in zip(string.ascii_lowercase, PANELS):
        arm = ARMS[net_id]
        ax = fig.add_axes([x, PANEL_BOTTOM, prov, h])
        x += prov + PANEL_GAP
        nodes, edges = nets[net_id]
        ax._net_id = net_id
        panel_labels[net_id] = draw_panel(
            ax, style, net_id, nodes, edges, rows[net_id], pos[net_id],
            shared[net_id], routes[net_id], xlims[net_id], ylim, key_nodes[net_id])
        ax.set_title(f"{arm.label}\n{len(nodes)} nodes, {len(edges)} edges",
                     fontsize=FONT_TITLE, fontweight="semibold", color=INK,
                     pad=14, linespacing=1.5)
        ax.text(0.0, 1.055, tag, transform=ax.transAxes, fontsize=FONT_TAG,
                fontweight="bold", color=INK, ha="left", va="bottom")
        axes.append(ax)

    for i, (name, _) in enumerate(TIERS):
        axes[0].text(-0.035, -i * Y_GAP, textwrap.fill(name, 13, break_long_words=False),
                     transform=axes[0].get_yaxis_transform(),
                     ha="right", va="center", fontsize=FONT_TIER, color=MUTED,
                     fontweight="semibold", linespacing=1.25)

    # ---- one legend for all three panels ----
    present = [t for t in LEGEND_TYPE_ORDER
               if any(m["type"] == t for nodes, _ in nets.values() for m in nodes.values())]
    handles = [Line2D([], [], linestyle="none",
                      marker=shape_of(style, t), markersize=15,
                      markerfacecolor=pick(style["node"]["fill"], t),
                      markeredgecolor=pick(style["node"]["border"], t),
                      markeredgewidth=1.6, label=TYPE_LABELS.get(t, t))
               for t in present]
    handles += [
        Line2D([], [], color=pick(style["edge"]["line"], "positive"), lw=3.0,
               marker=">", markersize=9, markevery=[-1],
               markerfacecolor=pick(style["edge"]["arrowColor"], "positive"),
               markeredgecolor="none", label="activates  (→)"),
        Line2D([], [], color=pick(style["edge"]["line"], "negative"), lw=3.0,
               marker="|", markersize=13, markeredgewidth=3.0, markevery=[-1],
               markeredgecolor=pick(style["edge"]["arrowColor"], "negative"),
               label="inhibits  (—|)"),
        Line2D([], [], linestyle="none", marker="o", markersize=15,
               markerfacecolor="none", markeredgecolor=KEY_RING_COLOR,
               markeredgewidth=3.0, label="sucrose / phenotype node"),
        Line2D([], [], linestyle="none", marker="o", markersize=15,
               markerfacecolor="none", markeredgecolor=RING_SHARED, markeredgewidth=2.0,
               label=f"concept {_all_phrase(len(PANELS))} ({len(present_in_all)})"),
        Line2D([], [], linestyle="none", marker="o", markersize=15,
               markerfacecolor="none", markeredgecolor=RING_UNIQUE, markeredgewidth=1.4,
               label=f"not {_all_phrase(len(PANELS))}"),
    ]
    leg = fig.legend(handles=handles, loc="lower center", ncol=5,
                     bbox_to_anchor=(0.5, 0.004), frameon=True, fontsize=FONT_LEGEND,
                     handletextpad=0.7, columnspacing=1.9, labelspacing=1.0,
                     borderpad=1.0)
    leg.get_frame().set_facecolor("#ffffff")
    leg.get_frame().set_edgecolor("#d7e0db")
    leg.get_frame().set_linewidth(1.0)

    missing_legend = {m["type"] for nodes, _ in nets.values()
                      for m in nodes.values()} - set(present)
    if missing_legend:
        raise SystemExit(f"node types drawn but absent from the legend: {sorted(missing_legend)}")

    # Labels overhang the node field by an amount that depends on the scale they
    # are setting, so measure, widen, then solve the geometry -- and repeat, since
    # a new scale changes what overhangs. Two passes is enough in practice.
    for _ in range(2):
        xlims = fit_labels(fig, axes, panel_labels, xlims)
        meas = place_axes(fig, axes, xlims, ylim)

    for net_id, hits in label_collisions(fig, panel_labels).items():
        print(f"  {net_id}: overlapping node labels: {len(hits)}"
              + ("" if not hits else "  -> " + ", ".join(
                  f"{a}/{b}" for a, b in hits[:8])
                 + (" ..." if len(hits) > 8 else "")))
    print(f"  scale: {meas:.3f} in per data unit")
    print(f"  node labels are {FONT_NODE:.0f} pt on a {FIG_W:.1f} in figure "
          f"-> {FONT_NODE * WORD_WIDTH_IN / FIG_W:.1f} pt when placed "
          f"{WORD_WIDTH_IN:.1f} in wide in Word")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = "FigM02_Networks_SideBySide"
    # FIG_W is chosen for the full registry, so a restricted arm set centres a
    # shorter row of panels in a canvas built for a longer one -- 14% dead margin
    # either side at two panels. Crop that back; the node glyphs and label point
    # sizes are untouched, the figure is simply placed narrower. Never taken at the
    # full registry, so FigM02 there is byte-for-byte what it always was.
    bbox = "tight" if RESTRICTED else None
    if bbox:
        print("  restricted arm set -- cropping the margin the shorter row leaves")
    for ext in ("png", "svg"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}", dpi=DPI, facecolor=BG,
                    bbox_inches=bbox)
    plt.close(fig)
    print(f"  Saved: {OUT_DIR / (stem + '.png')}")
    print(f"  Saved: {OUT_DIR / (stem + '.svg')}")


if __name__ == "__main__":
    main()
