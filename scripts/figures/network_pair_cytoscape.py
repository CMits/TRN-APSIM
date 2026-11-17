"""The hand-curated and Flash-P networks in the Cytoscape idiom, on ONE skeleton.

  04_figures/methods/FigM02b_Networks_Cytoscape.(png|svg)
      left  Hand-curated (32n/56e)      right  Flash-P (30n/44e)
  04_figures/methods/FigM02b_Networks_Cytoscape_Legend.(png|svg)
      the key, as its own file

WHY THIS FIGURE, GIVEN FigM02 ALREADY EXISTS
FigM02 (scripts/figures/network_diagram.py) answers "do the TIERS line up": eight
labelled bands, small nodes, every name dropped BELOW its node on a three-level
stagger so 69 of them fit across a landscape page. It is a tier map.

This one answers a different question -- "are these the same network?" -- and answers
it in the visual language the reader is already using on these graphs: Cytoscape's.
Node name INSIDE the node, one shape and fill per node type, a sparse top-to-bottom
hierarchy rather than dense bands.

THE ONE IDEA
Both panels are drawn on a SINGLE coordinate skeleton keyed on CANONICAL CONCEPT,
not on node id. `SL` in panel a and `Strigolactone` in panel b are the same point on
the page; so are `BRC1_2`/`BRC1`, `SMXL6_7_8`/`SMXL678`, `HB21_40_53`/`GT1` and
`PAT_Bud`/`PIN1` -- while each panel keeps its OWN network's spelling, because the
spelling is what the reader has in their own files. Co-position is imposed, not
discovered: there is no per-panel coordinate that could drift (assert_registered()
proves it). What is NOT imposed is which slots are filled and which edges exist, and
that is the comparison. Every difference therefore draws itself:

  * a concept one network lacks leaves a dashed EMPTY slot in the other panel,
    labelled with the missing node's name in grey;
  * an edge only one network has is drawn at reduced opacity, so the 31 shared
    edges stand out as the common skeleton and the 25 / 13 panel-only ones as the
    delta.

Measured from the networks themselves: 32 concepts, of which 30 are in both; 68
union edges, of which 31 are in both.

THREE THINGS THE CAPTION MUST CARRY
1. Node positions are identical between panels BY CONSTRUCTION -- the panels are not
   two independent layouts that happen to agree.
2. Panel b is drawn from network/network.json, the file the pipeline actually runs.
   The Cytoscape exports beside it (network/cytoscape/network.graphml) are STALE for
   the _suc arms: scripts/extend_sucrose.py rebuilds network.json and the equations
   and copies network/cytoscape/* through untouched, so that GraphML is the 40-edge
   PRE-sucrose network. The four edges it lacks --
       Sucrose -> Cytokinin, Sucrose -| MAX2, Sucrose -> SMXL678, Sucrose -> Tiller_Number
   -- are exactly the ones that define the arm, so a figure built from that GraphML
   would omit the sucrose wiring the whole analysis rests on. This one shows them.
3. Two nodes have DIFFERENT TYPES in the two networks, so the same slot is a
   different colour left and right. SMXL6,7,8 is a GENE in PSoup and a
   PROTEIN_COMPLEX in Flash-P; bud release is a PHENOTYPE in PSoup and a PROCESS in
   Flash-P. Each panel is drawn with its own network's type. That is information,
   not an error, and a caption that does not say so invites the reader to read it
   as one.

PAT_Bud <-> PIN1
01_audit/ortholog_map.json separates these because their node TYPES differ (PSoup
models auxin transport as a PROCESS, Flash-P as the PIN1 gene), but PSoup's PAT_Bud
carries gene modifier PIN_1_3_4_7_n -- it IS the PIN slot. CONCEPT_ALIAS joins them,
which is what makes the node counts reconcile: 32 - 30 = the two nodes only PSoup
has (Perception_SL, SMXL6_7_8_Transcription). Merging them does create a two-cycle
(PSoup has Aux_Bud -> PAT_Bud, Flash-P has PIN1 -> Auxin_Bud), which is why
AUXIN_BUD -> PIN_1_3_4_7 appears in FEEDBACK below.

WHY THE LAYERS ARE SOLVED, NOT DECLARED
FigM02's eight-band grammar is the right frame for a tier map and the wrong one
here: it puts seven biosynthesis genes in one row when the biology is a six-step
chain. This figure runs a longest-path layering over the union of both edge sets,
then a tighten sweep that pulls every node down to just above its earliest
successor. That is what produces the Cytoscape tree: 19 layers, never more than 4
nodes wide, with D27 -> MAX3 -> MAX4 -> MAX1 -> CLAMT -> LBO -> SL -> D14 -> MAX2
running dead straight down the centre of BOTH panels. Without the tighten sweep the
raw longest path parks every source on L00 and the widest row is 6.

Reads  02_networks/{Hand_Curated,Tillering_Sorghum_Benchmark_suc}
       scripts/figures/style/flashp_style.json
       01_audit/ortholog_map.json
Writes <FIG_ROOT>/methods/FigM02b_Networks_Cytoscape.{png,svg}
       <FIG_ROOT>/methods/FigM02b_Networks_Cytoscape_Legend.{png,svg}
       scripts/figures/network_pair_layout.json   (the shared skeleton, frozen)

THE LEGEND IS A SEPARATE FILE, AND SO ARE THE COUNTS
Panel tags, the key and the "30 of 32 concepts shared" line were all on the canvas
and were all costing page area that the drawing needed more: 19 layers of nodes
against a fixed page height means every millimetre spent on furniture comes out of
node size, label size and the air between rows that makes an arrowhead legible.
The key now prints beside the figure from its own file, at 10.5 pt instead of 7.6,
and the counts print to stdout for the caption. Only the two panel titles remain,
because without them the reader cannot tell which network is which.

  python scripts/figures/network_pair_cytoscape.py [--relayout]
"""
from __future__ import annotations

import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# Keep SVG labels as real <text> so the figure stays editable in Illustrator /
# Inkscape, and use the house Arial/Helvetica stack.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib import font_manager                  # noqa: E402
from matplotlib.font_manager import FontProperties   # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402
from matplotlib.patches import PathPatch, Polygon    # noqa: E402
from matplotlib.textpath import TextPath             # noqa: E402
from matplotlib.transforms import Affine2D           # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS, FIG_ROOT                      # noqa: E402

# Data-layer helpers only -- nothing that reads network_diagram's geometry
# constants, which are tuned for nodes a third of this size with labels outside
# them. Importing load_network and canonical_map (rather than re-deriving them)
# guarantees the two figures can never disagree about what a network contains or
# how an ortholog resolves. Importing the module runs only `from arms import ARMS`;
# its main() is __main__-guarded, so nothing is drawn or written.
from network_diagram import (                        # noqa: E402
    SHAPES, TYPE_LABELS, canonical_map, load_network, pick,
)

STYLE_PATH = HERE / "style" / "flashp_style.json"
LAYOUT_PATH = HERE / "network_pair_layout.json"
OUT_DIR = FIG_ROOT / "methods"
STEM = "FigM02b_Networks_Cytoscape"

# =========================================================================== #
# ===============================  EDIT ME  ================================= #
# =========================================================================== #

# The two arms, left to right. A hard pair, intersected with the registry: under
# the full three-arm registry this still draws these two, and under an ARMS_ONLY
# set that drops one it prints why and returns 0 rather than failing the run.
PAIR = ["hand_curated", "flashp_suc"]

# Built at the size it occupies on an A4 PORTRAIT page, so Word never rescales it
# and every point size below is the printed point size. Portrait because the
# skeleton is inherently tall and narrow: 19 layers against 4 columns is the
# network's own longest path, and no layered drawing can be shorter.
FIG_W, FIG_H = 6.9, 9.45      # inches -- A4 portrait, ~0.75 in side margins
DPI = 300

FONT_TITLE = 12.5             # panel titles
FONT_GHOST = 6.2              # the other network's name in an empty slot
FONT_NODE_MAX = 9.0           # node label, tried first
FONT_NODE_MIN = 5.6           # ... and never shrunk below this. Below ~5.5 pt a
                              # printed label stops being readable, so a name that
                              # still does not fit is drawn at the minimum and
                              # allowed to OVERFLOW its node, on a background
                              # plate. That is what Cytoscape itself does -- see
                              # Decapitation_Signal and Strigolactone in the
                              # reference export -- and keeping the node shapes is
                              # worth more here than forcing every name inside,
                              # because shape IS part of the style being matched.

# --- grid ---------------------------------------------------------------- #
# X_GAP is wide relative to the node because the skeleton is only 4 columns
# across: the spare width buys clean vertical corridors at the odd half-multiples
# of X_GAP, which is what keeps 56 edges legible over 19 layers. Y_GAP must leave
# Y_GAP - 2*R_BASE of clear air between rows for the horizontal lanes.
X_GAP = 4.60                  # data units between adjacent node columns
Y_GAP = 3.30                  # data units between adjacent layers
R_BASE = 1.20                 # node radius (GENE / default) in data units
BORDER_LW = 1.25              # node border, points
GHOST_LW = 0.85               # dashed placeholder for a concept this panel lacks
GHOST_COLOR = "#b9c6bf"

EDGE_LW = 0.90                # edge line, points
EDGE_ALPHA_SHARED = 0.95      # edge present in BOTH networks
EDGE_ALPHA_UNIQUE = 0.38      # edge present only in this one -- the delta
HEAD_SCALE = 1.25             # arrowhead / inhibition-bar size multiplier
EDGE_PAD = 0.09               # gap between node border and edge end, data units
CORNER_R = 0.22               # bend fillet radius, data units

# The lane band has to start below the SOURCE row's nodes and finish above the
# NEXT row's, so LANE_ZONE * Y_GAP > R_BASE and (1 - LANE_ZONE) * Y_GAP > R_BASE.
# lane_y() asserts the second; this is the reason both are tied to R_BASE rather
# than being free numbers.
LANE_ZONE = 0.44              # first horizontal lane, as a fraction of the gap
LANE_STEP = 0.175             # lane pitch, data units (auto-compressed to fit)
FAN_SPREAD = 76               # % of the node face that endpoints fan across
CORR_SUB = (0.0, 0.30, -0.30)  # sub-lanes within one vertical corridor
SIDE_PAD = 1.05               # outermost node -> first side corridor
SIDE_STEP = 0.26              # side-corridor pitch
PEEL_SPAN = 8                 # layers spanned before an edge peels out to the side
CLEARANCE = 0.08              # extra gap an edge must keep from a foreign node

PANEL_LEFT, PANEL_RIGHT = 0.015, 0.985
PANEL_BOTTOM, PANEL_TOP = 0.012, 0.952
PANEL_GAP = 0.028             # gutter between panels, figure fraction

# The legend is its own file. Keeping it off the network figure buys ~15% of the
# page height, which goes straight into node size, label size and the air between
# rows that the arrowheads need to be readable. Place the two side by side in Word.
LEGEND_STEM = "FigM02b_Networks_Cytoscape_Legend"
LEGEND_W, LEGEND_H = 6.9, 1.15   # inches -- same width as the figure
LEGEND_NCOL = 4
FONT_LEGEND = 10.5

# Draw the empty slot where a panel's network has no node for a shared-skeleton
# concept, named with the other network's id. This is the most informative element
# in the figure -- it lets the reader SEE a difference rather than infer it.
SHOW_GHOSTS = True

# Fade the edges only one network has, so the 31 shared edges read as the common
# skeleton. Off draws every edge at full strength, i.e. pure Cytoscape.
FADE_UNIQUE_EDGES = True

# ortholog_map.json keys on node type as well as name, so it splits these two.
# They are the same slot -- see PAT_Bud <-> PIN1 in the header.
CONCEPT_ALIAS = {"PAT_Bud": "PIN_1_3_4_7"}

# Edges that run UP the hierarchy. Curated rather than discovered, because a
# generic feedback-arc heuristic breaks whichever edge its traversal order happens
# to reach first -- in testing that was MAX4 -> MAX1 and MAX1 -> CLAMT, destroying
# the biosynthesis chain this figure exists to show. Naming them costs five lines
# and makes the choice reviewable. Any OTHER cycle is a hard failure, which is how
# a network changing underneath this figure gets caught.
FEEDBACK = {
    ("Perception_SL", "MAX3"),   # SL signalling feeding back onto SL biosynthesis
    ("Perception_SL", "MAX4"),   #   (hand-curated only; Flash-P has no such loop)
    ("Perception_SL", "MAX1"),
    ("SMXL6_7_8", "SMXL6_7_8_Transcription"),   # SMXL6,7,8 autorepression
    ("AUXIN_BUD", "PIN_1_3_4_7"),               # induced by the PAT_Bud/PIN1 merge
}

# Nudge a concept after the auto-layout: "STRIGOLACTONE": (x, layer)
POS_OVERRIDE: dict = {}

# Regression gates. Not derived from the files being drawn, so a network changing
# underneath this figure is caught rather than silently redrawn.
EXPECT = {"hand_curated": (32, 56), "flashp_suc": (30, 44)}
EXPECT_SKELETON = {"concepts": 32, "union_edges": 68, "layers": 19, "max_width": 4}

# =========================================================================== #
# ============================  END OF EDIT ME  ============================= #
# =========================================================================== #

BG = "#fbfdfc"
INK = "#13251c"
MUTED = "#5b6b63"
def _measuring_font():
    """The first installed family of the house stack, as a FontProperties.

    Resolved once, to the family that is actually present, for two reasons: a bare
    "sans-serif" string is parsed as a fontconfig pattern and chokes on the hyphen,
    and handing TextPath a stack containing an absent family makes it emit one
    findfont warning PER MEASUREMENT -- thousands of lines for one figure.
    """
    installed = {f.name for f in font_manager.fontManager.ttflist}
    for name in matplotlib.rcParams["font.sans-serif"]:
        if name in installed:
            return FontProperties(family=[name])
    return FontProperties()


FP = _measuring_font()


# --------------------------------------------------------------------------- #
# geometry -- local, because network_diagram's are sized for a different figure
# --------------------------------------------------------------------------- #

def radius_of(style, ntype):
    """Node radius in data units, keeping the house 65 / 80 / 30 size ratios."""
    base = float(style["node"]["size"]["default"])
    return R_BASE * float(pick(style["node"]["size"], ntype)) / base


def half_extent(style, ntype):
    """Half-width / half-height of the node's bounding box, in data units."""
    r = radius_of(style, ntype)
    if pick(style["node"]["shape"], ntype) == "round-rectangle":
        return 1.75 * r, 0.92 * r
    return r, r


def inner_box(style, ntype):
    """Usable width / height for a label drawn INSIDE the node, in data units.

    A diamond or a vee wastes most of its bounding box, so the inscribed box is far
    smaller than the extent. Getting this fraction wrong is what makes labels spill
    silently instead of being routed to the overflow path.
    """
    hw, hh = half_extent(style, ntype)
    k = {"diamond": 0.58, "vee": 0.44, "triangle": 0.44,
         "round-rectangle": 0.90}.get(pick(style["node"]["shape"], ntype), 0.72)
    return hw * k * 2.0, hh * k * 2.0


def unit_patch(style, ntype, xy, scale=1.0, **kw):
    r = radius_of(style, ntype) * scale
    p = SHAPES[pick(style["node"]["shape"], ntype)]
    return PathPatch(p.transformed(Affine2D().scale(r).translate(*xy)), **kw)


def _fillet(pts, r):
    """Round every interior corner of an axis-aligned polyline."""
    if len(pts) < 3 or r <= 0:
        return pts
    out = [pts[0]]
    for i in range(1, len(pts) - 1):
        a, b, c = pts[i - 1], pts[i], pts[i + 1]
        d0 = math.hypot(b[0] - a[0], b[1] - a[1])
        d1 = math.hypot(c[0] - b[0], c[1] - b[1])
        k = min(r, d0 / 2.0, d1 / 2.0)
        if k < 1e-9 or d0 < 1e-9 or d1 < 1e-9:
            out.append(b)
            continue
        u0 = ((b[0] - a[0]) / d0, (b[1] - a[1]) / d0)
        u1 = ((c[0] - b[0]) / d1, (c[1] - b[1]) / d1)
        p0 = (b[0] - u0[0] * k, b[1] - u0[1] * k)
        p1 = (b[0] + u1[0] * k, b[1] + u1[1] * k)
        out.append(p0)
        for t in (0.25, 0.5, 0.75):
            out.append(((1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * b[0] + t * t * p1[0],
                        (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * b[1] + t * t * p1[1]))
        out.append(p1)
    out.append(pts[-1])
    return out


def _pack(items):
    """Greedy interval colouring. items = [(lo, hi, key)] -> {key: lane}, n_lanes."""
    lanes: list[list[tuple[float, float]]] = []
    got = {}
    for lo, hi, key in sorted(items, key=lambda z: (z[0], z[1])):
        for k, used in enumerate(lanes):
            if all(hi <= a or lo >= b for a, b in used):
                used.append((lo, hi))
                got[key] = k
                break
        else:
            lanes.append([(lo, hi)])
            got[key] = len(lanes) - 1
    return got, len(lanes)


# --------------------------------------------------------------------------- #
# the shared skeleton
# --------------------------------------------------------------------------- #

def union_graph(nets, canon):
    """(concepts, union edge set) over BOTH networks, at the concept level."""
    V, U = set(), set()
    for a, (nodes, edges) in nets.items():
        V |= {canon[a][n] for n in nodes}
        U |= {(canon[a][s], canon[a][t]) for s, t, _g in edges}
    return V, U


def assign_layers(V, U):
    """Longest path over the DAG, then tighten. -> {concept: layer}, n_layers.

    The tighten sweep is what makes the drawing sparse. Raw longest path parks every
    node with no predecessor on layer 0, so CXE, SMXL6_7_8_Transcription and PIN1
    all start at the top with 6-to-12-layer edges hanging off them and the widest
    row is 6. Tightening slides each non-sink down to just above its earliest
    successor, which drops CXE beside LBO, SMXL mRNA beside SMXL6,7,8 and PIN1
    beside BRC1 -- 19 layers, never wider than 4.
    """
    D = [(s, t) for s, t in U if s != t and (s, t) not in FEEDBACK]

    adj = defaultdict(list)
    for s, t in D:
        adj[s].append(t)
    color, bad = {}, []

    def walk(u):
        color[u] = 1
        for v in adj[u]:
            if color.get(v, 0) == 0:
                walk(v)
            elif color[v] == 1:
                bad.append((u, v))
        color[u] = 2

    for u in sorted(V):
        if color.get(u, 0) == 0:
            walk(u)
    if bad:
        raise SystemExit(
            "the union graph still has a cycle after FEEDBACK was removed; add "
            "the edge that closes it (with a reason) to FEEDBACK:\n  "
            + "\n  ".join(f"{s} -> {t}" for s, t in sorted(set(bad))))

    pred, succ = defaultdict(list), defaultdict(list)
    for s, t in D:
        pred[t].append(s)
        succ[s].append(t)

    layer = {c: 0 for c in V}
    for _ in range(len(V)):
        for s, t in D:
            if layer[t] < layer[s] + 1:
                layer[t] = layer[s] + 1

    for _ in range(40):
        moved = False
        for c in sorted(V):
            if not succ[c]:
                continue
            lo = max([layer[p] + 1 for p in pred[c]], default=0)
            hi = min([layer[s] - 1 for s in succ[c]], default=layer[c])
            if hi < lo:
                continue
            best = min(range(lo, hi + 1),
                       key=lambda L: (sum(abs(L - layer[x])
                                          for x in pred[c] + succ[c]), -L))
            if best != layer[c]:
                layer[c] = best
                moved = True
        if not moved:
            break

    used = sorted(set(layer.values()))
    squash = {d: i for i, d in enumerate(used)}
    return {c: squash[v] for c, v in layer.items()}, len(used)


def order_columns(layer, n_layers, U):
    """Barycentre sweep over the UNION edge set -> an x per concept.

    Two deliberate choices. The sweep is driven by BOTH networks' edges at once, so
    neither panel is untangled at the other's expense -- which is what would happen
    if each were laid out from its own edges and the results averaged. And a
    neighbour's pull is weighted 1/span, so a fourteen-layer edge does not drag a
    node sideways as hard as a one-layer edge; that is what keeps the biosynthesis
    chain straight despite sucrose's five long-haul edges.
    """
    rows = [sorted(c for c in layer if layer[c] == r) for r in range(n_layers)]
    nb = defaultdict(list)
    for s, t in U:
        if s == t or layer[s] == layer[t]:
            continue
        w = 1.0 / abs(layer[s] - layer[t])
        nb[s].append((t, w))
        nb[t].append((s, w))

    def slots(rs):
        out = {}
        for members in rs:
            off = -(len(members) - 1) / 2.0
            for i, c in enumerate(members):
                out[c] = (off + i) * X_GAP
        return out

    for it in range(60):
        pos = slots(rows)
        order = range(n_layers) if it % 2 == 0 else range(n_layers - 1, -1, -1)
        for r in order:
            key = {}
            for c in rows[r]:
                near = nb[c]
                key[c] = (sum(pos[m] * w for m, w in near) / sum(w for _m, w in near)
                          if near else pos[c])
            rows[r] = sorted(rows[r], key=lambda c: (key[c], c))
    return slots(rows)


def build_skeleton(nets, canon, relayout):
    """The one coordinate skeleton both panels are drawn on, seeded then frozen."""
    V, U = union_graph(nets, canon)
    layer, n_layers = assign_layers(V, U)
    width = max(Counter(layer.values()).values())

    got = {"concepts": len(V), "union_edges": len(U),
           "layers": n_layers, "max_width": width}
    if got != EXPECT_SKELETON:
        raise SystemExit(
            f"skeleton changed: {got} != {EXPECT_SKELETON} -- "
            f"update EXPECT_SKELETON if intended")

    xs = order_columns(layer, n_layers, U)
    seeded = {c: (round(xs[c], 4), layer[c]) for c in V}

    saved = {}
    if LAYOUT_PATH.is_file() and not relayout:
        blob = json.loads(LAYOUT_PATH.read_text(encoding="utf-8"))
        saved = blob.get("skeleton", {})
    if set(saved) == set(seeded):
        skel = {c: (float(saved[c][0]), int(saved[c][1])) for c in seeded}
        print(f"  using the skeleton from {LAYOUT_PATH.name}")
    else:
        skel = seeded
        LAYOUT_PATH.write_text(json.dumps({
            "_note": "concept -> [x, layer]. ONE skeleton, both panels: editing an "
                     "entry moves that node in BOTH. --relayout to re-seed.",
            "feedback_edges": sorted(list(e) for e in FEEDBACK),
            "skeleton": {c: [v[0], v[1]] for c, v in sorted(skel.items())},
        }, indent=2) + "\n", encoding="utf-8")
        print(f"  wrote {LAYOUT_PATH.name} "
              f"(edit it to nudge a concept; --relayout to re-seed)")

    skel.update({c: tuple(v) for c, v in POS_OVERRIDE.items() if c in skel})
    pos = {c: (x, -L * Y_GAP) for c, (x, L) in skel.items()}
    row = {c: L for c, (_x, L) in skel.items()}
    return pos, row


def assert_registered(nets, canon, pos):
    """The figure's central claim, asserted: a concept is at ONE point on the page."""
    seen: dict = {}
    for a, (nodes, _e) in nets.items():
        for n in nodes:
            c = canon[a][n]
            if c in seen and seen[c] != pos[c]:
                raise SystemExit(f"concept {c!r} is not co-positioned")
            seen[c] = pos[c]
        slots = [pos[canon[a][n]] for n in nodes]
        if len(set(slots)) != len(slots):
            dup = [c for c, k in Counter(canon[a][n] for n in nodes).items() if k > 1]
            raise SystemExit(f"{a}: two nodes share the skeleton slot(s) {dup}")


# --------------------------------------------------------------------------- #
# edge routing -- axis-aligned, on the layer grid the skeleton defines
# --------------------------------------------------------------------------- #

def route_edges(style, nodes, edges, pos, row):
    """Fixed node positions -> one axis-aligned polyline per edge."""
    ext = {n: half_extent(style, nodes[n]["type"]) for n in nodes}
    xs = [pos[n][0] for n in nodes]
    x_lo, x_hi, x_mid = min(xs), max(xs), (min(xs) + max(xs)) / 2.0

    recs = []
    for i, (s, t, sign) in enumerate(edges):
        span = row[t] - row[s]
        recs.append({"i": i, "s": s, "t": t, "sign": sign, "span": span,
                     "kind": ("self" if s == t else
                              "flat" if span == 0 else
                              "direct" if abs(span) == 1 else
                              "peel" if abs(span) >= PEEL_SPAN else "corridor")})

    # --- 1. fan the endpoints across the node faces ------------------------ #
    groups: dict = {}
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
            off = 0.0 if len(grp) == 1 else \
                (j / (len(grp) - 1.0) - 0.5) * (FAN_SPREAD / 100.0) * 2 * hw
            r["off_src" if role == "src" else "off_tgt"] = off

    for r in recs:
        if r["kind"] == "self":
            continue
        (sx, sy), (tx, ty) = pos[r["s"]], pos[r["t"]]
        sh, th = ext[r["s"]][1], ext[r["t"]][1]
        down_out, down_in = r["span"] >= 0, r["span"] > 0
        r["p0"] = (sx + r.get("off_src", 0.0),
                   sy - sh - EDGE_PAD if down_out else sy + sh + EDGE_PAD)
        r["p1"] = (tx + r.get("off_tgt", 0.0),
                   ty + th + EDGE_PAD if down_in else ty - th - EDGE_PAD)
        rs, rt = row[r["s"]], row[r["t"]]
        if r["span"] >= 0:
            r["g1"], r["g2"] = rs, (rt - 1 if r["span"] > 1 else rs)
        else:
            r["g1"], r["g2"] = rs - 1, rt

    # --- 2. drop the multi-layer edges through free in-field corridors ----- #
    slots_ = sorted({round((2 * k + 1) * X_GAP / 2.0 + d, 4)
                     for k in range(int((x_lo - X_GAP) / X_GAP) - 1,
                                    int((x_hi + X_GAP) / X_GAP) + 1)
                     for d in CORR_SUB})
    busy = {c: [] for c in slots_}
    for r in sorted((r for r in recs if r["kind"] == "corridor"),
                    key=lambda r: (abs(r["span"]), r["i"])):
        want = r["p0"][0]
        lo, hi = min(r["g1"], r["g2"]), max(r["g1"], r["g2"])
        free = [c for c in slots_ if all(hi < a or lo > b for a, b in busy[c])]
        if not free:
            r["kind"] = "peel"
            continue
        c = min(free, key=lambda c: (abs(c - want), c))
        busy[c].append((lo, hi))
        r["mid_x"] = c

    # --- 3. side corridors for the long haul ------------------------------- #
    sides: dict = {"L": [], "R": []}
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

    # --- 4. horizontal lanes, packed per layer gap ------------------------- #
    runs: dict = {}
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

    def lane_y(g, k, n):
        """y of horizontal lane k of n, inside the gap below layer g."""
        top = -g * Y_GAP - Y_GAP * LANE_ZONE
        room = Y_GAP * (1.0 - LANE_ZONE) - R_BASE - EDGE_PAD
        return top - k * min(LANE_STEP, max(room, 0.02) / max(n, 1))

    # --- 5. build the polylines -------------------------------------------- #
    out = []
    for r in recs:
        if r["kind"] == "self":
            (x, y) = pos[r["s"]]
            hw, hh = ext[r["s"]]
            w, h = hw + 0.60, hh + 0.48
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
        clean = [pts[0]]
        for p in pts[1:]:
            if math.hypot(p[0] - clean[-1][0], p[1] - clean[-1][1]) > 1e-9:
                clean.append(p)
        out.append((_fillet(clean, CORNER_R), r["sign"]))
    return out


def foreign_hits(style, nodes, pos, edges, routes):
    """Routed edges that pass over a node that is not one of their endpoints."""
    ext = {n: half_extent(style, nodes[n]["type"]) for n in nodes}
    hits = []
    for (pts, _sg), (s, t, _x) in zip(routes, edges):
        for n, (x, y) in pos.items():
            if n in (s, t):
                continue
            hw, hh = ext[n][0] + CLEARANCE, ext[n][1] + CLEARANCE
            if any(abs(px - x) <= hw and abs(py - y) <= hh for px, py in pts):
                hits.append((s, t, n))
                break
    return hits


def draw_route(ax, style, pts, sign, alpha):
    key = "negative" if sign < 0 else "positive"
    ax.add_line(Line2D([q[0] for q in pts], [q[1] for q in pts],
                       color=pick(style["edge"]["line"], key), lw=EDGE_LW,
                       alpha=alpha, solid_capstyle="round",
                       solid_joinstyle="round", zorder=1))
    # Head geometry ported from flashp-prop-core.js::headD -- a filled triangle for
    # activation, a perpendicular bar for inhibition.
    head_c = pick(style["edge"]["arrowColor"], key)
    tip, back = pts[-1], pts[max(0, len(pts) - 4)]
    m = math.hypot(tip[0] - back[0], tip[1] - back[1]) or 1.0
    ca, sa = (tip[0] - back[0]) / m, (tip[1] - back[1]) / m
    u = R_BASE * HEAD_SCALE
    if sign > 0:
        length, half = 0.50 * u, 0.27 * u
        ax.add_patch(Polygon(
            [(tip[0], tip[1]),
             (tip[0] - length * ca + half * sa, tip[1] - length * sa - half * ca),
             (tip[0] - length * ca - half * sa, tip[1] - length * sa + half * ca)],
            closed=True, facecolor=head_c, edgecolor="none", alpha=alpha, zorder=1.5))
    else:
        bar = 0.39 * u
        ax.add_line(Line2D([tip[0] + bar * sa, tip[0] - bar * sa],
                           [tip[1] - bar * ca, tip[1] + bar * ca],
                           color=head_c, lw=EDGE_LW + 0.9, alpha=alpha,
                           solid_capstyle="butt", zorder=1.5))


# --------------------------------------------------------------------------- #
# labels inside the nodes
# --------------------------------------------------------------------------- #

def _text_w(s, size):
    """Width of one line in POINTS, measured from the glyph outlines."""
    return TextPath((0, 0), s, size=size, prop=FP).get_extents().width if s else 0.0


def _splits(name):
    """Candidate line-breakings of a node id, fewest lines first.

    Breaking on '_' only -- never renaming, never abbreviating, because the whole
    point of the figure is that each panel shows its own network's spelling. It is
    what keeps SMXL6_7_8_Transcription, Decapitation_Signal and Sustained_growth
    inside their nodes: the parts are meaningful on their own lines.
    """
    parts = [p for p in name.split("_") if p]
    out = [[name]]
    for k in (2, 3):
        if len(parts) < k:
            continue
        per, lines, cur = math.ceil(len(name) / k), [], ""
        for p in parts:
            if cur and len(cur) + 1 + len(p) > per and len(lines) < k - 1:
                lines.append(cur)
                cur = p
            else:
                cur = f"{cur}_{p}" if cur else p
        lines.append(cur)
        if len(lines) == k and lines not in out:
            out.append(lines)
    return out


def fit_label(name, box_w_pt, box_h_pt):
    """(lines, fontsize, overflowed) for `name` inside a box_w x box_h point box."""
    best = None
    for lines in _splits(name):
        size = FONT_NODE_MAX
        while size >= FONT_NODE_MIN - 1e-9:
            w = max(_text_w(ln, size) for ln in lines)
            h = len(lines) * size * 1.16
            if w <= box_w_pt and h <= box_h_pt:
                return lines, size, False
            over = max(w / box_w_pt, h / box_h_pt)
            if best is None or over < best[0]:
                best = (over, lines, size)
            size -= 0.25
    return best[1], best[2], True


# --------------------------------------------------------------------------- #
# panel
# --------------------------------------------------------------------------- #

def draw_panel(ax, style, nodes, edges, pos, row, shared_edges, ghosts,
               xlim, ylim, pt_per_unit):
    ax.set_facecolor(BG)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")

    if SHOW_GHOSTS:
        for ntype, xy, name in ghosts:
            ax.add_patch(unit_patch(style, ntype, xy, 1.0, facecolor="none",
                                    edgecolor=GHOST_COLOR, lw=GHOST_LW,
                                    linestyle=(0, (1.6, 2.0)), zorder=2))
            ax.text(xy[0], xy[1], "\n".join(_splits(name)[-1]), ha="center",
                    va="center", fontsize=FONT_GHOST, color=GHOST_COLOR,
                    style="italic", zorder=2.5, linespacing=1.05)

    routes = route_edges(style, nodes, edges, pos, row)
    for (pts, sign), e in zip(routes, edges):
        a = EDGE_ALPHA_SHARED
        if FADE_UNIQUE_EDGES and e not in shared_edges:
            a = EDGE_ALPHA_UNIQUE
        draw_route(ax, style, pts, sign, a)

    for n, meta in nodes.items():
        ty = meta["type"]
        xy = pos[n]
        ax.add_patch(unit_patch(style, ty, xy, 1.0,
                                facecolor=pick(style["node"]["fill"], ty),
                                edgecolor=pick(style["node"]["border"], ty),
                                lw=BORDER_LW, zorder=3))
        bw, bh = inner_box(style, ty)
        lines, size, spill = fit_label(n, bw * pt_per_unit, bh * pt_per_unit)
        # A label too long for its shape overflows on a plate, as Cytoscape does,
        # rather than being shrunk past legibility or renamed.
        box = dict(boxstyle="round,pad=0.12", facecolor=BG, edgecolor="none",
                   alpha=0.82) if spill else None
        ax.text(xy[0], xy[1], "\n".join(lines), ha="center", va="center",
                fontsize=size, color=style["node"].get("labelColor", INK),
                zorder=4, linespacing=1.05, bbox=box)

    return foreign_hits(style, nodes, pos, edges, routes)


def legend_handles(style, types):
    h, lab = [], []
    for ty in ("GENE", "HORMONE", "METABOLITE", "PROTEIN_COMPLEX", "PROCESS",
               "ENVIRONMENT", "PHENOTYPE"):
        if ty not in types:
            continue
        marker = {"ellipse": "o", "diamond": "D", "vee": "v",
                  "round-rectangle": "s",
                  "triangle": "^"}[pick(style["node"]["shape"], ty)]
        h.append(Line2D([], [], marker=marker, linestyle="none", markersize=7,
                        markerfacecolor=pick(style["node"]["fill"], ty),
                        markeredgecolor=pick(style["node"]["border"], ty),
                        markeredgewidth=1.1))
        lab.append(TYPE_LABELS[ty])
    missing = types - {"GENE", "HORMONE", "METABOLITE", "PROTEIN_COMPLEX",
                       "PROCESS", "ENVIRONMENT", "PHENOTYPE"}
    if missing:
        raise SystemExit(f"node type(s) drawn but absent from the legend: {missing}")
    h.append(Line2D([], [], color=pick(style["edge"]["line"], "positive"), lw=2.2))
    lab.append("activates  (→)")
    h.append(Line2D([], [], color=pick(style["edge"]["line"], "negative"), lw=2.2))
    lab.append("inhibits  (—|)")
    if FADE_UNIQUE_EDGES:
        h.append(Line2D([], [], color=MUTED, lw=2.2, alpha=EDGE_ALPHA_UNIQUE))
        lab.append("edge in this network only")
    if SHOW_GHOSTS:
        h.append(Line2D([], [], marker="o", linestyle="none", markersize=7,
                        markerfacecolor="none", markeredgecolor=GHOST_COLOR,
                        markeredgewidth=1.0))
        lab.append("slot the other network fills")
    return h, lab


def write_legend(style, types):
    """The key, as its own file, at a size meant to be read.

    Split out so the network figure gets the whole page. It is drawn from the same
    style file and the same legend_handles() the panels are built from, so it can
    never describe a shape or a colour the figure does not use.
    """
    fig = plt.figure(figsize=(LEGEND_W, LEGEND_H), dpi=DPI)
    fig.patch.set_facecolor("white")
    h, lab = legend_handles(style, types)
    fig.legend(h, lab, loc="center", bbox_to_anchor=(0.5, 0.5),
               ncol=LEGEND_NCOL, frameon=False, fontsize=FONT_LEGEND,
               handletextpad=0.7, columnspacing=1.9, labelspacing=0.95)
    for ext in ("png", "svg"):
        p = OUT_DIR / f"{LEGEND_STEM}.{ext}"
        fig.savefig(p, dpi=DPI, facecolor="white")
        print(f"  wrote {p.relative_to(ROOT).as_posix()}")
    plt.close(fig)


# --------------------------------------------------------------------------- #

def main():
    relayout = "--relayout" in sys.argv
    missing = [a for a in PAIR if a not in ARMS]
    if missing:
        print(f"network_pair_cytoscape: {', '.join(missing)} not in the registry "
              f"(ARMS_ONLY={','.join(ARMS)}) -- nothing to draw")
        return

    style = json.loads(STYLE_PATH.read_text(encoding="utf-8"))
    nets = {a: load_network(ARMS[a]) for a in PAIR}
    for a, (nodes, edges) in nets.items():
        want = EXPECT.get(a)
        if want and (len(nodes), len(edges)) != want:
            raise SystemExit(
                f"{a}: network is {len(nodes)} nodes / {len(edges)} edges, "
                f"expected {want[0]} / {want[1]} -- update EXPECT if intended")

    canon = canonical_map(nets)
    for a in PAIR:
        for n, c in list(canon[a].items()):
            canon[a][n] = CONCEPT_ALIAS.get(c, c)

    pos_c, row_c = build_skeleton(nets, canon, relayout)
    assert_registered(nets, canon, pos_c)

    # Concept-level edge sets, for the shared / panel-only split.
    csets = {a: {(canon[a][s], canon[a][t], g) for s, t, g in nets[a][1]}
             for a in PAIR}
    shared_c = csets[PAIR[0]] & csets[PAIR[1]]
    all_concepts = set(pos_c)
    have = {a: {canon[a][n] for n in nets[a][0]} for a in PAIR}
    both = have[PAIR[0]] & have[PAIR[1]]

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)
    fig.patch.set_facecolor("white")

    # One shared data window, so the two panels really are in register.
    xs = [xy[0] for xy in pos_c.values()]
    ys = [xy[1] for xy in pos_c.values()]
    xlim = (min(xs) - X_GAP * 0.55, max(xs) + X_GAP * 0.55)
    ylim = (min(ys) - Y_GAP * 0.95, max(ys) + Y_GAP * 0.95)

    span_w = (PANEL_RIGHT - PANEL_LEFT - PANEL_GAP) / 2.0
    pt_per_unit = (span_w * FIG_W * 72.0) / (xlim[1] - xlim[0])

    types_seen, all_hits = set(), []
    for k, a in enumerate(PAIR):
        nodes, edges = nets[a]
        other = PAIR[1 - k]
        pos = {n: pos_c[canon[a][n]] for n in nodes}
        row = {n: row_c[canon[a][n]] for n in nodes}
        ghosts = []
        for c in sorted(all_concepts - have[a]):
            src = next(n for n in nets[other][0] if canon[other][n] == c)
            ghosts.append((nets[other][0][src]["type"], pos_c[c], src))
        shared_here = {(s, t, g) for s, t, g in edges
                       if (canon[a][s], canon[a][t], g) in shared_c}
        types_seen |= {m["type"] for m in nodes.values()}

        left = PANEL_LEFT + k * (span_w + PANEL_GAP)
        ax = fig.add_axes([left, PANEL_BOTTOM, span_w, PANEL_TOP - PANEL_BOTTOM])
        hits = draw_panel(ax, style, nodes, edges, pos, row, shared_here, ghosts,
                          xlim, ylim, pt_per_unit)
        all_hits += [(a,) + h for h in hits]

        ax.set_title(ARMS[a].fig_label.replace(" (", "\n("), fontsize=FONT_TITLE,
                     fontweight="bold", color=INK, pad=7)
        print(f"  {a:<14}{len(nodes)} nodes, {len(edges)} edges, "
              f"{len(shared_here)} shared, {len(ghosts)} empty slot(s), "
              f"{len(hits)} edge(s) over a foreign node")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        p = OUT_DIR / f"{STEM}.{ext}"
        fig.savefig(p, dpi=DPI, facecolor="white")
        print(f"  wrote {p.relative_to(ROOT).as_posix()}")
    plt.close(fig)

    write_legend(style, types_seen)

    # The counts the figure used to print along its bottom edge. They belong in
    # the caption, not on the canvas, so they are reported here instead of being
    # paid for in page area.
    print(f"  caption numbers: node positions identical in both panels; "
          f"{len(both)} of {len(all_concepts)} concepts shared; {len(shared_c)} of "
          f"{len(csets[PAIR[0]])} / {len(csets[PAIR[1]])} edges shared")

    if all_hits:
        print(f"  NOTE {len(all_hits)} edge(s) pass over a foreign node; "
              f"nudge via POS_OVERRIDE if it reads badly")


if __name__ == "__main__":
    main()
