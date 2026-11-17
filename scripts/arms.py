"""The arm registry — one source of truth for which networks the pipeline runs.

Reads scripts/arms.json and exposes it as frozen dataclasses with ROOT-resolved
absolute paths. run_all.ps1 reads the same JSON directly, so the arm list is never
written down twice.

An "arm" is one network carried through the whole PSoup -> APSIM chain. Everything
downstream of the network is identical between arms by construction, so any
difference in the figures is attributable to the network.

TWO ENVIRONMENT HOOKS, both for building a figure set over a SUBSET of the arms
without editing arms.json (see scripts/make_2network_figures.py):

  ARMS_ONLY   comma-separated arm names. Restricts ARMS to those, in registry
              order. Unset means every registered arm, which is the normal case.
  FIGS_ROOT   figure output root, relative to ROOT. Defaults to "04_figures".
              Every figure script builds its own output directory under FIG_ROOT,
              so one variable redirects the whole set.

They are environment variables rather than arguments because ARMS is module-level
and read at import time: one process can only ever hold one configuration, so a
different arm set means a different subprocess. make_figures.py already documents
the same constraint for scripts/figures/_figure_common.py.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
_REGISTRY = ROOT / "scripts" / "arms.json"

# Figure output root. Scripts join their own subdirectory onto this rather than
# onto a hardcoded "04_figures", so a whole alternative figure set lands beside
# the default one instead of overwriting it.
FIG_ROOT = ROOT / os.environ.get("FIGS_ROOT", "04_figures")

# True when ARMS_ONLY actually dropped an arm. A figure whose canvas is sized for
# the full registry can use this to reclaim the margin a shorter set leaves behind,
# without changing anything about the full-registry output.
RESTRICTED = False


@dataclass(frozen=True)
class Arm:
    name: str
    label: str
    kind: str                    # "r_psoup" | "python_flashp"
    base_network_dir: Path | None  # pristine FlashP export the _suc copy is built FROM
    network_dir: Path | None
    equations: Path | None       # algebraic equations APSIM parameters derive from
    sweep_network_dir: Path | None
    response_csv: Path
    pattern_column: str
    solver_out_stem: Path | None
    phenotype: str
    sucrose_node: str
    gene_map: str | None         # variant key into gene_map.VARIANTS
    audit_key: str | None
    extend_sucrose: bool
    color: str
    nodes: int | None            # counted from the network, not declared
    edges: int | None

    @property
    def is_flashp(self) -> bool:
        return self.kind == "python_flashp"

    @property
    def size_label(self) -> str:
        """"(30 nodes, 44 edges)" -- empty if the network could not be counted."""
        if self.nodes is None or self.edges is None:
            return ""
        return f"({self.nodes} nodes, {self.edges} edges)"

    @property
    def fig_label(self) -> str:
        """Name plus size, for figure legends and column headers.

        Both FlashP arms are called "Flash-P"; the size is what tells them apart,
        which is the point -- the reader sees the topology the arm was built on.
        """
        return f"{self.label} {self.size_label}".strip()

    @property
    def fig_label_2line(self) -> str:
        """Same, broken so it fits a rotated row label or a narrow column head."""
        return "\n".join(p for p in (self.label, self.size_label) if p)


def _path(rel: str | None) -> Path | None:
    return None if rel is None else ROOT / rel


def _rows(csv: Path) -> int:
    """Data rows in a plain one-record-per-line CSV (header dropped)."""
    lines = csv.read_text(encoding="utf-8-sig").splitlines()
    return sum(1 for ln in lines[1:] if ln.strip())


def _size(net: Path | None) -> tuple[int | None, int | None]:
    """Node and edge counts, read from the network itself.

    Counted rather than declared in arms.json so a figure label can never drift
    from the topology it describes. FlashP exports carry the graph in
    network/network.json; the hand-curated PSoup model carries it as two CSVs.
    """
    if net is None:
        return None, None
    nj = net / "network" / "network.json"
    if nj.exists():
        j = json.loads(nj.read_text(encoding="utf-8"))
        return len(j["nodes"]), len(j["edges"])
    nodes, edges = net / "PSoup_Network_nodes.csv", net / "Edges.csv"
    if nodes.exists() and edges.exists():
        return _rows(nodes), _rows(edges)
    return None, None


def _restrict(arms: dict[str, Arm], reference: str) -> dict[str, Arm]:
    """ARMS_ONLY applied, or `arms` unchanged when it is unset.

    Registry order is preserved regardless of the order the names are given in, so
    a figure's panel order never depends on how the variable was typed. The
    reference arm has to survive: every pairwise figure compares against it, and
    losing it here would surface much later as a confusing failure inside
    others_present().
    """
    raw = os.environ.get("ARMS_ONLY", "").strip()
    if not raw:
        return arms
    keep = [n.strip() for n in raw.split(",") if n.strip()]
    bad = [n for n in keep if n not in arms]
    if bad:
        raise SystemExit(
            f"ARMS_ONLY names unknown arm(s): {', '.join(bad)} -- "
            f"choose from {', '.join(arms)}")
    if reference not in keep:
        raise SystemExit(
            f"ARMS_ONLY must keep the reference arm {reference!r} -- "
            f"every pairwise figure is drawn against it")
    global RESTRICTED
    RESTRICTED = len(keep) < len(arms)
    return {n: arms[n] for n in arms if n in keep}


def _load() -> tuple[dict[str, Arm], str]:
    raw = json.loads(_REGISTRY.read_text(encoding="utf-8"))
    arms = {}
    for a in raw["arms"]:
        nodes, edges = _size(_path(a["network_dir"]))
        arms[a["name"]] = Arm(
            name=a["name"], label=a["label"], kind=a["kind"],
            base_network_dir=_path(a.get("base_network_dir")),
            network_dir=_path(a["network_dir"]),
            equations=_path(a["equations"]),
            sweep_network_dir=_path(a["sweep_network_dir"]),
            response_csv=_path(a["response_csv"]),
            pattern_column=a["pattern_column"],
            solver_out_stem=_path(a["solver_out_stem"]),
            phenotype=a["phenotype"], sucrose_node=a["sucrose_node"],
            gene_map=a["gene_map"], audit_key=a["audit_key"],
            extend_sucrose=a["extend_sucrose"], color=a["color"],
            nodes=nodes, edges=edges,
        )
    return _restrict(arms, raw["reference"]), raw["reference"]


ARMS, REFERENCE = _load()


def arm_names() -> list[str]:
    return list(ARMS)


def arm(name: str) -> Arm:
    try:
        return ARMS[name]
    except KeyError:
        raise SystemExit(
            f"unknown arm {name!r} -- choose from {', '.join(ARMS)}") from None


def flashp_arms() -> list[Arm]:
    return [a for a in ARMS.values() if a.is_flashp]


def select(argv: list[str], kind: str | None = None) -> list[Arm]:
    """Arms named on the command line, else all of them (optionally filtered by kind).

    Every consumer takes an optional arm filter so that `run_all.ps1 -Arm <x>` never
    touches another arm's outputs -- build_apsim_runs.py in particular DELETES the .db
    beside each .apsimx it stages, so an unfiltered run there is destructive.
    """
    pool = [a for a in ARMS.values() if kind is None or a.kind == kind]
    if not argv:
        return pool
    bad = [n for n in argv if n not in ARMS]
    if bad:
        raise SystemExit(
            f"unknown arm(s): {', '.join(bad)} -- choose from {', '.join(ARMS)}")
    picked = [ARMS[n] for n in argv]
    wrong = [a.name for a in picked if kind is not None and a.kind != kind]
    if wrong:
        raise SystemExit(f"arm(s) {', '.join(wrong)} are not of kind {kind!r}")
    return picked
