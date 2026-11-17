"""Add the missing sucrose pathway to a LOCAL copy of a FlashP network.

Usage:  python extend_sucrose.py [arm ...] [--with-grass-twins]
        (default: every arm flagged extend_sucrose in scripts/arms.json)

Nothing in FlashP_Repo is touched, and neither is the pristine export this reads.
Each arm's `base_network_dir` is copied to its `network_dir` and the edges are added
to the COPY:

    02_networks/Tillering_Sorghum_Benchmark  ->  ..._Benchmark_suc   (flashp_suc)
    02_networks/Shoot_Branching_Sorghum      ->  ..._Sorghum_suc     (flashp_shoot)

The hand-curated network is never passed here. It is the reference arm, it already
carries all five of its sucrose edges natively, and rewiring it would move the
yardstick every other arm is measured against.

WHY
Both FlashP builds under-wire sucrose relative to the hand-curated PSoup model,
which has five outgoing SUC edges:

  benchmark (30 nodes)  2 edges: Sucrose -| BRC1, Sucrose -> Bud_Release
                        -> compressed the WT amplitude to 0.745..1.265
                           (hand-curated: 0.317..2.059)
  shoot (69 nodes)      3 edges: SUC -| SL, SUC -| BRC1_2, SUC -> Bud_release
                        -> and NONE of them reaches the phenotype directly, so
                           sucrose moved Shoot_Branching by only +-4.5% and all
                           32 genotypes collapsed onto near-identical flat lines

With so few routes the mutants cannot differentiate by sucrose at all. That is a
wiring gap, not a biological result, so it is closed here rather than reported.

Each added edge carries a literature DOI and the sentence that supports it, in
the same spirit as FlashP's own evidence grounding. Species is recorded so the
Arabidopsis-vs-grass provenance stays visible.

The benchmark arm additionally flips Cytokinin from a source node to a regulated
one: it had NO upstream regulator at all, which is why the whole cytokinin arm was
inert. The shoot network's CK is already regulated (IPT ->, LOG ->, CKX -|), so
nothing flips there -- the flip is driven by what actually gained a regulator, not
declared per arm.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from arms import Arm, select  # noqa: E402

# --------------------------------------------------------------------------
# The sucrose edges to add, per arm, each with its evidence.
#
# The two FlashP arms get the SAME four biological routes -- direct metabolic
# supply to the phenotype, the two SL-signalling arms, and cytokinin -- so that
# neither network is advantaged by a richer sucrose wiring than the other. Only
# the node NAMES differ, because the two builds use different label conventions.
# --------------------------------------------------------------------------
_SUPPLY = {
    "d": "10.1111/pce.12500",
    "species": "Sorghum bicolor",
    "title": "Photosynthetic leaf area modulates tiller bud outgrowth in sorghum",
    "authors": "Kebrom TH, Mullet JE",
    "year": 2015, "journal": "Plant, Cell & Environment 38(8):1471-1478",
    "evidence": ("The SWEET transporter may facilitate sucrose unloading from the phloem "
                 "to the apoplast where cell wall invertases generate monosaccharides for "
                 "uptake and utilization to sustain bud outgrowth."),
}
_MAX2 = {
    "d": "10.1111/nph.17834",
    "species": "Oryza sativa",
    "title": "Sucrose promotes D53 accumulation and tillering in rice",
    "authors": "Patil SB, Barbier FF, Zhao J, et al.",
    "year": 2022, "journal": "New Phytologist 234:122-136",
    "evidence": ("Sucrose inhibits the gene expression of D3, the orthologue of the "
                 "Arabidopsis F-box MAX2 required for SL signalling."),
    "note": "Grass-direct. D3 is the rice/sorghum orthologue of MAX2 (= SbMAX2).",
}
_SMXL = {
    "d": "10.1111/nph.17834",
    "species": "Oryza sativa",
    "title": "Sucrose promotes D53 accumulation and tillering in rice",
    "authors": "Patil SB, Barbier FF, Zhao J, et al.",
    "year": 2022, "journal": "New Phytologist 234:122-136",
    "evidence": ("sucrose antagonizes the suppressive action of SL on tillering in rice "
                 "and on the degradation of D53, a major target of SL signalling."),
    "note": ("Sucrose stabilises D53/SMXL6,7,8 by opposing SL-driven degradation. "
             "Grass-direct, same study."),
}
_CK = {
    "d": "10.1093/plphys/kiab003",
    "species": "Solanum tuberosum / Pisum sativum",
    "title": "Sucrose promotes stem branching through cytokinin",
    "authors": "Salam BB, Barbier F, Danieli R, et al.",
    "year": 2021, "journal": "Plant Physiology 185(4):1708-1721",
    "evidence": ("Sucrose feeding to detached stems promotes the accumulation of "
                 "cytokinin (CK), as well as the expression of vacuolar invertase (VInv)."),
}

EDGES: dict[str, list[dict]] = {
    # ---- 30-node tillering benchmark. FROZEN: these four entries reproduce
    # 02_networks/Tillering_Sorghum_Benchmark_suc byte-for-byte. Do not reword.
    "flashp_suc": [
        {"s": "Sucrose", "t": "MAX2", "x": -1, **_MAX2},
        {"s": "Sucrose", "t": "SMXL678", "x": 1, **_SMXL},
        {"s": "Sucrose", "t": "Tiller_Number", "x": 1, **_SUPPLY,
         "note": ("SORGHUM-DIRECT - the strongest species grounding available. Represents the "
                  "direct metabolic-supply role of sucrose in SUSTAINING outgrowth, which the "
                  "hand-curated model encodes by putting SUC inside the phenotype formula: "
                  "Sustained_growth = ((SUC + min(GA,PAT_Bud))/2) * Bud_release. "
                  "This is the single largest missing contributor to sucrose response amplitude.")},
        {"s": "Sucrose", "t": "Cytokinin", "x": 1, **_CK,
         "note": ("NOT grass-direct - dicot evidence. Flagged. Included because it also "
                  "converts Cytokinin from an unregulated source node into a regulated one, "
                  "activating the otherwise inert CK arm (Cytokinin -| BRC1, -> Auxin_Bud).")},
    ],

    # ---- 69-node shoot-branching network.
    "flashp_shoot": [
        # The supply edge targets Sustained_growth, NOT Shoot_Branching. This network
        # already carries Sustained_growth -> Shoot_Branching, and the hand-curated model
        # puts SUC literally inside that node:
        #     Sustained_growth = ((SUC + min(GA,PAT_Bud))/2) * Bud_release
        # so this reproduces the curated structure instead of inventing a shortcut edge
        # straight to the phenotype.
        {"s": "SUC", "t": "Sustained_growth", "x": 1, **_SUPPLY,
         "note": ("SORGHUM-DIRECT - the strongest species grounding available. Represents the "
                  "direct metabolic-supply role of sucrose in SUSTAINING outgrowth. Placed on "
                  "Sustained_growth rather than Shoot_Branching because this network already "
                  "has Sustained_growth -> Shoot_Branching, and the hand-curated model encodes "
                  "supply in exactly that node: "
                  "Sustained_growth = ((SUC + min(GA,PAT_Bud))/2) * Bud_release. "
                  "This is the single largest missing contributor to sucrose response amplitude.")},
        # MAX2 and SMXL6_7_8 are the nodes the max2 / smxl678 genotypes actually perturb
        # (gene_map.py variant psoup_style). The network also carries a parallel grass arm
        # (D3 -| D53 -| TB1/IPA1) that no genotype touches -- see OPTIONAL_GRASS_TWINS.
        {"s": "SUC", "t": "MAX2", "x": -1, **_MAX2},
        {"s": "SUC", "t": "SMXL6_7_8", "x": 1, **_SMXL},
        {"s": "SUC", "t": "CK", "x": 1, **_CK,
         "note": ("NOT grass-direct - dicot evidence. Flagged. Unlike the benchmark network, "
                  "CK here is already regulated (IPT ->, LOG ->, CKX -|), so this edge adds a "
                  "route without changing CK's source status.")},
    ],
}

# Optional, OFF by default: the literal rice genes in Patil et al. 2022. The shoot
# network carries the grass orthologue arm (D3 -| D53) alongside the Arabidopsis-named
# one (MAX2 -| SMXL6_7_8), so the same evidence can be placed on either. No genotype
# perturbs D3 or D53, so these edges raise WT amplitude WITHOUT letting the mutants
# differentiate -- which is why they are not the default. Enable with --with-grass-twins
# to see what they cost; the script reports the amplitude either way.
OPTIONAL_GRASS_TWINS: dict[str, list[dict]] = {
    "flashp_shoot": [
        {"s": "SUC", "t": "D3", "x": -1, **_MAX2,
         "note": ("Grass-orthologue twin of SUC -| MAX2, on the node the Patil 2022 paper "
                  "literally names. Off by default: no genotype perturbs D3.")},
        {"s": "SUC", "t": "D53", "x": 1, **_SMXL,
         "note": ("Grass-orthologue twin of SUC -> SMXL6_7_8, on the node the Patil 2022 paper "
                  "literally names. Off by default: no genotype perturbs D53.")},
    ],
}

REASON = ("FlashP wired sucrose with too few outgoing edges; the hand-curated PSoup "
          "benchmark has 5. Added the missing sucrose->SL-signalling and "
          "sucrose->cytokinin arms, each with a literature DOI.")

# The benchmark arm's reason string predates the second arm and is kept verbatim so
# Tillering_Sorghum_Benchmark_suc regenerates byte-for-byte.
REASON_BY_ARM = {
    "flashp_suc": ("FlashP wired sucrose with 2 outgoing edges; the hand-curated PSoup "
                   "benchmark has 5. Added the missing sucrose->SL-signalling and "
                   "sucrose->cytokinin arms, each with a literature DOI."),
}

FLOOR = 0.01
EPS = 0.1
K = 10.0


def formula(node: str, act: list[str], inh: list[str], src: bool) -> str:
    """Rebuild a FlashP algebraic formula string for a node."""
    if src and not act and not inh:
        return f"{node} = gene_modifier + exogenous_supply"
    parts = []
    if act:
        terms = " * ".join(f"max({a}, {FLOOR})" for a in act)
        parts.append(f"({terms})^(1/{len(act)})" if len(act) > 1 else f"{terms}^(1/1)")
    if inh:
        prod = " * ".join(inh)
        parts.append(f"min(1/max({prod}, {EPS}), {K})")
    body = " * ".join(parts) if parts else "1.0"
    return f"{node} = {body} * gene_modifier + exogenous_supply"


def extend(a: Arm, with_twins: bool) -> None:
    src_dir, dst_dir = a.base_network_dir, a.network_dir
    if src_dir is None:
        raise SystemExit(f"{a.name}: no base_network_dir in arms.json -- cannot extend")
    if not src_dir.exists():
        raise SystemExit(f"{a.name}: missing base network {src_dir}")

    specs = list(EDGES.get(a.name, []))
    if with_twins:
        specs += OPTIONAL_GRASS_TWINS.get(a.name, [])
    if not specs:
        raise SystemExit(f"{a.name}: flagged extend_sucrose but no edges defined for it")

    print(f"\n{'='*70}\n{a.name}\n  from {src_dir.name}\n  to   {dst_dir.name}\n{'='*70}")

    # Idempotent: copy over the top rather than deleting first (cloud-synced
    # folders can hold locks on subfolders, so rmtree is unreliable). Both JSONs are
    # always regenerated from the pristine SRC copies below, so re-running is safe.
    shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)

    # ---------------- network.json (always rebuilt from SRC) ----------------
    npath = dst_dir / "network" / "network.json"
    net = json.loads((src_dir / "network" / "network.json").read_text(encoding="utf-8"))

    # Fail closed on a typo'd target. Without this the edge is written anyway, the
    # equation it should have joined is untouched, and the arm silently produces the
    # same flat sucrose response the extension was meant to fix -- the same failure
    # mode build_apsim_parameters.py:118-127 guards against for the sucrose node.
    known = {n["id"] for n in net["nodes"]}
    unknown = sorted({s["t"] for s in specs if s["t"] not in known} |
                     {s["s"] for s in specs if s["s"] not in known})
    if unknown:
        raise SystemExit(
            f"{a.name}: edge endpoint(s) not in {src_dir.name}: {', '.join(unknown)}")

    existing = {(e["s"], e["t"]) for e in net["edges"]}
    next_id = max(int(e["eid"][1:]) for e in net["edges"]) + 1
    added = []
    for spec in specs:
        if (spec["s"], spec["t"]) in existing:
            print(f"  skip (already present): {spec['s']} -> {spec['t']}")
            continue
        net["edges"].append({
            "s": spec["s"], "t": spec["t"], "x": spec["x"],
            "eid": f"N{next_id:03d}", "d": spec["d"],
        })
        added.append(spec)
        next_id += 1

    # Any node that gained an upstream regulator is no longer a source. Driven by
    # what actually changed rather than named per arm: Cytokinin in the benchmark,
    # nothing in the shoot network (its CK already had IPT / LOG / CKX upstream).
    gained = {s["t"] for s in added}
    flipped = sorted(n["id"] for n in net["nodes"] if n["id"] in gained and n.get("src"))
    for n in net["nodes"]:
        if n["id"] in flipped:
            n["src"] = False

    net["metadata"]["total_edges"] = len(net["edges"])
    net["metadata"]["local_extension"] = {
        "applied_by": "FlashP_APSIM_Reproduction/scripts/extend_sucrose.py",
        "reason": REASON_BY_ARM.get(a.name, REASON),
        "edges_added": [
            {k: v for k, v in s.items() if k != "x"} | {"sign": s["x"]} for s in added
        ],
        "also": (" ".join(f"{n} changed from source to regulated node." for n in flipped)
                 if flipped else "No source node changed."),
        "note": "FlashP_Repo is unmodified; this extension lives only in this project.",
    }
    npath.write_text(json.dumps(net, indent=1), encoding="utf-8")

    # ---------------- algebraic_equations.json (always rebuilt from SRC) ----------------
    epath = dst_dir / "network" / "algebraic_equations.json"
    eqs = json.loads((src_dir / "network" / "algebraic_equations.json").read_text(encoding="utf-8"))
    by_node = {e["n"]: e for e in eqs["equations"]}

    for spec in added:
        tgt = by_node[spec["t"]]
        bucket = "a" if spec["x"] > 0 else "inh"
        tgt.setdefault(bucket, [])
        if spec["s"] not in tgt[bucket]:
            tgt[bucket].append(spec["s"])
        if spec["t"] in flipped:
            tgt["src"] = False
        tgt["f"] = formula(tgt["n"], tgt.get("a") or [], tgt.get("inh") or [], tgt.get("src", False))

    eqs.setdefault("metadata", {})["local_extension"] = net["metadata"]["local_extension"]
    epath.write_text(json.dumps(eqs, indent=1), encoding="utf-8")

    # ---------------- report ----------------
    out_edges = [e for e in net["edges"] if e["s"] == a.sucrose_node]
    print(f"\nwrote: {dst_dir}")
    print(f"edges: {len(net['edges'])} (was {len(net['edges']) - len(added)})")
    print(f"outgoing {a.sucrose_node} edges: {len(out_edges)} "
          f"(hand-curated has 5)")
    for s in added:
        arrow = "->" if s["x"] > 0 else "-|"
        print(f"\n  {s['s']} {arrow} {s['t']}   [{s['species']}]  {s['d']}")
        print(f"     {s['evidence']}")
        if "NOT grass-direct" in s.get("note", ""):
            print(f"     !! {s['note']}")
    print("\nupdated equations:")
    for s in added:
        e = by_node[s["t"]]
        print(f"  [{e['n']}] src={e.get('src')} a={e.get('a')} inh={e.get('inh')}")
        print(f"     {e['f']}")


def main() -> None:
    argv = [x for x in sys.argv[1:] if not x.startswith("--")]
    with_twins = "--with-grass-twins" in sys.argv[1:]
    bad = [x for x in sys.argv[1:] if x.startswith("--") and x != "--with-grass-twins"]
    if bad:
        raise SystemExit(f"unknown flag(s): {', '.join(bad)}")

    picked = [a for a in select(argv) if a.extend_sucrose]
    if not picked:
        raise SystemExit("no selected arm is flagged extend_sucrose in arms.json")
    if with_twins:
        print("!! --with-grass-twins: adding SUC -| D3 and SUC -> D53. These raise WT")
        print("   amplitude but no genotype perturbs them, so they do NOT help the")
        print("   mutants differentiate. Compare the amplitudes before keeping this.")
    for a in picked:
        extend(a, with_twins)


if __name__ == "__main__":
    main()
