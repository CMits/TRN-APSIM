"""Stage 1 - gene coverage audit.

Answers the question the whole project hinges on: do the FlashP networks contain
the genes the paper perturbed, and therefore which of its 34 genotypes can be
reproduced at all?

Read-only. Writes to ../01_audit/:
  gene_coverage.csv      one row per (network, paper gene token)
  genotype_coverage.csv  one row per (network, paper genotype) - constructible y/n
  ortholog_map.json      the frozen, human-signable mapping used above

The ortholog map is deliberately explicit rather than inferred by string
similarity: BRC1_2 -> TB1 and SMXL6_7_8 -> D53 are real orthology claims, not
name matches, and they need to be reviewable.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER_DATA = ROOT / "00_reference" / "paper_data"
sys.path.insert(0, str(ROOT / "scripts"))
import gene_map  # noqa: E402
from arms import ARMS  # noqa: E402

NETWORKS = ROOT / "02_networks"
OUT = ROOT / "01_audit"

# ----------------------------------------------------------------------
# The paper's gene tokens -> node id per network.
# None means "no node and no accepted ortholog in this network".
# ----------------------------------------------------------------------
ORTHOLOG_MAP: dict[str, dict[str, str | None]] = {
    #  paper token        Tillering_Sorghum_Benchmark(_suc)
    "BRC1_2":       {"sorghum": "TB1",    "arabidopsis": "BRC1",    "note": "BRC1/TB1 TCP orthologues"},
    "D14":          {"sorghum": "D14",    "arabidopsis": "D14",     "note": "SL receptor, same symbol"},
    "MAX2":         {"sorghum": "SBMAX2", "arabidopsis": "MAX2",    "note": "SbMAX2 = rice D3"},
    "MAX3":         {"sorghum": "D17",    "arabidopsis": "CCD7",    "note": "MAX3 = CCD7 = D17"},
    "MAX4":         {"sorghum": "D10",    "arabidopsis": "CCD8",    "note": "MAX4 = CCD8 = D10"},
    "SMXL6_7_8":    {"sorghum": "D53",    "arabidopsis": "SMXL678", "note": "D53/SMXL6,7,8 repressor family"},
    "D27":          {"sorghum": "D27",    "arabidopsis": "D27",     "note": "same symbol"},
    "PIN_1_3_4_7":  {"sorghum": "PIN1",   "arabidopsis": "PIN1",    "note": "auxin efflux; paper pools PIN1/3/4/7"},
    "NCED3":        {"sorghum": "NCED1",  "arabidopsis": "NCED3",   "note": "NCED paralogues"},
    "SPL9_15":      {"sorghum": "IPA1",   "arabidopsis": "SPL9",    "note": "IPA1 = OsSPL14; SPL family"},
    "MAX1":         {"sorghum": None,     "arabidopsis": "MAX1",    "note": "SL biosynthesis tail"},
    "LBO":          {"sorghum": None,     "arabidopsis": "LBO",     "note": "SL biosynthesis tail"},
    "CLAMT":        {"sorghum": None,     "arabidopsis": None,      "note": "SL biosynthesis tail - absent from both"},
    "CXE":          {"sorghum": None,     "arabidopsis": None,      "note": "SL catabolism - absent from both"},
    "FHY3_FAR1":    {"sorghum": None,     "arabidopsis": "FHY3",    "note": "shade -> ABA arm"},
    "HB21_40_53":   {"sorghum": None,     "arabidopsis": "HB21",    "note": "shade -> ABA arm, HB21/40/53"},
}

# Networks do not share a node-naming convention: the benchmarks label nodes
# benchmark-style (BRC1, SMXL678, ...), Shoot_Branching_Sorghum PSoup-style
# (BRC1_2, SMXL6_7_8, ...). Both mappings already exist as the single source of
# truth in gene_map.VARIANTS, so derive every convention column from there rather
# than hand-writing columns that could drift.
for _tok, _m in ORTHOLOG_MAP.items():
    for _variant, _map in gene_map.VARIANTS.items():
        _m[_variant] = _map.get(f"{_tok}_n")

NETWORK_DIRS = {
    # The two benchmark networks are audited directly (the base one is not itself an
    # arm). The earlier Tillering_Sorghum_v1 / ShootBranching_Arabidopsis networks
    # this script used to audit now sit in _archive_unused/.
    #
    # Order matters: these come first so the existing rows of gene_coverage.csv and
    # genotype_coverage.csv stay unchanged as arms are added.
    "Tillering_Sorghum_Benchmark": ("benchmark", NETWORKS / "Tillering_Sorghum_Benchmark" / "network" / "network.json"),
    "Tillering_Sorghum_Benchmark_suc": ("benchmark", NETWORKS / "Tillering_Sorghum_Benchmark_suc" / "network" / "network.json"),
}

# ...then every registered arm that names an ortholog convention and is not already
# covered above.
for _a in ARMS.values():
    if _a.audit_key is None or _a.network_dir is None:
        continue
    if _a.network_dir.name in NETWORK_DIRS:
        continue
    NETWORK_DIRS[_a.network_dir.name] = (
        _a.audit_key, _a.network_dir / "network" / "network.json")


def load_nodes(path: Path) -> set[str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {n["id"] if isinstance(n, dict) else n for n in data["nodes"]}


def tokens_of(genotype: str) -> list[str]:
    """Split a paper genotype string into its gene tokens.

    Handles 'BRC1_2_n, MAX_1_n' -> ['BRC1_2', 'MAX1'], the MAX_n_n typo form,
    the bare 'SPL9_15' (no _n), and the 'CXE_n_OE' overexpressor.
    """
    out = []
    for raw in genotype.split(","):
        t = raw.strip()
        if not t or t == "WT":
            continue
        t = re.sub(r"_OE$", "", t)          # CXE_n_OE -> CXE_n
        t = re.sub(r"_n$", "", t)           # BRC1_2_n -> BRC1_2
        t = re.sub(r"^MAX_(\d)$", r"MAX\1", t)  # MAX_1 -> MAX1
        out.append(t)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    genotypes = [r["Genotype"] for r in
                 csv.DictReader((PAPER_DATA / "Genotype_Linear_Equations.csv").open(encoding="utf-8-sig"))]

    nodes = {name: load_nodes(p) for name, (_, p) in NETWORK_DIRS.items()}

    # ---------------- gene-level coverage ----------------
    gene_rows = []
    for net_name, (key, _) in NETWORK_DIRS.items():
        for token, m in ORTHOLOG_MAP.items():
            mapped = m[key]
            present = mapped is not None and mapped in nodes[net_name]
            gene_rows.append({
                "network": net_name,
                "paper_gene": token,
                "mapped_node": mapped or "",
                "present": "yes" if present else "no",
                "relation": "identical" if mapped == token else ("absent" if mapped is None else "ortholog"),
                "note": m["note"],
            })
            # a mapping that names a node the network doesn't have is a bug in this table
            if mapped is not None and mapped not in nodes[net_name]:
                raise SystemExit(
                    f"ortholog_map error: {net_name} (convention '{key}') has no "
                    f"node '{mapped}' for {token}")

    with (OUT / "gene_coverage.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(gene_rows[0]))
        w.writeheader()
        w.writerows(gene_rows)

    # ---------------- genotype-level coverage ----------------
    geno_rows = []
    for net_name, (key, _) in NETWORK_DIRS.items():
        for g in genotypes:
            toks = tokens_of(g)
            missing = [t for t in toks
                       if t not in ORTHOLOG_MAP or ORTHOLOG_MAP[t][key] is None]
            unknown = [t for t in toks if t not in ORTHOLOG_MAP]
            geno_rows.append({
                "network": net_name,
                "genotype": g,
                "tokens": " + ".join(toks) or "WT",
                "mapped_nodes": " + ".join(
                    ORTHOLOG_MAP[t][key] or "?" for t in toks if t in ORTHOLOG_MAP),
                "constructible": "yes" if not missing else "no",
                "missing_tokens": " + ".join(missing),
                "unrecognised_tokens": " + ".join(unknown),
            })

    with (OUT / "genotype_coverage.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(geno_rows[0]))
        w.writeheader()
        w.writerows(geno_rows)

    (OUT / "ortholog_map.json").write_text(
        json.dumps(ORTHOLOG_MAP, indent=2), encoding="utf-8")

    # ---------------- summary ----------------
    print(f"paper genotypes: {len(genotypes)}\n")
    for net_name, (key, _) in NETWORK_DIRS.items():
        present = [t for t, m in ORTHOLOG_MAP.items() if m[key]]
        absent = [t for t, m in ORTHOLOG_MAP.items() if not m[key]]
        ok = [r for r in geno_rows if r["network"] == net_name and r["constructible"] == "yes"]
        print(f"{net_name}")
        print(f"  nodes            : {len(nodes[net_name])}")
        print(f"  genes covered    : {len(present)}/{len(ORTHOLOG_MAP)}")
        print(f"  missing          : {', '.join(absent) or '-'}")
        print(f"  genotypes usable : {len(ok)}/{len(genotypes)}")
        print()

    # the three cultivars the headline figure is built from
    print("headline figure representatives (_figure_common.py GENOTYPES):")
    for label, token in (("max2 (High)", "MAX2"), ("lbo (Medium)", "LBO"), ("fhy3far1 (Low)", "FHY3_FAR1")):
        marks = {n: ("ok" if ORTHOLOG_MAP[token][k] else "MISSING")
                 for n, (k, _) in NETWORK_DIRS.items()}
        print(f"  {label:<16} " + "  ".join(f"{n}={v}" for n, v in marks.items()))

    print(f"\nwrote: {OUT / 'gene_coverage.csv'}")
    print(f"wrote: {OUT / 'genotype_coverage.csv'}")
    print(f"wrote: {OUT / 'ortholog_map.json'}")


if __name__ == "__main__":
    main()
