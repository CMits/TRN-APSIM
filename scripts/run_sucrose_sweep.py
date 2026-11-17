"""Step 2 — sucrose sweep through a FlashP network.

Usage:  python run_sucrose_sweep.py [arm ...]      (default: every FlashP arm)

For each of the paper's genotypes x sucrose {0.5, 1.0, 1.5}, propagate through the
arm's network and record its phenotype node. The network, the phenotype, the sucrose
node and the gene-map variant all come from scripts/arms.json.

Sucrose is applied as a *gene_modifier*, not exogenous_supply: the paper drives it
as SUC_n in {0.5, 1.0, 1.5} multiplicatively, and every FlashP equation ends in
`* gene_modifier`, so the semantics are identical.

Genotype modifier VALUES are read from DataSetS1_Results_APSIM_Suc.csv rather than
re-derived, so arm B is perturbed exactly as arm A was. Only the node NAMES are
translated.

Fails closed: an unmapped genotype token raises rather than silently simulating
wild type.

Writes each arm's response_csv (01_audit/flashp_sucrose_response.csv for flashp_suc).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# Vendored into this project so the folder runs standalone (see 00_reference/flashp_source).
FLASHP_SHARED = ROOT / "00_reference" / "flashp_source"
sys.path.insert(0, str(FLASHP_SHARED))

from flashp_validator import FlashPSimulator, SimulationConfig, load_equations  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from gene_map import gene_to_node, no_node                                       # noqa: E402
from arms import Arm, select                                                     # noqa: E402

SRC = ROOT / "00_reference" / "paper_data" / "DataSetS1_Results_APSIM_Suc.csv"


def sweep(arm: Arm) -> None:
    gene_map = gene_to_node(arm.gene_map)
    unmapped = no_node(arm.gene_map)
    phenotype = arm.phenotype

    network = load_equations(str(arm.sweep_network_dir / "network" / "algebraic_equations.json"))
    sim = FlashPSimulator(network, SimulationConfig())

    rows = list(csv.DictReader(SRC.open(encoding="utf-8-sig")))
    out_rows = []
    unconverged = []

    for r in rows:
        genotype = r["genotype"]
        label = r["sucrose_label"].strip().lower()

        mods: dict[str, float] = {}
        for col, node in gene_map.items():
            if col not in r or r[col] == "":
                continue
            val = float(r[col])
            if val != 1.0:                     # only non-WT modifiers matter
                mods[node] = val
        # sucrose is always set explicitly, including the WT medium level
        mods[arm.sucrose_node] = float(r["SUC_n"])

        # fail closed: every non-WT paper modifier must map somewhere
        for col, val in r.items():
            if col.endswith("_n") and col not in gene_map and col not in unmapped:
                raise SystemExit(f"unmapped modifier column '{col}' — refusing to guess")

        values, converged, iters = sim.simulate(mods, {})
        if phenotype not in values:
            raise SystemExit(f"'{phenotype}' not produced. got: {sorted(values)[:12]}")
        if not converged:
            unconverged.append((genotype, label))

        out_rows.append({
            "genotype": genotype,
            "sucrose_label": label,
            "SUC_n": r["SUC_n"],
            phenotype: f"{values[phenotype]:.9f}",
            "converged": int(converged),
            "iterations": iters,
            "modifiers": ";".join(f"{k}={v}" for k, v in sorted(mods.items())),
        })

    out = arm.response_csv
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)

    genos = {r["genotype"] for r in out_rows}
    print(f"\n=== {arm.name} ({arm.sweep_network_dir.name}) ===")
    print(f"rows: {len(out_rows)}   genotypes: {len(genos)}")
    print(f"converged: {sum(int(r['converged']) for r in out_rows)}/{len(out_rows)}")
    if unconverged:
        print(f"  NOT converged: {unconverged[:6]}")

    wt = {r["sucrose_label"]: float(r[phenotype])
          for r in out_rows if r["genotype"] == "WT"}
    print(f"\nWT {phenotype}  low={wt.get('low'):.4f}  medium={wt.get('medium'):.4f}  "
          f"high={wt.get('high'):.4f}")
    print(f"\nwrote: {out}")


def main() -> None:
    for a in select(sys.argv[1:], kind="python_flashp"):
        sweep(a)


if __name__ == "__main__":
    main()
