"""Which update rule reproduces the hand-curated mutant ordering?

Usage:  python compare_solvers.py [arm ...]        (default: every FlashP arm)

Holds the network TOPOLOGY constant within an arm and varies only the node update
rule, comparing each against the hand-curated PSoup model:

  flashp_algebraic   geometric mean activators; inhibition = min(1/max(prod,0.1),10)
  flashp_ode         Hill functions, K and n read from the network's own metadata
  psoup_form         arithmetic mean activators; inhibition = 2/(1+mean) -- no clamp

Writes each arm's solver_out_stem .csv and .txt
(01_audit/solver_comparison.* for flashp_suc).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "00_reference" / "psoup_source"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "00_reference" / "flashp_source"))

from psoup_simulator import simulate as psoup_sim, GENE_NAMES          # noqa: E402
from psoup_form_sim import PSoupFormSimulator                          # noqa: E402
from compare_patterns import spearman, kendall_tau_b                   # noqa: E402
from gene_map import gene_to_node                                      # noqa: E402
from arms import Arm, select                                           # noqa: E402
from flashp_validator import (FlashPSimulator, SimulationConfig,       # noqa: E402
                              load_equations as load_alg)
from ode_validator import (ODESimulator, ODEConfig,                    # noqa: E402
                           load_equations as load_ode)

PAPER = ROOT / "00_reference" / "paper_data" / "DataSetS1_Results_APSIM_Suc.csv"
LEVELS = ["low", "medium", "high"]


def ode_tuning(path: Path) -> tuple[float, int]:
    """K and n as tuned for THIS network, from its own equations metadata.

    The benchmark is K=0.5 n=2; Shoot_Branching_Sorghum is K=1.0 n=3. Reading them
    rather than hardcoding keeps each arm on its own tuned setting.
    """
    meta = json.loads(path.read_text(encoding="utf-8")).get("metadata", {})
    params = {**meta, **meta.get("parameters", {})}
    return float(params.get("K", 0.5)), int(params.get("n", 2))


def build_solvers(net: Path, pheno: str):
    K, n = ode_tuning(net / "ode_equations.json")
    alg = FlashPSimulator(load_alg(str(net / "algebraic_equations.json")), SimulationConfig())
    ode = ODESimulator(load_ode(str(net / "ode_equations.json")), ODEConfig(K=K, n=n))
    psf = PSoupFormSimulator(net / "algebraic_equations.json")
    return {
        "flashp_algebraic": lambda m: alg.simulate(m, {})[0][pheno],
        "flashp_ode":       lambda m: ode.simulate(m, {})[0][pheno],
        "psoup_form":       lambda m: psf.simulate(m)[0][pheno],
    }, K, n


def compare(arm: Arm) -> None:
    NET = arm.sweep_network_dir / "network"
    PHENO = arm.phenotype
    GENE_TO_NODE = gene_to_node(arm.gene_map)
    solvers, K, n = build_solvers(NET, PHENO)
    rows = list(csv.DictReader(PAPER.open(encoding="utf-8-sig")))

    ref: dict[tuple[str, str], float] = {}
    got: dict[str, dict[tuple[str, str], float]] = {k: {} for k in solvers}

    for r in rows:
        key = (r["genotype"], r["sucrose_label"].strip().lower())
        ref[key] = psoup_sim(genotype={g: float(r[g]) for g in GENE_NAMES})["Sustained_growth"]
        mods = {node: float(r[col]) for col, node in GENE_TO_NODE.items()
                if col in r and r[col] != "" and float(r[col]) != 1.0}
        mods[arm.sucrose_node] = float(r["SUC_n"])
        for name, fn in solvers.items():
            got[name][key] = fn(mods)

    genos = sorted({g for g, _ in ref})
    lines: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        lines.append(s)

    say("REFERENCE: hand-curated PSoup (psoup_simulator.py), same 37 genotypes.")
    say(f"Topology held constant = {arm.sweep_network_dir.name} (phenotype {PHENO}).")
    say(f"Only the rule varies. ODE tuning for this network: K={K}, n={n}.")
    say()
    say(f"  {'rule':<20}{'level':<9}{'Spearman':>10}{'Kendall':>10}")
    results = []
    for name in solvers:
        for lv in LEVELS:
            a = [ref[(g, lv)] for g in genos]
            b = [got[name][(g, lv)] for g in genos]
            rho, tau = spearman(a, b), kendall_tau_b(a, b)
            results.append({"rule": name, "level": lv,
                            "spearman": f"{rho:.4f}", "kendall": f"{tau:.4f}"})
            say(f"  {name:<20}{lv:<9}{rho:>10.3f}{tau:>10.3f}")
        say()

    say("RESPONSE-SHAPE DIVERSITY (distinct normalised curves / 37 genotypes)")
    for name in solvers:
        shapes: dict[tuple, int] = {}
        for g in genos:
            m = got[name][(g, "medium")]
            k = tuple(round(got[name][(g, lv)] / m, 3) if m else 0.0 for lv in LEVELS)
            shapes[k] = shapes.get(k, 0) + 1
        say(f"  {name:<20}{len(shapes):>3} shapes, largest cluster {max(shapes.values())}")
    hand = {}
    for g in genos:
        m = ref[(g, "medium")]
        k = tuple(round(ref[(g, lv)] / m, 3) if m else 0.0 for lv in LEVELS)
        hand[k] = hand.get(k, 0) + 1
    say(f"  {'hand-curated':<20}{len(hand):>3} shapes, largest cluster {max(hand.values())}")
    say()

    say("WT SUCROSE RESPONSE (normalised to own medium)")
    for name in solvers:
        w = [got[name][("WT", lv)] for lv in LEVELS]
        say(f"  {name:<20}{w[0]/w[1]:>8.3f}{1.0:>8.3f}{w[2]/w[1]:>8.3f}")
    w = [ref[("WT", lv)] for lv in LEVELS]
    say(f"  {'hand-curated':<20}{w[0]/w[1]:>8.3f}{1.0:>8.3f}{w[2]/w[1]:>8.3f}")

    out = arm.solver_out_stem.with_suffix(".csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        wtr = csv.DictWriter(f, fieldnames=["rule", "level", "spearman", "kendall"])
        wtr.writeheader()
        wtr.writerows(results)
    arm.solver_out_stem.with_suffix(".txt").write_text("\n".join(lines), encoding="utf-8")
    say(f"\nwrote: {out}")


def main() -> None:
    for a in select(sys.argv[1:], kind="python_flashp"):
        compare(a)


if __name__ == "__main__":
    main()
