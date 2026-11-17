"""Step 3 — does the FlashP network reproduce the hand-curated mutant x sucrose patterns?

Compares three columns over the same 37 genotypes x 3 sucrose levels:

  A_published  hand-curated PSoup, the paper's own numbers (authoritative)
  A_canonical  hand-curated PSoup equations re-run here (psoup_simulator.py)
  B_flashp*    one column per FlashP arm (see scripts/arms.json), each that arm's
               network under FlashP's algebraic engine

Reports, per pair:
  - rank agreement across genotypes at each sucrose level (Spearman rho, Kendall tau)
  - sucrose-response shape per genotype, normalised to its own medium level
  - direction agreement (does sucrose move each genotype the same way?)

Rank statistics are implemented directly so the script has no scipy dependency.

Writes 01_audit/pattern_comparison.csv and 01_audit/pattern_summary.txt
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "00_reference" / "psoup_source"))
from psoup_simulator import simulate, GENE_NAMES  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from arms import flashp_arms  # noqa: E402

PAPER = ROOT / "00_reference" / "paper_data" / "DataSetS1_Results_APSIM_Suc.csv"
OUT_CSV = ROOT / "01_audit" / "pattern_comparison.csv"
OUT_TXT = ROOT / "01_audit" / "pattern_summary.txt"

LEVELS = ["low", "medium", "high"]


# ----------------------------- rank statistics -----------------------------

def _ranks(xs: list[float]) -> list[float]:
    """Average ranks, ties shared."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(a: list[float], b: list[float]) -> float:
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = sum((x - ma) ** 2 for x in a) ** 0.5
    db = sum((y - mb) ** 2 for y in b) ** 0.5
    return num / (da * db) if da and db else float("nan")


def spearman(a: list[float], b: list[float]) -> float:
    return pearson(_ranks(a), _ranks(b))


def kendall_tau_b(a: list[float], b: list[float]) -> float:
    n = len(a)
    conc = disc = ta = tb = 0
    for i in range(n):
        for j in range(i + 1, n):
            da, db = a[i] - a[j], b[i] - b[j]
            s = (da > 0) - (da < 0)
            t = (db > 0) - (db < 0)
            if s == 0 and t == 0:
                continue
            if s == 0:
                ta += 1
            elif t == 0:
                tb += 1
            elif s * t > 0:
                conc += 1
            else:
                disc += 1
    denom = ((conc + disc + ta) * (conc + disc + tb)) ** 0.5
    return (conc - disc) / denom if denom else float("nan")


# ----------------------------------- data -----------------------------------

def flashp_columns() -> list[tuple[str, object]]:
    """(column name, arm) for every FlashP arm whose sweep has actually been run.

    Skipping arms with no response_csv keeps `run_all.ps1 -Arm <x>` from demanding
    another arm's sweep output.
    """
    return [(a.pattern_column, a) for a in flashp_arms() if a.response_csv.exists()]


def load_columns(fcols) -> dict[str, dict[tuple[str, str], float]]:
    cols: dict[str, dict[tuple[str, str], float]] = {"A_published": {}, "A_canonical": {}}
    for name, _ in fcols:
        cols[name] = {}

    for r in csv.DictReader(PAPER.open(encoding="utf-8-sig")):
        key = (r["genotype"], r["sucrose_label"].strip().lower())
        cols["A_published"][key] = float(r["Sustained_growth"])
        geno = {g: float(r[g]) for g in GENE_NAMES}
        cols["A_canonical"][key] = simulate(genotype=geno)["Sustained_growth"]

    for name, a in fcols:
        for r in csv.DictReader(a.response_csv.open(encoding="utf-8-sig")):
            cols[name][(r["genotype"], r["sucrose_label"])] = float(r[a.phenotype])

    return cols


def main() -> None:
    fcols = flashp_columns()
    if not fcols:
        raise SystemExit("no FlashP sweep output found -- run run_sucrose_sweep.py first")
    FLASHP_COLS = [name for name, _ in fcols]
    ALL_COLS = ["A_published", "A_canonical", *FLASHP_COLS]
    cols = load_columns(fcols)
    genotypes = sorted({g for g, _ in cols["A_published"]})
    lines: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        lines.append(s)

    say(f"genotypes: {len(genotypes)}    sucrose levels: {LEVELS}")
    say()
    say("REFERENCE = A_canonical (hand-curated PSoup equations, re-run here).")
    say("A_published is retained only as provenance, not as the comparison target:")
    say("  - no simulator on disk reproduces its magnitudes (max delta 2.28, and the")
    say("    sucrose-extension variant is worse), so it cannot be regenerated;")
    say("  - 5 PIN genotypes are published as exactly 0.0, an artefact the paper's own")
    say("    run_psoup_34mutants_PIN075.py exists to work around;")
    say("  - it costs nothing to drop: published vs canonical ranks at rho 0.97-0.999,")
    say("    so the mutant ORDERING is the same either way.")
    say()

    # ---------- 1. rank agreement at each sucrose level ----------
    say("RANK AGREEMENT ACROSS GENOTYPES (does the model order mutants the same way?)")
    say(f"  {'pair':<34}{'level':<9}{'Spearman':>10}{'Kendall':>10}")
    pairs = ([("A_canonical", c) for c in FLASHP_COLS]
             + [("A_published", "A_canonical")]
             + [("A_published", c) for c in FLASHP_COLS])
    for x, y in pairs:
        for lv in LEVELS:
            a = [cols[x][(g, lv)] for g in genotypes]
            b = [cols[y][(g, lv)] for g in genotypes]
            say(f"  {x+' vs '+y:<34}{lv:<9}{spearman(a,b):>10.3f}{kendall_tau_b(a,b):>10.3f}")
        say()

    # ---------- 2. sucrose-response shape, normalised to medium ----------
    say("SUCROSE RESPONSE SHAPE  (value / own medium)  — WT and extremes")
    say(f"  {'genotype':<30}{'model':<14}{'low':>8}{'med':>8}{'high':>8}")
    for g in ["WT", "BRC1_2_n", "MAX2_n", "LBO_n", "FHY3_FAR1_n"]:
        if g not in genotypes:
            continue
        for c in ALL_COLS:
            med = cols[c][(g, "medium")]
            vals = [cols[c][(g, lv)] / med if med else float("nan") for lv in LEVELS]
            say(f"  {g:<30}{c:<14}" + "".join(f"{v:>8.3f}" for v in vals))
        say()

    # ---------- 3. direction agreement ----------
    say("SUCROSE DIRECTION AGREEMENT (sign of high - low, per genotype)")
    def direction(c: str, g: str) -> int:
        d = cols[c][(g, "high")] - cols[c][(g, "low")]
        return (d > 1e-9) - (d < -1e-9)
    for x, y in pairs:
        agree = sum(1 for g in genotypes if direction(x, g) == direction(y, g))
        say(f"  {x} vs {y}: {agree}/{len(genotypes)} genotypes move the same way")
    say()

    # flat genotypes (no sucrose response) per model
    for c in ["A_canonical", *FLASHP_COLS, "A_published"]:
        flat = [g for g in genotypes if direction(c, g) == 0]
        say(f"  {c}: {len(flat)} genotypes are sucrose-unresponsive" +
            (f" -> {flat}" if flat and len(flat) <= 8 else ""))
    say()

    # ---------- 4. response-shape diversity ----------
    say("RESPONSE-SHAPE DIVERSITY (distinct normalised low/med/high curves, 3dp)")
    for c in ["A_canonical", *FLASHP_COLS]:
        shapes: dict[tuple, list[str]] = {}
        for g in genotypes:
            med = cols[c][(g, "medium")]
            key = tuple(round(cols[c][(g, lv)] / med, 3) if med else float("nan")
                        for lv in LEVELS)
            shapes.setdefault(key, []).append(g)
        say(f"  {c}: {len(shapes)} distinct shapes across {len(genotypes)} genotypes")
        biggest = max(shapes.items(), key=lambda kv: len(kv[1]))
        if len(biggest[1]) > 1:
            say(f"     largest cluster {biggest[0]} shared by {len(biggest[1])} genotypes")
    say()

    # ---------- write per-genotype table ----------
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["genotype", "sucrose_label", *ALL_COLS,
                    *[f"{c}_norm" for c in ALL_COLS]])
        for g in genotypes:
            for lv in LEVELS:
                row = [g, lv]
                raw = [cols[c][(g, lv)] for c in ALL_COLS]
                norm = []
                for c in ALL_COLS:
                    med = cols[c][(g, "medium")]
                    norm.append(cols[c][(g, lv)] / med if med else "")
                w.writerow(row + [f"{v:.6f}" for v in raw] +
                           [f"{v:.6f}" if v != "" else "" for v in norm])

    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")
    say(f"wrote: {OUT_CSV}")
    say(f"wrote: {OUT_TXT}")


if __name__ == "__main__":
    main()
