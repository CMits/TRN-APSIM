"""Column 1 control — reproduce the published PSoup values with PSoup's own simulator.

This is the foundation check for the whole comparison. Before we can attribute any
difference to the rule form or to the network, we must show the harness reproduces
the numbers the paper actually published.

Reads the genotype definitions straight out of DataSetS1_Results_APSIM_Suc.csv (it
carries all 26 `_n` modifier columns per row), re-runs psoup_simulator, and diffs
Sustained_growth against the published value.

Read-only. Writes 01_audit/col1_reproduction_check.csv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "00_reference" / "psoup_source"))

from psoup_simulator import simulate, GENE_NAMES  # noqa: E402

SRC = ROOT / "00_reference" / "paper_data" / "DataSetS1_Results_APSIM_Suc.csv"
OUT = ROOT / "01_audit" / "col1_reproduction_check.csv"


def main() -> None:
    rows = list(csv.DictReader(SRC.open(encoding="utf-8-sig")))
    print(f"published rows: {len(rows)}")

    out_rows = []
    worst = 0.0
    for r in rows:
        # the CSV carries the full genotype: one column per `_n` modifier
        geno = {g: float(r[g]) for g in GENE_NAMES}
        published = float(r["Sustained_growth"])
        recomputed = simulate(genotype=geno)["Sustained_growth"]
        delta = abs(recomputed - published)
        worst = max(worst, delta)
        out_rows.append({
            "genotype": r["genotype"],
            "sucrose_label": r["sucrose_label"].strip().lower(),
            "SUC_n": r["SUC_n"],
            "published_Sustained_growth": f"{published:.9f}",
            "recomputed_Sustained_growth": f"{recomputed:.9f}",
            "abs_delta": f"{delta:.9f}",
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)

    # ---- report ----
    deltas = sorted(out_rows, key=lambda x: -float(x["abs_delta"]))
    print(f"max |delta| vs published : {worst:.6f}")
    n_close = sum(1 for r in out_rows if float(r["abs_delta"]) < 1e-6)
    print(f"rows matching to 1e-6    : {n_close}/{len(out_rows)}")
    print("\nworst 8 rows:")
    print(f"  {'genotype':<34}{'suc':<8}{'published':>12}{'recomputed':>12}{'delta':>10}")
    for r in deltas[:8]:
        print(f"  {r['genotype']:<34}{r['sucrose_label']:<8}"
              f"{float(r['published_Sustained_growth']):>12.4f}"
              f"{float(r['recomputed_Sustained_growth']):>12.4f}"
              f"{float(r['abs_delta']):>10.4f}")
    print(f"\nwrote: {OUT}")

    if worst < 1e-6:
        print("\nVERDICT: exact reproduction — harness is faithful, use it as column 1.")
    else:
        print(f"\nVERDICT: does NOT reproduce exactly (max delta {worst:.4f}).")
        print("The published values came from R; this Python port is known to diverge")
        print("at high sucrose. Column 1 must therefore be the PUBLISHED csv values,")
        print("and this simulator used only for genotypes absent from it.")


if __name__ == "__main__":
    main()
