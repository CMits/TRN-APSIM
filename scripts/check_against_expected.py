"""Compare a fresh run against the results shipped in expected_results/.

Checks the files that carry the study's numbers:

  03_runs/Genotype_APSIM_Parameters_<arm>.csv   tillerSdIntercept / tillerSdSlope per cultivar
  03_runs/Arm_Comparison_long.csv               parameters, TTN and rank per arm x cultivar

Numbers are compared with a small absolute tolerance (the files are written to 3
decimals); text cells must match exactly. Exits non-zero on any difference.

  python scripts/check_against_expected.py
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS  # noqa: E402

EXPECTED = ROOT / "expected_results"
TOL = 1.5e-3          # one unit in the third decimal, plus rounding slack
NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def cells_equal(a: str, b: str) -> bool:
    if a == b:
        return True
    na, nb = NUM.findall(a), NUM.findall(b)
    if not na or len(na) != len(nb) or NUM.sub("#", a) != NUM.sub("#", b):
        return False
    return all(abs(float(x) - float(y)) <= TOL for x, y in zip(na, nb))


def compare(name: str, key: list[str]) -> int:
    exp_p, got_p = EXPECTED / name, ROOT / "03_runs" / name
    if not got_p.exists():
        print(f"  MISSING  {got_p.relative_to(ROOT)}")
        return 1
    exp = {tuple(r[k] for k in key): r for r in rows(exp_p)}
    got = {tuple(r[k] for k in key): r for r in rows(got_p)}
    bad = 0
    for k in sorted(set(exp) | set(got)):
        if k not in got or k not in exp:
            print(f"  {name}: row {k} only in {'expected' if k in exp else 'new run'}")
            bad += 1
            continue
        for col, ev in exp[k].items():
            gv = got[k].get(col)
            if gv is None or not cells_equal(ev, gv):
                print(f"  {name}: {k} {col}: expected {ev!r}, got {gv!r}")
                bad += 1
    print(f"  {'OK  ' if not bad else 'DIFF'}  {name}  ({len(exp)} rows)")
    return bad


def main() -> None:
    bad = 0
    for arm in ARMS:
        bad += compare(f"Genotype_APSIM_Parameters_{arm}.csv", ["Genotype"])
    bad += compare("Arm_Comparison_long.csv", ["arm", "apsim_name"])
    if bad:
        raise SystemExit(f"\n{bad} difference(s) against expected_results/")
    print("\nall results match expected_results/")


if __name__ == "__main__":
    main()
