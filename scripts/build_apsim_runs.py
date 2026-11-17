"""Stage a self-contained APSIM run for each arm, patched with that arm's parameters.

Takes the PRISTINE templates (the .bak files, before the paper's own patcher touched
them) from 00_reference/apsim_inputs/ and rewrites exactly three things:

  1. every cultivar's tillerSdIntercept / tillerSdSlope, from
     03_runs/Genotype_APSIM_Parameters_<arm>.csv
  2. the Genotype factor, to the representatives we chose
  3. the weather FileName, to "Bullarah.sim". The weather file is copied next to
     each staged .apsimx, and APSIM resolves a bare file name relative to the
     .apsimx, so the staged runs hold no machine-specific absolute path.

Nothing else in the tree is touched -- Clock, Soil, SoilWater factors, Manager
scripts, reports and the other 36 cultivars are left exactly as they were. The
templates in 00_reference/ are opened read-only and never written.

Output: 03_runs/apsim/<arm>/*.apsimx, ready for Models.exe.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import select  # noqa: E402
PAPER = ROOT / "00_reference" / "apsim_inputs"
WEATHER = ROOT / "00_reference" / "apsim_inputs" / "Bullarah.sim"
OUT = ROOT / "03_runs" / "apsim"

TEMPLATES = ["Exploration_NoTiller.apsimx", "DensityManagement_NoTiller.apsimx"]

# Low / Medium / High. nced3 and max2 are the chosen pair; lbo is kept from the
# paper's own trio as the Medium so the figure scripts keep their 3-way structure.
REPRESENTATIVES = ["nced3", "lbo", "max2"]

INTERCEPT = "[Sorghum].Leaf.Parameters.tillerSdIntercept.FixedValue ="
SLOPE = "[Sorghum].Leaf.Parameters.tillerSdSlope.FixedValue ="

# Our CSV uses the supplementary-table names; the apsimx tree uses historical
# short names. Same mapping the paper's updater uses, so cultivars line up.
CSV_TO_APSIMX = {
    "brc12": "brc1", "brc12d14": "brc1d14", "brc12hb214053": "brc1hb21",
    "brc12max1": "brc1max1", "brc12max2": "brc1max2", "brc12max4": "brc1max4",
    "brc12pin1347": "brc1pin", "cxeoe": "CXEOE", "d14pin1347": "d14pin",
    "hb214053": "hb21", "pin1347": "pin", "pin1347max2": "pinmax2",
    "pin1347max4": "pinmax4", "smxl678": "smxl6", "smxl678max2": "smxl6max2",
    "smxl678max3": "smxl6max3", "smxl678spl915": "smxl6spl9", "spl915": "spl9",
    "spl915fhy3far1": "spl9fhy3", "fhy3far1": "fhy3",
}


def read_params(arm: str) -> dict[str, tuple[float, float]]:
    p = ROOT / "03_runs" / f"Genotype_APSIM_Parameters_{arm}.csv"
    out = {}
    for r in csv.DictReader(p.open(encoding="utf-8")):
        out[r["Genotype"].strip()] = (
            float(r["Intercept_Line"].rsplit("=", 1)[1]),
            float(r["Slope_Line"].rsplit("=", 1)[1]),
        )
    return out


def cultivars(tree: dict) -> list[dict]:
    for c in tree.get("Children", []):
        if c.get("Name", "").lower() == "replacements":
            for s in c.get("Children", []):
                if s.get("Name", "").lower() in ("cultivars", "cultivar"):
                    return s["Children"]
    raise RuntimeError("Replacements -> Cultivars not found")


def walk(node: dict):
    yield node
    for c in node.get("Children", []) or []:
        yield from walk(c)


def patch(tree: dict, params: dict[str, tuple[float, float]]) -> dict[str, int]:
    by_name = {c["Name"]: c for c in cultivars(tree)}
    counts = {"cultivars": 0, "missing": 0, "weather": 0, "factor": 0}

    for csv_name, (a, b) in params.items():
        node = by_name.get(csv_name) or by_name.get(CSV_TO_APSIMX.get(csv_name, ""))
        if node is None:
            print(f"    WARNING no apsimx cultivar for {csv_name!r}")
            counts["missing"] += 1
            continue
        node["Command"] = [f"{INTERCEPT} {a:.3f}", f"{SLOPE} {b:.3f}"]
        node["Name"] = csv_name          # rename to the supplementary-table name
        counts["cultivars"] += 1

    for n in walk(tree):
        # weather: a bare file name, resolved by APSIM relative to the .apsimx.
        # main() copies Bullarah.sim next to every staged file.
        if "FileName" in n and isinstance(n["FileName"], str) and \
                n["FileName"].lower().endswith((".sim", ".met")):
            n["FileName"] = WEATHER.name
            counts["weather"] += 1
        # genotype factor: run our representatives
        if n.get("Name") == "Genotype" and "Specification" in n:
            n["Specification"] = ("[SowingRule].Script.CultivarName = "
                                  + ", ".join(REPRESENTATIVES))
            counts["factor"] += 1
    return counts


def main() -> None:
    if not WEATHER.exists():
        raise SystemExit(f"weather file missing: {WEATHER}")

    # Arm-filtered on purpose: staging DELETES the .db beside each .apsimx it writes,
    # so an unfiltered run would discard another arm's completed APSIM output and
    # force a full re-run. run_all.ps1 passes the selected arms through.
    for a in select(sys.argv[1:]):
        arm = a.name
        params = read_params(arm)
        dest = OUT / arm
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(WEATHER, dest / WEATHER.name)
        print(f"\n=== {arm}  ({len(params)} genotypes in the parameter file) ===")

        for name in TEMPLATES:
            src = PAPER / (name + ".bak")
            if not src.exists():                       # no .bak -> use the live file
                src = PAPER / name
            tree = json.loads(src.read_text(encoding="utf-8"))
            n_before = len(cultivars(tree))
            counts = patch(tree, params)
            target = dest / name
            target.write_text(json.dumps(tree, indent=2) + "\n", encoding="utf-8")

            # a stale .db next to the apsimx would be appended to, not replaced
            db = target.with_suffix(".db")
            if db.exists():
                db.unlink()

            print(f"  {name}")
            print(f"    source        {src.name}")
            print(f"    cultivars     {counts['cultivars']}/{n_before} patched"
                  f"{', ' + str(counts['missing']) + ' missing' if counts['missing'] else ''}")
            print(f"    weather       {counts['weather']} node(s) -> Bullarah.sim (local)")
            print(f"    genotypes     {', '.join(REPRESENTATIVES)}"
                  f"  ({counts['factor']} factor node)")

    print(f"\nstaged under {OUT}")
    print("nothing in 00_reference was modified.")


if __name__ == "__main__":
    main()
