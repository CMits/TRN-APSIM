"""Step 4/5 — derive APSIM tillerSd parameters, and pick the new Low/High cultivars.

Ports the paper's math (Pipeline.ipynb cells 12/16/20) and VERIFIES the port by
reproducing the published Genotype_APSIM_Parameters.csv to machine precision before
applying it to anything new.

  b_hat[level]  = anchor_medium * response[level] / response[medium]
  slope_raw     = (b_hat_high - b_hat_low) / (SD_max - SD_min)
  intercept_raw = b_hat_low - slope_raw * SD_min
  intercept     = minmax(intercept_raw) -> [-1.363, 1.5]
  slope         = minmax(slope_raw)     -> [0, 0.4]

Note there is NO cross-genotype regression: each genotype is its own two-point line.
The rescale uses the observed min/max across the genotype set being scaled, so the
set matters -- dropping fhy3far1 moves the intercept floor, since fhy3far1 defines it
in the published set.

Writes 01_audit/apsim_params_<arm>.csv and 01_audit/cultivar_choice.txt
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "00_reference" / "psoup_source"))
sys.path.insert(0, str(ROOT / "scripts"))

from psoup_simulator import simulate as psoup_sim, GENE_NAMES   # noqa: E402
from psoup_form_sim import PSoupFormSimulator                   # noqa: E402
from gene_map import GENE_TO_NODE                               # noqa: E402

PD = ROOT / "00_reference" / "paper_data"
LEVELS = ["low", "medium", "high"]
INT_BAND = (-1.363, 1.5)
SLOPE_BAND = (0.0, 0.4)

# genotypes excluded: FHY3_FAR1 has no confirmed grass orthologue (FlashP's own
# crosswalk returns "unconfirmed"), so it cannot be tested in both models.
EXCLUDE_TOKEN = "FHY3_FAR1"


def minmax(v: float, lo_obs: float, hi_obs: float, lo: float, hi: float) -> float:
    if hi_obs == lo_obs:
        return lo
    return lo + (v - lo_obs) / (hi_obs - lo_obs) * (hi - lo)


def derive(bhat: dict[str, dict[str, float]], sd_min: float, sd_max: float
           ) -> dict[str, dict[str, float]]:
    """b_hat per genotype -> raw and rescaled APSIM parameters."""
    raw = {}
    for g, b in bhat.items():
        s = (b["high"] - b["low"]) / (sd_max - sd_min)
        raw[g] = {"slope_raw": s, "intercept_raw": b["low"] - s * sd_min}
    irs = [r["intercept_raw"] for r in raw.values()]
    srs = [r["slope_raw"] for r in raw.values()]
    for g, r in raw.items():
        r["intercept"] = minmax(r["intercept_raw"], min(irs), max(irs), *INT_BAND)
        r["slope"] = minmax(r["slope_raw"], min(srs), max(srs), *SLOPE_BAND)
    return raw


def load_paper_bhat() -> dict[str, dict[str, float]]:
    out = {}
    for r in csv.DictReader((PD / "weighted_arabidopsis_branches.csv").open(encoding="utf-8-sig")):
        out[r["genotype"]] = {lv: float(r[lv]) for lv in LEVELS}
    return out


def normalise_genotype(name: str) -> str:
    """Same normalisation Pipeline.ipynb cell 4 applies: MAX_4_n -> MAX4_n, tidy commas.

    DataSetS1 carries both spellings for some doubles; the paper collapses them and
    averages, so we must too or the genotype keys will not line up.
    """
    s = re.sub(r"MAX_(\d)_n", r"MAX\1_n", name)
    return re.sub(r"\s*,\s*", ", ", s).strip()


def flashp_response(eq_path: Path) -> dict[str, dict[str, float]]:
    sim = PSoupFormSimulator(eq_path)
    acc: dict[str, dict[str, list[float]]] = {}
    for r in csv.DictReader((PD / "DataSetS1_Results_APSIM_Suc.csv").open(encoding="utf-8-sig")):
        mods = {n: float(r[c]) for c, n in GENE_TO_NODE.items()
                if c in r and r[c] != "" and float(r[c]) != 1.0}
        mods["Sucrose"] = float(r["SUC_n"])
        v, _, _ = sim.simulate(mods)
        g = normalise_genotype(r["genotype"])
        acc.setdefault(g, {}).setdefault(r["sucrose_label"].strip().lower(), []).append(
            v["Tiller_Number"])
    # average across duplicate spellings, as the notebook does
    return {g: {lv: sum(xs) / len(xs) for lv, xs in d.items()} for g, d in acc.items()}


def main() -> None:
    lin = list(csv.DictReader((PD / "Genotype_Linear_Equations.csv").open(encoding="utf-8-sig")))
    sd_min, sd_max = float(lin[0]["SD_min"]), float(lin[0]["SD_max"])
    published = {r["Genotype"]: (float(r["Intercept"]), float(r["Slope"])) for r in lin}

    lines: list[str] = []
    def say(s: str = "") -> None:
        print(s); lines.append(s)

    # ---------- 1. validate the port against the published parameters ----------
    paper_bhat = {g: b for g, b in load_paper_bhat().items() if g in published}
    check = derive(paper_bhat, sd_min, sd_max)
    worst = max(max(abs(check[g]["intercept"] - published[g][0]),
                    abs(check[g]["slope"] - published[g][1])) for g in published)
    say(f"PORT VALIDATION: max deviation from published parameters = {worst:.2e}")
    if worst > 1e-9:
        raise SystemExit("port does not reproduce the paper — fix before proceeding")
    say("  -> exact. The notebook math is correctly reproduced.\n")

    # ---------- 2. anchor + FlashP sucrose shape ----------
    anchor = {g: b["medium"] for g, b in paper_bhat.items()}
    arms = {
        "flashp_base": ROOT / "02_networks" / "Tillering_Sorghum_Benchmark" / "network" / "algebraic_equations.json",
        "flashp_suc":  ROOT / "02_networks" / "Tillering_Sorghum_Benchmark_suc" / "network" / "algebraic_equations.json",
    }

    eligible = [g for g in published if EXCLUDE_TOKEN not in g]
    say(f"eligible genotypes: {len(eligible)}/{len(published)}  "
        f"(excluded {len(published)-len(eligible)} containing {EXCLUDE_TOKEN})")
    say(f"SD range: [{sd_min:.4f}, {sd_max:.4f}]\n")

    derived: dict[str, dict[str, dict[str, float]]] = {}
    for arm, path in arms.items():
        resp = flashp_response(path)
        bhat = {g: {lv: anchor[g] * resp[g][lv] / resp[g]["medium"] for lv in LEVELS}
                for g in eligible}
        derived[arm] = derive(bhat, sd_min, sd_max)
        with (ROOT / "01_audit" / f"apsim_params_{arm}.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["Genotype", "Intercept", "Slope", "Intercept_raw", "Slope_raw"])
            for g in sorted(derived[arm]):
                d = derived[arm][g]
                w.writerow([g, f"{d['intercept']:.4f}", f"{d['slope']:.4f}",
                            f"{d['intercept_raw']:.6f}", f"{d['slope_raw']:.6f}"])

    # hand-curated arm on the same eligible set, for a like-for-like ranking
    hand = derive({g: paper_bhat[g] for g in eligible}, sd_min, sd_max)

    # ---------- 3. pick Low / Medium / High from the parameters ----------
    say("CULTIVAR SELECTION — ranked by derived tillerSdIntercept")
    say("(both arms must agree on the ordering for a genotype to be usable)\n")
    arm = "flashp_suc"
    rank_f = sorted(eligible, key=lambda g: derived[arm][g]["intercept"])
    rank_h = sorted(eligible, key=lambda g: hand[g]["intercept"])
    pos_f = {g: i for i, g in enumerate(rank_f)}
    pos_h = {g: i for i, g in enumerate(rank_h)}

    say(f"  {'rank':<6}{'FlashP (suc)':<32}{'int':>8}{'slope':>8}   {'hand-curated':<32}{'int':>8}")
    for i in range(len(eligible)):
        gf, gh = rank_f[i], rank_h[i]
        say(f"  {i+1:<6}{gf:<32}{derived[arm][gf]['intercept']:>8.3f}"
            f"{derived[arm][gf]['slope']:>8.3f}   {gh:<32}{hand[gh]['intercept']:>8.3f}")

    say()
    # candidates where the two arms broadly agree (rank gap small)
    agree = [g for g in eligible if abs(pos_f[g] - pos_h[g]) <= 6]
    lo_c = [g for g in rank_f if g in agree][:3]
    hi_c = [g for g in rank_f if g in agree][-3:]
    mid_i = len(rank_f) // 2
    md_c = [g for g in rank_f[mid_i - 3:mid_i + 3] if g in agree][:3]

    say(f"genotypes where both arms agree on rank (within 6 places): {len(agree)}/{len(eligible)}")
    say(f"\n  LOW  candidates : {lo_c}")
    say(f"  MED  candidates : {md_c}")
    say(f"  HIGH candidates : {hi_c}")
    if lo_c and hi_c:
        say(f"\n  proposed trio: LOW={lo_c[0]}  MED={md_c[0] if md_c else '?'}  HIGH={hi_c[-1]}")
        for g in (lo_c[0], md_c[0] if md_c else None, hi_c[-1]):
            if g:
                say(f"     {g:<32} FlashP int={derived[arm][g]['intercept']:+.3f} "
                    f"slope={derived[arm][g]['slope']:.3f} | hand int={hand[g]['intercept']:+.3f}")

    (ROOT / "01_audit" / "cultivar_choice.txt").write_text("\n".join(lines), encoding="utf-8")
    say(f"\nwrote: 01_audit/apsim_params_*.csv and cultivar_choice.txt")


if __name__ == "__main__":
    main()
