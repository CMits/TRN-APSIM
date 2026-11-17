"""Compare the APSIM arms against each other and against the published run.

Question: does each network produce the same crop-level outcome as the hand-curated
one, once all have gone through the identical parameter chain?

Three checks:
  1. yield ordering of the three representatives -- must be Low < Medium < High in all
  2. per-cultivar yield agreement, each arm against the reference arm
  3. the reference arm vs the published .db, for the two cultivars they share
     (max2, lbo) -- a sanity check that our staging did not change the crop model

Writes 03_runs/APSIM_Arm_Comparison.txt
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from arms import ARMS as REGISTRY, REFERENCE  # noqa: E402

RUNS = ROOT / "03_runs" / "apsim"
PUB = ROOT / "00_reference" / "apsim_published"

# only arms whose APSIM run has actually produced output
ARMS = [a for a in REGISTRY
        if (RUNS / a / "DensityManagement_NoTiller.db").exists()]
FILES = ["Exploration_NoTiller.db", "DensityManagement_NoTiller.db"]
YIELD = '"Sorghum.Grain.Wt"'

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def stats(db: Path) -> dict[str, dict[str, float]]:
    c = sqlite3.connect(db)
    q = (f"select Genotype, count(*), avg({YIELD}), min({YIELD}), max({YIELD}), "
         f"avg(TTN), max(TTN), avg(FTN) from HarvestReport group by 1")
    return {r[0]: dict(n=r[1], mean=r[2], lo=r[3], hi=r[4],
                       ttn=r[5], ttn_max=r[6], ftn=r[7]) for r in c.execute(q)}


def main() -> None:
    if not ARMS:
        raise SystemExit(f"no APSIM output under {RUNS} -- run APSIM first")
    # fig_label, not label: both FlashP arms are called "Flash-P" and only the
    # node/edge counts distinguish them.
    say("APSIM: " + " vs ".join(REGISTRY[a].fig_label for a in ARMS))
    say("=" * 74)
    say("Same parameter chain, same APSIM files, same weather. Only the network differs.")

    for fname in FILES:
        got = {arm: stats(RUNS / arm / fname) for arm in ARMS}
        genos = sorted(set().union(*(got[a] for a in ARMS)),
                       key=lambda g: got[ARMS[0]].get(g, {}).get("mean", 0))

        say(f"\n\n### {fname}")
        say(f"  {'cultivar':<10}{'arm':<14}{'n':>5}{'yield':>9}{'range':>17}"
            f"{'meanTTN':>9}{'maxTTN':>8}{'meanFTN':>9}")
        for g in genos:
            for arm in ARMS:
                d = got[arm].get(g)
                if not d:
                    continue
                say(f"  {g:<10}{arm:<14}{d['n']:>5}{d['mean']:>9.1f}"
                    f"{d['lo']:>8.0f}..{d['hi']:<8.0f}"
                    f"{d['ttn']:>9.2f}{d['ttn_max']:>8.2f}{d['ftn']:>9.2f}")
            say("")

        # ---- 1. ordering ----
        say("  ORDERING (mean yield, low to high)")
        for arm in ARMS:
            order = sorted(got[arm], key=lambda g: got[arm][g]["mean"])
            say(f"    {arm:<14}{' < '.join(order)}")
        same = len({tuple(sorted(got[a], key=lambda g: got[a][g]["mean"]))
                    for a in ARMS}) == 1
        say(f"    -> {'SAME ordering' if same else 'ORDERING DIFFERS'}")

        # ---- 2. arm agreement, each arm against the reference ----
        for arm in [a for a in ARMS if a != REFERENCE]:
            if REFERENCE not in got:
                break
            say(f"\n  ARM AGREEMENT ({arm} vs {REFERENCE})")
            for g in genos:
                a, b = got[REFERENCE].get(g), got[arm].get(g)
                if not (a and b):
                    continue
                d = b["mean"] - a["mean"]
                say(f"    {g:<10}{a['mean']:>8.1f} -> {b['mean']:<8.1f}"
                    f"{d:>+8.1f} kg/ha  ({d / a['mean'] * 100:+.2f}%)")

        # ---- 3. hand-curated vs published ----
        pub_path = PUB / fname
        if pub_path.exists():
            pub = stats(pub_path)
            shared = sorted(set(pub) & set(got["hand_curated"]))
            say(f"\n  SANITY: hand-curated arm vs published .db  (shared: {', '.join(shared)})")
            for g in shared:
                a, p = got["hand_curated"][g], pub[g]
                d = (a["mean"] - p["mean"]) / p["mean"] * 100
                flag = "ok" if abs(d) <= 5 else "OUT OF BAND"
                say(f"    {g:<10}published {p['mean']:>7.1f}   ours {a['mean']:>7.1f}"
                    f"   {d:>+6.2f}%  {flag}")
        else:
            say(f"\n  SANITY vs published .db skipped: {pub_path.relative_to(ROOT)} not found"
                "\n  (download it from the GitHub release, see README)")

    (ROOT / "03_runs" / "APSIM_Arm_Comparison.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    say("\n\nwrote 03_runs/APSIM_Arm_Comparison.txt")


if __name__ == "__main__":
    main()
