"""Regenerate Figures 5-9 for both arms, using the paper's own figure code.

The figure scripts under scripts/figures/ are a COPY of SimulationsV2/scripts/,
changed only so that the .db, the output directory and the three representative
genotypes come from the environment instead of being hardcoded. The plotting
itself is the paper's, unmodified -- that is what makes the two sets comparable.

Each arm runs in its own subprocess because _figure_common reads the environment
at import time, so one process can only ever hold one configuration.

  04_figures/hand_curated/   Figures 5-9 from the hand-curated network
  04_figures/flashp_suc/     Figures 5-9 from the FlashP network
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
FIGDIR = ROOT / "scripts" / "figures"
RUNS = ROOT / "03_runs" / "apsim"

from arms import arm_names, FIG_ROOT  # noqa: E402

OUT = FIG_ROOT

ARMS = arm_names()


def _selected_arms() -> list[str]:
    """Optional arm filter, so `run_all.ps1 -Arm <x>` does not demand the other arm's .db."""
    if len(sys.argv) > 1:
        picked = [a for a in sys.argv[1:] if a in ARMS]
        bad = [a for a in sys.argv[1:] if a not in ARMS]
        if bad:
            raise SystemExit(f"unknown arm(s): {', '.join(bad)} -- choose from {', '.join(ARMS)}")
        return picked
    return ARMS

# HIGH -> LOW, matching the paper's ordering. nced3 replaces fhy3far1 as the Low
# representative (fhy3far1 has no confirmed grass orthologue). The Medium
# representative (lbo) is still staged into APSIM but no longer plotted -- the
# results show the high/low contrast only. Re-inserting "lbo" between the two
# restores the previous three-genotype figures without any other change.
GENOTYPES = "max2,nced3"


def main() -> None:
    arms = _selected_arms()
    failed = []
    for arm in arms:
        db = RUNS / arm / "DensityManagement_NoTiller.db"
        if not db.exists():
            raise SystemExit(f"missing {db} -- run build_apsim_runs.py and APSIM first")
        dest = OUT / arm
        dest.mkdir(parents=True, exist_ok=True)

        env = dict(os.environ,
                   FIGS_DB=str(db),
                   FIGS_OUT=str(dest),
                   FIGS_GENOTYPES=GENOTYPES)

        print(f"\n{'='*70}\n{arm}\n  db  {db}\n  out {dest}\n  {GENOTYPES}\n{'='*70}")
        r = subprocess.run([sys.executable, "make_all_figures.py"],
                           cwd=FIGDIR, env=env, text=True,
                           capture_output=True)
        print(r.stdout, end="")
        if r.returncode != 0:
            print(r.stderr, end="")
            failed.append(arm)

    print("\n" + "=" * 70)
    for arm in arms:
        pngs = sorted((OUT / arm).glob("*.png"))
        print(f"  {arm:<14}{len(pngs)} figures  {'FAILED' if arm in failed else ''}")
        for p in pngs:
            print(f"      {p.name}")
    if failed:
        raise SystemExit(f"failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
