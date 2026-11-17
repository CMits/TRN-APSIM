"""Draw every figure of the chapter for the two networks (hand_curated, flashp_suc).

  04_figures/methods/            FigM01 .. FigM06 and the network diagrams
  04_figures/figures_combined/   the cross-arm outcome figures
  04_figures/<arm>/              Figures 5-9, one directory per arm

Nothing here draws anything. It runs the figure scripts in order, one subprocess
each (scripts/arms.py builds the arm registry at import time), with ARMS_ONLY and
FIGS_ROOT set, then strips the _AllNetworks / _BothNetworks filename suffixes,
which carry no information when only one pair of networks is drawn.

  python scripts/make_2network_figures.py                 # all figures
  python scripts/make_2network_figures.py --only methods  # no APSIM output needed
  python scripts/make_2network_figures.py --only results  # needs the APSIM .db files
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"

# The two arms to keep, in registry order, and where their figures go. Re-pointing
# this at a different pair is the only edit a different comparison needs -- but the
# reference arm has to stay in, and arms.py says so if it is dropped.
KEEP = ["hand_curated", "flashp_suc"]
DEST = "04_figures"

# Mirrors run_all.ps1 steps 05c-05g and 09a-09f. make_figures.py takes the arm
# names on the command line as well, because it hard-fails on a missing .db rather
# than skipping, and its own filter is argv-based.
STEPS: list[tuple[str, str, list[str]]] = [
    ("05c methods sucrose (FigM01)",    "methods_figure_sucrose.py",    []),
    ("05d methods networks (FigM02)",   "figures/network_diagram.py",   []),
    ("05e methods arabidopsis (FigM03)", "methods_figure_arabidopsis.py", []),
    ("05f methods S/D (FigM04, M05)",   "methods_figure_sd.py",         []),
    ("05g methods pair (FigM02b)",      "figures/network_pair_cytoscape.py", []),
    ("09a figures 5-9 per arm",         "make_figures.py",              KEEP),
    ("09b combined crossover",          "combined_figure_05.py",        []),
    ("09c combined ranks",              "combined_figure_rank.py",      []),
    ("09d combined values",             "combined_figure_values.py",    []),
    ("09e all-network pairs",           "combined_figure_pairs.py",     []),
    ("09f combined win rate",           "combined_figure_06.py",        []),
]

# stem -> stem, relative to DEST, applied to both .png and .svg.
RENAMES = {
    "methods/FigM01_Sucrose_Response_AllNetworks":
        "methods/FigM01_Sucrose_Response",
    "figures_combined/Fig05_Crossover_AllNetworks":
        "figures_combined/Fig05_Crossover",
    "figures_combined/Fig06_LowTilleringWins_AllNetworks":
        "figures_combined/Fig06_LowTilleringWins",
    "figures_combined/FigRank_GenotypeRanks_BothNetworks":
        "figures_combined/FigRank_GenotypeRanks",
    "figures_combined/FigPTN_TillerNumber_BothNetworks":
        "figures_combined/FigPTN_TillerNumber",
    "figures_combined/FigSlope_SupplyDemandSlope_BothNetworks":
        "figures_combined/FigSlope_SupplyDemandSlope",
}


def run(label: str, script: str, argv: list[str], env: dict) -> bool:
    print(f"\n=== {label} ===")
    r = subprocess.run([sys.executable, str(SCRIPTS / script), *argv],
                       env=env, text=True, capture_output=True)
    print(r.stdout, end="")
    if r.returncode != 0:
        print(r.stderr, end="")
    return r.returncode == 0


def strip_suffixes(dest: Path) -> None:
    """Drop the _AllNetworks / _BothNetworks suffixes -- see FILENAMES above."""
    print(f"\n=== renaming ===")
    for old_stem, new_stem in RENAMES.items():
        for ext in ("png", "svg"):
            old, new = dest / f"{old_stem}.{ext}", dest / f"{new_stem}.{ext}"
            if not old.exists():
                continue
            old.replace(new)          # replace, so a re-run is idempotent
            print(f"  {old.name}  ->  {new.name}")


def main() -> None:
    # --only methods : steps 05c-05g (no APSIM output needed)
    # --only results : steps 09a-09f (need the APSIM .db files)
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else "all"
    if only not in ("all", "methods", "results"):
        raise SystemExit(f"--only must be methods, results or all, not {only!r}")
    steps = [s for s in STEPS
             if only == "all" or (only == "methods") == s[0].startswith("05")]
    dest = ROOT / DEST
    env = dict(os.environ, ARMS_ONLY=",".join(KEEP), FIGS_ROOT=DEST)

    print(f"arms   {', '.join(KEEP)}")
    print(f"out    {dest}")

    failed = [label for label, script, argv in steps
              if not run(label, script, argv, env)]
    if failed:
        raise SystemExit(f"\nfailed: {', '.join(failed)}")

    strip_suffixes(dest)

    print(f"\n=== {dest.name} ===")
    for p in sorted(dest.rglob("*.png")) + sorted(dest.rglob("*.svg")):
        print(f"  {p.relative_to(dest).as_posix()}")

    # The whole point of the folder: the third network must be gone from it.
    stray = [p for p in dest.rglob("*") if "flashp_shoot" in p.name]
    if stray:
        raise SystemExit("\nflashp_shoot leaked into the two-network set:\n  " +
                         "\n  ".join(str(p) for p in stray))
    print(f"\n{len(list(dest.rglob('*.png')))} png + "
          f"{len(list(dest.rglob('*.svg')))} svg, no flashp_shoot anywhere")


if __name__ == "__main__":
    main()
