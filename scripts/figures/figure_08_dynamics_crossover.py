"""Figure 8 — Daily dynamics in the CROSSOVER representative year.

Same 5×3 layout as Figures 7 and 9 (see figure_07_dynamics_dry.py for the
row / column / variable conventions). This figure captures the year where
the mean percent difference between Low and High across the three
sowing dates is closest to zero — the genotypes swap dominance depending on
sowing date, which is the regime change the paper highlights.

Data source
-----------
  SimulationsV2/outputs/DensityManagement_NoTiller.db
  Tables: HarvestReport + DailyReport (Stage > 0).

Usage
-----
    python figure_08_dynamics_crossover.py

Output
------
    Figures/Fig08_Dynamics_<year>_Crossover.png
"""
from __future__ import annotations

from _figure_common import (
    FIG_DIR, SOWING_ORDER, apply_rcparams, load_db, growing_season_rain,
    pick_rep_years_per_sowing, all_distinct_years,
    dynamics_ylims, draw_dynamics,
)


SUBTITLE = 'Crossover Years - Equal Performance per Sowing'


def main() -> None:
    apply_rcparams()
    harvest, daily = load_db()
    gs_rain = growing_season_rain(daily, harvest)
    per_sow = pick_rep_years_per_sowing(harvest)
    years   = {sow: per_sow[sow]['crossover'] for sow in SOWING_ORDER}
    ylims   = dynamics_ylims(daily, harvest, all_distinct_years(per_sow))
    print(f'Crossover per-sowing years = {years}')
    out = FIG_DIR / 'Fig08_Dynamics_Crossover.png'
    draw_dynamics(daily, harvest, years, SUBTITLE, gs_rain, ylims, out)


if __name__ == '__main__':
    main()
