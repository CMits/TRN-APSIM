"""Figure 9 — Daily dynamics in a WET (High-tillering-wins) representative year.

Same 5×3 layout as Figures 7 and 8 (see figure_07_dynamics_dry.py for the
row / column / variable conventions). This figure is drawn for the year with
the largest mean High advantage across the three sowing dates — a
high-yield, favourable season where producing the extra tillers pays off.

Data source
-----------
  SimulationsV2/outputs/DensityManagement_NoTiller.db
  Tables: HarvestReport + DailyReport (Stage > 0).

Usage
-----
    python figure_09_dynamics_wet.py

Output
------
    Figures/Fig09_Dynamics_<year>_Wet.png
"""
from __future__ import annotations

from _figure_common import (
    FIG_DIR, SOWING_ORDER, apply_rcparams, load_db, growing_season_rain,
    pick_rep_years_per_sowing, all_distinct_years,
    dynamics_ylims, draw_dynamics,
)


SUBTITLE = 'Largest High-Tillering Advantage'


def main() -> None:
    apply_rcparams()
    harvest, daily = load_db()
    gs_rain = growing_season_rain(daily, harvest)
    per_sow = pick_rep_years_per_sowing(harvest)
    years   = {sow: per_sow[sow]['wet'] for sow in SOWING_ORDER}
    ylims   = dynamics_ylims(daily, harvest, all_distinct_years(per_sow))
    print(f'Wet per-sowing years = {years}')
    out = FIG_DIR / 'Fig09_Dynamics_Wet.png'
    draw_dynamics(daily, harvest, years, SUBTITLE, gs_rain, ylims, out)


if __name__ == '__main__':
    main()
