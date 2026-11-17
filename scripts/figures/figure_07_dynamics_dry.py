"""Figure 7 — Daily dynamics in a DRY (Low-tillering-wins) representative year.

5 rows × 3 columns. Rows:
  1. Tiller number — the per-genotype TTN cap set by intercept + slope ×
     S/D ratio early in the season. (FTN is omitted: the APSIM 2026.1
     DailyReport exposes it as a transient emergence spike that drops to
     zero after flowering, which would be misleading alongside the cap.)
  2. Leaf Area Index (m²/m²).
  3. Above-ground biomass (solid) and grain yield (dashed), t/ha.
  4. Plant-available water (mm).
  5. Cumulative in-season rainfall (mm).

Columns are the three sowing dates (Oct 15 / Nov 15 / Dec 15). Y-axis limits
are shared across Figures 7/8/9 so cross-figure comparisons stay honest.

Year is chosen data-driven by ``pick_rep_years`` — the year with the
largest mean (Low − High)/High across the three sowings,
considering only years with complete 3-genotype × 3-sowing coverage.

Data source
-----------
  SimulationsV2/outputs/DensityManagement_NoTiller.db
  Tables: HarvestReport + DailyReport (the DailyReport rows are filtered to
  Stage > 0 to retain the crop cycle only).

Usage
-----
    python figure_07_dynamics_dry.py

Output
------
    Figures/Fig07_Dynamics_<year>_Dry.png
"""
from __future__ import annotations

from _figure_common import (
    FIG_DIR, SOWING_ORDER, apply_rcparams, load_db, growing_season_rain,
    pick_rep_years_per_sowing, all_distinct_years,
    dynamics_ylims, draw_dynamics,
)


SUBTITLE = 'Largest Low-Tillering Advantage'


def main() -> None:
    apply_rcparams()
    harvest, daily = load_db()
    gs_rain = growing_season_rain(daily, harvest)
    per_sow = pick_rep_years_per_sowing(harvest)
    years   = {sow: per_sow[sow]['dry'] for sow in SOWING_ORDER}
    ylims   = dynamics_ylims(daily, harvest, all_distinct_years(per_sow))
    print(f'Dry per-sowing years = {years}')
    out = FIG_DIR / 'Fig07_Dynamics_Dry.png'
    draw_dynamics(daily, harvest, years, SUBTITLE, gs_rain, ylims, out)


if __name__ == '__main__':
    main()
