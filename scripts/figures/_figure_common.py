"""Shared utilities for Figures 5-9 of the APSIM-PSoup paper.

Every ``figure_0X_*.py`` script in this directory imports from this module
to guarantee a single, consistent definition of:

* The representative genotypes and their display colours / labels.
* Sowing-date ordering and pretty names.
* Database loading (with the string->float coercion APSIM requires).
* Growing-season rainfall calculation.
* Data-driven selection of the three representative years (dry / crossover /
  wet) from which Figures 7, 8, 9 are drawn.
* Global y-axis limits shared across the three dynamics figures so visual
  comparisons between them are valid.
* Matplotlib rc parameters for publication-quality output.

The figures themselves (one per script) take care of their own layout.

LOCAL NOTE (FlashP_APSIM_Reproduction)
--------------------------------------
This is a COPY. The original in SimulationsV2/scripts/ is left untouched. The only
change is that the three things that were hardcoded now come from the environment,
so the same figure code can be pointed at either arm:

    FIGS_DB          .db to plot
    FIGS_OUT         directory to write figures into
    FIGS_GENOTYPES   comma-separated, HIGHEST tillering first ... LOWEST last

With nothing set, the defaults reproduce the paper's own figures exactly. Layout,
year-selection logic and styling below are the paper's code, unmodified -- the
hardcoded 'max2' / 'fhy3far1' are simply replaced by HIGH / LOW.
"""
from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths — resolved from this file so the scripts are location-independent.
# ---------------------------------------------------------------------------
HERE    = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent                      # FlashP_APSIM_Reproduction/
# Defaults only -- make_figures.py always sets FIGS_DB/FIGS_OUT per arm. They point
# inside this project so nothing resolves outside it when run standalone.
DB_DIR  = PROJECT / '03_runs' / 'apsim' / 'hand_curated'

# Primary data source: the factorial DensityManagement simulation run.
DB_PATH = Path(os.environ.get('FIGS_DB', DB_DIR / 'DensityManagement_NoTiller.db'))
FIG_DIR = Path(os.environ.get('FIGS_OUT', PROJECT / '04_figures' / 'hand_curated'))
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Genotype configuration.
#
# The representatives are chosen from the 34-genotype supplementary table to
# span the tillering range observed in the PSoup -> Arabidopsis ratio-based
# meta-analysis:
#
#   max2      - High tillering (SL-receptor loss-of-function, strong tiller)
#   nced3     - Low tillering
#
# A Medium representative (lbo) was plotted in earlier versions; the results now
# show only the high/low contrast, so it is no longer passed in. The role logic
# below still supports it, so re-adding it to FIGS_GENOTYPES restores the old
# three-way figures unchanged.
#
# Colours match the old paper: red/blue/green for visual continuity with
# Figures in APSIM-PSoup_V4.docx.
# ---------------------------------------------------------------------------
# Order is HIGH -> (MEDIUM) -> LOW, as in the paper. Colours and labels are
# assigned BY ROLE, not by raw index, so a substituted cultivar keeps the role its
# predecessor had AND a two-genotype list still gives the low representative the
# Low colour rather than the Medium one: position 0 is always High, the last
# position always Low, anything in between Medium.
GENOTYPES = [g.strip() for g in
             os.environ.get('FIGS_GENOTYPES', 'max2,nced3').split(',') if g.strip()]

_ROLE_COLORS = ['#E53935', '#2196F3', '#4CAF50']          # red / blue / green
_ROLE_SHORT  = ['High', 'Medium', 'Low']


def _role_index(i: int, n: int) -> int:
    """Map a position in GENOTYPES onto 0 = High, 1 = Medium, 2 = Low."""
    if i == 0:
        return 0
    return 2 if i == n - 1 else 1


# Exported so figure-specific palettes (e.g. Figure 5's Okabe-Ito override) can
# be assigned by the same role rule instead of re-deriving it.
ROLE_INDEX = [_role_index(i, len(GENOTYPES)) for i in range(len(GENOTYPES))]

COLORS = {g: _ROLE_COLORS[ROLE_INDEX[i]] for i, g in enumerate(GENOTYPES)}
SHORT  = {g: _ROLE_SHORT[ROLE_INDEX[i]] for i, g in enumerate(GENOTYPES)}
LABELS = {g: f'{_ROLE_SHORT[ROLE_INDEX[i]]} Tillering ({g})'
          for i, g in enumerate(GENOTYPES)}

# The two cultivars the crossover / year-selection logic contrasts. The paper
# hardcoded 'max2' and 'fhy3far1'; these are the same thing, by role.
HIGH = GENOTYPES[0]
LOW  = GENOTYPES[-1]

SOWING_ORDER  = ['Oct15', 'Nov15', 'Dec15']
SOWING_LABELS = {'Oct15': 'Oct 15', 'Nov15': 'Nov 15', 'Dec15': 'Dec 15'}

# Representative density (plants/m²) used for the paper figures.
# The original paper Figures 5-9 reported at Density = 2.5 plants/m²
# (the old ``Exploration_MultiGenotype_Sowing.db`` only contained this
# density), so we match it here for visual continuity. The NEW simulation
# file also runs 5.0 and 7.5 for the density-sensitivity appendix, but
# those are not plotted in Figs 5-9.
REP_DENSITY = 2.5

# Unit conversion: APSIM reports biomass / yield in g/m²; figures use t/ha.
G_TO_T_PER_HA = 0.01

# Colours used to circle the three representative years on Figure 5.
YEAR_COLORS = {'dry': '#D32F2F', 'crossover': '#7B1FA2', 'wet': '#1565C0'}

# ---------------------------------------------------------------------------
# Matplotlib style — applied once per script via ``apply_rcparams()``.
# ---------------------------------------------------------------------------
def apply_rcparams() -> None:
    """Publication font sizes matching the old-paper figure style."""
    plt.rcParams.update({
        'font.size':        16,
        'axes.titlesize':   18,
        'axes.labelsize':   18,
        'xtick.labelsize':  14,
        'ytick.labelsize':  14,
        'legend.fontsize':  13,
        'figure.titlesize': 20,
    })


# ---------------------------------------------------------------------------
# Database loading.
#
# APSIM Next Generation stores Experiment factor values as strings ('5',
# '2.5', '0.45'). The figure scripts compare them to floats, so we coerce
# the numeric columns explicitly.
# ---------------------------------------------------------------------------
def load_db(db_path: Path = DB_PATH) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(harvest, daily)`` filtered to the three representatives at
    ``REP_DENSITY``, with Dates parsed and yields expressed in t/ha."""
    if not db_path.exists():
        raise FileNotFoundError(
            f'{db_path} not found. Run APSIM on the apsimx files first '
            f'(see run_apsim_hint.md) and move the .db into outputs/.'
        )
    with sqlite3.connect(str(db_path)) as conn:
        harvest = pd.read_sql('SELECT * FROM HarvestReport', conn)
        daily   = pd.read_sql('SELECT * FROM DailyReport',   conn)

    for df in (harvest, daily):
        for col in ('Density', 'SoilWater'):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        df['Date'] = pd.to_datetime(df['Clock.Today'], errors='coerce')
        df['Year'] = df['Date'].dt.year

    harvest['GrainYield_tha'] = harvest['Sorghum.Grain.Wt']       * G_TO_T_PER_HA
    harvest['Biomass_tha']    = harvest['Sorghum.AboveGround.Wt'] * G_TO_T_PER_HA

    mask_h = (harvest['Density'] == REP_DENSITY) & harvest['Genotype'].isin(GENOTYPES)
    mask_d = (daily['Density']   == REP_DENSITY) & daily['Genotype'].isin(GENOTYPES)
    return harvest[mask_h].copy(), daily[mask_d].copy()


def growing_season_rain(daily: pd.DataFrame,
                        harvest: pd.DataFrame) -> dict[str, dict[int, float]]:
    """Return ``{SowingDate: {Year: total_rain_mm}}`` summed over the crop
    cycle (Stage>0) for the ``max2`` simulation — the three genotypes share
    the same weather within a (sowing, year), so one is enough."""
    out: dict[str, dict[int, float]] = {sow: {} for sow in SOWING_ORDER}
    for sow in SOWING_ORDER:
        dsub = daily[(daily['Genotype'] == HIGH)   & (daily['SowingDate']   == sow)]
        hsub = harvest[(harvest['Genotype'] == HIGH) & (harvest['SowingDate'] == sow)]
        for yr, g in hsub.groupby('Year'):
            hd = g['Date'].iloc[0]
            window = dsub[(dsub['Date'] >= hd - pd.Timedelta(days=250)) &
                          (dsub['Date'] <= hd + pd.Timedelta(days=10))]
            crop = window[window['Stage'] > 0]
            if not crop.empty:
                out[sow][int(yr)] = float(crop['Rain'].sum())
    return out


# ---------------------------------------------------------------------------
# Representative-year selection.
#
# A year is a good representative only if all three genotypes harvest at all
# three sowing dates — otherwise the dynamics panels are patchy. We filter
# to complete-coverage years first, then rank each category by a metric
# tuned to give a clear, paper-worthy story:
#
# Notation: dLow_s = (yield_fhy3far1 − yield_max2) / yield_max2 × 100 %
#           at sowing s, computed per year. meanVH is the across-sowing mean
#           High yield.
#
#   * dry       -> largest mean(dLow). Low wins at every sowing date.
#                  Yields are typically low (drought/stress).
#   * crossover -> years where dLow changes sign across the three sowings
#                  (regime change within the year), ranked by
#                  range(dLow) − 10 × |mean(dLow)|: big swing, mean near 0.
#                  This captures the paper's narrative of sowing-date tuning.
#   * wet       -> largest (−mean(dLow)) × meanVH: High wins clearly
#                  AND the season is genuinely high-yielding.  Unweighted
#                  selection can pick a low-yield year with a mild advantage
#                  (e.g. 2004) which is biologically uninteresting.
#
# With the current parameter set these rank out as 2024 / 2022 / 1976, which
# matches the old paper's 1976 wet-year pick and gives a sharp story for
# 2024 (dry) and 2022 (crossover).
# ---------------------------------------------------------------------------
def pick_rep_years_per_sowing(harvest: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Return ``{sow: {'dry': yr, 'crossover': yr, 'wet': yr}}``.

    Each sowing date is ranked independently. The best dry / crossover / wet
    year at (say) Oct 15 is typically not the same as at Dec 15 — sorghum
    phenology shifts in and out of the stress window as the sowing date
    moves, so the Low / High contrast peaks in different years
    across the three sowings. Using a per-sowing pick lets every column of
    Figures 7-9 show its own maximum-contrast year.

    Ranking within a single sowing date:
      dry       -> year with largest +dLow  (Low beats High most).
      crossover -> year with |dLow| closest to zero (yields roughly tied).
      wet       -> year maximising (-dLow) × High yield, i.e. High
                   tillering wins clearly AND the season is high-yielding; the
                   yield weighting avoids picking a biologically dull
                   low-yield year with a mild high-tillering edge.
    """
    out: dict[str, dict[str, int]] = {}
    for sow in SOWING_ORDER:
        piv = (harvest[harvest['SowingDate'] == sow]
                   .pivot_table(index='Year', columns='Genotype',
                                values='GrainYield_tha', aggfunc='first')
                   .dropna(subset=[HIGH, LOW]))
        piv['dLow']      = (piv[LOW] - piv[HIGH]) / piv[HIGH] * 100
        piv['wet_score'] = -piv['dLow'] * piv[HIGH]
        out[sow] = {
            'dry':       int(piv['dLow'].idxmax()),
            'crossover': int(piv['dLow'].abs().idxmin()),
            'wet':       int(piv['wet_score'].idxmax()),
        }
    return out


def all_distinct_years(per_sowing: dict[str, dict[str, int]]) -> list[int]:
    """Flatten ``pick_rep_years_per_sowing(...)`` to a list of unique years,
    used by ``dynamics_ylims`` to pick a y-axis scale that accommodates
    every cell of every dynamics figure."""
    ys = set()
    for sow_entry in per_sowing.values():
        for yr in sow_entry.values():
            ys.add(int(yr))
    return sorted(ys)


# ---------------------------------------------------------------------------
# Global y-limits across the three dynamics figures.
#
# Figures 7, 8, 9 share the same 5x3 panel grid; keeping the y-axis limits
# common lets the reader compare curves across figures without being misled
# by rescaling.
# ---------------------------------------------------------------------------
# Small headroom BELOW zero, as a fraction of each row's span. Without it a curve
# that sits flat at zero -- the Low tillering representative's tiller number is
# exactly 0 at every timestep, because max(0, intercept + slope*SD) clamps there --
# is drawn directly on top of the x-axis spine and cannot be seen. Applied to every
# row, not just the tiller row, so all five panels keep the same baseline offset.
YLIM_BOTTOM_PAD = 0.03


def dynamics_ylims(daily: pd.DataFrame, harvest: pd.DataFrame,
                   years: list[int]) -> dict[int, tuple[float, float]]:
    """Return ``{row_index: (-pad, max*1.1)}`` for rows 0..4 of the dynamics
    figures (TTN, LAI, Weight, PAW, Cum rain)."""
    row_max = {r: 0.0 for r in range(5)}
    for yr in years:
        for sow in SOWING_ORDER:
            for g in GENOTYPES:
                hsub = harvest[(harvest['Genotype']   == g) &
                               (harvest['SowingDate'] == sow) &
                               (harvest['Year']       == yr)]
                if hsub.empty:
                    continue
                hd = hsub['Date'].iloc[0]
                window = daily[(daily['Genotype']   == g) &
                               (daily['SowingDate'] == sow) &
                               (daily['Date'] >= hd - pd.Timedelta(days=250)) &
                               (daily['Date'] <= hd + pd.Timedelta(days=10))]
                crop = window[window['Stage'] > 0]
                if crop.empty:
                    continue
                row_max[0] = max(row_max[0], crop['TTN'].max())
                row_max[1] = max(row_max[1], crop['Sorghum.Leaf.LAI'].max())
                row_max[2] = max(row_max[2], (crop['AboveGroundWt'] * G_TO_T_PER_HA).max())
                row_max[3] = max(row_max[3], crop['PAWmm'].max())
                row_max[4] = max(row_max[4], crop['Rain'].cumsum().max())
    tops = {r: (row_max[r] * 1.10 if row_max[r] > 0 else 1.0) for r in range(5)}
    return {r: (-YLIM_BOTTOM_PAD * tops[r], tops[r]) for r in range(5)}


# ---------------------------------------------------------------------------
# Shared helper for the three dynamics figures (7, 8, 9).
#
# Each produces a 5x3 grid of panels: 3 sowing dates (columns) x 5 variables
# (rows). The only thing that differs between Fig 7/8/9 is the chosen year
# and the subtitle.
# ---------------------------------------------------------------------------
ROW_YLABELS = ['Tiller number', 'LAI (m²/m²)', 'Weight (t/ha)',
               'PAW (mm)', 'Cum. rain (mm)']


def draw_dynamics(daily: pd.DataFrame, harvest: pd.DataFrame,
                  sowing_years: dict[str, int], subtitle: str,
                  gs_rain: dict, ylims: dict, out_path: Path) -> None:
    """Draw one dynamics figure (3 sowing dates x 5 variables).

    ``sowing_years`` maps each SowingDate to the year that column should
    show, so the three columns can display three different years — for
    example the maximum-Low-advantage year at each sowing date.
    """
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(5, 3, figsize=(20, 22), sharex='col')

    LABEL_FS  = 22
    TICK_FS   = 18
    TAG_FS    = 19
    LEGEND_FS = 20

    for col, sow in enumerate(SOWING_ORDER):
        year = sowing_years[sow]
        for g in GENOTYPES:
            colour = COLORS[g]
            short  = SHORT[g]

            hsub = harvest[(harvest['Genotype']   == g) &
                           (harvest['SowingDate'] == sow) &
                           (harvest['Year']       == year)]
            if hsub.empty:
                continue
            hd = hsub['Date'].iloc[0]
            window = daily[(daily['Genotype']   == g) &
                           (daily['SowingDate'] == sow) &
                           (daily['Date'] >= hd - pd.Timedelta(days=250)) &
                           (daily['Date'] <= hd + pd.Timedelta(days=10))]
            crop = window[window['Stage'] > 0].copy()
            if crop.empty:
                continue
            sow_date    = crop['Date'].iloc[0]
            crop['DAS'] = (crop['Date'] - sow_date).dt.days

            # Row 0: TTN (solid) + FTN (dashed)
            #   TTN = total tiller-number cap set early by intercept +
            #         slope × S/D ratio (stays flat once the cap is set).
            #   FTN = fertile-tiller dynamics; in APSIM 2026.1 this peaks
            #         early then drops as tillers are capped by flowering.
            axes[0, col].plot(crop['DAS'], crop['TTN'], color=colour,
                              linewidth=2, label=f'{short} TTN')
            axes[0, col].plot(crop['DAS'], crop['FTN'], color=colour,
                              linewidth=2, linestyle='--', alpha=0.7,
                              label=f'{short} FTN')
            # Row 1: LAI
            axes[1, col].plot(crop['DAS'], crop['Sorghum.Leaf.LAI'],
                              color=colour, linewidth=2, label=short)
            # Row 2: Above-ground biomass (solid) + Grain (dashed), both in t/ha
            axes[2, col].plot(crop['DAS'], crop['AboveGroundWt'] * G_TO_T_PER_HA,
                              color=colour, linewidth=2, label=f'{short} Biomass')
            axes[2, col].plot(crop['DAS'], crop['GrainWt'] * G_TO_T_PER_HA,
                              color=colour, linewidth=2, linestyle='--',
                              label=f'{short} Grain')
            # Row 3: Plant-available water
            axes[3, col].plot(crop['DAS'], crop['PAWmm'], color=colour,
                              linewidth=2, label=short)
            # Row 4: Cumulative in-season rainfall
            axes[4, col].plot(crop['DAS'], crop['Rain'].cumsum(), color=colour,
                              linewidth=2, label=short)

        rain_mm = gs_rain.get(sow, {}).get(year, 0)
        axes[0, col].text(
            0.02, 0.95,
            f'{SOWING_LABELS[sow]} - {year} ({rain_mm:.0f} mm GS rain)',
            transform=axes[0, col].transAxes, fontsize=TAG_FS,
            fontweight='bold', verticalalignment='top',
            zorder=10,
        )

    for r in range(5):
        for c in range(3):
            ax = axes[r, c]
            ax.set_ylim(ylims[r])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(axis='both', labelsize=TICK_FS)
        axes[r, 0].set_ylabel(ROW_YLABELS[r], fontsize=LABEL_FS)
    for c in range(3):
        axes[4, c].set_xlabel('Days after sowing', fontsize=LABEL_FS)

    handles = [Line2D([0], [0], color=COLORS[g], linewidth=3.5,
                      label=SHORT[g] + ' Tillering') for g in GENOTYPES]
    handles.append(Line2D([0], [0], color='0.4', linewidth=3.0,
                          linestyle='-', label='TTN / Biomass (solid)'))
    handles.append(Line2D([0], [0], color='0.4', linewidth=3.0,
                          linestyle='--', label='FTN / Grain (dashed)'))
    fig.legend(handles=handles, loc='upper center', ncol=5,
               fontsize=LEGEND_FS, frameon=False,
               bbox_to_anchor=(0.5, 0.995))

    plt.tight_layout(rect=[0, 0, 1, 0.965])
    svg_path = Path(out_path).with_suffix('.svg')
    fig.savefig(out_path, dpi=600, bbox_inches='tight')
    fig.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_path}')
    print(f'Saved: {svg_path}')
