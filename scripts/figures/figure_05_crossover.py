"""Figure 5 — Crossover: Low vs High tillering yield difference by sowing date.

Three panels (Oct 15 / Nov 15 / Dec 15). For each sowing date we plot, year
by year:

  x = grain yield of the High tillering representative (max2), t/ha
  y = (yield_LOW - yield_max2) / yield_max2 × 100 %

Only the Low vs High contrast is shown, and it is the only contrast the results
report. A second-order polynomial is fit to the Low scatter to visualise the
crossover with High tillering yield.

Each panel independently picks its own dry / crossover / wet representative
year (via ``pick_rep_years_per_sowing``) and circles those three points,
cross-referencing Figures 7, 8, 9. Mean growing-season rainfall is printed
as a small italic caption in the top-left of each panel.

Data source
-----------
  SimulationsV2/outputs/DensityManagement_NoTiller.db
  Table: HarvestReport, filtered to Density=2.5, Genotype ∈ {max2, fhy3far1}.

Usage
-----
    python figure_05_crossover.py

Output
------
    Figures/Fig05_Crossover_Yield.png
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _figure_common import (
    FIG_DIR, GENOTYPES, ROLE_INDEX, HIGH, LOW, SOWING_ORDER, SOWING_LABELS,
    apply_rcparams, load_db, growing_season_rain, pick_rep_years_per_sowing,
)

# Nature Plants-style palette (Okabe-Ito, colour-blind safe) — overrides the
# default red/blue/green only inside this figure so the other panels keep
# their established colours. Indexed through ROLE_INDEX rather than zipped
# positionally, so the Low representative keeps the green even when no Medium
# genotype is plotted.
_ROLE_PALETTE = ['#D55E00',   # Vermillion   — High tillering
                 '#0072B2',   # Deep blue    — Medium tillering
                 '#009E73']   # Bluish green — Low tillering
COLORS = {g: _ROLE_PALETTE[ROLE_INDEX[i]] for i, g in enumerate(GENOTYPES)}
YEAR_COLORS = {
    'dry':       '#B2182B',  # Deep red
    'crossover': '#542788',  # Muted purple
    'wet':       '#1F4E79',  # Navy
}


def make_figure(harvest: pd.DataFrame,
                per_sowing: dict[str, dict[str, int]],
                gs_rain: dict,
                out_path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(22, 8), sharey=True)

    annot_offset = {'dry': (-35, 30), 'crossover': (35, -35), 'wet': (40, 25)}
    annot_label  = {'dry': 'Dry', 'crossover': 'Intermediate', 'wet': 'Wet'}

    for i, sow in enumerate(SOWING_ORDER):
        ax = axes[i]
        piv = (harvest[harvest['SowingDate'] == sow]
                   .pivot_table(index='Year', columns='Genotype',
                                values='GrainYield_tha', aggfunc='first'))
        common = piv.dropna(subset=[HIGH, LOW]).copy()
        if common.empty:
            ax.set_title(f'(no data for {SOWING_LABELS[sow]})')
            continue

        common['dLow'] = ((common[LOW] - common[HIGH])
                          / common[HIGH]) * 100

        # Visual anchors: zero line + faint fill for the winning region.
        ax.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.5)
        ax.axhspan(0, 100,  alpha=0.04, color='#2196F3')  # Low wins
        ax.axhspan(-100, 0, alpha=0.04, color='#E53935')  # High wins

        # Low vs High scatter + quadratic fit.
        x = common[HIGH].values
        y = common['dLow'].values
        ax.scatter(x, y, color=COLORS[LOW], s=50, alpha=0.65,
                   zorder=3, label='Low vs High')
        if len(x) > 3:
            coef = np.polyfit(x, y, 2)
            xs = np.linspace(x.min(), x.max(), 100)
            ax.plot(xs, np.polyval(coef, xs),
                    color=COLORS[LOW], linewidth=2.5, alpha=0.75)

        # Circle + annotate this panel's own three representative years.
        rep_years = per_sowing[sow]
        for key, yr in rep_years.items():
            if yr not in common.index:
                continue
            row = common.loc[yr]
            c   = YEAR_COLORS[key]
            ax.scatter([row[HIGH]], [row['dLow']], s=200,
                       facecolors='none', edgecolors=c,
                       linewidths=2.5, zorder=5)
            ax.annotate(f"{yr}\n({annot_label[key]})",
                        xy=(row[HIGH], row['dLow']),
                        xytext=annot_offset[key], textcoords='offset points',
                        fontsize=12, fontweight='bold', color=c,
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.5),
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  edgecolor=c, alpha=0.9))

        ax.set_xlabel('High Tillering yield (t/ha)')
        if i == 0:
            ax.set_ylabel('Yield difference from High Tillering (%)')
        ax.text(0.02, 0.97, SOWING_LABELS[sow],
                transform=ax.transAxes, fontsize=16, fontweight='bold',
                verticalalignment='top')
        ax.legend(fontsize=12, loc='upper right', frameon=False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=600, bbox_inches='tight')
    svg_path = out_path.with_suffix('.svg')
    fig.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_path}')
    print(f'Saved: {svg_path}')


def main() -> None:
    apply_rcparams()
    harvest, daily = load_db()
    gs_rain    = growing_season_rain(daily, harvest)
    per_sowing = pick_rep_years_per_sowing(harvest)
    print(f'Per-sowing representative years: {per_sowing}')
    make_figure(harvest, per_sowing, gs_rain,
                FIG_DIR / 'Fig05_Crossover_Yield.png')


if __name__ == '__main__':
    main()
