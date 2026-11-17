"""Figure 6 — Win rate: proportion of years each tillering strategy achieves
the highest grain yield, by sowing date.

Three panels (Oct 15 / Nov 15 / Dec 15). For each sowing date we take every
year in the long-term record, identify which of the three representatives
gave the maximum yield in that year, and draw a bar chart of the resulting
percentages. Under each bar we report the mean yield of that genotype across
all years at that sowing date (context for the "win" statistic — a low win
rate with a high mean is a reliable runner-up).

Data source
-----------
  SimulationsV2/outputs/DensityManagement_NoTiller.db
  Table: HarvestReport, filtered to Density=5, Genotype ∈ {max2, lbo, fhy3far1}.

Usage
-----
    python figure_06_win_rate.py

Output
------
    Figures/Fig06_Win_Rate.png
"""
from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt

from _figure_common import (
    FIG_DIR, GENOTYPES, COLORS, SHORT, SOWING_ORDER, SOWING_LABELS,
    apply_rcparams, load_db,
)


def make_figure(harvest: pd.DataFrame, out_path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    for i, sow in enumerate(SOWING_ORDER):
        ax = axes[i]
        piv = (harvest[harvest['SowingDate'] == sow]
                   .pivot_table(index='Year', columns='Genotype',
                                values='GrainYield_tha', aggfunc='first')
                   [GENOTYPES].dropna(how='all'))
        if piv.empty:
            ax.set_title(f'(no data {SOWING_LABELS[sow]})')
            continue

        winners = piv.idxmax(axis=1)
        counts  = winners.value_counts()
        n       = len(winners)
        pcts    = [100 * counts.get(g, 0) / n for g in GENOTYPES]

        bars = ax.bar(range(len(GENOTYPES)), pcts,
                      color=[COLORS[g] for g in GENOTYPES],
                      alpha=0.85, edgecolor='black', linewidth=0.5, width=0.6)

        for j, (bar, p) in enumerate(zip(bars, pcts)):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 1.5,
                    f'{p:.0f}%\n({counts.get(GENOTYPES[j], 0)}/{n})',
                    ha='center', va='bottom', fontsize=14, fontweight='bold')

        ax.set_xticks(range(len(GENOTYPES)))
        ax.set_xticklabels([SHORT[g] for g in GENOTYPES], fontsize=14)
        ax.text(0.5, 0.95,
                f'({chr(ord("a") + i)}) {SOWING_LABELS[sow]} ({n} years)',
                transform=ax.transAxes, fontsize=16, fontweight='bold',
                ha='center', va='top')
        if i == 0:
            ax.set_ylabel('Win rate (% of years)')
        ax.set_ylim(0, max(pcts) * 1.4 if max(pcts) > 0 else 100)
        ax.grid(True, alpha=0.3, axis='y')

        for j, g in enumerate(GENOTYPES):
            if g in piv.columns:
                ax.text(j, -0.10 * ax.get_ylim()[1],
                        f'mean: {piv[g].mean():.2f} t/ha',
                        ha='center', va='top', fontsize=12,
                        color='black', fontweight='bold')

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_path}')


def main() -> None:
    apply_rcparams()
    harvest, _ = load_db()
    make_figure(harvest, FIG_DIR / 'Fig06_Win_Rate.png')


if __name__ == '__main__':
    main()
