"""Convenience driver: regenerate Figures 5-9 in one command.

Calls each ``figure_0X_*.py`` script's ``main()`` in order. Equivalent to
running them one at a time — useful after a fresh APSIM re-run.

Usage
-----
    python make_all_figures.py
"""
from __future__ import annotations

import figure_05_crossover
import figure_06_win_rate
import figure_07_dynamics_dry
import figure_08_dynamics_crossover
import figure_09_dynamics_wet


def main() -> None:
    figure_05_crossover.main()
    figure_06_win_rate.main()
    figure_07_dynamics_dry.main()
    figure_08_dynamics_crossover.main()
    figure_09_dynamics_wet.main()


if __name__ == '__main__':
    main()
