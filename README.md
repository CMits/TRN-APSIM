# TRN-APSIM: from a trait regulatory network to sorghum yield

Code and data to reproduce the analysis linking a qualitative **trait regulatory network (TRN)** for
shoot branching to the **APSIM Next Generation** sorghum crop model. Gene perturbations are propagated
through the network, turned into the APSIM tillering parameters `tillerSdIntercept` and
`tillerSdSlope` for each genotype, and simulated at field scale (Bullarah, NSW, Australia).

The pipeline runs twice, with everything downstream held identical:

| Arm | Network | Size | Simulated with |
|---|---|---|---|
| `hand_curated` | Hand-curated shoot-branching network (PSoup) | 32 nodes, 56 edges | R, the PSoup model's own equations |
| `flashp_suc` | Sorghum tillering network built automatically with Flash-P, plus 4 sucrose edges | 30 nodes, 44 edges | Python, the Flash-P algebraic engine |

## The pipeline

1. **Network response.** Each of the 32 genotypes is simulated in each network at sucrose levels 0.5, 1.0 and 1.5.
2. **Arabidopsis anchor.** Responses are anchored to published *Arabidopsis* shoot-branch counts (fold change against each experiment's own wild type).
3. **Supply/demand mapping.** Low and high sucrose are placed on APSIM's supply/demand (S/D) axis.
4. **Linear tillering line.** An intercept and slope are fitted per genotype.
5. **Rescale to APSIM.** Values are min–max rescaled into APSIM's published ranges (intercept −1.363 to 1.5, slope 0 to 0.4) and written as APSIM cultivars.

APSIM is then run for a high-tillering (`max2`) and a low-tillering (`nced3`) genotype across
historical seasons at three sowing dates, and the yield difference is compared between seasons.

## Requirements

| Tool | Version | Needed for |
|---|---|---|
| Python | ≥ 3.10 | sweeps, parameters, comparisons, figures (`pip install -r environment/requirements.txt`) |
| R | ≥ 4.4, **base R only** | the hand-curated sucrose sweep (see `environment/R-packages.md`) |
| APSIM Next Generation | tested with 2026.1.7961.0 | the crop simulations (<https://www.apsim.info>) |
| PowerShell | 5.1 (Windows) or 7+ (`pwsh`, macOS/Linux) | `run_all.ps1` |

`run_all.ps1` finds Python, `Rscript` and APSIM `Models` on the PATH or in the usual install
folders, and stops with a clear message if one is missing.

## Running it

```pwsh
pip install -r environment/requirements.txt
pwsh ./run_all.ps1                     # everything, both networks, including APSIM
pwsh ./run_all.ps1 -SkipApsim          # phases 00-06 only: network sweeps, parameters, methods figures
pwsh ./run_all.ps1 -Arm hand_curated   # one network only
pwsh ./run_all.ps1 -Arm ?              # list the networks
```

| Phase | Does | Writes |
|---|---|---|
| 00 | check Python, R, APSIM | |
| 01 | gene/genotype coverage of each network; Python cross-check of the PSoup model | `01_audit/` |
| 02 | add the sucrose edges to the Flash-P network | `02_networks/Tillering_Sorghum_Benchmark_suc/` |
| 03 | sucrose sweeps, each network at SUC 0.5 / 1.0 / 1.5 | `01_audit/*_sucrose_response.csv` |
| 04 | pattern and solver comparisons | `01_audit/` |
| 05 | network response → APSIM `tillerSdIntercept` / `tillerSdSlope`; methods figures | `03_runs/Genotype_APSIM_Parameters_<arm>.csv`, `04_figures/methods/` |
| 06 | write a patched `.apsimx` per network from the original templates | `03_runs/apsim/<arm>/` |
| 07 | run APSIM | `03_runs/apsim/<arm>/*.db` |
| 08 | compare the networks' APSIM output, and the published run | `03_runs/*Comparison*.txt` |
| 09 | result figures | `04_figures/` |
| 10 | compare the new parameters with `expected_results/` | |

A full run takes a few minutes, most of it in APSIM. Phase 10 fails if your parameters differ from
the ones used in the paper.

## Figures and tables in the chapter

| Chapter item | File (in `04_figures/`, reference copy in `expected_results/figures/`) |
|---|---|
| Fig. 5.6 (the two networks) | `methods/FigM02_Networks_SideBySide`, `methods/FigM02b_Networks_Cytoscape` |
| Fig. 5.7 (mutant response to sucrose) | `methods/FigM01_Sucrose_Response` |
| Fig. 5.8 (shoot branches vs S/D ratio) | `methods/FigM04_05_SupplyDemand_Combined` (also `FigM04`, `FigM05`, `FigM06`) |
| Fig. 5.10 (LT vs HT yield difference) | `figures_combined/Fig05_Crossover` |
| Suppl. Fig. 5.1 (Arabidopsis raw data) | `methods/FigM03_Arabidopsis_Branches` |
| Suppl. Fig. 5.2 (potential tiller number) | `figures_combined/FigPTN_TillerNumber` |
| Suppl. Fig. 5.3 (genotype ranks) | `figures_combined/FigRank_GenotypeRanks` |
| Suppl. Fig. 5.4 (tillerSdSlope) | `figures_combined/FigSlope_SupplyDemandSlope` |
| Suppl. Fig. 5.5 (win rate) | `figures_combined/Fig06_LowTilleringWins` |
| Suppl. Figs 5.6–5.11 (crop dynamics) | `<arm>/Fig07_Dynamics_Dry`, `<arm>/Fig08_Dynamics_Crossover`, `<arm>/Fig09_Dynamics_Wet` |
| Suppl. Table 5.2 (intercepts and slopes) | `03_runs/Genotype_APSIM_Parameters_<arm>.csv`, `03_runs/Arm_Comparison_long.csv` |

## The published APSIM run (optional download)

Phase 08 also compares the hand-curated arm with the APSIM output of the original study. These
`.db` files (~90 MB) are too large for git and are attached to the
[GitHub release](https://github.com/CMits/TRN-APSIM/releases). To include that check, unzip
`apsim_published.zip` into `00_reference/apsim_published/`. Without it, that comparison is skipped
and everything else runs.

## Repository layout

| Folder | Holds |
|---|---|
| `run_all.ps1` | Runs the whole pipeline |
| `scripts/` | Every script the pipeline runs; `scripts/arms.json` defines the two networks |
| `00_reference/paper_data/` | Input data: PSoup sucrose results, Arabidopsis branch counts, transformation weights, S/D range, published APSIM parameter set |
| `00_reference/apsim_inputs/` | Original APSIM templates (`*.apsimx.bak`) and the Bullarah weather file |
| `00_reference/psoup_source/` | Python port of the PSoup update rule (cross-check) |
| `00_reference/flashp_source/` | The Flash-P simulation engine used for the Flash-P network |
| `02_networks/Hand_Curated/` | The hand-curated network (SBGN, edges, nodes) and its PSoup R model |
| `02_networks/Tillering_Sorghum_Benchmark/` | The Flash-P export as built (network, evidence, validation) |
| `02_networks/Tillering_Sorghum_Benchmark_suc/` | The same network with the sucrose edges added (regenerated by phase 02) |
| `expected_results/` | The parameters and figures used in the paper, for comparison |
| `environment/` | Python requirements and R notes |

Generated folders (`01_audit/`, `03_runs/`, `04_figures/`) are created by `run_all.ps1` and are not
under version control.

## Notes on interpretation

- **The hand-curated arm is a fresh R run.** `run_psoup_sweep.R` executes the PSoup model's own
  equations. `01_audit/handcurated_vs_published.csv` records how it compares with the earlier
  published sweep (`DataSetS1`).
- **Only the relative pattern reaches APSIM.** Step 5 min–max rescales within each network, so the
  absolute scale of a network's output and the wild-type reference in step 2 cancel out. What has to
  agree between networks is the genotype ordering (Spearman ρ = 0.998 between the two arms).
- **The crossover figure cannot tell the networks apart.** `nced3` is clamped to zero tillers and
  `max2` sits far above the fertile-tiller ceiling in both arms, so Fig. 5.10 is identical for both
  networks. The network comparison is in the rank and value figures (Suppl. Figs 5.2–5.4).
- `DatasetS2_Arabidopsis_branches.csv` is cp1252-encoded, not UTF-8.

## Data sources

- Weather: SILO daily climate data for Bullarah, NSW (−29.45, 149.05), 1972–2025, Queensland
  Government, <https://www.longpaddock.qld.gov.au/silo/>, licensed CC BY 4.0.
- Crop model: APSIM Next Generation (Holzworth et al., 2014), <https://www.apsim.info>. APSIM is not
  redistributed here; install it separately under its own licence.
- PSoup: Fortuna et al. (2026), <https://github.com/NicoleZFortuna/PSoup>.
- Arabidopsis shoot-branch counts and the genotype list come from the studies cited in the chapter
  (Supplementary Table 5.1 and Supplementary Fig. 5.1).
- `02_networks/Tillering_Sorghum_Benchmark*/data/evidence.json` holds the literature evidence
  (DOIs, titles, abstracts and short quoted sentences) that Flash-P used for each edge. The full
  texts of the articles are not included.

## Licence

Code: MIT. Data and results: CC BY 4.0. Third-party material keeps its own licence. See `LICENSE`.
