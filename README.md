# hSEM_simulations

A simulation study comparing methods for inferring species interactions from camera trap data.

This repository contains the code for a simulation study that evaluates eleven modelling approaches for estimating species interactions, and frames hierarchical structural equation models (hSEMs) as a general framework within which familiar approaches (GLMs, SEMs, and hierarchical models) arise as special cases. The models are compared on their ability to recover known interaction strengths from data generated under a three-species predator-mesopredator-prey process.

## The models

The eleven approaches span three model families, each fitted under a naive, two-step, or joint (hSEM) framework:

- **Abundance / N-mixture** formulations
- **Royle-Nichols** formulations
- **Occupancy** formulations

Each is evaluated across scenarios that vary species abundance (abundant / rare), the number of sites (50 / 500 / 1000), the presence of a mesopredator, and landscape covariates.

## Pipeline

The scripts run in order:

1. `hSEM/01_simulate_and_fit.R` — simulates data under the three-species data-generating process and fits the eleven models. This is the simulation driver, run once per scenario on an HPC. It reads its scenario from the environment variables `SCENARIO`, `LANDSCAPE`, and `SETTING` (see "Running on an HPC" below) rather than command-line arguments. For each model it records point estimates and 95% interval bounds (Wald confidence intervals for the frequentist models, posterior credible intervals for the hSEMs) for the three species interaction values.
2. `hSEM/02_extract_results.R` — extracts parameter estimates and interval bounds from the fitted model outputs and assembles them, with the known truth, into a single comparison table across scenarios.
3. `hSEM/03_make_figures.R` — produces the percent-bias recovery figures.
4. `hSEM/04_coverage_analysis.R` — computes interval coverage (the proportion of replicates whose 95% interval contains the truth) and mean interval width for each model, interaction pathway, and scenario.

Two further scripts turn the coverage output into paper-ready outputs:

- `hSEM/build_table_S3.R` — builds the coverage-and-width supplementary table (Table S3) as a formatted Word document.
- `hSEM/coverage_width_figure.R` — produces the coverage-against-width figure.

Both exclude three occupancy two-step replicates (nSites = 50, abundant scenario) that failed to converge, consistent with the bias analysis.

`hSEM/Chris_original_pred_prey_sim.R` is the original predator-prey simulation framework (developed by C. Sutherland) that this work adapts, retained for reference and attribution.

## Running on an HPC

The simulation runs as a set of SLURM array jobs, one per scenario. `SLURM/generate_slurm_v4rn.sh` writes those per-scenario job scripts: it loops over the twelve scenarios (abundant / rare, at 50 / 500 / 1000 sites, with and without a mesopredator) crossed with landscape on and off, and sets walltime and memory according to the number of sites. Each generated job is a 100-replicate array (`--array=1-100%10`) that exports `SCENARIO` and `LANDSCAPE` and calls `01_simulate_and_fit.R`.

The generator is specific to the cluster it was written for. To run it elsewhere, adapt the account (`--account`), the absolute paths, and the module load line to your own system.

## Data

This is a simulation study. All data is generated within the scripts from the data-generating process, so no external data is required to reproduce the results. Set the seed in `01_simulate_and_fit.R` to reproduce a specific run.

## Dependencies

- **R** 4.4.2
- **NIMBLE** for fitting the hierarchical (hSEM) models, with **coda** for posterior diagnostics
- **JAGS**, installed separately, if running the JAGS-based fits
- R packages: `nimble`, `coda`, `unmarked`, `glmmTMB`, `lme4`, `piecewiseSEM`, `tidyverse` (`dplyr`, `tidyr`, `ggplot2`), `patchwork`, `scales`, `paletteer`, `ggrepel`, `flextable`, `officer`

The packages map onto the model families compared in the study: `unmarked` for the occupancy, N-mixture, and Royle-Nichols fits, `piecewiseSEM` for the SEM approaches, `lme4` and `glmmTMB` for the GLM baselines, and `nimble` for the joint hSEMs. `flextable` and `officer` build the supplementary table; `ggrepel` labels the coverage figure.

## Citation

TBC

## License

Released under the MIT License. See `LICENSE`.