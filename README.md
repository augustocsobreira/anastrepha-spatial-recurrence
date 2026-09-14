# Spatial structure and trap recurrence in the monitoring of *Anastrepha fraterculus*

Reproducible analysis code for the manuscript *Spatial Structure and Trap Recurrence Influence the Monitoring Efficiency of Anastrepha fraterculus* (commercial apple orchard, Vacaria, Rio Grande do Sul, Brazil, 2019 to 2025).

## Contents

| File | Purpose |
|---|---|
| `spatial_recurrence_analysis.R` | Single script that runs the complete analysis, from the master monitoring workbook to every table and figure |
| `data/Final.xlsx` | Master monitoring workbook (see *Data availability*) |
| `data/trap_captures_long.csv` | Derived, de-identified trap-level dataset (see *Data availability*); used automatically when the workbook is absent |
| `output/tables.xlsx` | All tables, one sheet each (generated) |
| `output/figures/` | Figures at 300 dpi (generated) |
| `output/analysis_log.txt` | Full console log of the run (generated) |
| `output/session_info.txt` | R and package versions of the run (generated) |

## How to run

From the repository root, in R:

```r
source("spatial_recurrence_analysis.R")
```

or from a shell:

```bash
Rscript spatial_recurrence_analysis.R
```

On Google Colab, open `spatial_recurrence_analysis.ipynb` (File > Upload notebook), upload `Final.xlsx` to `/content` through the file panel and run all cells. The configuration cell detects the Colab environment and writes to `/content/output`.

Missing packages are installed automatically. The negative-binomial GLMMs (section 14) take a few minutes; set `config$run_glmm <- FALSE` to skip them.

## Input

`data/Final.xlsx` has two sheets:

- **Consolidado**: one row per monitoring date (375 dates, 17 Oct 2019 to 24 Apr 2025). Columns: date, climate from INMET station A880, orchard-level FTD (flies per trap per day), number of valid traps, lure replacement and insecticide application flags, and one column per McPhail trap (`N°1` to `N°146`; trap numbers 93 to 96 are not used) with the number of *A. fraterculus* captured. The value `-1` means no reading on that date.
- **Coordenadas Geográficas**: trap name, latitude and longitude (WGS 84).

The script converts the workbook to a long table (one row per trap and date), keeps readings with a non-negative capture as valid, and projects coordinates to SIRGAS 2000 / UTM 22S (EPSG:31982) for all distance-based operations.

## What the script does

| Section | Analysis | Output sheets |
|---|---|---|
| 1 | Data preparation, trap activity by season | `trap_activity`, `traps_per_season` |
| 2 to 3 | Trap-level summaries, k-nearest-neighbour weights (k = 4), global Moran's I and LISA for the full period | `main_text_numbers` |
| 4 to 5 | Critical events (orchard FTD > 0.5), recurrence classes, capture intensity, global Moran's I and LISA during critical events, sensitivity to k | `S1_moran_by_k`, `class_summary` |
| 6 | Threshold scenarios for the recurrence classes, Kruskal-Wallis and Dunn's tests | `S2_class_distribution`, `S3_kruskal_by_scenario`, `S4_dunn_S1_and_S4`, `S5_dunn_S2`, `S6_dunn_S3`, `S7_dunn_S5` |
| 7 | Critical events cross-referenced with recorded insecticide applications | `critical_events`, `insecticide_applications` |
| 8 | Monte Carlo permutation tests (999 permutations) for global and local Moran; FDR and Bonferroni corrections for the LISA; overlap between full-period and critical-event hotspots | `moran_global`, `lisa_critical_by_criterion`, `lisa_critical_traps`, `lisa_full_by_criterion`, `hotspot_overlap` |
| 9 | Hotspot stability: leave-one-event-out and k = 3 to 8 | `hotspot_by_k`, `hotspot_loo_by_event`, `hotspot_trap_stability` |
| 10 | Fixed-distance neighbourhoods (250 to 1000 m) | `distance_bands` |
| 11 | Sensitivity to the FTD threshold defining a critical event (0.2 to 0.7) | `ftd_threshold_sensitivity` |
| 12 | Leave-one-event-out validation: traps classified on n − 1 events, classes tested on the held-out event | `loo_validation`, `loo_pooled_by_class` |
| 13 | Descriptive assessment of orchard FTD around insecticide applications | `ftd_by_application_window`, `ftd_before_after_application` |
| 14 | Negative-binomial GLMMs with random intercepts for trap and date (glmmTMB) | `glmm_coefficients`, `glmm_random_effects`, `glmm_fit` |
| 15 | Figures | `output/figures/` |
| 16 | Per-trap table joining all results | `per_trap` |

Random permutations use `set.seed(123)`.

## Data availability

The raw monitoring workbook (`Final.xlsx`) was provided by Rasip Agro and is subject to data-ownership restrictions; it is not included here and is available from the corresponding author on reasonable request and with the owner's permission. The derived, de-identified analysis dataset needed to reproduce every result is provided as `data/trap_captures_long.csv` (one row per trap and inspection date: date, trap, coordinates, number of A. fraterculus captured, validity flag, exposure interval, orchard-level FTD, number of valid traps, climatic variables from INMET station A880 and the operational flags for lure replacement and insecticide application; the value -1 in the flags means "not recorded"). All outputs in `output/` were generated from the raw workbook; the script reads the workbook when present and otherwise reproduces the analysis from the CSV.

## Software

R 4.5.1 with readxl 1.4.5, dplyr 1.2.1, tidyr 1.3.1, purrr 1.1.0, ggplot2 4.0.3, sf 1.0.21, spdep 1.4.2, FSA 0.10.1, writexl 1.5.4 and glmmTMB 1.1.14. The exact environment of each run is written to `output/session_info.txt`.

## Release

The release corresponding to the revised manuscript is tagged `v1.0`.

## License

MIT, see `LICENSE`.