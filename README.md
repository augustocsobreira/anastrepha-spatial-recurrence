# Spatial structure and recurrence of *Anastrepha fraterculus* captures in a commercial apple orchard

Analysis code and derived dataset for the manuscript *Spatial Structure and Recurrence of Anastrepha fraterculus Captures in a Commercial Apple Orchard* (Vacaria, Rio Grande do Sul, Brazil; 142 McPhail traps, 2019 to 2025), submitted to the *International Journal of Tropical Insect Science*. The release corresponding to the revised manuscript is tagged `v1.0`.

## Contents

| File | Purpose |
|---|---|
| `spatial_recurrence_analysis.R` | Single script that runs the complete analysis and writes every table and every analytical figure |
| `spatial_recurrence_analysis.ipynb` | The same script as a Google Colab notebook |
| `data/trap_captures_long.csv` | Derived, de-identified trap-level dataset (one row per trap and inspection date) |
| `data/study_area_map.png` | Study-area map (Figure 1), produced in QGIS |
| `output/tables.xlsx` | All tables, one sheet each |
| `output/figures/` | Figures 1 to 5 and S1 to S3 of the manuscript (600 dpi); Figure 1 is a copy of the QGIS map |
| `output/analysis_log.txt` | Console log of the run that produced `output/` |
| `output/session_info.txt` | R and package versions of that run |

The raw monitoring workbook (`data/Final.xlsx`) is not included (see *Data availability*). When it is absent, the script reads `data/trap_captures_long.csv` and reproduces every output.

## How to run

From the repository root, in R:

```r
source("spatial_recurrence_analysis.R")
```

or from a shell:

```bash
Rscript spatial_recurrence_analysis.R
```

On Google Colab, open `spatial_recurrence_analysis.ipynb`, upload `data/trap_captures_long.csv` to `/content` through the file panel and run all cells; outputs are written to `/content/output`.

Missing packages are installed automatically. The mixed models (section 14) take a few minutes; set `config$run_glmm <- FALSE` to skip them. Random permutations use `set.seed(123)`.

## Dataset

`data/trap_captures_long.csv` has 53,250 rows (142 traps x 375 recorded dates, 17 Oct 2019 to 24 Apr 2025; 350 of the dates are trap inspections and 25 are insecticide applications recorded without an inspection) and the columns:

| Column | Content |
|---|---|
| `date` | Inspection date |
| `trap` | Trap identifier (`N°1` to `N°146`; numbers 93 to 96 do not exist) |
| `latitude`, `longitude` | Trap coordinates (WGS 84); projected by the script to SIRGAS 2000 / UTM 22S (EPSG:31982) |
| `capture` | Number of *A. fraterculus* captured; empty (NA) when the trap was not read on that date |
| `valid_reading` | 1 if `capture` is a valid reading, 0 otherwise (7,989 rows) |
| `exposure_days` | Interval in days recorded in the workbook; the script recomputes the exposure as the number of days between consecutive inspection dates (they differ only on 30 Mar 2020) |
| `orchard_FTD` | Orchard-level flies per trap per day recorded in the workbook; the script recomputes it from the readings (identical on all critical dates) |
| `valid_traps`, `total_capture` | Number of traps read and total flies on the date |
| `precipitation_mm`, `tmax_C`, `tmean_C`, `tmin_C`, `humidity_pct`, `wind_ms` | Daily records of INMET station A880; `-1` = not recorded (treated as missing by the script) |
| `lure_replacement` | 1 = attractant replaced on the date; -1 = not recorded |
| `insecticide_application` | 1 = application recorded on the date; -1 = none recorded |

## What the script does

| Section | Analysis | Output sheets |
|---|---|---|
| 1 | Data preparation, trap activity by season | `trap_activity`, `traps_per_season` |
| 2 to 3 | Trap-level summaries, k-nearest-neighbor weights (k = 4), global Moran's I and LISA for the full period | `main_text_numbers` |
| 3b | Global Moran's I and LISA on the non-critical dates and in each growing season, with FDR correction | `moran_by_season`, `lisa_by_season` |
| 4 to 5 | Critical events (orchard FTD > 0.5), recurrence classes, capture intensity, Moran's I and LISA during critical events, sensitivity to k, common trap sets, geometry of the High-High cluster | `S1_moran_by_k`, `class_summary`, `moran_common_sets`, `cluster_geometry` |
| 6 | Threshold scenarios for the recurrence classes, Kruskal-Wallis and Dunn's tests | `S2_class_distribution`, `S3_kruskal_by_scenario`, `S4_dunn_S1_and_S4` to `S7_dunn_S5` |
| 7 | Critical events cross-referenced with the insecticide applications | `critical_events`, `insecticide_applications` |
| 8 | Monte Carlo permutation tests (999 permutations); FDR, permutation + FDR and Bonferroni criteria with the High-High trap identifiers; overlap between full-period and critical-event clusters | `moran_global`, `lisa_critical_by_criterion`, `lisa_critical_traps`, `lisa_full_by_criterion`, `hotspot_overlap` |
| 9 | Cluster stability: leave-one-event-out and k = 3 to 8 | `hotspot_by_k`, `hotspot_loo_by_event`, `hotspot_trap_stability` |
| 10 | Fixed-distance neighborhoods (250 to 1000 m) | `distance_bands` |
| 11 | Sensitivity to the FTD threshold defining a critical event (0.2 to 0.7), with the High-High traps of each threshold and their overlap with the 0.5 cluster | `ftd_threshold_sensitivity` |
| 12 | Leave-one-event-out validation: classes defined on n - 1 events and tested on the held-out event | `loo_validation`, `loo_pooled_by_class` |
| 13 | Orchard FTD around insecticide applications | `ftd_by_application_window`, `ftd_before_after_application` |
| 14 | Negative-binomial GLMMs (glmmTMB) with random intercepts for trap and date, residual diagnostics and sensitivity analyses | `glmm_coefficients`, `glmm_random_effects`, `glmm_fit`, `glmm_diagnostics`, `glmm_sensitivity` |
| 15 | Figures 2 to 5 and S1 to S3 of the manuscript (Figure 1 copied from `data/`) | `output/figures/` |
| 16 | Per-trap table joining all results | `per_trap` |

## Data availability

The raw monitoring workbook was provided by Rasip Agro and is subject to data-ownership restrictions; requests should be addressed to the corresponding author and depend on the owner's permission. The derived dataset above is sufficient to reproduce every result: running the script from the CSV yields the same tables as running it from the workbook.

## Software

R 4.5.1 with readxl 1.4.5, dplyr 1.2.1, tidyr 1.3.1, purrr 1.1.0, ggplot2 4.0.3, ggrepel 0.9.8, sf 1.0.21, spdep 1.4.2, FSA 0.10.1, writexl 1.5.4 and glmmTMB 1.1.14 (`output/session_info.txt`).

## License

MIT, see `LICENSE`.