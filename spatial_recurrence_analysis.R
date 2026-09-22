# =============================================================================
# Spatial structure and recurrence of Anastrepha fraterculus captures
# in a commercial apple orchard (Vacaria, RS, Brazil)
#
# Reproducible analysis script accompanying the manuscript.
#
# Input   data/Final.xlsx  - master monitoring workbook with two sheets:
#           "Consolidado"             one row per monitoring date; one column per
#                                     McPhail trap (N°1 ... N°146) with the number
#                                     of A. fraterculus captured (-1 = no reading),
#                                     plus orchard-level FTD, climate and
#                                     management fields
#           "Coordenadas Geográficas" trap name, latitude, longitude (WGS 84)
# Output  output/tables.xlsx          one sheet per table (see README)
#         output/figures/*.png        300 dpi
#         output/analysis_log.txt     full console log
#         output/session_info.txt     R and package versions
#
# Sections
#    0  Configuration
#    1  Data preparation (wide to long format, validity flags)
#    2  Trap-level summaries and spatial weights
#    3  Global Moran's I and local indicators (LISA), full period
#    4  Critical events, recurrence classes and capture intensity
#    5  Global Moran's I and LISA, critical events (Table S1, hotspot)
#    6  Threshold scenarios for the recurrence classes (Tables S2 to S7)
#    7  Critical events and insecticide applications
#    8  Permutation inference and multiple-testing correction
#    9  Hotspot stability (leave-one-event-out; number of neighbors)
#   10  Fixed-distance neighborhoods
#   11  Sensitivity to the FTD threshold that defines a critical event
#   12  Leave-one-event-out validation of the recurrence classification
#   13  Insecticide applications: descriptive assessment
#   14  Negative-binomial generalized linear mixed models
#   15  Figures
#   16  Export
#
# Run from the repository root:  source("spatial_recurrence_analysis.R")
# On Google Colab, upload trap_captures_long.csv (or Final.xlsx) and study_area_map.png to /content and source this script.
# =============================================================================

# ---- 0. Configuration -------------------------------------------------------
config <- list(
  master_file    = "data/Final.xlsx",
  output_dir     = "output",
  seed           = 123,
  n_permutations = 999,     # Monte Carlo permutations (global and local Moran)
  ftd_threshold  = 0.5,     # orchard-level FTD defining a critical event
  k_neighbours   = 4,       # k nearest neighbors for the main analyses
  crs_utm        = 31982,   # SIRGAS 2000 / UTM zone 22S (metres)
  run_glmm       = TRUE     # section 14 takes a few minutes; set FALSE to skip
)
if (dir.exists("/content") && (file.exists("/content/Final.xlsx") || file.exists("/content/trap_captures_long.csv"))) {   # Google Colab
  config$master_file <- "/content/Final.xlsx"
  config$output_dir  <- "/content/output"
}

packages <- c("readxl", "dplyr", "tidyr", "purrr", "ggplot2", "ggrepel", "sf", "spdep", "FSA", "writexl", "glmmTMB")
for (p in packages) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))

fig_dir <- file.path(config$output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(config$seed)
tbl <- list()                                   # every exported table is stored here
sink(file.path(config$output_dir, "analysis_log.txt"), split = TRUE)
cat("Spatial recurrence analysis -", format(Sys.time()), "\nMaster file:", config$master_file, "\n")

# ---- Helper functions -------------------------------------------------------
section <- function(title) cat("\n\n==========", title, "==========\n")

format_p <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, digits = 3, format = "f"))

# Growing season label: a season starts in August (e.g. 2020-12-03 -> "2020/21")
season_label <- function(d) {
  y <- as.integer(format(d, "%Y")); m <- as.integer(format(d, "%m"))
  start <- ifelse(m >= 8, y, y - 1)
  paste0(start, "/", substr(start + 1, 3, 4))
}

# Recurrence class from the relative frequency of participation in critical events
classify_recurrence <- function(freq, high = 0.85, medium = 0.50) {
  factor(dplyr::case_when(freq >= high ~ "High", freq >= medium ~ "Medium", freq > 0 ~ "Low", TRUE ~ "None"),
         levels = c("None", "Low", "Medium", "High"))
}

# LISA cluster label from the Moran scatterplot quadrant and a p-value vector
lisa_clusters <- function(local_moran, p, alpha = 0.05) {
  quadrant <- as.character(attr(local_moran, "quadr")$mean)
  ifelse(p < alpha, quadrant, "Not significant")
}

knn_weights <- function(coords, k) nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")

utm_coords <- function(df) {
  st_coordinates(st_transform(st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326), config$crs_utm))
}

moran_row <- function(x, weights, label) {
  test <- moran.test(x, weights, zero.policy = TRUE)
  tibble(configuration = label, moran_I = unname(test$estimate[1]), expected = unname(test$estimate[2]),
         variance = unname(test$estimate[3]), z = unname(test$statistic), p_value = test$p.value)
}

days_since <- function(dates, reference) sapply(dates, function(d) { p <- reference[reference < d]; if (length(p)) as.integer(d - max(p)) else NA_integer_ })
days_until <- function(dates, reference) sapply(dates, function(d) { p <- reference[reference > d]; if (length(p)) as.integer(min(p) - d) else NA_integer_ })

# Per-trap recurrence: proportion of the critical dates on which the trap was READ that had at least one
# capture. Traps without any reading on the critical dates are not classified (absence of observation is
# not a zero capture), and traps read on fewer dates use their own number of readings as denominator.
recurrence_table <- function(obs_critical, n_events = NULL, all_traps = NULL) {
  obs_critical %>% group_by(Armadilha, Latitude, Longitude) %>%
    summarise(n_events_read = n(), n_events_positive = sum(Capturas > 0), .groups = "drop") %>%
    mutate(rel_frequency = n_events_positive / n_events_read)
}

# Per-trap capture intensity during critical events (all readings, including zeros)
intensity_table <- function(obs_critical) {
  obs_critical %>% group_by(Armadilha, Latitude, Longitude) %>%
    summarise(mean_capture_critical = mean(Capturas), total_capture_critical = sum(Capturas),
              n_readings_critical = n(), .groups = "drop")
}

# ---- 1. Data preparation ----------------------------------------------------
section("1. Data preparation")

build_long_table <- function(master_file) {
  wide <- read_excel(master_file, sheet = "Consolidado", .name_repair = "minimal")
  names(wide) <- gsub("[[:space:]]+", " ", names(wide))
  wide <- wide %>% filter(!is.na(Data)) %>% mutate(Data = as.Date(Data))
  trap_cols <- grep("^N°[0-9]+$", names(wide), value = TRUE)

  coords <- read_excel(master_file, sheet = "Coordenadas Geográficas", .name_repair = "minimal")[, 1:3]
  names(coords) <- c("name", "Latitude", "Longitude")
  coords <- coords %>% filter(!is.na(name)) %>%
    mutate(Armadilha = paste0("N°", as.integer(gsub("[^0-9]", "", name))))
  stopifnot(!any(duplicated(coords$Armadilha)))          # one coordinate pair per trap

  long <- wide %>%
    select(Data, all_of(trap_cols), `Diferença em dias`, FTD, `Número de armadilhas válidas`, `Total de capturas`,
           `Precipitação (mm)`, `Tmax (°C)`, `Tmed (°C)`, `Tmin (°C)`, `Umidade média (%)`, `Vento (m/s)`,
           `Troca de atrativos`, `Aplicação de inseticida`) %>%
    pivot_longer(all_of(trap_cols), names_to = "Armadilha", values_to = "capture_raw") %>%
    mutate(capture_raw = suppressWarnings(as.numeric(as.character(capture_raw)))) %>%
    left_join(coords %>% select(Armadilha, Latitude, Longitude), by = "Armadilha") %>%
    transmute(Data, Armadilha, Latitude, Longitude,
              Capturas      = ifelse(!is.na(capture_raw) & capture_raw >= 0, capture_raw, NA_real_),
              valid_reading = as.integer(!is.na(capture_raw) & capture_raw >= 0),   # -1 or blank = no reading
              interval_days = `Diferença em dias`,
              orchard_FTD   = FTD,
              valid_traps_sheet = `Número de armadilhas válidas`,
              total_capture_sheet = `Total de capturas`,
              precipitation_mm = `Precipitação (mm)`, tmax_C = `Tmax (°C)`, tmean_C = `Tmed (°C)`, tmin_C = `Tmin (°C)`,
              humidity_pct = `Umidade média (%)`, wind_ms = `Vento (m/s)`,
              lure_change = ifelse(`Troca de atrativos` == 1, 1L, -1L), insecticide = ifelse(`Aplicação de inseticida` %in% c(1, 2), 1L, -1L)) %>%   # 1 = recorded, -1 = not recorded
    arrange(Data, Armadilha)
  stopifnot(!any(is.na(long$Latitude)))                  # every trap has coordinates
  cat("Master sheet:", nrow(wide), "monitoring dates x", length(trap_cols), "traps\n")
  long
}

read_long_csv <- function(csv) {
  d <- read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8", fileEncoding = "UTF-8")
  data.frame(Data = as.Date(d$date), Armadilha = d$trap, Latitude = d$latitude, Longitude = d$longitude, Capturas = d$capture,
             valid_reading = d$valid_reading, interval_days = d$exposure_days, orchard_FTD = d$orchard_FTD,
             valid_traps_sheet = d$valid_traps, total_capture_sheet = d$total_capture, precipitation_mm = d$precipitation_mm,
             tmax_C = d$tmax_C, tmean_C = d$tmean_C, tmin_C = d$tmin_C, humidity_pct = d$humidity_pct, wind_ms = d$wind_ms,
             lure_change = d$lure_replacement, insecticide = d$insecticide_application, stringsAsFactors = FALSE)
}
csv_file <- file.path(dirname(config$master_file), "trap_captures_long.csv")
if (file.exists(config$master_file)) {
  records <- build_long_table(config$master_file)
} else if (file.exists(csv_file)) {
  cat("Master workbook not found; reading the derived dataset", csv_file, "\n")
  records <- read_long_csv(csv_file)
} else stop("Neither the master workbook nor data/trap_captures_long.csv was found")

# Exposure interval = days between consecutive inspection dates (dates with at least one valid reading);
# orchard FTD = total captures / (traps with a valid reading x exposure). Negative climate values are
# missing-value codes of the workbook and are treated as missing.
inspections <- records %>% filter(valid_reading == 1) %>% group_by(Data) %>%
  summarise(n_valid = n(), flies = sum(Capturas), .groups = "drop") %>% arrange(Data) %>%
  mutate(gap_days = as.integer(Data - lag(Data)))
records <- records %>% left_join(inspections, by = "Data") %>%
  mutate(interval_days = ifelse(!is.na(gap_days), gap_days, interval_days),
         orchard_FTD   = ifelse(!is.na(n_valid), flies / (n_valid * pmax(interval_days, 1)), orchard_FTD)) %>%
  select(-n_valid, -flies, -gap_days) %>%
  mutate(across(c(precipitation_mm, tmax_C, tmean_C, tmin_C, humidity_pct, wind_ms), ~ ifelse(.x < 0, NA_real_, .x)))
obs <- records %>% filter(valid_reading == 1)
stopifnot(all(table(obs$Armadilha, obs$Data) <= 1))
cat("Trap-date records:", nrow(records), "| valid readings:", nrow(obs),
    "| traps:", n_distinct(obs$Armadilha), "| dates:", n_distinct(obs$Data), "\n")

# Trap activity over the study (traps added or deactivated during monitoring)
tbl$trap_activity <- obs %>% group_by(trap = Armadilha) %>%
  summarise(first_date = min(Data), last_date = max(Data), n_dates = n(), .groups = "drop") %>%
  mutate(season_in = season_label(first_date), season_out = season_label(last_date))
cat("Traps by number of valid readings:\n"); print(table(tbl$trap_activity$n_dates))
tbl$traps_per_season <- obs %>% count(Data, name = "valid_traps") %>% mutate(season = season_label(Data)) %>%
  group_by(season) %>% summarise(dates = n(), min_traps = min(valid_traps), max_traps = max(valid_traps), .groups = "drop")
print(as.data.frame(tbl$traps_per_season))

# ---- 2. Trap-level summaries and spatial weights ----------------------------
section("2. Trap-level summaries and spatial weights")
trap_summary <- obs %>% group_by(Armadilha) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude),
            mean_capture = mean(Capturas), median_capture = median(Capturas), total_capture = sum(Capturas),
            sd_capture = sd(Capturas), n_readings = n(), prop_positive = mean(Capturas > 0), .groups = "drop") %>%
  mutate(cv_capture = ifelse(mean_capture > 0, sd_capture / mean_capture, NA))
trap_coords <- utm_coords(trap_summary)
w_full <- knn_weights(trap_coords, config$k_neighbours)
nn_distance <- apply(as.matrix(dist(trap_coords)), 1, function(x) mean(sort(x)[2:(config$k_neighbours + 1)]))
cat(sprintf("Mean distance to the %d nearest neighbors (m): mean %.1f, min %.1f, max %.1f\n",
            config$k_neighbours, mean(nn_distance), min(nn_distance), max(nn_distance)))

# ---- 3. Global Moran's I and LISA, full period ------------------------------
section("3. Global Moran's I and LISA - full period")
moran_full <- moran.test(trap_summary$mean_capture, w_full)
print(moran_full)
lisa_full_raw <- localmoran(trap_summary$mean_capture, w_full)
trap_summary <- trap_summary %>%
  mutate(Ii_full = lisa_full_raw[, 1], p_full = lisa_full_raw[, 5],
         cluster_full = lisa_clusters(lisa_full_raw, p_full))
print(table(trap_summary$cluster_full))

# ---- 3b. Spatial structure outside the critical dates and by season ---------
section("3b. Global Moran's I on non-critical dates and by season")
critical_dates_pre <- obs %>% filter(orchard_FTD > config$ftd_threshold) %>% distinct(Data) %>% pull(Data)
noncrit_summary <- obs %>% filter(!Data %in% critical_dates_pre) %>% group_by(Armadilha) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude), mean_capture = mean(Capturas), .groups = "drop")
w_nc <- knn_weights(utm_coords(noncrit_summary), config$k_neighbours)
moran_noncrit <- moran.test(noncrit_summary$mean_capture, w_nc)
mc_noncrit <- moran.mc(noncrit_summary$mean_capture, w_nc, nsim = config$n_permutations)
lm_nc <- localmoran(noncrit_summary$mean_capture, w_nc)
noncrit_summary$cluster <- lisa_clusters(lm_nc, lm_nc[, 5])
noncrit_summary$cluster_fdr <- ifelse(p.adjust(lm_nc[, 5], "BH") < 0.05, as.character(attr(lm_nc, "quadr")$mean), "Not significant")
noncrit_hh <- noncrit_summary$Armadilha[noncrit_summary$cluster == "High-High"]
cat("Non-critical High-High traps (analytical p < 0.05):", paste(noncrit_hh, collapse = ", "),
    "| surviving FDR:", sum(noncrit_summary$cluster_fdr == "High-High"), "\n")
cat("Non-critical dates only:", n_distinct(obs$Data) - length(critical_dates_pre), "dates,", nrow(noncrit_summary), "traps; Moran I =",
    round(moran_noncrit$estimate[1], 4), "analytical p =", signif(moran_noncrit$p.value, 3), "permutation p =", mc_noncrit$p.value, "\n")
share_critical <- sum(obs$Capturas[obs$Data %in% critical_dates_pre]) / sum(obs$Capturas)
cat("Share of all flies captured on the critical dates:", round(100 * share_critical, 1), "%\n")
season_lisa <- list()
tbl$moran_by_season <- map_dfr(sort(unique(season_label(obs$Data))), function(s) {
  d <- obs %>% filter(season_label(Data) == s) %>% group_by(Armadilha) %>%
    summarise(Latitude = first(Latitude), Longitude = first(Longitude), mean_capture = mean(Capturas), n = n(), .groups = "drop")
  w <- knn_weights(utm_coords(d), config$k_neighbours)
  t <- moran.test(d$mean_capture, w); mc <- moran.mc(d$mean_capture, w, nsim = config$n_permutations)
  lm <- localmoran(d$mean_capture, w); cl <- lisa_clusters(lm, lm[, 5]); hh <- sum(cl == "High-High")
  hh_fdr <- sum(cl == "High-High" & p.adjust(lm[, 5], "BH") < 0.05)
  season_lisa[[s]] <<- d %>% mutate(season = s, cluster = cl)
  tibble(season = s, dates = n_distinct(obs$Data[season_label(obs$Data) == s]), traps = nrow(d),
         critical_events = sum(critical_dates_pre %in% obs$Data[season_label(obs$Data) == s]),
         total_flies = sum(obs$Capturas[season_label(obs$Data) == s]),
         moran_I = round(unname(t$estimate[1]), 4), p_analytical = t$p.value, p_permutation = mc$p.value, high_high_traps = hh,
         high_high_traps_fdr = hh_fdr)
})
tbl$moran_by_season <- bind_rows(tbl$moran_by_season,
  tibble(season = "All non-critical dates", dates = n_distinct(obs$Data) - length(critical_dates_pre), traps = nrow(noncrit_summary),
         critical_events = 0L, total_flies = sum(obs$Capturas[!obs$Data %in% critical_dates_pre]),
         moran_I = round(unname(moran_noncrit$estimate[1]), 4), p_analytical = moran_noncrit$p.value, p_permutation = mc_noncrit$p.value,
         high_high_traps = { lm <- localmoran(noncrit_summary$mean_capture, w_nc); sum(lisa_clusters(lm, lm[, 5]) == "High-High") },
         high_high_traps_fdr = sum(noncrit_summary$cluster_fdr == "High-High")))
print(as.data.frame(tbl$moran_by_season))

# ---- 4. Critical events, recurrence classes and capture intensity -----------
section("4. Critical events and recurrence classes")
critical_dates <- obs %>% filter(orchard_FTD > config$ftd_threshold) %>% distinct(Data) %>% arrange(Data) %>% pull(Data)
n_events <- length(critical_dates)
cat("Critical events (orchard FTD >", config$ftd_threshold, "):", n_events, "\n"); print(format(critical_dates))
obs_critical <- obs %>% filter(Data %in% critical_dates)
all_traps <- obs %>% distinct(Armadilha, Latitude, Longitude)

recurrence <- recurrence_table(obs_critical) %>%
  mutate(recurrence_class = classify_recurrence(rel_frequency))
unobserved_traps <- setdiff(all_traps$Armadilha, recurrence$Armadilha)
cat("Traps read on at least one critical date (classified):", nrow(recurrence),
    "| traps without any reading on a critical date (not classified):", length(unobserved_traps), "\n")
cat("Traps by number of critical dates read:\n"); print(table(recurrence$n_events_read))
cat("Relative frequency distribution:\n"); print(table(round(recurrence$rel_frequency, 3)))
cat("Recurrence classes:\n"); print(table(recurrence$recurrence_class))

intensity_critical <- intensity_table(obs_critical)
cat("Traps with at least one reading during critical events:", nrow(intensity_critical), "\n")
trap_data <- recurrence %>%
  left_join(intensity_critical, by = c("Armadilha", "Latitude", "Longitude")) %>%
  mutate(across(c(mean_capture_critical, total_capture_critical, n_readings_critical), ~ replace_na(.x, 0)))

# ---- 5. Global Moran's I and LISA, critical events --------------------------
section("5. Global Moran's I and LISA - critical events")
crit_coords <- utm_coords(intensity_critical)
w_crit <- knn_weights(crit_coords, config$k_neighbours)
x_crit <- intensity_critical$mean_capture_critical
moran_crit <- moran.test(x_crit, w_crit)
print(moran_crit)

# Table S1: sensitivity of the global index to the number of neighbors
tbl$S1_moran_by_k <- map_dfr(3:8, function(k) moran_row(x_crit, knn_weights(crit_coords, k), paste0("k = ", k))) %>%
  mutate(across(c(moran_I, expected, variance, z), ~ round(.x, 4)))
cat("Table S1:\n"); print(as.data.frame(tbl$S1_moran_by_k))

lisa_crit_raw <- localmoran(x_crit, w_crit)
intensity_critical <- intensity_critical %>%
  mutate(Ii = lisa_crit_raw[, 1], p_analytical = lisa_crit_raw[, 5],
         cluster = lisa_clusters(lisa_crit_raw, p_analytical))
print(table(intensity_critical$cluster))
cat("High-High traps (analytical p < 0.05):",
    paste(intensity_critical$Armadilha[intensity_critical$cluster == "High-High"], collapse = ", "), "\n")

# Sensitivity to unequal sampling histories: traps observed over the whole period / on all seven critical dates
hh_crit_ids <- intensity_critical$Armadilha[intensity_critical$cluster == "High-High"]
cl_pts <- st_transform(st_as_sf(intensity_critical %>% filter(Armadilha %in% hh_crit_ids), coords = c("Longitude", "Latitude"), crs = 4326), config$crs_utm)
cl_xy <- st_coordinates(cl_pts); all_xy <- utm_coords(all_traps)
nn_cl <- apply(as.matrix(dist(cl_xy)), 1, function(d) min(d[d > 0]))
in_cl <- intensity_critical$Armadilha %in% hh_crit_ids
tbl$cluster_geometry <- tibble(
  item = c("High-High traps on the critical dates", "longitude, westernmost trap", "longitude, easternmost trap", "latitude, southernmost trap",
           "latitude, northernmost trap", "convex hull area (ha)", "distance from the westernmost cluster trap to the western limit of the network (m)",
           "nearest High-High neighbor, minimum distance (m)", "nearest High-High neighbor, maximum distance (m)", "nearest High-High neighbor, mean distance (m)",
           "mean of the per-trap mean captures on the critical dates, cluster traps", "mean of the per-trap mean captures on the critical dates, other traps",
           "other traps with a reading on the critical dates"),
  value = c(length(hh_crit_ids), min(intensity_critical$Longitude[in_cl]), max(intensity_critical$Longitude[in_cl]),
            min(intensity_critical$Latitude[in_cl]), max(intensity_critical$Latitude[in_cl]),
            as.numeric(st_area(st_convex_hull(st_union(cl_pts)))) / 1e4, min(cl_xy[, 1]) - min(all_xy[, 1]),
            min(nn_cl), max(nn_cl), mean(nn_cl),
            mean(intensity_critical$mean_capture_critical[in_cl]), mean(intensity_critical$mean_capture_critical[!in_cl]), sum(!in_cl)))
cat("Geometry of the critical-event cluster:\n"); print(as.data.frame(tbl$cluster_geometry), digits = 5)
seasons_per_trap <- obs %>% group_by(Armadilha) %>% summarise(n_seasons = n_distinct(season_label(Data)), .groups = "drop")
common_full <- trap_summary %>% filter(Armadilha %in% seasons_per_trap$Armadilha[seasons_per_trap$n_seasons == max(seasons_per_trap$n_seasons)])
common_crit <- intensity_critical %>% filter(n_readings_critical == n_events)
moran_common <- function(d, v, label) { w <- knn_weights(utm_coords(d), config$k_neighbours); t <- moran.test(d[[v]], w); mc <- moran.mc(d[[v]], w, nsim = config$n_permutations)
  lm <- localmoran(d[[v]], w); cl <- lisa_clusters(lm, lm[, 5])
  tibble(subset = label, traps = nrow(d), moran_I = round(unname(t$estimate[1]), 4), p_analytical = t$p.value, p_permutation = mc$p.value,
         high_high_traps = sum(cl == "High-High"), high_high_in_critical_cluster = sum(d$Armadilha[cl == "High-High"] %in% hh_crit_ids)) }
tbl$moran_common_sets <- bind_rows(
  moran_common(common_full, "mean_capture", "Full period, traps present in all six seasons"),
  moran_common(common_crit, "mean_capture_critical", "Critical events, traps read on all seven dates"))
print(as.data.frame(tbl$moran_common_sets))
season_overlap <- vapply(names(season_lisa), function(s) sum(season_lisa[[s]]$Armadilha[season_lisa[[s]]$cluster == "High-High"] %in% hh_crit_ids), 1L)
tbl$moran_by_season$high_high_in_critical_cluster <- c(season_overlap[tbl$moran_by_season$season[seq_along(season_overlap)]], sum(noncrit_hh %in% hh_crit_ids))
tbl$moran_by_season$high_high_traps_ids <- c(vapply(names(season_lisa), function(s) paste(season_lisa[[s]]$Armadilha[season_lisa[[s]]$cluster == "High-High"], collapse = ", "), ""), paste(noncrit_hh, collapse = ", "))
tbl$lisa_by_season <- bind_rows(season_lisa) %>% select(season, trap = Armadilha, Latitude, Longitude, mean_capture, cluster) %>% mutate(in_critical_cluster = trap %in% hh_crit_ids)
cat("Seasonal High-High traps inside the critical-event cluster:\n"); print(as.data.frame(tbl$moran_by_season %>% select(season, high_high_traps, high_high_in_critical_cluster, high_high_traps_ids)))

# Capture intensity by recurrence class (baseline thresholds 0.85 / 0.50)
section("5b. Capture intensity by recurrence class - baseline")
kw_baseline <- kruskal.test(mean_capture_critical ~ recurrence_class, data = trap_data)
print(kw_baseline)
dunn_baseline <- dunnTest(mean_capture_critical ~ recurrence_class, data = trap_data, method = "bonferroni")
print(dunn_baseline$res)
tbl$class_summary <- trap_data %>% group_by(recurrence_class) %>%
  summarise(n_traps = n(), median_capture = median(mean_capture_critical), mean_capture = mean(mean_capture_critical),
            sd_capture = sd(mean_capture_critical), median_rel_frequency = median(rel_frequency), .groups = "drop")
print(as.data.frame(tbl$class_summary))

# ---- 6. Threshold scenarios for the recurrence classes (Tables S2 to S7) ----
section("6. Threshold scenarios for the recurrence classes")
scenarios <- tribble(
  ~scenario, ~description, ~high_threshold, ~medium_threshold,
  "S1 - Baseline",                        "Original operational thresholds",                        0.85, 0.50,
  "S2 - More conservative",               "Stricter classification for high recurrence",            0.90, 0.60,
  "S3 - Less conservative",               "More permissive classification thresholds",              0.80, 0.40,
  "S4 - Moderate (lower high threshold)", "Reduced high threshold with same medium threshold",      0.75, 0.50,
  "S5 - Moderate (higher high threshold)","Increased high threshold with unchanged medium threshold", 0.90, 0.50)
scenario_data <- scenarios %>%
  mutate(data = map2(high_threshold, medium_threshold,
                     ~ trap_data %>% mutate(recurrence_class = classify_recurrence(rel_frequency, .x, .y))))

tbl$S2_class_distribution <- scenario_data %>%
  mutate(d = map(data, ~ as.data.frame(table(.x$recurrence_class)) %>% pivot_wider(names_from = Var1, values_from = Freq))) %>%
  select(scenario, d) %>% unnest(d)
print(as.data.frame(tbl$S2_class_distribution))

tbl$S3_kruskal_by_scenario <- scenario_data %>%
  mutate(kw = map(data, ~ kruskal.test(mean_capture_critical ~ recurrence_class, data = .x)),
         kruskal_wallis_chi2 = round(map_dbl(kw, ~ unname(.x$statistic)), 2),
         df = map_dbl(kw, ~ unname(.x$parameter)), p_value = map_dbl(kw, ~ .x$p.value)) %>%
  select(scenario, description, high_threshold, medium_threshold, kruskal_wallis_chi2, df, p_value)
print(as.data.frame(tbl$S3_kruskal_by_scenario))

dunn_by_scenario <- scenario_data %>%
  mutate(dunn = map(data, ~ dunnTest(mean_capture_critical ~ recurrence_class, data = .x, method = "bonferroni")$res)) %>%
  select(scenario, dunn) %>% unnest(dunn) %>%
  transmute(scenario, comparison = Comparison, Z = round(Z, 4), p_unadjusted = P.unadj,
            p_adjusted = P.adj, p_adjusted_text = format_p(P.adj))
tbl$S4_dunn_S1_and_S4 <- dunn_by_scenario %>% filter(scenario == "S1 - Baseline")
tbl$S5_dunn_S2 <- dunn_by_scenario %>% filter(scenario == "S2 - More conservative")
tbl$S6_dunn_S3 <- dunn_by_scenario %>% filter(scenario == "S3 - Less conservative")
tbl$S7_dunn_S5 <- dunn_by_scenario %>% filter(scenario == "S5 - Moderate (higher high threshold)")
cat("Dunn's test, all scenarios:\n"); print(as.data.frame(dunn_by_scenario %>% select(-p_adjusted)))

# ---- 7. Critical events and insecticide applications ------------------------
section("7. Critical events and insecticide applications")
date_level <- records %>% group_by(Data) %>%
  summarise(orchard_FTD = first(orchard_FTD), valid_traps_sheet = first(valid_traps_sheet),
            total_capture_sheet = first(total_capture_sheet), insecticide = first(insecticide),
            interval_days = first(interval_days), .groups = "drop") %>% arrange(Data) %>%
  mutate(has_reading = valid_traps_sheet > 0 & total_capture_sheet >= 0)   # -1 in the sheet = no inspection on that date
application_dates <- date_level %>% filter(insecticide == 1) %>% pull(Data)
cat("Recorded insecticide applications:", length(application_dates), "\n")

tbl$critical_events <- obs_critical %>% group_by(date = Data) %>%
  summarise(orchard_FTD = round(first(orchard_FTD), 3), exposure_days = first(interval_days),
            valid_traps_sheet = first(valid_traps_sheet), traps_with_reading = n(),
            traps_positive = sum(Capturas > 0), total_flies = sum(Capturas), max_trap_capture = max(Capturas),
            mean_temp_C = first(tmean_C), precipitation_mm = first(precipitation_mm), .groups = "drop") %>%
  mutate(season = season_label(date), application_same_day = date %in% application_dates,
         days_since_previous_application = days_since(date, application_dates),
         days_to_next_application = days_until(date, application_dates)) %>%
  select(date, season, everything())
print(as.data.frame(tbl$critical_events))

tbl$insecticide_applications <- date_level %>% filter(Data %in% application_dates) %>%
  transmute(date = Data, season = season_label(Data), code_in_sheet = insecticide,
            reading_same_day = has_reading,
            orchard_FTD = ifelse(has_reading, round(orchard_FTD, 3), NA),
            critical_event_same_day = Data %in% critical_dates,
            days_since_last_critical_event = days_since(Data, critical_dates))
print(as.data.frame(tbl$insecticide_applications))

# ---- 8. Permutation inference and multiple-testing correction ---------------
section("8. Permutation inference and multiple-testing correction")
mc_full <- moran.mc(trap_summary$mean_capture, w_full, nsim = config$n_permutations)
mc_crit <- moran.mc(x_crit, w_crit, nsim = config$n_permutations)
tbl$moran_global <- tibble(
  dataset = c("Full period (all dates)", "Critical events only"),
  n_traps = c(nrow(trap_summary), nrow(intensity_critical)),
  moran_I = round(c(moran_full$estimate[1], moran_crit$estimate[1]), 4),
  p_analytical = c(moran_full$p.value, moran_crit$p.value),
  p_permutation = c(mc_full$p.value, mc_crit$p.value),
  n_permutations = config$n_permutations)
print(as.data.frame(tbl$moran_global))

perm_column <- function(m) m[, grep("Sim", colnames(m), value = TRUE)[1]]
lisa_crit_perm <- localmoran_perm(x_crit, w_crit, nsim = config$n_permutations)
lisa_crit <- intensity_critical %>%
  mutate(p_permutation = perm_column(lisa_crit_perm),
         q_fdr = p.adjust(p_analytical, "BH"), p_bonferroni = p.adjust(p_analytical, "bonferroni"),
         q_fdr_permutation = p.adjust(p_permutation, "BH"),
         quadrant = as.character(attr(lisa_crit_raw, "quadr")$mean),
         cluster_permutation = ifelse(p_permutation < 0.05, quadrant, "Not significant"),
         cluster_fdr = ifelse(q_fdr < 0.05, quadrant, "Not significant"),
         cluster_permutation_fdr = ifelse(q_fdr_permutation < 0.05, quadrant, "Not significant"),
         cluster_bonferroni = ifelse(p_bonferroni < 0.05, quadrant, "Not significant"))
count_hh <- function(v) sum(v == "High-High"); count_sig <- function(v) sum(v != "Not significant")
tbl$lisa_critical_by_criterion <- tibble(
  criterion = c("Analytical p < 0.05", "Permutation p < 0.05", "FDR (Benjamini-Hochberg) q < 0.05",
                "Permutation + FDR q < 0.05", "Bonferroni p < 0.05"),
  high_high = c(count_hh(lisa_crit$cluster), count_hh(lisa_crit$cluster_permutation), count_hh(lisa_crit$cluster_fdr),
                count_hh(lisa_crit$cluster_permutation_fdr), count_hh(lisa_crit$cluster_bonferroni)),
  any_significant = c(count_sig(lisa_crit$cluster), count_sig(lisa_crit$cluster_permutation), count_sig(lisa_crit$cluster_fdr),
                      count_sig(lisa_crit$cluster_permutation_fdr), count_sig(lisa_crit$cluster_bonferroni)),
  high_high_ids = vapply(c("cluster", "cluster_permutation", "cluster_fdr", "cluster_permutation_fdr", "cluster_bonferroni"),
                         function(v) paste(lisa_crit$Armadilha[lisa_crit[[v]] == "High-High"], collapse = ", "), ""))
print(as.data.frame(tbl$lisa_critical_by_criterion))
tbl$lisa_critical_traps <- lisa_crit %>%
  filter(cluster != "Not significant" | cluster_permutation != "Not significant") %>%
  transmute(trap = Armadilha, mean_capture_critical = round(mean_capture_critical, 3), Ii = round(Ii, 4), quadrant,
            p_analytical = round(p_analytical, 4), p_permutation = round(p_permutation, 4),
            q_fdr = round(q_fdr, 4), p_bonferroni = round(p_bonferroni, 4)) %>% arrange(p_analytical)
print(as.data.frame(tbl$lisa_critical_traps))

lisa_full_perm <- localmoran_perm(trap_summary$mean_capture, w_full, nsim = config$n_permutations)
lisa_full <- trap_summary %>%
  mutate(p_permutation = perm_column(lisa_full_perm), q_fdr = p.adjust(p_full, "BH"), p_bonferroni = p.adjust(p_full, "bonferroni"),
         q_fdr_permutation = p.adjust(p_permutation, "BH"),
         quadrant = as.character(attr(lisa_full_raw, "quadr")$mean),
         cluster_permutation = ifelse(p_permutation < 0.05, quadrant, "Not significant"),
         cluster_fdr = ifelse(q_fdr < 0.05, quadrant, "Not significant"),
         cluster_permutation_fdr = ifelse(q_fdr_permutation < 0.05, quadrant, "Not significant"),
         cluster_bonferroni = ifelse(p_bonferroni < 0.05, quadrant, "Not significant"))
full_criteria <- c("cluster_full", "cluster_permutation", "cluster_fdr", "cluster_permutation_fdr", "cluster_bonferroni")
tbl$lisa_full_by_criterion <- tibble(
  criterion = c("Analytical p < 0.05", "Permutation p < 0.05", "FDR (Benjamini-Hochberg) q < 0.05", "Permutation + FDR q < 0.05", "Bonferroni p < 0.05"),
  high_high = vapply(full_criteria, function(v) count_hh(lisa_full[[v]]), 1L),
  low_low = vapply(full_criteria, function(v) sum(lisa_full[[v]] == "Low-Low"), 1L),
  any_significant = vapply(full_criteria, function(v) count_sig(lisa_full[[v]]), 1L),
  high_high_ids = vapply(full_criteria, function(v) paste(lisa_full$Armadilha[lisa_full[[v]] == "High-High"], collapse = ", "), ""))
cat("LISA, full period, by criterion:\n"); print(as.data.frame(tbl$lisa_full_by_criterion))

hh_full <- lisa_full$Armadilha[lisa_full$cluster_full == "High-High"]
hh_crit <- lisa_crit$Armadilha[lisa_crit$cluster == "High-High"]
tbl$hotspot_overlap <- lisa_full %>% select(trap = Armadilha, cluster_full, q_fdr_full = q_fdr) %>%
  full_join(lisa_crit %>% select(trap = Armadilha, cluster_critical = cluster, q_fdr_critical = q_fdr), by = "trap") %>%
  left_join(recurrence %>% select(trap = Armadilha, recurrence_class, rel_frequency), by = "trap") %>%
  filter(cluster_full == "High-High" | cluster_critical %in% "High-High") %>%
  mutate(across(c(q_fdr_full, q_fdr_critical, rel_frequency), ~ round(.x, 4))) %>%
  arrange(desc(cluster_full == "High-High" & cluster_critical %in% "High-High"))
cat("High-High traps in the full period and/or in critical events:\n"); print(as.data.frame(tbl$hotspot_overlap))
cat("High-High full period:", length(hh_full), "| critical events:", length(hh_crit), "| both:", length(intersect(hh_full, hh_crit)), "\n")

# ---- 9. Hotspot stability ---------------------------------------------------
section("9. Hotspot stability")
lisa_on_subset <- function(obs_subset, k = config$k_neighbours) {
  r <- intensity_table(obs_subset)
  lm <- localmoran(r$mean_capture_critical, knn_weights(utm_coords(r), k))
  r$cluster <- lisa_clusters(lm, lm[, 5]); r
}
hh_by_k <- map_dfr(3:8, function(k) { r <- lisa_on_subset(obs_critical, k); tibble(k, trap = r$Armadilha, cluster = r$cluster) })
tbl$hotspot_by_k <- hh_by_k %>% group_by(k) %>%
  summarise(high_high = sum(cluster == "High-High"), low_high = sum(cluster == "Low-High"),
            high_low = sum(cluster == "High-Low"), low_low = sum(cluster == "Low-Low"), .groups = "drop")
print(as.data.frame(tbl$hotspot_by_k))

hh_loo <- map_dfr(seq_len(n_events), function(i) {
  r <- lisa_on_subset(obs_critical %>% filter(Data != critical_dates[i]))
  tibble(event_removed = critical_dates[i], trap = r$Armadilha, cluster = r$cluster, n_high_high = sum(r$cluster == "High-High"))
})
tbl$hotspot_loo_by_event <- hh_loo %>% distinct(event_removed, n_high_high)
print(as.data.frame(tbl$hotspot_loo_by_event))
tbl$hotspot_trap_stability <- lisa_crit %>% select(trap = Armadilha, cluster_all_events = cluster) %>%
  left_join(hh_loo %>% filter(cluster == "High-High") %>% count(trap, name = "times_high_high_loo"), by = "trap") %>%
  left_join(hh_by_k %>% filter(cluster == "High-High") %>% count(trap, name = "times_high_high_k3_to_k8"), by = "trap") %>%
  mutate(across(c(times_high_high_loo, times_high_high_k3_to_k8), ~ replace_na(.x, 0L))) %>%
  filter(cluster_all_events == "High-High" | times_high_high_loo > 0 | times_high_high_k3_to_k8 > 0) %>%
  arrange(desc(times_high_high_loo), desc(times_high_high_k3_to_k8))
cat("High-High traps across leave-one-event-out runs (of", n_events, ") and k = 3..8 (of 6):\n")
print(as.data.frame(tbl$hotspot_trap_stability))

# ---- 10. Fixed-distance neighborhoods --------------------------------------
section("10. Fixed-distance neighborhoods")
radii_m <- c(250, 300, 400, 500, 750, 1000)
tbl$distance_bands <- map_dfr(radii_m, function(d) {
  nb_c <- dnearneigh(crit_coords, 0, d); w_c <- nb2listw(nb_c, style = "W", zero.policy = TRUE)
  t_c <- moran.test(x_crit, w_c, zero.policy = TRUE)
  mc_c <- moran.mc(x_crit, w_c, nsim = config$n_permutations, zero.policy = TRUE)
  nb_f <- dnearneigh(trap_coords, 0, d); w_f <- nb2listw(nb_f, style = "W", zero.policy = TRUE)
  t_f <- moran.test(trap_summary$mean_capture, w_f, zero.policy = TRUE)
  tibble(radius_m = d, islands = sum(card(nb_c) == 0), mean_links = round(mean(card(nb_c)), 2),
         moran_I_critical = round(unname(t_c$estimate[1]), 4), p_critical = t_c$p.value, p_permutation_critical = mc_c$p.value,
         moran_I_full = round(unname(t_f$estimate[1]), 4), p_full = t_f$p.value)
})
print(as.data.frame(tbl$distance_bands))

# ---- 11. Sensitivity to the FTD threshold -----------------------------------
section("11. Sensitivity to the FTD threshold defining critical events")
ftd_thresholds <- c(0.2, 0.3, 0.4, 0.5, 0.7)
tbl$ftd_threshold_sensitivity <- map_dfr(ftd_thresholds, function(th) {
  dates_th <- obs %>% filter(orchard_FTD > th) %>% distinct(Data) %>% pull(Data)
  obs_th <- obs %>% filter(Data %in% dates_th)
  r <- intensity_table(obs_th)
  w <- knn_weights(utm_coords(r), config$k_neighbours)
  t <- moran.test(r$mean_capture_critical, w); mc <- moran.mc(r$mean_capture_critical, w, nsim = config$n_permutations)
  lm <- localmoran(r$mean_capture_critical, w); cl_th <- lisa_clusters(lm, lm[, 5]); hh <- sum(cl_th == "High-High")
  hh_ids_th <- r$Armadilha[cl_th == "High-High"]
  cl <- recurrence_table(obs_th, length(dates_th), all_traps) %>% mutate(class = classify_recurrence(rel_frequency)) %>%
    left_join(r, by = c("Armadilha", "Latitude", "Longitude")) %>% mutate(mean_capture_critical = replace_na(mean_capture_critical, 0))
  kw <- kruskal.test(mean_capture_critical ~ class, data = cl)
  tibble(ftd_threshold = th, n_events = length(dates_th), n_traps = nrow(r), moran_I = round(unname(t$estimate[1]), 4),
         p_analytical = t$p.value, p_permutation = mc$p.value, high_high_traps = hh,
         high_high_in_baseline_cluster = sum(hh_ids_th %in% hh_crit_ids), high_high_ids = paste(hh_ids_th, collapse = ", "),
         class_none = sum(cl$class == "None"), class_low = sum(cl$class == "Low"),
         class_medium = sum(cl$class == "Medium"), class_high = sum(cl$class == "High"),
         kruskal_chi2 = round(unname(kw$statistic), 2), kruskal_p = kw$p.value)
})
print(as.data.frame(tbl$ftd_threshold_sensitivity))

# ---- 12. Leave-one-event-out validation of the classification ---------------
section("12. Leave-one-event-out validation of the recurrence classification")
# Traps are classified with n - 1 events and the classes are tested on the held-out event,
# so that classification and test never share observations.
loo_fold <- function(i) {
  held_out <- critical_dates[i]
  classes <- recurrence_table(obs_critical %>% filter(Data != held_out)) %>%
    mutate(class = classify_recurrence(rel_frequency))
  obs_critical %>% filter(Data == held_out) %>% select(Armadilha, Capturas) %>%
    inner_join(classes, by = "Armadilha") %>% mutate(held_out_event = held_out)
}
loo_folds <- map(seq_len(n_events), loo_fold)
tbl$loo_validation <- map_dfr(loo_folds, function(f) {
  kw <- kruskal.test(Capturas ~ class, data = f)
  rho <- suppressWarnings(cor.test(f$rel_frequency, f$Capturas, method = "spearman"))
  medians <- f %>% group_by(class) %>% summarise(m = median(Capturas), .groups = "drop") %>%
    pivot_wider(names_from = class, values_from = m, names_prefix = "median_capture_")
  bind_cols(tibble(held_out_event = f$held_out_event[1], n_traps_tested = nrow(f),
                   kruskal_chi2 = round(unname(kw$statistic), 2), kruskal_p = kw$p.value,
                   spearman_rho = round(unname(rho$estimate), 3), spearman_p = rho$p.value,
                   mean_capture_high = round(mean(f$Capturas[f$class == "High"]), 2),
                   mean_capture_none = round(mean(f$Capturas[f$class == "None"]), 2)), medians)
})
print(as.data.frame(tbl$loo_validation))
cat("Folds with Kruskal-Wallis p < 0.05:", sum(tbl$loo_validation$kruskal_p < 0.05), "of", n_events,
    "| Spearman rho: median", median(tbl$loo_validation$spearman_rho),
    "range", min(tbl$loo_validation$spearman_rho), "-", max(tbl$loo_validation$spearman_rho), "\n")
tbl$loo_pooled_by_class <- bind_rows(loo_folds) %>% group_by(class_from_remaining_events = class) %>%
  summarise(n_trap_events = n(), mean_capture = round(mean(Capturas), 2), median_capture = median(Capturas),
            prop_positive = round(mean(Capturas > 0), 3), .groups = "drop")
print(as.data.frame(tbl$loo_pooled_by_class))

# ---- 13. Insecticide applications: descriptive assessment -------------------
section("13. Insecticide applications - descriptive assessment")
readings <- date_level %>% filter(has_reading) %>%
  mutate(days_since_application = days_since(Data, application_dates),
         window = case_when(Data %in% application_dates ~ "application day",
                            !is.na(days_since_application) & days_since_application <= 7 ~ "1-7 days after",
                            !is.na(days_since_application) & days_since_application <= 14 ~ "8-14 days after",
                            !is.na(days_since_application) & days_since_application <= 30 ~ "15-30 days after",
                            TRUE ~ "> 30 days or none"),
         critical = Data %in% critical_dates)
window_levels <- c("application day", "1-7 days after", "8-14 days after", "15-30 days after", "> 30 days or none")
tbl$ftd_by_application_window <- readings %>% group_by(window) %>%
  summarise(n_dates = n(), mean_FTD = round(mean(orchard_FTD), 4), median_FTD = round(median(orchard_FTD), 4),
            n_critical_events = sum(critical), .groups = "drop") %>% arrange(factor(window, levels = window_levels))
print(as.data.frame(tbl$ftd_by_application_window))
tbl$ftd_before_after_application <- map_dfr(application_dates, function(a) {
  before <- readings %>% filter(Data <= a) %>% slice_max(Data, n = 1)
  after  <- readings %>% filter(Data > a + 3) %>% slice_min(Data, n = 1)
  tibble(application = a, date_before = before$Data[1], FTD_before = round(before$orchard_FTD[1], 3),
         date_after = after$Data[1], FTD_after = round(after$orchard_FTD[1], 3))
}) %>% mutate(change = FTD_after - FTD_before)
print(as.data.frame(tbl$ftd_before_after_application))
cat("Median FTD change after application:", median(tbl$ftd_before_after_application$change, na.rm = TRUE),
    "| decreases:", sum(tbl$ftd_before_after_application$change < 0, na.rm = TRUE), "of",
    sum(!is.na(tbl$ftd_before_after_application$change)), "\n")

# ---- 14. Negative-binomial GLMMs --------------------------------------------
if (config$run_glmm) {
  section("14. Negative-binomial GLMMs (glmmTMB)")
  glmm_data <- obs %>%
    mutate(season = factor(season_label(Data)), trap = factor(Armadilha), date = factor(Data),
           days_since_application = days_since(Data, application_dates),
           insecticide_7d = as.integer(Data %in% application_dates | (!is.na(days_since_application) & days_since_application <= 7)),
           critical = as.integer(Data %in% critical_dates), exposure = pmax(interval_days, 1)) %>%   # first inspection of the series: 1 day
    left_join(recurrence %>% select(Armadilha, recurrence_class), by = "Armadilha") %>%
    filter(!is.na(tmean_C), !is.na(precipitation_mm), !is.na(humidity_pct), !is.na(wind_ms))
  cat("Observations in the GLMMs:", nrow(glmm_data), "\n")

  # M1: climate, recent insecticide application and season; random intercepts for trap and date
  t0 <- Sys.time()
  m1 <- glmmTMB(Capturas ~ scale(tmean_C) + scale(precipitation_mm) + scale(humidity_pct) + scale(wind_ms) +
                  insecticide_7d + season + offset(log(exposure)) + (1 | trap) + (1 | date),
                family = nbinom2, data = glmm_data)
  cat("M1 fitted in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n"); print(summary(m1))

  # M2: does the recurrence class (defined on critical events) predict captures on NON-critical dates?
  # M2a is the same model without the class, on the same observations, so that the change in the
  # between-trap variance can be attributed to the class.
  glmm_noncritical <- glmm_data %>% filter(critical == 0, !is.na(recurrence_class))
  m2a <- glmmTMB(Capturas ~ scale(tmean_C) + scale(precipitation_mm) + scale(humidity_pct) + scale(wind_ms) +
                   insecticide_7d + season + offset(log(exposure)) + (1 | trap) + (1 | date),
                 family = nbinom2, data = glmm_noncritical)
  cat("M2a: same observations without the recurrence class\n"); print(summary(m2a))
  m2 <- glmmTMB(Capturas ~ recurrence_class + scale(tmean_C) + scale(precipitation_mm) + scale(humidity_pct) + scale(wind_ms) +
                  insecticide_7d + season + offset(log(exposure)) + (1 | trap) + (1 | date),
                family = nbinom2, data = glmm_noncritical)
  cat("M2: recurrence class on non-critical dates\n"); print(summary(m2))

  coef_table <- function(m, label) {
    s <- summary(m)$coefficients$cond
    tibble(model = label, term = rownames(s), estimate = round(s[, 1], 4), SE = round(s[, 2], 4),
           z = round(s[, 3], 2), p_value = s[, 4], rate_ratio = exp(s[, 1]),
           rr_lower95 = exp(s[, 1] - 1.96 * s[, 2]), rr_upper95 = exp(s[, 1] + 1.96 * s[, 2]))
  }
  # convergence and residual diagnostics: temporal autocorrelation of date-level mean Pearson residuals
  # (lag 1 in the sequence of inspection dates) and spatial autocorrelation of trap-level mean residuals
  diag_model <- function(m, d, label) {
    r <- residuals(m, type = "pearson")
    by_date <- tapply(r, as.character(d$Data), mean); by_date <- by_date[order(as.Date(names(by_date)))]
    acf1 <- cor(by_date[-1], by_date[-length(by_date)])
    by_trap <- tapply(r, as.character(d$Armadilha), mean)
    tr <- trap_summary %>% filter(Armadilha %in% names(by_trap)); w <- knn_weights(utm_coords(tr), config$k_neighbours)
    mi <- moran.test(as.numeric(by_trap[tr$Armadilha]), w)
    tibble(model = label, converged = is.null(m$fit$convergence) || m$fit$convergence == 0, n_obs = nrow(d),
           date_residual_lag1_autocorrelation = round(acf1, 3), trap_residual_moran_I = round(unname(mi$estimate[1]), 4),
           trap_residual_moran_p = mi$p.value)
  }
  tbl$glmm_coefficients <- bind_rows(coef_table(m1, "M1: climate + insecticide + season"),
                                     coef_table(m2a, "M2a: non-critical dates, without recurrence class"),
                                     coef_table(m2, "M2: non-critical dates, with recurrence class"))
  variance_table <- function(m, label) { v <- VarCorr(m)$cond
    tibble(model = label, group = names(v), variance = round(sapply(v, function(x) x[1, 1]), 4)) }
  tbl$glmm_random_effects <- bind_rows(variance_table(m1, "M1"), variance_table(m2a, "M2a"), variance_table(m2, "M2"))
  tbl$glmm_fit <- tibble(model = c("M1", "M2a", "M2"), n_obs = c(nrow(glmm_data), nrow(glmm_noncritical), nrow(glmm_noncritical)),
                         AIC = c(AIC(m1), AIC(m2a), AIC(m2)), dispersion_theta = c(sigma(m1), sigma(m2a), sigma(m2)))
  tbl$glmm_diagnostics <- bind_rows(diag_model(m1, glmm_data, "M1"), diag_model(m2a, glmm_noncritical, "M2a"), diag_model(m2, glmm_noncritical, "M2"))

  # Sensitivity analyses: (a) without the first inspection of each season, whose exposure interval spans the
  # inter-season gap; (b) with log(1 + mean capture of the previous inspection) as a covariate, to check whether
  # the associations change when the capture level of the preceding date is accounted for.
  first_dates <- obs %>% group_by(s = season_label(Data)) %>% summarise(d = min(Data), .groups = "drop") %>% pull(d)
  lag_tab <- obs %>% group_by(Data) %>% summarise(mc = mean(Capturas), .groups = "drop") %>% arrange(Data) %>% mutate(lag_mc = dplyr::lag(mc))
  f1 <- Capturas ~ scale(tmean_C) + scale(precipitation_mm) + scale(humidity_pct) + scale(wind_ms) + insecticide_7d + season + offset(log(exposure)) + (1 | trap) + (1 | date)
  f2 <- update(f1, . ~ . + recurrence_class)
  d1_nf <- glmm_data %>% filter(!Data %in% first_dates); d2_nf <- glmm_noncritical %>% filter(!Data %in% first_dates)
  d1_lag <- glmm_data %>% left_join(lag_tab %>% select(Data, lag_mc), by = "Data") %>% filter(!is.na(lag_mc)) %>% mutate(log_lag = log1p(lag_mc))
  d2_lag <- glmm_noncritical %>% left_join(lag_tab %>% select(Data, lag_mc), by = "Data") %>% filter(!is.na(lag_mc)) %>% mutate(log_lag = log1p(lag_mc))
  m1_nf <- glmmTMB(f1, family = nbinom2, data = d1_nf); m2_nf <- glmmTMB(f2, family = nbinom2, data = d2_nf)
  m1_lag <- glmmTMB(update(f1, . ~ . + log_lag), family = nbinom2, data = d1_lag); m2_lag <- glmmTMB(update(f2, . ~ . + log_lag), family = nbinom2, data = d2_lag)
  sens_row <- function(m, d, label) { s <- summary(m)$coefficients$cond; rr <- function(t) if (t %in% rownames(s)) sprintf("%.2f (%.2f–%.2f)", exp(s[t, 1]), exp(s[t, 1] - 1.96 * s[t, 2]), exp(s[t, 1] + 1.96 * s[t, 2])) else ""
    v <- VarCorr(m)$cond; r <- residuals(m, type = "pearson"); bd <- tapply(r, as.character(d$Data), mean); bd <- bd[order(as.Date(names(bd)))]
    tibble(model = label, n_obs = nrow(d), temperature = rr("scale(tmean_C)"), wind = rr("scale(wind_ms)"), humidity = rr("scale(humidity_pct)"),
           insecticide_7d = rr("insecticide_7d"), lag_mean_capture = rr("log_lag"), class_low = rr("recurrence_classLow"), class_medium = rr("recurrence_classMedium"),
           class_high = rr("recurrence_classHigh"), trap_variance = round(v$trap[1, 1], 3), date_residual_lag1 = round(cor(bd[-1], bd[-length(bd)]), 3)) }
  tbl$glmm_sensitivity <- bind_rows(sens_row(m1, glmm_data, "M1 (reference)"), sens_row(m1_nf, d1_nf, "M1 without the first inspection of each season"),
                                    sens_row(m1_lag, d1_lag, "M1 + log(1 + mean capture of the previous inspection)"),
                                    sens_row(m2, glmm_noncritical, "M2 (reference)"), sens_row(m2_nf, d2_nf, "M2 without the first inspection of each season"),
                                    sens_row(m2_lag, d2_lag, "M2 + log(1 + mean capture of the previous inspection)"))
  print(as.data.frame(tbl$glmm_sensitivity))
  print(as.data.frame(tbl$glmm_random_effects)); print(as.data.frame(tbl$glmm_fit)); print(as.data.frame(tbl$glmm_diagnostics))
}

# ---- 15. Figures ------------------------------------------------------------
section("15. Figures")
# Figures 1 to 5 and S1 to S3 of the manuscript (600 dpi). Figure 1 is the study-area map produced in QGIS
# (data/study_area_map.png); it is copied, not generated.
map_file <- file.path(dirname(config$master_file), "study_area_map.png")   # data/ locally, /content on Google Colab
if (file.exists(map_file)) {
  stopifnot(file.copy(map_file, file.path(fig_dir, "Figure 1.png"), overwrite = TRUE))
} else message("Figure 1 (study-area map) not found at ", map_file, "; it is produced in QGIS and is not generated by this script")
theme_paper <- theme_minimal(base_size = 14) + theme(legend.position = "right")
axis_labels <- labs(x = "Longitude (\u00b0W)", y = "Latitude (\u00b0S)")
lisa_levels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
lisa_colours <- c("High-High" = "#D55E00", "Low-Low" = "#0072B2", "High-Low" = "#E69F00", "Low-High" = "#F8766D", "Not significant" = "#00BFC4")
class_colours <- c("Not observed" = "grey40", "None" = "grey65", "Low" = "#56B4E9", "Medium" = "#E69F00", "High" = "#D55E00")
save_fig <- function(p, name, w = 8, h = 6) ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 600, bg = "white")

save_fig(ggplot(trap_summary, aes(Longitude, Latitude, size = mean_capture)) + geom_point(alpha = 0.7, colour = "#00BFC4") +
           theme_paper + axis_labels + labs(size = "Mean capture"), "Figure 2.png")
fig2 <- lisa_full %>% mutate(cl = factor(cluster_full, levels = lisa_levels))
save_fig(ggplot(fig2, aes(Longitude, Latitude, colour = cl)) + geom_point(size = 3) +
           scale_colour_manual(values = lisa_colours, drop = TRUE) + theme_paper + axis_labels + labs(colour = "Cluster Type"), "Figure 3.png")
fig3 <- lisa_crit %>% mutate(cl = factor(cluster, levels = lisa_levels), fdr = cluster_fdr == "High-High")
save_fig(ggplot(fig3, aes(Longitude, Latitude, colour = cl)) + geom_point(size = 3) +
           geom_point(data = filter(fig3, fdr), shape = 21, size = 5.5, stroke = 1.1, colour = "black", fill = NA) +
           ggrepel::geom_text_repel(data = filter(fig3, cluster == "High-High"), aes(label = sub("N\u00b0", "", Armadilha)), colour = "black", size = 3,
                                    min.segment.length = 0, segment.size = 0.3, box.padding = 0.45, point.padding = 0.25, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
           scale_colour_manual(values = lisa_colours, drop = TRUE) + theme_paper + axis_labels + labs(colour = "Cluster Type"), "Figure 4.png")
fig4 <- all_traps %>% left_join(recurrence %>% select(Armadilha, recurrence_class), by = "Armadilha") %>%
  mutate(cl = factor(ifelse(is.na(recurrence_class), "Not observed", as.character(recurrence_class)), levels = names(class_colours)))
save_fig(ggplot(fig4, aes(Longitude, Latitude, colour = cl, shape = cl)) + geom_point(size = 3) +
           scale_colour_manual(values = class_colours) + scale_shape_manual(values = c("Not observed" = 4, "None" = 16, "Low" = 16, "Medium" = 16, "High" = 16)) +
           theme_paper + axis_labels + labs(colour = "Recurrence class", shape = "Recurrence class"), "Figure 5.png")
figS1 <- lisa_crit %>% left_join(hh_loo %>% filter(cluster == "High-High") %>% count(Armadilha = trap, name = "n_loo"), by = "Armadilha") %>%
  mutate(n_loo = replace_na(n_loo, 0L))
save_fig(ggplot(figS1, aes(Longitude, Latitude, colour = n_loo)) + geom_point(size = 3) +
           scale_colour_gradient(low = "grey85", high = "#D55E00", breaks = 0:n_events) + theme_paper + axis_labels +
           labs(colour = paste0("High-High\n(of ", n_events, " runs)")), "Figure S1.png")
figS2 <- tbl$loo_validation %>% mutate(event = factor(format(held_out_event, "%d %b %Y"), levels = format(sort(held_out_event), "%d %b %Y")))
save_fig(ggplot(figS2, aes(event, spearman_rho)) + geom_col(fill = "#0072B2", width = 0.6) + geom_hline(yintercept = 0) + theme_paper +
           labs(x = "Held-out critical event", y = "Spearman's rho") + theme(axis.text.x = element_text(angle = 30, hjust = 1)), "Figure S2.png")
figS3 <- tbl$lisa_by_season %>% mutate(cl = factor(cluster, levels = lisa_levels))
save_fig(ggplot(figS3, aes(Longitude, Latitude, colour = cl)) + geom_point(size = 1.8) +
           geom_point(data = filter(figS3, in_critical_cluster), shape = 21, size = 3.2, stroke = 0.8, colour = "black", fill = NA) +
           scale_colour_manual(values = lisa_colours, drop = TRUE) + facet_wrap(~ season, ncol = 3) + theme_minimal(base_size = 12) +
           theme(legend.position = "bottom", axis.text = element_text(size = 7)) + axis_labels + labs(colour = "Cluster type"), "Figure S3.png")

# ---- 16. Export -------------------------------------------------------------
section("16. Export")
tbl$main_text_numbers <- tibble(
  item = c("valid observations", "traps", "monitoring dates", "critical events",
           "traps with reading in critical events", "Moran I full period (k=4)", "p full period (analytical)",
           "p full period (permutation)", "Moran I critical events (k=4)", "p critical (analytical)", "p critical (permutation)",
           "High-High traps critical (analytical p<0.05)", "High-High traps critical (FDR q<0.05)",
           "High-High traps full period (analytical p<0.05)", "High-High traps full period (FDR q<0.05)", "High-High traps in both analyses",
           "Kruskal-Wallis chi2 baseline", "Kruskal-Wallis p baseline",
           "mean 4-NN distance (m)", "LOO folds with KW p<0.05", "median Spearman rho (LOO)",
           "traps classified by recurrence", "traps not observed on any critical date"),
  value = c(nrow(obs), n_distinct(obs$Armadilha), n_distinct(obs$Data), n_events, nrow(intensity_critical),
            round(moran_full$estimate[1], 4), signif(moran_full$p.value, 3), mc_full$p.value,
            round(moran_crit$estimate[1], 4), signif(moran_crit$p.value, 3), mc_crit$p.value,
            length(hh_crit), count_hh(lisa_crit$cluster_fdr), length(hh_full), count_hh(lisa_full$cluster_fdr),
            length(intersect(hh_full, hh_crit)), round(unname(kw_baseline$statistic), 2), signif(kw_baseline$p.value, 3),
            round(mean(nn_distance), 1), sum(tbl$loo_validation$kruskal_p < 0.05), median(tbl$loo_validation$spearman_rho),
            nrow(recurrence), length(unobserved_traps)))
tbl$per_trap <- lisa_full %>% select(Armadilha, Latitude, Longitude, mean_capture_full = mean_capture, cluster_full, q_fdr_full = q_fdr) %>%
  left_join(trap_data %>% select(-Latitude, -Longitude), by = "Armadilha") %>%
  left_join(lisa_crit %>% select(Armadilha, Ii_critical = Ii, p_critical = p_analytical, q_fdr_critical = q_fdr, cluster_critical = cluster), by = "Armadilha") %>%
  rename(trap = Armadilha)

tbl <- c(tbl["main_text_numbers"], tbl[names(tbl) != "main_text_numbers"])
tbl <- lapply(tbl, function(d) { d <- as.data.frame(d); d[] <- lapply(d, function(x) if (inherits(x, "Date")) format(x) else x); d })
write_xlsx(tbl, file.path(config$output_dir, "tables.xlsx"))
cat("Tables written to", file.path(config$output_dir, "tables.xlsx"), "\n")
cat("\nNumbers for the main text:\n"); print(tbl$main_text_numbers)
sink()
writeLines(c(R.version.string, paste(Sys.info()[c("sysname", "release")], collapse = " "), "", capture.output(sessionInfo())),
           file.path(config$output_dir, "session_info.txt"))
cat("\nDone. Log:", file.path(config$output_dir, "analysis_log.txt"), "\n")