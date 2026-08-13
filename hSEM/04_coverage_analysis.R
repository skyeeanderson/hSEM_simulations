# ============================================================
#  Interval coverage analysis for the three species-interaction
#  values (SIVs): b_pred, c_pred, b_meso.
#
#  Coverage = fraction of replicates whose 95% interval contains
#  the true value. Two-step models report Wald confidence
#  intervals; the hSEMs report Bayesian credible intervals. A
#  well-calibrated interval covers near 0.95; a model that fails
#  to propagate first-stage uncertainty covers below 0.95 because
#  its intervals are too narrow.
#
#  Reads the compiled comparison table produced by
#  02_extract_results.R, which now carries the *_lo / *_hi bound
#  rows for each SIV alongside the point estimates.
#
#  Run as: Rscript 04_coverage_analysis.R
# ============================================================

rm(list = ls())

library(dplyr)
library(tidyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
results_root <- "/scratch/user/uqsand24/simulations/results/"
comps_file   <- file.path(results_root, "v4rn_all_scenarios_comparisons.rds")
out_file     <- file.path(results_root, "v4rn_coverage_table.rds")
csv_file     <- file.path(results_root, "v4rn_coverage_table.csv")

comps <- readRDS(comps_file)
cat("Loaded comparisons:", nrow(comps), "rows,", ncol(comps), "cols\n")

# ── Model columns (the eleven reported models) ────────────────────────────────
# coabund is excluded (disabled / NA). These names are the comparison-table
# column names, not the paper model numbers.
model_cols <- c(
  "psem_maxcount",   # 1  naive count-index SEM
  "psem_blup_null",  # 2  two-step null N-mixture SEM
  "psem_blup_full",  # 3  two-step N-mixture SEM with covariates
  "sem_nmix",        # 4  abundance hSEM
  "rn_blup_null",    # 5  two-step null RN SEM
  "rn_blup_full",    # 6  two-step RN SEM with covariates
  "rn_nmix",         # 7  RN hSEM
  "naive_sem",       # 8  naive detection-rate SEM
  "null_occ_sem",    # 9  two-step null occupancy SEM
  "full_occ_sem",    # 10 two-step occupancy SEM with covariates
  "occ_sem"          # 11 occupancy hSEM
)

model_labels <- c(
  psem_maxcount  = "1. Naive count-index SEM",
  psem_blup_null = "2. Two-step null N-mixture SEM",
  psem_blup_full = "3. Two-step N-mixture SEM (covariates)",
  sem_nmix       = "4. Abundance hSEM",
  rn_blup_null   = "5. Two-step null RN SEM",
  rn_blup_full   = "6. Two-step RN SEM (covariates)",
  rn_nmix        = "7. RN hSEM",
  naive_sem      = "8. Naive detection-rate SEM",
  null_occ_sem   = "9. Two-step null occupancy SEM",
  full_occ_sem   = "10. Two-step occupancy SEM (covariates)",
  occ_sem        = "11. Occupancy hSEM"
)

# Which models use credible (Bayesian) vs confidence (Wald) intervals — for the
# footnote in the supplementary table. hSEMs are credible; all others Wald.
interval_type <- ifelse(model_cols %in% c("sem_nmix", "rn_nmix", "occ_sem"),
                        "credible", "confidence")
names(interval_type) <- model_cols

# ── SIVs and their bound-row names ────────────────────────────────────────────
sivs <- c("b_pred", "c_pred", "b_meso")

# ── Regime from scenario name ─────────────────────────────────────────────────
comps <- comps %>%
  mutate(regime = ifelse(grepl("abundant", scenario), "Abundant", "Rare"))

# ── Reshape: one row per replicate, with est/lo/hi side by side ───────────────
# For each SIV we need, per replicate (identified by scenario + seed), the
# point estimate row, the _lo row and the _hi row, for every model column.
# We pull them out and join on the replicate id.

rep_id_cols <- c("scenario", "seed", "regime")

coverage_one_siv <- function(siv) {
  
  lo_name <- paste0(siv, "_lo")
  hi_name <- paste0(siv, "_hi")
  
  # Point estimates: rows where estimand == siv
  est <- comps %>%
    filter(estimand == siv) %>%
    select(all_of(c(rep_id_cols, "truth", model_cols))) %>%
    pivot_longer(all_of(model_cols), names_to = "model", values_to = "est")
  
  lo <- comps %>%
    filter(estimand == lo_name) %>%
    select(all_of(c(rep_id_cols, model_cols))) %>%
    pivot_longer(all_of(model_cols), names_to = "model", values_to = "lo")
  
  hi <- comps %>%
    filter(estimand == hi_name) %>%
    select(all_of(c(rep_id_cols, model_cols))) %>%
    pivot_longer(all_of(model_cols), names_to = "model", values_to = "hi")
  
  # Join est, lo, hi on replicate id + model
  joined <- est %>%
    left_join(lo, by = c(rep_id_cols, "model")) %>%
    left_join(hi, by = c(rep_id_cols, "model")) %>%
    mutate(siv = siv)
  
  joined
}

all_siv <- bind_rows(lapply(sivs, coverage_one_siv))

# ── Compute coverage and width per model x SIV x regime ───────────────────────
# covered: does the 95% interval contain the truth?
# A replicate contributes only if all of truth, lo, hi are non-missing.
coverage_table <- all_siv %>%
  mutate(covered = ifelse(!is.na(lo) & !is.na(hi) & !is.na(truth),
                          truth >= lo & truth <= hi, NA),
         width   = ifelse(!is.na(lo) & !is.na(hi), hi - lo, NA)) %>%
  group_by(siv, regime, model) %>%
  summarise(
    n_reps    = sum(!is.na(covered)),
    coverage  = mean(covered, na.rm = TRUE),
    mean_width = mean(width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model_label   = model_labels[model],
    interval_type = interval_type[model]
  ) %>%
  # Order rows by model number (via the factor level of model_cols)
  mutate(model = factor(model, levels = model_cols)) %>%
  arrange(siv, regime, model)

# ── Report ────────────────────────────────────────────────────────────────────
cat("\n===== Coverage table (fraction of 95% intervals containing truth) =====\n")
for (s in sivs) {
  cat(sprintf("\n--- SIV: %s ---\n", s))
  sub <- coverage_table %>% filter(siv == s)
  # Wide view: coverage by regime
  wide_cov <- sub %>%
    select(model_label, regime, coverage) %>%
    pivot_wider(names_from = regime, values_from = coverage)
  print(as.data.frame(wide_cov), digits = 2, row.names = FALSE)
}

cat("\n===== Mean interval width by SIV, regime, model =====\n")
for (s in sivs) {
  cat(sprintf("\n--- SIV: %s ---\n", s))
  sub <- coverage_table %>% filter(siv == s)
  wide_w <- sub %>%
    select(model_label, regime, mean_width) %>%
    pivot_wider(names_from = regime, values_from = mean_width)
  print(as.data.frame(wide_w), digits = 3, row.names = FALSE)
}

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(coverage_table, out_file)
write.csv(coverage_table, csv_file, row.names = FALSE)
cat(sprintf("\nSaved coverage table: %s\n", out_file))
cat(sprintf("Saved CSV: %s\n", csv_file))

# ── Footnote reminder ─────────────────────────────────────────────────────────
cat("\nNOTE for the supplement: hSEM intervals (models 4, 7, 11) are 95% Bayesian\n")
cat("credible intervals; all other models report 95% Wald confidence intervals.\n")
cat("Both are the model's stated 95% uncertainty region; coverage compares how\n")
cat("often each contains the truth. Nominal coverage is 0.95.\n")