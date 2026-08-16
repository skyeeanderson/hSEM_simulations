# ============================================================
#  Build Table S3: interval coverage and mean width for the
#  three SIVs, by model, regime, and scenario.
#
#  Reads the compiled comparison table, EXCLUDES the three
#  occupancy two-step convergence-failure replicates
#  (nSites=50, abundant, seeds 69/79/87), recomputes coverage
#  and width cleanly, and formats a flextable for the supplement.
#
#  Run as: Rscript build_table_S3.R
# ============================================================

library(dplyr)
library(tidyr)
library(flextable)
library(officer)

# ── Paths ─────────────────────────────────────────────────────────────────────
results_root <- "~/Dropbox/Skye PhD UQ projects/Skye Thesis/Ch3 Skye species interactions/Analyses- Skye species interaction sims/skye_ch3_analysis/simulations/results/integrated_SEM_v4/"
comps_file   <- file.path(results_root, "v4rn_all_scenarios_comparisons.rds")
out_docx     <- file.path(results_root, "Table_S3_coverage.docx")

comps <- readRDS(comps_file)

# ── Exclude convergence-failure replicates ────────────────────────────────────
# Three replicates (nSites = 50, abundant scenario, seeds 69/79/87) showed
# occupancy first-stage convergence failure (Wald SEs exploded, e.g. interval
# upper bounds > 600). Excluded from all summaries, consistent with the bias
# analysis.
excl_seeds    <- c(69, 79, 87)
excl_scenario <- "abundant_50_meso"
n_before <- nrow(comps)
comps <- comps[!(comps$scenario == excl_scenario & comps$seed %in% excl_seeds), ]
cat(sprintf("Excluded %d rows (%d convergence-failure replicates).\n",
            n_before - nrow(comps), length(excl_seeds)))

# ── Model columns and display labels (eleven reported models) ─────────────────
model_cols <- c(
  "psem_maxcount", "psem_blup_null", "psem_blup_full", "sem_nmix",
  "rn_blup_null", "rn_blup_full", "rn_nmix",
  "naive_sem", "null_occ_sem", "full_occ_sem", "occ_sem"
)
model_labels <- c(
  psem_maxcount  = "1. Naive count-index SEM",
  psem_blup_null = "2. Two-step null N-mixture SEM",
  psem_blup_full = "3. Two-step N-mixture SEM (cov.)",
  sem_nmix       = "4. Abundance hSEM",
  rn_blup_null   = "5. Two-step null RN SEM",
  rn_blup_full   = "6. Two-step RN SEM (cov.)",
  rn_nmix        = "7. RN hSEM",
  naive_sem      = "8. Naive detection-rate SEM",
  null_occ_sem   = "9. Two-step null occupancy SEM",
  full_occ_sem   = "10. Two-step occupancy SEM (cov.)",
  occ_sem        = "11. Occupancy hSEM"
)

sivs      <- c("b_pred", "c_pred", "b_meso")
siv_names <- c(b_pred = "Predator \u2192 prey",
               c_pred = "Predator \u2192 meso",
               b_meso = "Meso \u2192 prey")

comps <- comps %>% mutate(regime = ifelse(grepl("abundant", scenario),
                                          "Abundant", "Rare"))

# ── Compute coverage and mean width per model x SIV x regime ──────────────────
compute_siv <- function(siv) {
  lo_n <- paste0(siv, "_lo"); hi_n <- paste0(siv, "_hi")
  
  est <- comps %>% filter(estimand == siv) %>%
    select(scenario, seed, regime, truth, all_of(model_cols)) %>%
    pivot_longer(all_of(model_cols), names_to = "model", values_to = "est")
  lo <- comps %>% filter(estimand == lo_n) %>%
    select(scenario, seed, regime, all_of(model_cols)) %>%
    pivot_longer(all_of(model_cols), names_to = "model", values_to = "lo")
  hi <- comps %>% filter(estimand == hi_n) %>%
    select(scenario, seed, regime, all_of(model_cols)) %>%
    pivot_longer(all_of(model_cols), names_to = "model", values_to = "hi")
  
  est %>%
    left_join(lo, by = c("scenario","seed","regime","model")) %>%
    left_join(hi, by = c("scenario","seed","regime","model")) %>%
    mutate(siv = siv,
           covered = ifelse(!is.na(lo) & !is.na(hi) & !is.na(truth),
                            truth >= lo & truth <= hi, NA),
           width   = ifelse(!is.na(lo) & !is.na(hi), hi - lo, NA))
}

all_siv <- bind_rows(lapply(sivs, compute_siv))

summ <- all_siv %>%
  group_by(siv, model, regime) %>%
  summarise(coverage   = mean(covered, na.rm = TRUE),
            mean_width = mean(width,   na.rm = TRUE),
            .groups = "drop")

# ── Reshape to wide: one row per model, columns = SIV x regime x metric ───────
wide <- summ %>%
  mutate(cell = sprintf("%.2f (%.2f)", coverage, mean_width)) %>%   # "coverage (width)"
  select(siv, model, regime, cell) %>%
  pivot_wider(names_from = c(siv, regime), values_from = cell,
              names_sep = "_") %>%
  mutate(model = factor(model, levels = model_cols)) %>%
  arrange(model) %>%
  mutate(Model = model_labels[as.character(model)]) %>%
  select(Model,
         `b_pred_Abundant`, `b_pred_Rare`,
         `c_pred_Abundant`, `c_pred_Rare`,
         `b_meso_Abundant`, `b_meso_Rare`)

# ── Build flextable ───────────────────────────────────────────────────────────
ft <- flextable(wide)

# Two-row header: SIV group over Abundant/Rare
ft <- set_header_labels(ft,
                        Model = "Model",
                        b_pred_Abundant = "Abundant", b_pred_Rare = "Rare",
                        c_pred_Abundant = "Abundant", c_pred_Rare = "Rare",
                        b_meso_Abundant = "Abundant", b_meso_Rare = "Rare")

ft <- add_header_row(ft,
                     values = c("", siv_names["b_pred"], siv_names["c_pred"], siv_names["b_meso"]),
                     colwidths = c(1, 2, 2, 2))

ft <- theme_booktabs(ft)
ft <- align(ft, part = "all", align = "center")
ft <- align(ft, j = 1, part = "all", align = "left")
ft <- bold(ft, part = "header")
ft <- fontsize(ft, size = 8, part = "all")
ft <- autofit(ft)

# Highlight the well-calibrated hSEM rows (optional visual aid; comment out if
# the journal prefers no shading)
hsem_rows <- which(wide$Model %in% c("4. Abundance hSEM", "7. RN hSEM",
                                     "11. Occupancy hSEM"))
ft <- bg(ft, i = hsem_rows, bg = "#EFEFEF", part = "body")

# ── Save as Word doc (landscape) ──────────────────────────────────────────────
sect <- prop_section(page_size = page_size(orient = "landscape"))
doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_end_block_section(block_section(sect))
print(doc, target = out_docx)

cat("Saved:", out_docx, "\n")
cat("\nEach cell shows: coverage (mean interval width).\n")
cat("Nominal coverage = 0.95.\n")
