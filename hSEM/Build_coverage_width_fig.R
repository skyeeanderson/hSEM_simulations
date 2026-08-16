# ============================================================
#  Coverage vs interval width figure (main text).
#
#  One point per model x SIV x regime. Coverage on the y-axis,
#  mean interval width on the x-axis. The well-calibrated hSEMs
#  sit in the top-left "good" region (high coverage, narrow
#  width); two-step models fall low (under-covering); the
#  occupancy hSEM sits low and far right (wide but still
#  missing). A horizontal line marks nominal 0.95 coverage.
#
#  Excludes the three occupancy-convergence-failure replicates,
#  consistent with all other summaries.
#
#  Run as: Rscript coverage_width_figure.R
# ============================================================

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# ── Paths (EDIT for your machine) ─────────────────────────────────────────────
comps_file <- "~/Dropbox/Skye PhD UQ projects/Skye Thesis/Ch3 Skye species interactions/Analyses- Skye species interaction sims/skye_ch3_analysis/simulations/results/integrated_SEM_v4/v4rn_all_scenarios_comparisons.rds"
out_png    <- "~/Dropbox/Skye PhD UQ projects/Skye Thesis/Ch3 Skye species interactions/Analyses- Skye species interaction sims/skye_ch3_analysis/simulations/figures/integrated_SEM_v4/Fig_coverage_width.png"

comps <- readRDS(comps_file)

# ── Exclude convergence-failure replicates ────────────────────────────────────
excl_seeds    <- c(69, 79, 87)
excl_scenario <- "abundant_50_meso"
comps <- comps[!(comps$scenario == excl_scenario & comps$seed %in% excl_seeds), ]

# ── Models, labels, families ──────────────────────────────────────────────────
model_cols <- c(
  "psem_maxcount", "psem_blup_null", "psem_blup_full", "sem_nmix",
  "rn_blup_null", "rn_blup_full", "rn_nmix",
  "naive_sem", "null_occ_sem", "full_occ_sem", "occ_sem"
)
model_num <- c(
  psem_maxcount = "1", psem_blup_null = "2", psem_blup_full = "3",
  sem_nmix = "4", rn_blup_null = "5", rn_blup_full = "6", rn_nmix = "7",
  naive_sem = "8", null_occ_sem = "9", full_occ_sem = "10", occ_sem = "11"
)
# Family, for colour
model_family <- c(
  psem_maxcount = "Abundance", psem_blup_null = "Abundance",
  psem_blup_full = "Abundance", sem_nmix = "Abundance",
  rn_blup_null = "Royle-Nichols", rn_blup_full = "Royle-Nichols",
  rn_nmix = "Royle-Nichols",
  naive_sem = "Occupancy", null_occ_sem = "Occupancy",
  full_occ_sem = "Occupancy", occ_sem = "Occupancy"
)
# Estimation type, for shape (hSEMs vs the rest)
model_type <- case_when(model_cols %in% c("sem_nmix","rn_nmix","occ_sem") ~ "Joint hSEM",
                        model_cols %in% c("psem_maxcount", "naive_sem") ~ "Naive",
                        model_cols %in% c("psem_blup_null", "psem_blup_full",
                                          "rn_blup_null", "rn_blup_full",
                                          "null_occ_sem", "full_occ_sem") ~ "Two-step")
names(model_type) <- model_cols

sivs      <- c("b_pred", "c_pred", "b_meso")
siv_names <- c(b_pred = "Predator \u2192 prey",
               c_pred = "Predator \u2192 meso",
               b_meso = "Meso \u2192 prey")

comps <- comps %>% mutate(regime = ifelse(grepl("abundant", scenario),
                                          "Abundant", "Rare"))

# ── Compute coverage and width ────────────────────────────────────────────────
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

pts <- bind_rows(lapply(sivs, compute_siv)) %>%
  group_by(siv, model, regime) %>%
  summarise(coverage = mean(covered, na.rm = TRUE),
            width    = mean(width,   na.rm = TRUE), .groups = "drop") %>%
  mutate(
    #num    = model_num[model],
    family = model_family[model],
    type   = model_type[model],
    siv_lab = factor(siv_names[siv], levels = siv_names)
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
# Width is log-scaled because Model 11's widths are ~10-30x the others; on a
# linear axis every calibrated model collapses into a sliver at the left.
fam_cols <- c("Abundance" = "#78A641FF", "Royle-Nichols" = "#12A2A8FF",
              "Occupancy" = "#FFAA0EFF")

p <- ggplot(pts, aes(x = width, y = coverage,
                     colour = family, shape = type)) +
  # nominal coverage reference
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.90, ymax = 1.0,
           fill = "grey90", alpha = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_point(size = 4, stroke = 1) +
  ggrepel::geom_text_repel(aes(label = num), size = 4, seed = 1,
                           show.legend = FALSE, max.overlaps = 20,
                           colour = "black") +
  facet_grid(regime ~ siv_lab) +
  scale_x_continuous(trans = "log10",
                     breaks = c(0.1, 0.3, 1, 3),
                     labels = c("0.1", "0.3", "1", "3")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = fam_cols, name = "State formulation") +
  scale_shape_manual(values = c("Joint hSEM" = 17,
                                "Naive" = 16,
                                "Two-step" = 15),
                     name = "Approaches") +
  labs(
    x = "Mean 95% interval width (log scale)",
    y = "Coverage (proportion of intervals containing the truth)"
    ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold"),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(out_png, p, width = 11, height = 7, dpi = 300)
cat("Saved:", out_png, "\n")
