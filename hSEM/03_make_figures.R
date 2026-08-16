# ============================================================
#  Compile v4 simulation results — FA/FB/FC percent bias figures
#
#  Layout: facet_grid(abundance ~ estimand_label)
#    - Rows    = Abundant / Rare
#    - Columns = estimands (e.g. 3 SIVs, 3 env effects, 2-3 intercepts)
#    - x-axis  = model (method_label)
#    - Dodged within each model by nSites (50/500/1000)
#    - Coloured jitter = model colour
#    - Black CI bars + black summary points
#    - Shapes: circle=50, triangle=500, square=1000
#
#  Produces 12 figures (6 abundance + 6 occupancy), LF only.
# ============================================================

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(paletteer)
library(scales)

# ── Paths ─────────────────────────────────────────────────────────────────────
results_dir <- "~/Dropbox/Skye PhD UQ projects/Skye Thesis/Ch3 Skye species interactions/Analyses- Skye species interaction sims/skye_ch3_analysis/simulations/results/integrated_SEM_v4/"
fig_dir     <- "~/Dropbox/Skye PhD UQ projects/Skye Thesis/Ch3 Skye species interactions/Analyses- Skye species interaction sims/skye_ch3_analysis/simulations/figures/integrated_SEM_v4/"
# dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
comparisons <- readRDS(file.path(results_dir, "v4rn_all_scenarios_comparisons.rds"))

cat("Comparisons:", nrow(comparisons), "rows,", ncol(comparisons), "cols\n")

# ── Add scenario metadata ─────────────────────────────────────────────────────
comparisons <- comparisons %>%
  mutate(
    abundance = ifelse(grepl("^abundant", scenario), "Abundant", "Rare"),
    nSites    = case_when(
      grepl("_50",   scenario) & !grepl("_500",  scenario) ~   50L,
      grepl("_500",  scenario) & !grepl("_1000", scenario) ~  500L,
      grepl("_1000", scenario)                              ~ 1000L
    ),
    meso = grepl("meso", scenario)
  )

# ── Model groupings ───────────────────────────────────────────────────────────
all_model_cols <- c(
  "psem_maxcount", "psem_blup_null", "psem_blup_full",
  "coabund", "sem_nmix",
  "rn_blup_null", "rn_blup_full", "rn_nmix",
  "naive_sem", "null_occ_sem", "full_occ_sem", "occ_sem"
)

# coabund excluded from main comparison
abund_lmer_models   <- c("psem_maxcount", "psem_blup_null", "psem_blup_full")
abund_nimble_models <- c("sem_nmix")
abund_all_models    <- c(abund_lmer_models, abund_nimble_models)

rn_lmer_models   <- c("rn_blup_null", "rn_blup_full")
rn_nimble_models <- c("rn_nmix")
rn_all_models    <- c(rn_lmer_models, rn_nimble_models)

occ_lmer_models   <- c("naive_sem", "null_occ_sem", "full_occ_sem")
occ_nimble_models <- c("occ_sem")
occ_all_models    <- c(occ_lmer_models, occ_nimble_models)

# ── Method labels ─────────────────────────────────────────────────────────────
method_labels <- c(
  psem_maxcount  = "1: Naive count index\n SEM",
  psem_blup_null = "2: Two-step null\nN-mix SEM",
  psem_blup_full = "3: Two-step N-mix\nwith covariates SEM",
  sem_nmix       = "4: Abundance hSEM",
  rn_blup_null   = "5: Two-step null\nRN SEM",
  rn_blup_full   = "6: Two-step RN\nwith covariates SEM",
  rn_nmix        = "7: Royle-Nichols hSEM",
  naive_sem      = "8: Naive detection\nrate SEM",
  null_occ_sem   = "9: Two-step null\nOcc SEM",
  full_occ_sem   = "10: Two-step Occ\nwith covariates SEM",
  occ_sem        = "11: Occupancy hSEM"
)

# ── Colours ───────────────────────────────────────────────────────────────────
make_pal <- function(pal, df) {
  matching <- names(pal)[names(pal) %in% unique(df$model)]
  setNames(pal[matching], method_labels[matching])
}

full_pal <- as.character(paletteer::paletteer_d("ggthemes::Classic_Cyclic", n = 13))
pal_12   <- setNames(full_pal[c(1, 2, 3, 5, 6, 7, 8, 9, 11, 1, 2, 3)], all_model_cols)
pal_12["psem_maxcount"]  <- "#12A2A8FF"
pal_12["psem_blup_null"] <- "#2CA030FF"
pal_12["psem_blup_full"] <- "#78A641FF"
pal_12["sem_nmix"]       <- "#FFAA0EFF"
# RN mirrors abundance palette with transparency shift via alpha in plots
pal_12["rn_blup_null"]   <- "#2CA030FF"
pal_12["rn_blup_full"]   <- "#78A641FF"
pal_12["rn_nmix"]        <- "#FFAA0EFF"
# Occupancy mirrors abundance palette
pal_12["naive_sem"]      <- "#12A2A8FF"
pal_12["null_occ_sem"]   <- "#2CA030FF"
pal_12["full_occ_sem"]   <- "#78A641FF"
pal_12["occ_sem"]        <- "#FFAA0EFF"

# ── Pivot to long format, LF only ────────────────────────────────────────────
comp_long <- comparisons %>%
  pivot_longer(
    cols      = all_of(all_model_cols),
    names_to  = "model",
    values_to = "estimate"
  ) %>%
  filter(!is.na(truth)) %>%
  mutate(
    percent_bias = (estimate - truth) / abs(truth) * 100,
    method_label = factor(method_labels[model], levels = method_labels),
    nSites_label = factor(paste0("n = ", nSites),
                          levels = c("n = 50", "n = 500", "n = 1000"))
  ) %>%
  filter(!use_landscape)   # LF only

cat("Rows in comp_long:", nrow(comp_long), "\n")

comp_long_clean <- comp_long %>%
  filter(!(model == "null_occ_sem" & 
             seed %in% c(69, 79, 87) & 
             nSites == 50 & 
             abundance == "Abundant"))

# ── Core figure function ──────────────────────────────────────────────────────
#
#  Layout: facet_grid(abundance ~ estimand_label)
#    Rows = Abundant / Rare
#    Cols = estimands
#    x    = model (method_label)
#    Dodge within model by nSites
#    Coloured jitter, BLACK summary stats, shapes for nSites
#
make_bias_fig_FABC <- function(df, title, estimand_map,
                               subtitle   = "Percent bias = (estimate \u2212 truth) / |truth| \u00d7 100",
                               ylim       = c(100, 100),
                               fig_width  = 12,
                               fig_height = 7,
                               pal        = pal_12,
                               fname      = NULL) {
  
  df <- df %>%
    filter(estimand %in% names(estimand_map),
    is.finite(percent_bias)) %>%
    mutate(
      estimand_label = factor(estimand_map[estimand], levels = estimand_map),
      method_label   = droplevels(method_label)
    )
  
  if (nrow(df) == 0) {
    warning("Empty dataframe, skipping: ", fname)
    return(invisible(NULL))
  }
  
  pal_use <- make_pal(pal, df)
  
  p <- ggplot(df, aes(x = method_label, y = percent_bias, colour = method_label)) +
    
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.6, colour = "grey40") +
    
    # Coloured jitter — dodged by nSites within each model
    geom_jitter(aes(group = nSites_label, shape = nSites_label),
                position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width  = 0.75),
                alpha = 0.3, size = 1.2) +
    
    # Black CI bars — dodged by nSites
    stat_summary(aes(group = nSites_label, shape = nSites_label),
                 fun      = mean,
                 fun.min  = function(x) mean(x) - sd(x),
                 fun.max  = function(x) mean(x) + sd(x),
                 geom     = "linerange",
                 colour   = "black",
                 linewidth = 0.3,
                 position = position_dodge(0.75)) +

    # Black summary points
    stat_summary(aes(group = nSites_label, shape = nSites_label),
                 fun      = mean,
                 geom     = "point",
                 colour   = "black",
                 fill     = "black",
                 size     = 1.3,
                 position = position_dodge(0.75)) +
    
    facet_grid(abundance ~ estimand_label) +
    coord_cartesian(ylim = ylim) +
    
    scale_colour_manual(values = pal_use, name = "Model") +
    scale_shape_manual(
      values = c("n = 50" = 21, "n = 500" = 24, "n = 1000" = 22),
      name   = "Sample size"
    ) +
    
    labs(title    = title,
         subtitle = subtitle,
         x        = NULL,
         y        = "Percent bias (%)") +
    
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 11),
      legend.text      = element_text(size = 9),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y      = element_text(size = 11),
      axis.title.y     = element_text(size = 13, face = "bold"),
      strip.text       = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    )
  
  if (!is.null(fname)) {
    ggsave(fname, p, width = fig_width, height = fig_height, dpi = 300)
    cat("Saved:", basename(fname), "\n")
  }
  invisible(p)
}

# ── Estimand maps ─────────────────────────────────────────────────────────────
estimands_siv_nonmeso <- c(
  "b_pred" = "SIV: predator \u2192 prey"
)
estimands_siv_meso <- c(
  "b_pred" = "SIV: predator \u2192 prey",
  "c_pred" = "SIV: predator \u2192 meso",
  "b_meso" = "SIV: meso \u2192 prey"
)
estimands_env_nonmeso <- c(
  "a_env"         = "Env\u2081 \u2192 predator",
  "b_env"         = "Env\u2082 \u2192 prey",
  "indirect_Env1" = "Indirect: Env\u2081 \u2192 prey"
)
estimands_env_meso <- c(
  "a_env"                  = "Env\u2081 \u2192 predator",
  "b_env"                  = "Env\u2082 \u2192 prey",
  "indirect_Env1_via_pred" = "Indirect: Env\u2081 \u2192 prey"
)

estimands_int_nonmeso <- c(
  "a_int" = "Predator intercept",
  "b_int" = "Prey intercept"
)
estimands_int_meso <- c(
  "a_int" = "Predator intercept",
  "c_int" = "Meso intercept",
  "b_int" = "Prey intercept"
)

# ── ABUNDANCE figures: FA / FB / FC ──────────────────────────────────────────
cat("\n── Abundance figures ──\n")

for (is_meso in c(FALSE, TRUE)) {
  meso_tag   <- if (is_meso) "meso"  else "nonmeso"
  meso_label <- if (is_meso) "Meso"  else "Non-meso"
  
  df_base <- comp_long_clean %>%
    filter(model %in% abund_all_models,
           meso  == is_meso,
           !is.na(estimate), !is.na(truth))
  
  siv_map <- if (is_meso) estimands_siv_meso    else estimands_siv_nonmeso
  env_map <- if (is_meso) estimands_env_meso    else estimands_env_nonmeso
  int_map <- if (is_meso) estimands_int_meso    else estimands_int_nonmeso
  
  # FA — SIVs
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("SIV % bias \u2014 Abundance models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = siv_map,
    ylim         = c(-100, 100),
    fig_width    = if (is_meso) 11 else 8,
    fig_height   = 6,
    pal          = pal_12[abund_all_models],
    fname        = file.path(fig_dir, paste0("FA_SIV_abund_", meso_tag, "_LF.png"))
  )
  
  # FB — Environmental
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("Environmental % bias \u2014 Abundance models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = env_map,
    ylim         = c(-100, 100),
    fig_width    = 11,
    fig_height   = 6,
    pal          = pal_12[abund_all_models],
    fname        = file.path(fig_dir, paste0("FB_Env_abund_", meso_tag, "_LF.png"))
  )
  
  # FC — Intercepts
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("Intercept % bias \u2014 Abundance models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = int_map,
    ylim         = c(-100, 100),
    fig_width    = if (is_meso) 11 else 12,
    fig_height   = 6,
    pal          = pal_12[abund_all_models],
    fname        = file.path(fig_dir, paste0("FC_Int_abund_", meso_tag, "_LF.png"))
  )
}


# ── RN figures: FA / FB / FC ─────────────────────────────────────────────────
cat("\n── Royle-Nichols figures ──\n")

for (is_meso in c(FALSE, TRUE)) {
  meso_tag   <- if (is_meso) "meso"  else "nonmeso"
  meso_label <- if (is_meso) "Meso"  else "Non-meso"
  
  df_base <- comp_long_clean %>%
    filter(model %in% rn_all_models,
           meso  == is_meso,
           !is.na(estimate), !is.na(truth))
  
  siv_map <- if (is_meso) estimands_siv_meso    else estimands_siv_nonmeso
  env_map <- if (is_meso) estimands_env_meso    else estimands_env_nonmeso
  int_map <- if (is_meso) estimands_int_meso    else estimands_int_nonmeso
  
  # FA — SIVs
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("SIV % bias \u2014 Royle-Nichols models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = siv_map,
    ylim         = c(-100, 100),
    fig_width    = if (is_meso) 11 else 8,
    fig_height   = 6,
    pal          = pal_12[rn_all_models],
    fname        = file.path(fig_dir, paste0("FA_SIV_rn_", meso_tag, "_LF.png"))
  )
  
  # FB — Environmental
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("Environmental % bias \u2014 Royle-Nichols models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = env_map,
    ylim         = c(-100, 100),
    fig_width    = 11,
    fig_height   = 6,
    pal          = pal_12[rn_all_models],
    fname        = file.path(fig_dir, paste0("FB_Env_rn_", meso_tag, "_LF.png"))
  )
  
  # FC — Intercepts
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("Intercept % bias \u2014 Royle-Nichols models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = int_map,
    ylim         = c(-100, 100),
    fig_width    = if (is_meso) 11 else 12,
    fig_height   = 6,
    pal          = pal_12[rn_all_models],
    fname        = file.path(fig_dir, paste0("FC_Int_rn_", meso_tag, "_LF.png"))
  )
}


# ── OCCUPANCY figures: FA / FB / FC ──────────────────────────────────────────
cat("\n── Occupancy figures ──\n")

for (is_meso in c(FALSE, TRUE)) {
  meso_tag   <- if (is_meso) "meso"  else "nonmeso"
  meso_label <- if (is_meso) "Meso"  else "Non-meso"
  
  df_base <- comp_long_clean %>%
    filter(model %in% occ_all_models,
           meso  == is_meso,
           !is.na(estimate), !is.na(truth))
  
  siv_map <- if (is_meso) estimands_siv_meso    else estimands_siv_nonmeso
  env_map <- if (is_meso) estimands_env_meso    else estimands_env_nonmeso
  int_map <- if (is_meso) estimands_int_meso    else estimands_int_nonmeso
  
  # FA — SIVs
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("SIV % bias \u2014 Occupancy models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = siv_map,
    ylim         = c(-100, 100),
    fig_width    = if (is_meso) 12 else 7,
    fig_height   = 7,
    pal          = pal_12[occ_all_models],
    fname        = file.path(fig_dir, paste0("FA_SIV_occ_", meso_tag, "_LF.png"))
  )
  
  # FB — Environmental
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("Environmental % bias \u2014 Occupancy models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = env_map,
    ylim         = c(-100, 100),
    fig_width    = 12,
    fig_height   = 7,
    pal          = pal_12[occ_all_models],
    fname        = file.path(fig_dir, paste0("FB_Env_occ_", meso_tag, "_LF.png"))
  )
  
  # FC — Intercepts
  make_bias_fig_FABC(
    df           = df_base,
    title        = paste0("Intercept % bias \u2014 Occupancy models \u2014 ", meso_label, " \u2014 No landscape"),
    estimand_map = int_map,
    ylim         = c(-100, 100),
    fig_width    = if (is_meso) 12 else 12,
    fig_height   = 7,
    pal          = pal_12[occ_all_models],
    fname        = file.path(fig_dir, paste0("FC_Int_occ_", meso_tag, "_LF.png"))
  )
}

cat("\nDone. 12 figures saved to:\n", fig_dir, "\n")



# ── Summary-stats-only version of make_bias_fig_FABC ─────────────────────────
make_bias_fig_summary <- function(df, title, estimand_map,
                                  subtitle   = "Mean \u00b1 SD; shape = sample size",
                                  ylim       = c(-100, 100),
                                  fig_width  = 12,
                                  fig_height = 7,
                                  pal        = pal_12,
                                  fname      = NULL) {
  
  df <- df %>%
    filter(estimand %in% names(estimand_map),
           is.finite(percent_bias)) %>%
    mutate(
      estimand_label = factor(estimand_map[estimand], levels = estimand_map),
      method_label   = droplevels(method_label)
    )
  
  if (nrow(df) == 0) {
    warning("Empty dataframe, skipping: ", fname)
    return(invisible(NULL))
  }
  
  pal_use <- make_pal(pal, df)
  
  p <- ggplot(df, aes(x = method_label, y = percent_bias, colour = method_label)) +
    
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.6, colour = "grey40") +

    # Coloured CI bars — dodged by nSites
    stat_summary(aes(group = nSites_label),
                 fun      = mean,
                 fun.min  = function(x) mean(x) - sd(x),
                 fun.max  = function(x) mean(x) + sd(x),
                 geom     = "linerange",
                 linewidth = 0.8,
                 position = position_dodge(0.75)) +

    # Coloured summary points
    stat_summary(aes(group = nSites_label, shape = nSites_label),
                 fun      = mean,
                 geom     = "point",
                 size     = 3,
                 position = position_dodge(0.75)) +
    
    facet_grid(abundance ~ estimand_label) +
    coord_cartesian(ylim = ylim) +
    
    scale_colour_manual(values = pal_use, name = "Model") +
    scale_shape_manual(
      values = c("n = 50" = 16, "n = 500" = 17, "n = 1000" = 15),
      name   = "Sample size"
    ) + 
    
    labs(title    = title,
         subtitle = subtitle,
         x        = NULL,
         y        = "Percent bias (%)") +
    
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 11),
      legend.text      = element_text(size = 9),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y      = element_text(size = 11),
      axis.title.y     = element_text(size = 13, face = "bold"),
      strip.text       = element_text(face = "bold", size = 11),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    )
  
  if (!is.null(fname)) {
    ggsave(fname, p, width = fig_width, height = fig_height, dpi = 300)
    cat("Saved:", basename(fname), "\n")
  }
  invisible(p)
}

# ── FA SIV summary figures — Meso only ───────────────────────────────────────

# Abundance
df_abund <- comp_long_clean %>%
  filter(model %in% abund_all_models,
         meso  == TRUE,
         !is.na(estimate), !is.na(truth))

make_bias_fig_summary(
  df           = df_abund,
  title        = "SIV % bias summary \u2014 Abundance models \u2014 Meso \u2014 No landscape",
  estimand_map = estimands_siv_meso,
  ylim         = c(-100, 100),
  fig_width    = 11,
  fig_height   = 6,
  pal          = pal_12[abund_all_models],
  fname        = file.path(fig_dir, "FA_SIV_summary_abund_meso_LF.png")
)

# Royle Nichols
# RN
df_rn <- comp_long_clean %>%
  filter(model %in% rn_all_models,
         meso  == TRUE,
         !is.na(estimate), !is.na(truth))

make_bias_fig_summary(
  df           = df_rn,
  title        = "SIV % bias summary \u2014 Royle-Nichols models \u2014 Meso \u2014 No landscape",
  estimand_map = estimands_siv_meso,
  ylim         = c(-100, 100),
  fig_width    = 11,
  fig_height   = 6,
  pal          = pal_12[rn_all_models],
  fname        = file.path(fig_dir, "FA_SIV_summary_rn_meso_LF.png")
)

# Occupancy
df_occ <- comp_long_clean %>%
  filter(model %in% occ_all_models,
         meso  == TRUE,
         !is.na(estimate), !is.na(truth))

make_bias_fig_summary(
  df           = df_occ,
  title        = "SIV % bias summary \u2014 Occupancy models \u2014 Meso \u2014 No landscape",
  estimand_map = estimands_siv_meso,
  ylim         = c(-100, 100),
  fig_width    = 11,
  fig_height   = 6,
  pal          = pal_12[occ_all_models],
  fname        = file.path(fig_dir, "FA_SIV_summary_occ_meso_LF.png")
)


# ── Summary bias figure ───────────────────────────────────────────────────────
# For each rep × model × nSites: sum |percent_bias| across estimands
# Then average across nSites and reps → one value per model per abundance level

summary_estimands <- c(
  "b_pred", "b_int", "a_int", "a_env", "b_env",            # non-meso
  "c_pred", "b_meso", "c_int",                     # meso-only (ignored if absent)
  "indirect_Env1", "indirect_Env1_via_pred"
)

summary_df <- comp_long_clean %>%
  filter(model   %in% c(abund_all_models, rn_all_models, occ_all_models),
         meso    == TRUE,          # or FALSE — run separately
         estimand %in% summary_estimands,
         is.finite(percent_bias)) %>%
  group_by(seed, model, method_label, abundance, nSites) %>%
  dplyr::summarise(total_abs_bias = sum(abs(percent_bias), na.rm = TRUE) + 1,
            .groups = "drop") %>%
  mutate(
    method_label = droplevels(factor(method_labels[model], levels = method_labels)),
    model_group  = case_when(
      model %in% abund_all_models ~ "Abundance models",
      model %in% rn_all_models    ~ "Royle-Nichols models",
      model %in% occ_all_models   ~ "Occupancy models"
    ),
    model_group = factor(model_group,
                         levels = c("Abundance models",
                                    "Royle-Nichols models",
                                    "Occupancy models"))
  )

matching <- names(pal_12)[names(pal_12) %in% unique(summary_df$model)]
pal_use  <- setNames(pal_12[matching], method_labels[matching])

p_summary <- ggplot(summary_df,
                    aes(x = method_label, y = total_abs_bias,
                        colour = method_label)) +
  geom_jitter(alpha = 0.2, width = 0.15, size = 0.9) +
  stat_summary(fun      = mean,
               fun.min  = function(x) mean(x) - sd(x),
               fun.max  = function(x) mean(x) + sd(x),
               geom     = "linerange",
               colour   = "black", linewidth = 0.6) +
  stat_summary(fun    = mean, geom = "point",
               colour = "black", size = 2.5) +
  facet_grid(abundance ~ model_group, scales = "free_x", space = "free_x") + 
  scale_y_continuous(trans  = "log10",
                     breaks = c(50, 100, 500, 1000, 3000),
                     labels = scales::label_comma()) +
  coord_cartesian(ylim = c(50, 2900)) +
  scale_colour_manual(values = pal_use, name = "Model") +
  labs(title    = "Total absolute bias — all estimands — Meso — No landscape",
       subtitle = "Sum of |% bias| across SIVs, intercepts, and environmental effects per replicate\na_int excluded (truth = 0 for rare predator)",
       x = NULL, y = "Total absolute bias (%) — log10 scale") + 
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y      = element_text(size = 11),
        axis.title.y     = element_text(size = 13, face = "bold"),
        strip.text       = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(size = 9, colour = "grey40"))

ggsave(file.path(fig_dir, "F00_summary_bias_meso_LF.png"),
       p_summary, width = 12, height = 7, dpi = 300)



# ── Summary stats figure (no raw dots) ───────────────────────────────────────
p_summary_clean <- ggplot(summary_df,
                          aes(x      = method_label,
                              y      = total_abs_bias,
                              colour = method_label,
                              shape  = abundance)) +
  stat_summary(fun      = mean,
               fun.min  = function(x) mean(x) - sd(x),
               fun.max  = function(x) mean(x) + sd(x),
               geom     = "linerange",
               linewidth = 0.8,
               position = position_dodge(width = 0.6)) +
  stat_summary(fun      = mean,
               geom     = "point",
               size     = 4,
               position = position_dodge(width = 0.6)) +
  facet_wrap(~ model_group, nrow = 1, scales = "free_x") +
  coord_cartesian(ylim = c(50, 3000)) +
  scale_y_continuous(trans  = "log10",
                     breaks = c(50, 100, 500, 1000, 3000),
                     labels = scales::label_comma()) +
  scale_colour_manual(values = pal_use, name = "Model") +
  scale_shape_manual(values = c("Abundant" = 16, "Rare" = 17),
                     name   = "Scenario") +
  labs(title    = "Total absolute bias — all estimands — Meso — No landscape",
       subtitle = "Mean ± 95% CI of summed |% bias| across estimands per replicate\na_int excluded (truth = 0 for rare predator)",
       x        = NULL,
       y        = "Total absolute bias (%) — log10 scale") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "right",
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y      = element_text(size = 11),
        axis.title.y     = element_text(size = 13, face = "bold"),
        strip.text       = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(size = 9, colour = "grey40"))

ggsave(file.path(fig_dir, "F00_summary_bias_clean_meso_LF.png"),
       p_summary_clean, width = 12, height = 5, dpi = 300)


 # ── Model 6 diagnostic histogram ─────────────────────────────────────────────
p_mod6_hist <- summary_df %>%
  filter(model == "null_occ_sem") %>%
  ggplot(aes(x = total_abs_bias)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white", linewidth = 0.3) +
  geom_vline(aes(xintercept = median(total_abs_bias)),
             colour = "red", linetype = "dashed", linewidth = 0.7) +
  facet_grid(abundance ~ nSites, labeller = label_both) +
  labs(x = "Total absolute bias (%) + 1",
       y = "Count",
       title = "Model 6 (Two-step null Occ SEM) — total absolute bias distribution",
       subtitle = "Red dashed line = median; note +1 offset from log10 transformation") +
  theme_bw(base_size = 11) +
  theme(strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

plot(p_mod6_hist)

summary_df %>%
  filter(model == "null_occ_sem") %>%
  arrange(desc(total_abs_bias)) %>%
  head(10) %>%
  select(seed, abundance, nSites, total_abs_bias)


# ------------- table of mean ± SD per model × estimand × abundance × nSites ##
comp_long %>%
  filter(model %in% abund_all_models,
         meso == TRUE,
         estimand %in% c("b_pred", "c_pred", "b_meso"),
         is.finite(percent_bias)) %>%
  group_by(model, estimand, abundance, nSites) %>%
  dplyr::summarise(
    mean_bias = mean(percent_bias, na.rm = TRUE),
    sd_bias   = sd(percent_bias, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(model, estimand, abundance, nSites) %>%
  print(n = 100)

