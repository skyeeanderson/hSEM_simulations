# ============================================================
#  Extract comparison tables and Bayesian power from
#  v4 simulation replicates.
#  Run on Bunya before downloading to avoid transferring ~30GB
#
#  Output:
#    - v4_all_scenarios_comparisons.rds
#    - v4_all_scenarios_power.rds
#    - per-scenario _comparisons_only.rds and _power_only.rds
#
#  Run as: Rscript 02_extract_results.R
# ============================================================

rm(list = ls())

results_root <- "/scratch/user/uqsand24/simulations/results/"
save_path    <- "/scratch/user/uqsand24/simulations/results/"

scenarios <- c(
  "v4rn_abundant_50_meso_LT",
  "v4rn_abundant_500_meso_LT",
  "v4rn_abundant_1000_meso_LT",
  "v4rn_rare_50_meso_LT",
  "v4rn_rare_500_meso_LT",
  "v4rn_rare_1000_meso_LT",
  "v4rn_abundant_50_meso_LF",
  "v4rn_abundant_500_meso_LF",
  "v4rn_abundant_1000_meso_LF",
  "v4rn_rare_50_meso_LF",
  "v4rn_rare_500_meso_LF",
  "v4rn_rare_1000_meso_LF"
)

# NIMBLE model names as stored in fit object — must match results$ names
# 4a = coabund_nimble, 4b = occ_sem_nimble, 4c = sem_nimble
nimble_models <- c(
  coabund    = "4. Co-abundance (N-mixture)",
  sem_nimble = "5. Integrated N-mixture SEM",
  rn_nimble  = "8. Integrated RN SEM",
  occ_sem    = "12. Integrated occupancy SEM"
)

# ADD near top of extract script:
safe_val <- function(expr) {
  val <- tryCatch(expr, error = function(e) NA)
  if (is.null(val) || length(val) == 0) NA else val
}

# Map NIMBLE parameter names to comparison table estimand names
estimand_map <- c(
  a_int                  = "a_int",
  a_env                  = "a_env",
  p_pred                 = "p_pred",
  lp_pred                = "p_pred",      # logit scale — power evaluated on logit scale
  b_int                  = "b_int",
  b_env                  = "b_env",
  b_pred                 = "b_pred",
  p_prey                 = "p_prey",
  lp_prey                = "p_prey",
  c_int                  = "c_int",
  c_pred                 = "c_pred",
  b_meso                 = "b_meso",
  p_meso                 = "p_meso",
  lp_meso                = "p_meso",
  r_pred                 = "p_pred",
  r_meso                 = "p_meso",
  r_prey                 = "p_prey",
  indirect_Env1_on_prey  = "indirect_Env1",
  indirect_Env1_via_pred = "indirect_Env1_via_pred",
  indirect_Env1_via_meso = "indirect_Env1_via_meso",
  indirect_pred_via_meso = "indirect_pred_via_meso",
  total_Env1_on_prey     = "total_Env1_on_prey",
  total_pred_on_prey     = "total_pred_on_prey",
  sigma_land_pred        = "sigma_land_pred",
  sigma_land_meso        = "sigma_land_meso",
  sigma_land_prey        = "sigma_land_prey"
)

# For lp_ params, truth is on logit scale not probability scale
logit_params <- c("lp_pred", "lp_prey", "lp_meso")

all_comparisons <- list()
all_power       <- list()

for (sc in scenarios) {
  
  sc_path   <- file.path(results_root, sc)
  rds_files <- list.files(
    sc_path,
    pattern = "rep_[0-9]+\\.rds$",
    full.names = TRUE
  )
  
  cat(sprintf("\n%s: %d files found\n", sc, length(rds_files)))
  
  sc_comparisons <- list()
  sc_power       <- list()
  failed         <- 0
  
  for (f in rds_files) {
    
    rep_result <- tryCatch(readRDS(f), error = function(e) NULL)
    
    if (is.null(rep_result)) {
      failed <- failed + 1
      next
    }
    
    # ── Comparison table ─────────────────────────────────────
    comp <- tryCatch(rep_result$fit$comparison, error = function(e) NULL)
    if (is.null(comp)) {
      # Try alternative structure
      comp <- tryCatch(rep_result$comparison, error = function(e) NULL)
    }
    if (is.null(comp)) {
      failed <- failed + 1
      rm(rep_result); gc(verbose = FALSE)
      next
    }
    
    # Metadata — use tryCatch for all fields since structure may vary
    comp$seed     <- tryCatch(rep_result$seed,     error = function(e) NA)
    comp$slurm    <- tryCatch(rep_result$slurm,    error = function(e) NA)
    comp$scenario <- tryCatch(rep_result$scenario, error = function(e) sc)
    comp$nSites   <- tryCatch(rep_result$nSites,   error = function(e) NA)
    comp$meso     <- tryCatch(rep_result$meso,     error = function(e) NA)
    comp$use_landscape <- grepl("_LT$", sc)   # TRUE for LT, FALSE for LF
    
    # Species-specific mean lambda (new parameterisation)
    comp$meanLambda_pred <- tryCatch(rep_result$meanLambda_pred, error = function(e) NA)
    comp$meanLambda_prey <- tryCatch(rep_result$meanLambda_prey, error = function(e) NA)
    comp$meanLambda_meso <- tryCatch(rep_result$meanLambda_meso, error = function(e) NA)
    
    # Realised mean abundances from simulation
    comp$mean_N_pred <- tryCatch(rep_result$mean_N_pred, error=function(e) NA)
    comp$mean_N_prey <- tryCatch(rep_result$mean_N_prey, error=function(e) NA)
    comp$mean_N_meso <- tryCatch(rep_result$mean_N_meso, error=function(e) NA)
    
    sc_comparisons[[length(sc_comparisons) + 1]] <- comp
    
    # ── Bayesian power from NIMBLE summaries ─────────────────
    # Power = 95% CrI excludes zero in the correct direction
    # For lp_ params: truth converted to logit scale for sign reference
    
    truth_vec <- setNames(comp$truth, comp$estimand)
    
    for (nm in names(nimble_models)) {
      
      # Flexible access — try fit$ wrapper first, then direct
      summ <- tryCatch(
        rep_result$fit[[nm]]$summary,
        error = function(e) NULL
      )
      if (is.null(summ)) {
        summ <- tryCatch(
          rep_result[[nm]]$summary,
          error = function(e) NULL
        )
      }
      if (is.null(summ)) next
      
      power_rows <- lapply(seq_len(nrow(summ)), function(k) {
        
        param <- summ$param[k]
        q2.5  <- summ$q2.5[k]
        q97.5 <- summ$q97.5[k]
        
        estimand <- estimand_map[param]
        if (is.na(estimand)) return(NULL)
        
        truth <- truth_vec[estimand]
        if (is.na(truth)) return(NULL)
        
        # For logit-scale detection params, convert truth to logit scale
        if (param %in% logit_params) {
          truth <- qlogis(truth)   # e.g. p=0.35 → logit=-0.619
        }
        
        # Power: CrI excludes zero in correct direction
        if (truth < 0) {
          power <- as.integer(q97.5 < 0)
        } else if (truth > 0) {
          power <- as.integer(q2.5 > 0)
        } else {
          power <- NA_integer_
        }
        
        data.frame(
          seed           = safe_val(rep_result$seed),
          slurm          = safe_val(rep_result$slurm),
          scenario = if (!is.na(safe_val(rep_result$scenario))) rep_result$scenario else sc,
          nSites         = safe_val(rep_result$nSites),
          meanLambda_pred = safe_val(rep_result$meanLambda_pred),
          meanLambda_prey = safe_val(rep_result$meanLambda_prey),
          meso           = safe_val(rep_result$meso),
          model          = nm,
          model_label    = nimble_models[nm],
          param          = param,
          estimand       = estimand,
          truth          = truth,
          mean_est       = summ$mean[k],
          sd_est         = summ$sd[k],
          q2.5           = q2.5,
          q97.5          = q97.5,
          Rhat           = summ$Rhat[k],
          power          = power,
          use_landscape = grepl("_LT$", sc),
          stringsAsFactors = FALSE
        )
      })
      
      power_rows <- do.call(rbind, Filter(Negate(is.null), power_rows))
      if (!is.null(power_rows) && nrow(power_rows) > 0)
        sc_power[[length(sc_power) + 1]] <- power_rows
    }
    
    rm(rep_result)
    gc(verbose = FALSE)
  }
  
  cat(sprintf("  Successful: %d | Failed: %d\n",
              length(sc_comparisons), failed))
  
  if (length(sc_comparisons) > 0) {
    sc_df    <- do.call(rbind, sc_comparisons)
    out_file <- file.path(save_path, paste0(sc, "_comparisons_only.rds"))
    saveRDS(sc_df, out_file)
    cat(sprintf("  Saved: %s (%.1f KB)\n", basename(out_file),
                file.size(out_file) / 1024))
    all_comparisons[[sc]] <- sc_df
  }
  
  if (length(sc_power) > 0) {
    sc_power_df <- do.call(rbind, sc_power)
    out_file    <- file.path(save_path, paste0(sc, "_power_only.rds"))
    saveRDS(sc_power_df, out_file)
    cat(sprintf("  Saved: %s (%.1f KB)\n", basename(out_file),
                file.size(out_file) / 1024))
    all_power[[sc]] <- sc_power_df
  }
}

# ── Save combined files ───────────────────────────────────────────────────────
if (length(all_comparisons) > 0) {
  combined      <- do.call(rbind, all_comparisons)
  combined_file <- file.path(save_path, "v4rn_all_scenarios_comparisons.rds")
  saveRDS(combined, combined_file)
  cat(sprintf("\nComparisons saved: %s (%.1f KB)\n",
              combined_file, file.size(combined_file) / 1024))
  cat(sprintf("Total rows: %d\n", nrow(combined)))
  cat(sprintf("Scenarios: %s\n", paste(unique(combined$scenario), collapse = ", ")))
  cat(sprintf("\nReplicates per scenario (b_pred estimand):\n"))
  print(table(combined$scenario[combined$estimand == "b_pred"]))
}

if (length(all_power) > 0) {
  combined_power <- do.call(rbind, all_power)
  power_file    <- file.path(save_path, "v4rn_all_scenarios_power.rds")
  saveRDS(combined_power, power_file)
  cat(sprintf("\nPower saved: %s (%.1f KB)\n",
              power_file, file.size(power_file) / 1024))
  cat(sprintf("Total power rows: %d\n", nrow(combined_power)))
  
  cat(sprintf("\nPower summary (b_pred, by model and scenario):\n"))
  bp <- combined_power[combined_power$estimand == "b_pred", ]
  if (nrow(bp) > 0) {
    print(
      tapply(bp$power, list(bp$scenario, bp$model_label), mean, na.rm = TRUE),
      digits = 2
    )
  }
}
