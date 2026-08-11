# ============================================================
#  Predator-prey simulation — combined abundance + occupancy (v4)
#  (Skye Anderson, UQ PhD Ch3)
#
#  SAME SIMULATION DATASET applied to all 12 models.
#
#  DATA-GENERATING MODEL:
#    Abundance DGP: N_i ~ Poisson(lambda_i), log link.
#    SIV (b_pred) on LOG scale in DGP.
#    Binary detection histories Y_ij = ifelse(C_ij > 0, 1, 0)
#    shared by Royle-Nichols and occupancy models.
#    Count histories C_ij shared by N-mixture models.
#    Landscape structure toggled via LANDSCAPE env var.
#
#  FITTED MODELS:
#
#    Abundance — N-mixture (b_pred on log scale):
#    Model 1: Max count + SEM (lmer)
#    Model 2: Null N-mixture BLUPs + SEM (lmer)
#    Model 3: Full N-mixture BLUPs + SEM (lmer)
#    Model 4: Co-abundance N-mixture (NIMBLE)
#    Model 5: Integrated N-mixture SEM (NIMBLE)
#
#    Abundance — Royle-Nichols (b_pred on log scale):
#    Model 6: Null Royle-Nichols BLUPs + SEM (lmer)
#    Model 7: Full Royle-Nichols BLUPs + SEM (lmer)
#    Model 8: Integrated Royle-Nichols SEM (NIMBLE)
#
#    Occupancy applied to binarised Y_ij from abundance DGP:
#    Model 9:  Naïve detection rate + SEM (lmer, probability scale)
#    Model 10: Null occupancy model BLUPs + SEM (lmer, probability scale)
#    Model 11: Full occupancy model predicted psi + SEM (lmer, probability scale)
#    Model 12: Integrated Bayesian occupancy SEM (NIMBLE, logit scale)
#
#  Note: b_pred in Models 1–8 is on the log scale, matching the DGP.
#        Models 6–8 use binary Y_ij but maintain a Poisson abundance
#        process, recovering the same estimand as Models 1–5.
#        b_pred in Models 9–11 is on the probability scale (lmer Gaussian);
#        b_pred in Model 12 is on the logit scale. Both differ from the
#        log-scale truth — this is expected and part of the comparison.
# ============================================================

rm(list = ls())

# ── Settings ──────────────────────────────────────────────────────────────────
setting       <- Sys.getenv("SETTING");  if (setting == "")    setting  <- "LOCAL"
scenario      <- Sys.getenv("SCENARIO"); if (scenario == "")   scenario <- "rare_50_meso"
landscape_str <- Sys.getenv("LANDSCAPE"); if (landscape_str == "") landscape_str <- "FALSE"
use_landscape <- !tolower(landscape_str) %in% c("false", "0", "no")
rep_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); if (is.na(rep_id)) rep_id <- 1
seed          <- rep_id

cat(sprintf("Setting: %s | Scenario: %s | Landscape: %s\n",
            setting, scenario, use_landscape))

# ── Unique temp dir for NIMBLE compilation (HPC only) ────────────────────────
if (setting == "HPC") {
  task_dir <- file.path(tempdir(),
                        paste0("task_", Sys.getenv("SLURM_ARRAY_JOB_ID"),
                               "_", Sys.getenv("SLURM_ARRAY_TASK_ID")))
  dir.create(task_dir, recursive=TRUE, showWarnings=FALSE)
  setwd(task_dir)
}

# ── Helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ── Scenario parameters ───────────────────────────────────────────────────────
scenario_params <- list(
  
  # # ── Abundant (prey=10, meso=4, pred=2) ──────────────────────────────────
  # abundant_1000 = list(
  #   nSites=1000, nVisits=10,
  #   meanLambda_pred=2.0, meanLambda_prey=10.0, meanLambda_meso=4.0,
  #   run_meso=FALSE, meso_cpred=-0.2, meso_bmeso=-0.15),
  # abundant_500 = list(
  #   nSites=500,  nVisits=10,
  #   meanLambda_pred=2.0, meanLambda_prey=10.0, meanLambda_meso=4.0,
  #   run_meso=FALSE, meso_cpred=-0.2, meso_bmeso=-0.15),
  # abundant_50 = list(
  #   nSites=50,   nVisits=10,
  #   meanLambda_pred=2.0, meanLambda_prey=10.0, meanLambda_meso=4.0,
  #   run_meso=FALSE, meso_cpred=-0.2, meso_bmeso=-0.15),
  # 
  # # ── Rare (prey=5, meso=2, pred=1) ───────────────────────────────────────
  # rare_1000 = list(
  #   nSites=1000, nVisits=10,
  #   meanLambda_pred=0.5, meanLambda_prey=2, meanLambda_meso=1.5,
  #   run_meso=FALSE, meso_cpred=-0.2, meso_bmeso=-0.15),
  # rare_500 = list(
  #   nSites=500,  nVisits=10,
  #   meanLambda_pred=0.5, meanLambda_prey=2, meanLambda_meso=1.5,
  #   run_meso=FALSE, meso_cpred=-0.2, meso_bmeso=-0.15),
  # rare_50 = list(
  #   nSites=50,   nVisits=10,
  #   meanLambda_pred=0.5, meanLambda_prey=2, meanLambda_meso=1.5,
  #   run_meso=FALSE, meso_cpred=-0.2, meso_bmeso=-0.15),
  
  # ── Abundant + meso ──────────────────────────────────────────────────────
  abundant_1000_meso = list(
    nSites=1000, nVisits=10,
    meanLambda_pred=2.0, meanLambda_prey=10.0, meanLambda_meso=4.0,
    run_meso=TRUE, meso_cpred=-0.2, meso_bmeso=-0.15),
  abundant_500_meso = list(
    nSites=500,  nVisits=10,
    meanLambda_pred=2.0, meanLambda_prey=10.0, meanLambda_meso=4.0,
    run_meso=TRUE, meso_cpred=-0.2, meso_bmeso=-0.15),
  abundant_50_meso = list(
    nSites=50,   nVisits=10,
    meanLambda_pred=2.0, meanLambda_prey=10.0, meanLambda_meso=4.0,
    run_meso=TRUE, meso_cpred=-0.2, meso_bmeso=-0.15),
  
  # ── Rare + meso ──────────────────────────────────────────────────────────
  rare_1000_meso = list(
    nSites=1000, nVisits=10,
    meanLambda_pred=0.5, meanLambda_prey=2, meanLambda_meso=1.5,
    run_meso=TRUE, meso_cpred=-0.2, meso_bmeso=-0.15),
  rare_500_meso = list(
    nSites=500,  nVisits=10,
    meanLambda_pred=0.5, meanLambda_prey=2, meanLambda_meso=1.5,
    run_meso=TRUE, meso_cpred=-0.2, meso_bmeso=-0.15),
  rare_50_meso = list(
    nSites=50,   nVisits=10,
    meanLambda_pred=0.5, meanLambda_prey=2, meanLambda_meso=1.5,
    run_meso=TRUE, meso_cpred=-0.2, meso_bmeso=-0.15)
)

if (!scenario %in% names(scenario_params))
  stop(sprintf("Unknown scenario: '%s'", scenario))

sc               <- scenario_params[[scenario]]
nSites           <- sc$nSites
nVisits          <- sc$nVisits
meanLambda_pred  <- sc$meanLambda_pred
meanLambda_prey  <- sc$meanLambda_prey
meanLambda_meso  <- sc$meanLambda_meso
run_meso         <- sc$run_meso
meso_cpred       <- sc$meso_cpred
meso_bmeso       <- sc$meso_bmeso

n_landscapes <- 10

# ── Libraries ─────────────────────────────────────────────────────────────────
if (setting == "LOCAL") {
  library(lme4); library(unmarked); library(nimble); library(coda)
} else {
  for (pkg in c("lme4","unmarked","nimble","coda"))
    library(pkg, lib.loc="/home/uqsand24/R/x86_64-pc-linux-gnu-library/4.4",
            character.only=TRUE)
}

# ── MCMC settings ─────────────────────────────────────────────────────────────
if (setting == "LOCAL") {
  nimble_iter=5000; nimble_burnin=1000; nimble_chains=3; nimble_thin=2
} else {
  nimble_iter=30000; nimble_burnin=10000; nimble_chains=3; nimble_thin=5
}

# ── Seed and replicate ────────────────────────────────────────────────────────
if (setting == "LOCAL") { slurm <- 1; seed <- 42 } else {
  slurm <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if (is.na(slurm))
    stop("SLURM_ARRAY_TASK_ID not found")
  seed  <- slurm
}

# ── Save paths ────────────────────────────────────────────────────────────────
land_tag <- ifelse(use_landscape, "LT", "LF")

# ============================================================
#  SIMULATE DATA — ABUNDANCE DGP
#  Same as v3. Generates count data (C) and binarised
#  detection histories (Y) for both model families.
#  Landscape structure toggled by use_landscape.
# ============================================================

simulate_predprey <- function(
    nSites, nVisits=10, seed=NULL,
    n_landscapes=10, n_both_landscapes=4,
    n_pred_only=3, n_prey_only=3,
    a_int=log(2.0), a_env=0.8,
    b_int=log(10), b_env=0.6, b_pred=-0.25,
    meso, c_int=log(4.0), c_pred=-0.2, b_meso=-0.15,
    p_pred=0.35, p_prey=0.65, p_meso=0.50,
    sigma_land_pred=0.3, sigma_land_prey=0.3, sigma_land_meso=0.3,
    meanLambda_pred=NULL, meanLambda_prey=NULL, meanLambda_meso=NULL,
    use_landscape
) {
  if (!is.null(meanLambda_pred)) a_int <- log(meanLambda_pred)
  if (!is.null(meanLambda_prey)) b_int <- log(meanLambda_prey)
  if (!is.null(meanLambda_meso)) c_int <- log(meanLambda_meso)
  if (!is.null(seed)) set.seed(seed)
  
  if (use_landscape) {
    stopifnot(nSites %% n_landscapes == 0)
    stopifnot(n_both_landscapes + n_pred_only + n_prey_only == n_landscapes)
    sites_per_land <- nSites / n_landscapes
    landscape  <- rep(1:n_landscapes, each=sites_per_land)
    land_type  <- c(rep("both",n_both_landscapes),
                    rep("pred_only",n_pred_only),
                    rep("prey_only",n_prey_only))
    b_land_pred <- rnorm(n_landscapes, 0, sigma_land_pred)
    b_land_prey <- rnorm(n_landscapes, 0, sigma_land_prey)
    if (meso) b_land_meso <- rnorm(n_landscapes, 0, sigma_land_meso)
    z_land_pred <- ifelse(land_type[landscape] %in% c("both","pred_only"), 1L, 0L)
    z_land_prey <- ifelse(land_type[landscape] %in% c("both","prey_only"), 1L, 0L)
  } else {
    landscape   <- rep(1L, nSites)  # single pseudo-landscape
    land_type   <- "both"
    b_land_pred <- 0
    b_land_prey <- 0
    if (meso) b_land_meso <- 0
    z_land_pred <- rep(1L, nSites)
    z_land_prey <- rep(1L, nSites)
    n_landscapes <- 1L
  }
  
  Env1 <- rnorm(nSites); Env2 <- rnorm(nSites)
  
  # ── Predator ────────────────────────────────────────────────────────────────
  lambda_pred <- exp(a_int + a_env*Env1 + b_land_pred[landscape])
  N_pred      <- ifelse(z_land_pred==1L, rpois(nSites, lambda_pred), 0L)
  C_pred      <- matrix(rbinom(nSites*nVisits, rep(N_pred,nVisits), p_pred),
                        nrow=nSites, ncol=nVisits)
  Y_pred      <- ifelse(C_pred > 0L, 1L, 0L)
  
  # ── Mesopredator ─────────────────────────────────────────────────────────────
  if (meso) {
    lambda_meso <- exp(c_int + c_pred*N_pred + b_land_meso[landscape])
    N_meso      <- rpois(nSites, lambda_meso)
    C_meso      <- matrix(rbinom(nSites*nVisits, rep(N_meso,nVisits), p_meso),
                          nrow=nSites, ncol=nVisits)
    Y_meso      <- ifelse(C_meso > 0L, 1L, 0L)
  }
  
  # ── Prey ────────────────────────────────────────────────────────────────────
  if (meso) {
    lambda_prey <- exp(b_int + b_env*Env2 +
                         b_pred*N_pred + b_meso*N_meso +
                         b_land_prey[landscape])
  } else {
    lambda_prey <- exp(b_int + b_env*Env2 +
                         b_pred*N_pred +
                         b_land_prey[landscape])
  }
  N_prey <- ifelse(z_land_prey==1L, rpois(nSites, lambda_prey), 0L)
  C_prey <- matrix(rbinom(nSites*nVisits, rep(N_prey,nVisits), p_prey),
                   nrow=nSites, ncol=nVisits)
  Y_prey <- ifelse(C_prey > 0L, 1L, 0L)
  
  # ── Effect decomposition (log-scale path products) ───────────────────────────
  indirect_Env1_via_pred <- a_env * b_pred
  total_Env1_on_prey     <- indirect_Env1_via_pred
  if (meso) {
    indirect_Env1_via_meso <- a_env * c_pred * b_meso
    indirect_pred_via_meso <- c_pred * b_meso
    total_pred_on_prey     <- b_pred + indirect_pred_via_meso
    total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
  }
  
  sim_data <- list(
    nSites=nSites, nVisits=nVisits,
    n_landscapes=n_landscapes, landscape=landscape,
    land_type=land_type, use_landscape=use_landscape,
    z_land_pred=z_land_pred, z_land_prey=z_land_prey,
    N_pred=N_pred, N_prey=N_prey,
    C_pred=C_pred, C_prey=C_prey,
    Y_pred=Y_pred, Y_prey=Y_prey,
    Env1=Env1, Env2=Env2)
  if (meso) {
    sim_data$N_meso=N_meso; sim_data$C_meso=C_meso; sim_data$Y_meso=Y_meso }
  
  truth <- list(
    a_int=a_int, a_env=a_env, b_int=b_int, b_env=b_env, b_pred=b_pred,
    p_pred=p_pred, p_prey=p_prey,
    sigma_land_pred=sigma_land_pred, sigma_land_prey=sigma_land_prey,
    b_land_pred=b_land_pred, b_land_prey=b_land_prey,
    N_pred=N_pred, N_prey=N_prey,
    indirect_Env1_via_pred=indirect_Env1_via_pred,
    total_Env1_on_prey=total_Env1_on_prey)
  if (meso) {
    truth$c_int=c_int; truth$c_pred=c_pred; truth$b_meso=b_meso
    truth$p_meso=p_meso; truth$N_meso=N_meso
    truth$sigma_land_meso=sigma_land_meso; truth$b_land_meso=b_land_meso
    truth$indirect_pred_via_meso=indirect_pred_via_meso
    truth$total_pred_on_prey=total_pred_on_prey
    truth$indirect_Env1_via_meso=indirect_Env1_via_meso }
  
  list(sim_data=sim_data, truth=truth)
}

# ============================================================
#  FIT ALL MODELS
# ============================================================

fit_all_models <- function(sim_out,
                           nimble_iter=30000, nimble_burnin=10000,
                           nimble_chains=3, nimble_thin=5,
                           meso, use_landscape) {
  
  dat          <- sim_out$sim_data
  truth        <- sim_out$truth
  nSites       <- dat$nSites; nVisits <- dat$nVisits
  n_landscapes <- dat$n_landscapes
  landscape    <- dat$landscape; land_f <- factor(landscape)
  Env1         <- dat$Env1; Env2 <- dat$Env2
  C_pred       <- dat$C_pred; C_prey <- dat$C_prey
  Y_pred       <- dat$Y_pred; Y_prey <- dat$Y_prey
  z_land_pred  <- dat$z_land_pred; z_land_prey <- dat$z_land_prey
  if (meso) { C_meso <- dat$C_meso; Y_meso <- dat$Y_meso }
  
  # Landscape-type presence indices (used for subsetting)
  idx_pred <- which(z_land_pred==1); idx_prey <- which(z_land_prey==1)
  
  results <- list()
  
  # ── Shared helpers ───────────────────────────────────────────────────────────
  post_summary <- function(samples, params) {
    do.call(rbind, lapply(params, function(p) {
      draws <- as.vector(as.matrix(samples[,p]))
      data.frame(param=p, mean=mean(draws), sd=sd(draws),
                 q2.5=quantile(draws,.025), q97.5=quantile(draws,.975),
                 Rhat=tryCatch(gelman.diag(samples[,p])$psrf[1],
                               error=function(e) NA))
    }))
  }
  gpm <- function(samps,p) tryCatch(mean(as.vector(as.matrix(samps[,p]))),
                                    error=function(e) NA_real_)
  safe_get <- function(lst,nm,def=NA_real_) {
    v <- tryCatch(lst[[nm]], error=function(e) def)
    if (is.null(v)||length(v)==0) def else as.numeric(v) }
  
  # ── NIMBLE helper: compile and run ───────────────────────────────────────────
  run_nimble <- function(code, constants, data, inits, monitors,
                         iter, burnin, chains, thin) {
    m    <- suppressWarnings(nimbleModel(code, constants, data, inits))
    cm   <- compileNimble(m)
    conf <- configureMCMC(cm, monitors=monitors)
    mc   <- buildMCMC(conf)
    cmc  <- compileNimble(mc, project=m)
    runMCMC(cmc, niter=iter, nburnin=burnin, nchains=chains, thin=thin,
            samplesAsCodaMCMC=TRUE)
  }
  
  # ── Abundance lmer SEM helper ─────────────────────────────────────────────
  # Models 1–3: regresses count-scale N_hat on N_hat.
  # use_landscape controls (1|landscape) inclusion.
  # TO:
  fit_abund_sem <- function(N_hat_pred, N_hat_prey, label,
                            N_hat_meso=NULL) {
    # Round to integers, floor at 1 (required for Poisson GLM)
    N_hat_pred <- pmax(round(N_hat_pred), 1L)
    N_hat_prey <- pmax(round(N_hat_prey), 1L)
    if (!is.null(N_hat_meso)) N_hat_meso <- pmax(round(N_hat_meso), 1L)
    
    fit_glm <- function(formula_str, data) {
      if (use_landscape)
        tryCatch(lme4::glmer(as.formula(formula_str), data=data, family=poisson),
                 error=function(e) NULL)
      else
        tryCatch(glm(as.formula(sub("\\+ \\(1\\|landscape\\)","",formula_str)),
                     data=data, family=poisson), error=function(e) NULL)
    }
    if (!is.null(N_hat_meso)) {
      df_pred <- data.frame(Npred=N_hat_pred[idx_pred], Env1=Env1[idx_pred],
                            landscape=land_f[idx_pred])
      df_meso <- data.frame(Nmeso=N_hat_meso, Npred=N_hat_pred,
                            Env1=Env1, landscape=land_f)
      df_prey <- data.frame(Nprey=N_hat_prey[idx_prey], Npred=N_hat_pred[idx_prey],
                            Nmeso=N_hat_meso[idx_prey],
                            Env2=Env2[idx_prey], landscape=land_f[idx_prey])
      m1 <- fit_glm("Npred ~ Env1 + (1|landscape)", df_pred)
      m2 <- fit_glm("Nmeso ~ Env1 + Npred + (1|landscape)", df_meso)
      m3 <- fit_glm("Nprey ~ Env2 + Npred + Nmeso + (1|landscape)", df_prey)
      if (is.null(m1)||is.null(m2)||is.null(m3))
        return(list(label=label, a_int=NA,c_int=NA,b_int=NA,a_env=NA,b_env=NA,
                    b_pred=NA,c_pred=NA,b_meso=NA,indirect_Env1=NA,
                    indirect_pred_via_meso=NA,total_pred_on_prey=NA,
                    indirect_Env1_via_meso=NA,total_Env1_on_prey=NA))
      cf1 <- coef(m1); cf2 <- coef(m2); cf3 <- coef(m3)
      if (use_landscape) { cf1 <- lme4::fixef(m1); cf2 <- lme4::fixef(m2)
      cf3 <- lme4::fixef(m3) }
      bE <- cf1["Env1"]; bpp <- cf3["Npred"]; cpm <- cf2["Npred"]; bmp <- cf3["Nmeso"]
      return(list(label=label, a_int=cf1["(Intercept)"], c_int=cf2["(Intercept)"],
                  b_int=cf3["(Intercept)"], a_env=bE, b_env=cf3["Env2"],
                  b_pred=bpp, c_pred=cpm, b_meso=bmp,
                  indirect_Env1=bE*bpp, indirect_pred_via_meso=cpm*bmp,
                  total_pred_on_prey=bpp+cpm*bmp,
                  indirect_Env1_via_meso=bE*cpm*bmp,
                  total_Env1_on_prey=bE*bpp+bE*cpm*bmp))
    } else {
      df_pred <- data.frame(Npred=N_hat_pred[idx_pred], Env1=Env1[idx_pred],
                            landscape=land_f[idx_pred])
      df_prey <- data.frame(Nprey=N_hat_prey[idx_prey], Npred=N_hat_pred[idx_prey],
                            Env2=Env2[idx_prey], landscape=land_f[idx_prey])
      m1 <- fit_glm("Npred ~ Env1 + (1|landscape)", df_pred)
      m2 <- fit_glm("Nprey ~ Env2 + Npred + (1|landscape)", df_prey)
      if (is.null(m1)||is.null(m2))
        return(list(label=label,a_int=NA,b_int=NA,a_env=NA,b_env=NA,
                    b_pred=NA,indirect_Env1=NA))
      cf1 <- coef(m1); cf2 <- coef(m2)
      if (use_landscape) { cf1 <- lme4::fixef(m1); cf2 <- lme4::fixef(m2) }
      return(list(label=label, a_int=cf1["(Intercept)"],
                  b_int=cf2["(Intercept)"], a_env=cf1["Env1"],
                  b_env=cf2["Env2"], b_pred=cf2["Npred"],
                  indirect_Env1=cf1["Env1"]*cf2["Npred"]))
    }
  }
  
  # ── Occupancy lmer SEM helper ─────────────────────────────────────────────
  # Models 6–8: regresses psi estimates (probability scale).
  fit_occ_sem <- function(psi_pred, psi_prey, label, psi_meso=NULL) {
    fit_lm <- function(formula_str, data) {
      if (use_landscape)
        tryCatch(lme4::lmer(as.formula(formula_str), data=data),
                 error=function(e) NULL)
      else
        tryCatch(lm(as.formula(sub("\\+ \\(1\\|landscape\\)","",formula_str)),
                    data=data), error=function(e) NULL)
    }
    if (sd(psi_pred, na.rm=TRUE) < 0.01) {
      if (!is.null(psi_meso))
        return(list(label=label, a_int=NA, c_int=NA, b_int=NA, a_env=NA, b_env=NA,
                    b_pred=NA, c_pred=NA, b_meso=NA, indirect_Env1=NA,
                    indirect_pred_via_meso=NA, total_pred_on_prey=NA,
                    indirect_Env1_via_meso=NA, total_Env1_on_prey=NA))
      else
        return(list(label=label, a_int=NA, b_int=NA, a_env=NA,
                    b_env=NA, b_pred=NA, indirect_Env1=NA))
    }
    if (!is.null(psi_meso)) {
      df_pred <- data.frame(psi_pred=psi_pred[idx_pred], Env1=Env1[idx_pred],
                            landscape=land_f[idx_pred])
      df_meso <- data.frame(psi_meso=psi_meso, psi_pred=psi_pred,
                            Env1=Env1, landscape=land_f)
      df_prey <- data.frame(psi_prey=psi_prey[idx_prey], psi_pred=psi_pred[idx_prey],
                            psi_meso=psi_meso[idx_prey],
                            Env2=Env2[idx_prey], landscape=land_f[idx_prey])
      m1 <- fit_lm("psi_pred ~ Env1 + (1|landscape)", df_pred)
      m2 <- fit_lm("psi_meso ~ Env1 + psi_pred + (1|landscape)", df_meso)
      m3 <- fit_lm("psi_prey ~ Env2 + psi_pred + psi_meso + (1|landscape)", df_prey)
      if (is.null(m1)||is.null(m2)||is.null(m3))
        return(list(label=label,a_int=NA,c_int=NA,b_int=NA,a_env=NA,b_env=NA,
                    b_pred=NA,c_pred=NA,b_meso=NA,indirect_Env1=NA,
                    indirect_pred_via_meso=NA,total_pred_on_prey=NA,
                    indirect_Env1_via_meso=NA,total_Env1_on_prey=NA))
      cf1 <- coef(m1); cf2 <- coef(m2); cf3 <- coef(m3)
      if (use_landscape) { cf1 <- lme4::fixef(m1); cf2 <- lme4::fixef(m2)
      cf3 <- lme4::fixef(m3) }
      bE <- cf1["Env1"]; bpp <- cf3["psi_pred"]; cpm <- cf2["psi_pred"]
      bmp <- cf3["psi_meso"]
      return(list(label=label, a_int=cf1["(Intercept)"], c_int=cf2["(Intercept)"],
                  b_int=cf3["(Intercept)"], a_env=bE, b_env=cf3["Env2"],
                  b_pred=bpp, c_pred=cpm, b_meso=bmp,
                  indirect_Env1=bE*bpp, indirect_pred_via_meso=cpm*bmp,
                  total_pred_on_prey=bpp+cpm*bmp,
                  indirect_Env1_via_meso=bE*cpm*bmp,
                  total_Env1_on_prey=bE*bpp+bE*cpm*bmp))
    } else {
      df_pred <- data.frame(psi_pred=psi_pred[idx_pred], Env1=Env1[idx_pred],
                            landscape=land_f[idx_pred])
      df_prey <- data.frame(psi_prey=psi_prey[idx_prey], psi_pred=psi_pred[idx_prey],
                            Env2=Env2[idx_prey], landscape=land_f[idx_prey])
      m1 <- fit_lm("psi_pred ~ Env1 + (1|landscape)", df_pred)
      m2 <- fit_lm("psi_prey ~ Env2 + psi_pred + (1|landscape)", df_prey)
      cat("sd(psi_prey):", sd(df_prey$psi_prey, na.rm=TRUE), "\n")
      cat("sd(psi_pred in prey df):", sd(df_prey$psi_pred, na.rm=TRUE), "\n")
      cat("sd(psi_meso in prey df):", sd(df_prey$psi_meso, na.rm=TRUE), "\n")
      cat("any NA in df_prey:", anyNA(df_prey), "\n")
      if (is.null(m1)||is.null(m2))
        return(list(label=label,a_int=NA,b_int=NA,a_env=NA,
                    b_env=NA,b_pred=NA,indirect_Env1=NA))
      cf1 <- coef(m1); cf2 <- coef(m2)
      if (use_landscape) { cf1 <- lme4::fixef(m1); cf2 <- lme4::fixef(m2) }
      return(list(label=label, a_int=cf1["(Intercept)"],
                  b_int=cf2["(Intercept)"], a_env=cf1["Env1"],
                  b_env=cf2["Env2"], b_pred=cf2["psi_pred"],
                  indirect_Env1=cf1["Env1"]*cf2["psi_pred"]))
    }
  }
  
  # ── unmarked siteCovs ────────────────────────────────────────────────────────
  if (use_landscape) {
    sc_pred <- data.frame(Env1=Env1[idx_pred], landscape=droplevels(land_f[idx_pred]))
    sc_prey <- data.frame(Env2=Env2[idx_prey], landscape=droplevels(land_f[idx_prey]))
  } else {
    sc_pred <- data.frame(Env1=Env1[idx_pred])
    sc_prey <- data.frame(Env2=Env2[idx_prey])
  }
  
  # ============================================================
  #  MODELS 1–3: ABUNDANCE TWO-STAGE (lmer SEM)
  # ============================================================
  
  # ── Model 1: Max count SEM ────────────────────────────────────────────────
  message("\n--- Model 1: Max count SEM ---")
  maxC_pred <- numeric(nSites); maxC_pred[idx_pred] <- apply(C_pred[idx_pred,,drop=F],1,max)
  maxC_prey <- numeric(nSites); maxC_prey[idx_prey] <- apply(C_prey[idx_prey,,drop=F],1,max)
  if (meso) {
    maxC_meso <- apply(C_meso, 1, max)
    results$psem_maxcount <- fit_abund_sem(maxC_pred, maxC_prey, "max count SEM",
                                           N_hat_meso=maxC_meso)
  } else {
    results$psem_maxcount <- fit_abund_sem(maxC_pred, maxC_prey, "max count SEM")
  }
  
  # ── Models 2–3: N-mixture BLUPs / predicted SEM ───────────────────────────
  message("\n--- Models 2–3: N-mixture BLUPs + SEM ---")
  K_pred <- max(C_pred) * 2 + 10
  K_prey <- max(C_prey) * 2 + 10
  
  # Null N-mixture
  umf_pred_null <- if (use_landscape)
    unmarkedFramePCount(y=C_pred[idx_pred,], siteCovs=sc_pred) else
      unmarkedFramePCount(y=C_pred[idx_pred,], siteCovs=sc_pred)
  umf_prey_null <- unmarkedFramePCount(y=C_prey[idx_prey,], siteCovs=sc_prey)
  
  nm_pred_null <- suppressWarnings(tryCatch(
    pcount(~1~1, data=umf_pred_null, K=K_pred), error=function(e) NULL))
  nm_prey_null <- suppressWarnings(tryCatch(
    pcount(~1~1, data=umf_prey_null, K=K_prey), error=function(e) NULL))
  
  Nhat_pred_null <- numeric(nSites)
  Nhat_prey_null <- numeric(nSites)
  Nhat_pred_null[idx_pred] <- if (!is.null(nm_pred_null))
    bup(ranef(nm_pred_null), stat="mean") else rep(NA, length(idx_pred))
  Nhat_prey_null[idx_prey] <- if (!is.null(nm_prey_null))
    bup(ranef(nm_prey_null), stat="mean") else rep(NA, length(idx_prey))
  
  # Full N-mixture
  full_formula_pred <- if (use_landscape) ~1~Env1+landscape else ~1~Env1
  full_formula_prey <- if (use_landscape) ~1~Env2+landscape else ~1~Env2
  
  nm_pred_full <- suppressWarnings(
    tryCatch(
      pcount(full_formula_pred, data=umf_pred_null, K=K_pred),
      error=function(e) NULL))
  
  nm_prey_full <- suppressWarnings(
    tryCatch(
      pcount(full_formula_prey, data=umf_prey_null, K=K_prey),
      error=function(e) NULL))
  
  Nhat_pred_full <- numeric(nSites)
  Nhat_prey_full <- numeric(nSites)
  Nhat_pred_full[idx_pred] <- if (!is.null(nm_pred_full))
    bup(ranef(nm_pred_full), stat="mean") else rep(NA, length(idx_pred))
  Nhat_prey_full[idx_prey] <- if (!is.null(nm_prey_full))
    bup(ranef(nm_prey_full), stat="mean") else rep(NA, length(idx_prey))
  
  if (meso) {
    K_meso <- max(C_meso) * 2 + 10
    sc_meso <- if (use_landscape) data.frame(Env1=Env1, landscape=land_f) else
      data.frame(Env1=Env1)
    full_formula_meso <- if (use_landscape) ~1~Env1+landscape else ~1~Env1
    umf_meso_null <- unmarkedFramePCount(y=C_meso, siteCovs=sc_meso)
    nm_meso_null  <- suppressWarnings(tryCatch(
      pcount(~1~1, data=umf_meso_null, K=K_meso), error=function(e) NULL))
    nm_meso_full <- suppressWarnings(
      tryCatch(
        pcount(full_formula_meso, data=umf_meso_null, K=K_meso),
        error=function(e) NULL))
    Nhat_meso_null <- if (!is.null(nm_meso_null))
      bup(ranef(nm_meso_null), stat="mean") else rep(NA, nSites)
    Nhat_meso_full <- if (!is.null(nm_meso_full))
      bup(ranef(nm_meso_full), stat="mean") else rep(NA, nSites)
    
    results$psem_blup_null <- fit_abund_sem(Nhat_pred_null, Nhat_prey_null,
                                            "null BLUP SEM", N_hat_meso=Nhat_meso_null)
    results$psem_blup_full <- fit_abund_sem(Nhat_pred_full, Nhat_prey_full,
                                            "full BLUP SEM", N_hat_meso=Nhat_meso_full)
  } else {
    results$psem_blup_null <- fit_abund_sem(Nhat_pred_null, Nhat_prey_null, "null BLUP SEM")
    results$psem_blup_full <- fit_abund_sem(Nhat_pred_full, Nhat_prey_full, "full BLUP SEM")
  }
  
  # ============================================================
  #
  # MODELS 6–7: ROYLE-NICHOLS TWO-STAGE (BUPs → lmer SEM)
  #
  # ============================================================
  
  message("\n--- Models 6–7: Royle-Nichols BLUPs + SEM ---")
  
  umf_pred_rn_null <- unmarkedFrameOccu(y=Y_pred[idx_pred,,drop=F], siteCovs=sc_pred)
  umf_prey_rn_null <- unmarkedFrameOccu(y=Y_prey[idx_prey,,drop=F], siteCovs=sc_prey)
  
  rn_pred_null <- suppressWarnings(tryCatch(
    occuRN(~1~1, data=umf_pred_rn_null), error=function(e) NULL))
  rn_prey_null <- suppressWarnings(tryCatch(
    occuRN(~1~1, data=umf_prey_rn_null), error=function(e) NULL))
  
  Nhat_pred_rn_null <- numeric(nSites)
  Nhat_prey_rn_null <- numeric(nSites)
  Nhat_pred_rn_null[idx_pred] <- if (!is.null(rn_pred_null))
    bup(ranef(rn_pred_null), stat="mean") else rep(NA, length(idx_pred))
  Nhat_prey_rn_null[idx_prey] <- if (!is.null(rn_prey_null))
    bup(ranef(rn_prey_null), stat="mean") else rep(NA, length(idx_prey))
  
  full_formula_pred_rn <- if (use_landscape) ~1~Env1+landscape else ~1~Env1
  full_formula_prey_rn <- if (use_landscape) ~1~Env2+landscape else ~1~Env2
  
  rn_pred_full <- suppressWarnings(tryCatch(
    occuRN(full_formula_pred_rn, data=umf_pred_rn_null), error=function(e) NULL))
  rn_prey_full <- suppressWarnings(tryCatch(
    occuRN(full_formula_prey_rn, data=umf_prey_rn_null), error=function(e) NULL))
  
  Nhat_pred_rn_full <- numeric(nSites)
  Nhat_prey_rn_full <- numeric(nSites)
  Nhat_pred_rn_full[idx_pred] <- if (!is.null(rn_pred_full))
    bup(ranef(rn_pred_full), stat="mean") else rep(NA, length(idx_pred))
  Nhat_prey_rn_full[idx_prey] <- if (!is.null(rn_prey_full))
    bup(ranef(rn_prey_full), stat="mean") else rep(NA, length(idx_prey))
  
  if (meso) {
    sc_meso_rn <- if (use_landscape) data.frame(Env1=Env1, landscape=land_f) else
      data.frame(Env1=Env1)
    full_formula_meso_rn <- if (use_landscape) ~1~Env1+landscape else ~1~Env1
    umf_meso_rn_null <- unmarkedFrameOccu(y=Y_meso, siteCovs=sc_meso_rn)
    rn_meso_null <- suppressWarnings(tryCatch(
      occuRN(~1~1, data=umf_meso_rn_null), error=function(e) NULL))
    rn_meso_full <- suppressWarnings(tryCatch(
      occuRN(full_formula_meso_rn, data=umf_meso_rn_null), error=function(e) NULL))
    Nhat_meso_rn_null <- if (!is.null(rn_meso_null))
      bup(ranef(rn_meso_null), stat="mean") else rep(NA, nSites)
    Nhat_meso_rn_full <- if (!is.null(rn_meso_full))
      bup(ranef(rn_meso_full), stat="mean") else rep(NA, nSites)
    
    results$rn_blup_null <- fit_abund_sem(Nhat_pred_rn_null, Nhat_prey_rn_null,
                                          "RN null BLUP SEM", N_hat_meso=Nhat_meso_rn_null)
    results$rn_blup_full <- fit_abund_sem(Nhat_pred_rn_full, Nhat_prey_rn_full,
                                          "RN full BLUP SEM", N_hat_meso=Nhat_meso_rn_full)
  } else {
    results$rn_blup_null <- fit_abund_sem(Nhat_pred_rn_null, Nhat_prey_rn_null,
                                          "RN null BLUP SEM")
    results$rn_blup_full <- fit_abund_sem(Nhat_pred_rn_full, Nhat_prey_rn_full,
                                          "RN full BLUP SEM")
  }
  
  # ============================================================
  #  MODELS 4–5: ABUNDANCE NIMBLE
  #  Landscape toggle: with_land uses b_land/sigma_land;
  #  without uses no landscape terms.
  # ============================================================
  message("\n--- Models 4–5: Abundance NIMBLE ---")
  
  nim_const_abund <- if (use_landscape)
    list(nSites=nSites, nVisits=nVisits, n_landscapes=n_landscapes,
         landscape=landscape, Env1=Env1, Env2=Env2) else
           list(nSites=nSites, nVisits=nVisits, Env1=Env1, Env2=Env2)
  
  # ── Model 4: Co-abundance N-mixture ───────────────────────────────────────
  if (meso) {
    if (use_landscape) {
      coabund_code_meso <- nimbleCode({
        a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)
        c_int ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); p_meso ~ dbeta(1,1)
        b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
        b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); p_prey ~ dbeta(1,1)
        sigma_land_pred ~ dexp(1); sigma_land_meso ~ dexp(1); sigma_land_prey ~ dexp(1)
        for (l in 1:n_landscapes) {
          b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
          b_land_meso[l] ~ dnorm(0, sd=sigma_land_meso)
          b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
        }
        for (i in 1:nSites) {
          log(lambda_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
          N_pred[i] ~ dpois(lambda_pred[i] * z_land_pred[i])
          for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }
          log(lambda_meso[i]) <- c_int + c_pred*N_pred[i] + b_land_meso[landscape[i]]
          N_meso[i] ~ dpois(lambda_meso[i])
          for (j in 1:nVisits) { C_meso[i,j] ~ dbin(p_meso, N_meso[i]) }
          log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i] +
            b_meso*N_meso[i] + b_land_prey[landscape[i]]
          N_prey[i] ~ dpois(lambda_prey[i] * z_land_prey[i])
          for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
        }
        indirect_Env1_via_pred <- a_env * b_pred
        indirect_Env1_via_meso <- a_env * c_pred * b_meso
        indirect_pred_via_meso <- c_pred * b_meso
        total_pred_on_prey     <- b_pred + indirect_pred_via_meso
        total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
      })
      nd4 <- list(C_pred=C_pred, C_meso=C_meso, C_prey=C_prey,
                  z_land_pred=z_land_pred, z_land_prey=z_land_prey)
      ni4 <- list(a_int=0,a_env=0,p_pred=0.5,
                  c_int=0,c_pred=0,p_meso=0.5,
                  b_int=0,b_env=0,b_pred=0,b_meso=0,p_prey=0.5,
                  sigma_land_pred=0.3,sigma_land_meso=0.3,sigma_land_prey=0.3,
                  b_land_pred=rep(0,n_landscapes),
                  b_land_meso=rep(0,n_landscapes),
                  b_land_prey=rep(0,n_landscapes),
                  N_pred=ifelse(z_land_pred==1L,apply(C_pred,1,max)+1L,0L),
                  N_meso=apply(C_meso,1,max)+1L,
                  N_prey=ifelse(z_land_prey==1L,apply(C_prey,1,max)+1L,0L))
      params4 <- c("a_int","a_env","p_pred","c_int","c_pred","p_meso",
                   "b_int","b_env","b_pred","b_meso","p_prey",
                   "sigma_land_pred","sigma_land_meso","sigma_land_prey",
                   "indirect_Env1_via_pred","indirect_Env1_via_meso",
                   "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
    } else {
      coabund_code_meso <- nimbleCode({
        a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)
        c_int ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); p_meso ~ dbeta(1,1)
        b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
        b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); p_prey ~ dbeta(1,1)
        for (i in 1:nSites) {
          log(lambda_pred[i]) <- a_int + a_env*Env1[i]
          N_pred[i] ~ dpois(lambda_pred[i])
          for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }
          log(lambda_meso[i]) <- c_int + c_pred*N_pred[i]
          N_meso[i] ~ dpois(lambda_meso[i])
          for (j in 1:nVisits) { C_meso[i,j] ~ dbin(p_meso, N_meso[i]) }
          log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i] +
            b_meso*N_meso[i]
          N_prey[i] ~ dpois(lambda_prey[i])
          for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
        }
        indirect_Env1_via_pred <- a_env * b_pred
        indirect_Env1_via_meso <- a_env * c_pred * b_meso
        indirect_pred_via_meso <- c_pred * b_meso
        total_pred_on_prey     <- b_pred + indirect_pred_via_meso
        total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
      })
      nd4 <- list(C_pred=C_pred, C_meso=C_meso, C_prey=C_prey)
      ni4 <- list(a_int=0,a_env=0,p_pred=0.5,
                  c_int=0,c_pred=0,p_meso=0.5,
                  b_int=0,b_env=0,b_pred=0,b_meso=0,p_prey=0.5,
                  N_pred=apply(C_pred,1,max)+1L,
                  N_meso=apply(C_meso,1,max)+1L,
                  N_prey=apply(C_prey,1,max)+1L)
      params4 <- c("a_int","a_env","p_pred","c_int","c_pred","p_meso",
                   "b_int","b_env","b_pred","b_meso","p_prey",
                   "indirect_Env1_via_pred","indirect_Env1_via_meso",
                   "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
    }
  } else {
    
   if (use_landscape) {
    coabund_code <- nimbleCode({
      a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)  # Beta prior
      b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2); b_pred ~ dnorm(0,sd=2)
      p_prey ~ dbeta(1,1)
      sigma_land_pred ~ dexp(1); sigma_land_prey ~ dexp(1)
      for (l in 1:n_landscapes) {
        b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
        b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
      }
      for (i in 1:nSites) {
        log(lambda_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
        N_pred[i] ~ dpois(lambda_pred[i] * z_land_pred[i])         # z_land_pred as data
        for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }  # direct p
        log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i] +
          b_land_prey[landscape[i]]
        N_prey[i] ~ dpois(lambda_prey[i] * z_land_prey[i])
        for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
      }
      indirect_Env1 <- a_env * b_pred
    })
    nd4 <- list(C_pred=C_pred, C_prey=C_prey,
                z_land_pred=z_land_pred, z_land_prey=z_land_prey)   # map z_land → z_land_pred
    ni4 <- list(a_int=0, a_env=0, p_pred=0.5,            # p not lp
                b_int=0, b_env=0, b_pred=0, p_prey=0.5,
                sigma_land_pred=0.3, sigma_land_prey=0.3,
                b_land_pred=rep(0,n_landscapes), b_land_prey=rep(0,n_landscapes),
                N_pred=ifelse(z_land_pred==1L, apply(C_pred,1,max)+1L, 0L),
                N_prey=ifelse(z_land_prey==1L, apply(C_prey,1,max)+1L, 0L))
    params4 <- c("a_int","a_env","p_pred","b_int","b_env","b_pred","p_prey",
                 "sigma_land_pred","sigma_land_prey","indirect_Env1")
  } else {
    coabund_code <- nimbleCode({
      a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)
      b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2); b_pred ~ dnorm(0,sd=2)
      p_prey ~ dbeta(1,1)
      for (i in 1:nSites) {
        log(lambda_pred[i]) <- a_int + a_env*Env1[i]
        N_pred[i] ~ dpois(lambda_pred[i])
        for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }
        log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i]
        N_prey[i] ~ dpois(lambda_prey[i])
        for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
      }
      indirect_Env1 <- a_env * b_pred
    })
    nd4 <- list(C_pred=C_pred, C_prey=C_prey)
    ni4 <- list(a_int=0, a_env=0, p_pred=0.5,
                b_int=0, b_env=0, b_pred=0, p_prey=0.5,
                N_pred=apply(C_pred,1,max)+1L,
                N_prey=apply(C_prey,1,max)+1L)
    params4 <- c("a_int","a_env","p_pred","b_int","b_env","b_pred",
                 "p_prey","indirect_Env1")
  } 
    }
  
    # These run for BOTH meso and non-meso:
    m4 <- nimbleModel(if (meso) coabund_code_meso else coabund_code,
                      nim_const_abund, nd4, ni4)
    cm4    <- compileNimble(m4)
    conf4  <- configureMCMC(cm4, monitors=params4)
    mcmc4  <- buildMCMC(conf4)
    cmcmc4 <- compileNimble(mcmc4, project=m4)
    samps4 <- runMCMC(cmcmc4, niter=nimble_iter, nburnin=nimble_burnin,
                      nchains=nimble_chains, thin=nimble_thin,
                      samplesAsCodaMCMC=TRUE)
    results$coabund <- if (setting=="LOCAL") {
      list(samples=samps4, summary=post_summary(samps4,params4))
    } else list(summary=post_summary(samps4,params4))
  

  
  
  # ── Model 5: Integrated N-mixture SEM ────────────────────────────────────
    message("\n--- Model 4c: Integrated N-mixture SEM (NIMBLE) ---")
    
    if (meso) {
      if (use_landscape) {
        sem_code_meso <- nimbleCode({
          a_int  ~ dnorm(0,sd=2); a_env  ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)
          c_int  ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); p_meso ~ dbeta(1,1)
          b_int  ~ dnorm(0,sd=2); b_env  ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); p_prey ~ dbeta(1,1)
          sigma_land_pred ~ dexp(1); sigma_land_meso ~ dexp(1); sigma_land_prey ~ dexp(1)
          for (l in 1:n_landscapes) {
            b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
            b_land_meso[l] ~ dnorm(0, sd=sigma_land_meso)
            b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
          }
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
            N_pred[i] ~ dpois(lambda_pred[i] * z_land_pred[i])
            for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }
            log(lambda_meso[i]) <- c_int + c_pred*N_pred[i] + b_land_meso[landscape[i]]
            N_meso[i] ~ dpois(lambda_meso[i])
            for (j in 1:nVisits) { C_meso[i,j] ~ dbin(p_meso, N_meso[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] +
              b_pred*N_pred[i] + b_meso*N_meso[i] +
              b_land_prey[landscape[i]]
            N_prey[i] ~ dpois(lambda_prey[i] * z_land_prey[i])
            for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
          }
          indirect_Env1_via_pred <- a_env * b_pred
          indirect_Env1_via_meso <- a_env * c_pred * b_meso
          indirect_pred_via_meso <- c_pred * b_meso
          total_pred_on_prey     <- b_pred + indirect_pred_via_meso
          total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
        })
        nd5 <- list(C_pred=C_pred, C_meso=C_meso, C_prey=C_prey,
                    z_land_pred=z_land_pred, z_land_prey=z_land_prey)
        ni5 <- list(a_int=0,a_env=0,p_pred=0.5,
                    c_int=0,c_pred=0,p_meso=0.5,
                    b_int=0,b_env=0,b_pred=0,b_meso=0,p_prey=0.5,
                    sigma_land_pred=0.3,sigma_land_meso=0.3,sigma_land_prey=0.3,
                    b_land_pred=rep(0,n_landscapes),
                    b_land_meso=rep(0,n_landscapes),
                    b_land_prey=rep(0,n_landscapes),
                    N_pred=ifelse(z_land_pred==1L,apply(C_pred,1,max)+1L,0L),
                    N_meso=apply(C_meso,1,max)+1L,
                    N_prey=ifelse(z_land_prey==1L,apply(C_prey,1,max)+1L,0L))
        params5 <- c("a_int","a_env","p_pred","c_int","c_pred","p_meso",
                     "b_int","b_env","b_pred","b_meso","p_prey",
                     "sigma_land_pred","sigma_land_meso","sigma_land_prey",
                     "indirect_Env1_via_pred","indirect_Env1_via_meso",
                     "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
        
      } else {  # meso, no landscape
        sem_code_meso <- nimbleCode({
          a_int  ~ dnorm(0,sd=2); a_env  ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)
          c_int  ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); p_meso ~ dbeta(1,1)
          b_int  ~ dnorm(0,sd=2); b_env  ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); p_prey ~ dbeta(1,1)
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i]
            N_pred[i] ~ dpois(lambda_pred[i])
            for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }
            log(lambda_meso[i]) <- c_int + c_pred*N_pred[i]
            N_meso[i] ~ dpois(lambda_meso[i])
            for (j in 1:nVisits) { C_meso[i,j] ~ dbin(p_meso, N_meso[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] +
              b_pred*N_pred[i] + b_meso*N_meso[i]
            N_prey[i] ~ dpois(lambda_prey[i])
            for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
          }
          indirect_Env1_via_pred <- a_env * b_pred
          indirect_Env1_via_meso <- a_env * c_pred * b_meso
          indirect_pred_via_meso <- c_pred * b_meso
          total_pred_on_prey     <- b_pred + indirect_pred_via_meso
          total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
        })
        nd5 <- list(C_pred=C_pred, C_meso=C_meso, C_prey=C_prey)
        ni5 <- list(a_int=0,a_env=0,p_pred=0.5,
                    c_int=0,c_pred=0,p_meso=0.5,
                    b_int=0,b_env=0,b_pred=0,b_meso=0,p_prey=0.5,
                    N_pred=apply(C_pred,1,max)+1L,
                    N_meso=apply(C_meso,1,max)+1L,
                    N_prey=apply(C_prey,1,max)+1L)
        params5 <- c("a_int","a_env","p_pred","c_int","c_pred","p_meso",
                     "b_int","b_env","b_pred","b_meso","p_prey",
                     "indirect_Env1_via_pred","indirect_Env1_via_meso",
                     "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
      }
      
    } else {  # non-meso
      
      if (use_landscape) {
        sem_code <- nimbleCode({
          a_int  ~ dnorm(0,sd=2); a_env  ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)
          b_int  ~ dnorm(0,sd=2); b_env  ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); p_prey ~ dbeta(1,1)
          sigma_land_pred ~ dexp(1); sigma_land_prey ~ dexp(1)
          for (l in 1:n_landscapes) {
            b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
            b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
          }
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
            N_pred[i] ~ dpois(lambda_pred[i] * z_land_pred[i])
            for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i] +
              b_land_prey[landscape[i]]
            N_prey[i] ~ dpois(lambda_prey[i] * z_land_prey[i])
            for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
          }
          indirect_Env1_on_prey <- a_env * b_pred
          total_Env1_on_prey    <- indirect_Env1_on_prey
        })
        nd5 <- list(C_pred=C_pred, C_prey=C_prey,
                    z_land_pred=z_land_pred, z_land_prey=z_land_prey)
        ni5 <- list(a_int=0,a_env=0,p_pred=0.5,
                    b_int=0,b_env=0,b_pred=0,p_prey=0.5,
                    sigma_land_pred=0.3,sigma_land_prey=0.3,
                    b_land_pred=rep(0,n_landscapes),b_land_prey=rep(0,n_landscapes),
                    N_pred=ifelse(z_land_pred==1L,apply(C_pred,1,max)+1L,0L),
                    N_prey=ifelse(z_land_prey==1L,apply(C_prey,1,max)+1L,0L))
        params5 <- c("a_int","a_env","p_pred","b_int","b_env","b_pred","p_prey",
                     "sigma_land_pred","sigma_land_prey",
                     "indirect_Env1_on_prey","total_Env1_on_prey")
        
      } else {  # non-meso, no landscape
        sem_code <- nimbleCode({
          a_int  ~ dnorm(0,sd=2); a_env  ~ dnorm(0,sd=2); p_pred ~ dbeta(1,1)
          b_int  ~ dnorm(0,sd=2); b_env  ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); p_prey ~ dbeta(1,1)
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i]
            N_pred[i] ~ dpois(lambda_pred[i])
            for (j in 1:nVisits) { C_pred[i,j] ~ dbin(p_pred, N_pred[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i]
            N_prey[i] ~ dpois(lambda_prey[i])
            for (j in 1:nVisits) { C_prey[i,j] ~ dbin(p_prey, N_prey[i]) }
          }
          indirect_Env1_on_prey <- a_env * b_pred
          total_Env1_on_prey    <- indirect_Env1_on_prey
        })
        nd5 <- list(C_pred=C_pred, C_prey=C_prey)
        ni5 <- list(a_int=0,a_env=0,p_pred=0.5,
                    b_int=0,b_env=0,b_pred=0,p_prey=0.5,
                    N_pred=apply(C_pred,1,max)+1L,
                    N_prey=apply(C_prey,1,max)+1L)
        params5 <- c("a_int","a_env","p_pred","b_int","b_env","b_pred","p_prey",
                     "indirect_Env1_on_prey","total_Env1_on_prey")
      }
    }
    
    m5 <- nimbleModel(if (meso) sem_code_meso else sem_code,
                      nim_const_abund, nd5, ni5)
    cm5    <- compileNimble(m5)
    conf5  <- configureMCMC(cm5, monitors=params5)
    mcmc5  <- buildMCMC(conf5)
    cmcmc5 <- compileNimble(mcmc5, project=m5)
    samps5 <- runMCMC(cmcmc5, niter=nimble_iter, nburnin=nimble_burnin,
                      nchains=nimble_chains, thin=nimble_thin,
                      samplesAsCodaMCMC=TRUE)
    results$sem_nimble <- if (setting=="LOCAL") {
      list(samples=samps5, summary=post_summary(samps5,params5))
    } else list(summary=post_summary(samps5,params5))
  
    
    # ============================================================
    #  MODEL 8: INTEGRATED ROYLE-NICHOLS SEM (NIMBLE)
    #  Same abundance process as Model 5. Detection via
    #  Y[i,j] ~ dbern(1 - (1-r)^N[i]) on binary data.
    #  b_pred on log scale — same estimand as Model 5.
    # ============================================================
    message("\n--- Model 8: Integrated Royle-Nichols SEM (NIMBLE) ---")
    
    if (meso) {
      if (use_landscape) {
        rn_code_meso <- nimbleCode({
          a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); r_pred ~ dbeta(1,1)
          c_int ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); r_meso ~ dbeta(1,1)
          b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); r_prey ~ dbeta(1,1)
          sigma_land_pred ~ dexp(1); sigma_land_meso ~ dexp(1); sigma_land_prey ~ dexp(1)
          for (l in 1:n_landscapes) {
            b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
            b_land_meso[l] ~ dnorm(0, sd=sigma_land_meso)
            b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
          }
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
            N_pred[i] ~ dpois(lambda_pred[i] * z_land_pred[i])
            for (j in 1:nVisits) { Y_pred[i,j] ~ dbern(1 - (1-r_pred)^N_pred[i]) }
            log(lambda_meso[i]) <- c_int + c_pred*N_pred[i] + b_land_meso[landscape[i]]
            N_meso[i] ~ dpois(lambda_meso[i])
            for (j in 1:nVisits) { Y_meso[i,j] ~ dbern(1 - (1-r_meso)^N_meso[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] +
              b_pred*N_pred[i] + b_meso*N_meso[i] +
              b_land_prey[landscape[i]]
            N_prey[i] ~ dpois(lambda_prey[i] * z_land_prey[i])
            for (j in 1:nVisits) { Y_prey[i,j] ~ dbern(1 - (1-r_prey)^N_prey[i]) }
          }
          indirect_Env1_via_pred <- a_env * b_pred
          indirect_Env1_via_meso <- a_env * c_pred * b_meso
          indirect_pred_via_meso <- c_pred * b_meso
          total_pred_on_prey     <- b_pred + indirect_pred_via_meso
          total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
        })
        nd8 <- list(Y_pred=Y_pred, Y_meso=Y_meso, Y_prey=Y_prey,
                    z_land_pred=z_land_pred, z_land_prey=z_land_prey)
        ni8 <- list(a_int=0, a_env=0, r_pred=0.5,
                    c_int=0, c_pred=0, r_meso=0.5,
                    b_int=0, b_env=0, b_pred=0, b_meso=0, r_prey=0.5,
                    sigma_land_pred=0.3, sigma_land_meso=0.3, sigma_land_prey=0.3,
                    b_land_pred=rep(0,n_landscapes),
                    b_land_meso=rep(0,n_landscapes),
                    b_land_prey=rep(0,n_landscapes),
                    N_pred=ifelse(z_land_pred==1L, apply(Y_pred,1,max)+1L, 0L),
                    N_meso=apply(Y_meso,1,max)+1L,
                    N_prey=ifelse(z_land_prey==1L, apply(Y_prey,1,max)+1L, 0L))
        params8 <- c("a_int","a_env","r_pred","c_int","c_pred","r_meso",
                     "b_int","b_env","b_pred","b_meso","r_prey",
                     "sigma_land_pred","sigma_land_meso","sigma_land_prey",
                     "indirect_Env1_via_pred","indirect_Env1_via_meso",
                     "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
        
      } else {  # meso, no landscape
        rn_code_meso <- nimbleCode({
          a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); r_pred ~ dbeta(1,1)
          c_int ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); r_meso ~ dbeta(1,1)
          b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); r_prey ~ dbeta(1,1)
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i]
            N_pred[i] ~ dpois(lambda_pred[i])
            for (j in 1:nVisits) { Y_pred[i,j] ~ dbern(1 - (1-r_pred)^N_pred[i]) }
            log(lambda_meso[i]) <- c_int + c_pred*N_pred[i]
            N_meso[i] ~ dpois(lambda_meso[i])
            for (j in 1:nVisits) { Y_meso[i,j] ~ dbern(1 - (1-r_meso)^N_meso[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] +
              b_pred*N_pred[i] + b_meso*N_meso[i]
            N_prey[i] ~ dpois(lambda_prey[i])
            for (j in 1:nVisits) { Y_prey[i,j] ~ dbern(1 - (1-r_prey)^N_prey[i]) }
          }
          indirect_Env1_via_pred <- a_env * b_pred
          indirect_Env1_via_meso <- a_env * c_pred * b_meso
          indirect_pred_via_meso <- c_pred * b_meso
          total_pred_on_prey     <- b_pred + indirect_pred_via_meso
          total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
        })
        nd8 <- list(Y_pred=Y_pred, Y_meso=Y_meso, Y_prey=Y_prey)
        ni8 <- list(a_int=0, a_env=0, r_pred=0.5,
                    c_int=0, c_pred=0, r_meso=0.5,
                    b_int=0, b_env=0, b_pred=0, b_meso=0, r_prey=0.5,
                    N_pred=apply(Y_pred,1,max)+1L,
                    N_meso=apply(Y_meso,1,max)+1L,
                    N_prey=apply(Y_prey,1,max)+1L)
        params8 <- c("a_int","a_env","r_pred","c_int","c_pred","r_meso",
                     "b_int","b_env","b_pred","b_meso","r_prey",
                     "indirect_Env1_via_pred","indirect_Env1_via_meso",
                     "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
      }
      
    } else {  # non-meso
      
      if (use_landscape) {
        rn_code <- nimbleCode({
          a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); r_pred ~ dbeta(1,1)
          b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); r_prey ~ dbeta(1,1)
          sigma_land_pred ~ dexp(1); sigma_land_prey ~ dexp(1)
          for (l in 1:n_landscapes) {
            b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
            b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
          }
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
            N_pred[i] ~ dpois(lambda_pred[i] * z_land_pred[i])
            for (j in 1:nVisits) { Y_pred[i,j] ~ dbern(1 - (1-r_pred)^N_pred[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i] +
              b_land_prey[landscape[i]]
            N_prey[i] ~ dpois(lambda_prey[i] * z_land_prey[i])
            for (j in 1:nVisits) { Y_prey[i,j] ~ dbern(1 - (1-r_prey)^N_prey[i]) }
          }
          indirect_Env1_on_prey <- a_env * b_pred
          total_Env1_on_prey    <- indirect_Env1_on_prey
        })
        nd8 <- list(Y_pred=Y_pred, Y_prey=Y_prey,
                    z_land_pred=z_land_pred, z_land_prey=z_land_prey)
        ni8 <- list(a_int=0, a_env=0, r_pred=0.5,
                    b_int=0, b_env=0, b_pred=0, r_prey=0.5,
                    sigma_land_pred=0.3, sigma_land_prey=0.3,
                    b_land_pred=rep(0,n_landscapes), b_land_prey=rep(0,n_landscapes),
                    N_pred=ifelse(z_land_pred==1L, apply(Y_pred,1,max)+1L, 0L),
                    N_prey=ifelse(z_land_prey==1L, apply(Y_prey,1,max)+1L, 0L))
        params8 <- c("a_int","a_env","r_pred","b_int","b_env","b_pred","r_prey",
                     "sigma_land_pred","sigma_land_prey",
                     "indirect_Env1_on_prey","total_Env1_on_prey")
        
      } else {  # non-meso, no landscape
        rn_code <- nimbleCode({
          a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); r_pred ~ dbeta(1,1)
          b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
          b_pred ~ dnorm(0,sd=2); r_prey ~ dbeta(1,1)
          for (i in 1:nSites) {
            log(lambda_pred[i]) <- a_int + a_env*Env1[i]
            N_pred[i] ~ dpois(lambda_pred[i])
            for (j in 1:nVisits) { Y_pred[i,j] ~ dbern(1 - (1-r_pred)^N_pred[i]) }
            log(lambda_prey[i]) <- b_int + b_env*Env2[i] + b_pred*N_pred[i]
            N_prey[i] ~ dpois(lambda_prey[i])
            for (j in 1:nVisits) { Y_prey[i,j] ~ dbern(1 - (1-r_prey)^N_prey[i]) }
          }
          indirect_Env1_on_prey <- a_env * b_pred
          total_Env1_on_prey    <- indirect_Env1_on_prey
        })
        nd8 <- list(Y_pred=Y_pred, Y_prey=Y_prey)
        ni8 <- list(a_int=0, a_env=0, r_pred=0.5,
                    b_int=0, b_env=0, b_pred=0, r_prey=0.5,
                    N_pred=apply(Y_pred,1,max)+1L,
                    N_prey=apply(Y_prey,1,max)+1L)
        params8 <- c("a_int","a_env","r_pred","b_int","b_env","b_pred","r_prey",
                     "indirect_Env1_on_prey","total_Env1_on_prey")
      }
    }
    
    m8 <- nimbleModel(if (meso) rn_code_meso else rn_code,
                      nim_const_abund, nd8, ni8)
    cm8    <- compileNimble(m8)
    conf8  <- configureMCMC(cm8, monitors=params8)
    mcmc8  <- buildMCMC(conf8)
    cmcmc8 <- compileNimble(mcmc8, project=m8)
    samps8 <- runMCMC(cmcmc8, niter=nimble_iter, nburnin=nimble_burnin,
                      nchains=nimble_chains, thin=nimble_thin,
                      samplesAsCodaMCMC=TRUE)
    results$rn_nimble <- if (setting=="LOCAL") {
      list(samples=samps8, summary=post_summary(samps8,params8))
    } else list(summary=post_summary(samps8,params8))
    
    

    # ========================================================================
    # ── Model 9: Naïve detection rate SEM ──────────────────────────────────
    # =========================================================================
    
    message("\n--- Model 6: Naïve detection rate SEM ---")
    
    naive_pred <- numeric(nSites)
    naive_pred[idx_pred] <- rowMeans(Y_pred[idx_pred,,drop=F])
    naive_prey <- numeric(nSites)
    naive_prey[idx_prey] <- rowMeans(Y_prey[idx_prey,,drop=F])
    
    if (meso) {
      naive_meso <- rowMeans(Y_meso)
      results$naive_sem <- fit_occ_sem(naive_pred, naive_prey, "naive SEM",
                                       psi_meso=naive_meso)
    } else {
      results$naive_sem <- fit_occ_sem(naive_pred, naive_prey, "naive SEM")
    }
    
    # ── Models 10–11: Occupancy model BLUPs / predicted psi SEM ──────────────
    message("\n--- Models 7–8: Occupancy model + SEM ---")
    
    occ_pred_null <- NULL; occ_prey_null <- NULL
    occ_pred_full <- NULL; occ_prey_full <- NULL
    occ_meso_null <- NULL; occ_meso_full <- NULL
    
    umf_pred_occ_null <- unmarkedFrameOccu(y=Y_pred[idx_pred,,drop=F],
                                           siteCovs=sc_pred)
    umf_prey_occ_null <- unmarkedFrameOccu(y=Y_prey[idx_prey,,drop=F],
                                           siteCovs=sc_prey)
    
    occ_pred_null <- suppressWarnings(tryCatch(
      occu(~1~1, data=umf_pred_occ_null), error=function(e) NULL))
    occ_prey_null <- suppressWarnings(tryCatch(
      occu(~1~1, data=umf_prey_occ_null), error=function(e) NULL))
    
    occ_full_form_pred <- if (use_landscape) ~1~Env1+landscape else ~1~Env1
    occ_full_form_prey <- if (use_landscape) ~1~Env2+landscape else ~1~Env2
    
    occ_pred_full <- suppressWarnings(tryCatch(
      occu(occ_full_form_pred, data=umf_pred_occ_null), error=function(e) NULL))
    occ_prey_full <- suppressWarnings(tryCatch(
      occu(occ_full_form_prey, data=umf_prey_occ_null), error=function(e) NULL))
    
    psi_pred_null <- numeric(nSites)
    psi_prey_null <- numeric(nSites)
    psi_pred_null[idx_pred] <- if (!is.null(occ_pred_null))
      bup(ranef(occ_pred_null), stat="mean") else rep(NA, length(idx_pred))
    psi_prey_null[idx_prey] <- if (!is.null(occ_prey_null))
      bup(ranef(occ_prey_null), stat="mean") else rep(NA, length(idx_prey))
    
    psi_pred_full <- numeric(nSites)
    psi_prey_full <- numeric(nSites)
    psi_pred_full[idx_pred] <- if (!is.null(occ_pred_full))
      bup(ranef(occ_pred_full), stat="mean") else rep(NA, length(idx_pred))
    psi_prey_full[idx_prey] <- if (!is.null(occ_prey_full))
      bup(ranef(occ_prey_full), stat="mean") else rep(NA, length(idx_prey))
    
    if (meso) {
      sc_meso_occ <- if (use_landscape) data.frame(Env1=Env1, landscape=land_f) else
        data.frame(Env1=Env1)
      occ_full_form_meso <- if (use_landscape) ~1~Env1+landscape else ~1~Env1
      umf_meso_occ <- unmarkedFrameOccu(y=Y_meso, siteCovs=sc_meso_occ)
      occ_meso_null <- suppressWarnings(tryCatch(
        occu(~1~1, data=umf_meso_occ), error=function(e) NULL))
      occ_meso_full <- suppressWarnings(tryCatch(
        occu(occ_full_form_meso, data=umf_meso_occ), error=function(e) NULL))
      psi_meso_null <- if (!is.null(occ_meso_null))
        bup(ranef(occ_meso_null), stat="mean") else rep(NA, nSites)
      psi_meso_full <- if (!is.null(occ_meso_full))
        bup(ranef(occ_meso_full), stat="mean") else rep(NA, nSites)
      
      results$null_occ_sem <- fit_occ_sem(psi_pred_null, psi_prey_null,
                                          "null occ SEM", psi_meso=psi_meso_null)
      results$full_occ_sem <- fit_occ_sem(psi_pred_full, psi_prey_full,
                                          "full occ SEM", psi_meso=psi_meso_full)
    } else {
      results$null_occ_sem <- fit_occ_sem(psi_pred_null, psi_prey_null,
                                          "null occ SEM")
      results$full_occ_sem <- fit_occ_sem(psi_pred_full, psi_prey_full,
                                          "full occ SEM")
    }
    
  # ============================================================
  #  MODEL 12: INTEGRATED BAYESIAN OCCUPANCY SEM
  #  Applied to binarised Y_ij. SIV via z_land_pred (binary latent).
  #  b_pred on LOGIT scale.
  # ============================================================
  message("\n--- Model 9: Integrated Bayesian occupancy SEM ---")
  
  nim_const_occ <- if (use_landscape)
    list(nSites=nSites, nVisits=nVisits, n_landscapes=n_landscapes,
         landscape=landscape, Env1=Env1, Env2=Env2) else
           list(nSites=nSites, nVisits=nVisits, Env1=Env1, Env2=Env2)
  
  z_land_pred_init <- ifelse(z_land_pred==1L, apply(Y_pred,1,max), 0L)
  z_land_prey_init <- ifelse(z_land_prey==1L, apply(Y_prey,1,max), 0L)
  if (meso) z_meso_init <- apply(Y_meso,1,max)
  
  if (!meso) {
    if (use_landscape) {
      occ9_code <- nimbleCode({
        a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); lp_pred ~ dnorm(0,sd=2)
        b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
        b_pred ~ dnorm(0,sd=2); lp_prey ~ dnorm(0,sd=2)
        sigma_land_pred ~ dexp(1); sigma_land_prey ~ dexp(1)
        for (l in 1:n_landscapes) {
          b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
          b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
        }
        for (i in 1:nSites) {
          logit(psi_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
          z_pred[i] ~ dbern(psi_pred[i] * z_land_pred[i])
          for (j in 1:nVisits) {
            Y_pred[i,j] ~ dbern(z_pred[i] * ilogit(lp_pred))
          }
          logit(psi_prey[i]) <- b_int + b_env*Env2[i] + b_pred*z_pred[i] +
            b_land_prey[landscape[i]]
          z_prey[i] ~ dbern(psi_prey[i] * z_land_prey[i])
          for (j in 1:nVisits) {
            Y_prey[i,j] ~ dbern(z_prey[i] * ilogit(lp_prey))
          }
        }                                    # closes for (i in 1:nSites)
        indirect_Env1 <- a_env * b_pred
        total_Env1    <- indirect_Env1
      })
      nd9 <- list(Y_pred=Y_pred, Y_prey=Y_prey,
                  z_land_pred=z_land_pred, z_land_prey=z_land_prey)
      ni9 <- list(a_int=0,a_env=0,lp_pred=0,b_int=0,b_env=0,b_pred=0,lp_prey=0,
                  sigma_land_pred=0.3, sigma_land_prey=0.3,
                  b_land_pred=rep(0,n_landscapes), b_land_prey=rep(0,n_landscapes),
                  z_pred=z_land_pred_init, z_prey=z_land_prey_init)
      params9 <- c("a_int","a_env","lp_pred","b_int","b_env","b_pred","lp_prey",
                   "sigma_land_pred","sigma_land_prey","indirect_Env1","total_Env1")
    } else {
      occ9_code <- nimbleCode({
        a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); lp_pred ~ dnorm(0,sd=2)
        b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
        b_pred ~ dnorm(0,sd=2); lp_prey ~ dnorm(0,sd=2)
        for (i in 1:nSites) {
          logit(psi_pred[i]) <- a_int + a_env*Env1[i]
          z_pred[i] ~ dbern(psi_pred[i])
          for (j in 1:nVisits) { Y_pred[i,j] ~ dbern(z_pred[i] * ilogit(lp_pred)) }
          logit(psi_prey[i]) <- b_int + b_env*Env2[i] + b_pred*z_pred[i]
          z_prey[i] ~ dbern(psi_prey[i])
          for (j in 1:nVisits) { Y_prey[i,j] ~ dbern(z_prey[i] * ilogit(lp_prey)) }
        }
        indirect_Env1 <- a_env * b_pred
        total_Env1    <- indirect_Env1
      })
      nd9 <- list(Y_pred=Y_pred, Y_prey=Y_prey)
      ni9 <- list(a_int=0,a_env=0,lp_pred=0,b_int=0,b_env=0,b_pred=0,lp_prey=0,
                  z_pred=z_land_pred_init, z_prey=z_land_prey_init)
      params9 <- c("a_int","a_env","lp_pred","b_int","b_env","b_pred","lp_prey",
                   "indirect_Env1","total_Env1")
    }
  } else {  # meso
    if (use_landscape) {
      occ9_code <- nimbleCode({
        a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); lp_pred ~ dnorm(0,sd=2)
        c_int ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); lp_meso ~ dnorm(0,sd=2)
        b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
        b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); lp_prey ~ dnorm(0,sd=2)
        sigma_land_pred ~ dexp(1); sigma_land_meso ~ dexp(1); sigma_land_prey ~ dexp(1)
        for (l in 1:n_landscapes) {
          b_land_pred[l] ~ dnorm(0, sd=sigma_land_pred)
          b_land_meso[l] ~ dnorm(0, sd=sigma_land_meso)
          b_land_prey[l] ~ dnorm(0, sd=sigma_land_prey)
        }
        for (i in 1:nSites) {
          logit(psi_pred[i]) <- a_int + a_env*Env1[i] + b_land_pred[landscape[i]]
          z_pred[i] ~ dbern(psi_pred[i] * z_land_pred[i])
          for (j in 1:nVisits) {
            Y_pred[i,j] ~ dbern(z_pred[i] * ilogit(lp_pred))
          }
          logit(psi_meso[i]) <- c_int + c_pred*z_pred[i] + b_land_meso[landscape[i]]
          z_meso[i] ~ dbern(psi_meso[i])
          for (j in 1:nVisits) {
            Y_meso[i,j] ~ dbern(z_meso[i] * ilogit(lp_meso)) }
          logit(psi_prey[i]) <- b_int + b_env*Env2[i] + b_pred*z_pred[i] +
            b_meso*z_meso[i] + b_land_prey[landscape[i]]
          z_prey[i] ~ dbern(psi_prey[i] * z_land_prey[i])
          for (j in 1:nVisits) {
            Y_prey[i,j] ~ dbern(z_prey[i] * ilogit(lp_prey))
          }
        }
        indirect_Env1_via_pred <- a_env * b_pred
        indirect_Env1_via_meso <- a_env * c_pred * b_meso
        indirect_pred_via_meso <- c_pred * b_meso
        total_pred_on_prey     <- b_pred + indirect_pred_via_meso
        total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
      })
      nd9 <- list(Y_pred=Y_pred, Y_meso=Y_meso, Y_prey=Y_prey,
                  z_land_pred=z_land_pred, z_land_prey=z_land_prey)
      ni9 <- list(a_int=0,a_env=0,lp_pred=0,c_int=0,c_pred=0,lp_meso=0,
                  b_int=0,b_env=0,b_pred=0,b_meso=0,lp_prey=0,
                  sigma_land_pred=0.3,sigma_land_meso=0.3,sigma_land_prey=0.3,
                  b_land_pred=rep(0,n_landscapes),
                  b_land_meso=rep(0,n_landscapes),
                  b_land_prey=rep(0,n_landscapes),
                  z_pred=z_land_pred_init, z_meso=z_meso_init, z_prey=z_land_prey_init)
      params9 <- c("a_int","a_env","lp_pred","c_int","c_pred","lp_meso",
                   "b_int","b_env","b_pred","b_meso","lp_prey",
                   "sigma_land_pred","sigma_land_meso","sigma_land_prey",
                   "indirect_Env1_via_pred","indirect_Env1_via_meso",
                   "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
    } else {
      occ9_code <- nimbleCode({
        a_int ~ dnorm(0,sd=2); a_env ~ dnorm(0,sd=2); lp_pred ~ dnorm(0,sd=2)
        c_int ~ dnorm(0,sd=2); c_pred ~ dnorm(0,sd=2); lp_meso ~ dnorm(0,sd=2)
        b_int ~ dnorm(0,sd=2); b_env ~ dnorm(0,sd=2)
        b_pred ~ dnorm(0,sd=2); b_meso ~ dnorm(0,sd=2); lp_prey ~ dnorm(0,sd=2)
        for (i in 1:nSites) {
          logit(psi_pred[i]) <- a_int + a_env*Env1[i]
          z_pred[i] ~ dbern(psi_pred[i])
          for (j in 1:nVisits) { Y_pred[i,j] ~ dbern(z_pred[i] * ilogit(lp_pred)) }
          logit(psi_meso[i]) <- c_int + c_pred*z_pred[i]
          z_meso[i] ~ dbern(psi_meso[i])
          for (j in 1:nVisits) { Y_meso[i,j] ~ dbern(z_meso[i] * ilogit(lp_meso)) }
          logit(psi_prey[i]) <- b_int + b_env*Env2[i] +
            b_pred*z_pred[i] + b_meso*z_meso[i]
          z_prey[i] ~ dbern(psi_prey[i])
          for (j in 1:nVisits) { Y_prey[i,j] ~ dbern(z_prey[i] * ilogit(lp_prey)) }
        }
        indirect_Env1_via_pred <- a_env * b_pred
        indirect_Env1_via_meso <- a_env * c_pred * b_meso
        indirect_pred_via_meso <- c_pred * b_meso
        total_pred_on_prey     <- b_pred + indirect_pred_via_meso
        total_Env1_on_prey     <- indirect_Env1_via_pred + indirect_Env1_via_meso
      })
      nd9 <- list(Y_pred=Y_pred, Y_meso=Y_meso, Y_prey=Y_prey)
      ni9 <- list(a_int=0,a_env=0,lp_pred=0,c_int=0,c_pred=0,lp_meso=0,
                  b_int=0,b_env=0,b_pred=0,b_meso=0,lp_prey=0,
                  z_pred=z_land_pred_init, z_meso=z_meso_init, z_prey=z_land_prey_init)
      params9 <- c("a_int","a_env","lp_pred","c_int","c_pred","lp_meso",
                   "b_int","b_env","b_pred","b_meso","lp_prey",
                   "indirect_Env1_via_pred","indirect_Env1_via_meso",
                   "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey")
    }
  }
  
  samps9 <- run_nimble(occ9_code, nim_const_occ, nd9, ni9, params9,
                       nimble_iter, nimble_burnin, nimble_chains, nimble_thin)
  results$occ_sem <- if (setting=="LOCAL")
    list(samples=samps9, summary=post_summary(samps9,params9)) else
      list(summary=post_summary(samps9,params9))
  
  
  # ============================================================
  #  COMPARISON TABLE
  # ============================================================
  message("\n--- Building comparison table ---")
  
  get4 <- function(key, nm) safe_get(results[[key]], nm)
  gpm4 <- function(key, nm) {
    tryCatch({
      samps <- results[[key]]$samples
      summ  <- results[[key]]$summary
      if (!is.null(samps)) {
        as.numeric(gpm(samps, nm))
      } else {
        val <- summ$mean[summ$param == nm][1]
        if (is.null(val) || length(val) == 0) NA_real_ else as.numeric(val)
      }
    }, error = function(e) NA_real_)
  }
  
  if (!meso) {
    
    comparison <- tryCatch({
      data.frame(
        estimand = c("a_int","a_env","p_pred",
                     "b_int","b_env","b_pred","p_prey",
                     "sigma_land_pred","sigma_land_prey",
                     "indirect_Env1","total_Env1"),
        truth = c(truth$a_int, truth$a_env, truth$p_pred,
                  truth$b_int, truth$b_env, truth$b_pred, truth$p_prey,
                  truth$sigma_land_pred, truth$sigma_land_prey,
                  truth$indirect_Env1_via_pred, truth$total_Env1_on_prey),
        
        # Abundance lmer (count scale; no detection estimated)
        psem_maxcount = c(
          get4("psem_maxcount","a_int"), get4("psem_maxcount","a_env"), NA,
          get4("psem_maxcount","b_int"), get4("psem_maxcount","b_env"),
          get4("psem_maxcount","b_pred"), NA, NA, NA,
          get4("psem_maxcount","indirect_Env1"), NA),
        
        psem_blup_null = c(
          get4("psem_blup_null","a_int"), get4("psem_blup_null","a_env"), NA,
          get4("psem_blup_null","b_int"), get4("psem_blup_null","b_env"),
          get4("psem_blup_null","b_pred"), NA, NA, NA,
          get4("psem_blup_null","indirect_Env1"), NA),
        
        psem_blup_full = c(
          get4("psem_blup_full","a_int"), get4("psem_blup_full","a_env"), NA,
          get4("psem_blup_full","b_int"), get4("psem_blup_full","b_env"),
          get4("psem_blup_full","b_pred"), NA, NA, NA,
          get4("psem_blup_full","indirect_Env1"), NA),
        
        rn_blup_null = c(
          get4("rn_blup_null","a_int"), get4("rn_blup_null","a_env"), NA,
          get4("rn_blup_null","b_int"), get4("rn_blup_null","b_env"),
          get4("rn_blup_null","b_pred"), NA, NA, NA,
          get4("rn_blup_null","indirect_Env1"), NA),
        
        rn_blup_full = c(
          get4("rn_blup_full","a_int"), get4("rn_blup_full","a_env"), NA,
          get4("rn_blup_full","b_int"), get4("rn_blup_full","b_env"),
          get4("rn_blup_full","b_pred"), NA, NA, NA,
          get4("rn_blup_full","indirect_Env1"), NA),
        
        rn_nmix = c(
          gpm4("rn_nimble","a_int"), gpm4("rn_nimble","a_env"),
          gpm4("rn_nimble","r_pred"),
          gpm4("rn_nimble","b_int"), gpm4("rn_nimble","b_env"),
          gpm4("rn_nimble","b_pred"), gpm4("rn_nimble","r_prey"),
          if(use_landscape) gpm4("rn_nimble","sigma_land_pred") else NA,
          if(use_landscape) gpm4("rn_nimble","sigma_land_prey") else NA,
          gpm4("rn_nimble","indirect_Env1_on_prey"),
          gpm4("rn_nimble","total_Env1_on_prey")),
        
        # Co-abundance NIMBLE (log scale)
        coabund = c(
          gpm4("coabund","a_int"), gpm4("coabund","a_env"),
          gpm4("coabund","p_pred"),
          gpm4("coabund","b_int"), gpm4("coabund","b_env"),
          gpm4("coabund","b_pred"), gpm4("coabund","p_prey"),
          if(use_landscape) gpm4("coabund","sigma_land_pred") else NA,
          if(use_landscape) gpm4("coabund","sigma_land_prey") else NA,
          gpm4("coabund","indirect_Env1"), NA),
        
        # Integrated N-mix SEM NIMBLE (log scale)
        sem_nmix = c(
          gpm4("sem_nimble","a_int"), gpm4("sem_nimble","a_env"),
          gpm4("sem_nimble","p_pred"),
          gpm4("sem_nimble","b_int"), gpm4("sem_nimble","b_env"),
          gpm4("sem_nimble","b_pred"), gpm4("sem_nimble","p_prey"),
          if(use_landscape) gpm4("sem_nimble","sigma_land_pred") else NA,
          if(use_landscape) gpm4("sem_nimble","sigma_land_prey") else NA,
          gpm4("sem_nimble","indirect_Env1_on_prey"),
          gpm4("sem_nimble","total_Env1_on_prey")),
        
        # Naïve SEM (probability scale)
        naive_sem = c(
          get4("naive_sem","a_int"), get4("naive_sem","a_env"), NA,
          get4("naive_sem","b_int"), get4("naive_sem","b_env"),
          get4("naive_sem","b_pred"), NA, NA, NA,
          get4("naive_sem","indirect_Env1"), NA),
        
        # Null occ SEM (probability scale; p_* from occ sub-model)
        null_occ_sem = c(
          get4("null_occ_sem","a_int"), get4("null_occ_sem","a_env"),
          if(!is.null(occ_pred_null)) as.numeric(plogis(coef(occ_pred_null)["p(Int)"])) else NA_real_,
          get4("null_occ_sem","b_int"), get4("null_occ_sem","b_env"),
          get4("null_occ_sem","b_pred"),
          if(!is.null(occ_prey_null)) as.numeric(plogis(coef(occ_prey_null)["p(Int)"])) else NA_real_,
          NA, NA,
          get4("null_occ_sem","indirect_Env1"), NA),
        
        # Full occ SEM (probability scale; p_* from occ sub-model)
        full_occ_sem = c(
          get4("full_occ_sem","a_int"), get4("full_occ_sem","a_env"),
          if(!is.null(occ_pred_full)) as.numeric(plogis(coef(occ_pred_full)["p(Int)"])) else NA_real_,
          get4("full_occ_sem","b_int"), get4("full_occ_sem","b_env"),
          get4("full_occ_sem","b_pred"),
          if(!is.null(occ_prey_full)) as.numeric(plogis(coef(occ_prey_full)["p(Int)"])) else NA_real_,
          NA, NA,
          get4("full_occ_sem","indirect_Env1"), NA),
        
        # Integrated occ SEM NIMBLE (logit scale; p_* converted to probability)
        occ_sem = c(
          gpm4("occ_sem","a_int"), gpm4("occ_sem","a_env"),
          as.numeric(plogis(gpm4("occ_sem","lp_pred"))),
          gpm4("occ_sem","b_int"), gpm4("occ_sem","b_env"),
          gpm4("occ_sem","b_pred"),
          as.numeric(plogis(gpm4("occ_sem","lp_prey"))),
          if(use_landscape) gpm4("occ_sem","sigma_land_pred") else NA,
          if(use_landscape) gpm4("occ_sem","sigma_land_prey") else NA,
          gpm4("occ_sem","indirect_Env1"),
          gpm4("occ_sem","total_Env1"))
      )
    }, error = function(e) {
      message("Comparison error: ", e$message)
      NULL
    })
    
  } else {  # meso
    
    comparison <- tryCatch({
      data.frame(
        estimand = c(
          "a_int","a_env","p_pred",
          "c_int","c_pred","p_meso",
          "b_int","b_env","b_pred","b_meso","p_prey",
          "sigma_land_pred","sigma_land_meso","sigma_land_prey",
          "indirect_Env1_via_pred","indirect_Env1_via_meso",
          "indirect_pred_via_meso","total_pred_on_prey","total_Env1_on_prey"),
        truth = c(
          truth$a_int, truth$a_env, truth$p_pred,
          truth$c_int, truth$c_pred, truth$p_meso,
          truth$b_int, truth$b_env, truth$b_pred, truth$b_meso, truth$p_prey,
          truth$sigma_land_pred, truth$sigma_land_meso, truth$sigma_land_prey,
          truth$indirect_Env1_via_pred, truth$indirect_Env1_via_meso,
          truth$indirect_pred_via_meso, truth$total_pred_on_prey,
          truth$total_Env1_on_prey),
        
        # Models 1–3: Abundance lmer (count scale; no detection)
        psem_maxcount = c(
          get4("psem_maxcount","a_int"),  get4("psem_maxcount","a_env"),  NA,
          get4("psem_maxcount","c_int"),  get4("psem_maxcount","c_pred"), NA,
          get4("psem_maxcount","b_int"),  get4("psem_maxcount","b_env"),
          get4("psem_maxcount","b_pred"), get4("psem_maxcount","b_meso"), NA,
          NA, NA, NA,
          get4("psem_maxcount","indirect_Env1"),
          get4("psem_maxcount","indirect_Env1_via_meso"),
          get4("psem_maxcount","indirect_pred_via_meso"),
          get4("psem_maxcount","total_pred_on_prey"),
          get4("psem_maxcount","total_Env1_on_prey")),
        
        psem_blup_null = c(
          get4("psem_blup_null","a_int"),  get4("psem_blup_null","a_env"),  NA,
          get4("psem_blup_null","c_int"),  get4("psem_blup_null","c_pred"), NA,
          get4("psem_blup_null","b_int"),  get4("psem_blup_null","b_env"),
          get4("psem_blup_null","b_pred"), get4("psem_blup_null","b_meso"), NA,
          NA, NA, NA,
          get4("psem_blup_null","indirect_Env1"),
          get4("psem_blup_null","indirect_Env1_via_meso"),
          get4("psem_blup_null","indirect_pred_via_meso"),
          get4("psem_blup_null","total_pred_on_prey"),
          get4("psem_blup_null","total_Env1_on_prey")),
        
        psem_blup_full = c(
          get4("psem_blup_full","a_int"),  get4("psem_blup_full","a_env"),  NA,
          get4("psem_blup_full","c_int"),  get4("psem_blup_full","c_pred"), NA,
          get4("psem_blup_full","b_int"),  get4("psem_blup_full","b_env"),
          get4("psem_blup_full","b_pred"), get4("psem_blup_full","b_meso"), NA,
          NA, NA, NA,
          get4("psem_blup_full","indirect_Env1"),
          get4("psem_blup_full","indirect_Env1_via_meso"),
          get4("psem_blup_full","indirect_pred_via_meso"),
          get4("psem_blup_full","total_pred_on_prey"),
          get4("psem_blup_full","total_Env1_on_prey")),
        
        rn_blup_null = c(
          get4("rn_blup_null","a_int"),  get4("rn_blup_null","a_env"),  NA,
          get4("rn_blup_null","c_int"),  get4("rn_blup_null","c_pred"), NA,
          get4("rn_blup_null","b_int"),  get4("rn_blup_null","b_env"),
          get4("rn_blup_null","b_pred"), get4("rn_blup_null","b_meso"), NA,
          NA, NA, NA,
          get4("rn_blup_null","indirect_Env1"),
          get4("rn_blup_null","indirect_Env1_via_meso"),
          get4("rn_blup_null","indirect_pred_via_meso"),
          get4("rn_blup_null","total_pred_on_prey"),
          get4("rn_blup_null","total_Env1_on_prey")),
        
        rn_blup_full = c(
          get4("rn_blup_full","a_int"),  get4("rn_blup_full","a_env"),  NA,
          get4("rn_blup_full","c_int"),  get4("rn_blup_full","c_pred"), NA,
          get4("rn_blup_full","b_int"),  get4("rn_blup_full","b_env"),
          get4("rn_blup_full","b_pred"), get4("rn_blup_full","b_meso"), NA,
          NA, NA, NA,
          get4("rn_blup_full","indirect_Env1"),
          get4("rn_blup_full","indirect_Env1_via_meso"),
          get4("rn_blup_full","indirect_pred_via_meso"),
          get4("rn_blup_full","total_pred_on_prey"),
          get4("rn_blup_full","total_Env1_on_prey")),
        
        rn_nmix = c(
          gpm4("rn_nimble","a_int"),  gpm4("rn_nimble","a_env"),  gpm4("rn_nimble","r_pred"),
          gpm4("rn_nimble","c_int"),  gpm4("rn_nimble","c_pred"), gpm4("rn_nimble","r_meso"),
          gpm4("rn_nimble","b_int"),  gpm4("rn_nimble","b_env"),
          gpm4("rn_nimble","b_pred"), gpm4("rn_nimble","b_meso"), gpm4("rn_nimble","r_prey"),
          if(use_landscape) gpm4("rn_nimble","sigma_land_pred") else NA,
          if(use_landscape) gpm4("rn_nimble","sigma_land_meso") else NA,
          if(use_landscape) gpm4("rn_nimble","sigma_land_prey") else NA,
          gpm4("rn_nimble","indirect_Env1_via_pred"),
          gpm4("rn_nimble","indirect_Env1_via_meso"),
          gpm4("rn_nimble","indirect_pred_via_meso"),
          gpm4("rn_nimble","total_pred_on_prey"),
          gpm4("rn_nimble","total_Env1_on_prey")),
        
        
        # Models 4–5: Abundance NIMBLE (p_* probability scale)
        coabund = c(
          gpm4("coabund","a_int"),  gpm4("coabund","a_env"),  gpm4("coabund","p_pred"),
          gpm4("coabund","c_int"),  gpm4("coabund","c_pred"), gpm4("coabund","p_meso"),
          gpm4("coabund","b_int"),  gpm4("coabund","b_env"),
          gpm4("coabund","b_pred"), gpm4("coabund","b_meso"), gpm4("coabund","p_prey"),
          if(use_landscape) gpm4("coabund","sigma_land_pred") else NA,
          if(use_landscape) gpm4("coabund","sigma_land_meso") else NA,
          if(use_landscape) gpm4("coabund","sigma_land_prey") else NA,
          gpm4("coabund","indirect_Env1_via_pred"),
          gpm4("coabund","indirect_Env1_via_meso"),
          gpm4("coabund","indirect_pred_via_meso"),
          gpm4("coabund","total_pred_on_prey"),
          gpm4("coabund","total_Env1_on_prey")),
        
        sem_nmix = c(
          gpm4("sem_nimble","a_int"),  gpm4("sem_nimble","a_env"),  gpm4("sem_nimble","p_pred"),
          gpm4("sem_nimble","c_int"),  gpm4("sem_nimble","c_pred"), gpm4("sem_nimble","p_meso"),
          gpm4("sem_nimble","b_int"),  gpm4("sem_nimble","b_env"),
          gpm4("sem_nimble","b_pred"), gpm4("sem_nimble","b_meso"), gpm4("sem_nimble","p_prey"),
          if(use_landscape) gpm4("sem_nimble","sigma_land_pred") else NA,
          if(use_landscape) gpm4("sem_nimble","sigma_land_meso") else NA,
          if(use_landscape) gpm4("sem_nimble","sigma_land_prey") else NA,
          gpm4("sem_nimble","indirect_Env1_via_pred"),
          gpm4("sem_nimble","indirect_Env1_via_meso"),
          gpm4("sem_nimble","indirect_pred_via_meso"),
          gpm4("sem_nimble","total_pred_on_prey"),
          gpm4("sem_nimble","total_Env1_on_prey")),
        
        # Models 6–8: Occupancy lmer (probability scale; no sigma_land)
        naive_sem = c(
          get4("naive_sem","a_int"),  get4("naive_sem","a_env"),  NA,
          get4("naive_sem","c_int"),  get4("naive_sem","c_pred"), NA,
          get4("naive_sem","b_int"),  get4("naive_sem","b_env"),
          get4("naive_sem","b_pred"), get4("naive_sem","b_meso"), NA,
          NA, NA, NA,
          get4("naive_sem","indirect_Env1"),
          get4("naive_sem","indirect_Env1_via_meso"),
          get4("naive_sem","indirect_pred_via_meso"),
          get4("naive_sem","total_pred_on_prey"),
          get4("naive_sem","total_Env1_on_prey")),
        
        null_occ_sem = c(
          get4("null_occ_sem","a_int"),  get4("null_occ_sem","a_env"),
          if(!is.null(occ_pred_null)) as.numeric(plogis(coef(occ_pred_null)["p(Int)"])) else NA_real_,
          get4("null_occ_sem","c_int"),  get4("null_occ_sem","c_pred"),
          if(!is.null(occ_meso_null)) as.numeric(plogis(coef(occ_meso_null)["p(Int)"])) else NA_real_,
          get4("null_occ_sem","b_int"),  get4("null_occ_sem","b_env"),
          get4("null_occ_sem","b_pred"), get4("null_occ_sem","b_meso"),
          if(!is.null(occ_prey_null)) as.numeric(plogis(coef(occ_prey_null)["p(Int)"])) else NA_real_,
          NA, NA, NA,
          get4("null_occ_sem","indirect_Env1"),
          get4("null_occ_sem","indirect_Env1_via_meso"),
          get4("null_occ_sem","indirect_pred_via_meso"),
          get4("null_occ_sem","total_pred_on_prey"),
          get4("null_occ_sem","total_Env1_on_prey")),
        
        full_occ_sem = c(
          get4("full_occ_sem","a_int"),  get4("full_occ_sem","a_env"),
          if(!is.null(occ_pred_full)) as.numeric(plogis(coef(occ_pred_full)["p(Int)"])) else NA_real_,
          get4("full_occ_sem","c_int"),  get4("full_occ_sem","c_pred"),
          if(!is.null(occ_meso_full)) as.numeric(plogis(coef(occ_meso_full)["p(Int)"])) else NA_real_,
          get4("full_occ_sem","b_int"),  get4("full_occ_sem","b_env"),
          get4("full_occ_sem","b_pred"), get4("full_occ_sem","b_meso"),
          if(!is.null(occ_prey_full)) as.numeric(plogis(coef(occ_prey_full)["p(Int)"])) else NA_real_,
          NA, NA, NA,
          get4("full_occ_sem","indirect_Env1"),
          get4("full_occ_sem","indirect_Env1_via_meso"),
          get4("full_occ_sem","indirect_pred_via_meso"),
          get4("full_occ_sem","total_pred_on_prey"),
          get4("full_occ_sem","total_Env1_on_prey")),
        
        # Model 9: Integrated occ SEM NIMBLE (logit scale; p_* → probability)
        occ_sem = c(
          gpm4("occ_sem","a_int"),  gpm4("occ_sem","a_env"),
          as.numeric(plogis(gpm4("occ_sem","lp_pred"))),
          gpm4("occ_sem","c_int"),  gpm4("occ_sem","c_pred"),
          as.numeric(plogis(gpm4("occ_sem","lp_meso"))),
          gpm4("occ_sem","b_int"),  gpm4("occ_sem","b_env"),
          gpm4("occ_sem","b_pred"), gpm4("occ_sem","b_meso"),
          as.numeric(plogis(gpm4("occ_sem","lp_prey"))),
          if(use_landscape) gpm4("occ_sem","sigma_land_pred") else NA,
          if(use_landscape) gpm4("occ_sem","sigma_land_meso") else NA,
          if(use_landscape) gpm4("occ_sem","sigma_land_prey") else NA,
          gpm4("occ_sem","indirect_Env1_via_pred"),
          gpm4("occ_sem","indirect_Env1_via_meso"),
          gpm4("occ_sem","indirect_pred_via_meso"),
          gpm4("occ_sem","total_pred_on_prey"),
          gpm4("occ_sem","total_Env1_on_prey"))
      )
    }, error = function(e) {
      message("Comparison error: ", e$message)
      NULL
    })
    
  }  # end if/else meso
  
  results$comparison <- comparison
  message("\nDone.")
  results
}

# ============================================================
#  RUN
# ============================================================

cat(sprintf("\nSimulating data (seed=%d, landscape=%s)...\n", seed, use_landscape))
sim_out <- simulate_predprey(
  nSites=nSites, nVisits=nVisits, seed=seed,
  meanLambda_pred=meanLambda_pred,
  meanLambda_prey=meanLambda_prey,
  meanLambda_meso=meanLambda_meso,
  meso=run_meso, c_pred=meso_cpred, b_meso=meso_bmeso,
  use_landscape=use_landscape)

cat("\n── Simulation summary ──────────────────────────────────────\n")
cat("nSites:", nSites, "| landscape:", use_landscape, "\n")
cat("Mean N_pred:", round(mean(sim_out$sim_data$N_pred),2),
    "| N_prey:", round(mean(sim_out$sim_data$N_prey),2), "\n")
cat("True b_pred (log scale):", sim_out$truth$b_pred, "\n")

fit <- fit_all_models(
  sim_out,
  nimble_iter=nimble_iter, nimble_burnin=nimble_burnin,
  nimble_chains=nimble_chains, nimble_thin=nimble_thin,
  meso=run_meso, use_landscape=use_landscape)

# ── Print (LOCAL only) ────────────────────────────────────────────────────────
if (setting == "LOCAL") {
  print(fit$comparison, digits=3, row.names=FALSE)
}

# ── Save (HPC only) ──────────────────────────────────────────────────────────
if (setting == "HPC") {
  land_tag <- if (use_landscape) "LT" else "LF"
  out_dir <- file.path("/scratch/user/uqsand24/simulations/results",
                       paste0("v4rn_", scenario, "_", land_tag))
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  out_file <- file.path(out_dir, paste0("rep_", rep_id, ".rds"))
  
  saveRDS(list(
    comparison      = fit$comparison,
    sem_nimble      = list(summary = fit$sem_nimble$summary),
    coabund         = list(summary = fit$coabund$summary),
    rn_nimble       = list(summary = fit$rn_nimble$summary),
    occ_sem         = list(summary = fit$occ_sem$summary),
    seed            = rep_id,
    slurm           = Sys.getenv("SLURM_ARRAY_TASK_ID"),
    scenario        = scenario,
    nSites          = nSites,
    meso            = run_meso,
    meanLambda_pred = meanLambda_pred,
    meanLambda_prey = meanLambda_prey,
    meanLambda_meso = if (run_meso) meanLambda_meso else NA,
    mean_N_pred     = round(mean(sim_out$sim_data$N_pred), 2),
    mean_N_prey     = round(mean(sim_out$sim_data$N_prey), 2),
    mean_N_meso     = if (run_meso) round(mean(sim_out$sim_data$N_meso), 2) else NA
  ), file = out_file)
  
  message("Saved rep ", rep_id, " → ", out_file)
}


# 
# # ── Simulation data summaries for supplementary table ────────────────────────
# 
# scenarios <- list(
#   list(nSites = 50,   meanLambda_pred = 0.5, meanLambda_meso = 1.5, 
#        meanLambda_prey = 2,  label = "Rare, n=50"),
#   list(nSites = 500,  meanLambda_pred = 0.5, meanLambda_meso = 1.5, 
#        meanLambda_prey = 2,  label = "Rare, n=500"),
#   list(nSites = 1000, meanLambda_pred = 0.5, meanLambda_meso = 1.5, 
#        meanLambda_prey = 2,  label = "Rare, n=1000"),
#   list(nSites = 50,   meanLambda_pred = 2,   meanLambda_meso = 4,   
#        meanLambda_prey = 10, label = "Abundant, n=50"),
#   list(nSites = 500,  meanLambda_pred = 2,   meanLambda_meso = 4,   
#        meanLambda_prey = 10, label = "Abundant, n=500"),
#   list(nSites = 1000, meanLambda_pred = 2,   meanLambda_meso = 4,   
#        meanLambda_prey = 10, label = "Abundant, n=1000")
# )
# 
# summarise_sim <- function(sim_out, label) {
#   
#   dat <- sim_out$sim_data
#   
#   # ── Count data summaries (C matrices) ──────────────────────────────────────
#   count_summary <- function(C, species) {
#     site_totals <- rowSums(C)
#     data.frame(
#       scenario       = label,
#       species        = species,
#       mean_count     = round(mean(site_totals), 2),
#       min_count      = min(site_totals),
#       max_count      = max(site_totals),
#       prop_zero_sites = round(mean(site_totals == 0), 3)
#     )
#   }
#   
#   # ── Detection history summaries (Y matrices) ────────────────────────────────
#   det_summary <- function(Y, species) {
#     naive_rate     <- rowMeans(Y)
#     never_detected <- mean(rowSums(Y) == 0)
#     mean_dets_when_present <- mean(rowSums(Y)[rowSums(Y) > 0])
#     data.frame(
#       scenario                  = label,
#       species                   = species,
#       mean_naive_det_rate       = round(mean(naive_rate), 3),
#       prop_sites_never_detected = round(never_detected, 3),
#       mean_dets_if_detected     = round(mean_dets_when_present, 2)
#     )
#   }
#   
#   counts <- bind_rows(
#     count_summary(dat$C_pred, "Predator"),
#     count_summary(dat$C_meso, "Mesopredator"),
#     count_summary(dat$C_prey, "Prey")
#   )
#   
#   dets <- bind_rows(
#     det_summary(dat$Y_pred, "Predator"),
#     det_summary(dat$Y_meso, "Mesopredator"),
#     det_summary(dat$Y_prey, "Prey")
#   )
#   
#   list(counts = counts, detections = dets)
# }
# 
# # ── Run across all scenarios using seed = 1 as representative ────────────────
# all_counts <- list()
# all_dets   <- list()
# 
# for (sc in scenarios) {
#   sim_out <- simulate_predprey(
#     nSites          = sc$nSites,
#     nVisits         = nVisits,
#     seed            = 1,
#     meanLambda_pred = sc$meanLambda_pred,
#     meanLambda_prey = sc$meanLambda_prey,
#     meanLambda_meso = sc$meanLambda_meso,
#     meso            = run_meso,
#     use_landscape   = use_landscape
#   )
#   
#   result         <- summarise_sim(sim_out, sc$label)
#   all_counts[[sc$label]] <- result$counts
#   all_dets[[sc$label]]   <- result$detections
#   cat("Done:", sc$label, "\n")
# }
# 
# count_table <- bind_rows(all_counts)
# det_table   <- bind_rows(all_dets)
# 
# # ── View results ──────────────────────────────────────────────────────────────
# print(count_table)
# print(det_table)
# 
# # # ── Save as CSVs ──────────────────────────────────────────────────────────────
# # write.csv(count_table, file.path(fig_dir, "supp_count_summary.csv"),
# #           row.names = FALSE)
# # write.csv(det_table,   file.path(fig_dir, "supp_detection_summary.csv"),
# #           row.names = FALSE)
# 
# 
# 
# # ── Figure: Simulation structure — one realisation ─────────────────────────
# library(ggplot2)
# library(patchwork)
# 
# sim_fig <- simulate_predprey(
#   nSites          = 500,
#   nVisits         = 10,
#   seed            = 1,
#   meanLambda_pred = 2,
#   meanLambda_meso = 4,
#   meanLambda_prey = 10,
#   meso            = TRUE,
#   use_landscape   = FALSE
# )
# 
# dat   <- sim_fig$sim_data
# truth <- sim_fig$truth
# 
# # ── Prediction curves — use truth values directly ─────────────────────────────
# env1_seq <- seq(min(dat$Env1), max(dat$Env1), length.out = 200)
# env2_seq <- seq(min(dat$Env2), max(dat$Env2), length.out = 200)
# pred_seq <- seq(0, max(truth$N_pred), length.out = 200)
# 
# curve_env1_pred <- data.frame(
#   x = env1_seq,
#   y = exp(truth$a_int + truth$a_env * env1_seq)
# )
# curve_env2_prey <- data.frame(
#   x = env2_seq,
#   y = exp(truth$b_int + truth$b_env * env2_seq +
#             truth$b_pred * mean(truth$N_pred) +
#             truth$b_meso * mean(truth$N_meso))
# )
# curve_pred_meso <- data.frame(
#   x = pred_seq,
#   y = exp(truth$c_int + truth$c_pred * pred_seq)
# )
# curve_pred_prey <- data.frame(
#   x = pred_seq,
#   y = exp(truth$b_int + truth$b_pred * pred_seq +
#             truth$b_meso * mean(truth$N_meso))
# )
# 
# # ── Plot data frame ───────────────────────────────────────────────────────────
# df_main <- data.frame(
#   Env1      = dat$Env1,
#   Env2      = dat$Env2,
#   N_pred    = truth$N_pred,
#   N_meso    = truth$N_meso,
#   N_prey    = truth$N_prey,
#   maxC_pred = apply(dat$C_pred, 1, max),
#   maxC_meso = apply(dat$C_meso, 1, max),
#   maxC_prey = apply(dat$C_prey, 1, max),
#   det_pred  = rowMeans(dat$Y_pred),
#   det_meso  = rowMeans(dat$Y_meso),
#   det_prey  = rowMeans(dat$Y_prey)
# )
# 
# # ── Shared theme ──────────────────────────────────────────────────────────────
# sim_theme <- theme_bw(base_size = 11) +
#   theme(
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(colour = "grey92"),
#     plot.title       = element_text(face = "bold", size = 13),
#     axis.title       = element_text(size = 12),
#     axis.text        = element_text(size = 11)
#   )
# 
# col_pred <- "#D85A30AA"
# col_meso <- "#7F77DDAA"
# col_prey <- "#1D9E75AA"
# alpha_pt_state <- 0.25   # for panels A-C with many overlapping points
# alpha_pt_obs   <- 0.45   # for panels D-E where colour distinction matters
# 
# 
# 
# # ── Panel A: Env1 → predator ──────────────────────────────────────────────────
# pA <- ggplot(df_main, aes(x = Env1, y = N_pred)) +
#   geom_point(colour = col_pred, alpha = alpha_pt, size = 2) +
#   geom_line(data = curve_env1_pred, aes(x = x, y = y),
#             colour = col_pred, linewidth = 1.2) +
#   labs(title = "A. Env\u2081 \u2192 Predator abundance (\u03b2 = 0.8)",
#        x = "Environmental covariate 1",
#        y = "True predator abundance (N)") +
#   sim_theme
# 
# # ── Panel B: N_pred → mesopredator ───────────────────────────────────────────
# pB <- ggplot(df_main, aes(x = N_pred, y = N_meso)) +
#   geom_point(colour = col_meso, alpha = alpha_pt, size = 2) +
#   geom_line(data = curve_pred_meso, aes(x = x, y = y),
#             colour = col_meso, linewidth = 1.2) +
#   labs(title = "B. Predator \u2192 Mesopredator (c_pred = \u22120.2)",
#        x = "True predator abundance (N)",
#        y = "True mesopredator abundance (N)") +
#   sim_theme
# 
# # ── Panel C: N_pred → prey ────────────────────────────────────────────────────
# pC <- ggplot(df_main, aes(x = N_pred, y = N_prey)) +
#   geom_point(colour = col_prey, alpha = alpha_pt, size = 2) +
#   geom_line(data = curve_pred_prey, aes(x = x, y = y),
#             colour = col_prey, linewidth = 1.2) +
#   labs(title = "C. Predator \u2192 Prey (b_pred = \u22120.25)",
#        x = "True predator abundance (N)",
#        y = "True prey abundance (N)") +
#   sim_theme
# 
# # ── Panel D: True N vs max observed count ────────────────────────────────────
# df_obs <- rbind(
#   data.frame(N_true = df_main$N_pred, N_obs = df_main$maxC_pred,
#              Species = "Predator (p = 0.35)"),
#   data.frame(N_true = df_main$N_meso, N_obs = df_main$maxC_meso,
#              Species = "Mesopredator (p = 0.50)"),
#   data.frame(N_true = df_main$N_prey, N_obs = df_main$maxC_prey,
#              Species = "Prey (p = 0.65)")
# )
# df_obs$Species <- factor(df_obs$Species,
#                          levels = c("Predator (p = 0.35)",
#                                     "Mesopredator (p = 0.50)",
#                                     "Prey (p = 0.65)"))
# 
# pD <- ggplot(df_obs, aes(x = N_true, y = N_obs, colour = Species)) +
#   geom_point(alpha = alpha_pt_obs, size = 2) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed",
#               colour = "grey40", linewidth = 0.7) +
#   scale_colour_manual(values = c(
#     "Predator (p = 0.35)"     = col_pred,
#     "Mesopredator (p = 0.50)" = col_meso,
#     "Prey (p = 0.65)"         = col_prey)) +
#   labs(title = "D. Imperfect detection: true N vs. max observed count",
#        subtitle = "Dashed line = perfect detection",
#        x = "True abundance (N)",
#        y = "Maximum observed count",
#        colour = NULL) +
#   coord_cartesian(xlim = c(0, 30)) +
#   sim_theme +
#   theme(legend.position = "bottom",
#         plot.subtitle = element_text(size = 10, colour = "grey40"))
# 
# # ── Panel E: True N vs naive detection rate ───────────────────────────────────
# df_det <- rbind(
#   data.frame(N_true = df_main$N_pred, det_rate = df_main$det_pred,
#              Species = "Predator (p = 0.35)"),
#   data.frame(N_true = df_main$N_meso, det_rate = df_main$det_meso,
#              Species = "Mesopredator (p = 0.50)"),
#   data.frame(N_true = df_main$N_prey, det_rate = df_main$det_prey,
#              Species = "Prey (p = 0.65)")
# )
# df_det$Species <- factor(df_det$Species,
#                          levels = c("Predator (p = 0.35)",
#                                     "Mesopredator (p = 0.50)",
#                                     "Prey (p = 0.65)"))
# 
# pE <- ggplot(df_det, aes(x = N_true, y = det_rate, colour = Species)) +
#   geom_point(alpha = alpha_pt_obs, size = 2) + 
#   geom_hline(yintercept = 1, linetype = "dashed",
#              colour = "grey40", linewidth = 0.7) +
#   scale_colour_manual(values = c(
#     "Predator (p = 0.35)"     = col_pred,
#     "Mesopredator (p = 0.50)" = col_meso,
#     "Prey (p = 0.65)"         = col_prey)) +
#   labs(title = "E. Naive detection rate vs. true abundance",
#        subtitle = "Dashed line = ceiling (detection rate = 1)",
#        x = "True abundance (N)",
#        y = "Naive detection rate",
#        colour = NULL) +
#   coord_cartesian(xlim = c(0, 30)) +
#   sim_theme +
#   theme(legend.position = "bottom",
#         plot.subtitle = element_text(size = 10, colour = "grey40"))
# 
# # ── Assemble ──────────────────────────────────────────────────────────────────
# fig2 <- (pA | pB | pC) / (pD | pE) +
#   plot_layout(guides = "collect") + 
#   plot_annotation(
#     title    = "One realisation of the simulation (seed = 1, n = 500 sites, abundant regime)",
#     subtitle = "Points show one simulated dataset; curves show the expected value from the data-generating model.",
#     theme    = theme(
#       plot.title    = element_text(face = "bold", size = 14),
#       plot.subtitle = element_text(size = 11, colour = "grey40")
#     )
#   ) &
#   theme(legend.position = "bottom")
# 
# fig_dir <- "~/Dropbox/Skye PhD UQ projects/Skye Thesis/Ch3 Skye species interactions/skye_ch3_analysis/simulations/figures/integrated_SEM_v4/"
# 
# ggsave(file.path(fig_dir, "Fig2_simulation_realisation.png"),
#        fig2, width = 12, height = 8, dpi = 300)
# cat("Saved: Fig2_simulation_realisation.png\n")
