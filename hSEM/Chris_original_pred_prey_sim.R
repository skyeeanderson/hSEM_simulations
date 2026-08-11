# ============================================================
#  Predator-prey N-mixture SEM simulation and model comparison
#
#  Species:
#    - Predator:     abundance driven by Env1
#    - Prey:         abundance driven by Env2 + negative effect of predator
#    - Mesopredator: optional third species, abundance negatively affected
#                    by predator and negatively affecting prey.
#                    Adds two additional indirect pathways:
#                      Env1 → N_pred → N_meso → N_prey
#                      N_pred → N_meso → N_prey
#                    Activated via meso = TRUE in simulate_predprey()
#
#  Models fitted:
#    1.  Poisson GLM         (glmmTMB, site random effect)
#    2a. Poisson SEM         (piecewiseSEM, max count)
#    2b. Poisson SEM         (piecewiseSEM, null N-mix BLUPs)
#    2c. Poisson SEM         (piecewiseSEM, full N-mix BLUPs)
#    3.  N-mixture           (unmarked, separate species)
#    4.  Co-abundance N-mix  (NIMBLE, joint with interaction)
#    5.  Bayesian SEM N-mix  (NIMBLE, integrated observation model)
#
#  Note: models 2b, 2c, 4, and 5 account for imperfect detection.
#        Models 1 and 2a use raw or max counts and do not.
#        The mesopredator pathway is not currently implemented
#        in the model fitting functions — simulation only.
# ============================================================

library(glmmTMB)
library(piecewiseSEM)
library(unmarked)
library(nimble)
library(coda)

# ============================================================
#  SIMULATE DATA
# ============================================================

simulate_predprey <- function(
    
  # study design
  nSites  = 100,
  nVisits = 4,
  seed    = NULL,
  
  # predator parameters
  a_int  =  1.2,
  a_env  =  0.8,
  
  # prey parameters
  b_int  =  2.5,
  b_env  =  0.6,
  b_pred = -0.12,
  
  # mesopredator parameters (only used if meso = TRUE)
  meso   = FALSE,
  c_int  =  1.0,
  c_pred = -0.10,
  b_meso = -0.08,
  p_meso =  0.50,
  
  # detection probabilities
  p_pred = 0.35,
  p_prey = 0.65,
  
  # output control
  plot = TRUE
  
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # environmental covariates
  Env1 <- rnorm(nSites)
  Env2 <- rnorm(nSites)
  
  # predator
  lambda_pred <- exp(a_int + a_env * Env1)
  N_pred      <- rpois(nSites, lambda_pred)
  C_pred      <- matrix(
    rbinom(nSites * nVisits, size = rep(N_pred, nVisits), prob = p_pred),
    nrow = nSites, ncol = nVisits
  )
  
  # mesopredator (optional)
  if (meso) {
    lambda_meso <- exp(c_int + c_pred * N_pred)
    N_meso      <- rpois(nSites, lambda_meso)
    C_meso      <- matrix(
      rbinom(nSites * nVisits, size = rep(N_meso, nVisits), prob = p_meso),
      nrow = nSites, ncol = nVisits
    )
  }
  
  # prey
  if (meso) {
    lambda_prey <- exp(b_int + b_env * Env2 + b_pred * N_pred + b_meso * N_meso)
  } else {
    lambda_prey <- exp(b_int + b_env * Env2 + b_pred * N_pred)
  }
  N_prey <- rpois(nSites, lambda_prey)
  C_prey <- matrix(
    rbinom(nSites * nVisits, size = rep(N_prey, nVisits), prob = p_prey),
    nrow = nSites, ncol = nVisits
  )
  
  # effect decomposition
  indirect_Env1_via_pred <- a_env * b_pred
  total_Env1_on_N_prey   <- indirect_Env1_via_pred
  
  if (meso) {
    indirect_Env1_via_meso <- a_env * c_pred * b_meso
    indirect_pred_via_meso <- c_pred * b_meso
    total_Env1_on_N_prey   <- indirect_Env1_via_pred + indirect_Env1_via_meso
  }
  
  # bundle sim_data
  sim_data <- list(
    nSites  = nSites,
    nVisits = nVisits,
    C_pred  = C_pred,
    C_prey  = C_prey,
    Env1    = Env1,
    Env2    = Env2
  )
  if (meso) sim_data$C_meso <- C_meso
  
  # bundle truth
  truth <- list(
    a_int   = a_int,  a_env  = a_env,
    b_int   = b_int,  b_env  = b_env,  b_pred = b_pred,
    p_pred  = p_pred, p_prey = p_prey,
    N_pred  = N_pred, N_prey = N_prey,
    lambda_pred             = lambda_pred,
    lambda_prey             = lambda_prey,
    indirect_Env1_on_N_prey = indirect_Env1_via_pred,
    total_Env1_on_N_prey    = total_Env1_on_N_prey
  )
  if (meso) {
    truth$c_int                  <- c_int
    truth$c_pred                 <- c_pred
    truth$b_meso                 <- b_meso
    truth$p_meso                 <- p_meso
    truth$N_meso                 <- N_meso
    truth$lambda_meso            <- lambda_meso
    truth$indirect_Env1_via_meso <- indirect_Env1_via_meso
    truth$indirect_pred_via_meso <- indirect_pred_via_meso
  }
  
  # plots
  if (plot) {
    n_panels <- if (meso) 8 else 6
    par(mfrow = c(2, ceiling(n_panels / 2)))
    
    plot(Env1, N_pred, pch = 16, col = "#D85A30AA",
         xlab = "Env1", ylab = "True N_pred",
         main = "Env1 → predator abundance")
    lines(sort(Env1), exp(a_int + a_env * sort(Env1)),
          col = "#993C1D", lwd = 2)
    
    plot(Env2, N_prey, pch = 16, col = "#7F77DDAA",
         xlab = "Env2", ylab = "True N_prey",
         main = "Env2 → prey abundance")
    lines(sort(Env2),
          exp(b_int + b_env * sort(Env2) + b_pred * mean(N_pred)),
          col = "#534AB7", lwd = 2)
    
    plot(N_pred, N_prey, pch = 16, col = "#1D9E75AA",
         xlab = "True N_pred", ylab = "True N_prey",
         main = "N_pred → N_prey (direct)")
    lines(0:50, exp(b_int + b_pred * (0:50)),
          col = "#0F6E56", lwd = 2)
    
    if (meso) {
      plot(N_pred, N_meso, pch = 16, col = "#EF9F27AA",
           xlab = "True N_pred", ylab = "True N_meso",
           main = "N_pred → N_meso")
      lines(0:50, exp(c_int + c_pred * (0:50)),
            col = "#BA7517", lwd = 2)
      plot(N_meso, N_prey, pch = 16, col = "#5DCAA5AA",
           xlab = "True N_meso", ylab = "True N_prey",
           main = "N_meso → N_prey")
    }
    
    hist(N_pred, breaks = 20, col = "#F5C4B3", border = "white",
         main = "Distribution of N_pred", xlab = "True N_pred")
    hist(N_prey, breaks = 20, col = "#CECBF6", border = "white",
         main = "Distribution of N_prey", xlab = "True N_prey")
    
    plot(N_pred, apply(C_pred, 1, max), pch = 16, col = "#888780AA",
         xlab = "True N_pred", ylab = "Max observed C_pred",
         main = "Observation model (predator)")
    abline(0, 1, lty = 2)
    
    par(mfrow = c(1, 1))
  }
  
  list(sim_data = sim_data, truth = truth)
}

# ============================================================
#  FIT ALL MODELS
# ============================================================

fit_all_models <- function(sim_out,
                           nimble_iter   = 30000,
                           nimble_burnin = 10000,
                           nimble_chains = 3,
                           nimble_thin   = 5) {
  
  dat     <- sim_out$sim_data
  truth   <- sim_out$truth
  
  nSites  <- dat$nSites
  nVisits <- dat$nVisits
  Env1    <- dat$Env1
  Env2    <- dat$Env2
  C_pred  <- dat$C_pred
  C_prey  <- dat$C_prey
  
  max_pred <- apply(C_pred, 1, max)
  max_prey <- apply(C_prey, 1, max)
  
  results <- list()
  
  # ── helper: posterior summary ──────────────────────────────
  post_summary <- function(samples, params) {
    do.call(rbind, lapply(params, function(p) {
      draws <- as.vector(as.matrix(samples[, p]))
      data.frame(
        param = p,
        mean  = mean(draws),
        sd    = sd(draws),
        q2.5  = quantile(draws, 0.025),
        q97.5 = quantile(draws, 0.975),
        Rhat  = tryCatch(
          gelman.diag(samples[, p])$psrf[1],
          error = function(e) NA
        )
      )
    }))
  }
  
  get_post_mean <- function(samps, param) {
    mean(as.vector(as.matrix(samps[, param])))
  }
  
  # ── helper: fit psem from a pair of Nhat vectors ───────────
  fit_psem_blup <- function(Nhat_pred_in, Nhat_prey_in, label) {
    
    Nhat_pred_r <- pmax(round(Nhat_pred_in), 1L)
    Nhat_prey_r <- pmax(round(Nhat_prey_in), 1L)
    
    df <- data.frame(
      abund_pred = Nhat_pred_r,
      abund_prey = Nhat_prey_r,
      Env1       = Env1,
      Env2       = Env2
    )
    
    fit <- tryCatch(
      psem(
        glm(abund_pred ~ Env1,              data = df, family = poisson),
        glm(abund_prey ~ Env2 + abund_pred, data = df, family = poisson)
      ),
      error = function(e) {
        warning(sprintf("fit_psem_blup [%s]: psem failed — %s",
                        label, e$message))
        NULL
      }
    )
    
    if (is.null(fit)) {
      return(list(
        label              = label,
        model              = NULL,
        summary            = NULL,
        coef_table         = NULL,
        b_Env1_pred        = NA,
        b_Env2_prey        = NA,
        b_pred_prey        = NA,
        indirect_Env1_prey = NA
      ))
    }
    
    s <- summary(fit, .progressBar = FALSE)
    
    list(
      label              = label,
      model              = fit,
      summary            = s,
      coef_table         = s$coefficients,
      b_Env1_pred        = coef(fit[[1]])["Env1"],
      b_Env2_prey        = coef(fit[[2]])["Env2"],
      b_pred_prey        = coef(fit[[2]])["abund_pred"],
      indirect_Env1_prey = coef(fit[[1]])["Env1"] *
        coef(fit[[2]])["abund_pred"]
    )
  }
  
  # ==========================================================
  #  MODEL 1: Poisson GLM (glmmTMB, site random effect)
  #  Raw counts stacked long; site random intercept absorbs
  #  site-level abundance variation. Does not correct for
  #  detection — abundance index only, not absolute N.
  # ==========================================================
  message("\n--- Model 1: Poisson GLM (glmmTMB) ---")
  
  df_pred_long <- data.frame(
    count = as.vector(C_pred),
    Env1  = rep(Env1, nVisits),
    site  = rep(factor(1:nSites), nVisits)
  )
  df_prey_long <- data.frame(
    count = as.vector(C_prey),
    Env2  = rep(Env2, nVisits),
    site  = rep(factor(1:nSites), nVisits)
  )
  
  glmm_pred <- glmmTMB(count ~ Env1 + (1 | site),
                       data = df_pred_long, family = poisson)
  glmm_prey <- glmmTMB(count ~ Env2 + (1 | site),
                       data = df_prey_long, family = poisson)
  
  re_pred <- glmmTMB::ranef(glmm_pred)$cond$site[["(Intercept)"]]
  re_prey <- glmmTMB::ranef(glmm_prey)$cond$site[["(Intercept)"]]
  
  results$glmm <- list(
    model_pred       = glmm_pred,
    model_prey       = glmm_prey,
    abund_index_pred = exp(fixef(glmm_pred)$cond["(Intercept)"] + re_pred),
    abund_index_prey = exp(fixef(glmm_prey)$cond["(Intercept)"] + re_prey),
    coef_pred        = fixef(glmm_pred)$cond,
    coef_prey        = fixef(glmm_prey)$cond
  )
  
  # ==========================================================
  #  MODEL 3: N-mixture (unmarked) — fitted before SEMs
  #  because BLUP SEMs depend on these estimates.
  #
  #  3a. Null model — intercept only on occurrence and detection
  #  3b. Full model — env covariates on occurrence (data-generating
  #                   structure); used for parameter inference
  #
  #  BLUPs extracted via ranef() / bup() which condition on the
  #  observed count data at each site (empirical Bayes estimates).
  #  These vary across sites even under the null model because
  #  they are posterior means of N[i] | y[i,], not predict().
  # ==========================================================
  message("\n--- Model 3: N-mixture (unmarked) ---")
  
  K_pred <- max(C_pred) + 50
  K_prey <- max(C_prey) + 50
  
  # null
  umf_pred_null <- unmarkedFramePCount(y = C_pred)
  umf_prey_null <- unmarkedFramePCount(y = C_prey)
  
  nmix_pred_null <- pcount(~ 1 ~ 1, data = umf_pred_null, K = K_pred)
  nmix_prey_null <- pcount(~ 1 ~ 1, data = umf_prey_null, K = K_prey)
  
  # full
  umf_pred_full <- unmarkedFramePCount(
    y        = C_pred,
    siteCovs = data.frame(Env1 = Env1)
  )
  umf_prey_full <- unmarkedFramePCount(
    y        = C_prey,
    siteCovs = data.frame(Env2 = Env2)
  )
  
  nmix_pred_full <- pcount(~ 1 ~ Env1, data = umf_pred_full, K = K_pred)
  nmix_prey_full <- pcount(~ 1 ~ Env2, data = umf_prey_full, K = K_prey)
  
  # empirical Bayes BLUPs: posterior mean of N[i] | y[i,]
  # ranef() computes the full posterior; bup() extracts the mean
  Nhat_pred_null <- bup(ranef(nmix_pred_null), stat = "mean")
  Nhat_prey_null <- bup(ranef(nmix_prey_null), stat = "mean")
  Nhat_pred_full <- bup(ranef(nmix_pred_full), stat = "mean")
  Nhat_prey_full <- bup(ranef(nmix_prey_full), stat = "mean")
  
  results$nmixture <- list(
    # null
    model_pred_null = nmix_pred_null,
    model_prey_null = nmix_prey_null,
    Nhat_pred_null  = Nhat_pred_null,
    Nhat_prey_null  = Nhat_prey_null,
    # full
    model_pred_full = nmix_pred_full,
    model_prey_full = nmix_prey_full,
    Nhat_pred_full  = Nhat_pred_full,
    Nhat_prey_full  = Nhat_prey_full,
    coef_pred       = coef(nmix_pred_full),
    coef_prey       = coef(nmix_prey_full),
    p_pred_est      = plogis(coef(nmix_pred_full)["p(Int)"]),
    p_prey_est      = plogis(coef(nmix_prey_full)["p(Int)"])
  )
  
  # ==========================================================
  #  MODEL 2a: Poisson SEM — max count as abundance proxy
  #  No correction for detection. Expect attenuation in b_pred
  #  particularly for predator which has low detection.
  # ==========================================================
  message("\n--- Model 2a: Poisson SEM (max count) ---")
  
  df_sem_maxcount <- data.frame(
    abund_pred = max_pred,
    abund_prey = max_prey,
    Env1       = Env1,
    Env2       = Env2
  )
  
  psem_maxcount_fit <- psem(
    glm(abund_pred ~ Env1,              data = df_sem_maxcount, family = poisson),
    glm(abund_prey ~ Env2 + abund_pred, data = df_sem_maxcount, family = poisson)
  )
  
  s_mc <- summary(psem_maxcount_fit, .progressBar = FALSE)
  
  results$psem_maxcount <- list(
    label              = "max count",
    model              = psem_maxcount_fit,
    summary            = s_mc,
    coef_table         = s_mc$coefficients,
    b_Env1_pred        = coef(psem_maxcount_fit[[1]])["Env1"],
    b_Env2_prey        = coef(psem_maxcount_fit[[2]])["Env2"],
    b_pred_prey        = coef(psem_maxcount_fit[[2]])["abund_pred"],
    indirect_Env1_prey = coef(psem_maxcount_fit[[1]])["Env1"] *
      coef(psem_maxcount_fit[[2]])["abund_pred"]
  )
  
  # ==========================================================
  #  MODEL 2b: Poisson SEM — null N-mixture BLUPs
  #  Empirical Bayes Nhats from intercept-only N-mixture.
  #  Detection corrected but no env structure assumed.
  #  BLUPs vary across sites (conditional on observed counts)
  #  but shrink more toward global mean than full model BLUPs.
  # ==========================================================
  message("\n--- Model 2b: Poisson SEM (null N-mixture BLUPs) ---")
  
  results$psem_blup_null <- fit_psem_blup(
    Nhat_pred_null, Nhat_prey_null, "null N-mixture BLUPs"
  )
  
  # ==========================================================
  #  MODEL 2c: Poisson SEM — full N-mixture BLUPs
  #  Empirical Bayes Nhats from full covariate N-mixture.
  #  Best two-stage approach: detection corrected and env
  #  structure matches data-generating model. Still ignores
  #  uncertainty in Nhat (treats point estimates as known).
  # ==========================================================
  message("\n--- Model 2c: Poisson SEM (full N-mixture BLUPs) ---")
  
  results$psem_blup_full <- fit_psem_blup(
    Nhat_pred_full, Nhat_prey_full, "full N-mixture BLUPs"
  )
  
  # ==========================================================
  #  MODEL 4: Co-abundance N-mixture SEM (NIMBLE)
  #  Both species jointly modelled; N_pred enters prey lambda
  #  as a latent node — fully integrated, no two-stage bias.
  #  No shared site random effect (cf. Model 5).
  # ==========================================================
  message("\n--- Model 4: Co-abundance N-mixture SEM (NIMBLE) ---")
  
  coabund_code <- nimbleCode({
    
    # priors: predator
    a_int  ~ dnorm(0, sd = 2)
    a_env  ~ dnorm(0, sd = 2)
    p_pred ~ dbeta(1, 1)
    
    # priors: prey
    b_int  ~ dnorm(0, sd = 2)
    b_env  ~ dnorm(0, sd = 2)
    b_pred ~ dnorm(0, sd = 2)
    p_prey ~ dbeta(1, 1)
    
    for (i in 1:nSites) {
      
      # predator N-mixture
      log(lambda_pred[i]) <- a_int + a_env * Env1[i]
      N_pred[i] ~ dpois(lambda_pred[i])
      for (j in 1:nVisits) {
        C_pred[i, j] ~ dbin(p_pred, N_pred[i])
      }
      
      # prey N-mixture: N_pred[i] as latent covariate
      log(lambda_prey[i]) <- b_int + b_env * Env2[i] + b_pred * N_pred[i]
      N_prey[i] ~ dpois(lambda_prey[i])
      for (j in 1:nVisits) {
        C_prey[i, j] ~ dbin(p_prey, N_prey[i])
      }
    }
    
    # derived: effect decomposition
    indirect_Env1_on_prey <- a_env * b_pred
    total_Env1_on_prey    <- indirect_Env1_on_prey
  })
  
  nimble_constants <- list(
    nSites  = nSites,
    nVisits = nVisits,
    Env1    = Env1,
    Env2    = Env2
  )
  nimble_data <- list(
    C_pred = C_pred,
    C_prey = C_prey
  )
  nimble_inits4 <- list(
    a_int  = 0, a_env  = 0, p_pred = 0.5,
    b_int  = 0, b_env  = 0, b_pred = 0,  p_prey = 0.5,
    N_pred = max_pred + 1,
    N_prey = max_prey + 1
  )
  
  m4     <- nimbleModel(coabund_code, nimble_constants,
                        nimble_data,  nimble_inits4)
  cm4    <- compileNimble(m4)
  conf4  <- configureMCMC(cm4, monitors = c(
    "a_int", "a_env", "p_pred",
    "b_int", "b_env", "b_pred", "p_prey",
    "indirect_Env1_on_prey", "total_Env1_on_prey"
  ))
  mcmc4  <- buildMCMC(conf4)
  cmcmc4 <- compileNimble(mcmc4, project = m4)
  samps4 <- runMCMC(cmcmc4,
                    niter             = nimble_iter,
                    nburnin           = nimble_burnin,
                    nchains           = nimble_chains,
                    thin              = nimble_thin,
                    samplesAsCodaMCMC = TRUE)
  
  params4 <- c("a_int", "a_env", "p_pred",
               "b_int", "b_env", "b_pred", "p_prey",
               "indirect_Env1_on_prey", "total_Env1_on_prey")
  
  results$coabund_nimble <- list(
    samples = samps4,
    summary = post_summary(samps4, params4)
  )
  
  # ==========================================================
  #  MODEL 5: Bayesian integrated N-mixture SEM (NIMBLE)
  #  As Model 4 but adds a shared site-level random effect
  #  across both species equations, capturing unmeasured site
  #  heterogeneity that affects both species jointly.
  # ==========================================================
  message("\n--- Model 5: Bayesian integrated N-mixture SEM (NIMBLE) ---")
  
  sem_code <- nimbleCode({
    
    # priors: predator
    a_int  ~ dnorm(0, sd = 2)
    a_env  ~ dnorm(0, sd = 2)
    p_pred ~ dbeta(1, 1)
    
    # priors: prey
    b_int  ~ dnorm(0, sd = 2)
    b_env  ~ dnorm(0, sd = 2)
    b_pred ~ dnorm(0, sd = 2)
    p_prey ~ dbeta(1, 1)
    
    # shared site random effect
    sigma_re ~ dexp(1)
    for (i in 1:nSites) {
      re[i] ~ dnorm(0, sd = sigma_re)
    }
    
    for (i in 1:nSites) {
      
      # predator
      log(lambda_pred[i]) <- a_int + a_env * Env1[i] + re[i]
      N_pred[i] ~ dpois(lambda_pred[i])
      for (j in 1:nVisits) {
        C_pred[i, j] ~ dbin(p_pred, N_pred[i])
      }
      
      # prey: N_pred[i] is latent node sampled above
      log(lambda_prey[i]) <- b_int + b_env * Env2[i] +
        b_pred * N_pred[i] + re[i]
      N_prey[i] ~ dpois(lambda_prey[i])
      for (j in 1:nVisits) {
        C_prey[i, j] ~ dbin(p_prey, N_prey[i])
      }
    }
    
    # derived: effect decomposition
    indirect_Env1_on_prey <- a_env * b_pred
    total_Env1_on_prey    <- indirect_Env1_on_prey
  })
  
  nimble_inits5 <- list(
    a_int    = 0, a_env  = 0, p_pred = 0.5,
    b_int    = 0, b_env  = 0, b_pred = 0,  p_prey = 0.5,
    sigma_re = 0.3,
    re       = rep(0, nSites),
    N_pred   = max_pred + 1,
    N_prey   = max_prey + 1
  )
  
  m5     <- nimbleModel(sem_code, nimble_constants,
                        nimble_data, nimble_inits5)
  cm5    <- compileNimble(m5)
  conf5  <- configureMCMC(cm5, monitors = c(
    "a_int", "a_env", "p_pred",
    "b_int", "b_env", "b_pred", "p_prey",
    "sigma_re",
    "indirect_Env1_on_prey", "total_Env1_on_prey"
  ))
  mcmc5  <- buildMCMC(conf5)
  cmcmc5 <- compileNimble(mcmc5, project = m5)
  samps5 <- runMCMC(cmcmc5,
                    niter             = nimble_iter,
                    nburnin           = nimble_burnin,
                    nchains           = nimble_chains,
                    thin              = nimble_thin,
                    samplesAsCodaMCMC = TRUE)
  
  params5 <- c("a_int", "a_env", "p_pred",
               "b_int", "b_env", "b_pred", "p_prey",
               "sigma_re",
               "indirect_Env1_on_prey", "total_Env1_on_prey")
  
  results$sem_nimble <- list(
    samples = samps5,
    summary = post_summary(samps5, params5)
  )
  
  # ==========================================================
  #  COMPARISON TABLE
  # ==========================================================
  message("\n--- Building comparison table ---")
  
  comparison <- data.frame(
    estimand = c(
      "mean N_pred",
      "mean N_prey",
      "Env1 → N_pred (a_env)",
      "Env2 → N_prey (b_env)",
      "N_pred → N_prey (b_pred)",
      "indirect Env1 → N_prey",
      "p_pred",
      "p_prey"
    ),
    truth = c(
      mean(truth$N_pred),
      mean(truth$N_prey),
      truth$a_env,
      truth$b_env,
      truth$b_pred,
      truth$indirect_Env1_on_N_prey,
      truth$p_pred,
      truth$p_prey
    ),
    glmm = c(
      mean(results$glmm$abund_index_pred),
      mean(results$glmm$abund_index_prey),
      results$glmm$coef_pred["Env1"],
      results$glmm$coef_prey["Env2"],
      NA, NA, NA, NA
    ),
    psem_maxcount = c(
      mean(max_pred),
      mean(max_prey),
      results$psem_maxcount$b_Env1_pred,
      results$psem_maxcount$b_Env2_prey,
      results$psem_maxcount$b_pred_prey,
      results$psem_maxcount$indirect_Env1_prey,
      NA, NA
    ),
    psem_blup_null = c(
      mean(Nhat_pred_null),
      mean(Nhat_prey_null),
      results$psem_blup_null$b_Env1_pred,
      results$psem_blup_null$b_Env2_prey,
      results$psem_blup_null$b_pred_prey,
      results$psem_blup_null$indirect_Env1_prey,
      NA, NA
    ),
    psem_blup_full = c(
      mean(Nhat_pred_full),
      mean(Nhat_prey_full),
      results$psem_blup_full$b_Env1_pred,
      results$psem_blup_full$b_Env2_prey,
      results$psem_blup_full$b_pred_prey,
      results$psem_blup_full$indirect_Env1_prey,
      NA, NA
    ),
    nmixture = c(
      mean(Nhat_pred_full),
      mean(Nhat_prey_full),
      coef(results$nmixture$model_pred_full)["lam(Env1)"],
      coef(results$nmixture$model_prey_full)["lam(Env2)"],
      NA, NA,
      results$nmixture$p_pred_est,
      results$nmixture$p_prey_est
    ),
    coabund = c(
      get_post_mean(samps4, "a_int"),
      get_post_mean(samps4, "b_int"),
      get_post_mean(samps4, "a_env"),
      get_post_mean(samps4, "b_env"),
      get_post_mean(samps4, "b_pred"),
      get_post_mean(samps4, "indirect_Env1_on_prey"),
      get_post_mean(samps4, "p_pred"),
      get_post_mean(samps4, "p_prey")
    ),
    sem_bayes = c(
      get_post_mean(samps5, "a_int"),
      get_post_mean(samps5, "b_int"),
      get_post_mean(samps5, "a_env"),
      get_post_mean(samps5, "b_env"),
      get_post_mean(samps5, "b_pred"),
      get_post_mean(samps5, "indirect_Env1_on_prey"),
      get_post_mean(samps5, "p_pred"),
      get_post_mean(samps5, "p_prey")
    )
  )
  
  results$comparison <- comparison
  
  message("\nDone.")
  results
}

# ============================================================
#  RUN — 5 REPLICATES
# ============================================================

seeds <- c(1980, 1234, 5678, 4321, 9999)

all_reps <- lapply(seq_along(seeds), function(r) {
  
  message("\n============================================================")
  message(sprintf("  REPLICATE %d / %d  (seed = %d)", r, length(seeds), seeds[r]))
  message("============================================================")
  
  sim_out <- simulate_predprey(seed = seeds[r], plot = FALSE)
  fit     <- fit_all_models(sim_out)
  
  list(
    seed       = seeds[r],
    sim_out    = sim_out,
    fit        = fit,
    comparison = fit$comparison
  )
})

# ============================================================
#  SUMMARISE ACROSS REPLICATES
# ============================================================

model_cols <- c("glmm", "psem_maxcount", "psem_blup_null",
                "psem_blup_full", "nmixture", "coabund", "sem_bayes")

rep_tables <- do.call(rbind, lapply(seq_along(all_reps), function(r) {
  tab      <- all_reps[[r]]$comparison
  tab$rep  <- r
  tab$seed <- all_reps[[r]]$seed
  tab
}))

summary_across_reps <- do.call(rbind, lapply(
  unique(rep_tables$estimand), function(est) {
    sub <- rep_tables[rep_tables$estimand == est, ]
    row <- data.frame(estimand = est, truth = sub$truth[1])
    for (col in model_cols) {
      row[[paste0(col, "_mean")]] <- mean(sub[[col]], na.rm = TRUE)
      row[[paste0(col, "_sd")]]   <- sd(sub[[col]],   na.rm = TRUE)
    }
    row
  }
))

message("\n=== Mean estimates across replicates ===")
print(summary_across_reps, digits = 3, row.names = FALSE)

saveRDS(all_reps,            "predprey_all_reps.rds")
saveRDS(summary_across_reps, "predprey_summary.rds")

message("\nResults saved to predprey_all_reps.rds and predprey_summary.rds")
