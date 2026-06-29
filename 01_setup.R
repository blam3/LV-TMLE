# ==============================================================================
# FILE: 01_setup.R  (REVISED)
# PURPOSE: Estimator functions for the LV-TMLE simulation.
#
# WHAT CHANGED vs the original (and why):
#   1. Per-scenario true psi. The global TRUE_PSI = 2.0 is GONE. psi now varies
#      by (outcome x latent) and is computed by Monte Carlo in 00_dgp_variants.R,
#      cached to disk, and attached to sim_grid. Bias/coverage downstream use the
#      scenario's own true_psi.
#   2. One shared TMLE engine. `shift_tmle()` implements a single cross-fitted
#      shift-TMLE with the CORRECT efficient influence function (it keeps the
#      H*(Y-Q) score term the old code dropped). Naive, Oracle, and each
#      LV-TMLE imputation all call this same engine; they differ ONLY in which
#      L-values they receive. No more stacking -> no induced attenuation; and
#      because each fit sees an n-row (un-stacked) dataset, SuperLearner's inner
#      CV no longer leaks and the small-sample learner switch keys on the true N.
#   3. Oracle estimator: the true latent L_true plugged into `shift_tmle()`. This
#      is the performance ceiling (no measurement error).
#   4. Corrected LV-TMLE (`lv_tmle_mi`): proper multiple-imputation over the
#      latent posterior with Rubin's rules. Each of M draws is a COMPLETE n-row
#      dataset (one L per subject), fit separately; point estimates and
#      IF-variances are combined so the BETWEEN-imputation variance (the latent
#      uncertainty) enters the SE. The self-defeating naive-anchored truncation
#      is gone; truncation is a fixed, documented bound.
#
#   Requires 00_dgp_variants.R (DGP + true-psi machinery) and the SL stack.
# ==============================================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(data.table)
  library(SuperLearner)
  library(earth)
  library(rpart)
  library(dplyr)
})

# Pull in the modular DGPs, build_sim_grid(), build_psi_lookup(), attach_true_psi()
source("00_dgp_variants.R")

# ------------------------------------------------------------------------------
# 0. SIMULATION GRID WITH PER-SCENARIO TRUE PSI
# ------------------------------------------------------------------------------
# psi depends only on (outcome x latent); we cache the MC result so that the
# (potentially many) SLURM workers do not each recompute a 2e7-draw integral.
#
# RECOMMENDED: precompute the cache ONCE interactively to avoid a write race:
#     source("00_dgp_variants.R"); saveRDS(build_psi_lookup(n_mc = 2e7), "psi_lookup.rds")
# Then every worker just loads it below.
PSI_CACHE <- "psi_lookup.rds"
if (file.exists(PSI_CACHE)) {
  psi_lookup <- readRDS(PSI_CACHE)
} else {
  message("psi_lookup.rds not found -- computing true psi by Monte Carlo (slow). ",
          "Precompute once and save to avoid this on every worker.")
  psi_lookup <- build_psi_lookup(n_mc = 2e7)
  try(saveRDS(psi_lookup, PSI_CACHE), silent = TRUE)
}

# Factorial grid; each row carries its own true_psi.
sim_grid <- attach_true_psi(build_sim_grid(), psi_lookup)
sim_grid <- sim_grid[order(sim_grid$Scenario_ID), ]
rownames(sim_grid) <- NULL
# (Columns: outcome, latent, n_size, meas_error, indicators, Scenario_ID, true_psi)

# ------------------------------------------------------------------------------
# 1. SAFE SUPERLEARNER + SAFE PREDICTION  (unchanged in spirit; robust wrappers)
# ------------------------------------------------------------------------------
.lv_diag_env <- new.env(parent = emptyenv())

lv_diag_reset <- function() {
  .lv_diag_env$calls <- 0L
  .lv_diag_env$fallback_short_n <- 0L
  .lv_diag_env$fallback_low_variance <- 0L
  .lv_diag_env$superlearner_error <- 0L
  .lv_diag_env$glm_fallback <- 0L
  .lv_diag_env$mean_fallback <- 0L
  .lv_diag_env$discrete_sl <- 0L
  .lv_diag_env$predict_fallback <- 0L
  invisible(NULL)
}

lv_diag_bump <- function(name) {
  if (is.null(.lv_diag_env$calls)) lv_diag_reset()
  .lv_diag_env[[name]] <- as.integer(.lv_diag_env[[name]]) + 1L
  invisible(NULL)
}

lv_diag_snapshot <- function() {
  if (is.null(.lv_diag_env$calls)) lv_diag_reset()
  list(
    sl_calls = .lv_diag_env$calls,
    sl_fallback_short_n = .lv_diag_env$fallback_short_n,
    sl_fallback_low_variance = .lv_diag_env$fallback_low_variance,
    sl_superlearner_error = .lv_diag_env$superlearner_error,
    sl_glm_fallback = .lv_diag_env$glm_fallback,
    sl_mean_fallback = .lv_diag_env$mean_fallback,
    sl_discrete = .lv_diag_env$discrete_sl,
    sl_predict_fallback = .lv_diag_env$predict_fallback
  )
}

lv_diag_reset()

with_lavaan_core_fallback <- function(expr) {
  ns <- asNamespace("parallel")
  old_detect <- get("detectCores", envir = ns)
  detected <- try(c(old_detect(), old_detect(logical = FALSE)), silent = TRUE)
  needs_patch <- inherits(detected, "try-error") ||
    length(detected) == 0L || any(!is.finite(detected))

  if (needs_patch) {
    # lavaan 0.6-21 validates ncpus against parallel::detectCores(); in this
    # managed shell detectCores() can return NA, so force serial SEM fitting.
    was_locked <- bindingIsLocked("detectCores", ns)
    if (was_locked) unlockBinding("detectCores", ns)
    assign("detectCores", function(...) 1L, envir = ns)
    if (was_locked) lockBinding("detectCores", ns)
    on.exit({
      if (was_locked) unlockBinding("detectCores", ns)
      assign("detectCores", old_detect, envir = ns)
      if (was_locked) lockBinding("detectCores", ns)
    }, add = TRUE)
  }

  force(expr)
}

lv_sl_library <- function(n) {
  # SL.gam was dropped after pilot runs repeatedly removed it with singularity
  # errors. Keep the active library stable until a replacement is justified.
  if (n <= 100) {
    c("SL.glm", "SL.mean", "SL.rpart", "SL.earth")
  } else {
    c("SL.glm", "SL.earth", "SL.mean", "SL.rpart")
  }
}

robust_predict <- function(model, newdata, covars) {
  if (is.list(model) && !is.null(model$fallback_mean))
    return(rep(model$fallback_mean, nrow(newdata)))
  if (is.numeric(model) && length(model) == 1)
    return(rep(model, nrow(newdata)))

  out <- tryCatch({
    if (inherits(model, "SuperLearner")) {
      as.numeric(predict(model, newdata[, covars, drop = FALSE])$pred)
    } else {
      predict(model, newdata, type = "response")
    }
  }, error = function(e) NULL)

  if (is.null(out) || any(is.na(out)) || length(out) != nrow(newdata)) {
    lv_diag_bump("predict_fallback")
    fallback <- if (!is.null(model$Y)) mean(model$Y, na.rm = TRUE) else 0
    return(rep(fallback, nrow(newdata)))
  }
  as.numeric(out)
}

safe_SL <- function(Y, X, family, SL.library) {
  lv_diag_bump("calls")
  clean <- complete.cases(X, Y)
  Yc <- Y[clean]; Xc <- X[clean, , drop = FALSE]

  if (nrow(Xc) < 10) {
    lv_diag_bump("fallback_short_n")
    return(list(fallback_mean = mean(Yc, na.rm = TRUE), Y = Yc))
  }
  if (var(Yc) < 1e-8) {
    lv_diag_bump("fallback_low_variance")
    return(list(fallback_mean = mean(Yc, na.rm = TRUE), Y = Yc))
  }

  fit <- try(SuperLearner(Y = Yc, X = Xc, family = family,
                          SL.library = SL.library, verbose = FALSE), silent = TRUE)

  if (inherits(fit, "try-error") || all(fit$SL.predict == 0)) {
    lv_diag_bump("superlearner_error")
    fglm <- try(glm(Yc ~ ., data = Xc, family = family), silent = TRUE)
    if (!inherits(fglm, "try-error")) {
      lv_diag_bump("glm_fallback")
      return(fglm)
    }
    lv_diag_bump("mean_fallback")
    return(list(fallback_mean = mean(Yc, na.rm = TRUE), Y = Yc))
  }

  # Discrete SuperLearner in small samples (now keyed on the TRUE n, because the
  # datasets passed here are un-stacked, size n -- not n*M).
  if (nrow(Xc) < 1000) {
    best <- which.min(fit$cvRisk)
    fit$coef <- rep(0, length(fit$coef)); fit$coef[best] <- 1
    lv_diag_bump("discrete_sl")
  }
  fit
}

# ------------------------------------------------------------------------------
# 2. THE SHARED TMLE ENGINE  -- cross-fitted shift-TMLE with the correct EIF
# ------------------------------------------------------------------------------
# Estimand: the modified-treatment-policy contrast
#     psi = E[Y(A + delta)] - E[Y(A)] = E[Q(A+delta, X)] - E[Y].
# EIF (the thing the old SE was missing): for the shifted mean psi1 = E[Q(A+delta,X)],
#     D1(O) = H(A,X)*(Y - Q(A,X)) + Q(A+delta,X) - psi1,   H(a,x) = g(a-delta|x)/g(a|x),
# and for the contrast psi = psi1 - E[Y], subtract the (Y - E[Y]) component:
#     IC(O) = H(A,X)*(Y - Q*(A,X)) + Q*(A+delta,X) - psi1 - (Y - Ybar).
# SE = sd(IC)/sqrt(n). Note H*(Y - Q*) IS included (the original dropped it).
#
# Inputs:
#   Y  : outcome vector (length n)
#   A  : exposure values to treat as the (point-mass) latent for this fit
#   X  : data.frame/data.table of confounders (here W, Z)
# Returns: list(psi, se, ic, eps).
shift_tmle <- function(Y, A, X, sl_lib, delta = 1, V = 5,
                       gbound = 1e-3, hcap = 30, folds = NULL) {
  n <- length(Y)
  X <- as.data.frame(X)
  base <- data.frame(Y = Y, L_sim = A, X)
  covars_Q <- c("L_sim", colnames(X))
  covars_g <- colnames(X)

  if (is.null(folds)) {
    V <- max(2L, min(as.integer(V), n))
    folds <- sample(rep(1:V, length.out = n))
  } else {
    if (length(folds) != n) stop("folds must have length n")
    folds <- as.integer(folds)
    V <- length(unique(folds[is.finite(folds)]))
    if (V < 2L) stop("folds must contain at least two folds")
  }

  Q_A <- Q_As <- muG <- sigG <- rep(NA_real_, n)

  for (v in 1:V) {
    tr <- which(folds != v); va <- which(folds == v)
    if (length(tr) < 10 || length(va) == 0) next

    # Q-model: E[Y | L_sim, X]  (cross-fitted)
    qfit <- safe_SL(Y = base$Y[tr], X = base[tr, covars_Q, drop = FALSE],
                    family = gaussian(), SL.library = sl_lib)
    Q_A[va]  <- robust_predict(qfit, base[va, , drop = FALSE], covars_Q)
    bshift   <- base[va, , drop = FALSE]; bshift$L_sim <- bshift$L_sim + delta
    Q_As[va] <- robust_predict(qfit, bshift, covars_Q)

    # g-model: E[L_sim | X]; Gaussian working density g(a|x) = dnorm(a; muG, sigG)
    gfit <- safe_SL(Y = base$L_sim[tr], X = base[tr, covars_g, drop = FALSE],
                    family = gaussian(), SL.library = sl_lib)
    muG[va] <- robust_predict(gfit, base[va, , drop = FALSE], covars_g)
    mu_tr   <- robust_predict(gfit, base[tr, , drop = FALSE], covars_g)
    s <- sd(base$L_sim[tr] - mu_tr, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-6) s <- 1
    sigG[va] <- s
  }

  Q_A[is.na(Q_A)]   <- mean(Y, na.rm = TRUE)
  Q_As[is.na(Q_As)] <- mean(Y, na.rm = TRUE)
  muG[is.na(muG)]   <- mean(A, na.rm = TRUE)
  sigG[is.na(sigG)] <- { s <- sd(A, na.rm = TRUE); if (!is.finite(s) || s <= 0) 1 else s }

  g_A  <- pmax(dnorm(A,         muG, sigG), gbound)
  g_Am <- pmax(dnorm(A - delta, muG, sigG), gbound)
  g_Ap <- pmax(dnorm(A + delta, muG, sigG), gbound)
  H_A  <- pmin(g_Am / g_A, hcap)     # clever covariate at the observed A
  H_As <- pmin(g_A  / g_Ap, hcap)    # clever covariate at A + delta

  # One-step fluctuation (identity link; Y unbounded): Y ~ -1 + H, offset = Q
  eps <- tryCatch({
    cf <- coef(lm(Y ~ -1 + H_A, offset = Q_A))
    if (!is.finite(cf)) 0 else as.numeric(cf)
  }, error = function(e) 0)

  Qstar_A  <- Q_A  + eps * H_A
  Qstar_As <- Q_As + eps * H_As

  psi1 <- mean(Qstar_As, na.rm = TRUE)
  Ybar <- mean(Y, na.rm = TRUE)
  psi  <- psi1 - Ybar

  ic <- (H_A * (Y - Qstar_A) + Qstar_As - psi1) - (Y - Ybar)   # full EIF
  se <- sd(ic, na.rm = TRUE) / sqrt(n)

  list(psi = psi, se = se, ic = ic, eps = eps)
}

# ------------------------------------------------------------------------------
# 3. SEM FIT + LATENT POSTERIOR  (one fit serves the SEM estimator AND the draws)
# ------------------------------------------------------------------------------
# Returns the SEM's L->Y coefficient (its estimand for psi), Bartlett factor
# scores (the posterior mean mu_i), and the posterior SD of L given indicators.
# For this homoskedastic, complete-data measurement model the posterior SD is
# constant across subjects, so a scalar post_sd is correct here.
# NOTE: the ordinal branch fits an ordered CFA (WLSMV). Its post_sd is an
# APPROXIMATION on the latent-response scale and should be validated before
# publication; the continuous branch is exact for this model.
fit_sem_posterior <- function(dt, indicators_type) {
  dt <- copy(dt)
  dt[, Z_sq := Z^2]
  ordered_vars <- if (indicators_type == "ordinal") c("I1", "I2", "I3") else NULL

  cfa_model <- '
    L =~ I1 + I2 + I3
    L ~ W + Z + Z_sq
    Y ~ L + W + Z + Z_sq
  '
  fit <- try(with_lavaan_core_fallback(suppressWarnings(
    sem(cfa_model, data = as.data.frame(dt), ordered = ordered_vars,
        std.lv = FALSE, auto.fix.first = TRUE, ncpus = 1L)
  )), silent = TRUE)

  if (inherits(fit, "try-error") || !lavInspect(fit, "converged"))
    return(list(converged = FALSE, reason = "SEM_no_converge"))

  pe  <- parameterEstimates(fit)
  rec <- pe[pe$lhs == "Y" & pe$op == "~" & pe$rhs == "L", ]
  psi_sem <- if (nrow(rec)) rec$est else NA_real_
  se_sem  <- if (nrow(rec)) rec$se  else NA_real_
  if (is.na(psi_sem)) return(list(converged = FALSE, reason = "SEM_coef_NA"))

  mu <- try(lavPredict(fit, method = "bartlett")[, "L"], silent = TRUE)
  if (inherits(mu, "try-error") || any(!is.finite(mu)))
    mu <- try(lavPredict(fit, method = "regression")[, "L"], silent = TRUE)
  if (inherits(mu, "try-error") || any(!is.finite(mu)))
    return(list(converged = FALSE, reason = "Factor_scores_NA"))

  post_sd <- tryCatch({
    em    <- lavInspect(fit, "est")
    lam   <- em$lambda[c("I1", "I2", "I3"), "L"]
    th    <- pmax(diag(em$theta)[c("I1", "I2", "I3")], 1e-3)
    psiLL <- pmax(em$psi["L", "L"], 1e-3)
    sqrt(1 / (1 / psiLL + sum(lam^2 / th)))
  }, error = function(e) NA_real_)
  if (!is.finite(post_sd) || post_sd <= 0) post_sd <- 0.5 * sd(mu, na.rm = TRUE)

  list(converged = TRUE, psi_sem = psi_sem, se_sem = se_sem,
       mu = as.numeric(mu), post_sd = post_sd)
}

# ------------------------------------------------------------------------------
# 4. CORRECTED LV-TMLE  -- per-draw fits combined by Rubin's rules
# ------------------------------------------------------------------------------
# For m = 1..M: draw ONE L per subject from its posterior N(mu_i, post_sd),
# forming a complete n-row dataset; run the shared shift-TMLE. Then combine:
#     psi_bar = mean(psi^(m))                       (MI point estimate)
#     W_bar   = mean(var^(m))                        (within-imputation variance)
#     B       = var(psi^(m))                         (between-imputation variance)
#     T       = W_bar + (1 + 1/M) * B                (Rubin total variance)
#     SE      = sqrt(T)
# B is exactly the latent-measurement uncertainty the old SE discarded; W_bar is
# the proper EIF variance from `shift_tmle`. No stacking anywhere.
lv_tmle_mi <- function(dt, mu, post_sd, sl_lib, delta = 1, M = 25, V = 5,
                       gbound = 1e-3, hcap = 30, folds = NULL) {
  n <- nrow(dt)
  X <- as.data.frame(dt[, .(W, Z)])
  psis <- rep(NA_real_, M); vars <- rep(NA_real_, M)
  if (is.null(folds)) {
    V_eff <- max(2L, min(as.integer(V), n))
    folds <- sample(rep(1:V_eff, length.out = n))
  }

  for (m in 1:M) {
    A_m <- rnorm(n, mean = mu, sd = post_sd)        # one draw per subject
    fit <- tryCatch(
      shift_tmle(Y = dt$Y, A = A_m, X = X, sl_lib = sl_lib,
                 delta = delta, V = V, gbound = gbound, hcap = hcap,
                 folds = folds),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$psi) && is.finite(fit$se)) {
      psis[m] <- fit$psi
      vars[m] <- fit$se^2
    }
  }

  ok <- is.finite(psis) & is.finite(vars)
  if (sum(ok) < 2) return(list(psi = NA_real_, se = NA_real_, W_bar = NA, B = NA, df = NA,
                               m_used = sum(ok)))
  psis <- psis[ok]; vars <- vars[ok]; Mu <- length(psis)

  psi_bar <- mean(psis)
  W_bar   <- mean(vars)
  B       <- var(psis)                              # divides by Mu - 1 (Rubin)
  T_var   <- W_bar + (1 + 1 / Mu) * B
  se      <- sqrt(T_var)
  df      <- if (B > 0) (Mu - 1) * (1 + W_bar / ((1 + 1 / Mu) * B))^2 else Inf

  list(psi = psi_bar, se = se, W_bar = W_bar, B = B, df = df, m_used = Mu)
}

# ------------------------------------------------------------------------------
# 5. DRIVER: run all four estimators on one dataset
# ------------------------------------------------------------------------------
# Returns point estimates and SEs for Naive, SEM, Robust (LV-TMLE), and Oracle,
# plus light diagnostics (incl. the LV-TMLE between-imputation variance B, so the
# analysis can report how much of the SE comes from latent uncertainty).
run_comparison <- function(dt, indicators_type, M = 25, V = 5, delta = 1) {
  na_out <- function(reason) list(
    est_naive = NA, se_naive = NA, est_sem = NA, se_sem = NA,
    est_rob = NA, se_rob = NA, est_oracle = NA, se_oracle = NA,
    diagnostics = list(fail_location = "SEM", fail_reason = reason,
                       lv_B = NA, lv_Wbar = NA, lv_m = NA))

  sl_lib <- lv_sl_library(nrow(dt))

  # --- SEM fit + posterior (shared) ---
  post <- fit_sem_posterior(dt, indicators_type)
  if (!isTRUE(post$converged)) return(na_out(post$reason))

  X <- as.data.frame(dt[, .(W, Z)])

  # --- Naive TMLE: Bartlett factor scores as a point mass ---
  folds <- sample(rep(1:max(2L, min(as.integer(V), nrow(dt))), length.out = nrow(dt)))

  naive <- tryCatch(shift_tmle(dt$Y, post$mu, X, sl_lib, delta, V, folds = folds),
                    error = function(e) list(psi = NA, se = NA))

  # --- Oracle TMLE: the TRUE latent (ceiling; uses L_true) ---
  oracle <- tryCatch(shift_tmle(dt$Y, dt$L_true, X, sl_lib, delta, V, folds = folds),
                     error = function(e) list(psi = NA, se = NA))

  # --- Corrected LV-TMLE: MI over the posterior + Rubin SE ---
  lv <- lv_tmle_mi(dt, post$mu, post$post_sd, sl_lib, delta, M, V, folds = folds)

  list(
    est_naive  = naive$psi,  se_naive  = naive$se,
    est_sem    = post$psi_sem, se_sem  = post$se_sem,
    est_rob    = lv$psi,     se_rob    = lv$se,
    est_oracle = oracle$psi, se_oracle = oracle$se,
    diagnostics = list(fail_location = "None", fail_reason = "None",
                       lv_B = lv$B, lv_Wbar = lv$W_bar, lv_m = lv$m_used)
  )
}
