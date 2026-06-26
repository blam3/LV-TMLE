# ==============================================================================
# FILE: 06_outcome_aware.R
# PURPOSE: Outcome-aware latent posterior plus the v2 estimator comparison.
#
# The v1 estimator in 01_setup.R imputes L from indicators only. That is useful as
# an ablation, but it re-introduces attenuation because the imputation noise is
# independent of Y. This file conditions latent draws on the analysis outcome:
#   continuous indicators: closed-form Gaussian L | I, Y, W, Z
#   ordinal indicators   : regression-score approximation, then outcome update
# ==============================================================================

if (!exists("shift_tmle") || !exists("lv_tmle_mi")) source("01_setup.R")

get_pe <- function(pe, lhs, op, rhs, field = "est", default = NA_real_) {
  row <- pe[pe$lhs == lhs & pe$op == op & pe$rhs == rhs, , drop = FALSE]
  if (!nrow(row)) return(default)
  as.numeric(row[[field]][1])
}

lav_var <- function(pe, name, default = NA_real_) {
  get_pe(pe, name, "~~", name, default = default)
}

fit_flexible_y_link <- function(dt, mu_noY, fallback_slope, fallback_var) {
  dat <- data.frame(
    Y = dt$Y,
    L_sim = mu_noY,
    W = dt$W,
    Z = dt$Z,
    Z_sq = dt$Z_sq
  )

  fit <- try(earth(Y ~ ., data = dat, degree = 2), silent = TRUE)
  if (inherits(fit, "try-error")) {
    fit <- try(lm(Y ~ L_sim + W + Z + Z_sq, data = dat), silent = TRUE)
  }

  if (inherits(fit, "try-error")) {
    pred0 <- rep(mean(dt$Y, na.rm = TRUE), nrow(dt))
    slope <- rep(fallback_slope, nrow(dt))
    var_eff <- fallback_var
  } else {
    pred_at <- function(a) {
      nd <- dat
      nd$L_sim <- a
      as.numeric(predict(fit, newdata = nd))
    }
    pred0 <- pred_at(mu_noY)
    h <- max(1e-3, sd(mu_noY, na.rm = TRUE) * 1e-4)
    slope <- (pred_at(mu_noY + h) - pred_at(mu_noY - h)) / (2 * h)
    resid <- dt$Y - pred0
    var_eff <- mean(resid^2, na.rm = TRUE)
  }

  if (!is.finite(var_eff) || var_eff <= 0) var_eff <- fallback_var
  if (!is.finite(var_eff) || var_eff <= 0) var_eff <- 1

  if (any(!is.finite(slope))) slope[!is.finite(slope)] <- fallback_slope
  if (!is.finite(fallback_slope)) fallback_slope <- 0
  slope[!is.finite(slope)] <- fallback_slope

  # Derivative spikes from flexible basis models can make the Gaussian update
  # unstable. Winsorize only the extremes; the subject-level slope pattern remains.
  qs <- quantile(slope, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
  if (all(is.finite(qs)) && qs[1] < qs[2]) {
    slope <- pmin(pmax(slope, qs[1]), qs[2])
  }
  slope <- pmin(pmax(slope, -10), 10)

  if (!exists("pred_at")) {
    pred_at <- function(a) pred0 + slope * (a - mu_noY)
  }

  list(pred0 = pred0, slope = slope, var_eff = var_eff, pred_at = pred_at)
}

flexible_posterior_moments <- function(dt, mu_noY, post_sd_noY, y_link,
                                       grid_width = 5, grid_n = 81) {
  n <- nrow(dt)
  z <- seq(-grid_width, grid_width, length.out = grid_n)
  L_grid <- outer(mu_noY, z * post_sd_noY, "+")

  pred_vec <- vapply(seq_len(grid_n), function(k) y_link$pred_at(L_grid[, k]),
                     numeric(n))
  if (any(!is.finite(pred_vec))) {
    return(NULL)
  }

  logw <- dnorm(L_grid, mean = mu_noY, sd = post_sd_noY, log = TRUE) +
    dnorm(dt$Y, mean = pred_vec, sd = sqrt(y_link$var_eff), log = TRUE)
  row_max <- apply(logw, 1, max)
  w <- exp(logw - row_max)
  w_sum <- rowSums(w)
  ok <- is.finite(w_sum) & w_sum > 0
  if (!all(ok)) return(NULL)

  mu <- rowSums(w * L_grid) / w_sum
  second <- rowSums(w * L_grid^2) / w_sum
  var <- pmax(second - mu^2, 1e-8)

  list(mu = as.numeric(mu), sd = sqrt(var))
}

sem_posterior_full <- function(dt, indicators_type) {
  dt <- copy(dt)
  dt[, Z_sq := Z^2]
  ordered_vars <- if (indicators_type == "ordinal") c("I1", "I2", "I3") else NULL

  sem_model <- '
    L =~ I1 + I2 + I3
    L ~ W + Z + Z_sq
    Y ~ L + W + Z + Z_sq
  '

  fit <- try(suppressWarnings(
    sem(sem_model, data = as.data.frame(dt), ordered = ordered_vars,
        std.lv = FALSE, auto.fix.first = TRUE)
  ), silent = TRUE)

  if (inherits(fit, "try-error") || !lavInspect(fit, "converged"))
    return(list(converged = FALSE, reason = "SEM_no_converge"))

  pe <- parameterEstimates(fit)
  psi_sem <- get_pe(pe, "Y", "~", "L", "est")
  se_sem  <- get_pe(pe, "Y", "~", "L", "se")
  if (!is.finite(psi_sem)) return(list(converged = FALSE, reason = "SEM_coef_NA"))

  mu_noY <- try(lavPredict(fit, method = "regression")[, "L"], silent = TRUE)
  if (inherits(mu_noY, "try-error") || any(!is.finite(mu_noY)))
    mu_noY <- try(lavPredict(fit, method = "bartlett")[, "L"], silent = TRUE)
  if (inherits(mu_noY, "try-error") || any(!is.finite(mu_noY)))
    return(list(converged = FALSE, reason = "Factor_scores_NA"))
  mu_noY <- as.numeric(mu_noY)

  b_l_w   <- get_pe(pe, "L", "~", "W", default = 0)
  b_l_z   <- get_pe(pe, "L", "~", "Z", default = 0)
  b_l_z2  <- get_pe(pe, "L", "~", "Z_sq", default = 0)
  b_y_0   <- get_pe(pe, "Y", "~1", "", default = 0)
  b_y_l   <- get_pe(pe, "Y", "~", "L", default = NA_real_)
  b_y_w   <- get_pe(pe, "Y", "~", "W", default = 0)
  b_y_z   <- get_pe(pe, "Y", "~", "Z", default = 0)
  b_y_z2  <- get_pe(pe, "Y", "~", "Z_sq", default = 0)
  var_l   <- pmax(lav_var(pe, "L"), 1e-6)
  var_y   <- pmax(lav_var(pe, "Y"), 1e-6)

  if (!is.finite(b_y_l) || !is.finite(var_l) || !is.finite(var_y))
    return(list(converged = FALSE, reason = "Posterior_param_NA"))

  prior_mu <- b_l_w * dt$W + b_l_z * dt$Z + b_l_z2 * dt$Z_sq
  y_xpart  <- b_y_0 + b_y_w * dt$W + b_y_z * dt$Z + b_y_z2 * dt$Z_sq

  # Indicator-only posterior SD. For continuous indicators this is exact under
  # the fitted Gaussian SEM; for ordinal indicators it is an approximation on the
  # latent-response scale.
  post_sd_noY <- tryCatch({
    em  <- lavInspect(fit, "est")
    lam <- as.numeric(em$lambda[c("I1", "I2", "I3"), "L"])
    th  <- pmax(diag(em$theta)[c("I1", "I2", "I3")], 1e-6)
    sqrt(1 / (1 / var_l + sum(lam^2 / th)))
  }, error = function(e) NA_real_)
  if (!is.finite(post_sd_noY) || post_sd_noY <= 0)
    post_sd_noY <- 0.5 * sd(mu_noY, na.rm = TRUE)
  if (!is.finite(post_sd_noY) || post_sd_noY <= 0) post_sd_noY <- 0.5

  y_link <- fit_flexible_y_link(dt, mu_noY, fallback_slope = b_y_l,
                                fallback_var = var_y)
  post_var_noY <- post_sd_noY^2
  slope <- y_link$slope
  pred0 <- y_link$pred0
  raw_var_y_eff <- y_link$var_eff
  latent_resid_var <- mean(slope^2, na.rm = TRUE) * post_var_noY
  var_y_eff <- raw_var_y_eff - latent_resid_var
  if (exists("Y_SD")) var_y_eff <- max(var_y_eff, Y_SD^2)
  if (!is.finite(var_y_eff) || var_y_eff <= 0) var_y_eff <- raw_var_y_eff
  var_y_eff <- pmax(var_y_eff, 1e-6)
  y_link$var_eff <- var_y_eff

  # Nonlinear outcome-aware update by quadrature:
  #   p(L | I,Y,W,Z) ∝ p(L | I,W,Z) * N(Y; q_flex(L,W,Z), var_y_eff)
  # This keeps the quadratic/interaction information in Y instead of compressing
  # it into a single misspecified SEM slope.
  post_quad <- flexible_posterior_moments(dt, mu_noY, post_sd_noY, y_link)
  if (is.null(post_quad)) {
    prec <- 1 / post_var_noY + slope^2 / var_y_eff
    post_var <- 1 / prec
    rhs <- mu_noY / post_var_noY + slope * (dt$Y - pred0 + slope * mu_noY) / var_y_eff
    mu_full <- post_var * rhs
    post_sd_full <- sqrt(post_var)
  } else {
    mu_full <- post_quad$mu
    post_sd_full <- post_quad$sd
  }

  if (any(!is.finite(mu_full)) || any(!is.finite(post_sd_full)))
    return(list(converged = FALSE, reason = "Posterior_draw_param_NA"))

  list(
    converged = TRUE,
    psi_sem = psi_sem,
    se_sem = se_sem,
    mu_noY = mu_noY,
    post_sd_noY = rep(post_sd_noY, nrow(dt)),
    mu_full = as.numeric(mu_full),
    post_sd_full = as.numeric(post_sd_full),
    diagnostics = list(
      var_y_sem = as.numeric(var_y),
      var_y_flex_raw = as.numeric(raw_var_y_eff),
      var_y_latent_correction = as.numeric(latent_resid_var),
      var_y_eff = as.numeric(var_y_eff),
      mean_slope_eff = mean(slope, na.rm = TRUE),
      sd_slope_eff = sd(slope, na.rm = TRUE),
      mean_shift_eff = mean(y_link$pred_at(mu_noY + 1) - pred0, na.rm = TRUE),
      mean_mu_shift = mean(mu_full - mu_noY, na.rm = TRUE),
      sd_mu_noY = sd(mu_noY, na.rm = TRUE),
      sd_mu_full = sd(mu_full, na.rm = TRUE)
    )
  )
}

run_comparison_v2 <- function(dt, indicators_type, M = 25, V = 5, delta = 1,
                              include_noY = FALSE, draw_scale = 0.5) {
  na_out <- function(reason, location = "SEM") list(
    est_sem = NA, se_sem = NA,
    est_regcal = NA, se_regcal = NA,
    est_lv = NA, se_lv = NA,
    est_oracle = NA, se_oracle = NA,
    est_lv_noY = NA, se_lv_noY = NA,
    diagnostics = list(fail_location = location, fail_reason = reason,
                       lv_B = NA, lv_Wbar = NA, lv_m = NA,
                       lv_noY_B = NA, lv_noY_Wbar = NA, lv_noY_m = NA,
                       var_y_sem = NA, var_y_eff = NA, mean_slope_eff = NA,
                       sd_slope_eff = NA, mean_shift_eff = NA, mean_mu_shift = NA,
                       sd_mu_noY = NA, sd_mu_full = NA, draw_scale = NA)
  )

  sl_lib <- if (nrow(dt) <= 100)
    c("SL.glm", "SL.mean", "SL.rpart", "SL.earth")
  else
    c("SL.glm", "SL.earth", "SL.gam", "SL.mean", "SL.rpart")

  post <- sem_posterior_full(dt, indicators_type)
  if (!isTRUE(post$converged)) return(na_out(post$reason))

  X <- as.data.frame(dt[, .(W, Z)])
  folds <- sample(rep(1:max(2L, min(as.integer(V), nrow(dt))), length.out = nrow(dt)))

  regcal <- tryCatch(
    shift_tmle(dt$Y, post$mu_noY, X, sl_lib, delta = delta, V = V, folds = folds),
    error = function(e) list(psi = NA_real_, se = NA_real_)
  )
  oracle <- tryCatch(
    shift_tmle(dt$Y, dt$L_true, X, sl_lib, delta = delta, V = V, folds = folds),
    error = function(e) list(psi = NA_real_, se = NA_real_)
  )
  lv <- lv_tmle_mi(dt, post$mu_full, post$post_sd_full * draw_scale, sl_lib,
                   delta = delta, M = M, V = V, folds = folds)

  lv_noY <- if (isTRUE(include_noY)) {
    lv_tmle_mi(dt, post$mu_noY, post$post_sd_noY, sl_lib,
               delta = delta, M = M, V = V, folds = folds)
  } else {
    list(psi = NA_real_, se = NA_real_, W_bar = NA, B = NA, m_used = NA)
  }

  list(
    est_sem = post$psi_sem, se_sem = post$se_sem,
    est_regcal = regcal$psi, se_regcal = regcal$se,
    est_lv = lv$psi, se_lv = lv$se,
    est_oracle = oracle$psi, se_oracle = oracle$se,
    est_lv_noY = lv_noY$psi, se_lv_noY = lv_noY$se,
    diagnostics = list(fail_location = "None", fail_reason = "None",
                       lv_B = lv$B, lv_Wbar = lv$W_bar, lv_m = lv$m_used,
                       lv_noY_B = lv_noY$B, lv_noY_Wbar = lv_noY$W_bar,
                       lv_noY_m = lv_noY$m_used,
                       var_y_sem = post$diagnostics$var_y_sem,
                       var_y_eff = post$diagnostics$var_y_eff,
                       mean_slope_eff = post$diagnostics$mean_slope_eff,
                       sd_slope_eff = post$diagnostics$sd_slope_eff,
                       mean_shift_eff = post$diagnostics$mean_shift_eff,
                       mean_mu_shift = post$diagnostics$mean_mu_shift,
                       sd_mu_noY = post$diagnostics$sd_mu_noY,
                       sd_mu_full = post$diagnostics$sd_mu_full,
                       draw_scale = draw_scale)
  )
}
