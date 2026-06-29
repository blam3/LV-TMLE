# ==============================================================================
# FILE: 05d_posterior_gap_diagnostics.R
# PURPOSE: Diagnose the outcome-aware posterior/congeniality gap in
#          sem_posterior_full() before scaling run_comparison_v2().
#
# This script does NOT propose a production estimator. It compares the current
# outcome-aware posterior against oracle-only falsifiers to locate the bias source:
#   1. Does mu_full improve recovery of L_true over mu_noY?
#   2. Is the fitted outcome link's one-unit shift too low?
#   3. Does the residual-variance correction materially change the
#      posterior update?
#
# Fast run:
#   Rscript --vanilla -e 'options(lv_tmle.diag_reps=20); source("05d_posterior_gap_diagnostics.R")'
#
# Multi-scenario run:
#   Rscript --vanilla -e 'options(lv_tmle.diag_reps=20, lv_tmle.diag_scenarios=c(14,17,32), lv_tmle.diag_output="posterior_gap_diagnostics_scen_14_17_32.rds"); source("05d_posterior_gap_diagnostics.R")'
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lavaan)
  library(earth)
})

source("00_dgp_variants.R")
source("01_setup.R")
source("06_outcome_aware.R")

REPS <- as.integer(getOption("lv_tmle.diag_reps", 30))
if (!is.finite(REPS) || REPS <= 0) REPS <- 30
SEED <- as.integer(getOption("lv_tmle.diag_seed", 20260627))
if (!is.finite(SEED)) SEED <- 20260627
OUTPUT_FILE <- as.character(getOption("lv_tmle.diag_output",
                                      "posterior_gap_diagnostics.rds"))
if (length(OUTPUT_FILE) != 1L || !nzchar(OUTPUT_FILE)) {
  OUTPUT_FILE <- "posterior_gap_diagnostics.rds"
}

SCEN_OUTCOME    <- "nonlinear"
SCEN_LATENT     <- "normal"
SCEN_INDICATORS <- "continuous"
SCEN_N          <- 250
SCEN_ERROR      <- "medium"
DELTA           <- 1

requested_scenarios <- getOption("lv_tmle.diag_scenarios", NULL)
if (!is.null(requested_scenarios)) {
  requested_scenarios <- as.integer(requested_scenarios)
  requested_scenarios <- requested_scenarios[is.finite(requested_scenarios)]
}

if (length(requested_scenarios)) {
  diag_grid <- sim_grid[sim_grid$Scenario_ID %in% requested_scenarios, ]
  if (nrow(diag_grid) != length(unique(requested_scenarios))) {
    missing_ids <- setdiff(unique(requested_scenarios), diag_grid$Scenario_ID)
    stop(sprintf("Unknown diagnostic Scenario_ID(s): %s",
                 paste(missing_ids, collapse = ", ")))
  }
  diag_grid <- diag_grid[match(unique(requested_scenarios),
                               diag_grid$Scenario_ID), ]
} else {
  diag_grid <- data.frame(
    Scenario_ID = NA_integer_,
    outcome = SCEN_OUTCOME,
    latent = SCEN_LATENT,
    n_size = SCEN_N,
    meas_error = SCEN_ERROR,
    indicators = SCEN_INDICATORS,
    stringsAsFactors = FALSE
  )
}

scenario_cols <- c("Scenario_ID", "n_size", "meas_error", "outcome",
                   "latent", "indicators")

scenario_meta <- function(scen) {
  out <- as.data.frame(scen[, scenario_cols, drop = FALSE])
  rownames(out) <- NULL
  out
}

calibration_metrics <- function(mu, post_sd, L_true, label) {
  ok <- is.finite(mu) & is.finite(L_true)
  mu <- mu[ok]; L_true <- L_true[ok]
  psd <- if (length(post_sd) == 1) rep(post_sd, length(mu)) else post_sd[ok]
  z <- (L_true - mu) / psd
  fit <- lm(L_true ~ mu)
  data.frame(
    posterior = label,
    corr = cor(mu, L_true),
    rmse = sqrt(mean((mu - L_true)^2)),
    bias = mean(mu - L_true),
    sd_mu = sd(mu),
    sd_L = sd(L_true),
    cal_intercept = unname(coef(fit)[1]),
    cal_slope = unname(coef(fit)[2]),
    mean_post_sd = mean(psd),
    z_mean = mean(z),
    z_sd = sd(z),
    stringsAsFactors = FALSE
  )
}

make_current_link <- function(dt, post) {
  y_link <- fit_me_corrected_y_link(
    dt,
    post$mu_noY,
    post$post_sd_noY,
    fallback_slope = post$psi_sem,
    fallback_var = post$diagnostics$var_y_sem
  )
  raw_var <- y_link$var_eff
  latent_correction <- if (!is.null(y_link$latent_resid_var) &&
                           is.finite(y_link$latent_resid_var)) {
    y_link$latent_resid_var
  } else {
    mean(y_link$slope^2, na.rm = TRUE) * mean(post$post_sd_noY^2, na.rm = TRUE)
  }
  corrected_var <- raw_var - latent_correction
  if (exists("Y_SD")) corrected_var <- max(corrected_var, Y_SD^2)
  if (!is.finite(corrected_var) || corrected_var <= 0) corrected_var <- raw_var
  y_link$var_eff <- pmax(corrected_var, 1e-6)
  y_link
}

make_me_rawvar_link <- function(dt, post) {
  fit_me_corrected_y_link(
    dt,
    post$mu_noY,
    post$post_sd_noY,
    fallback_slope = post$psi_sem,
    fallback_var = post$diagnostics$var_y_sem
  )
}

make_rawvar_link <- function(dt, post) {
  fit_flexible_y_link(
    dt,
    post$mu_noY,
    fallback_slope = post$psi_sem,
    fallback_var = post$diagnostics$var_y_sem
  )
}

me_link_coefficients <- function(dt, mu_noY, post_sd_noY) {
  post_var <- if (length(post_sd_noY) == 1) {
    rep(post_sd_noY^2, nrow(dt))
  } else {
    post_sd_noY^2
  }
  dat <- data.frame(
    Y = dt$Y,
    L_m1 = mu_noY,
    L_m2 = mu_noY^2 + post_var,
    L_mW = mu_noY * dt$W,
    W = dt$W,
    Z = dt$Z,
    Z_sq = dt$Z_sq
  )
  fit <- try(lm(Y ~ L_m1 + L_m2 + L_mW + W + Z + Z_sq, data = dat),
             silent = TRUE)
  if (inherits(fit, "try-error") || any(!is.finite(coef(fit)))) return(NULL)
  cf <- coef(fit)
  cf[is.na(cf)] <- 0
  get_cf <- function(nm) if (nm %in% names(cf)) unname(cf[[nm]]) else 0
  c(
    b0 = get_cf("(Intercept)"),
    b1 = get_cf("L_m1"),
    b2 = get_cf("L_m2"),
    b3 = get_cf("L_mW"),
    bW = get_cf("W"),
    bZ = get_cf("Z"),
    bZ2 = get_cf("Z_sq")
  )
}

make_crossfit_me_link <- function(dt, post, V = 5, correct_var = TRUE) {
  n <- nrow(dt)
  V <- max(2L, min(as.integer(V), n))
  folds <- rep(seq_len(V), length.out = n)
  post_sd <- mean(post$post_sd_noY, na.rm = TRUE)
  post_var <- post_sd^2
  b0 <- b1 <- b2 <- b3 <- bW <- bZ <- bZ2 <- rep(NA_real_, n)

  for (v in seq_len(V)) {
    tr <- folds != v
    va <- folds == v
    cf <- me_link_coefficients(dt[tr], post$mu_noY[tr], post_sd)
    if (is.null(cf)) next
    b0[va] <- cf["b0"]; b1[va] <- cf["b1"]; b2[va] <- cf["b2"]
    b3[va] <- cf["b3"]; bW[va] <- cf["bW"]; bZ[va] <- cf["bZ"]
    bZ2[va] <- cf["bZ2"]
  }

  missing <- !is.finite(b0 + b1 + b2 + b3 + bW + bZ + bZ2)
  if (any(missing)) {
    fallback <- me_link_coefficients(dt, post$mu_noY, post_sd)
    if (is.null(fallback)) {
      return(make_current_link(dt, post))
    }
    b0[missing] <- fallback["b0"]; b1[missing] <- fallback["b1"]
    b2[missing] <- fallback["b2"]; b3[missing] <- fallback["b3"]
    bW[missing] <- fallback["bW"]; bZ[missing] <- fallback["bZ"]
    bZ2[missing] <- fallback["bZ2"]
  }

  pred_at <- function(a) {
    b0 + b1 * a + b2 * a^2 + b3 * a * dt$W +
      bW * dt$W + bZ * dt$Z + bZ2 * dt$Z_sq
  }
  pred0 <- pred_at(post$mu_noY)
  pred_obs <- b0 + b1 * post$mu_noY + b2 * (post$mu_noY^2 + post_var) +
    b3 * post$mu_noY * dt$W + bW * dt$W + bZ * dt$Z + bZ2 * dt$Z_sq
  slope <- b1 + 2 * b2 * post$mu_noY + b3 * dt$W
  raw_var <- mean((dt$Y - pred_obs)^2, na.rm = TRUE)
  if (!is.finite(raw_var) || raw_var <= 0) raw_var <- post$diagnostics$var_y_sem
  if (!is.finite(raw_var) || raw_var <= 0) raw_var <- 1

  c_lin <- b1 + b3 * dt$W
  latent_var <- c_lin^2 * post_var +
    4 * c_lin * b2 * post$mu_noY * post_var +
    b2^2 * (2 * post_var^2 + 4 * post$mu_noY^2 * post_var)
  latent_var <- pmax(latent_var, 0)

  if (any(!is.finite(slope))) slope[!is.finite(slope)] <- post$psi_sem
  if (!is.finite(post$psi_sem)) slope[!is.finite(slope)] <- 0
  slope[!is.finite(slope)] <- post$psi_sem
  qs <- quantile(slope, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  if (all(is.finite(qs)) && qs[1] < qs[2]) {
    slope <- pmin(pmax(slope, qs[1]), qs[2])
  }
  slope <- pmin(pmax(slope, -12), 12)

  var_eff <- raw_var
  if (isTRUE(correct_var)) {
    var_eff <- raw_var - mean(latent_var, na.rm = TRUE)
    if (exists("Y_SD")) var_eff <- max(var_eff, Y_SD^2)
    if (!is.finite(var_eff) || var_eff <= 0) var_eff <- raw_var
  }

  list(
    pred0 = pred0,
    slope = slope,
    var_eff = pmax(var_eff, 1e-6),
    pred_at = pred_at,
    latent_resid_var = mean(latent_var, na.rm = TRUE),
    link_type = if (isTRUE(correct_var)) "cf_me_quadratic" else "cf_me_raw"
  )
}

temper_y_link <- function(link, power) {
  out <- link
  if (!is.finite(power) || power <= 0) stop("power must be positive")
  out$var_eff <- out$var_eff / power
  out$link_type <- paste0(ifelse(is.null(out$link_type), "link", out$link_type),
                          "_power_", power)
  out
}

make_score_quadratic_link <- function(dt, post) {
  dat <- data.frame(
    Y = dt$Y,
    L_sim = post$mu_noY,
    W = dt$W,
    Z = dt$Z,
    Z_sq = dt$Z_sq
  )
  fit <- lm(Y ~ L_sim + I(L_sim^2) + L_sim:W + W + Z + Z_sq, data = dat)
  pred_at <- function(a) {
    nd <- dat
    nd$L_sim <- a
    as.numeric(predict(fit, newdata = nd))
  }
  pred0 <- pred_at(post$mu_noY)
  h <- max(1e-3, sd(post$mu_noY, na.rm = TRUE) * 1e-4)
  slope <- (pred_at(post$mu_noY + h) - pred_at(post$mu_noY - h)) / (2 * h)
  var_eff <- mean((dt$Y - pred0)^2, na.rm = TRUE)
  if (!is.finite(var_eff) || var_eff <= 0) {
    var_eff <- post$diagnostics$var_y_sem
  }
  list(pred0 = pred0, slope = slope, var_eff = var_eff, pred_at = pred_at)
}

correct_link_var <- function(link, post, floor_var = if (exists("Y_SD")) Y_SD^2 else 1) {
  out <- link
  latent_correction <- mean(out$slope^2, na.rm = TRUE) *
    mean(post$post_sd_noY^2, na.rm = TRUE)
  var_eff <- out$var_eff - latent_correction
  if (is.finite(floor_var)) var_eff <- max(var_eff, floor_var)
  if (!is.finite(var_eff) || var_eff <= 0) var_eff <- out$var_eff
  out$var_eff <- pmax(var_eff, 1e-6)
  out
}

make_oracle_true_link <- function(dt, mu_noY, outcome) {
  outcome_fun <- OUTCOME_FUNCS[[outcome]]
  if (is.null(outcome_fun)) stop(sprintf("Unknown outcome type: %s", outcome))
  pred_at <- function(a) outcome_fun(a, dt$W, dt$Z)
  pred0 <- pred_at(mu_noY)
  h <- max(1e-3, sd(mu_noY, na.rm = TRUE) * 1e-4)
  slope <- (pred_at(mu_noY + h) - pred_at(mu_noY - h)) / (2 * h)
  list(
    pred0 = pred0,
    slope = slope,
    var_eff = Y_SD^2,
    pred_at = pred_at
  )
}

make_oracle_earth_link <- function(dt, mu_noY) {
  dat <- data.frame(
    Y = dt$Y,
    L_sim = dt$L_true,
    W = dt$W,
    Z = dt$Z,
    Z_sq = dt$Z_sq
  )
  fit <- earth(Y ~ ., data = dat, degree = 2)
  pred_at <- function(a) {
    nd <- dat
    nd$L_sim <- a
    as.numeric(predict(fit, newdata = nd))
  }
  pred0 <- pred_at(mu_noY)
  h <- max(1e-3, sd(mu_noY, na.rm = TRUE) * 1e-4)
  slope <- (pred_at(mu_noY + h) - pred_at(mu_noY - h)) / (2 * h)
  var_eff <- mean((dt$Y - pred_at(dt$L_true))^2, na.rm = TRUE)
  if (!is.finite(var_eff) || var_eff <= 0) var_eff <- Y_SD^2
  list(pred0 = pred0, slope = slope, var_eff = var_eff, pred_at = pred_at)
}

posterior_from_link <- function(dt, post, link, label) {
  post_sd_noY <- mean(post$post_sd_noY, na.rm = TRUE)
  pm <- flexible_posterior_moments(dt, post$mu_noY, post_sd_noY, link)
  if (is.null(pm)) {
    return(list(metrics = NULL, link = NULL))
  }
  metrics <- calibration_metrics(pm$mu, pm$sd, dt$L_true, label)
  move <- pm$mu - post$mu_noY
  latent_error <- dt$L_true - post$mu_noY
  metrics$mean_move <- mean(move)
  metrics$move_error_corr <- cor(move, latent_error)
  metrics$mean_shift_eff <- mean(link$pred_at(post$mu_noY + DELTA) - link$pred0)
  metrics$mean_slope <- mean(link$slope, na.rm = TRUE)
  metrics$sd_slope <- sd(link$slope, na.rm = TRUE)
  metrics$var_eff <- link$var_eff
  list(metrics = metrics, link = link)
}

set.seed(SEED)
cat(sprintf("\nPosterior gap diagnostics: reps=%d scenario_count=%d seed=%d\n",
            REPS, nrow(diag_grid), SEED))
cat(sprintf("Writing diagnostics to %s\n", OUTPUT_FILE))

rows <- list()
link_rows <- list()
failures <- list()

for (s in seq_len(nrow(diag_grid))) {
  scen <- diag_grid[s, , drop = FALSE]
  meta <- scenario_meta(scen)
  scen_label <- if (is.finite(scen$Scenario_ID)) {
    sprintf("Scenario %d", scen$Scenario_ID)
  } else {
    "Default diagnostic scenario"
  }
  cat(sprintf(
    "\n%s: n=%d error=%s outcome=%s latent=%s indicators=%s\n",
    scen_label, scen$n_size, scen$meas_error, scen$outcome, scen$latent,
    scen$indicators
  ))

  for (b in seq_len(REPS)) {
    dt <- generate_data(scen$n_size, scen$meas_error, scen$outcome,
                        scen$latent, scen$indicators)
    dt[, Z_sq := Z^2]
    post <- sem_posterior_full(dt, scen$indicators)
    if (!isTRUE(post$converged)) {
      failures[[length(failures) + 1]] <- cbind(
        meta, data.frame(rep = b, reason = post$reason,
                         stringsAsFactors = FALSE)
      )
      next
    }

    base <- calibration_metrics(post$mu_noY, post$post_sd_noY, dt$L_true, "noY")
    base$mean_move <- 0
    base$move_error_corr <- NA_real_
    base$mean_shift_eff <- NA_real_
    base$mean_slope <- NA_real_
    base$sd_slope <- NA_real_
    base$var_eff <- NA_real_
    rows[[length(rows) + 1]] <- cbind(meta, rep = b, base)

    current <- calibration_metrics(post$mu_full, post$post_sd_full, dt$L_true,
                                   "current_full")
    current$mean_move <- mean(post$mu_full - post$mu_noY)
    current$move_error_corr <- cor(post$mu_full - post$mu_noY,
                                   dt$L_true - post$mu_noY)
    current$mean_shift_eff <- post$diagnostics$mean_shift_eff
    current$mean_slope <- post$diagnostics$mean_slope_eff
    current$sd_slope <- post$diagnostics$sd_slope_eff
    current$var_eff <- post$diagnostics$var_y_eff
    rows[[length(rows) + 1]] <- cbind(meta, rep = b, current)

    raw <- posterior_from_link(dt, post, make_rawvar_link(dt, post), "raw_var")
    if (!is.null(raw$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, raw$metrics)
    }

    me_raw <- posterior_from_link(dt, post, make_me_rawvar_link(dt, post),
                                  "me_raw_var")
    if (!is.null(me_raw$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, me_raw$metrics)
    }

    current_half <- posterior_from_link(
      dt, post, temper_y_link(make_current_link(dt, post), 0.5),
      "current_y_power_0.5"
    )
    if (!is.null(current_half$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, current_half$metrics)
    }

    current_quarter <- posterior_from_link(
      dt, post, temper_y_link(make_current_link(dt, post), 0.25),
      "current_y_power_0.25"
    )
    if (!is.null(current_quarter$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b,
                                        current_quarter$metrics)
    }

    me_raw_half <- posterior_from_link(
      dt, post, temper_y_link(make_me_rawvar_link(dt, post), 0.5),
      "me_raw_y_power_0.5"
    )
    if (!is.null(me_raw_half$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, me_raw_half$metrics)
    }

    cf_me <- posterior_from_link(
      dt, post, make_crossfit_me_link(dt, post, correct_var = TRUE),
      "cf_me_corr"
    )
    if (!is.null(cf_me$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, cf_me$metrics)
    }

    cf_me_raw <- posterior_from_link(
      dt, post, make_crossfit_me_link(dt, post, correct_var = FALSE),
      "cf_me_raw"
    )
    if (!is.null(cf_me_raw$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, cf_me_raw$metrics)
    }

    score_quad <- posterior_from_link(
      dt, post, make_score_quadratic_link(dt, post), "score_quadratic"
    )
    if (!is.null(score_quad$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, score_quad$metrics)
    }

    score_quad_corr <- posterior_from_link(
      dt, post,
      correct_link_var(make_score_quadratic_link(dt, post), post),
      "score_quad_corr"
    )
    if (!is.null(score_quad_corr$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b,
                                        score_quad_corr$metrics)
    }

    score_quad_y1 <- make_score_quadratic_link(dt, post)
    score_quad_y1$var_eff <- Y_SD^2
    score_quad_y1 <- posterior_from_link(dt, post, score_quad_y1,
                                         "score_quad_y1")
    if (!is.null(score_quad_y1$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b,
                                        score_quad_y1$metrics)
    }

    oracle_true <- posterior_from_link(
      dt, post, make_oracle_true_link(dt, post$mu_noY, scen$outcome),
      "oracle_true_q"
    )
    if (!is.null(oracle_true$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, oracle_true$metrics)
    }

    oracle_earth <- posterior_from_link(
      dt, post, make_oracle_earth_link(dt, post$mu_noY), "oracle_earth_trueL"
    )
    if (!is.null(oracle_earth$metrics)) {
      rows[[length(rows) + 1]] <- cbind(meta, rep = b, oracle_earth$metrics)
    }

    link_current <- make_current_link(dt, post)
    link_raw <- make_rawvar_link(dt, post)
    link_true <- make_oracle_true_link(dt, post$mu_noY, scen$outcome)
    link_rows[[length(link_rows) + 1]] <- cbind(
      meta,
      data.frame(
        rep = b,
        raw_var = link_raw$var_eff,
        current_var = link_current$var_eff,
        latent_correction = link_raw$var_eff - link_current$var_eff,
        current_shift = mean(link_current$pred_at(post$mu_noY + DELTA) -
                               link_current$pred0),
        true_shift_at_mu_noY = mean(link_true$pred_at(post$mu_noY + DELTA) -
                                      link_true$pred0),
        true_shift_at_L = mean(link_true$pred_at(dt$L_true + DELTA) -
                                 link_true$pred_at(dt$L_true)),
        stringsAsFactors = FALSE
      )
    )

    if (b %% 5 == 0 || b == 1) cat(sprintf("  rep %d/%d\n", b, REPS))
  }
}

diag <- rbindlist(rows, fill = TRUE)
links <- rbindlist(link_rows, fill = TRUE)
failure_dt <- rbindlist(failures, fill = TRUE)
saveRDS(list(posterior = diag, links = links, failures = failure_dt),
        OUTPUT_FILE)

summ_cols <- c("corr", "rmse", "bias", "sd_mu", "cal_slope", "mean_post_sd",
               "z_mean", "z_sd", "mean_move", "move_error_corr",
               "mean_shift_eff", "mean_slope", "sd_slope", "var_eff")
summary_tab <- diag[, lapply(.SD, mean, na.rm = TRUE),
                    by = c(scenario_cols, "posterior"), .SDcols = summ_cols]
summary_tab[, order_id := match(posterior, c("noY", "current_full", "raw_var",
                                             "me_raw_var",
                                             "current_y_power_0.5",
                                             "current_y_power_0.25",
                                             "me_raw_y_power_0.5",
                                             "cf_me_corr",
                                             "cf_me_raw",
                                             "score_quadratic",
                                             "score_quad_corr",
                                             "score_quad_y1",
                                             "oracle_true_q",
                                             "oracle_earth_trueL"))]
setorder(summary_tab, order_id)
summary_tab[, order_id := NULL]

cat("\n=== Posterior Recovery Summary ===\n")
print(format(summary_tab, digits = 3), row.names = FALSE)

cat("\n=== Outcome-Link Diagnostics ===\n")
link_summary <- links[, lapply(.SD, mean, na.rm = TRUE),
                      by = scenario_cols,
                      .SDcols = setdiff(names(links),
                                        c(scenario_cols, "rep"))]
print(format(link_summary, digits = 3), row.names = FALSE)

cat("\nReading guide:\n")
cat("  - If oracle_true_q improves RMSE/calibration while current_full does not,\n")
cat("    the posterior machinery is capable; the fitted outcome link is the gap.\n")
cat("  - If current_shift is below true_shift_at_mu_noY, the fitted outcome link\n")
cat("    is attenuating the outcome update before TMLE sees it.\n")
cat("  - raw_var vs current_full separates residual-variance weighting from the\n")
cat("    outcome-link congeniality problem.\n")

if (nrow(failure_dt)) {
  cat("\nFailures:\n")
  print(failure_dt[, .N, by = c(scenario_cols, "reason")][order(-N)])
}
