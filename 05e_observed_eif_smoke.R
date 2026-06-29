# ==============================================================================
# FILE: 05e_observed_eif_smoke.R
# PURPOSE: Smoke harness for the finite-dimensional observed-data EIF in
#          OBSERVED_DATA_EIF_DERIVATION.md.
#
# This is diagnostic infrastructure only. It does not change the production
# LV-TMLE estimator. It computes, on one continuous-indicator dataset:
#   1. posterior moments E[L^k | O], k = 1..4, by quadrature;
#   2. projected observed scores S_eta^O = E[S_eta^F | O];
#   3. empirical observed information I_eta;
#   4. direct target derivative dot_psi_eta;
#   5. D_eff = xi_eta(X) - psi + dot_psi_eta' I_eta^{-1} S_eta^O.
#
# Fast run:
#   Rscript --vanilla 05e_observed_eif_smoke.R
# Tuned run:
#   Rscript --vanilla -e 'options(lv_tmle.eif_n=1000, lv_tmle.eif_scenario=14); source("05e_observed_eif_smoke.R")'
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lavaan)
})

source("00_dgp_variants.R")
source("01_setup.R")

N <- as.integer(getOption("lv_tmle.eif_n", 750))
if (!is.finite(N) || N <= 0) N <- 750
SCENARIO_ID <- as.integer(getOption("lv_tmle.eif_scenario", 14))
if (!is.finite(SCENARIO_ID)) SCENARIO_ID <- 14
SEED <- as.integer(getOption("lv_tmle.eif_seed", 20260629))
if (!is.finite(SEED)) SEED <- 20260629
DELTA <- as.numeric(getOption("lv_tmle.eif_delta", 1))
OUT_FILE <- as.character(getOption("lv_tmle.eif_output",
                                   "observed_eif_smoke.rds"))
if (length(OUT_FILE) != 1L || !nzchar(OUT_FILE)) {
  OUT_FILE <- "observed_eif_smoke.rds"
}

scen <- sim_grid[sim_grid$Scenario_ID == SCENARIO_ID, ]
if (nrow(scen) != 1L) stop("Scenario_ID not found exactly once.")
if (scen$indicators != "continuous") {
  stop("05e currently supports continuous indicators only.")
}

get_pe_local <- function(pe, lhs, op, rhs, default = NA_real_) {
  row <- pe[pe$lhs == lhs & pe$op == op & pe$rhs == rhs, , drop = FALSE]
  if (!nrow(row)) return(default)
  as.numeric(row$est[1])
}

fit_measurement_eta <- function(dt) {
  d <- copy(dt)
  d[, Z_sq := Z^2]
  model <- '
    L =~ I1 + I2 + I3
    L ~ W + Z + Z_sq
  '
  fit <- try(with_lavaan_core_fallback(suppressWarnings(
    sem(model, data = as.data.frame(d), meanstructure = TRUE,
        std.lv = FALSE, auto.fix.first = TRUE, ncpus = 1L)
  )), silent = TRUE)
  if (inherits(fit, "try-error") || !lavInspect(fit, "converged")) {
    return(NULL)
  }

  pe <- parameterEstimates(fit)
  lambda <- c(
    get_pe_local(pe, "L", "=~", "I1", default = 1),
    get_pe_local(pe, "L", "=~", "I2"),
    get_pe_local(pe, "L", "=~", "I3")
  )
  theta <- pmax(c(
    get_pe_local(pe, "I1", "~~", "I1"),
    get_pe_local(pe, "I2", "~~", "I2"),
    get_pe_local(pe, "I3", "~~", "I3")
  ), 1e-6)
  nu <- c(
    get_pe_local(pe, "I1", "~1", "", default = 0),
    get_pe_local(pe, "I2", "~1", "", default = 0),
    get_pe_local(pe, "I3", "~1", "", default = 0)
  )
  alpha <- c(
    `(Intercept)` = get_pe_local(pe, "L", "~1", "", default = 0),
    W = get_pe_local(pe, "L", "~", "W", default = 0),
    Z = get_pe_local(pe, "L", "~", "Z", default = 0),
    Z_sq = get_pe_local(pe, "L", "~", "Z_sq", default = 0)
  )
  tau2 <- pmax(get_pe_local(pe, "L", "~~", "L"), 1e-6)
  if (any(!is.finite(lambda)) || any(!is.finite(theta)) ||
      any(!is.finite(nu)) || any(!is.finite(alpha)) ||
      !is.finite(tau2)) {
    return(NULL)
  }

  list(lambda = lambda, theta = theta, nu = nu, alpha = alpha,
       tau2 = tau2, fit = fit)
}

normal_moments4 <- function(mu, sig2) {
  cbind(
    M0 = rep(1, length(mu)),
    M1 = mu,
    M2 = mu^2 + sig2,
    M3 = mu^3 + 3 * mu * sig2,
    M4 = mu^4 + 6 * mu^2 * sig2 + 3 * sig2^2
  )
}

r_expectations <- function(M, W, Z, Z_sq) {
  coeff <- cbind(
    b0 = 1,
    bL = 1,
    bL2 = 1,
    bLW = W,
    bW = W,
    bZ = Z,
    bZ2 = Z_sq
  )
  power <- c(0, 1, 2, 1, 0, 0, 0)
  n <- nrow(M)
  p <- length(power)
  Er <- matrix(NA_real_, n, p)
  colnames(Er) <- colnames(coeff)
  Err <- array(NA_real_, dim = c(n, p, p),
               dimnames = list(NULL, colnames(coeff), colnames(coeff)))
  for (j in seq_len(p)) {
    Er[, j] <- coeff[, j] * M[, power[j] + 1L]
    for (k in seq_len(p)) {
      Err[, j, k] <- coeff[, j] * coeff[, k] *
        M[, power[j] + power[k] + 1L]
    }
  }
  list(Er = Er, Err = Err)
}

quad_form_rows <- function(Err, beta) {
  vapply(seq_len(dim(Err)[1]), function(i) {
    as.numeric(crossprod(beta, Err[i, , ] %*% beta))
  }, numeric(1))
}

indicator_posterior <- function(dt, eta) {
  Xmat <- cbind(1, dt$W, dt$Z, dt$Z_sq)
  prior_mu <- as.numeric(Xmat %*% eta$alpha)
  Icen <- sweep(as.matrix(dt[, .(I1, I2, I3)]), 2, eta$nu, "-")
  prec <- 1 / eta$tau2 + sum(eta$lambda^2 / eta$theta)
  sig2 <- 1 / prec
  num <- prior_mu / eta$tau2 +
    as.numeric(Icen %*% (eta$lambda / eta$theta))
  mu <- sig2 * num
  list(mu = mu, sig2 = rep(sig2, nrow(dt)), prior_mu = prior_mu)
}

fit_quadratic_outcome <- function(dt, post_i) {
  M <- normal_moments4(post_i$mu, post_i$sig2)
  rex <- r_expectations(M, dt$W, dt$Z, dt$Z_sq)
  beta <- as.numeric(qr.solve(rex$Er, dt$Y))
  names(beta) <- colnames(rex$Er)
  fitted_obs <- as.numeric(rex$Er %*% beta)
  resid_second <- dt$Y^2 - 2 * dt$Y * fitted_obs +
    quad_form_rows(rex$Err, beta)
  vy <- mean(resid_second, na.rm = TRUE)
  if (!is.finite(vy) || vy <= 0) vy <- mean((dt$Y - fitted_obs)^2)
  if (!is.finite(vy) || vy <= 0) vy <- 1
  list(beta = beta, vy = pmax(vy, 1e-6),
       fitted_obs = fitted_obs, resid_second = resid_second)
}

q_beta <- function(L, dt, beta) {
  beta["b0"] + beta["bL"] * L + beta["bL2"] * L^2 +
    beta["bLW"] * L * dt$W + beta["bW"] * dt$W +
    beta["bZ"] * dt$Z + beta["bZ2"] * dt$Z_sq
}

posterior_moments_y <- function(dt, post_i, beta, vy,
                                grid_width = 6, grid_n = 121) {
  n <- nrow(dt)
  z <- seq(-grid_width, grid_width, length.out = grid_n)
  L_grid <- matrix(post_i$mu, nrow = n, ncol = grid_n) +
    outer(sqrt(post_i$sig2), z, "*")
  q_mat <- vapply(seq_len(grid_n), function(k) q_beta(L_grid[, k], dt, beta),
                  numeric(n))
  mu_mat <- matrix(post_i$mu, nrow = n, ncol = grid_n)
  sig2_mat <- matrix(post_i$sig2, nrow = n, ncol = grid_n)
  y_mat <- matrix(dt$Y, nrow = n, ncol = grid_n)
  log_prior <- -0.5 * log(2 * pi * sig2_mat) -
    0.5 * (L_grid - mu_mat)^2 / sig2_mat
  log_y <- -0.5 * log(2 * pi * vy) - 0.5 * (y_mat - q_mat)^2 / vy
  logw <- log_prior + log_y
  row_max <- apply(logw, 1, max)
  w <- exp(logw - row_max)
  w_sum <- rowSums(w)
  if (any(!is.finite(w_sum)) || any(w_sum <= 0)) {
    stop("Posterior quadrature produced invalid weights.")
  }
  w <- w / w_sum
  M <- sapply(0:4, function(k) rowSums(w * L_grid^k))
  colnames(M) <- paste0("M", 0:4)
  list(M = M, mean = M[, "M1"], var = pmax(M[, "M2"] - M[, "M1"]^2, 0))
}

projected_scores <- function(dt, eta, beta, vy, post_y) {
  M <- post_y$M
  Xmat <- cbind(1, dt$W, dt$Z, dt$Z_sq)
  colnames(Xmat) <- c("alpha0", "alpha_W", "alpha_Z", "alpha_Z2")
  m <- as.numeric(Xmat %*% eta$alpha)
  M1 <- M[, "M1"]; M2 <- M[, "M2"]

  S_alpha <- Xmat * as.numeric((M1 - m) / eta$tau2)
  S_tau2 <- -1 / (2 * eta$tau2) +
    (M2 - 2 * m * M1 + m^2) / (2 * eta$tau2^2)

  I_mat <- as.matrix(dt[, .(I1, I2, I3)])
  S_nu <- matrix(NA_real_, nrow(dt), 3)
  S_lambda <- matrix(NA_real_, nrow(dt), 2)
  S_theta <- matrix(NA_real_, nrow(dt), 3)
  colnames(S_nu) <- paste0("nu", 1:3)
  colnames(S_lambda) <- paste0("lambda", 2:3)
  colnames(S_theta) <- paste0("theta", 1:3)
  for (j in 1:3) {
    I_center <- I_mat[, j] - eta$nu[j]
    S_nu[, j] <- (I_center - eta$lambda[j] * M1) / eta$theta[j]
    e2 <- I_center^2 - 2 * eta$lambda[j] * I_center * M1 +
      eta$lambda[j]^2 * M2
    S_theta[, j] <- -1 / (2 * eta$theta[j]) +
      e2 / (2 * eta$theta[j]^2)
    if (j >= 2) {
      S_lambda[, j - 1] <- (I_center * M1 - eta$lambda[j] * M2) /
        eta$theta[j]
    }
  }

  rex <- r_expectations(M, dt$W, dt$Z, dt$Z_sq)
  S_beta <- (dt$Y * rex$Er - t(apply(rex$Err, 1, function(Ei) {
    as.numeric(Ei %*% beta)
  }))) / vy
  colnames(S_beta) <- paste0("beta_", colnames(rex$Er))
  fitted_obs <- as.numeric(rex$Er %*% beta)
  S_vy <- -1 / (2 * vy) +
    (dt$Y^2 - 2 * dt$Y * fitted_obs + quad_form_rows(rex$Err, beta)) /
    (2 * vy^2)

  S <- cbind(S_alpha, tau2 = S_tau2, S_nu, S_lambda, S_theta,
             S_beta, vy = S_vy)
  storage.mode(S) <- "double"
  S
}

dot_psi <- function(dt, eta, beta, delta = DELTA) {
  Xmat <- cbind(1, dt$W, dt$Z, dt$Z_sq)
  colnames(Xmat) <- c("alpha0", "alpha_W", "alpha_Z", "alpha_Z2")
  m <- as.numeric(Xmat %*% eta$alpha)
  pnames <- c(colnames(Xmat), "tau2", paste0("nu", 1:3),
              paste0("lambda", 2:3), paste0("theta", 1:3),
              paste0("beta_", names(beta)), "vy")
  out <- setNames(rep(0, length(pnames)), pnames)

  dxi_dm <- 2 * beta["bL2"] * delta
  out[colnames(Xmat)] <- colMeans(Xmat) * dxi_dm
  out["beta_bL"] <- delta
  out["beta_bL2"] <- mean(2 * delta * m + delta^2)
  out["beta_bLW"] <- mean(delta * dt$W)
  out
}

xi_values <- function(dt, eta, beta, delta = DELTA) {
  Xmat <- cbind(1, dt$W, dt$Z, dt$Z_sq)
  m <- as.numeric(Xmat %*% eta$alpha)
  beta["bL"] * delta + beta["bL2"] * (2 * delta * m + delta^2) +
    beta["bLW"] * delta * dt$W
}

set.seed(SEED)
cat(sprintf(
  "\nObserved EIF smoke: scenario=%d n=%d seed=%d outcome=%s latent=%s error=%s\n",
  SCENARIO_ID, N, SEED, scen$outcome, scen$latent, scen$meas_error
))

dt <- generate_data(N, scen$meas_error, scen$outcome, scen$latent,
                    scen$indicators)
dt[, Z_sq := Z^2]

eta <- fit_measurement_eta(dt)
if (is.null(eta)) stop("Measurement/latent SEM failed.")
post_i <- indicator_posterior(dt, eta)
outcome <- fit_quadratic_outcome(dt, post_i)
post_y <- posterior_moments_y(dt, post_i, outcome$beta, outcome$vy)
S <- projected_scores(dt, eta, outcome$beta, outcome$vy, post_y)
dot <- dot_psi(dt, eta, outcome$beta)
dot <- dot[colnames(S)]

xi <- xi_values(dt, eta, outcome$beta)
psi_hat <- mean(xi)
D_X <- xi - psi_hat
S_centered <- scale(S, center = TRUE, scale = FALSE)
info <- crossprod(S_centered) / nrow(S_centered)
ridge <- max(1e-8, mean(diag(info), na.rm = TRUE) * 1e-6)
info_ridge <- info + diag(ridge, ncol(info))
coef_eta <- as.numeric(qr.solve(info_ridge, dot))
names(coef_eta) <- names(dot)
D_eta <- as.numeric(S_centered %*% coef_eta)
D_eff <- D_X + D_eta

cov_identity <- colMeans(D_eff * S_centered)
cov_gap <- cov_identity - dot
score_means <- colMeans(S)
eig <- eigen(info, symmetric = TRUE, only.values = TRUE)$values
finite_eig <- eig[is.finite(eig)]
cond <- if (length(finite_eig) && min(abs(finite_eig)) > 0) {
  max(abs(finite_eig)) / min(abs(finite_eig))
} else {
  Inf
}

z_cal <- (dt$L_true - post_y$mean) / sqrt(post_y$var)
summary_table <- data.frame(
  check = c("psi_hat_working", "true_psi_scenario", "posterior_z_mean",
            "posterior_z_sd", "max_abs_score_mean",
            "info_min_eigen", "info_condition", "ridge",
            "D_eff_mean", "D_eff_sd", "max_abs_cov_identity_gap"),
  value = c(
    psi_hat,
    scen$true_psi,
    mean(z_cal, na.rm = TRUE),
    sd(z_cal, na.rm = TRUE),
    max(abs(score_means), na.rm = TRUE),
    min(finite_eig, na.rm = TRUE),
    cond,
    ridge,
    mean(D_eff),
    sd(D_eff),
    max(abs(cov_gap), na.rm = TRUE)
  )
)

cat("\n=== Smoke Summary ===\n")
print(format(summary_table, digits = 4), row.names = FALSE)

cat("\n=== Nonzero dot_psi_eta Entries ===\n")
dot_nonzero <- dot[abs(dot) > 1e-10]
print(format(data.frame(parameter = names(dot_nonzero),
                        dot_psi = as.numeric(dot_nonzero)),
             digits = 4), row.names = FALSE)

cat("\n=== Largest Score Means ===\n")
ord_score <- order(abs(score_means), decreasing = TRUE)
print(format(data.frame(parameter = names(score_means)[head(ord_score, 8)],
                        mean_score = score_means[head(ord_score, 8)]),
             digits = 4), row.names = FALSE)

cat("\n=== Largest Covariance Identity Gaps ===\n")
ord_gap <- order(abs(cov_gap), decreasing = TRUE)
print(format(data.frame(parameter = names(cov_gap)[head(ord_gap, 8)],
                        cov_D_S = cov_identity[head(ord_gap, 8)],
                        dot_psi = dot[head(ord_gap, 8)],
                        gap = cov_gap[head(ord_gap, 8)]),
             digits = 4), row.names = FALSE)

saveRDS(
  list(
    config = list(n = N, scenario = scen, seed = SEED, delta = DELTA),
    eta = eta[setdiff(names(eta), "fit")],
    beta = outcome$beta,
    vy = outcome$vy,
    posterior = list(indicator_only = post_i, with_y = post_y),
    scores = S,
    score_means = score_means,
    observed_information = info,
    observed_information_ridge = info_ridge,
    dot_psi_eta = dot,
    D_X = D_X,
    D_eta = D_eta,
    D_eff = D_eff,
    covariance_identity = cov_identity,
    covariance_identity_gap = cov_gap,
    summary = summary_table
  ),
  OUT_FILE
)
cat(sprintf("\nSaved smoke object to %s\n", OUT_FILE))
