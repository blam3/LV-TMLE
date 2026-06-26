# ==============================================================================
# FILE: 00_dgp_variants.R
# PURPOSE: Modular data-generating processes for the LV-TMLE simulation, plus a
#          Monte Carlo routine that computes the TRUE target parameter (psi) for
#          each DGP once you leave the linear-Gaussian world.
#
# DESIGN PRINCIPLE (this is the thing that is easy to get wrong):
#   A DGP is composed of three swappable pieces --
#       (1) outcome_mean(L, W, Z)  : the DETERMINISTIC mean of Y  (E[Y|L,W,Z])
#       (2) gen_latent(n)          : draws (W, Z, L)   -- covariates + true latent
#       (3) gen_indicators(L,e_sd) : draws (I1, I2, I3) from L     -- measurement
#   The SAME outcome_mean() and gen_latent() functions are used BOTH to simulate
#   data AND to compute the true psi. There is therefore a single source of truth
#   per piece, so the estimand can never drift away from the data you generated.
#
#   KEY FACT: the true psi depends ONLY on (outcome_mean, gen_latent, delta).
#   It does NOT depend on the indicator model. Switching continuous <-> ordinal
#   indicators changes ESTIMATION difficulty, not the TRUTH. The non-normal
#   LATENT, by contrast, CAN change psi (through the outcome's curvature), so it
#   always gets its own MC evaluation.
#
#   ESTIMAND: psi = E[Y(L + delta)] - E[Y(L)]
#           = E_{(L,W,Z)}[ outcome_mean(L+delta, W, Z) - outcome_mean(L, W, Z) ]
#   (the Y-noise is mean-zero and cancels in expectation, so the MC uses the
#    noiseless mean function -- this is exact in the limit and low-variance.)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------------------------
# 0. SHARED CONSTANTS
# ------------------------------------------------------------------------------
DELTA  <- 1            # size of the latent shift defining the MTP
Y_SD   <- 1            # SD of the outcome noise
LAMBDA <- c(1.0, 0.9, 1.1)   # factor loadings for I1, I2, I3

# Theoretical Var(L) for the standard structural part used below:
#   L = 0.5*W + 0.5*Z + shock,  shock has mean 0, var 1,
#   Var(0.5W) = 0.25, Var(0.5Z) = 0.25 * Var(Unif(-2,2)) = 0.25 * (16/12) = 0.3333
#   => Var(L) = 0.25 + 0.3333 + 1 = 1.5833.
# All three latent generators below are variance-matched to this value so that
# the signal-to-noise ratio is comparable across latent shapes (a fair-comparison
# choice). If you change the structural coefficients or the shock variance, update
# VAR_L_THEORY or pass var_L explicitly to the ordinal generator.
VAR_L_THEORY <- 0.25 + 0.25 * (16 / 12) + 1   # = 1.583333...

# ------------------------------------------------------------------------------
# 1. OUTCOME MEAN FUNCTIONS  -- E[Y | L, W, Z]
# ------------------------------------------------------------------------------
# These return the deterministic conditional mean of Y. Noise is added only in
# the data generator, never here.

# (a) ORIGINAL: linear in L. The MTP shift effect equals the L coefficient (2.0)
#     regardless of the latent distribution, because E[(L+d) - L] = d.
#     This is the SEM's home turf -- include it as the honest negative control
#     where SEM SHOULD win.
outcome_mean_linear <- function(L, W, Z) {
  1 + 2.0 * L + 0.5 * W + sin(2 * Z)
}

# (b) NONLINEAR in L: quadratic curvature + a latent x covariate interaction.
#     - The L^2 term makes the true shift effect (psi) DIFFER from the linear
#       L->Y coefficient. A correctly fit *linear* SEM converges to the linear
#       coefficient (~2.0) and is therefore STRUCTURALLY BIASED for psi (~2.7).
#       This is the cleanest condition under which LV-TMLE can beat SEM on bias.
#     - The L*W term adds conditional heterogeneity that a linear SEM cannot
#       represent. (Its *marginal* contribution to psi is ~0 here because E[W]=0;
#       it breaks the conditional mean, not the average shift. The MC below
#       computes whatever it actually is -- never assume.)
outcome_mean_nonlinear <- function(L, W, Z) {
  1 + 2.0 * L + 0.7 * L^2 + 0.5 * L * W + 0.5 * W + sin(2 * Z)
}

OUTCOME_FUNCS <- list(
  linear    = outcome_mean_linear,
  nonlinear = outcome_mean_nonlinear
)

# ------------------------------------------------------------------------------
# 2. LATENT + COVARIATE GENERATORS  -- return data.frame(W, Z, L)
# ------------------------------------------------------------------------------
# All have E[L] = 0 and (by construction) matched Var(L) ~ 1.5833, so any
# difference in estimator performance across latent shapes is an ESTIMATION
# effect, not an estimand artifact.

# (a) ORIGINAL: normal latent.
gen_latent_normal <- function(n) {
  W <- rnorm(n, 0, 1)
  Z <- runif(n, -2, 2)
  L <- 0.5 * W + 0.5 * Z + rnorm(n, 0, 1)
  data.frame(W = W, Z = Z, L = L)
}

# (b) SKEWED latent: centered, scaled chi-square shock (mean 0, var 1, right-skew).
#     Breaks the normal-theory assumption underlying SEM's ML estimator (mainly
#     its standard errors / coverage; ML point estimates of linear structure are
#     partially robust to non-normality).
gen_latent_skew <- function(n) {
  W <- rnorm(n, 0, 1)
  Z <- runif(n, -2, 2)
  shock <- (rchisq(n, df = 3) - 3) / sqrt(6)   # mean 0, var 1, skewed
  L <- 0.5 * W + 0.5 * Z + shock
  data.frame(W = W, Z = Z, L = L)
}

# (c) MIXTURE (bimodal) latent: equal-weight two-component normal shock,
#     mean 0, var 1 (m^2 + s^2 = 0.75 + 0.25 = 1). Strongly non-normal.
gen_latent_mixture <- function(n) {
  W <- rnorm(n, 0, 1)
  Z <- runif(n, -2, 2)
  comp  <- rbinom(n, 1, 0.5)
  m <- sqrt(0.75); s <- 0.5
  shock <- ifelse(comp == 1, rnorm(n, m, s), rnorm(n, -m, s))
  L <- 0.5 * W + 0.5 * Z + shock
  data.frame(W = W, Z = Z, L = L)
}

LATENT_FUNCS <- list(
  normal  = gen_latent_normal,
  skew    = gen_latent_skew,
  mixture = gen_latent_mixture
)

# ------------------------------------------------------------------------------
# 3. INDICATOR (MEASUREMENT) MODELS  -- return data.frame(I1, I2, I3)
# ------------------------------------------------------------------------------
# These do NOT affect the true psi. They only change how hard L is to recover.

# (a) ORIGINAL: continuous Gaussian indicators. lavaan's linear-normal CFA is
#     correctly specified for the measurement model here.
gen_indicators_continuous <- function(L, e_sd, lambda = LAMBDA) {
  n <- length(L)
  data.frame(
    I1 = lambda[1] * L + rnorm(n, 0, e_sd),
    I2 = lambda[2] * L + rnorm(n, 0, e_sd),
    I3 = lambda[3] * L + rnorm(n, 0, e_sd)
  )
}

# (b) ORDINAL (Likert) indicators: the continuous indicator is thresholded into
#     n_cat ordered categories using FIXED (population) thresholds derived from
#     the theoretical marginal SD of each continuous indicator. Treating these as
#     continuous in a standard CFA (lavaan default, Pearson covariances) is
#     misspecified; a correct analysis needs ordered = TRUE / WLSMV.
#     Thresholds are fixed (not data-dependent) so they are identical across reps.
gen_indicators_ordinal <- function(L, e_sd, lambda = LAMBDA,
                                    n_cat = 5, var_L = VAR_L_THEORY) {
  n <- length(L)
  probs    <- seq(0, 1, length.out = n_cat + 1)
  interior <- probs[-c(1, length(probs))]          # interior cut probabilities
  make <- function(lam) {
    star <- lam * L + rnorm(n, 0, e_sd)
    sd_j <- sqrt(lam^2 * var_L + e_sd^2)            # theoretical marginal SD (E[I*]=0)
    thr  <- qnorm(interior, mean = 0, sd = sd_j)    # fixed, ~balanced thresholds
    as.integer(cut(star, breaks = c(-Inf, thr, Inf), labels = FALSE))
  }
  data.frame(I1 = make(lambda[1]), I2 = make(lambda[2]), I3 = make(lambda[3]))
}

INDIC_FUNCS <- list(
  continuous = gen_indicators_continuous,
  ordinal    = gen_indicators_ordinal
)

# ------------------------------------------------------------------------------
# 4. ASSEMBLER + DISPATCHER
# ------------------------------------------------------------------------------
# Build a generate_data() from the three pieces. Returns the same columns your
# pipeline expects (ID, W, Z, I1, I2, I3, Y) PLUS L_true.
#
#   *** L_true is included ONLY so you can build an oracle estimator (true L
#       plugged into the same TMLE) as a performance ceiling. Every non-oracle
#       estimator MUST ignore L_true. ***
make_generate_data <- function(outcome_mean, gen_latent, gen_indicators,
                               y_sd = Y_SD) {
  function(n, error_level) {
    e_sd <- switch(error_level, "small" = 0.2, "medium" = 0.6, "large" = 1.2,
                   stop("error_level must be 'small', 'medium', or 'large'"))
    d   <- gen_latent(n)
    ind <- gen_indicators(d$L, e_sd)
    Y   <- outcome_mean(d$L, d$W, d$Z) + rnorm(n, 0, y_sd)
    data.table(
      ID = seq_len(n), W = d$W, Z = d$Z,
      I1 = ind$I1, I2 = ind$I2, I3 = ind$I3,
      Y  = Y,
      L_true = d$L            # oracle only -- never use in real estimators
    )
  }
}

# One dispatcher that selects pieces by name. Note the new signature: you must
# pass outcome / latent / indicators in addition to n and error_level. Update the
# call in 02_run_sim.R accordingly (see "WIRING" note at the bottom).
generate_data <- function(n, error_level, outcome, latent, indicators) {
  make_generate_data(
    OUTCOME_FUNCS[[outcome]],
    LATENT_FUNCS[[latent]],
    INDIC_FUNCS[[indicators]]
  )(n, error_level)
}

# Convenience pre-assembled variants (each pairs a SEM-breaking feature with the
# nonlinear outcome, which is what gives LV-TMLE a chance to win on bias/RMSE):
generate_data_nonlinear <- make_generate_data(outcome_mean_nonlinear,
                                              gen_latent_normal,
                                              gen_indicators_continuous)
generate_data_ordinal   <- make_generate_data(outcome_mean_nonlinear,
                                              gen_latent_normal,
                                              gen_indicators_ordinal)
generate_data_nonnormal <- make_generate_data(outcome_mean_nonlinear,
                                              gen_latent_skew,
                                              gen_indicators_continuous)

# ------------------------------------------------------------------------------
# 5. THE MONTE CARLO TRUE-PSI ROUTINE  (this keeps the estimand correct)
# ------------------------------------------------------------------------------
# Computes psi = E[outcome_mean(L+delta, W, Z) - outcome_mean(L, W, Z)] by
# averaging the NOISELESS mean function over a huge draw from the true
# (latent, covariate) distribution. Chunked so memory stays bounded; returns a
# Monte Carlo standard error so you can confirm the precision is far below your
# simulation's resolution. Uses its own seed for reproducibility.
compute_true_psi <- function(outcome_mean, gen_latent, delta = DELTA,
                             n_mc = 2e7, chunk = 1e6, seed = 20240101) {
  set.seed(seed)
  n_done <- 0; s1 <- 0; s2 <- 0
  while (n_done < n_mc) {
    m <- min(chunk, n_mc - n_done)
    d <- gen_latent(m)
    diff <- outcome_mean(d$L + delta, d$W, d$Z) - outcome_mean(d$L, d$W, d$Z)
    s1 <- s1 + sum(diff)
    s2 <- s2 + sum(diff^2)
    n_done <- n_done + m
  }
  psi  <- s1 / n_done
  varc <- (s2 - s1^2 / n_done) / (n_done - 1)
  list(psi = psi, mc_se = sqrt(varc / n_done), n_mc = n_done)
}

# Build a lookup table of true psi for every (outcome x latent) combination.
# (Indicators are intentionally excluded -- they do not affect psi.)
build_psi_lookup <- function(outcomes = names(OUTCOME_FUNCS),
                             latents  = names(LATENT_FUNCS),
                             delta = DELTA, n_mc = 2e7) {
  rows <- list()
  for (oc in outcomes) for (lt in latents) {
    tp <- compute_true_psi(OUTCOME_FUNCS[[oc]], LATENT_FUNCS[[lt]],
                           delta = delta, n_mc = n_mc)
    rows[[length(rows) + 1]] <- data.frame(
      outcome = oc, latent = lt,
      true_psi = tp$psi, mc_se = tp$mc_se
    )
  }
  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------
# 6. DIAGNOSTIC: what a correctly-fit LINEAR SEM actually targets
# ------------------------------------------------------------------------------
# Under the nonlinear DGP, the SEM fits Y ~ L + W + Z + Z^2 (no L^2) and so
# converges to the population linear-projection coefficient on L. This function
# returns that coefficient (computed on a large ORACLE sample with true L), so
# you can SEE the structural gap between the SEM's estimand and the true psi.
# For a symmetric latent this comes out ~2.0 while true psi ~2.7 -- i.e. the SEM
# is biased for psi by ~the curvature contribution, by construction.
sem_linear_target <- function(outcome_mean, gen_latent, n_mc = 2e6,
                              seed = 7) {
  set.seed(seed)
  d <- gen_latent(n_mc)
  Y <- outcome_mean(d$L, d$W, d$Z)            # noiseless mean; coefficient is the same
  fit <- lm(Y ~ d$L + d$W + d$Z + I(d$Z^2))
  unname(coef(fit)["d$L"])
}

# ------------------------------------------------------------------------------
# 7. EXTEND THE SIMULATION GRID (replaces the old expand.grid in 01_setup.R)
# ------------------------------------------------------------------------------
# A full factorial. The honest "negative control" cell is
# outcome=linear, latent=normal, indicators=continuous: SEM should WIN there,
# and showing that is evidence your study is trustworthy.
build_sim_grid <- function() {
  g <- expand.grid(
    n_size     = c(50, 250, 1000),
    meas_error = c("small", "medium", "large"),
    outcome    = c("linear", "nonlinear"),
    latent     = c("normal", "skew"),          # add "mixture" if desired
    indicators = c("continuous", "ordinal"),
    stringsAsFactors = FALSE
  )
  g$Scenario_ID <- seq_len(nrow(g))
  g
}

# Merge the true psi onto the grid. Each scenario then carries its OWN true_psi,
# which downstream code MUST use in place of the old global TRUE_PSI = 2.0.
attach_true_psi <- function(sim_grid, psi_lookup) {
  merge(sim_grid, psi_lookup[, c("outcome", "latent", "true_psi")],
        by = c("outcome", "latent"), all.x = TRUE, sort = FALSE)
}

# ==============================================================================
# 8. DEMO / SANITY CHECK  (run this block to verify the estimand is correct)
# ==============================================================================
# Set RUN_DGP_DEMO <- TRUE before sourcing, or run interactively. Uses a smaller
# n_mc for speed; bump to 2e7 for publication-grade precision.
if (isTRUE(getOption("run_dgp_demo", FALSE)) ||
    (exists("RUN_DGP_DEMO") && isTRUE(RUN_DGP_DEMO))) {

  cat("\n--- TRUE psi by (outcome x latent), delta =", DELTA, "---\n")
  psi_tab <- build_psi_lookup(n_mc = 5e6)          # demo precision
  print(psi_tab, row.names = FALSE)

  cat("\n--- Sanity checks ---\n")
  lin_normal <- psi_tab$true_psi[psi_tab$outcome == "linear"   & psi_tab$latent == "normal"]
  non_normal <- psi_tab$true_psi[psi_tab$outcome == "nonlinear" & psi_tab$latent == "normal"]
  cat(sprintf("linear/normal     true psi ~ %.3f   (expected ~2.000: linear shift = coefficient)\n", lin_normal))
  cat(sprintf("nonlinear/normal  true psi ~ %.3f   (expected ~2.700: 2*delta + 0.7*delta^2, E[L]=0)\n", non_normal))

  cat("\n--- The SEM gap under the nonlinear DGP ---\n")
  sem_t <- sem_linear_target(outcome_mean_nonlinear, gen_latent_normal)
  cat(sprintf("Linear-SEM estimand (coef on L)  ~ %.3f\n", sem_t))
  cat(sprintf("True MTP shift psi               ~ %.3f\n", non_normal))
  cat(sprintf("=> SEM is structurally biased for psi by ~ %.3f (the curvature it cannot see)\n",
              non_normal - sem_t))

  cat("\n--- Quick structural-invariance check (psi vs latent shape, nonlinear outcome) ---\n")
  for (lt in names(LATENT_FUNCS)) {
    tp <- compute_true_psi(outcome_mean_nonlinear, LATENT_FUNCS[[lt]], n_mc = 5e6)
    cat(sprintf("  latent = %-8s  true psi = %.3f (+/- %.4f)\n", lt, tp$psi, tp$mc_se))
  }
  cat("(With E[L]=0 and matched Var(L), these are ~equal: latent shape changes\n",
      " estimation difficulty, not the estimand. Differences would appear only if\n",
      " you also change E[L] or Var(L).)\n")
}
