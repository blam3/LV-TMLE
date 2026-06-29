# ==============================================================================
# FILE: 05b_smoke_test_oaware.R
# PURPOSE: Outcome-aware smoke test. Same single cell as 05_smoke_test.R, but it
#          exercises the OUTCOME-AWARE estimator set from 06_outcome_aware.R via
#          run_comparison_v2(..., include_noY = TRUE):
#              SEM, RegCal (no-Y EAP), LV-TMLE (outcome-aware), Oracle, LV_noY (ablation)
#
#          This is the validation the handoff (§10, immediate step) calls for:
#          confirm the outcome-aware imputation fixes the attenuation that the
#          MI-without-Y "Robust" showed in 05_smoke_test.R.
#
# WHAT TO LOOK FOR  (true psi ~ 2.700 for this cell):
#   Oracle        ~ 2.70   bias ~ 0     engine ceiling                  [HARD]
#   LV (oaware)   ~ 2.70   bias ~ 0, NOT overshooting (+)               [HARD]
#   B > 0                  Rubin SE carries latent uncertainty          [HARD]
#   SEM           ~ 2.02   curvature-biased (~ -0.68), DGP wired right  [HARD]
#   RegCal(no-Y)  slightly < 2.70   mildly attenuated                   [DIAGNOSTIC]
#   LV_noY        clearly  < 2.70   attenuated ablation                 [DIAGNOSTIC]
#
# RUNTIME: noticeably heavier than 05 -- include_noY adds a SECOND block of M
#          shift-TMLE fits per rep (the ablation). For a fast sanity pass set
#          options(lv_tmle.reps = 15, lv_tmle.M = 10) before sourcing.
#
# ENVIRONMENT: R with lavaan, data.table, SuperLearner, earth, rpart, gam, dplyr,
#              and 00_dgp_variants.R / 01_setup.R / 06_outcome_aware.R in the wd.
# ==============================================================================

# ---- Field map from run_comparison_v2() --------------------------------------
F_EST <- list(SEM = "est_sem", RegCal = "est_regcal", LV = "est_lv",
              Oracle = "est_oracle", LV_noY = "est_lv_noY")
F_SE  <- list(SEM = "se_sem",  RegCal = "se_regcal",  LV = "se_lv",
              Oracle = "se_oracle",  LV_noY = "se_lv_noY")
F_B    <- "lv_B"       # diagnostics$  between-imputation variance (outcome-aware LV)
F_WBAR <- "lv_Wbar"    # diagnostics$  within-imputation variance
F_DIAG <- c("var_y_sem", "var_y_eff", "mean_slope_eff", "sd_slope_eff",
            "mean_shift_eff", "mean_mu_shift", "sd_mu_noY", "sd_mu_full",
            "draw_scale")
# ------------------------------------------------------------------------------

# ---- Config (mirror 05_smoke_test.R so the cells are comparable) -------------
read_int_knob <- function(option_name, env_name, default) {
  val <- getOption(option_name, NULL)
  if (is.null(val) || length(val) == 0 || !nzchar(as.character(val)[1])) {
    val <- Sys.getenv(env_name, unset = as.character(default))
  }
  out <- suppressWarnings(as.integer(val[1]))
  if (!is.finite(out) || out <= 0) default else out
}
read_num_knob <- function(option_name, env_name, default) {
  val <- getOption(option_name, NULL)
  if (is.null(val) || length(val) == 0 || !nzchar(as.character(val)[1])) {
    val <- Sys.getenv(env_name, unset = as.character(default))
  }
  out <- suppressWarnings(as.numeric(val[1]))
  if (!is.finite(out) || out <= 0) default else out
}

REPS <- read_int_knob("lv_tmle.reps", "LVTMLE_REPS", 40)
M    <- read_int_knob("lv_tmle.M",    "LVTMLE_M",    20)
V    <- read_int_knob("lv_tmle.V",    "LVTMLE_V",    5)
DRAW_SCALE <- read_num_knob("lv_tmle.draw_scale", "LVTMLE_DRAW_SCALE", 0.5)
SEED <- 12345

SCEN_OUTCOME    <- "nonlinear"
SCEN_LATENT     <- "normal"
SCEN_INDICATORS <- "continuous"
SCEN_N          <- 250
SCEN_ERROR      <- "medium"   # enough post_sd that B is clearly > 0
DELTA           <- 1

suppressPackageStartupMessages({ library(data.table); library(lavaan) })

# ---- Load DGPs + estimators (v1 engine) + outcome-aware layer (v2) -----------
source("00_dgp_variants.R")
source("01_setup.R")
source("06_outcome_aware.R")   # run_comparison_v2(), sem_posterior_full()

set.seed(SEED)
cat(sprintf("\nConfig: REPS=%d  M=%d  V=%d  draw_scale=%.2f  seed=%d\n",
            REPS, M, V, DRAW_SCALE, SEED))

# ---- True psi for this cell: cache AND an independent recompute ---------------
row <- subset(sim_grid,
              outcome == SCEN_OUTCOME & latent == SCEN_LATENT &
              indicators == SCEN_INDICATORS & n_size == SCEN_N &
              meas_error == SCEN_ERROR)
psi_cache <- row$true_psi[1]
psi_fresh <- compute_true_psi(OUTCOME_FUNCS[[SCEN_OUTCOME]],
                              LATENT_FUNCS[[SCEN_LATENT]],
                              delta = DELTA, n_mc = 2e6)$psi
PSI <- psi_cache
cat(sprintf("\nTrue psi  (cache) = %.4f   (fresh 2e6-draw recompute) = %.4f\n",
            psi_cache, psi_fresh))
if (abs(psi_cache - psi_fresh) > 0.02)
  cat("  !! cache and fresh psi disagree -- check the psi_lookup cache.\n")
set.seed(SEED)

# ---- Replication loop --------------------------------------------------------
cols <- c("SEM", "RegCal", "LV", "Oracle", "LV_noY")
est <- matrix(NA_real_, REPS, length(cols), dimnames = list(NULL, cols))
se_mat <- matrix(NA_real_, REPS, length(cols), dimnames = list(NULL, cols))
lvB <- lvW <- rep(NA_real_, REPS)
diag_mat <- matrix(NA_real_, REPS, length(F_DIAG), dimnames = list(NULL, F_DIAG))
rep_errors <- rep(NA_character_, REPS)
fail_location <- fail_reason <- rep(NA_character_, REPS)

cov_flag <- function(e, s) if (is.na(e) || is.na(s)) NA else
  ((e - 1.96 * s) <= PSI & PSI <= (e + 1.96 * s))
covm <- matrix(NA, REPS, length(cols), dimnames = list(NULL, cols))

# Tolerant field accessor (NA if the name is absent -- makes a wrong F_* loud
# rather than a silent crash).
get_field <- function(x, nm) { v <- x[[nm]]; if (is.null(v)) NA_real_ else as.numeric(v) }

t0 <- Sys.time()
for (b in 1:REPS) {
  dt <- generate_data(SCEN_N, SCEN_ERROR, SCEN_OUTCOME, SCEN_LATENT, SCEN_INDICATORS)

  r <- tryCatch(
    run_comparison_v2(dt, indicators_type = SCEN_INDICATORS,
                      M = M, V = V, include_noY = TRUE,
                      draw_scale = DRAW_SCALE),
    error = function(e) {
      rep_errors[b] <<- conditionMessage(e)
      NULL
    })

  if (!is.null(r)) {
    for (j in cols) {
      est[b, j] <- get_field(r, F_EST[[j]])
      se_mat[b, j]  <- get_field(r, F_SE[[j]])
      covm[b, j] <- cov_flag(est[b, j], se_mat[b, j])
    }
    lvB[b] <- get_field(r$diagnostics, F_B)
    lvW[b] <- get_field(r$diagnostics, F_WBAR)
    for (nm in F_DIAG) diag_mat[b, nm] <- get_field(r$diagnostics, nm)
    fail_location[b] <- r$diagnostics$fail_location
    fail_reason[b] <- r$diagnostics$fail_reason
  }
  if (b %% 5 == 0 || b == 1)
    cat(sprintf("  rep %d/%d  (elapsed %.1f min)\n", b, REPS,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(est = est, se = se_mat, lvB = lvB, lvW = lvW,
             posterior_diagnostics = diag_mat, errors = rep_errors,
             fail_location = fail_location, fail_reason = fail_reason, PSI = PSI),
        "smoke_test_oaware_raw.rds")

# ---- Summary table -----------------------------------------------------------
summ <- function(j) {
  e <- est[, j]; s <- se_mat[, j]; ok <- is.finite(e)
  n <- sum(ok)
  emp_sd <- sd(e[ok])
  if (n == 0) {
    return(data.frame(
      Estimator = j, n = 0,
      Mean = NA_real_, Bias = NA_real_, MCSE_bias = NA_real_,
      Emp_SD = NA_real_, Mean_SE = NA_real_, SE_ratio = NA_real_,
      Coverage = NA_real_
    ))
  }
  data.frame(
    Estimator = j, n = n,
    Mean = mean(e[ok]),
    Bias = mean(e[ok]) - PSI,
    MCSE_bias = emp_sd / sqrt(n),     # MC error on the mean -> how trustworthy 'Bias' is
    Emp_SD = emp_sd,
    Mean_SE = mean(s[ok], na.rm = TRUE),
    SE_ratio = mean(s[ok], na.rm = TRUE) / emp_sd,   # ~1 = well-calibrated SE
    Coverage = mean(covm[ok, j], na.rm = TRUE)
  )
}
tab <- do.call(rbind, lapply(cols, summ))

cat("\n============== OUTCOME-AWARE SMOKE TEST SUMMARY ==============\n")
cat(sprintf("Cell: outcome=%s latent=%s indic=%s  N=%d error=%s  | true psi=%.3f | REPS=%d \n\n",
            SCEN_OUTCOME, SCEN_LATENT, SCEN_INDICATORS, SCEN_N, SCEN_ERROR, PSI, REPS))
print(format(tab, digits = 3), row.names = FALSE)

mB <- mean(lvB, na.rm = TRUE); mW <- mean(lvW, na.rm = TRUE)
share <- mB / (mB + mW)
cat(sprintf("\nLV-TMLE (outcome-aware) variance decomposition:  within (W_bar)=%.4g  between (B)=%.4g  latent share=%.1f%%\n",
            mW, mB, 100 * share))

diag_avg <- colMeans(diag_mat, na.rm = TRUE)
cat(sprintf("\nOutcome-aware posterior diagnostics: SEM var_y=%.3f  outcome-link var_y=%.3f  mean dY/dL=%.3f  mean shift=%.3f  sd(mu_noY)=%.3f  sd(mu_full)=%.3f  draw_scale=%.2f\n",
            diag_avg["var_y_sem"], diag_avg["var_y_eff"],
            diag_avg["mean_slope_eff"], diag_avg["mean_shift_eff"],
            diag_avg["sd_mu_noY"],
            diag_avg["sd_mu_full"], diag_avg["draw_scale"]))

# ---- Pass / fail (infrastructure + the key method check) + diagnostics --------
get <- function(j, col) tab[tab$Estimator == j, col]
checks <- list()
add <- function(name, pass, msg) checks[[length(checks) + 1]] <<-
  list(name = name, pass = pass, msg = msg)

# Infra (HARD)
ok_pipe <- mean(is.finite(est[, "Oracle"]) & is.finite(est[, "LV"]) &
                is.finite(est[, "SEM"]) & is.finite(est[, "RegCal"]))
add("Pipeline completes", ok_pipe >= 0.8,
    sprintf("%.0f%% of reps produced SEM/RegCal/LV/Oracle", 100 * ok_pipe))

oracle_tol <- max(0.10, 2.5 * get("Oracle", "MCSE_bias"))
add("Oracle recovers psi (engine OK)", abs(get("Oracle", "Bias")) < oracle_tol,
    sprintf("Oracle bias = %+.3f (tol %.3f)", get("Oracle", "Bias"), oracle_tol))

add("B > 0 and propagates (Rubin SE OK)", is.finite(mB) && mB > 0 && share > 0.01,
    sprintf("B = %.4g, latent share = %.1f%%", mB, 100 * share))

add("Y update uses reduced residual variance", is.finite(diag_avg["var_y_eff"]) &&
      is.finite(diag_avg["var_y_sem"]) && diag_avg["var_y_eff"] < diag_avg["var_y_sem"],
    sprintf("var_y_eff = %.3f vs SEM var_y = %.3f",
            diag_avg["var_y_eff"], diag_avg["var_y_sem"]))

add("SEM is curvature-biased (DGP wired right)", get("SEM", "Bias") < -0.30,
    sprintf("SEM bias = %+.3f (expected ~ -0.68: coef ~2.0 vs psi ~2.7)", get("SEM", "Bias")))

# THE method check: outcome-aware LV recovers psi without overshooting (handoff §9).
lv_tol  <- max(0.10, 2.5 * get("LV", "MCSE_bias"))
lv_bias <- get("LV", "Bias")
if (!is.finite(lv_bias) || !is.finite(lv_tol)) {
  lv_status <- "LV estimate unavailable"
} else if (lv_bias > lv_tol) {
  lv_status <- "OVERSHOOT (+): outcome feedback too strong (§9b)"
} else if (lv_bias < -lv_tol) {
  lv_status <- "still attenuated (-): congeniality residual? (§9a)"
} else {
  lv_status <- "on target"
}
add("Outcome-aware LV recovers psi (no overshoot)", is.finite(lv_bias) &&
      is.finite(lv_tol) && abs(lv_bias) < lv_tol,
    sprintf("LV bias = %+.3f (tol %.3f); %s", lv_bias, lv_tol,
            lv_status))

# Method diagnostics (reported, not enforced)
regcal_b <- get("RegCal", "Bias")
lvnoY_b  <- get("LV_noY", "Bias")
regcal_mild <- is.finite(regcal_b) && regcal_b < 0 && regcal_b > -0.20
ablation_ok <- is.finite(lvnoY_b) && lvnoY_b < regcal_b && lvnoY_b < lv_bias

cat("\n--- Infrastructure / key checks ---\n")
for (c in checks)
  cat(sprintf("  [%s] %-46s %s\n", if (isTRUE(c$pass)) "PASS" else "FAIL", c$name, c$msg))

cat("\n--- Method diagnostics (not pass/fail) ---\n")
cat(sprintf("  RegCal(no-Y)  bias = %+.3f  -> %s\n", regcal_b,
            if (regcal_mild) "mildly attenuated, as expected (EAP plug-in; SE does not carry latent uncertainty)."
            else "NOT mildly attenuated -- inspect."))
cat(sprintf("  LV_noY (abl.) bias = %+.3f  -> %s\n", lvnoY_b,
            if (ablation_ok) "more attenuated than RegCal and outcome-aware LV -- ablation behaves as predicted."
            else "NOT the most attenuated -- ablation contrast weaker than expected; inspect."))

if (any(!is.na(rep_errors))) {
  cat("\n--- Replication errors ---\n")
  print(sort(table(rep_errors), decreasing = TRUE))
}
if (any(!is.na(fail_reason) & fail_reason != "None")) {
  cat("\n--- Structured fit failures ---\n")
  fail_tab <- as.data.frame(table(fail_location, fail_reason),
                            responseName = "n")
  fail_tab <- fail_tab[fail_tab$n > 0, , drop = FALSE]
  fail_tab <- fail_tab[order(fail_tab$n, decreasing = TRUE), , drop = FALSE]
  print(fail_tab, row.names = FALSE)
}

infra_pass <- all(vapply(checks, function(c) isTRUE(c$pass), logical(1)))
cat(sprintf("\n==> RESULT: %s\n",
    if (infra_pass) "PASS -- engine, Rubin SE, DGP, and the outcome-aware fix all behave as the handoff predicts. OK to wire run_comparison_v2 into 02."
    else "PROBLEM -- a hard check failed; do NOT wire run_comparison_v2 into 02 yet."))
