# ==============================================================================
# FILE: 05_smoke_test.R
# PURPOSE: Fast single-scenario diagnostic BEFORE committing the 72-scenario array.
#
# It runs the REAL pipeline (run_comparison: Naive, SEM, Robust/LV-TMLE, Oracle)
# on the (outcome=nonlinear, latent=normal, indicators=continuous, n=250) cell,
# and adds an inline RegCal estimator (regression/EAP factor scores -> same TMLE
# engine) so you can see the full de-attenuation spectrum in one run.
#
# WHAT TO LOOK FOR
#   true psi for this cell ~ 2.700.
#   Oracle  ~ 2.70  (engine works; this is the ceiling)            [HARD check]
#   SEM     ~ 2.00  (linear coef misses the curvature; biased low) [HARD check]
#   B > 0   (latent uncertainty enters the Rubin SE)               [HARD check]
#   Naive(Bartlett)  : attenuated (< 2.7)
#   Robust(MI no-Y)  : likely ALSO attenuated -- the issue to catch [DIAGNOSTIC]
#   RegCal(EAP score): closer to 2.7 than Naive                     [DIAGNOSTIC]
#
# RUNTIME: ~10-20 min single core at defaults. For a 2-5 min sanity pass set
#          REPS <- 15 and M <- 10 below.
# ==============================================================================

# ---- Config (the speed knobs) -----------------------------------------------
REPS <- 40          # replications
M    <- 20          # posterior draws for LV-TMLE
V    <- 5           # cross-fitting folds
SEED <- 12345

SCEN_OUTCOME    <- "nonlinear"
SCEN_LATENT     <- "normal"
SCEN_INDICATORS <- "continuous"
SCEN_N          <- 250
SCEN_ERROR      <- "medium"   # enough post_sd that B is clearly > 0
DELTA           <- 1

suppressPackageStartupMessages({ library(data.table); library(lavaan) })

# ---- Load DGPs + estimators -------------------------------------------------
# NOTE: sourcing 01_setup.R will build psi_lookup.rds on first run (~1 min,
# 6 (outcome x latent) combos at 2e7 draws). That cache is also what production
# needs, so this is a legitimate one-time cost, not throwaway work.
source("00_dgp_variants.R")
source("01_setup.R")

set.seed(SEED)

# ---- True psi for this cell: from the cache AND an independent recompute -----
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

# ---- Inline RegCal diagnostic: EAP/regression factor scores -> shift_tmle ----
# Demonstrates the de-attenuating plug-in (regression calibration). One extra
# SEM fit per rep; cheap relative to the SuperLearner cost.
cfa_model_str <- '
  L =~ I1 + I2 + I3
  L ~ W + Z + Z_sq
  Y ~ L + W + Z + Z_sq
'
regcal_estimate <- function(dt, sl_lib, V = 5, delta = 1) {
  d <- copy(dt); d[, Z_sq := Z^2]
  fit <- try(suppressWarnings(
    sem(cfa_model_str, data = as.data.frame(d),
        std.lv = FALSE, auto.fix.first = TRUE, ncpus = 1L)), silent = TRUE)
  if (inherits(fit, "try-error") || !lavInspect(fit, "converged"))
    return(list(psi = NA_real_, se = NA_real_))
  mu_reg <- try(lavPredict(fit, method = "regression")[, "L"], silent = TRUE)
  if (inherits(mu_reg, "try-error") || any(!is.finite(mu_reg)))
    return(list(psi = NA_real_, se = NA_real_))
  shift_tmle(d$Y, as.numeric(mu_reg), as.data.frame(d[, .(W, Z)]),
             sl_lib, delta = delta, V = V)
}

# Match run_comparison's learner library for this n.
sl_lib <- lv_sl_library(SCEN_N)

# ---- Replication loop -------------------------------------------------------
cols <- c("Naive", "SEM", "Robust", "Oracle", "RegCal")
est <- se <- matrix(NA_real_, REPS, length(cols), dimnames = list(NULL, cols))
lvB <- lvW <- rep(NA_real_, REPS)

cov_flag <- function(e, s) if (is.na(e) || is.na(s)) NA else
  ((e - 1.96 * s) <= PSI & PSI <= (e + 1.96 * s))
covm <- matrix(NA, REPS, length(cols), dimnames = list(NULL, cols))

t0 <- Sys.time()
for (b in 1:REPS) {
  dt <- generate_data(SCEN_N, SCEN_ERROR, SCEN_OUTCOME, SCEN_LATENT, SCEN_INDICATORS)

  r <- tryCatch(run_comparison(dt, indicators_type = SCEN_INDICATORS, M = M, V = V),
                error = function(e) NULL)
  rc <- tryCatch(regcal_estimate(dt, sl_lib, V = V, delta = DELTA),
                 error = function(e) list(psi = NA, se = NA))

  if (!is.null(r)) {
    est[b, ] <- c(r$est_naive, r$est_sem, r$est_rob, r$est_oracle, rc$psi)
    se[b, ]  <- c(r$se_naive,  r$se_sem,  r$se_rob,  r$se_oracle,  rc$se)
    lvB[b] <- r$diagnostics$lv_B; lvW[b] <- r$diagnostics$lv_Wbar
    for (j in cols) covm[b, j] <- cov_flag(est[b, j], se[b, j])
  }
  if (b %% 5 == 0 || b == 1)
    cat(sprintf("  rep %d/%d  (elapsed %.1f min)\n", b, REPS,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(est = est, se = se, lvB = lvB, lvW = lvW, PSI = PSI),
        "smoke_test_raw.rds")

# ---- Summary table ----------------------------------------------------------
summ <- function(j) {
  e <- est[, j]; s <- se[, j]; ok <- is.finite(e)
  n <- sum(ok)
  bias <- mean(e[ok]) - PSI
  emp_sd <- sd(e[ok])
  data.frame(
    Estimator = j, n = n,
    Mean = mean(e[ok]),
    Bias = bias,
    MCSE_bias = emp_sd / sqrt(n),     # MC error on the mean -> how trustworthy 'Bias' is
    Emp_SD = emp_sd,
    Mean_SE = mean(s[ok], na.rm = TRUE),
    SE_ratio = mean(s[ok], na.rm = TRUE) / emp_sd,   # ~1 = well-calibrated SE
    Coverage = mean(covm[ok, j], na.rm = TRUE)
  )
}
tab <- do.call(rbind, lapply(cols, summ))

cat("\n================ SMOKE TEST SUMMARY ================\n")
cat(sprintf("Cell: outcome=%s latent=%s indic=%s  N=%d error=%s  | true psi=%.3f | REPS=%d \n\n",
            SCEN_OUTCOME, SCEN_LATENT, SCEN_INDICATORS, SCEN_N, SCEN_ERROR, PSI, REPS))
print(format(tab, digits = 3), row.names = FALSE)

mB <- mean(lvB, na.rm = TRUE); mW <- mean(lvW, na.rm = TRUE)
share <- mB / (mB + mW)
cat(sprintf("\nLV-TMLE variance decomposition:  within (W_bar)=%.4g  between (B)=%.4g  latent share=%.1f%%\n",
            mW, mB, 100 * share))

# ---- Pass / fail (infrastructure) + diagnostics -----------------------------
get <- function(j, col) tab[tab$Estimator == j, col]
checks <- list()
add <- function(name, pass, msg) checks[[length(checks) + 1]] <<-
  list(name = name, pass = pass, msg = msg)

# Infrastructure (HARD)
ok_pipe <- mean(is.finite(est[, "Oracle"]) & is.finite(est[, "Robust"]) &
                is.finite(est[, "SEM"]) & is.finite(est[, "Naive"]))
add("Pipeline completes", ok_pipe >= 0.8,
    sprintf("%.0f%% of reps produced all four estimates", 100 * ok_pipe))

oracle_tol <- max(0.10, 2.5 * get("Oracle", "MCSE_bias"))
add("Oracle recovers psi (engine OK)", abs(get("Oracle", "Bias")) < oracle_tol,
    sprintf("Oracle bias = %+.3f (tol %.3f)", get("Oracle", "Bias"), oracle_tol))

add("B > 0 and propagates (Rubin SE OK)", is.finite(mB) && mB > 0 && share > 0.01,
    sprintf("B = %.4g, latent share = %.1f%%", mB, 100 * share))

add("SEM is curvature-biased (DGP wired right)", get("SEM", "Bias") < -0.30,
    sprintf("SEM bias = %+.3f (expected ~ -0.70: coef ~2.0 vs psi ~2.7)", get("SEM", "Bias")))

# Method-level (DIAGNOSTIC, reported not enforced)
rob_att   <- get("Robust", "Bias") < -0.15
regcal_ok <- abs(get("RegCal", "Bias")) < abs(get("Naive", "Bias"))

cat("\n--- Infrastructure checks ---\n")
for (c in checks)
  cat(sprintf("  [%s] %-40s %s\n", if (c$pass) "PASS" else "FAIL", c$name, c$msg))

cat("\n--- Method diagnostics (not pass/fail) ---\n")
cat(sprintf("  Robust(MI no-Y) bias = %+.3f  -> %s\n", get("Robust", "Bias"),
            if (rob_att) "ATTENUATED, as predicted: impute conditions on indicators only, not Y."
            else "not strongly attenuated (better than expected -- inspect)."))
cat(sprintf("  RegCal(EAP)     bias = %+.3f  -> %s than Naive(Bartlett) bias %+.3f\n",
            get("RegCal", "Bias"),
            if (regcal_ok) "CLOSER to psi" else "NOT closer", get("Naive", "Bias")))

infra_pass <- all(vapply(checks, function(c) c$pass, logical(1)))
cat(sprintf("\n==> INFRASTRUCTURE: %s\n", if (infra_pass) "OK -- engine, Rubin SE, and DGP are sound."
                                          else "PROBLEM -- fix before scaling up."))
cat("==> METHOD: the point-estimate fix (regression scores + outcome-aware draws)\n",
    "    still needs to be settled before the full array. See notes.\n")
