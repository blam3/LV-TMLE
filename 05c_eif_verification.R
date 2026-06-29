# ==============================================================================
# FILE: 05c_eif_verification.R
# PURPOSE: Pressure-test the STRUCTURE of the observed-data EIF for the latent
#          MTP-shift parameter BEFORE the algebra goes in the manuscript.
#
#   The full-data (oracle) EIF is settled:
#       D_F = H(L,X){Y - Q(L,X)} + Q(L+d,X) - psi1 - (Y - E[Y]).
#   The observed data is O = (W,Z,I,Y); only L is integrated out via the
#   posterior pi(l|O) ∝ g(l|X) m(I|l) p(Y|l,X)  -- note the Y factor: the EIF
#   FORCES outcome-aware conditioning, it is not a heuristic.
#
#   We cannot yet write a closed-form observed EIF D* (the latent-reconstruction
#   inflation term has no closed form under a nonlinear Q). So this harness does
#   NOT check D* directly; it checks the PROPERTIES any correct D* must satisfy.
#   It can FALSIFY the structure; it cannot prove it.
#
#   FOUR BLOCKS  (HARD = non-circular falsifier; DIAG = quantification/consistency)
#     1. Posterior calibration   [HARD-ish] : is pi(l|O) (incl. the Y term) right?
#                                 z = (L_true - Lbar)/sd_post must be ~ N(0,1).
#     2. Coarsening inequality   [HARD]     : Var(E_pi[D_F|O]) <= Var(D_F), and the
#                                 law-of-total-variance identity must close.
#     3. Gradient / eta identity [DIAG]     : E_n[Dbar * S_eta] vs dpsi/deta.
#                                 Size of the gap = size of the eta-estimation
#                                 correction => evidence for/against "estimated eta".
#     4. SE calibration          [HARD]     : with-Y Rubin SE vs Monte-Carlo SD of
#                                 psi-hat across reps; and de-attenuation vs no-Y.
#
#   SIMPLIFICATIONS (deliberate; see chat):
#     - Blocks 1-3 fit Q on the TRUE L (oracle Q) to isolate EIF algebra from
#       Q-estimation noise.
#     - The importance-sampled congeniality probe (conjugate-Q vs flexible-Q
#       posterior) is DEFERRED to the estimator build; Block 4 uses the existing
#       no-Y MI as the ablation instead.
#
#   RUNTIME: a few minutes single core at defaults.
# ==============================================================================

# ---- Config -----------------------------------------------------------------
N_BIG <- 2000;  M_BIG <- 200            # Blocks 1-3: one large sample, many draws
REPS  <- 40;    M_REP <- 25;  N_REP <- 250;  V <- 5   # Block 4: variance calibration
ERROR   <- "medium"
OUTCOME <- "nonlinear"; LATENT <- "normal"; INDIC <- "continuous"
DELTA <- 1; SEED <- 4242

suppressPackageStartupMessages({ library(data.table); library(lavaan) })
source("00_dgp_variants.R")
source("01_setup.R")          # safe_SL, robust_predict, shift_tmle,
                              # fit_sem_posterior, lv_tmle_mi, sim_grid
set.seed(SEED)

PSI <- subset(sim_grid, outcome == OUTCOME & latent == LATENT &
              indicators == INDIC & n_size == N_REP & meas_error == ERROR)$true_psi[1]

# dpsi/d(mu_L) for this DGP. psi = 2d + 0.7 d^2 + 1.4 d * E[L]  (see appendix in
# PROJECT_HANDOFF), so perturbing the latent mean moves psi at rate 1.4*delta.
DPSI_DMU <- 1.4 * DELTA

sl_lib <- lv_sl_library(N_REP)

# ==============================================================================
# HELPERS
# ==============================================================================

# Full SEM fit with meanstructure so we recover indicator/outcome intercepts that
# the with-Y posterior needs. Returns the parameter bundle, or NULL on failure.
fit_sem_full <- function(dt) {
  d <- copy(dt); d[, Z_sq := Z^2]
  model <- 'L =~ I1 + I2 + I3
            L ~ W + Z + Z_sq
            Y ~ L + W + Z + Z_sq'
  fit <- try(suppressWarnings(
    sem(model, data = as.data.frame(d), meanstructure = TRUE,
        std.lv = FALSE, auto.fix.first = TRUE, ncpus = 1L)), silent = TRUE)
  if (inherits(fit, "try-error") || !lavInspect(fit, "converged")) return(NULL)
  pe <- parameterEstimates(fit)
  g <- function(lhs, op, rhs) {
    r <- pe$est[pe$lhs == lhs & pe$op == op & pe$rhs == rhs]
    if (length(r)) r[1] else NA_real_
  }
  lam <- c(g("L","=~","I1"), g("L","=~","I2"), g("L","=~","I3"))
  th  <- pmax(c(g("I1","~~","I1"), g("I2","~~","I2"), g("I3","~~","I3")), 1e-4)
  nu  <- c(g("I1","~1",""), g("I2","~1",""), g("I3","~1","")); nu[is.na(nu)] <- 0
  list(
    lam = lam, th = th, nu = nu,
    psiLL  = max(g("L","~~","L"), 1e-4),
    aL     = ifelse(is.na(g("L","~1","")), 0, g("L","~1","")),
    bL     = c(W = g("L","~","W"), Z = g("L","~","Z"), Zsq = g("L","~","Z_sq")),
    beta_LY= g("Y","~","L"),
    bY     = c(W = g("Y","~","W"), Z = g("Y","~","Z"), Zsq = g("Y","~","Z_sq")),
    iY     = ifelse(is.na(g("Y","~1","")), 0, g("Y","~1","")),
    sY2    = max(g("Y","~~","Y"), 1e-4)
  )
}

# Posterior moments of L given O. withY = TRUE adds the outcome factor (the
# EIF-prescribed posterior); withY = FALSE is the indicator-only posterior (what
# 01_setup currently draws from). Homoskedastic complete-data model -> sd is a
# scalar (constant across subjects), as expected.
post_moments <- function(d2, P, withY = TRUE) {
  mu_g <- P$aL + P$bL["W"]*d2$W + P$bL["Z"]*d2$Z + P$bL["Zsq"]*d2$Z_sq
  Icen <- sweep(as.matrix(d2[, .(I1, I2, I3)]), 2, P$nu, "-")
  prec <- 1/P$psiLL + sum(P$lam^2 / P$th)
  num  <- as.numeric(mu_g)/P$psiLL + as.numeric(Icen %*% (P$lam / P$th))
  if (withY) {
    resY <- d2$Y - P$iY - P$bY["W"]*d2$W - P$bY["Z"]*d2$Z - P$bY["Zsq"]*d2$Z_sq
    prec <- prec + P$beta_LY^2 / P$sY2
    num  <- num  + P$beta_LY * as.numeric(resY) / P$sY2
  }
  list(Lbar = num / prec, psd = sqrt(1 / prec), mu_g = as.numeric(mu_g))
}

# Oracle outcome/density fit (on TRUE L) -> callable Q(.) and a Gaussian g(.|X).
fit_oracle_QG <- function(dt, sl_lib) {
  X    <- as.data.frame(dt[, .(W, Z)])
  base <- data.frame(Y = dt$Y, L_sim = dt$L_true, X)
  qfit <- safe_SL(base$Y,     base[, c("L_sim","W","Z")], gaussian(), sl_lib)
  gfit <- safe_SL(base$L_sim, base[, c("W","Z")],         gaussian(), sl_lib)
  muG  <- robust_predict(gfit, base, c("W","Z"))
  sigG <- sd(base$L_sim - muG); if (!is.finite(sigG) || sigG < 1e-6) sigG <- 1
  list(
    Qf   = function(ell, Xdf) robust_predict(qfit, data.frame(L_sim = ell, Xdf),
                                             c("L_sim","W","Z")),
    muG = muG, sigG = sigG
  )
}

# Full-data EIF evaluated at an arbitrary latent vector `ell`, using the oracle QG.
DF_at <- function(ell, dt, QG, psi1, Ybar, delta = DELTA,
                  gbound = 1e-3, hcap = 30) {
  X    <- as.data.frame(dt[, .(W, Z)])
  Q_l  <- QG$Qf(ell,         X)
  Q_lp <- QG$Qf(ell + delta, X)
  g_l  <- pmax(dnorm(ell,         QG$muG, QG$sigG), gbound)
  g_lm <- pmax(dnorm(ell - delta, QG$muG, QG$sigG), gbound)
  H_l  <- pmin(g_lm / g_l, hcap)
  H_l * (dt$Y - Q_l) + Q_lp - psi1 - (dt$Y - Ybar)
}

# With-Y multiple-imputation TMLE + Rubin SE (the EIF-prescribed estimator).
lv_mi_withY <- function(dt, P, sl_lib, M, V, delta = DELTA) {
  d2 <- copy(dt); d2[, Z_sq := Z^2]
  pm <- post_moments(d2, P, withY = TRUE)
  X  <- as.data.frame(dt[, .(W, Z)]); n <- nrow(dt)
  psis <- vars <- rep(NA_real_, M)
  for (m in 1:M) {
    A <- rnorm(n, pm$Lbar, pm$psd)
    f <- tryCatch(shift_tmle(dt$Y, A, X, sl_lib, delta, V), error = function(e) NULL)
    if (!is.null(f) && is.finite(f$psi) && is.finite(f$se)) {
      psis[m] <- f$psi; vars[m] <- f$se^2
    }
  }
  ok <- is.finite(psis) & is.finite(vars)
  if (sum(ok) < 2) return(list(psi = NA, se = NA))
  psis <- psis[ok]; vars <- vars[ok]; Mu <- length(psis)
  B <- var(psis); Tv <- mean(vars) + (1 + 1/Mu) * B
  list(psi = mean(psis), se = sqrt(Tv), B = B)
}

# ==============================================================================
# BLOCKS 1-3 : one large dataset, oracle Q
# ==============================================================================
cat(sprintf("\nCell: outcome=%s latent=%s indic=%s | N_big=%d M_big=%d | true psi=%.3f\n",
            OUTCOME, LATENT, INDIC, N_BIG, M_BIG, PSI))

dt_big <- generate_data(N_BIG, ERROR, OUTCOME, LATENT, INDIC)
d2_big <- copy(dt_big); d2_big[, Z_sq := Z^2]
P  <- fit_sem_full(dt_big)
if (is.null(P)) stop("SEM failed on the big dataset; rerun with a new seed.")

QG    <- fit_oracle_QG(dt_big, sl_lib)
Xbig  <- as.data.frame(dt_big[, .(W, Z)])
Ybar  <- mean(dt_big$Y)
psi1  <- mean(QG$Qf(dt_big$L_true + DELTA, Xbig))      # oracle E[Q(L+d,X)]

pmY  <- post_moments(d2_big, P, withY = TRUE)
pm0  <- post_moments(d2_big, P, withY = FALSE)

# ---- BLOCK 1: posterior calibration -----------------------------------------
z  <- (dt_big$L_true - pmY$Lbar) / pmY$psd
b1_mean <- mean(z); b1_sd <- sd(z)
b1_pass <- abs(b1_mean) < 0.10 && abs(b1_sd - 1) < 0.15

cat("\n=== BLOCK 1  Posterior calibration (validates pi(l|O), incl. Y term) ===\n")
cat(sprintf("  post sd: no-Y = %.4f   with-Y = %.4f   (Y must shrink it: info gain)\n",
            pm0$psd, pmY$psd))
cat(sprintf("  z = (L_true - Lbar)/sd_post :  mean = %+.3f (~0)   sd = %.3f (~1)\n",
            b1_mean, b1_sd))
cat(sprintf("  [%s] posterior is calibrated\n", if (b1_pass) "PASS" else "FAIL"))

# ---- BLOCK 2: coarsening inequality + law of total variance -----------------
DF_true <- DF_at(dt_big$L_true, dt_big, QG, psi1, Ybar)

DFm <- matrix(NA_real_, N_BIG, M_BIG)
for (m in 1:M_BIG) {
  ell <- rnorm(N_BIG, pmY$Lbar, pmY$psd)
  DFm[, m] <- DF_at(ell, dt_big, QG, psi1, Ybar)
}
Dbar   <- rowMeans(DFm)               # E_pi[D_F | O]
within <- apply(DFm, 1, var)          # Var_pi(D_F | O)

V_full   <- var(DF_true)
V_bar    <- var(Dbar)
E_within <- mean(within)
ltv_resid <- V_full - (V_bar + E_within)          # law of total variance -> ~0

b2_ineq <- V_bar < V_full                          # coarsening contraction
b2_ltv  <- abs(ltv_resid) / V_full < 0.15
b2_pass <- b2_ineq && b2_ltv

cat("\n=== BLOCK 2  Coarsening inequality + LTV identity [HARD falsifier] ===\n")
cat(sprintf("  Var(D_F oracle)      = %.4f\n", V_full))
cat(sprintf("  Var(E_pi[D_F|O])     = %.4f   (must be < oracle)\n", V_bar))
cat(sprintf("  E[Var_pi(D_F|O)]     = %.4f   (latent info loss = %.1f%% of oracle var)\n",
            E_within, 100 * E_within / V_full))
cat(sprintf("  LTV residual         = %+.4f  (%.1f%% of oracle var; should be ~0)\n",
            ltv_resid, 100 * ltv_resid / V_full))
cat(sprintf("  [%s] contraction holds   [%s] LTV closes\n",
            if (b2_ineq) "PASS" else "FAIL", if (b2_ltv) "PASS" else "FAIL"))

# ---- BLOCK 3: gradient / eta identity ---------------------------------------
# Observed-data eta-scores (posterior means of the full-data CFA scores).
L1 <- pmY$Lbar; L2 <- pmY$psd^2 + pmY$Lbar^2     # E_pi[L], E_pi[L^2]
S_lam2  <- ((d2_big$I2 - P$nu[2]) * L1 - P$lam[2] * L2) / P$th[2]   # pure nuisance
S_sigL2 <- ((pmY$psd^2 + (L1 - pmY$mu_g)^2) / (2*P$psiLL^2)) - 1/(2*P$psiLL)  # nuisance
S_alpha <- (L1 - pmY$mu_g) / P$psiLL                                 # moves psi

scores  <- list(lambda2 = S_lam2, sigmaL2 = S_sigL2, alpha = S_alpha)
targets <- c(lambda2 = 0, sigmaL2 = 0, alpha = DPSI_DMU)   # dpsi/deta

cat("\n=== BLOCK 3  Gradient / eta identity [DIAG: sizes the eta correction] ===\n")
cat("  E_n[Dbar * S_eta] vs dpsi/deta. Gap = magnitude of the eta-estimation\n")
cat("  correction along that direction (nonzero => 'estimated eta' is required).\n")
b3 <- data.frame()
for (k in names(scores)) {
  Chat <- mean(Dbar * scores[[k]])
  gap  <- Chat - targets[k]
  b3 <- rbind(b3, data.frame(direction = k, E_Dbar_S = Chat,
                             dpsi_deta = targets[k], gap = gap))
}
print(format(b3, digits = 3), row.names = FALSE)
alpha_ok <- abs(b3$gap[b3$direction == "alpha"]) <
            abs(b3$gap[b3$direction == "lambda2"]) + abs(b3$gap[b3$direction == "sigmaL2"]) + 1e-8
cat(sprintf("  Reading: gaps on lambda2/sigmaL2 are the nuisance corrections (expect > 0,\n",
            ""))
cat("           else fixing eta would be defensible). The alpha row tests whether\n")
cat("           Dbar already encodes the target derivative ~1.4*delta.\n")

# ==============================================================================
# BLOCK 4 : SE calibration across reps (with-Y vs no-Y ablation)
# ==============================================================================
cat(sprintf("\n=== BLOCK 4  SE calibration [HARD] : %d reps, N=%d, M=%d ===\n",
            REPS, N_REP, M_REP))

est_y <- se_y <- est_n <- se_n <- rep(NA_real_, REPS)
slr <- if (N_REP <= 100) c("SL.glm","SL.mean","SL.rpart","SL.earth") else sl_lib
t0 <- Sys.time()
for (b in 1:REPS) {
  dtb <- generate_data(N_REP, ERROR, OUTCOME, LATENT, INDIC)
  Pb  <- tryCatch(fit_sem_full(dtb), error = function(e) NULL)
  if (is.null(Pb)) next

  ry <- tryCatch(lv_mi_withY(dtb, Pb, slr, M_REP, V), error = function(e) NULL)
  if (!is.null(ry)) { est_y[b] <- ry$psi; se_y[b] <- ry$se }

  # no-Y ablation: the indicator-only posterior via 01's existing path.
  pf <- tryCatch(fit_sem_posterior(dtb, INDIC), error = function(e) NULL)
  if (!is.null(pf) && isTRUE(pf$converged)) {
    rn <- tryCatch(lv_tmle_mi(dtb, pf$mu, pf$post_sd, slr, DELTA, M_REP, V),
                   error = function(e) NULL)
    if (!is.null(rn)) { est_n[b] <- rn$psi; se_n[b] <- rn$se }
  }
  if (b %% 10 == 0 || b == 1)
    cat(sprintf("  rep %d/%d (%.1f min)\n", b, REPS,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cov_of <- function(e, s) mean(((e - 1.96*s) <= PSI) & (PSI <= (e + 1.96*s)), na.rm = TRUE)
report <- function(tag, e, s) {
  ok <- is.finite(e) & is.finite(s)
  mc_sd  <- sd(e[ok]); mean_se <- mean(s[ok])
  cat(sprintf("  %-7s  bias=%+.3f  MC-SD=%.3f  mean-SE=%.3f  SE/SD=%.2f  cov=%.2f  (n=%d)\n",
              tag, mean(e[ok]) - PSI, mc_sd, mean_se, mean_se / mc_sd,
              cov_of(e[ok], s[ok]), sum(ok)))
  c(bias = mean(e[ok]) - PSI, ratio = mean_se / mc_sd, cov = cov_of(e[ok], s[ok]))
}
ry_sum <- report("with-Y", est_y, se_y)
rn_sum <- report("no-Y",   est_n, se_n)

b4_se_ok  <- is.finite(ry_sum["ratio"]) && abs(ry_sum["ratio"] - 1) < 0.15
b4_cov_ok <- is.finite(ry_sum["cov"])   && ry_sum["cov"] > 0.90
b4_deatt  <- is.finite(ry_sum["bias"]) && is.finite(rn_sum["bias"]) &&
             abs(ry_sum["bias"]) < abs(rn_sum["bias"])
cat(sprintf("  [%s] with-Y Rubin SE tracks the MC-SD   [%s] coverage >= 0.90\n",
            if (b4_se_ok) "PASS" else "FAIL", if (b4_cov_ok) "PASS" else "FAIL"))
cat(sprintf("  [%s] Y-conditioning de-attenuates vs no-Y ablation\n",
            if (b4_deatt) "PASS" else "FAIL"))

# ==============================================================================
# VERDICT
# ==============================================================================
cat("\n================ EIF STRUCTURE VERDICT ================\n")
hard_pass <- b2_pass && b4_se_ok && b4_cov_ok
cat(sprintf("  Block 1 posterior calibration : %s\n", if (b1_pass) "PASS" else "FAIL"))
cat(sprintf("  Block 2 coarsening + LTV      : %s   [HARD]\n", if (b2_pass) "PASS" else "FAIL"))
cat(sprintf("  Block 4 SE / coverage         : %s   [HARD]\n",
            if (b4_se_ok && b4_cov_ok) "PASS" else "FAIL"))
cat(sprintf("\n  ==> %s\n", if (hard_pass && b1_pass)
    "STRUCTURE SURVIVES the falsifiers. Safe to derive the closed forms."
  else
    "STRUCTURE FAILED a falsifier. Do NOT proceed to the closed form -- diagnose first."))
cat("  Block 3 is diagnostic: the eta-gap sizes the estimation correction and so\n")
cat("  justifies (or undercuts) the 'estimated eta' modeling choice.\n")
