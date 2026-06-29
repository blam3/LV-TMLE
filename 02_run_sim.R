# ==============================================================================
# FILE: 02_run_sim.R  (REVISED)
# PURPOSE: Worker. Reads a scenario, runs N_REPS replications, records the v2
#          estimators (SEM, RegCal, outcome-aware LV-TMLE, Oracle, LV_noY) and computes coverage
#          against the SCENARIO'S OWN true_psi (no global TRUE_PSI any more).
# ==============================================================================

source("01_setup.R")
source("06_outcome_aware.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 02_run_sim.R [Scenario_ID] [N_Reps] [primary|exploratory] [Output_Dir]")
}

SCENARIO_ID <- as.integer(args[1])
N_REPS      <- if (length(args) > 1) as.integer(args[2]) else 1000
RUN_SCOPE   <- if (length(args) > 2) args[3] else "primary"
OUTPUT_DIR  <- if (length(args) > 3) args[4] else "output"
if (!RUN_SCOPE %in% c("primary", "exploratory")) {
  stop("Run scope must be 'primary' or 'exploratory'.")
}

scen     <- sim_grid[sim_grid$Scenario_ID == SCENARIO_ID, ]
if (nrow(scen) != 1) stop(sprintf("Scenario_ID %s not found exactly once.", SCENARIO_ID))
if (RUN_SCOPE == "primary" && scen$indicators != "continuous") {
  stop(sprintf(
    paste0("Scenario_ID %d uses ordinal indicators. Ordinal cells are currently ",
           "exploratory because pilot runs showed severe overshoot. Re-run with ",
           "third argument 'exploratory' only for diagnostic ordinal work."),
    SCENARIO_ID
  ))
}
TRUE_PSI <- scen$true_psi          # per-scenario truth

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
output_file <- file.path(OUTPUT_DIR, sprintf("res_scen_%02d.rds", SCENARIO_ID))

cat(sprintf("Scenario %d: N=%d | err=%s | outcome=%s | latent=%s | indic=%s | true_psi=%.3f\n",
            SCENARIO_ID, scen$n_size, scen$meas_error, scen$outcome,
            scen$latent, scen$indicators, TRUE_PSI))
cat(sprintf("Writing results to %s\n", output_file))

# Coverage helper now takes the scenario truth explicitly.
safe_cov <- function(est, se, truth) {
  if (is.na(est) || is.na(se)) return(NA)
  (est - 1.96 * se) <= truth & truth <= (est + 1.96 * se)
}

diag_value <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  as.numeric(x[1])
}

results_list <- vector("list", N_REPS)

for (b in 1:N_REPS) {
  if (b %% 50 == 0 || b == 1) cat(sprintf("[Scen %d] rep %d / %d\n", SCENARIO_ID, b, N_REPS))

  dt_sim <- generate_data(scen$n_size, scen$meas_error,
                          scen$outcome, scen$latent, scen$indicators)

  if (exists("lv_diag_reset")) lv_diag_reset()
  warning_messages <- character()
  res <- tryCatch(
    withCallingHandlers(
      run_comparison_v2(dt_sim, indicators_type = scen$indicators, include_noY = TRUE),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(
        est_sem = NA, se_sem = NA, est_regcal = NA, se_regcal = NA,
        est_lv = NA, se_lv = NA, est_oracle = NA, se_oracle = NA,
        est_lv_noY = NA, se_lv_noY = NA,
        diagnostics = list(fail_location = "CRASH", fail_reason = conditionMessage(e),
                           lv_B = NA, lv_Wbar = NA, lv_m = NA,
                           lv_noY_B = NA, lv_noY_Wbar = NA, lv_noY_m = NA,
                           var_y_sem = NA, var_y_eff = NA,
                           mean_slope_eff = NA, sd_slope_eff = NA,
                           mean_shift_eff = NA, mean_mu_shift = NA,
                           sd_mu_noY = NA, sd_mu_full = NA,
                           draw_scale = NA)
      )
    })
  learner_diag <- if (exists("lv_diag_snapshot")) lv_diag_snapshot() else list()
  warn_unique <- unique(warning_messages)

  d <- res$diagnostics
  results_list[[b]] <- data.frame(
    Scenario_ID = SCENARIO_ID, Rep = b, True_Psi = TRUE_PSI,
    Run_Scope = RUN_SCOPE,
    Primary_Analysis = scen$indicators == "continuous",

    SEM_Est    = res$est_sem,    SEM_Cov    = safe_cov(res$est_sem,    res$se_sem,    TRUE_PSI),
    RegCal_Est = res$est_regcal, RegCal_Cov = safe_cov(res$est_regcal, res$se_regcal, TRUE_PSI),
    LV_Est     = res$est_lv,     LV_Cov     = safe_cov(res$est_lv,     res$se_lv,     TRUE_PSI),
    Oracle_Est = res$est_oracle, Oracle_Cov = safe_cov(res$est_oracle, res$se_oracle, TRUE_PSI),
    LV_noY_Est = res$est_lv_noY, LV_noY_Cov = safe_cov(res$est_lv_noY, res$se_lv_noY, TRUE_PSI),

    # diagnostics
    Fail_Location = d$fail_location,
    Fail_Reason   = d$fail_reason,
    LV_B          = d$lv_B,        # between-imputation variance (latent uncertainty)
    LV_Wbar       = d$lv_Wbar,     # within-imputation variance (TMLE EIF variance)
    LV_M_used     = d$lv_m,
    LV_noY_B      = d$lv_noY_B,
    LV_noY_Wbar   = d$lv_noY_Wbar,
    LV_noY_M_used = d$lv_noY_m,
    Var_Y_SEM      = diag_value(d$var_y_sem),
    Var_Y_Eff      = diag_value(d$var_y_eff),
    Mean_Slope_Eff = diag_value(d$mean_slope_eff),
    SD_Slope_Eff   = diag_value(d$sd_slope_eff),
    Mean_Shift_Eff = diag_value(d$mean_shift_eff),
    Mean_Mu_Shift  = diag_value(d$mean_mu_shift),
    SD_Mu_noY      = diag_value(d$sd_mu_noY),
    SD_Mu_Full     = diag_value(d$sd_mu_full),
    Draw_Scale     = diag_value(d$draw_scale),
    Warning_Count = length(warning_messages),
    SL_Warning_Count = sum(grepl("SuperLearner|algorithm SL\\.", warning_messages)),
    SL_gam_Warning_Count = sum(grepl("SL\\.gam", warning_messages)),
    Warning_Examples = substr(paste(head(warn_unique, 3), collapse = " || "), 1, 500),
    SL_Calls = learner_diag$sl_calls,
    SL_Fallback_Short_N = learner_diag$sl_fallback_short_n,
    SL_Fallback_Low_Variance = learner_diag$sl_fallback_low_variance,
    SL_SuperLearner_Error = learner_diag$sl_superlearner_error,
    SL_GLM_Fallback = learner_diag$sl_glm_fallback,
    SL_Mean_Fallback = learner_diag$sl_mean_fallback,
    SL_Discrete = learner_diag$sl_discrete,
    SL_Predict_Fallback = learner_diag$sl_predict_fallback
  )
}

final_df <- do.call(rbind, results_list)
saveRDS(final_df, file = output_file)
cat("Finished.\n")

# NOTE for the primary full run: use continuous-indicator scenarios only
# (Scenario_ID 1-36). Ordinal scenarios remain exploratory until the ordinal
# posterior approximation is fixed or validated. To run an ordinal diagnostic:
#     Rscript --vanilla 02_run_sim.R 50 5 exploratory
