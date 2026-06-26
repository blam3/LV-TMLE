# ==============================================================================
# FILE: 02_run_sim.R  (REVISED)
# PURPOSE: Worker. Reads a scenario, runs N_REPS replications, records the v2
#          estimators (SEM, RegCal, outcome-aware LV-TMLE, Oracle, LV_noY) and computes coverage
#          against the SCENARIO'S OWN true_psi (no global TRUE_PSI any more).
# ==============================================================================

source("01_setup.R")
source("06_outcome_aware.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 02_run_sim.R [Scenario_ID] [N_Reps]")

SCENARIO_ID <- as.integer(args[1])
N_REPS      <- if (length(args) > 1) as.integer(args[2]) else 1000

scen     <- sim_grid[sim_grid$Scenario_ID == SCENARIO_ID, ]
if (nrow(scen) != 1) stop(sprintf("Scenario_ID %s not found exactly once.", SCENARIO_ID))
TRUE_PSI <- scen$true_psi          # per-scenario truth

if (!dir.exists("output")) dir.create("output")
output_file <- sprintf("output/res_scen_%02d.rds", SCENARIO_ID)

cat(sprintf("Scenario %d: N=%d | err=%s | outcome=%s | latent=%s | indic=%s | true_psi=%.3f\n",
            SCENARIO_ID, scen$n_size, scen$meas_error, scen$outcome,
            scen$latent, scen$indicators, TRUE_PSI))

# Coverage helper now takes the scenario truth explicitly.
safe_cov <- function(est, se, truth) {
  if (is.na(est) || is.na(se)) return(NA)
  (est - 1.96 * se) <= truth & truth <= (est + 1.96 * se)
}

results_list <- vector("list", N_REPS)

for (b in 1:N_REPS) {
  if (b %% 50 == 0 || b == 1) cat(sprintf("[Scen %d] rep %d / %d\n", SCENARIO_ID, b, N_REPS))

  dt_sim <- generate_data(scen$n_size, scen$meas_error,
                          scen$outcome, scen$latent, scen$indicators)

  res <- tryCatch(
    run_comparison_v2(dt_sim, indicators_type = scen$indicators, include_noY = TRUE),
    error = function(e) list(
      est_sem = NA, se_sem = NA, est_regcal = NA, se_regcal = NA,
      est_lv = NA, se_lv = NA, est_oracle = NA, se_oracle = NA,
      est_lv_noY = NA, se_lv_noY = NA,
      diagnostics = list(fail_location = "CRASH", fail_reason = conditionMessage(e),
                         lv_B = NA, lv_Wbar = NA, lv_m = NA,
                         lv_noY_B = NA, lv_noY_Wbar = NA, lv_noY_m = NA)))

  d <- res$diagnostics
  results_list[[b]] <- data.frame(
    Scenario_ID = SCENARIO_ID, Rep = b, True_Psi = TRUE_PSI,

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
    LV_noY_M_used = d$lv_noY_m
  )
}

final_df <- do.call(rbind, results_list)
saveRDS(final_df, file = output_file)
cat("Finished.\n")

# NOTE for 03_submit_job.sh: the grid now has nrow(sim_grid) scenarios
# (3 n x 3 error x 2 outcome x 2 latent x 2 indicators = 72). Update the array:
#     #SBATCH --array=1-72
# and bump --time / --mem: the LV-TMLE now runs M complete-data TMLE fits per rep.
