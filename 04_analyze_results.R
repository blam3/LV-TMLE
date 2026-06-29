# ==============================================================================
# FILE: 04_analyze_results.R  (REVISED)
# PURPOSE: Aggregate v2 results using each scenario's OWN true_psi, include the
#          Oracle ceiling, and summarize how much of the LV-TMLE SE comes from
#          latent uncertainty (between-imputation variance).
# ==============================================================================

# setwd("~/project_pi_ab2498/bnl24/lv_tmle_sim")   # adjust as needed

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
INPUT_DIR <- if (length(args) > 0) args[1] else "output"
ANALYSIS_DIR <- if (length(args) > 1) args[2] else "."
if (!dir.exists(ANALYSIS_DIR)) dir.create(ANALYSIS_DIR, recursive = TRUE)
analysis_file <- function(...) file.path(ANALYSIS_DIR, ...)

# 1. READ DATA -----------------------------------------------------------------
files <- list.files(INPUT_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No output files found.")
df <- do.call(rbind, lapply(files, readRDS))
cat("Reading", length(files), "result files from", INPUT_DIR, "\n")
cat("Writing analysis artifacts to", ANALYSIS_DIR, "\n")

# Re-attach grid factors (True_Psi already travels with each row).
if (file.exists("01_setup.R")) source("01_setup.R")  # gives sim_grid w/ factors
grid_cols <- c("Scenario_ID", "n_size", "meas_error", "outcome", "latent", "indicators")
df <- df %>% left_join(sim_grid[, grid_cols], by = "Scenario_ID")

# 2. FAILURE / DEGENERACY DIAGNOSTICS ------------------------------------------
fail_table <- df %>%
  filter(!is.na(Fail_Location) & Fail_Location != "None") %>%
  count(n_size, meas_error, outcome, latent, indicators,
        Fail_Location, Fail_Reason, name = "Count") %>%
  arrange(desc(Count))
write.csv(fail_table, analysis_file("failure_diagnostics.csv"), row.names = FALSE)
cat("Crash/non-convergence rows:", sum(fail_table$Count), "\n")

if (all(c("Warning_Count", "SL_Warning_Count", "SL_gam_Warning_Count",
          "SL_Calls", "SL_SuperLearner_Error", "SL_GLM_Fallback",
          "SL_Mean_Fallback", "SL_Predict_Fallback") %in% names(df))) {
  first_nonempty <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1] else NA_character_
  }
  learner_diag <- df %>%
    group_by(n_size, meas_error, outcome, latent, indicators) %>%
    summarise(
      Reps = n(),
      Reps_With_Warnings = sum(Warning_Count > 0, na.rm = TRUE),
      Total_Warnings = sum(Warning_Count, na.rm = TRUE),
      Total_SL_Warnings = sum(SL_Warning_Count, na.rm = TRUE),
      Total_SL_gam_Warnings = sum(SL_gam_Warning_Count, na.rm = TRUE),
      Mean_SL_Calls = mean(SL_Calls, na.rm = TRUE),
      Total_SL_SuperLearner_Errors = sum(SL_SuperLearner_Error, na.rm = TRUE),
      Total_SL_GLM_Fallbacks = sum(SL_GLM_Fallback, na.rm = TRUE),
      Total_SL_Mean_Fallbacks = sum(SL_Mean_Fallback, na.rm = TRUE),
      Total_SL_Predict_Fallbacks = sum(SL_Predict_Fallback, na.rm = TRUE),
      Warning_Example = first_nonempty(Warning_Examples),
      .groups = "drop"
    ) %>%
    arrange(desc(Total_SL_gam_Warnings), desc(Total_SL_Warnings))
  write.csv(learner_diag, analysis_file("learner_diagnostics.csv"), row.names = FALSE)
  cat("Learner warnings:", sum(learner_diag$Total_Warnings, na.rm = TRUE),
      "| SL.gam warnings:", sum(learner_diag$Total_SL_gam_Warnings, na.rm = TRUE),
      "\n")
} else {
  message("Skipping learner diagnostics: worker output predates warning/fallback columns.")
}

if (all(c("Var_Y_SEM", "Var_Y_Eff", "Mean_Slope_Eff", "SD_Slope_Eff",
          "Mean_Shift_Eff", "Mean_Mu_Shift", "SD_Mu_noY", "SD_Mu_Full",
          "Draw_Scale") %in% names(df))) {
  posterior_diag <- df %>%
    group_by(n_size, meas_error, outcome, latent, indicators) %>%
    summarise(
      Reps = n(),
      Mean_Var_Y_SEM = mean(Var_Y_SEM, na.rm = TRUE),
      Mean_Var_Y_Eff = mean(Var_Y_Eff, na.rm = TRUE),
      Mean_Slope_Eff = mean(Mean_Slope_Eff, na.rm = TRUE),
      Mean_SD_Slope_Eff = mean(SD_Slope_Eff, na.rm = TRUE),
      Mean_Shift_Eff = mean(Mean_Shift_Eff, na.rm = TRUE),
      Mean_Mu_Shift = mean(Mean_Mu_Shift, na.rm = TRUE),
      Mean_SD_Mu_noY = mean(SD_Mu_noY, na.rm = TRUE),
      Mean_SD_Mu_Full = mean(SD_Mu_Full, na.rm = TRUE),
      Mean_Draw_Scale = mean(Draw_Scale, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(n_size, outcome, latent, indicators, meas_error)
  write.csv(posterior_diag, analysis_file("posterior_update_diagnostics.csv"),
            row.names = FALSE)
  cat("Wrote posterior update diagnostics.\n")
} else {
  message("Skipping posterior update diagnostics: worker output predates columns.")
}

# 3. PERFORMANCE METRICS (vs per-row True_Psi) ---------------------------------
ok <- df %>% filter(is.na(Fail_Location) | Fail_Location == "None")

df_long <- ok %>%
  pivot_longer(c(SEM_Est, RegCal_Est, LV_Est, Oracle_Est, LV_noY_Est),
               names_to = "Estimator", values_to = "Estimate") %>%
  mutate(
    Estimator = sub("_Est$", "", Estimator),
    Coverage_Flag = case_when(             # renamed from the original 'Covariance'
      Estimator == "SEM"    ~ SEM_Cov,
      Estimator == "RegCal" ~ RegCal_Cov,
      Estimator == "LV"     ~ LV_Cov,
      Estimator == "Oracle" ~ Oracle_Cov,
      Estimator == "LV_noY" ~ LV_noY_Cov
    )
  )

df_summary <- df_long %>%
  group_by(n_size, meas_error, outcome, latent, indicators, Estimator) %>%
  summarise(
    Bias     = mean(Estimate - True_Psi, na.rm = TRUE),
    Abs_Bias = mean(abs(Estimate - True_Psi), na.rm = TRUE),
    RMSE     = sqrt(mean((Estimate - True_Psi)^2, na.rm = TRUE)),
    Coverage = mean(Coverage_Flag, na.rm = TRUE),
    N_reps   = sum(!is.na(Estimate)),       # actual successful reps (not /1000)
    .groups  = "drop"
  )
write.csv(df_summary, analysis_file("final_results_summary.csv"), row.names = FALSE)

# 3b. LV-TMLE variance decomposition: how much SE comes from latent uncertainty?
lv_var <- ok %>%
  group_by(n_size, meas_error, outcome, latent, indicators) %>%
  summarise(
    mean_within  = mean(LV_Wbar, na.rm = TRUE),   # TMLE EIF variance
    mean_between = mean(LV_B,    na.rm = TRUE),    # latent (between-imputation) variance
    latent_share = mean_between / (mean_within + mean_between),
    mean_within_noY  = mean(LV_noY_Wbar, na.rm = TRUE),
    mean_between_noY = mean(LV_noY_B,    na.rm = TRUE),
    latent_share_noY = mean_between_noY / (mean_within_noY + mean_between_noY),
    .groups = "drop"
  )
write.csv(lv_var, analysis_file("lv_variance_decomposition.csv"), row.names = FALSE)
# latent_share near 0 would mean the corrected SE still isn't propagating latent
# uncertainty (a red flag); a healthy nonzero share is the whole point of Goal 1.

# 4. PLOTS ---------------------------------------------------------------------
# Five factors don't fit one facet grid, so the default plots fix latent/indicators
# to the simplest cell and vary outcome x n_size. Change these filters to inspect
# other slices.
plot_slice <- df_summary %>% filter(latent == "normal", indicators == "continuous")

est_cols <- c(Oracle = "black", SEM = "#1b9e77", RegCal = "#d95f02",
              LV = "#7570b3", LV_noY = "#666666")

if (nrow(plot_slice) > 0) {
  p_bias <- ggplot(plot_slice,
                   aes(meas_error, Abs_Bias, color = Estimator, group = Estimator)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    facet_grid(outcome ~ n_size, labeller = label_both) +
    scale_color_manual(values = est_cols) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() + labs(title = "Absolute Bias (latent=normal, indicators=continuous)",
                      x = "Measurement error", y = "|Bias|")
  ggsave(analysis_file("plot_abs_bias.png"), p_bias, width = 10, height = 7)

  p_rmse <- ggplot(plot_slice,
                   aes(meas_error, RMSE, color = Estimator, group = Estimator)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    facet_grid(outcome ~ n_size, labeller = label_both) +
    scale_color_manual(values = est_cols) +
    theme_bw() + labs(title = "RMSE (latent=normal, indicators=continuous)",
                      x = "Measurement error", y = "RMSE")
  ggsave(analysis_file("plot_rmse.png"), p_rmse, width = 10, height = 7)

  p_cov <- ggplot(plot_slice,
                  aes(meas_error, Coverage, color = Estimator, group = Estimator)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    facet_grid(outcome ~ n_size, labeller = label_both) +
    scale_color_manual(values = est_cols) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() + labs(title = "95% CI Coverage (latent=normal, indicators=continuous)",
                      x = "Measurement error", y = "Coverage")
  ggsave(analysis_file("plot_coverage.png"), p_cov, width = 10, height = 7)
} else {
  message("Skipping default plots: no latent=normal, indicators=continuous rows.")
}

# A second coverage plot isolating the non-normality story (where the corrected
# Rubin SE is expected to help SEM-violating cells):
plot_latent <- df_summary %>% filter(outcome == "nonlinear", indicators == "continuous")
if (nrow(plot_latent) > 0) {
  p_cov_latent <- ggplot(plot_latent,
                         aes(meas_error, Coverage, color = Estimator, group = Estimator)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    facet_grid(latent ~ n_size, labeller = label_both) +
    scale_color_manual(values = est_cols) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() + labs(title = "Coverage by latent shape (outcome=nonlinear)",
                      x = "Measurement error", y = "Coverage")
  ggsave(analysis_file("plot_coverage_by_latent.png"), p_cov_latent, width = 10, height = 7)
} else {
  message("Skipping latent-shape coverage plot: no nonlinear continuous-indicator rows.")
}

cat("Done. Wrote final_results_summary.csv, lv_variance_decomposition.csv, and plots.\n")
