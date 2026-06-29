# LV-TMLE Project — Handoff Document

**Project:** Latent-Variable Targeted Maximum Likelihood Estimation (LV-TMLE) — propagating latent-variable (SEM) measurement uncertainty through a targeted causal estimator for a continuous latent exposure.
**Estimand:** Modified Treatment Policy (MTP) shift contrast, ψ = E[Y(L+δ)] − E[Y(L)], with δ = 1.
**Status:** Infrastructure rebuilt and validated by smoke test. The outcome-aware estimator (the statistical fix) is written and the `05b` smoke-test harness now runs. The outcome-link step in `sem_posterior_full()` has been replaced with a measurement-error-corrected quadratic moment link. A larger 40-rep `05b` validation pilot passed the hard checks. `run_comparison_v2` is wired into the worker/analysis path with learner warning/fallback diagnostics. A broader six-cell worker pilot completed with no crashes, but it surfaced repeated `SL.gam` warnings and severe ordinal-cell overshoot. Decision made: primary scaling is continuous-indicator only; ordinal cells remain exploratory/diagnostic until fixed or validated; `SL.gam` has been dropped from the active learner library.
**Audience:** A collaborator or a future session picking this up cold.

---

## 1. Executive summary

The original implementation demonstrated the *opposite* of its goal: the proposed LV-TMLE lost to both competitors on bias, RMSE, and coverage in essentially every cell. We traced this to specific, fixable causes — not bad luck. The project has since been rebuilt: a fair, SEM-adversarial set of data-generating processes; a single correct cross-fitted shift-TMLE engine shared by all estimators; per-scenario truth computed by Monte Carlo; an oracle ceiling; and a corrected variance via Rubin's rules. A single-cell smoke test confirmed the **infrastructure is sound** (oracle recovers the truth, the between-imputation variance is positive, the SEM shows the expected curvature bias) and confirmed the **last remaining method-level bug**: imputing the latent from the indicators alone re-introduces attenuation. The fix — imputing the latent conditional on the outcome — is implemented in `06_outcome_aware.R` and awaits validation.

The most important strategic update: a simple regression-calibration plug-in already nearly recovers the point estimate in the original smoke test (bias −0.09, coverage 0.93 in the test cell). So **the paper's contribution should be reframed around principled, robust *inference* (coverage that holds where simpler plug-ins' standard errors break), not bias/RMSE.** The first outcome-aware checks improve the no-Y MI ablation but do not yet eliminate attenuation.

---

## 2. Goals and the (reframed) success criterion

From `LV-TMLE-goals.md`:

1. Integrate SEM latent variables with TMLE.
2. Conduct a simulation comparing LV-TMLE with suitable alternatives.
3. Demonstrate LV-TMLE is viable — better than the alternatives in *most* (not all) conditions — in realistic conditions for psychologists.

**Honest reframing of Goal 3.** Against a correctly specified SEM on its home turf (Gaussian indicators, linear L→Y, normal latent), no semiparametric method should beat the SEM on efficiency — that is a theorem, not a tuning failure. So "viability" must be demonstrated where the SEM's assumptions break, and the natural axis of advantage is **valid inference (coverage)**, since a regression-calibration plug-in already handles the point estimate. The simulation now includes both SEM-violating cells (where LV-TMLE should win) and an honest negative-control cell (linear/normal/continuous, where the SEM should win) — showing the latter is evidence the study is trustworthy.

---

## 3. Current state at a glance

| Component | Status |
|---|---|
| Critique of original project | Complete (`lv_tmle_critique.md`) |
| Adversarial-to-SEM DGPs + Monte Carlo true-ψ | Complete (`00_dgp_variants.R`) |
| Shared cross-fitted shift-TMLE engine + correct EIF | Complete (`01_setup.R`) |
| Naive + Oracle estimators on shared engine | Complete (`01_setup.R`) |
| Per-scenario true_psi threaded through pipeline | Complete (`01/02/04`) |
| Corrected LV-TMLE via Rubin's rules (engine) | Complete (`01_setup.R`, `lv_tmle_mi`) |
| Single-cell smoke test | Run; results in hand (`05_smoke_test.R`) |
| **Outcome-aware imputation (the point-estimate fix)** | **Written; measurement-error-corrected outcome link passes the 40-rep `05b` pilot (`06_outcome_aware.R`)** |
| Outcome-aware smoke test (`05b`) | **Written and runnable; latest larger pilot passes hard checks** |
| `run_comparison_v2` wired into worker/analysis | **Complete; broader worker pilot passed execution checks and includes warning/fallback diagnostics** |
| Primary continuous-indicator array run | Not started |
| Ordinal-indicator scenarios | Exploratory only; blocked from default primary worker runs |
| Results report | Stale — describes the OLD broken pipeline |

---

## 4. File inventory

Two locations matter. The **original** project files live in `/mnt/project/` (the pre-critique versions). The **current working set** is in the outputs directory and supersedes them.

| File | Role | Status / notes |
|---|---|---|
| `LV-TMLE-goals.md` | The three project goals | Unchanged source of truth |
| `lv_tmle_critique.md` | Brutally honest critique mapping the work to the goals | Complete; the diagnostic that drove the rebuild |
| `00_dgp_variants.R` | Modular DGPs (outcome × latent × indicators) + `compute_true_psi`, `build_psi_lookup`, `build_sim_grid`, `attach_true_psi`, `sem_linear_target` | Complete. Single source of truth per DGP piece prevents estimand drift |
| `01_setup.R` (revised) | `shift_tmle` engine (correct EIF), `safe_SL`/`robust_predict`, `fit_sem_posterior`, `lv_tmle_mi` (Rubin), `run_comparison`; loads cached `psi_lookup.rds`, builds `sim_grid` with per-row `true_psi` | Complete; **infra validated**. Its `lv_tmle_mi` here draws from the indicator-only posterior (Bartlett mean) — the **attenuated** version, superseded for the point estimate by `06` |
| `02_run_sim.R` (revised) | Worker; passes DGP factors to `generate_data`, records SEM/RegCal/outcome-aware LV-TMLE/Oracle/LV_noY, computes coverage vs per-scenario `true_psi`, preserves SuperLearner warning/fallback diagnostics, and defaults to continuous-only `primary` scope | Complete for v2 primary runs; ordinal diagnostics require explicit `exploratory` scope |
| `04_analyze_results.R` (revised) | Aggregates vs per-row `True_Psi`, includes Oracle and LV_noY, writes LV variance decomposition plus learner diagnostics | Complete for v2 columns |
| `05_smoke_test.R` | Single-cell diagnostic harness | Run; results captured below |
| `05b_smoke_test_aware.R` | Outcome-aware single-cell diagnostic harness | Written and runnable. Use `Rscript --vanilla -e 'options(lv_tmle.reps=10, lv_tmle.M=10, lv_tmle.V=3); source("05b_smoke_test_aware.R")'` for quick checks; avoid env-prefixed `Rscript` invocations in this managed shell because they can break `parallel::detectCores()` and therefore lavaan. |
| `05e_observed_eif_smoke.R` | Observed-data latent-variable EIF diagnostic harness | Written and rerun on 2026-06-29 for scenario 14. Computes posterior moments, projected observed scores, empirical observed information, `dot_psi_eta`, and `D_eff`; diagnostic only, not production estimation |
| `06_outcome_aware.R` | `sem_posterior_full` (quadrature L\|I,Y,W,Z posterior with measurement-error-corrected quadratic outcome link) + `run_comparison_v2` (SEM, RegCal, outcome-aware LV-TMLE, Oracle, optional MI-no-Y ablation) | **Operational; passed 40-rep `05b` pilot, ready for worker wiring/pilot array** |
| `report_lv_tmle.md` | Original results report | **STALE.** Describes the broken pipeline's numbers and "LV-TMLE loses" conclusion. Do not cite; rewrite after the corrected run |
| `03_submit_job.sh` | SLURM array submitter | Needs `--array=1-72` and larger `--time`/`--mem` (see §9) |

---

## 5. The estimators (what each one is for)

The corrected design compares, on a single shared TMLE engine:

- **SEM** — `lavaan` CFA with structural paths; the L→Y coefficient is its estimand. The parametric benchmark. Under the nonlinear DGP it targets the linear coefficient (≈ 2.0) and is **structurally blind to the curvature** (true ψ ≈ 2.7).
- **Regression calibration (RegCal)** — plug E[L | I, W, Z] (the no-outcome posterior mean / EAP scores) into the shift-TMLE. The standard, honest "ignore-uncertainty" baseline. De-attenuates the point estimate but its SE does not formally carry latent uncertainty.
- **LV-TMLE (proposed)** — draw L from its posterior given **I, Y, W, Z**, fit the shift-TMLE on each draw, combine by Rubin's rules. De-attenuates *and* propagates latent uncertainty into the SE via the between-imputation variance.
- **Oracle** — true `L_true` plugged into the same engine. The performance ceiling (no measurement error); included only as a benchmark and ignored by all real estimators.
- **MI-without-Y (ablation)** — the cautionary version that imputes from indicators only. Kept for one figure to show *why* the outcome must enter the imputation.

---

## 6. Completed work

The project moved through a critique-then-rebuild arc:

The **critique** established that the original LV-TMLE inverted its own idea at three points and was tested on a field tilted toward the SEM. Building on that, the **DGP layer** was rebuilt so the estimand is computed by Monte Carlo from a single shared mean function and cached, guaranteeing the truth cannot drift when the model leaves the linear world; the nonlinear-in-L outcome creates a deliberate ≈ 0.7 gap between the true shift effect (≈ 2.7) and the SEM's linear target (≈ 2.0), which is the scientific signal. The **estimation layer** was rebuilt around one cross-fitted shift-TMLE engine carrying the full efficient influence function (including the previously dropped `H·(Y−Q)` score term), with the naive and oracle estimators differing from each other only in which L-values they receive. The broken "stack 100 draws and fit once" step was replaced by proper multiple imputation: M complete n-row datasets, each fit separately, combined by Rubin's rules so the **between-imputation variance B** carries the latent uncertainty the old SE discarded. Per-scenario truth was threaded through the worker and the analysis, the `medium ≡ none` duplicate was removed, and the grid became a clean factorial with an oracle ceiling and a negative-control cell.

The **smoke test** then validated the rebuilt infrastructure and surfaced the final method-level bug, after which the **outcome-aware imputation** was implemented to address it.

The first `05b` outcome-aware runs show that the v2 path is operational but not yet solved. With REPS=10, M=10, V=3, draw_scale=0.5, the hard infrastructure checks passed, but LV-TMLE bias was −0.106 with MCSE 0.041, narrowly failing the hard bias check (tolerance 0.102). The no-Y ablation remained much more attenuated (bias −0.407), so conditioning on Y helps. Increasing draw_scale to 1.0 raised the latent variance share (about 25%) but worsened LV bias to −0.193, so the residual bias is not simply caused by posterior draw shrinkage.

A follow-up posterior-only diagnostic (`05d_posterior_gap_diagnostics.R`) localized the gap. In the nonlinear/normal/continuous, N=250, medium-error cell over 100 reps, the current outcome-aware posterior improved latent recovery only modestly (`RMSE(L_hat, L_true)` about 0.333 for no-Y scores vs 0.309 for current full posterior). Oracle outcome-link falsifiers improved to about 0.28, so the posterior machinery can do better when the outcome likelihood is congenial. The fitted outcome link is attenuated: its average one-unit shift was about 2.58 while the true shift at the factor-score scale was about 2.70. The flexible residual-variance correction is not the primary problem: using raw residual variance instead of the corrected/floored variance gave essentially the same posterior RMSE. Explicit quadratic score-link variants raised the shift modestly (about 2.63) but did not improve posterior recovery. Diagnosis: the current `sem_posterior_full()` still uses an outcome likelihood learned on noisy factor scores, so Y-conditioning moves the posterior in the right direction but not enough; this is a congeniality / measurement-error-in-outcome-link issue rather than a Rubin variance or draw-scale issue.

The outcome-link step was then replaced with a measurement-error-corrected quadratic moment model: fit the outcome regression using `E[L|I]`, `E[L^2|I]`, and `E[LW|I]`, then evaluate the likelihood at candidate latent values `L`. After this replacement, `05d` over 100 reps showed mean shift about 2.63 and current posterior RMSE about 0.311 (still not close to the oracle-link 0.28 bound), but the estimator-level `05b` pilot improved materially. With REPS=10, M=10, V=3, draw_scale=0.5, LV-TMLE bias was −0.042 with MCSE 0.037 and passed the hard no-overshoot/on-target check; LV_noY remained attenuated at about −0.407. Interpretation: the replacement is useful at the estimator level, but the posterior-only evidence is not strong enough to justify full scaling without a larger pilot.

The larger `05b` validation pilot (REPS=40, M=20, V=3, draw_scale=0.5) passed all hard checks. True ψ was 2.6997. LV-TMLE mean was 2.65 (bias −0.052, MCSE 0.033, coverage 0.95, SE/SD 0.87); Oracle mean was 2.65 (bias −0.051, coverage 0.875); RegCal remained attenuated (bias −0.195, coverage 0.825); LV_noY remained strongly attenuated (bias −0.430, coverage 0.70); SEM remained curvature-biased (bias −0.774, coverage 0.075). Rubin between-imputation variance was positive (B ≈ 0.00518; latent share ≈ 14.4%). The pilot produced repeated `SL.gam` warnings (`Error in algorithm SL.gam; removed from SuperLearner`), but the wrapper continued and 100% of reps produced SEM/RegCal/LV/Oracle estimates. Before the full array, either log learner fallback rates explicitly or consider dropping/replacing `SL.gam` if warnings become too noisy.

A tiny end-to-end worker pilot was then run after archiving existing outputs to `output/archive_20260627_0506_pre_pilot/`. The fresh pilot used scenario 1 for 2 reps (linear/normal/N=50/small error/continuous), scenario 14 for 1 rep (nonlinear/normal/N=250/medium error/continuous), and scenario 50 for 1 rep (nonlinear/normal/N=250/medium error/ordinal), all through `02_run_sim.R` with the production default `M=25` and `include_noY=TRUE`. `04_analyze_results.R` completed, wrote `learner_diagnostics.csv`, and reported 0 crash/non-convergence rows, 0 learner warnings, 0 `SL.gam` warnings, and 0 SuperLearner/GLM/mean/prediction fallbacks. Mean `SL_Calls` was 520 per replication, confirming the diagnostic counters are active. Rubin between-imputation variance was nonzero in all three pilot cells. The one-rep ordinal cell showed outcome-aware LV overshoot in that draw, so it should be treated only as a plumbing result, not as evidence of ordinal performance.

A broader pilot was then run after archiving the tiny-pilot outputs to `output/archive_20260627_0515_tiny_pilot/`. The fresh run used 5 reps each for scenarios 14, 17, 32, 50, 53, and 68: nonlinear/N=250 cells crossing continuous vs ordinal indicators, normal vs skew latents, and medium/large measurement error. All six workers completed and `04_analyze_results.R` reported 0 crash/non-convergence rows. The new learner diagnostics captured 102 total warnings, 100 of them `SL.gam` removals, all in scenario 14 (medium-error/normal/continuous); there were still 0 counted SuperLearner/GLM/mean/prediction fallbacks. Continuous cells were roughly plausible at this small `N_reps=5` scale: LV bias was about -0.09 in scenario 14, +0.10 in scenario 17, and +0.15 in scenario 32, with nonzero Rubin between-imputation variance. Ordinal nonlinear cells were not ready for scaling: LV bias was about +0.35 in scenario 50, +1.81 in scenario 53, and +1.12 in scenario 68; RegCal and SEM also overshot badly in large-error/skew ordinal cells. Treat this as a blocking pilot result for ordinal-branch claims and full-array launch, not as a final performance estimate.

Follow-up decision: ordinal cells stay in the repository for diagnostics, but they are excluded from the primary scaling path and should not support manuscript claims until the ordinal posterior approximation is fixed or externally validated. `02_run_sim.R` now defaults to `primary` scope and refuses ordinal scenarios unless called with a third argument, e.g. `Rscript --vanilla 02_run_sim.R 50 5 exploratory`. The active learner library is centralized in `lv_sl_library()` and excludes `SL.gam`; a quick in-memory nonlinear/continuous check after the change produced 0 warnings and 0 `SL.gam` references.

The post-`SL.gam` continuous-only pilot was rerun on 2026-06-28 in an isolated output directory, `output/pilot_20260628_no_gam_continuous/`, using 5 reps each for scenarios 14, 17, and 32. Before the rerun, `02_run_sim.R` and `04_analyze_results.R` were made path-aware so pilots can write to separate directories without overwriting the main `output/` files or root-level summaries. The rerun also exposed a managed-shell issue with lavaan 0.6-21: `parallel::detectCores()` returns `NA`, and lavaan errors while validating `ncpus`. `01_setup.R` now wraps lavaan SEM calls with `with_lavaan_core_fallback()`, temporarily forcing `detectCores()` to `1L` only around the SEM fit when needed. With that fix, the 15-rep pilot completed with 0 crashes, 0 learner warnings, 0 `SL.gam` warnings, and 0 SuperLearner/GLM/mean/prediction fallbacks. Mean LV estimates were close to truth in the tiny pilot: scenario 14 bias about -0.020, scenario 17 bias about +0.043, scenario 32 bias about +0.038. Rubin between-imputation variance remained positive in all three cells, with latent shares around 8.6%-12.5%. Treat this as a plumbing and warning-cleanliness check, not a performance estimate.

A medium continuous-only pilot was then run on 2026-06-28 in `output/pilot_20260628_medium_continuous_50reps/`, using 50 reps each for scenarios 14, 17, and 32. It completed cleanly: 150 rows, 0 crashes/non-convergence rows, 0 learner warnings, 0 `SL.gam` warnings, and 0 SuperLearner/GLM/mean/prediction fallbacks. This validates the worker/analysis plumbing for a medium isolated pilot, but the estimator results are mixed enough to diagnose before primary scaling. Scenario 14 (normal latent, medium error) looked stable: LV bias about -0.017, RMSE 0.184, coverage 0.88; Oracle bias about -0.043, coverage 0.94; SEM retained the expected curvature bias (-0.685, coverage 0.06). Scenario 17 (normal latent, large error) showed LV positive bias about +0.171 and low coverage 0.50, while Oracle was near unbiased (bias +0.006, coverage 0.92) and LV_noY was strongly attenuated (bias -0.693). Scenario 32 (skew latent, medium error) also showed LV positive bias about +0.158 and coverage 0.82, with Oracle bias -0.053 and SEM overshoot/bad coverage. Rubin between-imputation variance was positive in all cells, with latent shares about 8.6%-15.9%. Interpretation: the continuous pipeline is now operational and warning-clean, but the outcome-aware posterior/congeniality gap still appears under large measurement error and skew latent stress cells; do not launch the full primary array until this is understood.

On 2026-06-29, `05d_posterior_gap_diagnostics.R` was made scenario-aware and rerun for scenarios 14, 17, and 32 with 15 reps each, writing `posterior_gap_diagnostics_scen_14_17_32.rds`. The diagnostic reproduced the pattern that scenario 14 is comparatively stable and localized different stress modes for 17 and 32. In scenario 17 (normal latent, large measurement error), the current outcome link under-shifted relative to the oracle link (`current_shift` about 2.56 vs true shift about 2.64 at the factor-score scale), and the oracle outcome link improved posterior RMSE materially (current full posterior RMSE about 0.517 vs oracle true-link RMSE about 0.433). Raw-variance/flexible links were slightly better than the current corrected-variance link in this cell, so residual-variance weighting may contribute, but the larger issue remains outcome-link congeniality under noisy factor scores. In scenario 32 (skew latent, medium error), the current link over-shifted (`current_shift` about 2.86 vs true shift about 2.70), and oracle-link posterior recovery again improved RMSE (about 0.279 vs current about 0.313). Interpretation: the overshoot is not a single generic "too much Y-conditioning" problem; large measurement error and latent skew push the fitted outcome link in opposite directions. `02_run_sim.R` now records posterior-update diagnostics (`Var_Y_SEM`, `Var_Y_Eff`, `Mean_Shift_Eff`, `Mean_Mu_Shift`, etc.), and `04_analyze_results.R` writes `posterior_update_diagnostics.csv` when those columns are present. A one-rep isolated smoke run in `output/pilot_20260629_diag_column_smoke/` verified the new worker/analysis columns.

Later on 2026-06-29, conservative posterior-update variants were added to `05d_posterior_gap_diagnostics.R` only: ME raw-variance weighting, partial Y-likelihood powers (`0.5`, `0.25`), and deterministic-fold cross-fitted ME links with corrected or raw variance. The 15-rep scenario 14/17/32 run wrote `posterior_gap_diagnostics_variants_scen_14_17_32.rds`. The best conservative variants improved posterior RMSE modestly without closing the oracle-link gap: scenario 14 current RMSE 0.320 vs best conservative 0.316 (`me_raw_var`), scenario 17 current 0.517 vs best 0.497 (tempered current link), and scenario 32 current 0.313 vs best 0.305 (tempered current link). Cross-fitting alone did not solve the issue. The practical next candidate is a half-strength outcome likelihood (`current_y_power_0.5`), but it should be tested at the estimator level in an isolated pilot before changing `06_outcome_aware.R`.

`EIF_AUDIT_NOTES.md` was also added. Current conclusion: the complete-data shift-TMLE EIF in `shift_tmle()` is structurally correct for the oracle/full-data shifted contrast, including the `H * (Y - Q)` term and the subtraction of the `E[Y]` EIF. However, the project does not yet have a complete observed-data EIF for the latent-variable estimator that accounts for all estimated posterior/link nuisance parameters. Mean and variance parameters are not passive: `post_sd_noY`, indicator residual variances, latent variance, `var_y_eff`, and outcome residual variance all need nuisance-score or cross-fitting arguments before making manuscript-level "efficient EIF" claims. Rubin's rules carry imputation variability, but they are not a derivation of the influence contribution of estimated posterior variance parameters.

`OBSERVED_DATA_EIF_DERIVATION.md` now gives the finite-dimensional observed-data derivation for the continuous Gaussian latent model. Key correction: the efficient observed-data IF is not simply a posterior average `E[D_full | O]` once posterior/model parameters are estimated. With nonparametric `P_X` and finite-dimensional conditional model `eta`, it is

`D_eff(O) = xi_eta(X) - psi + dot_psi_eta^T I_eta^{-1} S_eta^O(O)`,

where `S_eta^O(O) = E[S_eta^F | O]`. The note derives projected scores for latent mean, latent variance, loadings, indicator residual variances, outcome mean parameters, and outcome residual variance. It also treats `post_sd_noY` by delta method and flags the nonregular floor in `var_y_eff`. Cross-fitting helps with empirical-process control and EIF evaluation, but it does not remove the first-order nuisance-score term from a plug-in estimator when the target depends on estimated latent-model parameters.

`05e_observed_eif_smoke.R` now implements a finite-sample smoke harness for that derivation in the continuous Gaussian working model. The default rerun on 2026-06-29 used scenario 14 with N=750 and seed 20260629, writing `observed_eif_smoke.rds`. The run produced `psi_hat_working` 2.654 versus true ψ 2.700, posterior z mean 0.178 and posterior z SD 0.920, centered `D_eff` mean about 4.7e-16, and `D_eff` SD about 3.84. The observed information was near-singular (`min eigen` about 1.7e-16; condition number about 2.8e16), so the harness used a small ridge (`1.375e-06`) before solving. The largest covariance-identity gap was about 0.110, led by `beta_bZ2`, `alpha_Z`, and `nu1`; treat this as a diagnostic stability flag rather than production inference.

---

## 7. Crucial decisions (with rationale)

1. **Estimand is the MTP shift contrast**, not a regression coefficient. It stays well-defined under nonlinearity and is comparable across all estimators.
2. **Compute true ψ by Monte Carlo from one shared mean function, cached to `psi_lookup.rds`.** Prevents the estimand from silently diverging from the simulated data once the outcome is nonlinear. (Old code hard-coded `TRUE_PSI = 2.0` in two places — fatal once ψ varies by cell.)
3. **Per-draw Rubin combination, never stacking.** Stacking M draws into one dataset induces regression dilution on the L→Y slope (attenuation that *grows* with measurement error). Each draw must be a complete n-row dataset.
4. **One shared cross-fitted engine with the complete EIF.** Guarantees apples-to-apples comparison and a correct variance; restoring the `H·(Y−Q)` term was necessary for valid coverage. Un-stacking also auto-fixed the inner-CV leakage and the small-sample learner switch.
5. **Add an Oracle (true L) ceiling.** Distinguishes irreducible measurement error from method inefficiency.
6. **DGP must be adversarial to SEM.** Nonlinear-in-L outcome (curvature gap), non-normal latent, and ordinal indicators are the SEM-breaking levers; keep linear/normal/continuous as the honest negative control where the SEM should win.
7. **Use regression/EAP factor scores, not Bartlett, for the plug-in.** Bartlett scores are unbiased *for L* but inflate the predictor's variance, so they **attenuate** when L is used as a regressor. The smoke test made this stark (Bartlett-naive bias −0.37 vs EAP-RegCal −0.09).
8. **Impute the latent conditional on the outcome.** Imputing the exposure from indicators only makes the draw noise independent of Y, re-introducing attenuation (the "plausible values must include the analysis variables" rule from educational measurement). This is the core of `06`.
9. **Drop the `medium` misspecification duplicate; separate functional-form from confounder issues.** The old "misspecification" factor was one-third duplicate and conflated two distinct things; the DGP features now carry the SEM stress cleanly.
10. **Reframe the contribution toward robust inference.** Since regression calibration already nearly fixes the point estimate, LV-TMLE's distinctive value is the principled, latent-uncertainty-propagating SE (the B term) that should keep coverage where simpler plug-ins' informal SEs fail.

---

## 8. Key empirical findings (smoke test)

Cell: outcome = nonlinear, latent = normal, indicators = continuous, N = 250, error = medium, 40 reps. **True ψ = 2.700.**

| Estimator | Mean | Bias | SE/SD ratio | Coverage |
|---|---:|---:|---:|---:|
| Oracle (true L) | 2.69 | −0.008 | 1.04 | 0.90 |
| RegCal (EAP scores) | 2.61 | −0.090 | 1.03 | 0.93 |
| Naive (Bartlett) | 2.33 | −0.374 | 0.84 | 0.48 |
| Robust (MI **without** Y) | 2.12 | −0.585 | 1.47 | 0.38 |
| SEM | 2.02 | −0.683 | 0.71 | 0.08 |

LV-TMLE variance decomposition: within-imputation 0.0537, **between-imputation B = 0.0210 (latent share ≈ 28%)** — the SE is carrying real latent uncertainty.

**Reading:** Infrastructure is sound (Oracle recovers ψ; B > 0; SEM shows the predicted −0.68 curvature gap). The MI-without-Y "Robust" is *worse* than naive — exactly the attenuation predicted — which is what `06` exists to fix. RegCal shows the point estimate is largely recoverable with the right factor scores.

---

## 9. Open risks and questions

- **`06` is operational and passed a 40-rep validation pilot.** The measurement-error-corrected outcome link fixed the `05b` hard-check failure, but `05d` still shows a posterior-recovery gap relative to oracle outcome links. Overshoot has not appeared in validated continuous cells. `SL.gam` has been dropped from the active learner library after repeated pilot removals.
- **Robustness tension.** Conditioning the imputation on Y *via the parametric SEM* reintroduces dependence on the SEM's outcome model — partly in tension with the semiparametric-robustness selling point. The principled endpoint is the efficient influence function that integrates over the latent without leaning on the SEM's Y-model; that derivation is **not done** and is the deeper "technical heart" of the paper.
- **Ordinal cells are exploratory only.** The Gaussian conjugate posterior in `06` is exact only for continuous indicators; the ordinal branch falls back to a `lavPredict`-based approximation and the broader pilot showed severe overshoot. Validate or replace this branch before using ordinal cells for primary claims.
- **Compute cost.** LV-TMLE is now M complete-data TMLE fits per replication; the primary continuous grid is 36 scenarios × ~1000 reps × M. Budget runtime/memory before launching.
- **Degenerate-fit reporting.** "Zero crashes" can still hide runs where every learner fell back to the mean. Logging fallback activation rates (critique §2.9) is still outstanding.

---

## 10. Ongoing tasks (prioritized)

**Immediate**
1. **Estimator-level test of the tempered posterior candidate.** Keep production unchanged until tested. In a diagnostic branch or temporary harness, compare the current outcome-aware posterior against `current_y_power_0.5` for scenarios 14/17/32 using the full `run_comparison_v2` estimator path, then decide whether a production option is justified.
2. **EIF implementation work.** The first observed-data EIF smoke harness now exists in `05e_observed_eif_smoke.R`. Next, stress-test the harness across scenarios 14/17/32 and larger N, inspect the near-singular observed-information directions, and decide whether the parameterization needs constraints/reparameterization before any manuscript-level efficiency claim.
3. After the posterior-update diagnosis, run a broader continuous-only pilot across additional representative scenarios before the full array. Keep writing to isolated output directories, and use the new `posterior_update_diagnostics.csv` output.
4. Prepare the primary continuous-only scaling path: scenario IDs 1-36 in `primary` scope. Do not include ordinal IDs 37-72 in the primary array.
5. Keep the posterior-only `05d` gap and ordinal overshoot in method notes: the measurement-error-corrected link improves estimator behavior in continuous cells but does not validate the ordinal approximation.

**Before the full run**
3. **Precompute `psi_lookup.rds` once** (`saveRDS(build_psi_lookup(n_mc = 2e7), "psi_lookup.rds")`) to avoid a write race across workers.
4. **Update `03_submit_job.sh`**: `--array=1-72`, larger `--time`/`--mem`.
5. **Partial continuous-only pilot run** (a handful of scenarios, ~200 reps) to confirm runtime and cross-cell behavior before committing the primary array.

**Full study and write-up**
6. Run the primary 36 × 1000 continuous-indicator array; regenerate plots and the variance decomposition.
7. **Rewrite `report_lv_tmle.md`** around the corrected results and the reframed (inference-centric) contribution.

**Open / research**
8. Quantify the congeniality residual bias vs the Oracle across measurement-error levels.
9. Validate or bound the ordinal-indicator approximation.
10. Consider deriving the latent-shift EIF (the principled inference endpoint).
11. Log and report fallback/degenerate-fit rates.

---

## 11. How to reproduce / run

```r
# 0. One-time: precompute the true-psi cache (publication precision)
source("00_dgp_variants.R")
saveRDS(build_psi_lookup(n_mc = 2e7), "psi_lookup.rds")

# 1. Sanity-check the DGPs (optional)
RUN_DGP_DEMO <- TRUE; source("00_dgp_variants.R")   # prints true psi by cell

# 2. Infrastructure smoke test (already run; ~10–20 min)
Rscript 05_smoke_test.R

# 3. Outcome-aware smoke test
Rscript --vanilla -e 'options(lv_tmle.reps=10, lv_tmle.M=10, lv_tmle.V=3); source("05b_smoke_test_aware.R")'

# 4. Small primary continuous-only pilot before the full array
Rscript --vanilla 02_run_sim.R 14 5 primary
Rscript 04_analyze_results.R

# 5. Optional ordinal diagnostic only; not part of the primary array
Rscript --vanilla 02_run_sim.R 50 5 exploratory
```

Dependencies: R with `lavaan`, `data.table`, `SuperLearner`, `earth`, `rpart`, `dplyr`, `ggplot2`, `tidyr`.

---

## 12. Strategic framing for the paper

The story is no longer "our estimator beats naive on bias." Regression calibration with the right factor scores already does that. The defensible, novel contribution is:

> A targeted, latent-uncertainty-aware estimator that delivers **valid inference** for the causal effect of a latent exposure under realistic, SEM-violating conditions (nonlinear effects, non-normal/ordinal measurement), where parametric SEM intervals collapse and naive plug-in intervals undercover — demonstrated against an oracle ceiling and an honest negative control.

The between-imputation variance `B` is the visible mechanism of that contribution; coverage across the SEM-violating cells (vs the negative control where the SEM rightly wins) is the result that satisfies Goal 3.

---

## 13. Appendix — estimand reference

- ψ = E[Y(L+δ)] − E[Y(L)], δ = 1.
- Linear outcome: ψ = 2.0 (any latent shape; the shift equals the coefficient).
- Nonlinear/normal outcome (`1 + 2L + 0.7L² + 0.5LW + 0.5W + sin(2Z)`): ψ ≈ **2.700** = 2δ + 0.7δ² with E[L] = 0; the interaction's marginal contribution is ≈ 0 because E[W] = 0.
- SEM's implied (linear) target under that DGP: ≈ **2.0** — the ≈ 0.7 gap is the curvature it cannot see.
- With mean-zero, variance-matched latents, ψ is invariant to latent shape; latent shape changes *estimation difficulty*, not the estimand.
