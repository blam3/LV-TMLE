# LV-TMLE Project — Handoff Document

**Project:** Latent-Variable Targeted Maximum Likelihood Estimation (LV-TMLE) — propagating latent-variable (SEM) measurement uncertainty through a targeted causal estimator for a continuous latent exposure.
**Estimand:** Modified Treatment Policy (MTP) shift contrast, ψ = E[Y(L+δ)] − E[Y(L)], with δ = 1.
**Status:** Infrastructure rebuilt and validated by smoke test. The outcome-aware estimator (the statistical fix) is **written but not yet empirically validated**. Next concrete action: write/run the outcome-aware smoke test (`05b`).
**Audience:** A collaborator or a future session picking this up cold.

---

## 1. Executive summary

The original implementation demonstrated the *opposite* of its goal: the proposed LV-TMLE lost to both competitors on bias, RMSE, and coverage in essentially every cell. We traced this to specific, fixable causes — not bad luck. The project has since been rebuilt: a fair, SEM-adversarial set of data-generating processes; a single correct cross-fitted shift-TMLE engine shared by all estimators; per-scenario truth computed by Monte Carlo; an oracle ceiling; and a corrected variance via Rubin's rules. A single-cell smoke test confirmed the **infrastructure is sound** (oracle recovers the truth, the between-imputation variance is positive, the SEM shows the expected curvature bias) and confirmed the **last remaining method-level bug**: imputing the latent from the indicators alone re-introduces attenuation. The fix — imputing the latent conditional on the outcome — is implemented in `06_outcome_aware.R` and awaits validation.

The most important strategic update: a simple regression-calibration plug-in already nearly recovers the point estimate (bias −0.09, coverage 0.93 in the test cell). So **the paper's contribution should be reframed around principled, robust *inference* (coverage that holds where simpler plug-ins' standard errors break), not bias/RMSE.**

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
| **Outcome-aware imputation (the point-estimate fix)** | **Written, NOT validated (`06_outcome_aware.R`)** |
| Outcome-aware smoke test (`05b`) | **Not written — immediate next step** |
| `run_comparison_v2` wired into worker/analysis | **Not done** |
| Full 72-scenario array run | Not started |
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
| `02_run_sim.R` (revised) | SLURM worker; passes DGP factors to `generate_data`, records Naive/SEM/Robust/Oracle, computes coverage vs per-scenario `true_psi` | Complete for v1 estimators; **must be repointed to `run_comparison_v2`** |
| `04_analyze_results.R` (revised) | Aggregates vs per-row `True_Psi`, includes Oracle, fixes Covariance→Coverage naming, writes LV variance decomposition | Complete for v1 columns; update for v2 estimator set |
| `05_smoke_test.R` | Single-cell diagnostic harness | Run; results captured below |
| `06_outcome_aware.R` | `sem_posterior_full` (closed-form L\|I,Y,W,Z posterior) + `run_comparison_v2` (SEM, RegCal, outcome-aware LV-TMLE, Oracle, optional MI-no-Y ablation) | **Delivered, not yet validated** |
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

- **`06` is unvalidated.** Predicted to recover ψ with an honest Rubin SE. Two failure modes to watch in `05b`: (a) a small residual bias toward the SEM target from **congeniality** — the imputation model is linear in L while the analysis is nonlinear, so a linear conditioning model cannot fully feed the curvature; (b) **overshoot** (bias turning positive) if the outcome dominates the conditioning and feeds back. With three informative indicators, (b) is unlikely, but verify.
- **Robustness tension.** Conditioning the imputation on Y *via the parametric SEM* reintroduces dependence on the SEM's outcome model — partly in tension with the semiparametric-robustness selling point. The principled endpoint is the efficient influence function that integrates over the latent without leaning on the SEM's Y-model; that derivation is **not done** and is the deeper "technical heart" of the paper.
- **Ordinal cells are approximate.** The Gaussian conjugate posterior in `06` is exact only for continuous indicators; the ordinal branch falls back to a `lavPredict`-based approximation. Either validate it (e.g., against a WLSMV-based posterior) or restrict strong claims to continuous indicators.
- **Compute cost.** LV-TMLE is now M complete-data TMLE fits per replication; the full grid is 72 scenarios × ~1000 reps × M. Budget runtime/memory and bump the SBATCH array to `1-72` before launching.
- **Degenerate-fit reporting.** "Zero crashes" can still hide runs where every learner fell back to the mean. Logging fallback activation rates (critique §2.9) is still outstanding.

---

## 10. Ongoing tasks (prioritized)

**Immediate**
1. **Write and run `05b`** — outcome-aware smoke test on the same cell, calling `run_comparison_v2(..., include_noY = TRUE)`. Confirm: Oracle ≈ 2.70; outcome-aware LV-TMLE recovers ψ (bias near 0, not overshooting); B > 0; RegCal (no-Y) mildly attenuated; MI-no-Y attenuated (ablation). *(This was the agreed next step and is not yet done.)*
2. If `05b` passes, **wire `run_comparison_v2` into `02_run_sim.R`** (new columns: SEM, RegCal, LV, Oracle [+ LV_noY]) and update `04_analyze_results.R` to summarize the new estimator set.

**Before the full run**
3. **Precompute `psi_lookup.rds` once** (`saveRDS(build_psi_lookup(n_mc = 2e7), "psi_lookup.rds")`) to avoid a write race across workers.
4. **Update `03_submit_job.sh`**: `--array=1-72`, larger `--time`/`--mem`.
5. **Partial pilot run** (a handful of scenarios, ~200 reps) to confirm runtime and cross-cell behavior before committing the full array.

**Full study and write-up**
6. Run the full 72 × 1000 array; regenerate plots and the variance decomposition.
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

# 3. NEXT: outcome-aware smoke test (to be written)
#    sources 00, 01, 06; calls run_comparison_v2(..., include_noY = TRUE)
Rscript 05b_smoke_test_oaware.R

# 4. Full array (after wiring v2 and updating the .sh)
sbatch 03_submit_job.sh          # --array=1-72
Rscript 04_analyze_results.R
```

Dependencies: R with `lavaan`, `data.table`, `SuperLearner`, `earth`, `rpart`, `gam`, `dplyr`, `ggplot2`, `tidyr`.

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
