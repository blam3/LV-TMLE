# EIF and Posterior-Update Audit Notes

Date: 2026-06-29

Purpose: record the current status of the efficient influence function (EIF) logic and the diagnostic-only posterior-update variants tested before any production estimator change.

Update: the full observed-data finite-dimensional derivation is now in `OBSERVED_DATA_EIF_DERIVATION.md`. That document supersedes the earlier shortcut of treating `E[D_full | O]` as the complete observed-data EIF when posterior/model parameters are estimated.

## 1. Full-data shifted-contrast EIF

For complete data `O = (Y, A, X)`, where `A` is the latent exposure value treated as observed and `X = (W, Z)`, the target is

`psi = E[Q(A + delta, X)] - E[Y]`.

For the shifted mean `psi1 = E[Q(A + delta, X)]`, the standard modified treatment policy EIF is

`D1(O) = H(A, X) * {Y - Q(A, X)} + Q(A + delta, X) - psi1`,

with

`H(A, X) = g(A - delta | X) / g(A | X)`.

For the contrast, subtract the EIF for `E[Y]`:

`Dpsi(O) = D1(O) - {Y - E[Y]}`.

This is the structure implemented in `shift_tmle()` in `01_setup.R`, including the previously missing `H * (Y - Q)` term. For the complete-data/oracle estimator, this part is conceptually correct.

## 2. What Is Not Yet Fully Derived

The current LV-TMLE procedure uses posterior draws of the latent `L` and combines complete-data TMLE fits with Rubin's rules. That is a practical multiple-imputation estimator, but it is not yet a closed-form observed-data EIF for

`O_obs = (Y, I1, I2, I3, W, Z)`.

If the manuscript claims an observed-data EIF, the derivation must account for all nuisance parameters that define the posterior update, including at least:

- factor-score mean parameters and loadings,
- indicator residual variances,
- latent structural mean and variance,
- outcome-link parameters,
- outcome residual variance / likelihood weighting,
- any residual-variance correction or likelihood tempering used in the posterior.

In other words, if `eta` denotes the fitted posterior/link parameters, the finite-dimensional observed-data EIF cannot stop at a posterior-averaged full-data EIF. It requires the projected observed scores and observed information:

`D_eff(O) = xi_eta(X) - psi + dot_psi_eta^T I_eta^{-1} S_eta^O(O)`,

where `S_eta^O(O) = E[S_eta^F | O]`. Cross-fitting helps evaluate an EIF-based estimator, but it does not remove first-order nuisance effects from a plug-in estimator when the scientific target depends on estimated latent-model parameters. The current code does not yet implement this correction.

## 3. Mean and Variance Parameters

The user concern about "all parameters, e.g. the mean and the variance" is valid. The posterior mean and variance are functions of estimated SEM and outcome-link quantities. The variance pieces are not passive constants:

- `post_sd_noY` depends on latent variance, loadings, and indicator residual variances.
- `var_y_eff` depends on outcome residual variance and the latent-residual correction.
- Rubin's `B` captures between-draw variation, but it is not a substitute for deriving the influence contribution of the estimated posterior variance parameters.

Therefore, current simulation inference can be evaluated empirically, but a manuscript-level EIF claim should be framed cautiously until the nuisance-parameter correction is derived or a defensible cross-fitted asymptotic argument is written.

## 4. Diagnostic-Only Posterior Variants Tested

`05d_posterior_gap_diagnostics.R` now tests conservative variants without changing `06_outcome_aware.R`:

- `me_raw_var`: measurement-error-corrected quadratic link with raw residual variance, no latent-residual variance subtraction.
- `current_y_power_0.5`: current link with half-strength outcome likelihood.
- `current_y_power_0.25`: current link with quarter-strength outcome likelihood.
- `me_raw_y_power_0.5`: raw-variance ME link with half-strength outcome likelihood.
- `cf_me_corr`: deterministic-fold cross-fitted ME link with corrected variance.
- `cf_me_raw`: deterministic-fold cross-fitted ME link with raw variance.

Run:

```r
Rscript --vanilla -e 'options(lv_tmle.diag_reps=15, lv_tmle.diag_scenarios=c(14,17,32), lv_tmle.diag_output="posterior_gap_diagnostics_variants_scen_14_17_32.rds"); source("05d_posterior_gap_diagnostics.R")'
```

Results from the 15-rep diagnostic:

| Scenario | Current RMSE | Best conservative RMSE | Oracle-link RMSE | Reading |
|---:|---:|---:|---:|---|
| 14 | 0.320 | 0.316 (`me_raw_var`) | 0.288 | Small improvement; cell already stable. |
| 17 | 0.517 | 0.497 (`current_y_power_0.5` / `0.25`) | 0.433 | Conservative outcome weighting helps, but oracle-link gap remains. |
| 32 | 0.313 | 0.305 (`current_y_power_0.5`) | 0.279 | Tempering reduces over-pull under skew latent shape. |

The conservative variants also generally moved posterior calibration closer to `z_sd = 1`. Cross-fitting alone did not solve the gap; raw or tempered likelihood weighting was more useful.

## 5. Current Recommendation

Do not change the production estimator yet. The next diagnostic step should evaluate whether a tempered posterior update improves estimator-level bias/coverage, not only posterior recovery. The lowest-risk candidate to test first is `current_y_power_0.5`, because it improves scenarios 17 and 32 while leaving scenario 14 essentially unchanged.

Any production change should be gated by an isolated pilot across scenarios 14/17/32 and then a broader continuous-only pilot. The full primary array should remain blocked until this is done.
