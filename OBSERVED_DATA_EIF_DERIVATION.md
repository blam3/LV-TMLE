# Observed-Data EIF Derivation for LV-TMLE

Date: 2026-06-29

This note derives the observed-data influence function for the latent-variable shifted-exposure target under the continuous-indicator latent Gaussian working model used by the simulation. It also states how the posterior mean and variance components enter through nuisance scores.

Important distinction: the current production estimator draws latent values from an estimated posterior and then combines complete-data TMLE fits with Rubin's rules. That is not yet an implementation of the observed-data EIF below. This document is the derivation needed before making manuscript-level efficiency claims or rewriting the estimator.

## 1. Observed Data, Latent Model, and Target

Let

`O = (X, I, Y)`,

where `X = (W, Z)`, `I = (I1, I2, I3)`, and `Y` is the outcome. Let `L` denote the latent exposure.

For the continuous-indicator working model:

`L | X ~ N(m_alpha(X), tau2)`,

`Ij | L ~ N(nu_j + lambda_j L, theta_j)`, independently over `j = 1,2,3` conditional on `L`,

`Y | L, X ~ N(q_gamma(L, X), v_y)`.

Let `eta` collect all conditional-model parameters:

`eta = (alpha, tau2, nu, lambda, theta, gamma, v_y)`.

The covariate distribution `P_X` is left nonparametric. The scientific target is

`psi(P_X, eta) = E_X [ xi_eta(X) ]`,

where

`xi_eta(X) = E_eta[q_gamma(L + delta, X) - q_gamma(L, X) | X]`.

This is the latent modified-treatment-policy contrast. The measurement parameters identify the latent distribution from observed indicators, but in the correctly specified scientific target they do not directly change `psi` except through their role in estimating `eta`.

## 2. Observed Likelihood and Posterior Moments

The observed conditional likelihood is

`p_eta(I, Y | X) = integral p_eta(L | X) p_eta(I | L) p_eta(Y | L, X) dL`.

The posterior density is

`pi_eta(l | O) = p_eta(l | X) p_eta(I | l) p_eta(Y | l, X) / p_eta(I, Y | X)`.

Define posterior moments

`M_k(O) = E_eta[L^k | O]`.

For nonlinear `q_gamma`, these moments and related expectations may need numerical quadrature. For the quadratic working outcome link, moments through order four are enough for outcome mean and variance scores.

## 3. General Observed-Data IF

Because `P_X` is nonparametric and observed, the covariate-distribution part of the influence function is

`D_X(O) = xi_eta(X) - psi`.

For the conditional latent/measurement/outcome model, the observed score is

`S_eta^O(O) = partial log p_eta(I, Y | X) / partial eta`.

By Fisher's identity,

`S_eta^O(O) = E_eta[S_eta^F(X, L, I, Y) | O]`,

where `S_eta^F` is the complete-data score.

Let

`I_eta = E[S_eta^O(O) S_eta^O(O)^T]`

be the observed information matrix for the conditional model, and let

`dot_psi_eta = partial psi(P_X, eta) / partial eta`.

For a correctly specified finite-dimensional conditional model estimated by MLE, the efficient observed-data influence function is

`D_eff(O) = D_X(O) + dot_psi_eta^T I_eta^{-1} S_eta^O(O)`.

This is the key result. It is generally not enough to use `E[D_full | O]` when the posterior/model parameters are estimated. The estimated mean and variance components enter through `S_eta^O`, `I_eta^{-1}`, and `dot_psi_eta`.

## 4. Direct Derivatives of the Target

Let

`Delta q_eta(L, X) = q_gamma(L + delta, X) - q_gamma(L, X)`.

Then

`psi = E_X E_eta[Delta q_eta(L, X) | X]`.

For any conditional-model parameter `eta_j`,

`dot_psi_j = partial psi / partial eta_j`.

Equivalently, for latent-distribution parameters,

`dot_psi_j = E[Delta q_eta(L, X) S_j^L(L | X)]`,

where `S_j^L` is the score of `L | X`.

For outcome-regression parameters,

`dot_psi_gamma = E_X E_eta[partial q_gamma(L + delta, X)/partial gamma - partial q_gamma(L, X)/partial gamma | X]`.

For pure measurement parameters `(nu, lambda, theta)`, `dot_psi_j = 0` for the scientific target if the latent scale and structural model are separately parameterized. These parameters still affect the estimator's IF through `I_eta^{-1} S_eta^O`, because they are statistically coupled with latent structural and outcome parameters in the observed likelihood.

For the outcome residual variance `v_y`, `dot_psi_vy = 0` for the scientific target when `v_y` changes only the conditional variance of `Y | L, X`, not the conditional mean `q_gamma`.

## 5. Projected Scores for Latent Mean and Variance

Let

`m = m_alpha(X)`.

If `m_alpha(X) = a^T f(X)`, then the complete-data score for `a` is

`S_a^F = f(X) {L - m} / tau2`.

The observed score is

`S_a^O = f(X) {M_1 - m} / tau2`.

For the latent variance:

`S_tau2^F = -1/(2 tau2) + {L - m}^2 / (2 tau2^2)`.

Thus

`S_tau2^O = -1/(2 tau2) + {M_2 - 2m M_1 + m^2} / (2 tau2^2)`.

These are the terms that account for uncertainty in the latent mean and latent variance.

## 6. Projected Scores for Indicator Parameters

Define

`e_j(L) = Ij - nu_j - lambda_j L`.

Indicator intercept:

`S_nu_j^F = e_j(L) / theta_j`,

`S_nu_j^O = {Ij - nu_j - lambda_j M_1} / theta_j`.

Factor loading:

`S_lambda_j^F = L e_j(L) / theta_j`,

`S_lambda_j^O = {(Ij - nu_j) M_1 - lambda_j M_2} / theta_j`.

Indicator residual variance:

`S_theta_j^F = -1/(2 theta_j) + e_j(L)^2 / (2 theta_j^2)`.

Using

`E[e_j(L)^2 | O] = (Ij - nu_j)^2 - 2 lambda_j (Ij - nu_j) M_1 + lambda_j^2 M_2`,

the observed score is

`S_theta_j^O = -1/(2 theta_j) + {(Ij - nu_j)^2 - 2 lambda_j (Ij - nu_j) M_1 + lambda_j^2 M_2} / (2 theta_j^2)`.

These `theta_j` scores are the formal variance-parameter contributions for indicator residual variances.

## 7. Projected Scores for Outcome Mean and Outcome Variance

For a differentiable outcome mean `q_gamma(L, X)`:

`S_gamma^F = {Y - q_gamma(L, X)} dot_q_gamma(L, X) / v_y`,

where `dot_q_gamma = partial q_gamma / partial gamma`.

The observed score is

`S_gamma^O = E[{Y - q_gamma(L, X)} dot_q_gamma(L, X) | O] / v_y`.

For the outcome residual variance:

`S_vy^F = -1/(2 v_y) + {Y - q_gamma(L, X)}^2 / (2 v_y^2)`,

so

`S_vy^O = -1/(2 v_y) + E[{Y - q_gamma(L, X)}^2 | O] / (2 v_y^2)`.

For the quadratic working link

`q_gamma(L, X) = beta^T r(L, X)`,

with

`r(L, X) = (1, L, L^2, L W, W, Z, Z^2)^T`,

the outcome mean score is

`S_beta^O = {Y E[r | O] - E[r r^T | O] beta} / v_y`.

The outcome variance score is

`S_vy^O = -1/(2 v_y) + {Y^2 - 2Y beta^T E[r | O] + beta^T E[r r^T | O] beta} / (2 v_y^2)`.

Because `r r^T` contains powers up to `L^4`, this score needs posterior moments through at least the fourth order.

## 8. `post_sd_noY` as a Delta-Method Nuisance Function

The code's `post_sd_noY` is not a primitive parameter. Under the continuous Gaussian indicator model, the indicator-only posterior is

`L | I, X ~ N(mu_I, sigma_I2)`,

where

`C = 1/tau2 + sum_j lambda_j^2 / theta_j`,

`sigma_I2 = C^{-1}`,

`mu_I = sigma_I2 {m_alpha(X)/tau2 + sum_j lambda_j (Ij - nu_j) / theta_j}`.

Thus

`post_sd_noY = sigma_I = sqrt(sigma_I2)`.

Its IF is obtained by the delta method from the IF of `eta`:

`IF_sigma_I2 = grad_eta sigma_I2^T IF_eta`,

`IF_post_sd_noY = IF_sigma_I2 / (2 sigma_I)`.

Key gradients are

`partial sigma_I2 / partial tau2 = sigma_I2^2 / tau2^2`,

`partial sigma_I2 / partial theta_j = sigma_I2^2 lambda_j^2 / theta_j^2`,

`partial sigma_I2 / partial lambda_j = -2 sigma_I2^2 lambda_j / theta_j`.

The posterior mean also has an IF:

`IF_mu_I = grad_eta mu_I^T IF_eta`.

Since `mu_I = sigma_I2 B`, with

`B = m_alpha(X)/tau2 + sum_j lambda_j (Ij - nu_j) / theta_j`,

`d mu_I = B d sigma_I2 + sigma_I2 dB`.

This explicitly shows how latent mean, latent variance, loadings, intercepts, and residual variances enter the posterior mean and variance.

## 9. `var_y_eff` and Outcome Residual Variance

If `v_y` is the actual Gaussian outcome variance, it enters through `S_vy^O` and through the observed information matrix.

The code's `var_y_eff` is a working likelihood variance:

`var_y_eff = raw_var_y - latent_resid_correction`,

with a lower floor at `Y_SD^2`.

If `var_y_eff` is treated as an estimated nuisance parameter in an algorithmic posterior, define a smooth population version

`v_eff(eta) = R(eta) - C_latent(eta)`.

Away from the floor,

`IF_veff = grad_eta v_eff^T IF_eta`.

If it is estimated by an empirical residual moment, then

`IF_R(O) = r_eff(O) - R + grad_eta R^T IF_eta(O)`,

and similarly for `C_latent`.

At the floor,

`var_y_eff = max(v_eff, Y_SD^2)`.

This is non-differentiable at `v_eff = Y_SD^2`. Away from the kink:

- if `v_eff > Y_SD^2`, use `IF_veff`;
- if `v_eff < Y_SD^2`, the floored parameter is locally constant and has IF zero;
- if `v_eff = Y_SD^2`, the parameter is nonregular.

Therefore a smooth observed-data EIF cannot rely on the floored `var_y_eff` unless the analysis establishes that the true value is away from the floor or replaces the floor with a smooth positive parameterization.

## 10. Nuisance-Score Correction

Let `S_eta^O` stack all projected scores:

`S_eta^O = (S_alpha^O, S_tau2^O, S_nu^O, S_lambda^O, S_theta^O, S_gamma^O, S_vy^O)`.

Let `dot_psi_eta` stack the direct target derivatives from Section 4. The parametric conditional-model IF is

`IF_eta(O) = I_eta^{-1} S_eta^O(O)`.

The conditional-model contribution to the target IF is

`D_eta(O) = dot_psi_eta^T IF_eta(O)`.

Equivalently,

`D_eta(O) = dot_psi_eta^T I_eta^{-1} S_eta^O(O)`.

This formula is where all mean and variance nuisance components enter. Even when a component has zero direct derivative, it can still contribute through the inverse observed information if it is statistically coupled with components that have nonzero derivatives.

The full efficient IF under the finite-dimensional conditional model and nonparametric `P_X` is

`D_eff(O) = xi_eta(X) - psi + dot_psi_eta^T I_eta^{-1} S_eta^O(O)`.

This expression has mean zero and satisfies

`d/dt psi(P_t, eta_t)|0 = E[D_eff(O) S_t(O)]`

for regular submodels of the observed data law.

## 11. Cross-Fitting: What It Does and Does Not Justify

Cross-fitting is useful, but it does not make first-order nuisance effects disappear when the target depends on finite-dimensional nuisance parameters.

Cross-fitting justifies replacing unknown nuisance functions by estimates in the EIF when:

1. nuisance estimates are trained on separate folds from the observations used to evaluate the EIF;
2. the estimator solves an EIF estimating equation, such as `P_n D_eff(eta_hat, psi_hat) = o_p(n^{-1/2})`;
3. empirical-process terms are controlled by sample splitting;
4. second-order remainders are `o_p(n^{-1/2})`;
5. non-smooth nuisance transformations, such as the floor in `var_y_eff`, are avoided or are inactive with probability tending to one.

Cross-fitting alone is not a substitute for `dot_psi_eta^T I_eta^{-1} S_eta^O` in a plug-in estimator. It removes overfitting bias and Donsker restrictions; it does not erase the influence of estimating latent means, latent variances, residual variances, loadings, or outcome variances when those parameters are part of the identifying model.

## 12. Consequences for the Current Code

1. `shift_tmle()` has the correct complete-data EIF for the shifted contrast when `L` is observed or imputed as fixed complete data.
2. The observed-data latent-variable EIF under the finite-dimensional Gaussian working model is `D_eff(O)` above, not Rubin's rules alone.
3. `post_sd_noY` must be treated as a delta-method function of `(tau2, lambda, theta)`.
4. Latent mean and latent variance enter through `S_alpha^O` and `S_tau2^O`.
5. Indicator residual variances enter through `S_theta_j^O`.
6. Outcome residual variance enters through `S_vy^O`.
7. `var_y_eff` is a working nuisance with a nonregular floor; it needs either a smooth replacement or an explicit nonregularity caveat.
8. A future EIF implementation should compute posterior moments, projected scores, observed information, `dot_psi_eta`, and then use `D_eff` for one-step/TMLE inference.
9. Until that implementation exists, simulation coverage from the MI/Rubin estimator should be described empirically, not as guaranteed by the fully derived observed-data EIF.

