# Ordinal Posterior Branch: Problems, Theory, and Possible Solutions

## Purpose

This note explains why the ordinal-indicator branch of the LV-TMLE estimator is currently exploratory only. It is written as a technical planning document, not manuscript prose. It does not add citations. Before using any of these claims in a manuscript, the relevant literature should be checked carefully.

The short version:

- The continuous-indicator branch has a coherent Gaussian posterior approximation for the latent variable.
- The ordinal-indicator branch currently uses factor-score machinery plus a continuous-style posterior standard deviation. That is not a valid posterior for an ordinal measurement model.
- The broader pilot showed severe positive bias or overshoot in ordinal nonlinear cells.
- Ordinal scenarios should stay out of the primary simulation until the posterior is fixed or externally validated.

## Background in Plain Language

The project studies a latent exposure `L`: something not directly observed, measured through indicators `I1`, `I2`, and `I3`.

For continuous indicators, the measurement model is roughly:

```text
I_j = loading_j * L + measurement error
```

If the errors are Gaussian, then the posterior distribution of `L` given the indicators is also approximately Gaussian. That makes it reasonable to summarize each subject's latent uncertainty by:

```text
mean of L | indicators
standard deviation of L | indicators
```

For ordinal indicators, the observed data are categories, such as 1 through 5. The model is different. The usual latent-response story is:

```text
I*_j = loading_j * L + measurement error
I_j = category k if I*_j falls between thresholds k-1 and k
```

We do not observe `I*_j`; we only observe which interval it fell into. This matters because a category is not a noisy continuous number. A response of `4` means "the latent response landed between two thresholds," not "the continuous indicator equals 4."

## The Correct Posterior Target

For LV-TMLE, the object we need is the subject-level posterior distribution:

```text
p(L | I1, I2, I3, Y, W, Z)
```

This can be decomposed as:

```text
p(L | I, Y, W, Z) proportional to
  p(L | W, Z) * p(I1, I2, I3 | L) * p(Y | L, W, Z)
```

For continuous indicators, `p(I | L)` is a product of Gaussian densities.

For ordinal indicators, `p(I | L)` is a product of category probabilities:

```text
P(I_j = k | L) =
  Phi(threshold_k - loading_j * L) -
  Phi(threshold_{k-1} - loading_j * L)
```

where `Phi` is the normal cumulative distribution function on the latent-response scale.

This is the key theoretical point: an ordinal posterior must use probabilities of observed categories, not a Gaussian density for the category labels.

## What the Current Code Does

The current ordinal path in `sem_posterior_full()` fits a lavaan SEM with `ordered = c("I1", "I2", "I3")`, gets regression factor scores with `lavPredict()`, and then approximates:

```text
L | I, W, Z as Normal(mu_noY, post_sd_noY^2)
```

The posterior standard deviation is computed using the same loading/residual-variance style formula used for continuous indicators:

```text
sqrt(1 / (1 / var_l + sum(lambda^2 / theta)))
```

The code itself already notes that this is only an approximation on the latent-response scale.

The outcome-aware step then multiplies that approximate prior/posterior by an outcome likelihood learned from noisy factor scores.

## Why This Is Statistically Fragile

### 1. Ordinal categories are interval information, not continuous measurements

For an ordinal indicator, observing category `k` tells us that an unobserved continuous response fell between two thresholds. The likelihood contribution is a probability mass over an interval.

Treating the result like a continuous Gaussian measurement loses that interval structure.

### 2. The posterior for `L | ordinal indicators` is generally not Gaussian

Even if `L` and the latent responses are Gaussian before thresholding, conditioning on categories creates a truncated, non-Gaussian posterior. With several indicators it may look roughly bell-shaped in some cases, but that is an approximation that must be checked.

The posterior variance can also be subject-specific. A person with all middle-category responses may have different information about `L` than a person with all extreme-category responses.

The current branch uses a single posterior standard deviation for everyone. That is especially risky for ordinal data.

### 3. The latent scale may not match the intervention scale

The target parameter shifts the true latent exposure by `delta = 1`:

```text
psi = E[Y(L + 1)] - E[Y(L)]
```

This only makes sense if the estimated latent variable is on the same scale as the data-generating `L`.

For continuous indicators, fixing the first loading and using the continuous measurement model tends to keep the estimated latent scale closer to the DGP scale.

For ordinal models, the latent-response scale is identified through thresholds, loadings, and residual variance conventions. Factor scores may be on a scale that is not directly comparable to the true `L` scale. If a "one-unit shift" in the estimated ordinal factor is not the same as a one-unit shift in true `L`, then the TMLE is estimating the wrong intervention.

This scale mismatch is a plausible explanation for why RegCal, SEM, and LV-TMLE all overshot badly in some ordinal pilot cells.

### 4. Outcome conditioning can amplify a bad measurement posterior

Including `Y` in the latent posterior is conceptually correct. Imputation of a latent exposure should include the analysis outcome, otherwise imputation noise is independent of `Y` and can attenuate the effect.

But this only works if the measurement part is calibrated. If `p(L | I, W, Z)` is too diffuse, miscentered, or on the wrong scale, the outcome likelihood can dominate the posterior and pull `L` too aggressively toward values that explain `Y`. That can create positive bias or overshoot.

This is not an argument against outcome-aware imputation. It is an argument for using a coherent measurement likelihood before conditioning on the outcome.

### 5. Rubin's rules cannot fix a wrong posterior

Rubin's multiple-imputation variance combines:

```text
within-imputation variance + between-imputation variance
```

That is useful only if the imputations are draws from a reasonable posterior distribution. If the posterior is centered on the wrong scale or has the wrong variance, Rubin's rules faithfully propagate the wrong uncertainty.

In other words, better variance accounting does not repair a biased latent reconstruction model.

## Evidence From the Current Pilot

The broader six-cell pilot used 5 reps each for nonlinear `N = 250` cells.

Continuous-indicator cells were roughly plausible at this small pilot scale:

- Scenario 14: medium error, normal latent, continuous indicators; LV bias about `-0.09`.
- Scenario 17: large error, normal latent, continuous indicators; LV bias about `+0.10`.
- Scenario 32: medium error, skew latent, continuous indicators; LV bias about `+0.15`.

Ordinal-indicator cells were not plausible:

- Scenario 50: medium error, normal latent, ordinal indicators; LV bias about `+0.35`.
- Scenario 53: large error, normal latent, ordinal indicators; LV bias about `+1.81`.
- Scenario 68: medium error, skew latent, ordinal indicators; LV bias about `+1.12`.

The important pattern is not any single 5-rep estimate. The important pattern is that ordinal cells overshot strongly, and the overshoot appeared across different ordinal stress cells. That is enough to block scaling and manuscript claims for ordinal scenarios.

## Possible Solutions

### Option A: Exclude ordinal cells from the primary study for now

This is the current safest decision.

Use continuous indicators for the primary simulation and treat ordinal cells as exploratory. This allows the project to move forward on the core LV-TMLE contribution without making unsupported claims about ordinal measurement.

Pros:

- Low implementation risk.
- Avoids contaminating the primary results with an unvalidated approximation.
- Keeps the simulation study publishable around a defensible continuous-indicator scope.

Cons:

- The study is less broad.
- Psychology often uses ordinal items, so this limits immediate applied scope.

Recommendation: use this as the primary path unless there is time to implement and validate a coherent ordinal posterior.

### Option B: Implement a quadrature posterior for ordinal indicators

For each subject, evaluate a grid of candidate latent values `L_grid`. For each candidate value, compute:

```text
log p(L | W, Z)
+ sum_j log P(I_j = observed category | L)
+ log p(Y | L, W, Z)
```

Then normalize the grid weights to obtain:

```text
E[L | I, Y, W, Z]
SD[L | I, Y, W, Z]
draws from L | I, Y, W, Z
```

This would directly respect the ordinal threshold likelihood.

Pros:

- Statistically principled.
- Fits naturally with the current quadrature approach used for the outcome-aware update.
- Avoids pretending that ordinal category labels are continuous measurements.

Cons:

- Requires extracting thresholds, loadings, and the latent structural model correctly from lavaan.
- Requires careful scale alignment so `delta = 1` has the intended meaning.
- More code and more diagnostics.

Recommendation: this is the best technical fix if ordinal cells are important for the paper.

### Option C: Calibrate the estimated ordinal latent scale to the true DGP scale in simulations

Because this is a simulation, we observe `L_true`. We could estimate a calibration map:

```text
L_true approx a + b * L_score
```

and transform ordinal factor scores/posterior draws onto the true latent scale before applying TMLE.

Pros:

- May quickly diagnose how much of the ordinal overshoot is scale mismatch.
- Useful as a falsifier or diagnostic.

Cons:

- Not a real estimator for applied data, because applied users do not observe `L_true`.
- Should not be used as the proposed method unless reframed explicitly as an oracle diagnostic.

Recommendation: use only as a diagnostic, not as the production estimator.

### Option D: Use a full Bayesian latent-variable model for ordinal indicators

Fit a model that samples the latent `L` and threshold parameters jointly, then pass posterior draws of `L` into the TMLE.

Pros:

- Conceptually coherent.
- Naturally handles non-Gaussian and non-linear posterior shapes.

Cons:

- Adds a large dependency and major compute cost.
- Changes the project architecture.
- Harder to integrate with the current simulation worker.

Recommendation: not the immediate path unless ordinal indicators become central to the paper.

### Option E: Treat ordinal indicators as a future extension

The paper can focus on continuous or approximately continuous indicators and state that ordinal measurement requires a different posterior construction.

Pros:

- Honest and tractable.
- Preserves a clean methodological story.

Cons:

- Reviewers may ask about Likert-type items.

Recommendation: reasonable if the continuous-indicator results are strong enough.

## Recommended Decision

Keep ordinal cells exploratory for now. Do not include them in the primary simulation array or manuscript claims.

The primary study should proceed with continuous indicators while documenting that ordinal indicators require a threshold-likelihood posterior. The ordinal branch should be developed as a follow-up extension or a separate validation task.

This is not a retreat from the LV-TMLE idea. It is a recognition that the continuous and ordinal measurement models imply different likelihoods, and a valid latent posterior must respect that difference.

## Next Immediate Steps for the Project

1. Run a short continuous-only pilot after dropping `SL.gam`.

   Recommended cells:

   ```text
   Scenario 14: nonlinear, normal latent, N=250, medium error, continuous
   Scenario 17: nonlinear, normal latent, N=250, large error, continuous
   Scenario 32: nonlinear, skew latent, N=250, medium error, continuous
   ```

   Use 5-10 reps first. Confirm that learner warnings remain at zero and that the LV estimates remain broadly plausible.

2. Update analysis scripts to make primary vs exploratory results explicit.

   The worker now writes `Primary_Analysis`. The analysis should use that field to support two modes:

   ```text
   primary: continuous indicators only
   exploratory: include ordinal diagnostics
   ```

   This prevents accidental mixing of ordinal exploratory results into the primary summary plots.

3. Prepare a primary continuous-only array plan.

   The primary scenario IDs are 1-36. Before running 1000 reps per cell, run a medium pilot, for example 50-200 reps across a representative subset.

4. Add an ordinal posterior diagnostic script only if ordinal cells remain strategically important.

   The first diagnostic should compare, in ordinal simulation cells:

   ```text
   correlation(mu_noY, L_true)
   calibration slope of L_true ~ mu_noY
   RMSE(mu_noY, L_true)
   posterior z-score calibration
   effect of rescaling factor scores to L_true
   ```

   If simple scale calibration removes much of the overshoot, the main problem is likely latent scale mismatch. If not, the problem is deeper posterior-shape or likelihood misspecification.

5. If ordinal cells must be included, implement Option B.

   Build a threshold-likelihood quadrature posterior for:

   ```text
   p(L | I, Y, W, Z)
   ```

   Validate it against `L_true` before using it in `run_comparison_v2()`.

