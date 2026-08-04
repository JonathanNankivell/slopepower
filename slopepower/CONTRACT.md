# slopepower R port — interface contract

**Binding on all layers.** Parallel agents implement against this document, not against each
other's code. Do not change a signature or field name here without saying so explicitly in your
final report.

Reference: Nash, Morgan, Frost & Mulick (2021), *Stata Journal* 21(3): 575–601, and the original
`slopepower.ado` v2.1 in the parent directory.

---

## 0. File ownership

| File | Owner | Contents |
|---|---|---|
| `R/utils.R` | (foundation, already written) | shared validators, `%||%`, `qnorm` helpers |
| `R/design.R` | Layer 2 | `trial_design()`, `print.trial_design()` |
| `R/params.R` | Layer 1 | `slope_params()`, `slope_params_manual()`, `print.slope_params()` |
| `R/power.R` | Layer 3 | `slope_sigma()`, `slope_var()`, `slope_effect_size()`, `slope_power()`, `print.slope_power()`, `as.data.frame.slope_power()` |
| `R/grid.R`, `R/bootstrap.R`, `R/compat.R` | Layer 4 | `slope_power_grid()`, `slope_bootstrap()`, `slopepower()` |

**Touch only your own files.** `R/utils.R` is read-only for all layers.

---

## 1. Notation

Mapping from the paper to code, fixed for all layers:

| Paper | Code field | Meaning |
|---|---|---|
| σ²_a | `sigma2_intercept` | between-subject variance of random intercepts |
| σ²_b | `sigma2_slope` | between-subject variance of random slopes |
| σ_ab | `sigma_cov` | covariance of random intercept and slope |
| σ²_ε | `sigma2_residual` | within-subject residual variance |
| β′₁ | `slope` | slope of the untreated / case group |
| β′₁,hc or β′₂-related | `slope_comparator` | slope of healthy controls, or of the treated arm |
| β₂ | `tte` | target treatment effect for the future trial |
| s*² | `var_tte` | treatment-effect variance for a notional two-person trial |

---

## 2. `slope_params` object (Layer 1)

S3 class `"slope_params"`, a list with **exactly** these fields:

```r
list(
  slope             = <dbl>,   # untreated / case slope, per unit of `time`
  slope_comparator  = <dbl>,   # healthy-control or treated-arm slope; NA_real_ if none
  comparator        = <chr>,   # one of "none", "healthy", "treated"
  sigma2_intercept  = <dbl>,   # > 0
  sigma2_slope      = <dbl>,   # > 0
  sigma_cov         = <dbl>,   # unconstrained sign
  sigma2_residual   = <dbl>,   # > 0
  n_obs             = <int>,   # observations used in the fit (post NA-removal); NA for manual
  n_subjects        = <int>,   # subjects used in the fit; NA for manual
  common_variance   = <lgl>,   # TRUE if the comparator RE block was reduced (Stata `nocontvar`)
  time_shifted      = <lgl>,   # TRUE if any subject's first visit was moved to 0
  fit               = <model or NULL>,
  call              = <call>
)
```

Variance components are always those of the **untreated / case** group. When
`comparator == "healthy"` the healthy controls contribute **only** `slope_comparator`; their
variance components are estimated (to avoid contaminating the cases' estimates) and then
discarded. This mirrors the Stata behaviour and paper §2.3.

---

## 3. `trial_design` object (Layer 2)

S3 class `"trial_design"`, a list with **exactly** these fields:

```r
list(
  visits       = <dbl vector>,  # sorted, strictly increasing, visits[1] == 0, length >= 2
  dropout      = <dbl vector>,  # length == length(visits) - 1, INCREMENTAL, each >= 0
  has_dropout  = <lgl>,         # TRUE iff any(dropout > 0)
  dropout_type = <chr>          # "incremental" or "cumulative", as supplied by the user
)
```

`dropout[j]` is the proportion of randomised participants whose **last attended visit is
`visits[j]`** — i.e. they attend `visits[1:j]` and miss `visits[(j+1):K]`. Therefore:

- `dropout[1]` = attend baseline only (contributes **no** slope information)
- `sum(dropout) <= 1`, with tolerance `1e-8`
- completers proportion = `1 - sum(dropout)`

`dropout_type = "cumulative"` input is converted to incremental on construction via
`diff(c(0, cumulative))`; the cumulative vector must be non-decreasing and bounded by 1.
`dropout = NULL` means no dropout and yields a zero vector.

---

## 4. `slope_power` result object (Layer 3)

S3 class `"slope_power"`, a list with **exactly** these fields:

```r
list(
  n, n_per_arm, n_requested,     # n_requested = user's n before the even-number adjustment; NA when solving for n
  power, alpha,
  effectiveness,                 # NA_real_ when target == "observed"
  target,                        # "effectiveness" or "observed"
  tte,                           # target treatment effect, beta_2
  var_tte,                       # s*^2, or the effective (dropout-weighted) equivalent
  effect_size,                   # dropout-weighted standardised effect, signed
  slope_difference,              # slope - reference_slope
  reference_slope,
  solve_for,                     # "n" or "power"
  params, design                 # the input objects, retained
)
```

---

## 5. The mathematics — implement exactly as written

### 5.1 Covariance matrix at arbitrary visit times

For visit times `t = visits` (length `m`, including baseline 0), the marginal covariance is

```
Sigma[i, j] = sigma2_intercept
            + t[i] * t[j] * sigma2_slope
            + (t[i] + t[j]) * sigma_cov
            + (i == j) * sigma2_residual
```

Vectorised:

```r
Sigma <- sigma2_intercept +
         outer(t, t) * sigma2_slope +
         outer(t, t, "+") * sigma_cov +
         diag(sigma2_residual, length(t))
```

This is Σ* from paper p.579, generalised to arbitrary real `t`. **It is a deliberate divergence
from the Stata code**, which builds Σ on a unit-integer grid and selects rows with a matrix
product; the two agree exactly at integer times. There is no `scale()` argument in this port —
the caller expresses `visits` in the same units as the fitted `time` variable.

### 5.2 Two-person treatment-effect variance

```
Sigma_star = blockdiag(Sigma, Sigma)                       # 2m x 2m, one person per arm
X_star     = rbind(cbind(1, t, 0), cbind(1, t, t))         # 2m x 3, columns (1, t, g*t)
F_star     = solve(t(X_star) %*% solve(Sigma_star) %*% X_star)
s_squared  = F_star[3, 3]
```

`slope_var(params, visits)` returns `s_squared`. This is eq. (5) of the paper; element [3,3] is
V(β̂₂) because column 3 of `X_star` is the group-by-time interaction.

### 5.3 Reference slope and target treatment effect

```
comparator == "none"                             -> reference_slope = 0
comparator == "healthy"                          -> reference_slope = slope_comparator
comparator == "treated" && target == "effectiveness" -> reference_slope = 0
comparator == "treated" && target == "observed"  -> reference_slope = slope_comparator,
                                                    effectiveness forced to 1
```

```
slope_difference = slope - reference_slope
tte              = -effectiveness * slope_difference
```

Note the third row: with `target = "effectiveness"` on RCT data the treated-arm slope is
**ignored** and the effect is measured toward zero. That is Stata's model-3 default and it is
intentional. `target = "observed"` is Stata's `usetrt`.

### 5.4 Effect size with dropout (Dawson–Lagakos pattern mixture)

Let `K = length(visits) - 1`. Strata are indexed `j = 1..K`, where stratum `j` attends
`visits[1:j]`, plus a completer stratum attending all visits.

```
es_full = slope_difference / sqrt(slope_var(params, visits))
eff2    = (1 - sum(dropout)) * es_full^2
for (j in 2:K) {                              # j = 1 is skipped: baseline only carries no slope
  es_j <- slope_difference / sqrt(slope_var(params, visits[1:j]))
  eff2 <- eff2 + dropout[j] * es_j^2
}
effect_size = sign(slope_difference) * sqrt(eff2)
```

Skipping `j = 1` is correct: a participant with only a baseline measurement yields a singular
design matrix and contributes zero information, i.e. an infinite stratum-specific sample size and
a zero contribution to the sum. Layer 2 warns when `dropout[1] > 0`.

### 5.5 Sample size and power

Both directions, matching Stata exactly:

```r
z_a <- qnorm(1 - alpha / 2)

# solving for n:
n_per_arm <- ceiling((z_a + qnorm(power))^2 / (effect_size * effectiveness)^2)
n         <- 2 * n_per_arm

# solving for power:
n_actual  <- 2 * floor(n_requested / 2)        # force even, split 1:1
n_per_arm <- n_actual / 2
power     <- pnorm(abs(effect_size) * effectiveness * sqrt(n_per_arm) - z_a)
```

The `effectiveness` factor appears here **in addition to** its appearance inside `tte`. This is
**not** double counting: `effect_size` is on the `slope_difference` scale while `tte` is on the
β₂ scale, and dividing by `effectiveness^2` converts between them, reducing exactly to eq. (6),
`N = {(z_{1-α/2} + z_{1-β}) s* / β₂}²`. Do not "simplify" this away. Under
`target = "observed"`, `effectiveness == 1` and both formulas reduce correctly.

### 5.6 Reported `var_tte`

```
no dropout : var_tte = slope_var(params, visits)
dropout    : var_tte = n_per_arm * tte^2 / (z_a + qnorm(power))^2
```

The dropout branch back-solves the effective s*² because no single value applies across strata.
It is self-consistent in both the sample-size and the power branch (the algebra cancels).

---

## 6. Required guards

Errors (not warnings, not silent `NA`):

- `slope_difference == 0`, or `effect_size == 0` — sample size is undefined. Stata returns `N = .`;
  we error with a clear message.
- `alpha` outside (0, 1); `power` outside (0, 1); `effectiveness` outside (0, 1]
- `n < 2`, or non-integer `n`
- both or neither of `n` and `power` supplied (default: `power = 0.8` when both are `NULL`)
- `length(dropout) != length(visits) - 1`
- `sum(dropout) > 1 + 1e-8`
- non-positive variance components in `slope_params_manual()`
- a covariance matrix that is not positive definite

Warnings:

- `tte` points away from benefit (i.e. treatment would need to make the slope more extreme)
- `dropout[1] > 0` (those participants contribute nothing)
- any subject's time origin had to be shifted (Layer 1)
- `abs(slope) / se(slope) < 2.5` in `slope_bootstrap()` (paper §2.6)

**Floating point:** compare dropout sums with a tolerance. `dropout = c(0.3, 0.3, 0.4)` sums to
exactly 1 in decimal and must be accepted; naive accumulation makes it `-6.7e-17`.

---

## 7. Published results the test suite must reproduce

All from Nash et al. (2021). Data files are `../slpower1.dta`, `../slpower2.dta`, `../slpower3.dta`
relative to the package root, read with `haven::read_dta()`.

`slpower1` — 200 subjects, `visit` 0..3, single group, `effectiveness = 0.33`:

| Paper | Design | Expected |
|---|---|---|
| p.588 | `visits = c(0,1,2)` | slope −1.672, tte 0.552, **N = 712**, per arm 356 |
| p.588 | `visits = c(0,1,2,5)`, dropout `c(0,0,0.1)` | **N = 328**, per arm 164 |
| p.589 | `visits = c(0,0.5,1,1.5,2)` | **N = 620**, per arm 310 |

The third replaces Stata's `scale(0.5) schedule(1 2 3 4)`. Under this port's real-valued time it
must be expressible **either** as half-year visits on the original year scale (above) **or** by
rescaling time and using `visits = 0:4`; both must give N = 620. Stata reports slope −0.836 and
tte 0.276 on the rescaled axis, which are the year-scale values halved.

`slpower2` — 500 subjects (250 cases, 250 healthy controls), `vdate` a date, `case` 0/1.
Convert time to **years** (`as.numeric(vdate) / 365`), `effectiveness = 0.33`:

| Paper | Design | Expected |
|---|---|---|
| p.590 | `visits = c(0,1,2)` | slope cases −1.715, controls 0.975, difference −2.690, tte 0.888, **N = 296**, per arm 148 |
| p.593 | `visits = c(0,1,2)`, dropout `c(0.05,0.05)`, `n = 200` | **power = 0.597** |

`slpower3` — 150 subjects, `visit` in {0, 0.5, 2}, `treat` 0/1, `target = "observed"`:

| Paper | Design | Expected |
|---|---|---|
| p.594 | `visits = c(0,2,3)`, dropout `c(0.2,0.1)` | slope control −1.852, treated −1.104, difference −0.747, tte 0.747, **N = 318**, per arm 159 |

Table 1 (p.595), `slpower1`, `n = 450`, `effectiveness = 0.33`, powers as percentages:

| visits | dropout | power |
|---|---|---|
| `c(0,3)` | none | 0.798 |
| `0:3` | none | 0.817 |
| `seq(0,3,0.5)` | none | 0.865 |
| `c(0,3)` | `0.15` | 0.732 |
| `0:3` | `rep(0.05,3)` | 0.771 |
| `seq(0,3,0.5)` | `rep(0.025,6)` | 0.828 |
| `c(0,3)` | `0.3` | 0.648 |
| `0:3` | `rep(0.1,3)` | 0.716 |
| `seq(0,3,0.5)` | `rep(0.05,6)` | 0.783 |

**Tolerances.** Slopes and variance components come from a REML fit and will differ from Stata in
the last digits: compare to the paper's 3 d.p. with `tolerance = 5e-4` on the printed value.
Sample sizes are integers after `ceiling()` and should match **exactly**; if one is off by one,
report it rather than loosening the tolerance — an off-by-one means the underlying real-valued N
sits within rounding distance of an integer and that fact is worth knowing. Powers are printed to
3 d.p. in the paper; use `tolerance = 1e-3`.

---

## 8. Style

- Base R plus `nlme` and `stats` only in `R/`. No tidyverse inside the package (it is available
  for tests and exploration).
- `snake_case` functions, S3 classes, no R6/S4.
- Every exported function gets roxygen2 comments (`#'`) even though we are not running roxygen —
  they document intent and can generate `man/` later. Hand-write `NAMESPACE` entries in your
  report; the integrator maintains the file.
- Errors via `stop()` with a leading context, e.g. `stop("slope_power(): ...")`.
- No `library()` or `require()` calls inside `R/`. Use `nlme::lme`, `stats::qnorm`, etc.
