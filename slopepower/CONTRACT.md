# slopepower R port — porting notes and invariants

What the port guarantees, and — more usefully — which of its oddities are
deliberate. Several things here look like defects and are not: the notes saying
so are the point of the document, because the evidence behind them took a Stata
licence to obtain and cannot be re-derived by reading `slopepower.ado`.

Not shipped to users (see `.Rbuildignore`); this is for whoever maintains the
port. The user-facing documentation is the roxygen help pages. Code and tests
cite these sections by number, so **sections may be revised but not renumbered**.

Reference: Nash, Morgan, Frost & Mulick (2021), *Stata Journal* 21(3): 575–601, and the original
`slopepower.ado` v2.1 in the parent directory.

---

## 1. Notation

Mapping from the paper to code, fixed throughout:

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

## 2. `slope_params` object (`R/params.R`)

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

## 3. `trial_design` object (`R/design.R`)

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

A `dropout_rate(rate, per)` object may be supplied in place of the vector; `trial_design()`
expands it to `(rate / per) * diff(visits)` before the checks above. The expansion is the
constructor's, not the grid's, so `trial_design()` and `slope_*_grid()` cannot expand the same
rate differently. It yields incremental proportions by construction and is therefore rejected
alongside `dropout_type = "cumulative"`. The **stored field is always the expanded numeric
vector** — a `dropout_rate` is never a legal value of `dropout` in the object above.

---

## 4. Stage-two result objects (`R/power.R`)

Sample size and power are **two exported functions, not one function with a mode switch**. They
take different inputs, answer different questions, and return differently shaped objects:

```r
slope_sample_size(params, design, effectiveness = 0.25, target, power = 0.8, alpha = 0.05)
slope_power(params, design, n, effectiveness = 0.25, target, alpha = 0.05)
```

`n` is `slope_power()`'s third argument and is **required** — there is no default sample size to
fall back on. `power` is an ordinary argument of `slope_sample_size()` with its default in the
signature. Neither accepts the other's input: `slope_sample_size(..., n = 450)` and
`slope_power(..., power = 0.8)` are "unused argument" errors from R's own argument matching, so
the mismatch that used to be silently resolved by which argument was left `NULL` can no longer be
expressed. **Do not reintroduce a combined entry point**, and do not give `n` a default.

Both classes inherit from `"slope_result"`, which is what `as.data.frame()` dispatches on. So does
`slope_sample_size_floor` (section 4.3), which answers neither question but bounds the first.

### 4.1 `slope_sample_size`

S3 class `c("slope_sample_size", "slope_result")`, a list with **exactly** these 13 fields:

```r
list(
  n, n_per_arm,
  power, alpha,                  # power is the input here
  effectiveness,                 # NA_real_ when target == "observed"
  target,                        # "effectiveness" or "observed"
  tte,                           # target treatment effect, beta_2
  var_tte,                       # s*^2, or the effective (dropout-weighted) equivalent
  effect_size,                   # dropout-weighted standardised effect, signed
  slope_difference,              # slope - reference_slope
  reference_slope,
  params, design                 # the input objects, retained
)
```

### 4.2 `slope_power`

S3 class `c("slope_power", "slope_result")`. The same 13 fields, plus `n_requested` after
`n_per_arm` — the user's `n` before the even-number adjustment — for 14 in total. `n` is the input
here and `power` the output.

There is **no `solve_for` field**: the class carries that information. `as.data.frame()` still
emits a `solve_for` column, derived from the class, so that rows from both functions bind into one
interpretable table. Column names are identical across both classes and every comparator.

### 4.3 `slope_sample_size_floor`

S3 class `c("slope_sample_size_floor", "slope_result")`, from `R/floor.R`. The same fields as
`slope_sample_size` **minus `design`**, for 12:

```r
list(n, n_per_arm, power, alpha, effectiveness, target, tte, var_tte, effect_size,
     slope_difference, reference_slope, params)
```

There is no `design` because there is no design: the value bounds every schedule at once (section
5.7). **Do not add one**, and do not let it acquire a `visits` or `dropout` argument — the same
rule, for the same reason, as `slope_effect_size()`'s refusal of `effectiveness`. An argument the
answer does not depend on can only be ignored.

It is a **generic**, unlike the two stage-two entry points, with methods on `slope_params` (which
takes `effectiveness`, `target`, `power` and `alpha`, exactly as `slope_sample_size()` does) and on
`slope_result` (which reads all four off the object, so the floor and the design it bounds are
guaranteed comparable). `...` exists only because the generic has it and is rejected by
`reject_dots()`; a silently ignored `design =` would make the result look design-specific.

The `slope_result` method refuses a `slope_power` object whose `power` has saturated at exactly 1,
where `qnorm(1)` is `Inf` — the same double-precision edge that section 5.6 handles for `var_tte`.

`as.data.frame()` gains two consequences: `solve_for` is `"n_floor"`, and `n_follow_up` is
`NA_integer_` rather than a count, because the row belongs to no schedule. The 4.1/4.2 guarantee
that column *names* are identical across classes is unaffected.

### 4.4 `slope_sample_size_grid_boot` (`R/grid_boot.R`)

Not a `slope_result` subclass — a grid was always a data frame rather than a `slope_result`, and
this is a grid. S3 class `c("slope_sample_size_grid_boot", "data.frame")`: the fifteen columns
`grid_evaluate()` computes for a grid cell — both `n`/`n_per_arm` and both
`visits_total`/`visits_per_arm`, unreduced — plus eleven more: `n_mean`, `n_sd`, `n_lower`,
`n_upper`, `tte_mean`, `tte_sd`, `tte_lower`, `tte_upper`, `tte_ci_type`, `ci_type`, `n_failed`.
(`slope_sample_size_grid()` itself returns thirteen of those fifteen — see "Display basis" below —
so this object's own data frame carries two more columns than that function's return value, not
the same fifteen.) `ci_type` and `tte_ci_type` are separate because the two intervals have separate
bias corrections and separate jackknife columns: either can fall back from BCa to percentile
without the other, and the printed table marks each column with the method it actually used.

`tte` depends on `effectiveness` alone — not the visit schedule, the dropout pattern, `power` or
`alpha` — so it is resampled once per distinct `effectiveness` level rather than once per cell,
and cells sharing a level share one interval exactly rather than merely agreeing to rounding.
The printed target-effect frame is keyed by `effectiveness` for the same reason: one row per
level, not one per cell.

What every cell shares — `R`, `type`, `level`, `se`, `n_refit_failed`, `straddle`, and the same six
`slope_*` fields `slope_bootstrap()` itself returns — is carried as **attributes**, not columns: a
value that never varies across rows does not belong in one. `[.slope_sample_size_grid_boot()`
strips them (and the class) on every subset, since a value describing the whole table becomes
false, not merely stale, the moment the table it describes changes shape.

Every cell shares one set of resampled replicates: the resampling scheme (`R/bootstrap.R`'s
`boot_setup()`, `boot_replicate_matrix()`, `jackknife_values()`) depends only on `params`, never on
the design being priced, so `R` replicates refit once price every cell rather than `R` refits *per
cell*. A cell with fewer than two surviving replicates reports `NA` in its own row rather than
aborting the whole grid the way `slope_bootstrap()` aborts a single result (section 6) — a grid
that took several minutes to resample must not be discarded over one bad cell.

### 4.5 Display basis: `per_arm`

Every function in this section that reports a count accepts `per_arm` (default `TRUE`), which
picks whether a *printed* `n` — and, on a grid, `n`'s companion `visits` — reads as participants
per arm or as the trial total. It is **display only** and is carried as an attribute
(`attr(x, "per_arm")`), never as a field of a `slope_result` list and never as a column of a grid
`data.frame` — a value that changes how a number is shown, not what was computed, belongs beside
the object rather than inside it, the same reasoning section 4.4 gives for the bootstrap grid's
shared attributes. `print()` methods take their own `per_arm` argument, defaulting to `NULL`
("whatever the object's attribute says, or per arm if it carries none"), so the printed basis can
be overridden after the fact without recomputing anything: `print(x, per_arm = FALSE)`.

The two plain grids, `slope_sample_size_grid()` and `slope_power_grid()`, are the one exception.
They are undecorated data frames with no print method to consult an attribute at print time, so
for them `per_arm` instead shapes the *returned* data: the `n`/`n_per_arm` and
`visits_total`/`visits_per_arm` pairs each collapse to a single `n`/`visits` column on the chosen
basis (`basis_columns()`, `R/grid.R`), and the object carries thirteen columns, not fifteen. The
attribute is still recorded, so a table already in hand can be identified, but there is no
surviving twin column to fall back on — unlike every other function in this section, where both
bases are always retrievable from the object itself.

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
a zero contribution to the sum. `trial_design()` warns when `dropout[1] > 0`.

### 5.5 Sample size and power

Both directions, matching Stata exactly:

```r
z_a <- qnorm(1 - alpha / 2)

# slope_sample_size():
n_per_arm <- ceiling((z_a + qnorm(power))^2 / (effect_size * effectiveness)^2)
n         <- 2 * n_per_arm

# slope_power():
n_actual  <- 2 * floor(n_requested / 2)        # force even, split 1:1
n_per_arm <- n_actual / 2
power     <- pnorm(abs(effect_size) * effectiveness * sqrt(n_per_arm) - z_a)
```

Both live in one internal `solve_slope()`, called by the two exported functions. The **interfaces**
are separate; the **algebra** is shared, so the sample size a design needs and the power it
achieves cannot drift apart.

The `effectiveness` factor appears here **in addition to** its appearance inside `tte`. This is
**not** double counting: `effect_size` is on the `slope_difference` scale while `tte` is on the
β₂ scale, and dividing by `effectiveness^2` converts between them, reducing exactly to eq. (6),
`N = {(z_{1-α/2} + z_{1-β}) s* / β₂}²`. Do not "simplify" this away. Under
`target = "observed"`, `effectiveness == 1` and both formulas reduce correctly.

### 5.6 Reported `var_tte`

```
no dropout                     : var_tte = slope_var(params, visits)
dropout, slope_sample_size()   : var_tte = n_per_arm * tte^2 / (z_a + qnorm(power))^2
dropout, slope_power()         : var_tte = tte^2 / (effect_size * effectiveness)^2
```

The dropout branch back-solves the effective s*² because no single value applies across strata.

The two dropout forms are the same algebra but **not** the same number, and the difference is
real rather than a rounding artefact to be papered over. Solving for power, `qnorm` inverts the
`pnorm` above and `n_per_arm` cancels exactly, leaving `tte^2 / scaled_effect^2` — the true
effective s*², independent of `n`. That closed form is used directly. It is not merely a
simplification: past a few thousand per arm the power saturates at exactly 1 in double precision,
`qnorm(1)` is `Inf`, and the round trip silently reported `var_tte = 0`.

Solving for n there is no such cancellation, because `n_per_arm` has been rounded up to a whole
participant. The reported value carries that rounding and is slightly larger than the exact one,
bounded above by `n_per_arm / (n_per_arm - 1)` times it. This matches the Stata original, which
reports the same rounded quantity. Do not "fix" the two branches into agreement — assert the
inequality instead.

### 5.7 The floor over all designs (`R/floor.R`)

No counterpart in Stata or in the paper; derived in the "What s\* is" vignette, section 6.

```
slope_var_floor(params) = 2 * (sigma2_slope - sigma_cov^2 / sigma2_intercept)
```

Twice the Schur complement of the random-effects covariance matrix, i.e. twice `Var(b_i | a_i)`.
Positive whenever that matrix is positive definite, which `check_re_covariance()` already
guarantees on every route into a `slope_params` object — so there is no degenerate branch to
write.

Three facts, each load-bearing for how it is documented:

- It is the **infimum of `slope_var(params, visits)` over all `visits`, and is not attained**.
  Since `Sigma = Sigma_0 + sigma2_residual * I`, every contrast has `c'Sigma c > c'Sigma_0 c`
  strictly, so `t' Sigma^-1 t` stays strictly below its limit at every finite schedule. The
  documentation says "infimum, not minimum", and `test-floor.R` asserts the strict inequality over
  random schedules rather than an approximate equality.
- **Dropout can only raise it**, so the bound holds over designs *and* dropout patterns. That is
  why `slope_sample_size_floor()` needs no `design`, not merely why it has none.
- **Length alone does not reach it.** Two visits a distance `t` apart converge as `t -> Inf` on
  `2 * (sigma2_slope - sigma_cov^2 / (sigma2_intercept + sigma2_residual))`, strictly larger
  whenever `sigma_cov != 0`: only repeated measurement recovers the whole baseline correction. Do
  not shorten this to "as the trial gets longer" — a reader who lengthens a two-visit trial aims
  at the wrong number, and `test-floor.R` pins the two limits apart.

The sample size follows through `size_per_arm()`, the single implementation of equation (6) that
`solve_slope()` also calls (section 5.5). The sharing is the point: `slope_sample_size(...)$n >=
slope_sample_size_floor(...)$n` holds by construction rather than by two implementations happening
to agree. After `ceiling()` the bound is attained — on `slpower1` at `effectiveness = 0.33` both
the floor and a fifty-year schedule of 501 visits give N = 236, against 712 for the paper's p.588
design.

---

## 6. Required guards

Errors (not warnings, not silent `NA`):

- `slope_difference == 0`, or `effect_size == 0` — sample size is undefined. Stata returns `N = .`;
  we error with a clear message.
- `alpha` outside (0, 1); `power` outside (0, 1); `effectiveness` outside (0, 1]
- `n < 2`, or non-integer `n`
- `n` missing from `slope_power()` (there is no default sample size), and likewise from
  `slope_power_grid()`
- `slopepower()` only: both `n` and `power` supplied. The Stata command's single bimodal
  interface is mirrored **in the compatibility wrapper alone**; `power` defaults to 0.8 there when
  neither is given.
- `length(dropout) != length(visits) - 1`
- `sum(dropout) > 1 + 1e-8`
- non-positive variance components in `slope_params_manual()`
- a covariance matrix that is not positive definite

Warnings:

- `tte` points away from benefit (i.e. treatment would need to make the slope more extreme)
- `dropout[1] > 0` (those participants contribute nothing)
- any subject's time origin had to be shifted (`slope_params()`)
- `abs(slope) / se(slope) < 2.5` in `slope_bootstrap()` (paper §2.6)
- fewer than two replicates succeed for one cell of `slope_sample_size_grid_boot()` — that cell's
  interval columns are `NA` (collected into one warning naming every such cell); the call errors
  only if every cell is starved

**Floating point:** sum the dropout vector once and compare the total against 1 with a tolerance;
do not accumulate by repeated subtraction. `1 - 0.3 - 0.3 - 0.4` is `-5.551e-17`, so a bare `< 0`
guard would reject `dropout = c(0.3, 0.3, 0.4)`, which is legal and sums to exactly 1 in decimal.

This restates Stata rather than repairing it. `slopepower.ado:265` *does* accumulate by
subtraction, but a Stata local round-trips through a decimal string (`.7 - 0.3` stores as `.4`),
so its residue is exactly 0 and it accepts the list too — verified, see DIVERGENCES.md "Claims
checked and rejected". R has no such rounding, so the check has to be written this way here to
reach the same answer.

---

## 7. Published results the test suite must reproduce

All from Nash et al. (2021). The test suite reads the packaged datasets `slpower1`, `slpower2`
and `slpower3` through `load_paper_data()`, which adds the `time` column in the paper's units.
Because they ship inside the package, the parity tests run from a built tarball and under
`R CMD check` rather than skipping.

The packaged copies are built from `../slpower1.dta`, `../slpower2.dta`, `../slpower3.dta` by
`data-raw/make-data.R`. That makes them a second copy of the same numbers, so
`test-packaged-data.R` rebuilds the conversion from the `.dta` files and demands an exact match.
It is the only file that still needs `haven` and the `.dta` files, and so the only one that
skips outside a source checkout. If the `.dta` files ever change, re-run `data-raw/make-data.R`
or that test will fail.

Rows with an expected **N** are `slope_sample_size()` calls; the one row with an expected **power**
is a `slope_power()` call.

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
the last digits, so the rule is "still prints as the paper's figure": half a unit in the last
printed place. The paper prints slopes and powers to 3 d.p., giving `5e-4` and `5e-4`
respectively — but see the rounding-boundary exception below before applying that literally.
Sample sizes are integers after `ceiling()` and should match **exactly**; if one is off by one,
report it rather than loosening the tolerance — an off-by-one means the underlying real-valued N
sits within rounding distance of an integer and that fact is worth knowing.

**The slpower1 rounding boundary.** `slpower1`'s slope fits to −1.6724999999999988, which clears a
`5e-4` bound against the paper's −1.672 by 1.2e-15 — about one ulp. Two consequences, and both are
load-bearing:

- The value sits just *below* the tie, so nothing actually disagrees about it: Stata's `%5.3f` and
  R's `formatC`, `sprintf` and `round` all print −1.672. `test-paper-parity.R` nevertheless pins the
  value at 4 d.p. (`round(slope, 4) == -1.6725`) rather than comparing at 3, because agreement that
  rests on one ulp is not agreement worth asserting.
- Any check that compares this slope at 3 d.p. is one ulp from failing, and would fail for reasons
  having nothing to do with the port — an `nlme` or BLAS change is enough. `data-raw/make-data.R`'s
  regeneration guard therefore carries a **per-dataset** precision (`slpower1` at 4 d.p., the other
  two at 3) rather than one loosened bound for all three.

Do not "simplify" that guard back to a single tolerance. Loosening the shared bound to clear
`slpower1` also loosens `slpower2` and `slpower3`, which sit nowhere near their boundaries, and
costs the guard its meaning: at `1e-3`, `slpower2` could drift far enough to print −1.716 and still
pass, while the script's own success line printed a number contradicting the paper it claims to
have reproduced.

### 7.1 The Stata reference suite

Beyond the paper's 15 numbers, `../stata-reference/` holds 584 rows generated by running
`slopepower.ado` itself over the option space the paper never visits, plus 37 rows carrying the
fitted variance components. Copies live in `tests/testthat/fixtures/` so the parity tests need no
Stata licence. See that directory's README for how to regenerate.

**Two facts about the original, established under Stata on 2026-08-04, that this port must not
"fix":**

- The dropout-length guard at `slopepower.ado:279` is **live**. `` `sched_length ' `` has a space
  inside the macro name but Stata trims it, so both counters are correct and a wrong-length
  `dropouts()` list returns `_rc = 198`. This port's length validation restates Stata.
- `slopepower.ado:371` is **correct**. `mixed` stores the two residual parameters under
  `residuals(independent, by())` as `lnsig_e:_cons` — the base log-SD of the *reference* group,
  i.e. the controls — and `r_lns2ose:_cons`, a log-**ratio** offset. So `exp(b13+b14)^2` is the
  cases' variance, 10.354254, against 10.698895 for the controls.

**Tolerances against Stata.** Three strengths, and the difference between them is the whole point:

| Comparison | Tolerance | What it tests |
|---|---|---|
| Stata's variance components in via `slope_params_manual()` | 1e-12 | Σ\*, X\*, F\*, dropout weighting, `effectiveness` rescaling — **all 542 rows match exactly** (348 + 88 + 75 + 31, i.e. every comparable row of all four design grids) |
| Fitted components vs Stata's | 1e-3 to 5e-2 relative | how closely `nlme::lme` and `mixed` converge |
| End to end | dominated by the above | nothing extra; it cannot be tighter than the fits agree |

The middle row is calibrated to the port **as it ships** — through `slope_lme_control()`, the
tightened control `slope_params()` passes to every `lme()` call — not to `nlme`'s stock defaults.
Tightening *beyond* `slope_lme_control()` closes most of the remaining gap — `var_tte` on slpower1
at visits `0:3` moves from 8.36388645 to 8.36394022 against Stata's 8.36394172 — but that changes
the package's behaviour for every user to make a test look better, so it is deliberately not done.
(This section used to say the tolerances were calibrated to `nlme`'s defaults and cited 8.36388645
as the stock figure; both were wrong — 8.36388645 is what `slope_lme_control()` gives, and `nlme`'s
untouched defaults give 8.36399564 on that fit, about as far from Stata on the other side. Nothing
here is calibrated to stock `nlme`. Corrected 2026-08-10.) If `slope_lme_control()` is ever
tightened further, retighten these tolerances with it.

Each quantity is checked with a loose bound on the worst row and a tight bound on the median: the
maxima are dominated by fits that barely identify a random slope (two timepoints per subject, or
80 subjects at three visits), and a bound loose enough to admit those would hide a uniform drift.

---

## 8. Style

- Base R plus `nlme` and `stats` only in `R/`. No tidyverse inside the package (it is available
  for tests and exploration).
- `snake_case` functions, S3 classes, no R6/S4.
- Every exported function gets roxygen2 comments (`#'`). These are the source of truth:
  `NAMESPACE` and `man/*.Rd` are generated from them by `roxygen2::roxygenise()` and must not
  be edited by hand. Mark exports with `@export` and internal helpers with `@noRd`; package-level
  `@importFrom` tags live in `R/slopepower-package.R`, one line per tag.
- Errors via `stop()` with a leading context, e.g. `stop("slope_power(): ...")`.
- No `library()` or `require()` calls inside `R/`. Use `nlme::lme`, `stats::qnorm`, etc.
