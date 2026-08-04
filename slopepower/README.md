# slopepower

An R port of the Stata `slopepower` command of Nash, Morgan, Frost & Mulick (2021),
*Stata Journal* 21(3): 575–601. Sample size and power for two-arm parallel trials
whose outcome is a rate of change (slope) estimated by a linear mixed model.

The calculation has two stages: estimate slope and variance components from
previously collected longitudinal data (or supply them directly), then combine
them with the visit schedule of a proposed trial.

## Installing

The package source lives in this directory. From any other directory:

```r
install.packages("remotes")   # once
remotes::install_local("~/Documents/Masters Documents/Summer project/transpilation/slopepower")
```

or from a shell:

```sh
R CMD INSTALL "$HOME/Documents/Masters Documents/Summer project/transpilation/slopepower"
```

Then in any project:

```r
library(slopepower)
```

Requirements: R >= 4.1 and `nlme` (which ships with R). `haven` is only needed if
you read Stata `.dta` files; `lme4` and `testthat` only for the test suite.

Re-run the install command after pulling changes — there is no auto-reload. If you
are editing the package itself, `devtools::load_all("<that path>")` is faster.

## Basic use

```r
library(slopepower)

# Stage 1 — fit the mixed model to existing longitudinal data.
#   formula is  outcome ~ time | subject
pars <- slope_params(sdmt ~ visit | id, data = mydata)

# Stage 2 — describe the proposed trial.
design <- trial_design(visits = c(0, 1, 2))          # baseline + two follow-ups

# Then ask one of the two questions.
slope_sample_size(pars, design, effectiveness = 0.33)             # how many people?
slope_power(pars, design, n = 200, effectiveness = 0.33)          # what power at N = 200?
```

The two are separate functions on purpose. They take different inputs, answer
different questions, and return differently shaped objects, so there is no
argument you can leave out to switch between them:

| | input | output | class |
|---|---|---|---|
| `slope_sample_size()` | `power` (default 0.8) | `n`, `n_per_arm` | `slope_sample_size` |
| `slope_power()` | `n` (required) | `power` | `slope_power` |

Both print a formatted summary, and both have an `as.data.frame()` method giving
one tidy row with the same columns, so results from either can be bound into one
table — a `solve_for` column records which question produced each row. `alpha`
defaults to 0.05, two-sided.

`design` may be a bare numeric vector of visit times instead of a `trial_design`
object, so `slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)` also
works.

### Units of time

There is no `scale()` argument as there is in Stata. Whatever units the `time`
term in the formula is in, `visits` must use the same ones. Pick units on the
order of the study duration — fitting years-long follow-up on a day axis leaves
the random-slope variance near zero and converges poorly. If visits are recorded
as dates:

```r
pars <- slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, data = mydata)  # years
```

### Comparison groups

Three scenarios, matching §2.3 of the paper:

```r
slope_params(y ~ t | id, data = d)                    # single untreated group; effect measured toward slope 0
slope_params(y ~ t | id, data = d, healthy = case)    # observational: cases (1) vs healthy controls (0)
slope_params(y ~ t | id, data = d, treated = treat)   # previous RCT: treated (1) vs control (0)
```

`healthy` and `treated` take bare column names and are mutually exclusive. With
`treated`, `target = "observed"` targets the treatment effect actually observed
in that trial (Stata's `usetrt`) instead of an `effectiveness` fraction of the
slope difference.

### Dropout

`dropout[j]` is the proportion of participants whose **last attended visit is
`visits[j]`**, one entry per follow-up visit:

```r
trial_design(visits = c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
trial_design(visits = c(0, 1, 2, 5), dropout = c(0.05, 0.1, 0.2),
             dropout_type = "cumulative")   # converted to incremental
```

Participants who attend baseline only carry no slope information; a non-zero
first element warns.

## Parameters without data

If you have published estimates rather than raw data:

```r
pars <- slope_params_manual(
  slope            = -1.672,
  sigma2_intercept = 100, sigma2_slope = 2,
  sigma_cov        = 5,   sigma2_residual = 10,
  slope_comparator = 0.975, comparator = "healthy"   # optional
)
```

## Comparing many designs

The grid functions split the same way:

```r
# power of each design at a fixed N — this is table 1 of the paper
slope_power_grid(
  pars, n = 450, effectiveness = 0.33,
  visits  = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
  dropout = list(none = NULL, `5pc` = dropout_rate(0.05), `10pc` = dropout_rate(0.10))
)

# N each design needs for a fixed power
slope_sample_size_grid(
  pars, power = 0.8, effectiveness = 0.33,
  visits  = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
  dropout = list(none = NULL, `5pc` = dropout_rate(0.05))
)
```

Both return one row per visits × dropout combination, with identical columns. Use
`dropout_rate(rate, per = 1)` rather than a fixed vector when schedules have
different numbers of visits, so the same withdrawal rate is expanded correctly for
each.

## Uncertainty in the stage-one estimates

```r
slope_bootstrap(pars, R = 200, design = c(0, 1, 2), effectiveness = 0.33,
                statistic = "n", type = "bca", seed = 1)
```

Resamples subjects (stratified by group where relevant) and refits the mixed model
each time, so it is slow — start with a small `R`. Needs a `slope_params` object
from `slope_params()`; manually supplied parameters carry no data to resample.

`statistic` picks the entry point: `"power"` goes through `slope_power()` and so
needs an `n` in `...`, while `"n"` and `"tte"` go through `slope_sample_size()`
and take a target `power` instead. Passing the statistic you are bootstrapping as
an input is an error rather than a zero-width interval.

## Porting existing Stata scripts

`slopepower()` mirrors the Stata command's interface argument for argument —
`schedule` lists follow-up visits only, `dropouts` aligns with it, `scale` divides
time, the data type is declared with `obs`/`rct`, and a single call computes
either a sample size or a power depending on whether you pass `power` or `n`:

```r
d <- haven::read_dta("slpower1.dta")
slopepower(d, "sdmt", "id", "visit", schedule = c(1, 2),
           obs = TRUE, nocontrols = TRUE, effectiveness = 0.33)
```

That last point is the one thing the wrapper keeps that new code should not: it
is Stata's bimodal interface, preserved deliberately for parity. Use `slopepower()`
for mechanical translation and parity checks; prefer the pipeline above —
`slope_params()`, `trial_design()`, then `slope_sample_size()` or `slope_power()` —
for new work.

## Notes

- The package has no generated `man/` pages, so `?slope_sample_size` will not work after
  installing. The roxygen comments above each function in `R/` are the reference,
  and `CONTRACT.md` documents the mathematics and the deliberate divergences from
  the Stata original.
- `R CMD INSTALL` is fine; `R CMD check` will complain about the missing
  documentation.
