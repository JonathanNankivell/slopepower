# slopepower

An R port of the Stata `slopepower` command of Nash, Morgan, Frost & Mulick (2021),
*Stata Journal* 21(3): 575–601. Sample size and power for two-arm parallel trials
whose outcome is a rate of change (slope) estimated by a linear mixed model.

The calculation has two stages: estimate slope and variance components from
previously collected longitudinal data (or supply them directly), then combine
them with the visit schedule of a proposed trial.

## The model

One continuous outcome, measured repeatedly on **independent participants**,
whose mean is a straight line in time:

```
y[ij] = mean(t[ij]) + a[i] + b[i] * t[ij] + e[ij]
```

Each participant gets a random intercept `a[i]` and a random slope `b[i]`
(unstructured 2×2 covariance), with independent residuals `e[ij]`. Treatment
appears only as a difference in slopes, and the two arms share a single
intercept — baseline is modelled as part of the outcome vector, not entered as
a covariate. Nothing else is in the mean: no covariates, and no second level of
clustering above the participant.

The design matrices are fixed by the method, and they are not what you would
write by hand in `lme4`. For a previous trial the arms share one intercept and
one variance structure and differ only in slope; for cases against healthy
controls each group gets its own random-effects covariance *and* its own
residual variance, which `lme4` cannot fit at all. See `?slope_params` for the
three models written out.

### What this covers, and what it doesn't

Supported:

- two-arm parallel trials, 1:1 allocation, effect measured as a difference in slopes;
- stage-one data from a single untreated cohort, from cases plus healthy controls, or from a previous two-arm RCT;
- arbitrary, unequally spaced visit times, on a schedule common to all participants;
- monotone dropout, via the Dawson–Lagakos pattern mixture;
- parameters taken from the literature instead of fitted (`slope_params_manual()`).

Not supported:

- **covariate adjustment of any kind** — age, sex, disease duration, centre. The stage-one formula takes a single time term and nothing else, and rejects anything more;
- **multi-level / nested clustering** — visits within participants within sites, clinics or families. There is exactly one random-effects level, the participant;
- **baseline as an ANCOVA covariate** — baseline is a correlated outcome with a common intercept across arms (paper §2.1). If your planned analysis conditions on baseline instead, this is the wrong model;
- more than two arms, unequal allocation, cluster-randomised, crossover or stepped-wedge designs;
- non-linear trajectories — quadratic time, splines, change points — or any estimand that is not a slope difference;
- non-Gaussian outcomes: binary, ordinal, count or time-to-event endpoints;
- residual structures beyond independent errors, e.g. AR(1) or other serial correlation;
- intermittent or non-monotone missingness, and participant-specific visit schedules in the planned trial: dropout is assumed to truncate a common schedule.

## Installing

The package lives on GitHub at
[JonathanNankivell/slopepower](https://github.com/JonathanNankivell/slopepower),
in the `slopepower` subdirectory of the repository (the repository root also
holds the original Stata `.ado`/`.sthlp` files and the reference `.dta` data
this package was ported from). Install straight from there:

```r
install.packages("remotes")   # once
remotes::install_github("JonathanNankivell/slopepower", subdir = "slopepower",
                         build_vignettes = TRUE)
```

If you already have a local clone, install from that copy instead of fetching
from GitHub again:

```r
remotes::install_local("<path to your clone>/slopepower", build_vignettes = TRUE)
```

or from a shell, which knits the vignettes by default (this is what CRAN does):

```sh
R CMD build "<path to your clone>/slopepower" && R CMD INSTALL slopepower_*.tar.gz
```

Building the vignettes needs **pandoc** on `PATH` — see [Vignettes](#vignettes)
below — and adds a little time to the install. If you don't want them, drop
`build_vignettes` (or set it to `FALSE`), or pass `--no-build-vignettes` to
`R CMD build`; either way `vignette("...")` will then report that none exist.
Once installed the built HTML is self-contained; reading it needs neither
pandoc nor a network connection.

Then in any project:

```r
library(slopepower)
```

Requirements: R >= 4.1 and `nlme` (which ships with R). `haven` is only needed if
you read Stata `.dta` files; `testthat` only for the test suite.

Re-run the install command after pulling changes — there is no auto-reload. If you
are editing the package itself, `devtools::load_all("<path to your clone>/slopepower")`
is faster.

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

Both objects carry `n` and `n_per_arm` regardless; the printed summary shows one
of them at a time, per arm by default (`per_arm = TRUE`), and `print(x, per_arm =
FALSE)` shows the trial total instead without recomputing anything.

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
trial_design(visits = c(0, 1, 2, 3), dropout = dropout_rate(0.05))  # 5% per unit time
```

`dropout_rate(rate, per = 1)` states a constant withdrawal rate once and lets
`trial_design()` expand it to `(rate / per) * diff(visits)` — one proportion per
interval, whatever the schedule. It yields incremental proportions by
construction, so it cannot be combined with `dropout_type = "cumulative"`.

Participants who attend baseline only carry no slope information; a non-zero
first element warns.

Whatever the proportions, they are handled by the pattern-mixture method of
Dawson and Lagakos (1991, 1993), as in §2.5 of the paper: participants are
stratified by the visits they attend, each stratum is sized as if the whole trial
followed that pattern, and the strata are combined as the reciprocal of the
weighted mean of the reciprocals of those sizes. Withdrawers therefore still
contribute the visits they attended, which is less conservative than dividing a
completers-only sample size by the completion rate, and is the right adjustment
when the trial is to be analysed by a mixed model on all observed measurements.
Dropout is assumed monotone and unrelated to a participant's own trajectory —
every stratum shares the same variance components. `?trial_design` spells this
out.

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

Both return one row per visits × dropout combination, with identical columns,
including `n` and the anticipated visit burden `visits` — per arm by default
(`per_arm = TRUE`), or `per_arm = FALSE` for the trial total. Use `dropout_rate()`
rather than a fixed vector when the schedules being compared have different
numbers of visits: each cell expands it against its own schedule, so the same
withdrawal rate means the same thing in every row.

The assumptions are axes too. `effectiveness`, `alpha` and whichever of `n` and
`power` the grid is not solving for each take several values as readily as one, and
every value supplied is priced against every design — a sensitivity analysis in the
same table:

```r
slope_sample_size_grid(
  pars, visits = 0:3,
  dropout       = list(`0pc` = dropout_rate(0), `20pc` = dropout_rate(2 / 30)),
  power         = c(`80pc` = 0.8, `90pc` = 0.9),
  effectiveness = list(`20pc` = 0.2, `30pc` = 0.3)
)
```

Each of these axes reports its value in the column of the same name (`power`,
`alpha`, `effectiveness`, `n`); the names are what identify a cell in any message.

### The floor no design can beat

Before searching, it is worth knowing what the search can possibly achieve. The
treatment-effect variance has a greatest lower bound over *all* visit schedules,
`2 * (sigma2_slope - sigma_cov^2 / sigma2_intercept)`, so the sample size does
too:

```r
slope_sample_size_floor(pars, effectiveness = 0.33)   # the smallest N any design could need
slope_var_floor(pars)                                 # the limiting s*^2 behind it

# Or, for the settings of a result you already have:
ss <- slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
slope_sample_size_floor(ss)
```

Neither takes a visit schedule, because neither depends on one — and dropout can
only push the requirement up, so the bound covers designs with withdrawal too. If
the floor is already unaffordable, the visit schedule is not what needs to change.
The result inherits from `slope_result`, so its `as.data.frame()` row binds
together with the rows `as.data.frame()` gives for `slope_sample_size()` and
`slope_power()` results, marked `solve_for = "n_floor"`. (The grid functions
return a differently shaped table — one row per design, with the schedule
columns — so the floor row does not `rbind()` onto that.)

The bound is approached only as the schedule becomes both long *and* dense;
lengthening a two-visit trial only converges on the higher value
`2 * (sigma2_slope - sigma_cov^2 / (sigma2_intercept + sigma2_residual))`. The
`what-is-s-star` vignette derives all of this.

## Uncertainty in the stage-one estimates

```r
ss <- slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
slope_bootstrap(ss, R = 200, type = "bca", seed = 1)
```

Resamples subjects (stratified by group where relevant) and refits the mixed model
each time, so it is slow — start with a small `R`. Needs parameters from
`slope_params()`; manually supplied parameters carry no data to resample.

Hand it the result you want an interval around, not a fresh specification of the
calculation: the result already carries the design, effectiveness, target and
significance level it was solved with, and each replicate re-solves exactly that.
`slope_bootstrap()` dispatches on what it is given — a `slope_sample_size` object
for the required `n`, a `slope_power` object for the power achieved, or a
`slope_params` object for the fitted slope itself. Use `statistic = "tte"` on
either stage-two result for the target treatment effect behind it. Because each
method offers only quantities its object solved for or derived, bootstrapping one
of the calculation's own inputs — a zero-width interval dressed up as a result —
cannot be expressed.

## Example data

The three simulated datasets from Nash et al. (2021) ship with the package and
load lazily, so they are available by name after `library(slopepower)` with no
call to `data()` and no need for `haven`:

| | subjects | visits | grouping | paper |
|---|---|---|---|---|
| `slpower1` | 200 | 0, 1, 2, 3 years | none — decline against zero | pp. 588–589, Table 1 |
| `slpower2` | 500 | four calendar dates | `case` vs healthy control | pp. 590–593 |
| `slpower3` | 150 | 0, 0.5, 2 years | `treat` vs control arm | p. 594 |

See `?slpower1` and friends. `slpower2` records visit times as dates, so build the
time variable before fitting — `d$time <- as.numeric(d$vdate) / 365` for the
paper's units. They are rebuilt from the `.dta` files in the repository root by
`data-raw/make-data.R`, which asserts that the packaged copies still reproduce the
published slopes.

## Vignettes

```r
browseVignettes("slopepower")
vignette("introduction", package = "slopepower")
```

Both need the package to have been installed *with* its vignettes — see
[Installing the vignettes too](#installing-the-vignettes-too).

| | topic |
|---|---|
| `introduction` | A tour of the three stage-one situations `slope_params()` handles, worked through the three example datasets, ending in a full sample-size and power calculation |
| `from-stata` | A migration guide for existing Stata `slopepower` users: what maps directly, what changed, and what to watch for |
| `harmful-previous-trial` | What to reuse from a trial in which the treatment made things worse, why `target = "observed"` is the wrong tool there, and how much of the resulting sample size is bias rather than signal |
| `what-is-s-star` | What the "two-person trial" standard error actually is, why two is the minimum and not an approximation, the `s*/sqrt(N)` scaling derived rather than asserted, and a closed form that yields a floor on sample size no visit schedule can beat |

Vignettes are `rmarkdown::html_vignette`, so building them needs **pandoc** on
`PATH` — the flake's dev shell provides it. They set `math_method: mathml`, which
writes the equations into the file as MathML. The default would instead leave the
LaTeX as raw text and inject a loader that fetches MathJax from
`mathjax.rstudio.com` when the page is opened, so the equations would only render
online.
`harmful-previous-trial` runs a small Monte Carlo study and a bootstrap, and
takes around a minute to knit; `what-is-s-star` fits 200 mixed models to check
the `s*/sqrt(N)` scaling empirically, and takes around 20 seconds.

## Porting existing Stata scripts

`slopepower()` mirrors the Stata command's interface argument for argument —
`schedule` lists follow-up visits only, `dropouts` aligns with it, `scale` divides
time, the data type is declared with `obs`/`rct`, and a single call computes
either a sample size or a power depending on whether you pass `power` or `n`:

```r
slopepower(slpower1, "sdmt", "id", "visit", schedule = c(1, 2),
           obs = TRUE, nocontrols = TRUE, effectiveness = 0.33)
```

That last point is the one thing the wrapper keeps that new code should not: it
is Stata's bimodal interface, preserved deliberately for parity. Use `slopepower()`
for mechanical translation and parity checks; prefer the pipeline above —
`slope_params()`, `trial_design()`, then `slope_sample_size()` or `slope_power()` —
for new work.

## Notes

- `man/` and `NAMESPACE` are generated by roxygen2 from the `#'` comments in `R/` —
  do not edit them by hand. After changing those comments, re-run
  `roxygen2::roxygenise("<package path>")`. `CONTRACT.md` documents the
  mathematics and the deliberate divergences from the Stata original.
- The test suite reads the packaged datasets, so it runs in full under
  `R CMD check` from a built tarball. The one exception is
  `test-packaged-data.R`, which checks those datasets against the `.dta` files
  they were converted from and so skips outside a source checkout.
- `R CMD check` passes clean. `--run-donttest` additionally runs the
  `slope_bootstrap()` example, which fits a mixed model once per replicate and so
  takes appreciably longer than the rest.
