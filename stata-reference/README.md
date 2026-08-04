# Golden-reference tables from `slopepower.ado`

The published paper gives 15 numbers, and the R port reproduces all 15. That is
a strong check on the parts of the method the authors chose to demonstrate and
no check at all on the rest of the option space — irregular visit schedules,
dropout concentrated at one visit, `alpha` and `effectiveness` away from their
defaults, saturating power, `nocontvar`, `usetrt`, subsets that re-origin time.
Those are exactly the places a transcription slip survives.

These do-files run the original Stata command over that space and write CSVs the
R test suite reads back.

The CSVs need a Stata licence to produce, so a validated copy of each is
committed into the package as a test fixture:

```
Rscript make-fixtures.R      # stata-reference/*.csv -> ../slopepower/tests/testthat/fixtures/
```

The suite reads `tests/testthat/fixtures/` in preference to this directory, so
the parity tests run on any machine, with or without Stata. The copy is exact —
`export delimited` writes 17 significant digits and `make-fixtures.R` does not
reformat, because rounding would put the exact-arithmetic test out of reach. It
also refuses to write a partial set, checks every grid-5 fit converged, and
records a `MANIFEST` of row and column counts that the suite verifies, so a
regenerated `.do` file and a stale fixture cannot pass each other unnoticed.

Regenerating: run the `.do` files, then `make-fixtures.R`. Grid 5 must be
regenerated with the rest — the exact test joins its fits to the design grids on
dataset, model, scale, subset and `nocontvar`, so a stale grid 5 silently
compares against the wrong fit. If a row count changes on purpose, update
`EXPECTED_ROWS` in `make-fixtures.R`.

## Running

From this directory, with the three `.dta` files in the parent:

```stata
do 00_open_questions.do
do 01_grid_arithmetic.do
do 02_grid_comparators.do
do 03_grid_fits.do
do 04_edge_cases.do
do 05_variance_components.do
```

Only `00` opens a log, and it now starts with `capture log close _all`. On the
first run it aborted on a stray backtick before reaching `log close`, so the
other four files appended to its still-open log and question 2 never ran — if
`00_open_questions.log` comes out at hundreds of KB, that is what happened
again.

`00` writes a log to read by eye. The other five write `<stem>.dta` and
`<stem>.csv` beside themselves; the R suite reads the CSVs. Grid 5 must be
regenerated alongside the others: the exact test joins its fits to the design
grids on dataset, model, scale, subset and `nocontvar`, so a stale grid 5 either
drops rows from the join or compares against the wrong fit.

Stata refits the mixed model on **every** call — there is no way to reuse a fit
through the `slopepower` interface — so wall-clock time is dominated by the
number of rows times the cost of one fit. Rough shape:

| File | rows | dataset | cost per fit |
|---|---|---|---|
| `01_grid_arithmetic` | 348 | slpower1, 200 subjects, single group | cheap |
| `02_grid_comparators` | 88 | slpower2 (500 subjects, four variance components) and slpower3 | the 40 slpower2 rows dominate |
| `03_grid_fits` | 77 | all three | mixed |
| `04_edge_cases` | 71 | mostly slpower1; many rows error before fitting | cheap |
| `05_variance_components` | 37 | all three, one row per distinct fit | one fit per row, no repeats |

Around 580 calls in total. If that is too slow to iterate on, cut sweep B in
`01` first — it is the largest block and the most redundant, since alpha and
power enter only through `z_{1-α/2} + z_{1-β}`.

Run `01` first: if its CSV lands and the R suite is happy with it, the harness
is working end to end and the rest is just more of the same.

## What each file is for

**`00_open_questions.do`** — not a grid. Two questions about `slopepower.ado`
that were left open because there was no Stata to answer them with. **Both are
now answered, and both came back in the `.ado`'s favour** — see "What the run
found" below. The file is kept because it is the evidence, and because it
re-checks in a few seconds if either question is ever reopened.

1. `slopepower.ado:239` and `:276` both write `` `sched_length ' `` with a space
   inside the macro name. Probed three ways; Stata trims the name, so the guard
   at `:279` is live.
2. `slopepower.ado:371` computes the cases' residual variance as
   `exp(b[13]+b[14])^2`, which is right only if `residuals(independent, by())`
   stores a base log-SD plus an offset. The coefficient names in `e(b)` settle
   it outright, and they say base+offset.

**`01_grid_arithmetic.do`** — everything downstream of the fit, swept on the
cheapest dataset. Two sweeps rather than one cross: schedule × dropout pattern
at fixed scaling (where the matrix algebra lives), then effectiveness × alpha ×
(power | n) at two fixed designs (where the scaling lives, including the
`effectiveness` factor that appears both inside `tte` and again in the sample
size formula).

**`02_grid_comparators.do`** — the reference-slope rules. The one worth the most
is RCT data *without* `usetrt`: Stata estimates the treated arm's slope, prints
it, and then ignores it, measuring the effect toward zero instead. That is the
default, it is deliberate, and it is precisely the sort of thing a port
"corrects" on the way through. `effectiveness(1)` and `usetrt` differ only in
the reference slope, so any disagreement between that pair isolates cleanly.

**`03_grid_fits.do`** — everything that changes the fitted parameters. `scale()`
is the interesting lever: Stata divides time by `scale` before fitting and then
works on a unit integer grid, so `scale(s) schedule(j1 … jK)` is the only way
Stata can put a visit at a non-integer time, and it is the exact counterpart of
the R port's `visits = s * c(0, j1, …, jK)`. Since the port dropped `scale()`
deliberately, this is the sweep that tests that decision. Subsets are a free
extra dataset: each gives fresh variance components at no cost in new data, and
`if visit>0` also exercises the per-subject re-origining of time.

**`05_variance_components.do`** — added after the first run, and now the most
important file after `00`. Every other grid compares two REML fits and therefore
cannot be tighter than the two optimisers agree, which turned out to be ~1e-5 on
the variance components. This one posts the components themselves, so the R side
can feed Stata's own numbers into `slope_params_manual()` and require the
arithmetic to agree to machine precision — a test of Σ\*, X\* and F\* with the fit
removed from the comparison entirely. It also posts the whole of `e(b)` with its
column names, which settles open question 2 for every model at once: it computes
the cases' residual variance both the way `slopepower.ado:371` does and the way
the competing reading would, and prints them side by side against the R port's
10.354.

**`04_edge_cases.do`** — boundaries and guards, where `_rc` matters as much as
the numbers. The rows to read first are tagged `SAT-` (power saturation, where
Stata's back-solved `var_tte` goes **missing** — not zero — from `n = 10000`
upward while the R port's closed form stays finite), `SUM1-` (a dropout list
summing to exactly 1 in decimal, which both implementations accept), `ODD-`
(`n(451)` must come back as 450, and does), and `IGN-` (options Stata warns about
and then ignores, which the R port has no way to express).

## CSV layout

One row per call. Inputs first, so the R side can reconstruct the design without
parsing Stata syntax; then `rc`; then the ten columns of `r(table)` and the two
count scalars.

| Column | Meaning |
|---|---|
| `tag` | unique label within the file |
| `dataset` | `slpower1` / `slpower2` / `slpower3` |
| `model` | `obs_nocont` / `obs_cases` / `rct` (Stata models 2, 1, 3) |
| `scale` | the `scale()` argument |
| `subset` | the `if` condition, written so it is valid R as well as Stata |
| `sched` | the `schedule()` numlist |
| `drops` | the `dropouts()` numlist, empty for none |
| `effin` | the `effectiveness()` argument, empty when omitted or under `usetrt` |
| `usetrt`, `contvar` | 0/1 flags for `usetrt` and `nocontvar` |
| `alpha`, `mode`, `pin`, `nin` | `mode` is `power` (solve for N) or `n` (solve for power) |
| `rc` | `_rc`: 0 on success |
| `o_alpha` … `o_slope_comp` | `r(table)` columns 1–10 |
| `o_obs`, `o_subjects` | `r(obs_in_model)`, `r(subjects_in_model)` |

`o_slope` is the case / untreated slope; `o_slope_comp` is the control / treated
slope, and is **missing by construction** for model 2 and for model 3 without
`usetrt` even though a comparator slope was estimated.

## What the run found

All six files have been run (584 grid rows plus 37 fits). Summary, so nobody
re-derives it:

- **Given identical inputs the two implementations agree exactly.** Feeding
  Stata's own variance components from grid 5 into `slope_params_manual()`,
  all 542 comparable rows match: max |ΔN| = 0, power to 3.4e-15, `var_tte` to
  1.3e-14. Σ\*, X\*, F\*, the Dawson–Lagakos weighting, the `effectiveness`
  rescaling, both `var_tte` branches and all three comparator rules are exact.
- **Through the two fits, 534 of 542 match N exactly**, and every remaining
  difference is REML convergence rather than algebra — which is what the line
  above proves rather than assumes. On slpower1 the fixed effects agree to 1e-15
  while the variance components differ by ~1e-5, and that is not a coincidence:
  slpower1 is balanced, so the GLS slope is the OLS slope whatever the covariance
  structure, and the slope is insensitive to precisely the quantity the two
  optimisers disagree about. Tightening `nlme::lmeControl` moves R's `var_tte`
  from 8.36388645 to 8.36394022 against Stata's 8.36394172.
- **Grid 2 is a clean sweep**: 88/88 exact, including every `usetrt` and
  `nocontvar` row and every "ignore the treated slope" default.
- **Worst row throughout is `U3-visit<2`** — two timepoints, so the random-slope
  model has one degree of freedom per subject and the likelihood is nearly flat.
  It dominates the maximum of every quantity.
- **Open question 1 is answered, negatively.** Both a short and a long
  `dropouts()` list return `_rc = 198`, and the direct probe evaluates
  `` `counter ' + 1 `` to 8: Stata trims the trailing space in
  `` `sched_length ' ``, the counters are right, and the guard at `:279` is live.
  The R port's length check restates Stata rather than fixing it.
- **Open question 2 is answered: `slopepower.ado:371` is correct.** `mixed`
  stores the two residual parameters under `residuals(independent, by())` as
  `lnsig_e:_cons`, the base log-SD of the *reference* group (the controls), and
  `r_lns2ose:_cons`, a log-**ratio** offset — the `r_lns` prefix is what makes
  base+offset the right reading. So `exp(b13+b14)^2 = 10.354254`, matching
  `mixed`'s own printed `1: var(e) = 10.35425` for the cases, while
  `exp(b13)^2 = 10.698895` matches `0: var(e) = 10.6989` for the controls. The
  competing two-free-log-SDs reading would give 0.968, which is nonsense. Same
  result at `b[11]`, `b[12]` in the `nocontvar` branch.
- `dropouts(0.3 0.3 0.4)` is accepted by Stata (`N = 1418`), so the
  sum-to-one-in-decimal case is not a divergence either.
- Under dropout, solving for power, Stata's back-solved `var_tte` goes **missing**
  once power saturates at 1 — observed from `n = 10000` upward, finite at 2000.

Two things this exposed about the harness itself, both since fixed: Stata reports
the slope on the axis it fitted on, so the R side must multiply by
`scale / time_divisor` before comparing (`stata_axis_factor()`); and `treat` in
slpower3 is assigned by id block, so `if id<=75` leaves a single arm and
slopepower refuses it.

## Reading the results

Two test files consume this work.

`slopepower/tests/testthat/test-stata-reference.R` replays every row through the
R port and compares. It reports one failure per quantity per grid, listing the
offending rows, rather than several hundred separate expectations. It runs three
comparisons at three different strengths:

1. **Exact** (1e-12): Stata's own variance components in, so only the arithmetic
   is under test. All 542 comparable rows match — ΔN = 0, power to 3.4e-15.
2. **Fitted components** against grid 5, at tolerances calibrated to nlme's
   stock settings.
3. **End to end**, through both fits, which cannot be tighter than (2).

`slopepower/tests/testthat/test-stata-behaviour.R` states in R what the original
was found to do — the residual-variance parameterisation, the live dropout-length
guard, `2 * floor(n/2)`, the reference-slope rules, the boundary messages — with
the observed Stata values quoted inline. It needs no fixtures, so the knowledge
survives even if the tables are regenerated.

Each quantity is checked twice, with a loose bound on the worst row and a tight
bound on the median. A maximum alone would have to admit the ill-conditioned
fits and would then hide a uniform drift; a median alone would miss a single
broken case. The tolerances are calibrated to nlme's **defaults**, deliberately:
tightening `lmeControl` closes most of the gap, but changing what the package
does for every user to make a test look better is backwards.

Two things it does *not* do, on purpose:

- It does not assert agreement on `var_tte` for dropout rows solved for power.
  The two implementations define that number differently and the difference is
  real — see `CONTRACT.md` §5.6.
- It does not require the R port to accept everything Stata accepts. The
  `KNOWN_DIVERGENCES` vector in that file lists the exceptions, each with its
  reason. An unexplained divergence is a bug; an unexplained *entry* in that
  vector is worse.

When a row disagrees, the question is which implementation is right — not how to
widen a tolerance. The porting rule for this project is to follow the intent of
`slopepower.ado` rather than its letter, so a disagreement that traces to a
defect in the `.ado` gets recorded here and the R behaviour stands.

## What these tables cannot cover

- Non-integer visit times except through `scale()`. Stata's `schedule()` takes
  ascending integers ≥ 1 and nothing else, so a schedule like `c(0, 0.3, 1.7)`
  has no Stata counterpart at any single scale. The port supports it; only
  internal consistency checks can defend it.
- Unequal allocation. Neither implementation has it, and the R port left it out
  deliberately because it would break the two-person `s*` construction.
- Intermittent (non-monotone) missingness. Out of scope per the paper §2.5.
- `slope_bootstrap()`. There is no Stata counterpart at all.
