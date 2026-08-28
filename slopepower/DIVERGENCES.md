# Where the R port diverges from `slopepower.ado` v2.1

Every item below is a deliberate divergence, made under the rule in
`CONTRACT.md` (and the project memory): port Nash et al.'s *mathematics*
line-for-line, but implement the *architecture* the authors would have
written with more time — fixing Stata-idiom artifacts and validation
defects, not the statistical method. Two Stata-only claims (the
`sched_length ` macro-name space, and the `mixed`/`nlme` reference-group
convention) were checked against a real Stata licence and found *not* to be
bugs; those are noted where relevant but are not divergences.

Claims that turn out **not** to be divergences are recorded here too rather
than quietly dropped, because the reason they look like one is exactly what a
future reader will re-derive otherwise: see §10 on the input guards Stata
already has, and "Claims checked and rejected" at the end for the four that
were probed under a licence and came back in the `.ado`'s favour. Nothing in
this document now rests on unverified behaviour.

Each entry gives the Stata original, what the R port does instead, and why.

---

## 1. Arbitrary real-valued visit times, no `scale()`

**Stata** restricts `schedule()` to ascending integers ≥ 1 and builds the
marginal covariance matrix `V` by looping over integer positions, then adds a
`scale()` option so a fractional visit spacing (e.g. six-monthly) can be
expressed as an integer schedule on a rescaled time axis:

```stata
syntax varname(max=1) [if], ... SCHEDule(numlist ascending integer >=1) ...
    [... SCAle(real 1) ...]
...
replace `time' = `time' / `scale'
...
forvalues i=1/`tpoints' {
    local v_`i'_`i' = (`var_int')+(`var_res')+(((`i'-1)^2)*(`var_slope')) + (2*(`i'-1)*`cov_slopeint')
}
```

**R** builds Σ directly at whatever real visit times are supplied, so
`scale()` has no equivalent:

```r
Sigma <- sigma2_intercept +
         outer(t, t) * sigma2_slope +
         outer(t, t, "+") * sigma_cov +
         diag(sigma2_residual, length(t))
```

`trial_design(visits = c(0, 0.5, 1, 1.5, 2))` is now legal directly; Stata
requires `schedule(1 2 3 4) scale(0.5)`. The two are proven exactly
equivalent against Stata's own output for the paper's six-month design
(N = 620 either way, `stata-reference` grid 6). `slopepower()`, the
backward-compatible wrapper, still accepts `scale` for callers porting old
Stata calls, but only divides `time` by it before fitting and emits a
`message()` suggesting the caller drop it in new code — it is not enforced
as an integer/scale pair.

---

## 2. Baseline visit is explicit, not implicit

**Stata**'s `schedule()` lists only the *follow-up* visits; visit 0 is
assumed silently:

```stata
* Make the first row 1, 0, 0...
matrix `visit_matrix'[1,1] = 1
```

**R** requires the caller to write baseline into `visits` explicitly
(`trial_design(c(0, 1, 2))`, not `trial_design(c(1, 2))`), and rejects a
vector that doesn't start at 0 with a suggested correction, rather than
inserting it silently. Rationale in `design.R`: the implicit convention is a
common source of off-by-one design errors, and Stata's own users have to
remember it rather than see it.

---

## 3. Sample size and power are two functions, not one bimodal one

**Stata** has a single command that switches mode on which of `power()` /
`n()` was supplied, with four cases to check by hand:

```stata
if "`power'"!="" & "`n'"!="" {
    display as error "Only one of power and n may be specified, not both"
    exit 184
}
...
if "`power'"=="" & "`n'"=="" {
    local given_power = 0.8 // Default power is 80%
    local power_or_n power
}
```

**R** splits this into two exported functions with disjoint signatures:

```r
slope_sample_size(params, design, effectiveness = 0.25, target, power = 0.8, alpha = 0.05)
slope_power(params, design, n, effectiveness = 0.25, target, alpha = 0.05)
```

`slope_power(..., power = 0.8)` is an "unused argument" error from R's own
argument matching — the ambiguity Stata resolves with a run-time branch is a
type error at the call site instead. The bimodal interface is reintroduced
**only** in the compatibility wrapper `slopepower()`, which mimics the
original single-command surface (and defaults `power` to 0.8 there,
matching Stata) for users porting existing Stata calls.

---

## 4. The model is inferred from the data supplied, not declared alongside it

**Stata** requires the scenario to be *declared* with `obs` / `rct` /
`nocontrols`, **and** the matching variable supplied separately, and then
hand-checks the ways the two can disagree:

```stata
if "`obs'"!="" & "`rct'"!="" { ... exit 198 }
if "`obs'"=="" & "`rct'"=="" { ... exit 198 }
if "`obs'"!="" {
    if "`nocontrols'" == "" local model 1
        else local model 2
}
else if "`rct'"!="" & "`nocontrols'"=="" local model 3
    else { dis as error "You cannot specify nocontrols and rct " ... }
```

with two further checks (`.ado:135`, `:198`) that the declared model actually
has the variable it needs, and two more (`:131`, `:195`) warning that a
supplied variable does not apply to the declared model.

**R** has no model declaration at all. The scenario *is* which grouping
argument was given:

```r
slope_params(y ~ t | id, data)                    # Stata's obs nocontrols  (model 2)
slope_params(y ~ t | id, data, healthy = case)    # Stata's obs case()      (model 1)
slope_params(y ~ t | id, data, treated = treat)   # Stata's rct treat()     (model 3)
```

`healthy` and `treated` are mutually exclusive, and that single check is all
that survives of the six above: the remaining inconsistent states cannot be
written down. The declaration lives on only in `slopepower()`, which
reproduces every one of Stata's checks because its job is to accept Stata
calls verbatim.

---

## 5. Cumulative dropout is a first-class input option

**Stata** only accepts incremental dropout fractions (proportion whose last
attended visit is visit `i`). Converting from a cumulative attrition curve
is left to the user.

**R**'s `trial_design()` accepts `dropout_type = "cumulative"` and converts
internally:

```r
trial_design(c(0, 1, 2, 3), dropout = c(0.05, 0.10, 0.15),
             dropout_type = "cumulative")
```

converted via `diff(c(0, cumulative))`, with validation that the cumulative
vector is non-decreasing and bounded by 1. The stored `dropout` field is
always incremental regardless of which form was supplied.

---

## 6. `dropout_rate()` — a per-time-unit rate, expanded per schedule

**Stata** has no equivalent: a "5% per year over 3 years" dropout has to be
manually converted to `dropouts(0.15)`, `dropouts(.05 .05 .05)`, or
`dropouts(.025 .025 .025 .025 .025 .025)` depending on which visit schedule
is being compared, with the arithmetic redone by hand for each.

**R** adds `dropout_rate(rate, per = 1, type = c("linear", "cumulative"))`,
which `trial_design()` expands for whatever schedule it is given. `type`
chooses which of two withdrawal patterns `rate` describes:

- `"linear"` (the default) applies `rate` to the *original* cohort —
  `(rate / per) * (visits[j+1] - visits[j])` per stratum:

  ```r
  trial_design(c(0, 1, 2, 3),        dropout = dropout_rate(0.05))  # 0.05 x 3
  trial_design(seq(0, 3, by = 0.5),  dropout = dropout_rate(0.05))  # 0.025 x 6
  trial_design(c(0, 3),              dropout = dropout_rate(0.05))  # 0.15
  ```

- `"cumulative"` applies `rate` as the proportion of *whoever is still in
  follow-up* who withdraws per `per` units of time, compounding
  geometrically rather than shrinking the original cohort by equal amounts:
  survival at time `t` is `(1 - rate) ^ (t / per)`.

  ```r
  trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.05, type = "cumulative"))
  # 0.05, 0.0475, 0.045125 -- 5% of whoever remains, each interval
  ```

Because the expansion happens in the constructor, the same object drives a
single design and every row of a Table-1 style comparison alike: the grid
functions expand it per cell through the same code path. Whichever `type` is
used, a rate produces incremental proportions by construction, so it cannot be
paired with `dropout_type = "cumulative"` (section 5) — an unrelated choice,
about how the vector itself was written down rather than how withdrawal
behaves over time — and that combination is rejected.

This is new surface area, but it computes exactly what a Stata user would
otherwise compute by hand for `dropouts()`, so it is not a change to the
statistical method.

---

## 7. Stage one (model fit) and stage two (power/N) are separated; grids fit once

**Stata** re-runs the whole `mixed` model from scratch for every row of a
comparison table — Table 1 of the paper (nine designs) means nine REML fits
of the *same* data.

**R** splits parameter estimation (`slope_params()`, one `nlme::lme` fit)
from the sample-size/power algebra (`slope_sample_size()`/`slope_power()`,
closed-form given the fitted variance components), so
`slope_power_grid()`/`slope_sample_size_grid()` fit once and evaluate many
designs against the same estimates:

> Because this port separates parameter estimation from the sample-size
> calculation, the mixed model is fitted once and every design is evaluated
> against the same estimates — the Stata original refits the model for each
> row of that table.

Because the fit is reused, a row of a grid costs closed-form algebra rather
than a REML fit, and the grids take that further than the paper's table
does: `effectiveness`, `alpha` and whichever of `n` and `power` the grid is
not solving for each accept several values as readily as one, and every
value is priced against every design. A sensitivity analysis over the
assumptions is then one call and one table rather than one table per
assumption. Nothing about a cell differs from the single-design call it
stands for — each is exactly `slope_power()`/`slope_sample_size()` at those
arguments — so this is interface, not method.

`slope_params_manual()` additionally lets a user hand in variance
components with no data at all (e.g. from a published paper, or from
Stata's own fitted numbers) — Stata has no such entry point; the data-in-memory
model fit is mandatory there.

---

## 8. Variance components extracted by name, not by matrix position

**Stata** pulls fixed and random effects out of `e(b)`, a 1×`k` matrix, by
hard-coded column index, and the index shifts between the three models and
the `nocontvar` variant:

```stata
scalar `var_slope'=(exp(`mbeta'[1,7]))^2      // Model 1
scalar `var_res'  =(exp(`mbeta'[1,13]+`mbeta'[1,14]))^2   // Model 1, with control variance
scalar `var_res'  =(exp(`mbeta'[1,11]+`mbeta'[1,12]))^2   // Model 1, nocontvar
```

Confirmed correct under real Stata: `mixed`'s `residuals(independent, by())`
stores `lnsig_e:_cons` as the *reference* group's (controls') log-SD and
`r_lns2ose:_cons` as a log-**ratio** offset, so `exp(b13+b14)^2` genuinely is
the cases' variance. But it is correct only because the column numbering is
memorised, not derived — get a model variant wrong and the position is
silently wrong too.

**R** extracts every component from the fitted object by name: the
random-effects block out of `nlme::getVarCov()`'s dimnames (`extract_re()` in
`params.R`), and the residual out of
`coef(fit$modelStruct$varStruct, allCoef = TRUE)` looked up by level name
(`extract_residual()`), never by positional index. Adding or reordering a
term therefore cannot silently mis-map a variance component. This also
sidesteps a live trap the R side of the same convention creates:
`nlme::varIdent` puts the *reference* level's SD in `sigma(fit)`, and because
the reference here is `control`, naive `sigma(fit)^2` would return the
controls' residual variance rather than the cases' — a 3.3% error the R port
avoids by extracting by name against the known reference level, not by
position.

---

## 9. `slope_params`/`trial_design` invariants are re-checked at every use, not only at construction

**Stata** validates once, in the `syntax`/quietly block at the top of the
command, coupled to that specific call's arguments. There is no persisted
object to re-validate later.

**R** makes `slope_params` and `trial_design` real objects, and re-checks
their invariants (positive-definite Σ and G, dropout length/sum, finite
visits, etc.) at the boundary of *every* stage-two call — `check_params()`
and `as_trial_design()` — not only where the constructors run. That is what
closes the gap for objects the constructors never saw: a hand-assembled
`slope_params` with a non-PD random-effects G matrix but a PD marginal Σ
(large `sigma2_residual` masks it) would otherwise slip through and silently
return a wrong large N, and a `trial_design` whose `$dropout` was edited after
construction would escape the baseline-dropout warning entirely.

---

## 10. Errors instead of Stata's `N = .` / silent missing result

**Stata** returns a missing value and prints nothing further when the slope
difference is zero (undefined sample size) — the command completes as if
nothing went wrong:

```stata
* (no explicit guard; `effsize` becomes 0, sampsize_arm divides by 0^2 -> missing)
```

**R** raises an explicit, named error:

```r
if (slope_difference == 0) {
  stop(sprintf(paste0("%s: the slope difference is exactly zero, so the target treatment ",
                      ...)), call. = FALSE)
}
```

The same holds for a singular information matrix (`treatment_effect_var()`),
a non-positive-definite implied covariance (`sigma_at()`), a non-PD
random-effects G matrix on a hand-built `slope_params` (`check_params()`, see
§9), an `alpha` so small that `1 - alpha/2` is not distinguishable from 1 in
double precision (`z_alpha()`, below), and a missing `n` in
`slope_power()`/`slope_power_grid()` — where Stata does not fail at all but
silently switches to solving for sample size at 80% power.

The `alpha` case is worth spelling out because the degeneracy is purely
arithmetic and the true answer is unremarkable. Both languages compute the
two-sided critical value as `invnormal(1 - 0.5*alpha)` / `qnorm(1 - alpha/2)`
(`.ado:634`); below about `2.2e-16` that argument rounds to exactly 1, and
`invnormal(1)`/`qnorm(1)` is infinite, so the sample size divides out to a
missing value in Stata and would be `Inf` in R. The requirement being solved
for is finite — 7,584 at `alpha = 1e-16` on the reference parameters, against
712 at `alpha = 0.05` — so reporting "no sample size is enough" is wrong in
both. R says why instead. The expression itself is *not* changed to the
numerically preferable `qnorm(alpha/2, lower.tail = FALSE)`: the two differ in
the last ulp for a large fraction of alphas — 22% sampled over (1e-6, 0.5),
rising to 72% over the conventional (0.001, 0.1), and 0.05 itself is one of
them — and with `ceiling()` downstream that is enough to move a published N by
one. Parity wins over the last bit of accuracy, per §1.

Related, and also an error rather than a silently different answer: the
`healthy`/`treated` group indicator must be constant within a participant
(`slope_params()`). Stata has no such guard, and neither language would
complain — every model reads the indicator row by row, so a participant whose
coding changes part-way through follow-up is simply fitted as a case for some
visits and a control for others. The fit converges and the result looks
ordinary; it just answers a different question. Group membership is a property
of the participant, not of the visit.

**Not** divergences, though they read like them, and listed here so the
distinction survives: Stata already errors on `alpha` outside (0,1)
(`.ado:67`), `power` outside (0,1) (`:112-119`), `effectiveness` outside
(0,1] (`:183-190`), `n < 2` (`:98`) and non-integer `n` (`:94`). The R port
restates those checks; it does not add them.

---

## 11. `common_variance` is a first-class three-way argument, not a two-state flag

**Stata**'s `nocontvar` is a bare flag: absent means "fit the full
uns-covariance random-effects block for controls", present means "reduce
it":

```stata
noCONTVar
...
if "`contvar'"=="" { // Controls are allowed a variance parameter
    capture mixed ... res(ind, by(`case')) reml iterate(`iterate')
    if _rc!=0 {
        dis as error "Model failure. Try using the nocontvar option, or see error code (below)."
        exit _rc
    }
```

**R**'s `common_variance` takes `NULL` / `TRUE` / `FALSE`, of which Stata's
flag expresses two. `TRUE` forces the reduction — Stata's `nocontvar`
present. `FALSE` forces the full structure and errors if it does not
converge — Stata's `nocontvar` absent, which is the same contract, since the
original prints `"Model failure. Try using the nocontvar option..."` and
exits with the underlying return code rather than degrading silently.

What Stata has no way to say is the R default, `NULL`: try the full
structure, fall back to the reduced one with a `message()` if it fails to
converge, and record which was used on the returned object (the
`common_variance` field, shown by `print.slope_params()`). Stata's user has
to read the error and retype the command with the flag.

---

## 12. `target = "observed"` vs. Stata's `usetrt`, kept as a named enum

**Stata** encodes "which slope is the reference" as a bare flag (`usetrt`)
whose meaning is documented only in a comment, and which silently sets
`effectiveness` to 1 and then to missing as a side effect:

```stata
scalar `slope0' = 0 // Default option is to ignore the slope for treated
if "`usetrt'"!="" {
    scalar `slope0'=`mbeta'[1,1]
    local effectiveness = 1
    ...
}
...
if "`usetrt'"!="" local effectiveness = .
```

**R** models the same choice as an explicit `target` argument to
`slope_sample_size()`/`slope_power()` (`"effectiveness"` or `"observed"`),
recorded on the *result* object and documented in `CONTRACT.md` §5.3 as a
truth table:

```
comparator == "treated" && target == "effectiveness" -> reference_slope = 0
comparator == "treated" && target == "observed"      -> reference_slope = slope_comparator,
                                                          effectiveness forced to 1
```

It is deliberately **not** a field of `slope_params`, which carries only
`comparator`: which study the parameters came from is a property of the data,
what you are targeting is a property of the question, and one fit can answer
both.

The divergence is in the *input*: supplying `effectiveness` alongside
`target = "observed"` is an error (`check_target_effectiveness()`), where
`.ado:454` quietly overwrites it to 1. What the two report afterwards agrees
on purpose — the R result's `effectiveness` is `NA_real_` under that target,
which is `.ado:648`'s missing value faithfully preserved, and
`test-stata-behaviour.R` pins it as agreement rather than divergence.

---

## 13. Fields named for what they are, not their column position or Stata option name

**Stata** names its return scalars after the model-specific role
(`slope_controls`, `slope_cases`, `slope_untreated`, `slope_treated`), so the
same underlying "reference slope vs. index slope" concept has four different
names across the three models, and the on-screen legend has to be
cross-referenced to know which is which for a given `model`. E.g. Model 1's
report is:

```stata
return scalar slope_controls = `slope0'
return scalar slope_cases = `slope2'
```

while Model 3 reports `slope_untreated`/`slope_treated` for the same
positions.

**R** uses one pair of field names throughout — `slope` (the
untreated/case/index group) and `slope_comparator` (healthy control or
treated arm, `NA_real_` if there is none) — tagged by a separate
`comparator` field (`"none"`/`"healthy"`/`"treated"`) instead of the naming
itself carrying the model identity.

---

## 14. Warnings are `warning()`/`message()` conditions, classed and deduplicated, not `display as error`

**Stata** prints every advisory the same way (`dis as error "WARNING: ..."`),
indistinguishable at the call site from a genuine error except by reading the
text, and a warning fires again on every single row of a hand-rolled grid
loop with no way to collapse repeats.

**R** raises real R condition objects with a class —
`slopepower_tte_direction` (`power.R`) and `slopepower_baseline_dropout`
(`design.R`) — so callers can `tryCatch`/`withCallingHandlers` on the
specific condition rather than string-matching, and the grid functions
(`slope_power_grid()`, `slope_sample_size_grid()`) collect both by class and
report each once per grid instead of once per design, the way a naive
per-row Stata loop would.

Only those two carry a class, and deliberately: they are the ones a
legitimate grid sweep is *expected* to trip on most cells. The rest of the
port's guards — including the `dropout` length check, which is a hard error
rather than a warning — are ordinary `stop()`s and `warning()`s.

---

## 15. Time-origin shift is a `message()`, is optional, and is recorded on the object

**Stata** always shifts each subject's time to start at 0 and reports it as
a same-severity `WARNING` line mixed in with genuine errors:

```stata
sum `first_time'
if (r(mean)!=0) | (r(sd)!=0) dis as error "WARNING: time variable did not start at zero for all participants. Times have been adjusted such that the first visit for each person is treated as time zero."
```

**R** does the identical shift by default, but reports it as a `message()`
(not `warning()`), and records that it happened on the returned object
(`time_shifted = TRUE`) so a caller can detect it programmatically instead of
grepping console text:

```r
if (any(abs(first) > 1e-12)) {
  time_shifted <- TRUE
  message(sprintf(paste0("%s: time did not start at zero for all subjects. ...")))
}
```

R also lets a caller decline the shift altogether with `origin = "none"`,
which Stata has no equivalent for — there the rebasing is unconditional. The
default, `origin = "subject"`, is Stata's behaviour.

---

## 16. Bootstrap resampling scheme is fixed by construction, not four independent options

**Stata**'s bootstrap recipe (documented in the paper, §2.6 and §4.1.2, not
built into the command) requires the user to correctly assemble `cluster()`,
`idcluster()`, `strata()` and `jack()` options to `bootstrap:`, and warns
that it silently assumes no observations were excluded — four places to get
wrong with no structural guard.

**R**'s `slope_bootstrap()` has no equivalent options: subjects are always
the resampling unit, replicate identifiers are always freshly generated (so
a subject drawn twice counts as two people rather than colliding), and
stratification by group is automatic whenever `comparator != "none"`. There
is nothing to misconfigure because there is nothing to configure.

Every replicate also refits the random-effects structure the observed fit
ended up with. Stata gets this for free: `nocontvar` is part of the command
string the `bootstrap:` prefix re-runs, and there is no automatic fallback
for it to disagree with. R has one (divergence 11), so the structure a fit
ended up with is not always the one that was asked for, and it has to be
read off `params$common_variance` and pinned — otherwise a replicate that
could not fit the full structure would quietly fall back to the reduced one
and land in the same interval as a point estimate fitted the other way. A
replicate that cannot fit the pinned structure is counted in `n_failed`
instead. This is visible only under `healthy`, and only through the
controls' slope: the model factorises per group, so the case estimates are
invariant, but `slope_comparator` — what `slope_difference` is measured
against — is not.

`slope_bootstrap()` also warns when the fitted slope is less than 2.5× its
standard error, per the paper's own §2.6 recommendation — a check the Stata
side leaves to the user to remember to apply by eye.

---

## 17. `slope_var()`/Σ helpers are pure functions of `(params, visits)`, not tied to data in memory

**Stata**'s whole calculation is a `program` operating on "the data
currently in memory", per the file's own header comment (*"The data in
memory will not be altered by this program"* — `preserve`/`restore` is used
precisely because the command's working style is to mutate the loaded
dataset in place and then undo it).

**R** has no notion of "current data" at all in the stage-two layer:
`slope_var()`, `slope_effect_size()`, `slope_sample_size()`, `slope_power()`
and the grid functions take a `slope_params` object and a `trial_design`
object as plain arguments and return a plain value or list — they can be
called from a script with no dataset loaded, called from `slope_params_manual()`-built
inputs, or vectorised over many designs without any `preserve`/`restore`
bookkeeping.

---

## 18. `effectiveness` / `usetrt` mutual exclusivity is a hard error, not a warning that then proceeds

**Stata** warns and continues when incompatible options are combined at the
wrong model (e.g. `casecon` given for a non-Model-1 call), but for the
genuinely contradictory `usetrt` + `effectiveness` combination it does
error:

```stata
if "`usetrt'"!="" & "`effectiveness'"!="" {
    dis as error "Only one of usetrt and effectiveness may be specified, not both"
    exit 184
}
```

**R** keeps this as an error in `slopepower()` too, but generalises the
*warn-and-reset* half of the same pattern (`casecon`/`treat` supplied for the
wrong model, `usetrt` requested for non-RCT data) into one shared helper
(`warn_unused_arg()`) instead of three hand-written duplicate blocks, so the
"does this option apply to this model" logic and its wording live in one
place rather than being copy-pasted per option as in the `.ado` file.

One member of that family carries a typo of its own. `.ado:196`, the warning
for `treat()` supplied with observational data, reads

```stata
dis as error "Treatment can only be specified with RCT data`. Treatment variable will be ignored."
```

— with an unbalanced macro quote in the middle of the string, where its three
siblings have none. `slopepower()` warns and continues, which is what the
author meant.

**Settled, 2026-08-10** (`stata-reference/07_open_questions_2.do`, question
4): Stata prints the backtick as a literal character and carries straight on
— `_rc = 0`, N = 712, identical to the well-formed sibling at `.ado:132`. So
the blemish is cosmetic, nothing downstream is consumed, and `slopepower()`'s
warn-and-continue is a faithful reading of the original rather than a repair
of it. The only difference left is the stray character in the message text.

---

## 19. `var_tte` under dropout is back-solved in a form that survives saturated power

**Stata** reports `var_tte` as `FSTAR[3,3]` when there is no dropout, and
otherwise back-solves it by inverting the sample-size formula, because with
dropout no single s*² applies across strata:

```stata
mat `RESULTS'[1,8] = `sampsize_arm' * (`tte')^2 / (invnormal(1-`alpha'/2)+invnormal(`power'))^2
```

It uses that one expression in both directions — when `power` is the user's
input, and when `power` is the number it has just computed itself at `:645`.

**R** keeps Stata's form when solving for N, and uses the algebraically
identical but numerically stable form when solving for power:

```r
var_tte <- if (!comp$design$has_dropout) comp$var_full
           else if (solving_for_n)       n_per_arm * comp$tte^2 / z_sum_sq
           else                          comp$tte^2 / scaled_effect^2
```

The two agree exactly. `power = pnorm(scaled_effect * sqrt(n_per_arm) - z_a)`
makes `z_a + qnorm(power)` equal `scaled_effect * sqrt(n_per_arm)`, so
`n_per_arm` cancels out of the ratio entirely. Round-tripping through the
power anyway costs precision, and past a few thousand per arm costs
correctness: the power saturates at exactly 1, and the inverse normal of 1 is
not a number in either language, so Stata's expression yields a missing
`var_tte` for a calculation that otherwise succeeded (R's earlier literal
transcription of it silently reported 0, `qnorm(1)` being `Inf`).

Solving for N there is no such cancellation — `n_per_arm` has been rounded up
to a whole participant — so that branch keeps Stata's literal form, rounding
included, and reports the same very slightly inflated number Stata does.

---

## 20. Convergence control is fixed and inspectable, not an option

**Stata** exposes `ITERate(integer 16000)` and threads it into every `mixed`
call, leaving the optimiser otherwise at its defaults.

**R** has no such argument. The settings behind every fit live in one
exported, callable function, `slope_lme_control()` — more iterations than
`nlme`'s defaults, `opt = "optim"`, and tighter tolerances, chosen because
the untightened defaults converged less precisely on the two-block
random-effects structure fitted for `healthy`. It is exported so that the
settings behind a published number can be inspected and reproduced outside
the package, **not** so they can be overridden inside it: the model this
package fits is fixed (see "What these models do and do not include" in
`?slope_params`), and a tunable optimiser would let two runs of the same call
disagree with nothing on the returned object to say why.

This is the one place where the port's numbers can differ from Stata's at
all. `stata-reference/` records the measured size of that difference — every
end-to-end disagreement traces to REML convergence, and feeding Stata's own
variance components in through `slope_params_manual()` reproduces all 542
comparable rows exactly.

---

## 21. Two features Stata doesn't have were deliberately **not** added

Not divergences by omission — checked against the paper and left out on
purpose, because adding them would go beyond what Nash et al. specify:

- **Unequal treatment allocation.** The two-person treatment-effect-variance
  trick (`Sigma_star = blockdiag(Sigma, Sigma)`, one notional person per arm)
  depends on 1:1 allocation to work; generalising it to unequal allocation
  is out of scope per the transpilation's own ground rule and is not implied
  by anything the paper computes.
- **Intermittent (non-monotone) missingness.** Both Stata and R assume
  dropout is monotone — once a participant misses a visit they are treated
  as gone for good. The Dawson–Lagakos pattern-mixture approach in §2.5 of
  the paper only covers monotone dropout; a participant who misses one visit
  and returns has no representation in either implementation, and adding one
  is explicitly out of scope.

---

## 22. The published "×4 shortcut" figure (1,272) is reproduced *and* corrected, with both numbers kept

Not a code-level Stata-vs-R difference (this is about a number quoted in the
paper's prose, not one of the paper's 15 pinned results), but worth
recording because it is a place the port's output disagrees with the paper
on purpose. §4.1.3 says powering for a fraction `p` of a previously observed
effect scales `N` by `p^-2`, and quotes **1,272** for `p = 0.5` (318 × 4).
That shortcut multiplies the paper's *already-`ceiling()`-rounded* N by 4,
which magnifies the rounding. Rounding last (or halving the effect at the
source variance components) gives **1,266**. The R port's test suite pins
*both* numbers and asserts the gap is always upward and bounded by `p^-2` per
arm — it does not silently correct the paper's 1,272, and does not silently
accept it as the "right" answer either.

---

## 23. `slope_sample_size_floor()` — a bound the paper does not state

**Stata** has nothing like it, and neither does the paper: §2.2's `s*` is
always evaluated at a stated `schedule()`, and Table 1 explores schedules
one at a time.

**R** exports `slope_var_floor()` and `slope_sample_size_floor()`, the
greatest lower bound of `s*^2` — and hence of `N` — over *every* visit
schedule. It falls out of the closed form `s*^2 = 2 / (t' Sigma^-1 t)`
derived in the `what-is-s-star` vignette, and is
`2 * (sigma2_slope - sigma_cov^2 / sigma2_intercept)`, twice the variance of
a participant's random slope given their random intercept.

This is an *addition*, not a change: no existing number moves, and the port
still computes everything the paper computes the way the paper computes it.
It is included because the Table 1 exercise — try nine schedules, see which
is affordable — has no stopping rule without it. The floor says when the
search is pointless: on `slpower1` at 33% effectiveness it is 236, against
712 for the paper's own two-year design, so a trial that cannot afford 236
cannot be rescued by any visit schedule.

Two things kept it honest rather than clever. The sample size goes through
`size_per_arm()`, the same equation (6) `solve_slope()` calls, so the bound
and the thing it bounds cannot drift apart. And it takes no `design`
argument in either method (CONTRACT.md §4.3) — the value does not depend on
one, and dropout can only raise it, so the bound covers designs with
withdrawal too.

---

## 24. `slope_sample_size_grid_boot()` bootstraps a whole grid from one resampling pass

**Stata** has no equivalent of the grid functions at all (divergence 7 —
Table 1 is built by hand, one `slopepower` call per row), so it has no
equivalent of bootstrapping one either. The paper's own worked bootstrap
(§4.1.2) is a `bootstrap:` prefix around a single `slopepower` call, refit
completely for that one design; doing the same nine times over, once per row
of Table 1, is the only way the Stata recipe could be extended to a grid,
and it refits the stage-one model `n_cells * R` times.

**R** shares one set of resampled replicates across every cell instead of
bootstrapping each design independently. The resampling scheme —
`boot_setup()`, `boot_replicate_matrix()`, `jackknife_values()`, all
`R/bootstrap.R` — depends only on the stage-one `params`, never on the
design being priced, so refitting once prices every cell: `R` replicates
plus one leave-one-subject-out jackknife, not `n_cells` times that. A
nine-cell grid at the default `R = 999` costs about what one
`slope_bootstrap()` call does, not nine times it.

This is not only cheaper — it changes what the intervals mean together.
Every cell's interval is built from the *same* draws, so two designs'
intervals move together with whatever the resampling happened to do to the
slope on that draw. The comparison between two rows of the table is
therefore paired, the way it would be if the same nine hundred simulated
datasets were run through nine designs, rather than confounded with two
independent samples' worth of Monte Carlo noise the way `n_cells` separate
`bootstrap:` prefixes would be. Choosing between designs is exactly when
this matters: Table 1 exists to compare rows, not to read one at a time.

---

## Claims checked and rejected

Things that look like divergences in the `.ado` source and are not, recorded
so they are not "rediscovered" a third time. Every one was settled by running
Stata, not by reading it harder; the evidence is in `stata-reference/`.

- **The dropout-total guard does not reject a legal list.** `.ado:265`
  accumulates by repeated subtraction, and `1 - 0.3 - 0.3 - 0.4` is
  `-5.551e-17` as a double, so `dropouts(0.3 0.3 0.4)` — legal, and exactly 1
  in decimal — looked certain to trip the `< 0` guard at `:268`. It does not.
  A Stata local round-trips through a decimal *string* rather than carrying
  the double: `1 - 0.3` stores as `.7`, then `.7 - 0.3` stores as `.4`, not
  `0.39999999999999997`, and the residue lands on exactly 0. Stata accepts
  the list and returns N = 1418, which is what this port returns.
  (`07_open_questions_2.do` Q3; the two `%21x` lines it prints differ, which
  is the whole proof.) The port's `1e-8` tolerance therefore *restates*
  Stata's behaviour rather than repairing it. It is still the right way to
  write the check in R, where `sum()` does carry the full double — but it
  fixes nothing, and the guard at `:268` still fires correctly on a total
  genuinely over 1 (`dropouts(0.5 0.3 0.3)` → `_rc = 198`).
- **The unbalanced macro quote at `.ado:196` is cosmetic.** It prints the
  backtick literally and carries on; see §18.

- **`dropouts()` validation happens before the fit in Stata too.** The
  length and total guards are at `.ado:268` and `:279`, inside the
  syntax-section `quietly` block that closes at `:294` — well before the DATA
  section (`:301`) and the MODEL section's `mixed` call (`:340`). The R port
  briefly had this backwards, running the REML fit in `slopepower()` before
  validating `dropouts`; that was fixed, and the fix restores Stata's
  ordering rather than improving on it.
- **The `sched_length ` macro-name space is harmless.** `.ado:239`/`:276`
  write `` `sched_length ' `` with a trailing space inside the macro name.
  Stata trims it, both counters are correct, and the guard at `:279` fires on
  a mismatched list in both directions (`stata-reference/00_open_questions.log`).
- **`exp(b13+b14)^2` really is the cases' residual variance.** See §8.
