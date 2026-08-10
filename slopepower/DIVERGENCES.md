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
future reader will re-derive otherwise: see §11 on the input guards Stata
already has, and the preamble to §5 on the one unverified claim left in this
document.

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

## 5. Dropout sum tolerance

**Stata** checks only that dropouts don't exceed 100%, via subtraction with
no tolerance:

```stata
local dfrac_complete = `dfrac_complete' - `fr'
if `dfrac_complete' < 0 {
    display as error _n "Dropouts cannot exceed 100%"
    exit 198
}
```

**R** sums first and compares the total against 1 with an explicit `1e-8`
tolerance (`DROPOUT_TOL`, `utils.R`), rather than accumulating by repeated
subtraction. Repeated subtraction is what makes the Stata form fragile:
`1 - 0.3 - 0.3 - 0.4` evaluates to `-5.551e-17` in double precision and would
trip a bare `< 0` guard, though `c(0.3, 0.3, 0.4)` is a legal dropout vector
summing to exactly 1 in decimal.

Two caveats, both worth stating because this is the one entry in the document
resting on unverified behaviour. R's `sum()` happens to return exactly 1 for
that particular vector, so the tolerance is not what rescues *this* example —
it is what stops the guard depending on the order the elements arrive in.
And whether Stata's local-macro arithmetic really rejects
`dropouts(.3 .3 .4)` was never put to the licence, unlike the two claims in
the preamble; what is recorded here is the *form* of the check, which is
plain in the source either way.

---

## 6. Cumulative dropout is a first-class input option

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

## 7. `dropout_rate()` — a per-time-unit rate, expanded per schedule

**Stata** has no equivalent: a "5% per year over 3 years" dropout has to be
manually converted to `dropouts(0.15)`, `dropouts(.05 .05 .05)`, or
`dropouts(.025 .025 .025 .025 .025 .025)` depending on which visit schedule
is being compared, with the arithmetic redone by hand for each.

**R** adds `dropout_rate(rate, per = 1)`, expanded automatically for each
candidate schedule when building a design grid — `(rate / per) *
(visits[j+1] - visits[j])` per stratum — so the same object drives every row
of a Table-1 style comparison without hand recomputation. This is new surface
area, but it computes exactly what a Stata user would otherwise compute by
hand for `dropouts()`, so it is not a change to the statistical method.

---

## 8. Stage one (model fit) and stage two (power/N) are separated; grids fit once

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

`slope_params_manual()` additionally lets a user hand in variance
components with no data at all (e.g. from a published paper, or from
Stata's own fitted numbers) — Stata has no such entry point; the data-in-memory
model fit is mandatory there.

---

## 9. Variance components extracted by name, not by matrix position

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

## 10. `slope_params`/`trial_design` invariants are re-checked at every use, not only at construction

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

## 11. Errors instead of Stata's `N = .` / silent missing result

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
§10), and a missing `n` in `slope_power()`/`slope_power_grid()` — where Stata
does not fail at all but silently switches to solving for sample size at 80%
power.

**Not** divergences, though they read like them, and listed here so the
distinction survives: Stata already errors on `alpha` outside (0,1)
(`.ado:67`), `power` outside (0,1) (`:112-119`), `effectiveness` outside
(0,1] (`:183-190`), `n < 2` (`:98`) and non-integer `n` (`:94`). The R port
restates those checks; it does not add them.

---

## 12. `common_variance` is a first-class three-way argument, not a two-state flag

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

## 13. `target = "observed"` vs. Stata's `usetrt`, kept as a named enum

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

## 14. Fields named for what they are, not their column position or Stata option name

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

## 15. Warnings are `warning()`/`message()` conditions, classed and deduplicated, not `display as error`

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

## 16. Time-origin shift is a `message()`, is optional, and is recorded on the object

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

## 17. Bootstrap resampling scheme is fixed by construction, not four independent options

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

`slope_bootstrap()` also warns when the fitted slope is less than 2.5× its
standard error, per the paper's own §2.6 recommendation — a check the Stata
side leaves to the user to remember to apply by eye.

---

## 18. `slope_var()`/Σ helpers are pure functions of `(params, visits)`, not tied to data in memory

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

## 19. `effectiveness` / `usetrt` mutual exclusivity is a hard error, not a warning that then proceeds

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

One member of that family is not reproduced faithfully, on purpose.
`.ado:196`, the warning for `treat()` supplied with observational data,
reads

```stata
dis as error "Treatment can only be specified with RCT data`. Treatment variable will be ignored."
```

— with an unbalanced macro quote in the middle of the string. What Stata
actually does with that has not been put to the licence; it is at best a
mangled message and at worst an abort where the surrounding code plainly
intends a warning. `slopepower()` warns and continues, which is what the
author meant.

---

## 20. `var_tte` under dropout is back-solved in a form that survives saturated power

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

## 21. Convergence control is fixed and inspectable, not an option

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

## 22. Two features Stata doesn't have were deliberately **not** added

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

## 23. The published "×4 shortcut" figure (1,272) is reproduced *and* corrected, with both numbers kept

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

## Claims checked and rejected

Things that look like divergences in the `.ado` source and are not, recorded
so they are not "rediscovered" a third time:

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
- **`exp(b13+b14)^2` really is the cases' residual variance.** See §9.
