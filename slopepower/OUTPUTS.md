# Function outputs: structures and printed forms

What every exported function returns and prints, after the output
harmonisation on the `consistent-outputs` branch (2026-10-21). The rule the
names follow is CONTRACT.md section 4.0; `tests/testthat/test-output-names.R`
checks it across every output.

## The naming rule, in brief

- A count on a planned trial always comes as a pair, `<q>_total` and
  `<q>_per_arm`, total = 2 × per arm, both always stored: `n_*`,
  `n_requested_*`, `visits_*`. No bare `n`, `n_requested` or `visits`.
- A bootstrapped quantity `q` has six fields: `q`, `q_replicates`, `q_mean`,
  `q_sd`, `q_lower`, `q_upper` (`q` ∈ `n_total`, `n_per_arm`, `power`, `tte`,
  `slope`), and `<statistic>_ci_method`. No bare `ci`, `observed`,
  `boot_mean`, `boot_sd`, `se`.
- A scalar result's fields are named as its grid's columns.
- `_total` means the whole trial only (`dropout_cumulative`, not
  `dropout_total`); `n_` means participants (`replicates_failed`,
  `refits_failed`).
- `per_arm` is display only: it picks which twin `print()` shows, and the
  printed label or column always names the basis.
- Every print opens with a title line and a blank line.

## 1. Return types at a glance

| Function | Returns | Class | Print opens with |
|---|---|---|---|
| `slope_params()`, `slope_params_manual()` | list, 16 fields | `slope_params` | `Slope parameters: <kind>` |
| `slope_lme_control()` | `lmeControl` list | `list` | default |
| `fixef()`, `sigma()`, `slope_se()`, `slope_var()`, `slope_var_floor()`, `slope_effect_size()` | number(s) | – | default |
| `getVarCov()`, `vcov()`, `slope_sigma()` | matrix | – | default |
| `dropout_rate()` | list `rate`, `per`, `pattern` | `dropout_rate` | `Dropout rate: …` |
| `slope_sample_size()` | list, 15 fields | `slope_sample_size`, `slope_result` | `Sample size` |
| `slope_power()` | list, 17 fields | `slope_power`, `slope_result` | `Power` |
| `slope_sample_size_floor()` | list, 12 fields | `slope_sample_size_floor`, `slope_result` | `Sample size: lower bound over all visit schedules` |
| `slope_power_ceiling()` | list, 14 fields | `slope_power_ceiling`, `slope_result` | `Power: upper bound over all visit schedules` |
| `as.data.frame(<slope_result>)` | 1-row data.frame, 22 cols | `data.frame` | default |
| `slope_sample_size_grid()` | data.frame, 15 cols | `slope_sample_size_grid` | `Sample size grid (k rows)` |
| `slope_power_grid()` | data.frame, 17 cols | `slope_power_grid` | `Power grid (k rows)` |
| `slope_*_boot()` (5 functions) | list | `<fn>`, `slope_bootstrap` | `Bootstrapped <statistic>` |
| `confint(<slope_bootstrap>)` | matrix, rows by quantity | – | default |
| `slope_sample_size_grid_boot()` | data.frame | `…_grid_boot`, `…_grid`, `data.frame` | `Bootstrapped sample size grid (k rows)` |
| `slope_power_grid_boot()` | data.frame | as above, power | `Bootstrapped power grid (k rows)` |
| `confint(<grid_boot>)` | matrix, a row per cell | – | default |
| `[` on any grid | plain data.frame | `data.frame` | default |
| `slopepower()` | `slope_sample_size` or `slope_power` | as those | as those |

`result$design` is a `trial_design` (print: `Trial design`); `result$params`
is the full `slope_params`.

## 2. Structures

### `slope_result`
```
slope_sample_size:  n_total, n_per_arm, power, alpha, effectiveness, target,
                    tte, var_tte, effect_size, slope_difference, reference_slope,
                    visits_total, visits_per_arm, params, design
slope_power:        the same, with n_requested_total, n_requested_per_arm
                    after n_per_arm
floor / ceiling:    as their counterpart, minus visits_* and design
attr per_arm
```
`as.data.frame()`: `alpha, power, n_total, n_per_arm, n_requested_total,
n_requested_per_arm, effectiveness, target, tte, var_tte, effect_size, slope,
slope_comparator, comparator, reference_slope, slope_difference, n_follow_up,
visits_total, visits_per_arm, n_obs, n_subjects, solve_for`.

### Grids
```
design, dropout, scheduled_visits, last_visit, dropout_cumulative,
n_total, n_per_arm, n_requested_total, n_requested_per_arm, power, alpha,
effectiveness, tte, var_tte, effect_size, visits_total, visits_per_arm
attr per_arm
```
`n_requested_*` are in power grids only.
Bootstrap grids add `n_total_{mean,sd,lower,upper}`,
`n_per_arm_{mean,sd,lower,upper}` (or `power_*`), `n_ci_method` (or
`power_ci_method`), `tte_{mean,sd,lower,upper}`, `tte_ci_method`,
`replicates_failed`; attributes `R, requested_ci_method, statistic, level,
slope_se, refits_failed, named, per_arm, straddle, straddle_of, slope,
slope_replicates, slope_mean, slope_sd, slope_lower, slope_upper,
slope_ci_method`.

### Scalar bootstraps (here `statistic = "n"`)
```
statistic,
n_total,   n_total_replicates,   n_total_mean,   n_total_sd,   n_total_lower,   n_total_upper,
n_per_arm, n_per_arm_replicates, n_per_arm_mean, n_per_arm_sd, n_per_arm_lower, n_per_arm_upper,
n_ci_method,
slope,     slope_replicates,     slope_mean,     slope_sd,     slope_lower,     slope_upper,
slope_ci_method,
slope_se, R, level, replicates_failed, refits_failed, straddle, straddle_of, lattice
attr per_arm   (absent for slope_params_boot)
```
`power`/`tte` have one block in place of the two `n` blocks;
`slope_params_boot()` has the slope block only.

## 3. Printed forms

**`slope_result`**
```
Sample size

Data characteristics:
        number of observations in model = 800
        number of participants in model = 200
                         slope of cases = -1.672

Parameters for planned study:
                                  alpha = 0.050
                                  power = 0.800        | requested N per arm = 200
                          effectiveness = 0.330
  target treatment difference in slopes = 0.552
             number of follow-up visits = 2
                schedule (and dropouts) = 1 (0.05), 2 (0.05)

Estimated sample size:                                 | Estimated power:
                              N per arm = 388          |   power = 0.556
[Note: … the parameters' notes, as print.slope_params() gives them]
```
`per_arm = FALSE`: `N total = 776`, `requested N total = 400`.

**Scalar bootstrap**
```
Bootstrapped sample size

          estimate     mean       sd    lower   upper
n_per_arm      356   354.15  67.1556      258     513
slope      -1.6725 -1.69914 0.160655 -1.98959 -1.3985

  Bootstrap:   R = 40 replicates; 95% BCa intervals.
  Note:        0/40 (0.0%) bootstrap samples failed to refit.
  Note:        0/40 (0.0%) of replicates refit a slope on the opposite side
               of zero from the fitted one.
  Mean, SD:    each replicate is rounded up to a whole participant per arm
               before averaging, so the mean is not a runnable arm size.
```

**Grids**: title, then the table on the chosen basis under the stored names
(`n_per_arm`, `n_requested_per_arm` in a power grid, `visits_per_arm`, and for a bootstrap
grid `n_per_arm_mean` … `n_per_arm_upper`). A bootstrap grid hides the `tte_*`
columns (shown in a frame of their own when `effectiveness` varies), hides
`replicates_failed`, and shows `<statistic>_ci_method` only where a cell fell
back from the method asked for; its notes follow the table.

## 4. Status of the earlier inconsistencies

Resolved on this branch: the three meanings of `n`; bootstraps storing totals
but printing per arm; the four interval shapes (now `_lower`/`_upper`
everywhere, with `confint()` for scalar and grid bootstraps); the two field-name
conventions; the four header styles and the plain grids' missing basis
information; and results dropping the time-shift, common-variance and corSymm
notes. Fixed in passing: N labels ("N total", "requested N …"), heading
indentation, the nested parentheses in the `slope_params` title, the exposed
`added_cols` attribute,
`slope_params_boot()`'s meaningless `per_arm`, and `dropout_total`.

Still open, by choice: Stata's `.` for a missing value in the Stata-style
printouts; the order of the slope lines (`slope_params` lists the slopes then
the difference, a result the difference first); sign conventions (`tte`
positive, `effect_size` and `slope_difference` negative for a decline);
`as.data.frame()` for a bootstrap.
