# Layer 3 (continued) -- the limit of the layer-3 calculation over all designs:
# the sample-size floor and the power ceiling.
#
# Everything in power.R answers "what does *this* schedule cost?". This file
# answers "what does the best conceivable schedule cost?", which turns out to
# have a closed form that mentions no visit times at all. The two are the same
# arithmetic downstream of the variance: `size_per_arm()` in power.R is called
# from both, so the bound and the thing it bounds cannot drift apart.

# ---------------------------------------------------------------------------
# the limiting treatment-effect variance
# ---------------------------------------------------------------------------

#' The closed form, given already-validated parameters
#'
#' Positive whenever the random-effects covariance matrix is positive definite,
#' since `sigma2_slope - cov_intercept_slope^2 / sigma2_intercept` is its Schur
#' complement and positive definiteness makes the determinant
#' `sigma2_intercept * sigma2_slope - cov_intercept_slope^2` positive. That is enforced
#' by `check_re_covariance()`, which every route into a `slope_params` object
#' runs and which `check_params()` re-runs on hand-built ones -- so there is no
#' zero or negative branch to guard here.
#' @noRd
var_floor <- function(params, context) {
  if (residual_on_grid(params$residual)) {
    unstructured <- identical(params$residual$correlation, "corSymm")
    stop(sprintf(paste0(
      "%s: the residual structure (%s) is defined only at the visit times %s, so there is ",
      "no floor over all schedules to report: the schedules it can price are subsets of ",
      "those times%s. The smallest treatment-effect variance it allows is at all of them: ",
      "slope_var(params, c(%s))."),
      context,
      paste(c(if (unstructured) "corSymm()",
              if (!is.null(params$residual$sd_ratio)) "varIdent()"), collapse = " with "),
      label_numeric(params$residual$times),
      if (unstructured) paste0(", and under corSymm() the variance components the floor ",
                               "is built from are not separately identified") else "",
      label_numeric(params$residual$times)),
      call. = FALSE)
  }
  2 * (params$sigma2_slope - params$cov_intercept_slope^2 / params$sigma2_intercept)
}

#' Smallest treatment-effect variance any visit schedule can achieve
#'
#' The greatest lower bound of [slope_var()] over every possible set of visit
#' times:
#'
#' \deqn{\inf_t s^{*2} = 2\left(\sigma^2_b - \frac{\sigma^2_{ab}}{\sigma^2_a}\right)
#'       = 2\,\mathrm{Var}(b_i \mid a_i).}
#'
#' Twice the variance of a participant's own random slope given their random
#' intercept: what is left of the between-person spread in slopes once the
#' baseline measurement has told you everything it can. Since sample size is
#' proportional to \eqn{s^{*2}}, this is where the sample-size floor
#' [slope_sample_size_floor()] comes from.
#'
#' @details
#' # No design argument
#'
#' There is deliberately none. The bound does not depend on the visit schedule
#' --- that is the whole content of it --- and it does not depend on the dropout
#' pattern either, since dropout can only ever raise the variance. So the value
#' returned bounds every design the stage-two functions can express, not just
#' the ones without withdrawal.
#'
#' # It is an infimum, not a minimum
#'
#' No finite schedule attains it. Writing \eqn{\Sigma = \Sigma_0 + R}, with
#' \eqn{R} the residual covariance (\eqn{\sigma^2_\epsilon I} for independent
#' residuals), every contrast \eqn{c} has
#' \eqn{c^{\mathsf T}\Sigma c > c^{\mathsf T}\Sigma_0 c}, so
#' `slope_var(params, visits)` is strictly greater than this value for any
#' `visits`, however long or however dense. The gap closes as the number of
#' visits grows without bound.
#'
#' # Residual structures
#'
#' The value is the same with a serially correlated residual -- `corAR1()`,
#' `corCAR1()`, `corExp()` or `corGaus()` in [slope_params()] -- because each
#' of those correlations dies away with the time between visits, so a long
#' enough schedule still averages the residual out. It is approached more
#' slowly: visits packed close together are correlated and repeat each other,
#' so density alone gains less than it would with independent residuals.
#'
#' A structure with a parameter per visit, `corSymm()` or `varIdent()`, has
#' no floor and is refused: it is defined only at the visit times it was
#' fitted at, so the schedules it can price are subsets of those, and the
#' smallest variance among them is simply `slope_var()` at all of them.
#'
#' Lengthening a schedule is not enough on its own. Two visits a distance
#' \eqn{t} apart converge, as \eqn{t \to \infty}, on
#' \eqn{2(\sigma^2_b - \sigma^2_{ab}/(\sigma^2_a + \sigma^2_\epsilon))} --- the
#' same expression with the measurement error added to the baseline variance,
#' which is strictly larger than the floor whenever
#' \eqn{\sigma_{ab} \neq 0}. Only repeated measurement recovers the full
#' baseline correction, so reaching the floor takes a schedule that is both
#' long and dense.
#'
#' @param params A `slope_params` object.
#'
#' @return A single positive number: the infimum of \eqn{s^{*2}}.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_var_floor(pars)
#'
#' # Every schedule is above it, and a long dense one gets close.
#' slope_var(pars, c(0, 1, 2))
#' slope_var(pars, seq(0, 50, length.out = 501))
#'
#' @seealso [slope_var()], the same quantity at a stated schedule;
#'   [slope_sample_size_floor()], the sample size this implies.
#' @inherit stage_two references
#' @export
slope_var_floor <- function(params) {
  context <- "slope_var_floor()"
  check_params(params, context)
  var_floor(params, context)
}

# ---------------------------------------------------------------------------
# the sample-size floor and the power ceiling
# ---------------------------------------------------------------------------

#' Assemble a floor result from validated inputs
#'
#' The `effect_size` is the no-dropout `effect_components()` formula
#' -- `sign(d) * sqrt((d / sqrt(var))^2)`, i.e. `d / sqrt(var)` -- evaluated at
#' the limiting variance. Exactly one of `n` and `power` is supplied, as for
#' `solve_slope()`: the sample size comes from `size_per_arm()`, equation (6),
#' and the power from the same normal tail `solve_slope()` uses, both at the
#' limiting variance. Nothing here is a second implementation of anything.
#' @noRd
floor_result <- function(params, effectiveness, target, alpha, per_arm, context,
                         n = NULL, power = NULL) {
  check_params(params, context)
  solving_for_n <- check_n_or_power(alpha, n, power, context)
  per_arm <- check_per_arm(per_arm, context)

  comp <- target_components(params, target, effectiveness, context)
  var_tte <- var_floor(params, context)
  effect_size <- comp$slope_difference / sqrt(var_tte)
  scaled_effect <- scale_effect(effect_size, comp$effectiveness)
  z_a <- z_alpha(alpha, context)

  if (solving_for_n) {
    n_per_arm <- size_per_arm(scaled_effect, z_a, power)$n_per_arm
    cls <- "slope_sample_size_floor"
  } else {
    at_n <- power_at_n(scaled_effect, z_a, n)
    n_per_arm <- at_n$n_per_arm
    power <- at_n$power
    cls <- "slope_power_ceiling"
  }

  # No `design` argument, so the field is absent rather than NULL -- the
  # fields of CONTRACT.md section 4.3, from the same assembler that gives
  # solve_slope() its own.
  res <- stage_two_result(comp, n_per_arm = n_per_arm, power = power,
                          alpha = alpha, var_tte = var_tte,
                          effect_size = effect_size, params = params)
  if (!solving_for_n) {
    res <- add_n_requested(res, n)
  }
  structure(res, class = c(cls, "slope_result"), per_arm = per_arm)
}

#' Smallest sample size any trial design can need
#'
#' "Is this trial affordable at all?" --- the question worth asking before
#' searching over visit schedules. Because the treatment-effect variance
#' \eqn{s^{*2}} has a greatest lower bound that no schedule can beat (see
#' [slope_var_floor()]), so does the sample size. If the floor is already
#' unaffordable, no amount of design work will help and the target effect size
#' is what has to change.
#'
#' The bound is a property of the disease and the target effect, through
#' \eqn{\sigma^2_b} and \eqn{\sigma_{ab}/\sigma^2_a} --- not of the trial.
#'
#' @details
#' # No design argument
#'
#' Deliberately none. The bound holds for every visit
#' schedule, and for every dropout pattern too, since dropout only ever raises
#' the sample size. Passing a design would be passing something the answer does
#' not depend on --- the same reason [slope_effect_size()] refuses an
#' `effectiveness` argument.
#'
#' # How tight it is
#'
#' `slope_sample_size(params, visits, dropout, ...)$n` is greater than or equal
#' to this for every `visits` and `dropout`, strictly so before rounding. The bound is approached as
#' the schedule becomes long and dense, and after `ceiling()` a long enough
#' schedule reaches it exactly: on `slpower1` at the paper's 33% effectiveness
#' the floor is 236, and a fifty-year trial with visits every five weeks needs
#' 236. Ordinary schedules are nowhere near --- the paper's own two-year
#' three-visit design needs 712, three times the floor.
#'
#' `effectiveness` is not a detail here. Sample size scales as
#' `effectiveness^-2`, so a floor computed at the default 0.25 says nothing
#' about a trial powered for a 33% effect.
#'
#' @param params A `slope_params` object, from [slope_params()] or
#'   [slope_params_manual()].
#' @param effectiveness Proportion of the slope difference the treatment is
#'   expected to remove, in (0, 1]. Must not be supplied when
#'   `target = "observed"`, which fixes it at 1.
#' @param target `"effectiveness"` (the default) or `"observed"`. As in
#'   [slope_sample_size()]; see its "The reference slope" section.
#' @param power Desired power, between `alpha / 2` and 1. Defaults to 0.8.
#' @param alpha Two-sided significance level. Defaults to 0.05.
#' @param per_arm Which basis to print `N` on: `TRUE` (the default) for
#'   participants per arm, `FALSE` for the trial total. Display only, as in
#'   [slope_sample_size()]: the result carries both `n` and `n_per_arm`
#'   regardless, and `per_arm` is recorded as an attribute that `print()`
#'   consults and that can be overridden afterwards with
#'   `print(x, per_arm = FALSE)`.
#'
#' @return An object of class `slope_sample_size_floor`, a list with the same
#'   elements as a [slope_sample_size()] result **except `design`**, of which
#'   there is none: `n`, `n_per_arm`, `power`, `alpha`, `effectiveness`,
#'   `target`, `tte`, `var_tte`, `effect_size`, `slope_difference`,
#'   `reference_slope` and `params`. `var_tte` is [slope_var_floor()], the
#'   limiting \eqn{s^{*2}}. `per_arm` is recorded as an attribute rather than
#'   a list element, as for the other two entry points.
#'
#'   It inherits from `slope_result`, so [as.data.frame()][as.data.frame.slope_result]
#'   works and the row binds together with rows from the other two entry points.
#'   That row has `n_follow_up = NA` and `solve_for = "n_floor"`.
#'
#' @inheritSection stage_two The reference slope
#' @inherit stage_two references
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sample_size_floor(pars, effectiveness = 0.33)
#'
#' # How much of a design's sample size is the design's doing rather than the
#' # disease's: the same settings, with and without a schedule.
#' ss <- slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#' flr <- slope_sample_size_floor(pars, effectiveness = 0.33)
#' c(design = ss$n, floor = flr$n, ratio = ss$n / flr$n)
#'
#' @seealso [slope_var_floor()] for the variance behind it,
#'   [slope_sample_size()] for the sample size a stated design needs,
#'   [slope_sample_size_grid()] to search the designs that remain worth
#'   searching.
#' @export
slope_sample_size_floor <- function(params, power = 0.8, effectiveness = 0.25,
                                    target = c("effectiveness", "observed"),
                                    alpha = 0.05, per_arm = TRUE) {
  context <- "slope_sample_size_floor()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  floor_result(params, effectiveness, target, alpha, per_arm, context, power = power)
}

#' Highest power any trial design can achieve
#'
#' The power counterpart of [slope_sample_size_floor()]: at a fixed total
#' sample size `n`, the power of a trial whose treatment-effect variance is
#' the limiting [slope_var_floor()]. Because power rises as that variance
#' falls, no visit schedule and no dropout pattern can reach a higher power
#' with `n` participants: the variance's floor is the power's ceiling. If it is
#' already too low, more or better-placed visits will not rescue the trial: `n`
#' or the target effect has to change.
#'
#' Like [slope_sample_size_floor()] it is a bound, not an attainable value ---
#' `slope_power(params, visits, dropout, n = n, ...)$power` is strictly lower
#' for every finite schedule, approaching it only as the schedule becomes long
#' and dense. See [slope_sample_size_floor()] for how tight the bound is and
#' why it takes no visit schedule.
#'
#' @inheritParams slope_sample_size_floor
#' @param n Total number of participants across both arms, as in
#'   [slope_power()]. Required. Odd values are reduced by one so that the arms
#'   are equal, with the value as supplied kept in `n_requested`.
#' @param per_arm Which basis to print the counts on: `TRUE` (the default) for
#'   participants per arm, `FALSE` for the trial total. Display only, as in
#'   [slope_power()].
#'
#' @return An object of class `slope_power_ceiling`, a list with the same elements
#'   as a [slope_power()] result **except `design`**: `n`, `n_per_arm`,
#'   `n_requested`, `power`, `alpha`, `effectiveness`, `target`, `tte`,
#'   `var_tte`, `effect_size`, `slope_difference`, `reference_slope` and
#'   `params`. `power` is the upper bound and `var_tte` is [slope_var_floor()].
#'   It inherits from `slope_result`, so [as.data.frame()][as.data.frame.slope_result]
#'   gives a row with `n_follow_up = NA` and `solve_for = "power_ceiling"`.
#'
#' @inheritSection stage_two The reference slope
#' @inherit stage_two references
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_power_ceiling(pars, n = 200, effectiveness = 0.33)
#'
#' # How much of the achievable power a two-year, three-visit design delivers.
#' c(design  = slope_power(pars, c(0, 1, 2), n = 200, effectiveness = 0.33)$power,
#'   ceiling = slope_power_ceiling(pars, n = 200, effectiveness = 0.33)$power)
#'
#' @seealso [slope_sample_size_floor()], the sample-size bound from the same
#'   variance; [slope_power()] for the power of a stated design;
#'   [slope_var_floor()].
#' @export
slope_power_ceiling <- function(params, n, effectiveness = 0.25,
                              target = c("effectiveness", "observed"),
                              alpha = 0.05, per_arm = TRUE) {
  context <- "slope_power_ceiling()"
  require_n(missing(n) || is.null(n), paste0(
    "it is the sample size whose highest achievable\n",
    "  power is being bounded. For the smallest sample size any design could\n",
    "  need, use slope_sample_size_floor()."), context)
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  floor_result(params, effectiveness, target, alpha, per_arm, context, n = n)
}

#' Print a sample-size floor or power ceiling
#'
#' @param x A `slope_sample_size_floor` or `slope_power_ceiling` object.
#' @param ... Ignored.
#' @param per_arm Which basis to print `N` on: `TRUE` for participants per
#'   arm, `FALSE` for the trial total. Defaults to `NULL`, meaning "whatever
#'   the function that built `x` was called with" -- read from `x`'s
#'   `per_arm` attribute, or per arm if that is absent.
#' @return `x`, invisibly.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sample_size_floor(pars, effectiveness = 0.33)
#'
#' @export
print.slope_sample_size_floor <- function(x, ..., per_arm = NULL) {
  print_n_result(x, per_arm, "print.slope_sample_size_floor()")
}

#' @rdname print.slope_sample_size_floor
#' @export
print.slope_power_ceiling <- function(x, ..., per_arm = NULL) {
  print_power_result(x, per_arm, "print.slope_power_ceiling()")
}
