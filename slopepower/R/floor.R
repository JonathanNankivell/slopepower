# Layer 3 (continued) -- the limit of the layer-3 calculation over all designs.
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
#' since `sigma2_slope - sigma_cov^2 / sigma2_intercept` is its Schur
#' complement and positive definiteness makes the determinant
#' `sigma2_intercept * sigma2_slope - sigma_cov^2` positive. That is enforced
#' by `check_re_covariance()`, which every route into a `slope_params` object
#' runs and which `check_params()` re-runs on hand-built ones -- so there is no
#' zero or negative branch to guard here.
#' @noRd
var_floor <- function(params) {
  2 * (params$sigma2_slope - params$sigma_cov^2 / params$sigma2_intercept)
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
#' returned bounds every design [trial_design()] can express, not just the ones
#' without withdrawal.
#'
#' # It is an infimum, not a minimum
#'
#' No finite schedule attains it. Writing \eqn{\Sigma = \Sigma_0 + \sigma^2_\epsilon I},
#' every contrast \eqn{c} has \eqn{c^{\mathsf T}\Sigma c > c^{\mathsf T}\Sigma_0 c},
#' so `slope_var(params, visits)` is strictly greater than this value for any
#' `visits`, however long or however dense. The gap closes as the number of
#' visits grows without bound.
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
  var_floor(params)
}

# ---------------------------------------------------------------------------
# the sample-size floor
# ---------------------------------------------------------------------------

#' Assemble a floor result from validated inputs
#'
#' The `effect_size` is the no-dropout `effect_components()` formula
#' -- `sign(d) * sqrt((d / sqrt(var))^2)`, i.e. `d / sqrt(var)` -- evaluated at
#' the limiting variance, and the sample size comes from `size_per_arm()`, the
#' same equation (6) `solve_slope()` uses. Nothing here is a second
#' implementation of anything.
#' @noRd
floor_result <- function(params, effectiveness, target, power, alpha, context) {
  check_params(params, context)
  check_probability(alpha, "alpha", context)
  check_probability(power, "power", context)

  comp <- target_components(params, target, effectiveness, context)
  var_tte <- var_floor(params)
  effect_size <- comp$slope_difference / sqrt(var_tte)
  sized <- size_per_arm(abs(effect_size) * comp$effectiveness,
                        z_alpha(alpha, context), power)

  structure(
    list(
      n                = 2 * sized$n_per_arm,
      n_per_arm        = sized$n_per_arm,
      power            = power,
      alpha            = alpha,
      effectiveness    = if (identical(comp$target, "observed")) NA_real_ else comp$effectiveness,
      target           = comp$target,
      tte              = comp$tte,
      var_tte          = var_tte,
      effect_size      = effect_size,
      slope_difference = comp$slope_difference,
      reference_slope  = comp$reference_slope,
      params           = params
    ),
    class = c("slope_sample_size_floor", "slope_result")
  )
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
#' Deliberately none, in either method. The bound holds for every visit
#' schedule, and for every dropout pattern too, since dropout only ever raises
#' the sample size. Passing a design would be passing something the answer does
#' not depend on --- the same reason [slope_effect_size()] refuses an
#' `effectiveness` argument.
#'
#' # How tight it is
#'
#' `slope_sample_size(params, design, ...)$n` is greater than or equal to this
#' for every `design`, strictly so before rounding. The bound is approached as
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
#' @param x A `slope_params` object, or a result from [slope_sample_size()] or
#'   [slope_power()] whose settings should be reused.
#' @param effectiveness Proportion of the slope difference the treatment is
#'   expected to remove, in (0, 1]. Must not be supplied when
#'   `target = "observed"`, which fixes it at 1.
#' @param target `"effectiveness"` (the default) or `"observed"`. As in
#'   [slope_sample_size()]; see its "The reference slope" section.
#' @param power Desired power, in (0, 1). Defaults to 0.8.
#' @param alpha Two-sided significance level. Defaults to 0.05.
#' @param ... Not used; passing anything here is an error rather than being
#'   silently ignored.
#'
#' @return An object of class `slope_sample_size_floor`, a list with the same
#'   elements as a [slope_sample_size()] result **except `design`**, of which
#'   there is none: `n`, `n_per_arm`, `power`, `alpha`, `effectiveness`,
#'   `target`, `tte`, `var_tte`, `effect_size`, `slope_difference`,
#'   `reference_slope` and `params`. `var_tte` is [slope_var_floor()], the
#'   limiting \eqn{s^{*2}}.
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
#' # Reuse a result's own effectiveness, power and alpha, and see how much
#' # of the sample size is the design's doing rather than the disease's.
#' ss <- slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#' flr <- slope_sample_size_floor(ss)
#' c(design = ss$n, floor = flr$n, ratio = ss$n / flr$n)
#'
#' @seealso [slope_var_floor()] for the variance behind it,
#'   [slope_sample_size()] for the sample size a stated design needs,
#'   [slope_sample_size_grid()] to search the designs that remain worth
#'   searching.
#' @export
slope_sample_size_floor <- function(x, ...) {
  UseMethod("slope_sample_size_floor")
}

#' @describeIn slope_sample_size_floor Compute the floor from fitted or
#'   supplied parameters.
#' @export
slope_sample_size_floor.slope_params <- function(x,
                                                 effectiveness = 0.25,
                                                 target = c("effectiveness", "observed"),
                                                 power = 0.8, alpha = 0.05, ...) {
  context <- "slope_sample_size_floor()"
  reject_dots(list(...), paste0(
    "The floor holds for every visit schedule, so there is no `design`, no\n",
    "  `dropout` and no `n` to supply -- that is what makes it a floor. For the\n",
    "  sample size a particular design needs, use slope_sample_size()."),
    context)
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  floor_result(x, effectiveness, target, power, alpha, context)
}

#' @describeIn slope_sample_size_floor Compute the floor for a result already
#'   in hand, reusing its `effectiveness`, `target`, `power` and `alpha` so
#'   that the two numbers are comparable. For a [slope_power()] result the
#'   power reused is the one that design achieves, so the answer is the
#'   smallest sample size that could reach the same power.
#' @export
slope_sample_size_floor.slope_result <- function(x, ...) {
  context <- "slope_sample_size_floor()"
  reject_dots(list(...), paste0(
    "The calculation comes from the object, so `effectiveness`, `target`,\n",
    "  `power` and `alpha` belong to the call that produced it. To vary them,\n",
    "  take the floor of the parameters instead:\n",
    "    slope_sample_size_floor(x$params, power = 0.9)"),
    context)
  # A slope_power result whose design is comfortably over-powered reports
  # power == 1 in double precision; qnorm(1) is Inf and the floor would come
  # back as Inf, which reads as "no sample size is enough" -- the opposite of
  # the truth. check_probability()'s own message would name `power` as if the
  # caller had typed it, so the diagnosis is given here instead.
  if (isTRUE(x$power >= 1)) {
    stop(sprintf(paste0(
      "%s: this result's power is 1 to within double precision, so there is no\n",
      "  finite sample size that reaches it and the floor is undefined. Take the\n",
      "  floor at a stated power instead:\n",
      "    slope_sample_size_floor(x$params, power = 0.9)"), context), call. = FALSE)
  }
  # NA_real_ is what a target = "observed" result stores for `effectiveness`
  # (CONTRACT.md section 4.1); target_components() sets it to 1 on that branch
  # anyway, but passing the NA through would trip check_scalar() first if the
  # branch were ever reordered.
  effectiveness <- if (identical(x$target, "observed")) 1 else x$effectiveness
  floor_result(x$params, effectiveness, x$target, x$power, x$alpha, context)
}

#' @describeIn slope_sample_size_floor Reject anything else, with a pointer to
#'   what is accepted.
#' @export
slope_sample_size_floor.default <- function(x, ...) {
  stop(sprintf(paste0(
    "slope_sample_size_floor(): cannot compute a floor from an object of class %s.\n",
    "  The bound depends only on the variance components, so pass the parameters --\n",
    "  a slope_params object from slope_params() or slope_params_manual() -- or a\n",
    "  slope_sample_size or slope_power result to reuse the settings of."),
    paste(sQuote(class(x)), collapse = "/")), call. = FALSE)
}

#' Print a sample-size floor
#'
#' @param x A `slope_sample_size_floor` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sample_size_floor(pars, effectiveness = 0.33)
#'
#' @export
print.slope_sample_size_floor <- function(x, ...) {
  print_data_block(x)
  cat("\nParameters for planned study:\n")
  cat(fmt_line("alpha", x$alpha), "\n", sep = "")
  cat(fmt_line("power", x$power), "\n", sep = "")
  cat(fmt_line("effectiveness", x$effectiveness), "\n", sep = "")
  cat(fmt_line("target treatment difference in slopes", x$tte), "\n", sep = "")
  # Where print.slope_sample_size() shows the schedule and its dropouts. Saying
  # so explicitly, rather than omitting the line, is the point of the object:
  # the reader should not have to wonder which schedule produced the number.
  cat(fmt_line("visit schedule", "any (the bound holds for all)"), "\n", sep = "")
  cat("\n  Lower bound on sample size:\n")
  cat(fmt_line("N", x$n, digits = 0L), "\n", sep = "")
  cat(fmt_line("N per arm", x$n_per_arm, digits = 0L), "\n", sep = "")
  cat(fmt_line("limiting s*^2", x$var_tte), "\n", sep = "")
  cat("\n")
  invisible(x)
}
