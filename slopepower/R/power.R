# Layer 3 -- the mathematical core.
#
# Implements sections 5.1-5.6 of CONTRACT.md: the marginal covariance of the
# outcome at a set of visit times, the treatment-effect variance for a notional
# two-person trial, the Dawson-Lagakos pattern-mixture effect size, and the two
# closed-form solutions for sample size and power.

# ---------------------------------------------------------------------------
# internal validation
# ---------------------------------------------------------------------------

#' Required fields of a `slope_params` object
#' @noRd
PARAM_FIELDS <- c("slope", "slope_comparator", "comparator",
                  "sigma2_intercept", "sigma2_slope", "sigma_cov",
                  "sigma2_residual")

#' Validate a `slope_params` object well enough to compute with it
#' @noRd
check_params <- function(params, context) {
  if (!inherits(params, "slope_params")) {
    stop(sprintf("%s: `params` must be a `slope_params` object, as returned by ",
                 context),
         "`slope_params()` or `slope_params_manual()`.", call. = FALSE)
  }
  missing_fields <- setdiff(PARAM_FIELDS, names(params))
  if (length(missing_fields) > 0L) {
    stop(sprintf("%s: `params` is missing the field(s) %s.",
                 context, paste0("`", missing_fields, "`", collapse = ", ")),
         call. = FALSE)
  }
  check_scalar(params$slope, "params$slope", context)
  check_variance(params$sigma2_intercept, "params$sigma2_intercept", context)
  check_variance(params$sigma2_slope, "params$sigma2_slope", context)
  check_variance(params$sigma2_residual, "params$sigma2_residual", context)
  check_scalar(params$sigma_cov, "params$sigma_cov", context)

  if (!is.character(params$comparator) || length(params$comparator) != 1L ||
      !params$comparator %in% c("none", "healthy", "treated")) {
    stop(sprintf('%s: `params$comparator` must be one of "none", "healthy", "treated".',
                 context), call. = FALSE)
  }
  if (params$comparator != "none") {
    if (!is.numeric(params$slope_comparator) ||
        length(params$slope_comparator) != 1L ||
        !is.finite(params$slope_comparator)) {
      stop(sprintf('%s: `params$comparator` is "%s" but `params$slope_comparator` is not a finite number.',
                   context, params$comparator), call. = FALSE)
    }
  }
  invisible(params)
}

#' Validate a vector of visit times for the covariance construction
#'
#' Deliberately permissive: any strictly increasing finite vector of length >= 2
#' defines a valid covariance matrix. The stricter trial constraints (baseline at
#' zero) belong to `trial_design()` and are enforced in `slope_power()`.
#' @noRd
check_visits <- function(visits, context) {
  if (!is.numeric(visits) || length(visits) < 2L || any(!is.finite(visits))) {
    stop(sprintf("%s: `visits` must be a numeric vector of at least two finite times.",
                 context), call. = FALSE)
  }
  if (any(diff(visits) <= 0)) {
    stop(sprintf("%s: `visits` must be strictly increasing; got %s.",
                 context, paste(format(visits), collapse = ", ")), call. = FALSE)
  }
  invisible(as.numeric(visits))
}

#' Coerce and validate the `design` argument
#'
#' Accepts a `trial_design` object, or a bare numeric vector of visit times which
#' is passed to `trial_design()`.
#' @noRd
as_trial_design <- function(design, context) {
  if (is.numeric(design)) design <- trial_design(visits = design)
  if (!inherits(design, "trial_design")) {
    stop(sprintf("%s: `design` must be a `trial_design` object or a numeric vector of visit times.",
                 context), call. = FALSE)
  }
  # The rules for `visits` live in one place, `validate_visits()` in design.R, so
  # that a hand-built `trial_design` is held to exactly the same standard as one
  # from the constructor. The dropout checks below are invariant assertions only:
  # normalisation and the baseline-only warning belong to `trial_design()` and
  # must not be repeated here, or a constructed design would warn twice.
  visits <- validate_visits(design$visits, context)
  dropout <- design$dropout
  if (!is.numeric(dropout) || length(dropout) != length(visits) - 1L) {
    stop(sprintf("%s: `design$dropout` must have one entry per follow-up visit (%d), got %d.",
                 context, length(visits) - 1L, length(dropout)), call. = FALSE)
  }
  if (any(!is.finite(dropout)) || any(dropout < 0)) {
    stop(sprintf("%s: `design$dropout` entries must be finite and non-negative.",
                 context), call. = FALSE)
  }
  if (sum(dropout) > 1 + DROPOUT_TOL) {
    stop(sprintf("%s: dropout proportions sum to %g, which exceeds 1.",
                 context, sum(dropout)), call. = FALSE)
  }
  # `has_dropout` is what the calculation branches on, so it is derived here
  # rather than trusted. A hand-built design that omits it would otherwise fail
  # with "argument is of length zero", and one that sets it FALSE alongside a
  # non-zero `dropout` would silently report the unweighted s*^2.
  if (!identical(design$dropout_type, "incremental") &&
      !identical(design$dropout_type, "cumulative")) {
    stop(sprintf("%s: `design$dropout_type` must be \"incremental\" or \"cumulative\".",
                 context), call. = FALSE)
  }
  # `dropout_type` records only how the user supplied the vector; `dropout`
  # itself is always incremental, converted once by the constructor. So a value
  # of "cumulative" is provenance, not an instruction, and must not be acted on
  # here -- doing so rejected every design built with dropout_type =
  # "cumulative", including the one in trial_design()'s own examples, with an
  # error telling the user to do what they had already done. A hand-built list
  # that puts cumulative values in `dropout` cannot be detected anyway: the
  # values are indistinguishable from valid incremental ones.
  design$has_dropout <- any(dropout > 0)
  design$dropout <- dropout
  design$visits <- visits
  invisible(design)
}

# ---------------------------------------------------------------------------
# covariance and treatment-effect variance
# ---------------------------------------------------------------------------

#' Marginal covariance of the outcome at a set of visit times
#'
#' Builds the variance--covariance matrix of a single participant's repeated
#' measurements under the random intercept and slope model of Nash et al. (2021),
#' equation (1). This is the matrix printed in closed form on page 579 of the
#' paper, generalised to arbitrary real-valued visit times:
#'
#' \deqn{\Sigma_{ij} = \sigma^2_a + t_i t_j \sigma^2_b + (t_i + t_j)\sigma_{ab}
#'                     + [i = j]\sigma^2_\epsilon}
#'
#' The Stata implementation builds this on the unit-integer grid
#' \code{0:max(schedule)} and then selects the scheduled rows with a selection
#' matrix; the two agree exactly at integer times, and this form additionally
#' supports fractional visit spacing without the \code{scale()} option.
#'
#' @param params A `slope_params` object.
#' @param visits Numeric vector of visit times, strictly increasing, expressed in
#'   the same units as the `time` variable used to estimate `params`.
#'
#' @return A symmetric `length(visits)` by `length(visits)` matrix.
#' @export
slope_sigma <- function(params, visits) {
  context <- "slope_sigma()"
  check_params(params, context)
  t <- check_visits(visits, context)

  sigma <- params$sigma2_intercept +
    outer(t, t) * params$sigma2_slope +
    outer(t, t, "+") * params$sigma_cov +
    diag(params$sigma2_residual, length(t))

  if (!is_positive_definite(sigma)) {
    stop(sprintf("%s: the implied covariance matrix is not positive definite. ",
                 context),
         "Check the variance components in `params`.", call. = FALSE)
  }
  dimnames(sigma) <- list(format(t), format(t))
  sigma
}

#' Treatment-effect variance for a notional two-person trial
#'
#' Computes \eqn{s^{*2}}, the variance of the estimated difference in slopes in a
#' trial with one participant in each arm, following equation (5) of Nash et al.
#' (2021). Because the standard error of the treatment effect in a trial with
#' \eqn{N} participants per arm is \eqn{s^*/\sqrt{N}}, this single quantity
#' carries all of the design information needed for the sample-size formula.
#'
#' @param params A `slope_params` object.
#' @param visits Numeric vector of visit times.
#'
#' @return A single positive number, \eqn{s^{*2}}.
#' @export
slope_var <- function(params, visits) {
  context <- "slope_var()"
  check_params(params, context)
  t <- check_visits(visits, context)

  sigma <- slope_sigma(params, t)
  m <- length(t)
  zero <- matrix(0, m, m)

  # one person per arm: group indicator 0 for the first, 1 for the second
  sigma_star <- rbind(cbind(sigma, zero), cbind(zero, sigma))
  x_star <- rbind(cbind(1, t, 0), cbind(1, t, t))

  info <- t(x_star) %*% solve(sigma_star) %*% x_star
  f_star <- tryCatch(solve(info), error = function(e) {
    stop(sprintf("%s: the information matrix is singular at visit times %s. ",
                 context, paste(format(t), collapse = ", ")),
         "At least two distinct follow-up times are needed to identify a slope difference.",
         call. = FALSE)
  })
  f_star[3L, 3L]
}

# ---------------------------------------------------------------------------
# effect size
# ---------------------------------------------------------------------------

#' Resolve the reference slope, effectiveness and dropout-weighted effect size
#' @noRd
effect_components <- function(params, design, target, effectiveness,
                              effectiveness_supplied, context) {
  check_params(params, context)
  design <- as_trial_design(design, context)
  target <- match.arg(target, c("effectiveness", "observed"))

  comparator <- params$comparator

  if (identical(target, "observed")) {
    if (!identical(comparator, "treated")) {
      stop(sprintf('%s: target = "observed" reuses the treatment effect from a previous trial ',
                   context),
           sprintf('and requires `params` with comparator = "treated"; got "%s".',
                   comparator), call. = FALSE)
    }
    if (isTRUE(effectiveness_supplied)) {
      stop(sprintf('%s: supply only one of `effectiveness` and target = "observed". ',
                   context),
           'target = "observed" targets the previously observed effect in full.',
           call. = FALSE)
    }
    reference_slope <- params$slope_comparator
    effectiveness <- 1
  } else {
    check_scalar(effectiveness, "effectiveness", context,
                 lower = 0, upper = 1, lower_open = TRUE, upper_open = FALSE)
    reference_slope <- switch(comparator,
      none    = 0,
      healthy = params$slope_comparator,
      treated = 0
    )
  }

  slope_difference <- params$slope - reference_slope
  if (slope_difference == 0) {
    stop(sprintf("%s: the slope difference is exactly zero, so the target treatment ",
                 context),
         "effect is zero and no finite sample size gives power against it.",
         call. = FALSE)
  }
  tte <- -effectiveness * slope_difference

  visits <- design$visits
  dropout <- design$dropout
  n_follow_up <- length(visits) - 1L

  var_full <- slope_var(params, visits)
  es_full <- slope_difference / sqrt(var_full)

  # Dawson-Lagakos pattern mixture. Stratum j attends visits[1:j]; stratum 1 sees
  # baseline only and carries no slope information, so it is skipped -- an
  # infinite stratum-specific sample size contributes nothing to the sum.
  eff2 <- (1 - sum(dropout)) * es_full^2
  for (j in seq_len(n_follow_up)[-1L]) {
    if (dropout[j] == 0) next
    var_j <- slope_var(params, visits[seq_len(j)])
    eff2 <- eff2 + dropout[j] * (slope_difference / sqrt(var_j))^2
  }

  if (eff2 <= 0) {
    stop(sprintf("%s: every participant is expected to drop out before contributing ",
                 context),
         "slope information, so the effect size is zero. Check `design$dropout`.",
         call. = FALSE)
  }
  effect_size <- sign(slope_difference) * sqrt(eff2)

  # The target effect should move the untreated slope toward the reference, not
  # further from zero. It moves further out only when the reference is itself
  # more extreme than the group being treated.
  if (abs(params$slope + tte) > abs(params$slope)) {
    warning(sprintf(paste0("%s: the target treatment effect (%.4g) makes the slope more extreme ",
                           "(%.4g -> %.4g). The comparator slope is further from zero than the ",
                           "group being treated; check that `slope_comparator` is the intended target."),
                    context, tte, params$slope, params$slope + tte), call. = FALSE)
  }

  list(
    target           = target,
    effectiveness    = effectiveness,
    reference_slope  = reference_slope,
    slope_difference = slope_difference,
    tte              = tte,
    effect_size      = effect_size,
    var_full         = var_full,
    design           = design
  )
}

#' Dropout-weighted standardised effect size
#'
#' Combines the treatment-effect variance across dropout strata using the
#' pattern-mixture approach of Dawson and Lagakos (1991, 1993), as described in
#' section 2.5 of Nash et al. (2021). Each stratum's squared standardised effect
#' is weighted by the proportion of participants expected to follow that visit
#' pattern. Because sample size is inversely proportional to the squared effect
#' size, this is equivalent to taking the reciprocal of the weighted mean of the
#' reciprocals of the stratum-specific sample sizes.
#'
#' Participants who attend only the baseline visit contribute nothing, since a
#' single measurement carries no information about a slope.
#'
#' @param params A `slope_params` object.
#' @param design A `trial_design` object, or a numeric vector of visit times.
#' @param target `"effectiveness"` to measure the slope difference toward zero
#'   (or toward the healthy-control slope), or `"observed"` to measure it
#'   against the treated arm of a previous trial. See [slope_power()].
#'
#' @return A single signed number, on the **slope-difference scale**: it is
#'   \eqn{\Delta / s^*}, not \eqn{\beta_2 / s^*}. There is deliberately no
#'   `effectiveness` argument, because the returned value does not depend on
#'   one. Sample size is therefore
#'   \code{2 * ceiling((qnorm(1 - alpha / 2) + qnorm(power))^2 /
#'   (effect_size * effectiveness)^2)} --- omitting the `effectiveness` factor
#'   understates it by \eqn{1/e^2}. Use [slope_power()] unless you specifically
#'   want the unscaled quantity. The sign follows the slope difference; only the
#'   magnitude affects the sample size.
#' @export
slope_effect_size <- function(params, design,
                              target = c("effectiveness", "observed")) {
  # No `effectiveness` argument, deliberately. The returned effect size is on
  # the slope-difference scale (CONTRACT.md section 5.4) and does not depend on
  # effectiveness, so accepting the argument would silently ignore it: a caller
  # reconstructing N = 2 * ceiling((z + z)^2 / es^2) from this value would be
  # wrong by a factor of effectiveness^-2. See the note in @return.
  comp <- effect_components(params, design, target, effectiveness = 1,
                            effectiveness_supplied = FALSE,
                            context = "slope_effect_size()")
  comp$effect_size
}

# ---------------------------------------------------------------------------
# the main entry point
# ---------------------------------------------------------------------------

#' Sample size or power for a trial comparing slopes
#'
#' The second stage of the two-stage approach of Nash et al. (2021). Combines
#' slope and variance estimates -- from [slope_params()] fitted to previously
#' collected longitudinal data, or from [slope_params_manual()] -- with the
#' design of a proposed two-arm parallel trial, and returns either the required
#' sample size or the power of a given sample size.
#'
#' Supply exactly one of `n` and `power` and leave the other `NULL`; the missing
#' one is solved for. Supplying neither is equivalent to `power = 0.8`.
#'
#' The target treatment effect is defined relative to a reference slope, which
#' depends on the data the parameters came from:
#'
#' \describe{
#'   \item{`comparator = "none"`}{A treatment that is 100% effective halts change
#'     entirely; the reference slope is zero.}
#'   \item{`comparator = "healthy"`}{A treatment that is 100% effective reduces
#'     change to the rate seen in a disease-free population; the reference slope
#'     is that of the healthy controls.}
#'   \item{`comparator = "treated"` with `target = "observed"`}{The planned trial
#'     targets the same treatment effect as the previous trial.}
#'   \item{`comparator = "treated"` with `target = "effectiveness"`}{The previous
#'     trial's treated arm is ignored and the reference slope is zero, as in the
#'     Stata command's default behaviour for trial data.}
#' }
#'
#' To power a trial for an effect that is a fraction `p` of a previously observed
#' one, obtain the sample size with `target = "observed"` and multiply it by
#' `p^-2` (Nash et al. 2021, section 4.1.3).
#'
#' @param params A `slope_params` object.
#' @param design A `trial_design` object, or a numeric vector of visit times
#'   beginning at 0.
#' @param effectiveness Proportion of the slope difference the treatment is
#'   expected to remove, in (0, 1]. Ignored, and must not be supplied, when
#'   `target = "observed"`.
#' @param target `"effectiveness"` (the default) or `"observed"`. The latter is
#'   the equivalent of the Stata command's `usetrt` option and requires
#'   parameters estimated from a previous trial.
#' @param n Total number of participants across both arms. Odd values are reduced
#'   by one so that the arms are equal. Supply to compute power.
#' @param power Desired power. Supply to compute sample size.
#' @param alpha Two-sided significance level. Defaults to 0.05.
#'
#' @return An object of class `slope_power`. See [as.data.frame.slope_power()]
#'   for a tabular form with column names that are stable across designs.
#'
#' @references
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#'
#' @examples
#' pars <- slope_params_manual(
#'   slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
#'   sigma_cov = 5, sigma2_residual = 10
#' )
#' slope_power(pars, c(0, 1, 2), effectiveness = 0.33)
#' slope_power(pars, c(0, 1, 2), effectiveness = 0.33, n = 712)
#' @export
slope_power <- function(params, design,
                        effectiveness = 0.25,
                        target = c("effectiveness", "observed"),
                        n = NULL, power = NULL, alpha = 0.05) {
  context <- "slope_power()"

  if (!is.null(n) && !is.null(power)) {
    stop(sprintf("%s: supply only one of `n` and `power`; the other is solved for.",
                 context), call. = FALSE)
  }
  check_probability(alpha, "alpha", context)

  solve_for <- if (is.null(n)) "n" else "power"
  n_requested <- NA_real_

  if (identical(solve_for, "n")) {
    if (is.null(power)) power <- 0.8
    check_probability(power, "power", context)
  } else {
    check_scalar(n, "n", context, lower = 2, upper = Inf, lower_open = FALSE)
    if (n != floor(n)) {
      stop(sprintf("%s: `n` must be a whole number of participants; got %g.",
                   context, n), call. = FALSE)
    }
    n_requested <- n
  }

  comp <- effect_components(params, design, target, effectiveness,
                            effectiveness_supplied = !missing(effectiveness),
                            context = context)

  z_a <- stats::qnorm(1 - alpha / 2)
  scaled_effect <- abs(comp$effect_size) * comp$effectiveness

  if (identical(solve_for, "n")) {
    n_per_arm <- ceiling((z_a + stats::qnorm(power))^2 / scaled_effect^2)
    n_total <- 2 * n_per_arm
  } else {
    n_total <- 2 * floor(n_requested / 2)
    n_per_arm <- n_total / 2
    power <- stats::pnorm(scaled_effect * sqrt(n_per_arm) - z_a)
  }

  # With dropout no single s*^2 applies across strata, so report the effective
  # value obtained by inverting the sample-size formula. The algebra cancels in
  # both directions, recovering the dropout-weighted variance exactly.
  var_tte <- if (comp$design$has_dropout) {
    n_per_arm * comp$tte^2 / (z_a + stats::qnorm(power))^2
  } else {
    comp$var_full
  }

  structure(
    list(
      n                = n_total,
      n_per_arm        = n_per_arm,
      n_requested      = n_requested,
      power            = power,
      alpha            = alpha,
      effectiveness    = if (identical(comp$target, "observed")) NA_real_ else comp$effectiveness,
      target           = comp$target,
      tte              = comp$tte,
      var_tte          = var_tte,
      effect_size      = comp$effect_size,
      slope_difference = comp$slope_difference,
      reference_slope  = comp$reference_slope,
      solve_for        = solve_for,
      params           = params,
      design           = comp$design
    ),
    class = "slope_power"
  )
}

# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

#' Compact numeric formatting for visit times
#' @noRd
fmt_compact <- function(x) format(x, trim = TRUE, scientific = FALSE, digits = 7L)

#' Render the follow-up schedule the way the Stata command does
#' @noRd
schedule_string <- function(design) {
  follow_up <- design$visits[-1L]
  if (!design$has_dropout) return(paste(fmt_compact(follow_up), collapse = ", "))
  paste(sprintf("%s (%s)", fmt_compact(follow_up), fmt_compact(design$dropout)),
        collapse = ", ")
}

#' Print a slope sample-size or power calculation
#'
#' The layout follows the Stata command's output closely, which makes results
#' directly comparable with the worked examples in Nash et al. (2021).
#'
#' @param x A `slope_power` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.slope_power <- function(x, ...) {
  params <- x$params
  design <- x$design
  comparator <- params$comparator

  cat("\nData characteristics:\n")
  cat(fmt_line("number of observations in model", params$n_obs %||% NA, digits = 0L), "\n", sep = "")
  cat(fmt_line("number of participants in model", params$n_subjects %||% NA, digits = 0L), "\n", sep = "")

  show_difference <- comparator == "healthy" ||
    (comparator == "treated" && identical(x$target, "observed"))
  if (show_difference) {
    cat(fmt_line("observed difference in slopes", x$slope_difference), "\n", sep = "")
  }
  if (comparator == "treated") {
    cat(fmt_line("slope of control arm", params$slope), "\n", sep = "")
    if (identical(x$target, "observed")) {
      cat(fmt_line("slope of experimental arm", params$slope_comparator), "\n", sep = "")
    }
  } else {
    cat(fmt_line("slope of cases", params$slope), "\n", sep = "")
    if (comparator == "healthy") {
      cat(fmt_line("slope of healthy controls", params$slope_comparator), "\n", sep = "")
    }
  }

  cat("\nParameters for planned study:\n")
  cat(fmt_line("alpha", x$alpha), "\n", sep = "")
  if (identical(x$solve_for, "power")) {
    cat(fmt_line("specified N", x$n_requested, digits = 0L), "\n", sep = "")
    cat(fmt_line("actual N", x$n, digits = 0L), "\n", sep = "")
    cat(fmt_line("N per arm", x$n_per_arm, digits = 0L), "\n", sep = "")
  } else {
    cat(fmt_line("power", x$power), "\n", sep = "")
  }
  cat(fmt_line("effectiveness", x$effectiveness), "\n", sep = "")
  cat(fmt_line("target treatment difference in slopes", x$tte), "\n", sep = "")
  cat(fmt_line("number of follow-up visits", length(design$visits) - 1L, digits = 0L), "\n", sep = "")
  cat(fmt_line("schedule (and dropouts)", schedule_string(design)), "\n", sep = "")

  if (identical(x$solve_for, "n")) {
    cat("\n  Estimated sample size:\n")
    cat(fmt_line("N", x$n, digits = 0L), "\n", sep = "")
    cat(fmt_line("N per arm", x$n_per_arm, digits = 0L), "\n", sep = "")
  } else {
    cat("\nEstimated power:\n")
    cat(fmt_line("power", x$power), "\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' Coerce a slope power calculation to a one-row data frame
#'
#' Column names are identical for every combination of data type, comparator and
#' target, which makes results from different designs safe to bind together. This
#' differs deliberately from the Stata command, whose returned matrix renames its
#' slope columns depending on which model was fitted.
#'
#' @param x A `slope_power` object.
#' @param row.names,optional Passed on for consistency with the generic; ignored.
#' @param ... Ignored.
#' @return A one-row data frame.
#' @export
as.data.frame.slope_power <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    alpha            = x$alpha,
    power            = x$power,
    n                = x$n,
    n_per_arm        = x$n_per_arm,
    effectiveness    = x$effectiveness,
    target           = x$target,
    tte              = x$tte,
    var_tte          = x$var_tte,
    effect_size      = x$effect_size,
    slope            = x$params$slope,
    slope_comparator = as.numeric(x$params$slope_comparator %||% NA_real_),
    comparator       = x$params$comparator,
    reference_slope  = x$reference_slope,
    slope_difference = x$slope_difference,
    n_follow_up      = length(x$design$visits) - 1L,
    n_obs            = as.numeric(x$params$n_obs %||% NA_real_),
    n_subjects       = as.numeric(x$params$n_subjects %||% NA_real_),
    solve_for        = x$solve_for,
    row.names        = row.names,
    stringsAsFactors = FALSE
  )
}
