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
    stop(sprintf(paste0("%s: `params` must be a `slope_params` object, as returned by ",
                        "`slope_params()` or `slope_params_manual()`."), context), call. = FALSE)
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

  # new_slope_params() checks this at construction time for both routes into
  # the class, but a hand-built `slope_params` object bypasses it entirely --
  # and the marginal covariance built in sigma_at() can stay positive definite
  # even when the random-effects covariance itself is not, because a large
  # enough sigma2_residual masks it. So it is re-checked here, the one gate
  # every stage-two calculation funnels through, via the same helper
  # new_slope_params() uses -- see check_re_covariance() in utils.R.
  check_re_covariance(params$sigma2_intercept, params$sigma2_slope, params$sigma_cov, context)

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
#' zero) belong to `trial_design()` and are enforced by the stage-two entry
#' points, `slope_sample_size()` and `slope_power()`.
#' @noRd
check_visits <- function(visits, context) {
  if (!is.numeric(visits) || length(visits) < 2L || any(!is.finite(visits))) {
    stop(sprintf("%s: `visits` must be a numeric vector of at least two finite times.",
                 context), call. = FALSE)
  }
  if (is.unsorted(visits, strictly = TRUE)) {
    stop(sprintf("%s: `visits` must be strictly increasing; got %s.",
                 context, label_numeric(visits)), call. = FALSE)
  }
  invisible(as.numeric(visits))
}

# `as_trial_design()`, which coerces and revalidates the `design` argument of
# every stage-two call, lives in design.R beside the class it validates.

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
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sigma(pars, c(0, 1, 2))
#'
#' @seealso [slope_var()], the two-person treatment-effect variance built from
#'   this matrix; [slope_se()], the standard error of the fitted slope itself.
#' @export
slope_sigma <- function(params, visits) {
  context <- "slope_sigma()"
  check_params(params, context)
  t <- check_visits(visits, context)

  sigma <- sigma_at(params, t, context)
  # Only here, not in sigma_at(): the labels are for a human reading the printed
  # matrix, and formatting them costs more than the mathematics does (fmt_num()
  # is one format() call per element, against an outer() and an eigendecom-
  # position for the matrix itself). The internal callers index the matrix
  # numerically and would pay that on every grid cell and bootstrap replicate
  # for labels nothing reads.
  nm <- fmt_num(t)
  dimnames(sigma) <- list(nm, nm)
  sigma
}

#' Build the covariance matrix, given already-validated arguments
#'
#' The construction half of [slope_sigma()], split from it for the same reason
#' [treatment_effect_var()] is split from [slope_var()]: the callers on the hot
#' path have validated `params` and `visits` already, and re-validating is not
#' free -- `check_visits()` and `validate_visits()` would both run over the same
#' vector in a single stage-two call, which is the "two standards" problem those
#' validators exist to avoid.
#' @noRd
sigma_at <- function(params, t, context) {
  sigma <- params$sigma2_intercept +
    outer(t, t) * params$sigma2_slope +
    outer(t, t, "+") * params$sigma_cov +
    diag(params$sigma2_residual, length(t))

  if (!is_positive_definite(sigma)) {
    stop(sprintf(paste0("%s: the implied covariance matrix is not positive definite. ",
                        "Check the variance components in `params`."), context), call. = FALSE)
  }
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
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_var(pars, c(0, 1, 2))
#'
#' # A third follow-up visit adds information about the slope, so the
#' # treatment-effect variance shrinks.
#' slope_var(pars, c(0, 1, 2, 3))
#'
#' @seealso [slope_sigma()], the marginal covariance matrix this is built
#'   from; [slope_se()], the standard error of the fitted slope itself.
#' @export
slope_var <- function(params, visits) {
  context <- "slope_var()"
  check_params(params, context)
  t <- check_visits(visits, context)
  treatment_effect_var(sigma_at(params, t, context), t, context)
}

#' The GLS half of `slope_var()`: everything after `slope_sigma()` has built
#' and validated the covariance matrix
#'
#' Split out so [effect_components()] can pass a leading submatrix of an
#' already-validated covariance matrix for each dropout stratum, rather than
#' rebuilding it from `outer()` and re-running `slope_sigma()`'s
#' eigendecomposition-based positive-definite check on it. That re-check can
#' only ever confirm what is already known: a leading principal submatrix of a
#' positive definite matrix is itself always positive definite (Sylvester's
#' criterion), since its leading principal minors are a subset of the full
#' matrix's, which are all positive.
#'
#' @param sigma A covariance matrix already validated positive definite, e.g.
#'   by [slope_sigma()] or by being a leading submatrix of one.
#' @param t The visit times `sigma` was built at, already validated.
#' @noRd
treatment_effect_var <- function(sigma, t, context) {
  # One person per arm: group indicator 0 for the first, 1 for the second, so
  # equation (5)'s design matrix is X* = rbind(cbind(1, t, 0), cbind(1, t, t))
  # and its covariance Sigma* = diag(sigma, sigma) -- the two participants are
  # independent and identically scheduled, so both blocks are the same `sigma`.
  #
  # Neither 2m x 2m matrix is formed. A block diagonal matrix inverts blockwise,
  # so Sigma*^-1 = diag(sigma^-1, sigma^-1), and X*' Sigma*^-1 X* then separates
  # into the sum of the two participants' own contributions. That leaves one
  # m x m solve where the literal transcription does one of twice the order, for
  # about an eighth of the work; at thirteen visits it is a third of the cost of
  # this function. The result is bit-identical to the literal form at the visit
  # counts in Nash et al. (2021), and agrees to within rounding beyond them --
  # the summation order in the matrix products differs, nothing else.
  si <- solve(sigma)
  x1 <- cbind(1, t, 0)
  x2 <- cbind(1, t, t)

  info <- crossprod(x1, si %*% x1) + crossprod(x2, si %*% x2)
  f_star <- tryCatch(solve(info), error = function(e) {
    stop(sprintf(paste0("%s: the information matrix is singular at visit times %s. ",
                        "At least two distinct follow-up times are needed to identify ",
                        "a slope difference."),
                 context, label_numeric(t)), call. = FALSE)
  })
  f_star[3L, 3L]
}

# ---------------------------------------------------------------------------
# effect size
# ---------------------------------------------------------------------------

#' Resolve the reference slope, the effectiveness and the target treatment
#' effect --- everything that does not involve the visit schedule
#'
#' Split from [effect_components()] because the design-free half is the whole of
#' what [slope_sample_size_floor()] needs: the floor is a bound over every
#' schedule, so it has no `design` to weight by, but it targets exactly the same
#' \eqn{\beta_2} and must reach it by exactly the same route. Duplicating this
#' resolution there would let the two drift into targeting different effects
#' from the same `params`.
#'
#' `params` arrives already validated by [check_params()] -- both callers are
#' entry points that run it themselves, and re-running it here would be the
#' "two standards" duplication that [sigma_at()] and [treatment_effect_var()]
#' are likewise split to avoid.
#'
#' `effectiveness` and `target` arrive already reconciled: whether the caller
#' typed an `effectiveness` is knowable only at the exported boundary, so
#' [check_target_effectiveness()] is applied there -- by the stage-two entry
#' points, the grid wrappers and the floor alike -- rather than being plumbed
#' down here as a `missing()` flag on every intervening signature.
#' @noRd
target_components <- function(params, target, effectiveness, context) {
  target <- match.arg(target, c("effectiveness", "observed"))

  comparator <- params$comparator

  if (identical(target, "observed")) {
    if (!identical(comparator, "treated")) {
      stop(sprintf(paste0('%s: target = "observed" reuses the treatment effect from a previous ',
                          'trial and requires `params` with comparator = "treated"; got "%s".'),
                   context, comparator), call. = FALSE)
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
    stop(sprintf(paste0("%s: the slope difference is exactly zero, so the target treatment ",
                        "effect is zero and no finite sample size gives power against it."),
                 context), call. = FALSE)
  }
  tte <- -effectiveness * slope_difference

  # The target effect should move the untreated slope toward the reference, not
  # further from zero. It moves further out only when the reference is itself
  # more extreme than the group being treated.
  if (abs(params$slope + tte) > abs(params$slope)) {
    # Classed like the baseline-dropout warning in design.R, and for the same
    # reason: grid_impl() collects this one specifically, by class, to report
    # it once per grid rather than once per cell.
    warning(warningCondition(
      sprintf(paste0("%s: the target treatment effect (%.4g) makes the slope more extreme ",
                     "(%.4g -> %.4g). The comparator slope is further from zero than the ",
                     "group being treated; check that `slope_comparator` is the intended target."),
              context, tte, params$slope, params$slope + tte),
      class = "slopepower_tte_direction"))
  }

  list(
    target           = target,
    effectiveness    = effectiveness,
    reference_slope  = reference_slope,
    slope_difference = slope_difference,
    tte              = tte
  )
}

#' Resolve the reference slope, effectiveness and dropout-weighted effect size
#'
#' The design-dependent half; [target_components()] is the rest, and is called
#' *after* `as_trial_design()` so that a design warning still precedes a target
#' warning in the order they always did.
#' @noRd
effect_components <- function(params, design, target, effectiveness, context) {
  check_params(params, context)
  design <- as_trial_design(design, context)
  comp <- target_components(params, target, effectiveness, context)
  slope_difference <- comp$slope_difference

  visits <- design$visits
  dropout <- design$dropout
  n_follow_up <- length(visits) - 1L

  sigma_full <- sigma_at(params, visits, context)
  var_full <- treatment_effect_var(sigma_full, visits, context)
  es_full <- slope_difference / sqrt(var_full)

  # Dawson-Lagakos pattern mixture. Stratum j attends visits[1:j]; stratum 1 sees
  # baseline only and carries no slope information, so it is skipped -- an
  # infinite stratum-specific sample size contributes nothing to the sum.
  #
  # Stratum j's covariance is exactly the leading j x j submatrix of
  # `sigma_full` (Sigma_ik depends only on visits i and k, not on what other
  # visits exist), so it is sliced out rather than rebuilt and re-validated
  # from scratch -- see treatment_effect_var()'s note on why that re-check
  # would be redundant.
  eff2 <- (1 - sum(dropout)) * es_full^2
  for (j in seq_len(n_follow_up)[-1L]) {
    if (dropout[j] == 0) next
    idx <- seq_len(j)
    var_j <- treatment_effect_var(sigma_full[idx, idx, drop = FALSE], visits[idx], context)
    eff2 <- eff2 + dropout[j] * (slope_difference / sqrt(var_j))^2
  }

  if (eff2 <= 0) {
    stop(sprintf(paste0("%s: every participant is expected to drop out before contributing ",
                        "slope information, so the effect size is zero. Check `design$dropout`."),
                 context), call. = FALSE)
  }
  effect_size <- sign(slope_difference) * sqrt(eff2)

  c(comp, list(effect_size = effect_size, var_full = var_full, design = design))
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
#' single measurement carries no information about a slope. Every stratum uses
#' the same slope difference and the same variance components, so withdrawal is
#' assumed unrelated to a participant's own trajectory, and only monotone dropout
#' from a common schedule can be expressed. See [trial_design()] for what that
#' does and does not cover.
#'
#' @param params A `slope_params` object.
#' @param design A `trial_design` object, or a numeric vector of visit times.
#' @param target `"effectiveness"` to measure the slope difference toward zero
#'   (or toward the healthy-control slope), or `"observed"` to measure it
#'   against the treated arm of a previous trial. See [slope_sample_size()].
#'
#' @return A single signed number, on the **slope-difference scale**: it is
#'   \eqn{\Delta / s^*}, not \eqn{\beta_2 / s^*}. There is deliberately no
#'   `effectiveness` argument, because the returned value does not depend on
#'   one. Sample size is therefore
#'   \code{2 * ceiling((qnorm(1 - alpha / 2) + qnorm(power))^2 /
#'   (effect_size * effectiveness)^2)} --- omitting the `effectiveness` factor
#'   understates it by \eqn{1/e^2}. Use [slope_sample_size()] unless you
#'   specifically want the unscaled quantity. The sign follows the slope
#'   difference; only the magnitude affects the sample size.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_effect_size(pars, c(0, 1, 2))
#'
#' # A quarter of participants are expected to miss the final visit: the
#' # effect size shrinks because they contribute less slope information.
#' design <- trial_design(c(0, 1, 2), dropout = c(0, 0.25))
#' slope_effect_size(pars, design)
#'
#' @inherit stage_two references
#' @seealso [trial_design()] for how dropout is specified and what the
#'   pattern-mixture weighting assumes.
#' @export
slope_effect_size <- function(params, design,
                              target = c("effectiveness", "observed")) {
  # No `effectiveness` argument, deliberately. The returned effect size is on
  # the slope-difference scale (CONTRACT.md section 5.4) and does not depend on
  # effectiveness, so accepting the argument would silently ignore it: a caller
  # reconstructing N = 2 * ceiling((z + z)^2 / es^2) from this value would be
  # wrong by a factor of effectiveness^-2. See the note in @return.
  comp <- effect_components(params, design, target, effectiveness = 1,
                            context = "slope_effect_size()")
  comp$effect_size
}

# ---------------------------------------------------------------------------
# the shared solver
# ---------------------------------------------------------------------------

#' Equation (6): the per-arm sample size a scaled effect needs
#'
#' The one place `ceiling((z_{1-alpha/2} + z_{power})^2 / effect^2)` is written.
#' [solve_slope()] uses it for a stated design and
#' [slope_sample_size_floor()] for the limit over all of them; sharing it is
#' what makes "no design can beat this" a statement about the same arithmetic
#' rather than about a second implementation of it.
#'
#' `scaled_effect` is `abs(effect_size) * effectiveness` --- both factors, per
#' CONTRACT.md section 5.5.
#' @noRd
size_per_arm <- function(scaled_effect, z_a, power) {
  z_sum_sq <- (z_a + stats::qnorm(power))^2
  list(z_sum_sq = z_sum_sq, n_per_arm = ceiling(z_sum_sq / scaled_effect^2))
}

#' Everything `slope_sample_size()` and `slope_power()` have in common
#'
#' Exactly one of `n` and `power` is supplied; the other is solved for. This is
#' the shared machinery, not an interface: the two questions take different
#' inputs and produce differently shaped answers, so they are separate exported
#' functions rather than one function with a mode switch. Keeping the algebra
#' here means the sample size a design needs and the power it achieves can never
#' drift apart.
#'
#' Returns the fields common to both results, in their canonical order.
#' `n_requested` belongs only to the power branch and is added by its caller.
#' @noRd
solve_slope <- function(params, design, effectiveness,
                        target, alpha, n, power, context) {
  check_probability(alpha, "alpha", context)

  solving_for_n <- is.null(n)
  if (solving_for_n) {
    check_probability(power, "power", context)
  } else {
    check_whole_number(n, "n", "participants", context, lower = 2)
  }

  comp <- effect_components(params, design, target, effectiveness, context)

  z_a <- z_alpha(alpha, context)
  scaled_effect <- abs(comp$effect_size) * comp$effectiveness

  if (solving_for_n) {
    sized <- size_per_arm(scaled_effect, z_a, power)
    # (z_alpha + z_beta)^2, the numerator of the sample-size formula. Kept
    # because the var_tte back-solve below divides by this same quantity, and
    # the two cannot be allowed to drift apart.
    z_sum_sq <- sized$z_sum_sq
    n_per_arm <- sized$n_per_arm
    n_total <- 2 * n_per_arm
  } else {
    n_total <- 2 * floor(n / 2)
    n_per_arm <- n_total / 2
    power <- stats::pnorm(scaled_effect * sqrt(n_per_arm) - z_a)
  }

  # With dropout no single s*^2 applies across strata, so report the effective
  # value obtained by inverting the sample-size formula (CONTRACT.md 5.6).
  #
  # The two branches are the same algebra but not the same number. Solving for
  # power, qnorm inverts the pnorm above and n_per_arm cancels exactly, leaving
  # tte^2 / scaled_effect^2 -- so that form is used directly rather than round
  # -tripping through the power. The round trip is not merely redundant: past
  # n_per_arm of a few thousand the power saturates at exactly 1, qnorm(1) is
  # Inf, and the back-solve silently reported an effective variance of 0.
  # Solving for n there is no such cancellation, because n_per_arm has been
  # rounded up to a whole participant; the reported value carries that rounding
  # and is very slightly larger, which is what the Stata original reports too.
  var_tte <- if (!comp$design$has_dropout) {
    comp$var_full
  } else if (solving_for_n) {
    n_per_arm * comp$tte^2 / z_sum_sq
  } else {
    comp$tte^2 / scaled_effect^2
  }

  list(
    n                = n_total,
    n_per_arm        = n_per_arm,
    power            = power,
    alpha            = alpha,
    effectiveness    = if (identical(comp$target, "observed")) NA_real_ else comp$effectiveness,
    target           = comp$target,
    tte              = comp$tte,
    var_tte          = var_tte,
    effect_size      = comp$effect_size,
    slope_difference = comp$slope_difference,
    reference_slope  = comp$reference_slope,
    params           = params,
    design           = comp$design
  )
}

#' Shared documentation for the two stage-two entry points
#'
#' @param params A `slope_params` object, from [slope_params()] fitted to
#'   previously collected longitudinal data or from [slope_params_manual()].
#' @param design A `trial_design` object, from [trial_design()], or a numeric
#'   vector of visit times beginning at 0.
#' @param effectiveness Proportion of the slope difference the treatment is
#'   expected to remove, in (0, 1] (see "The reference slope" below). Must not
#'   be supplied when `target = "observed"`, which fixes it at 1.
#' @param target `"effectiveness"` (the default) or `"observed"`. The latter is
#'   the equivalent of the Stata command's `usetrt` option and requires
#'   parameters estimated from a previous trial.
#' @param alpha Two-sided significance level. Defaults to 0.05.
#'
#' @section The reference slope:
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
#' @section Dropout:
#'
#' When `design` carries a dropout pattern, the calculation is adjusted by the
#' pattern-mixture method of Dawson and Lagakos (1991, 1993), following section
#' 2.5 of Nash et al. (2021). Participants are stratified by the visits they
#' attend, each stratum is sized as though the whole trial followed that pattern,
#' and the strata are combined as the reciprocal of the weighted mean of the
#' reciprocals of those sizes. Withdrawers thus still contribute the visits they
#' attended, which is less conservative than dividing a completers-only sample
#' size by the completion rate, and is the right adjustment when the trial will
#' be analysed with a mixed model on all observed measurements. Only monotone
#' dropout from a schedule common to all participants is modelled, and every
#' stratum is assumed to share the same variance components, so withdrawal is
#' taken to be unrelated to a participant's own trajectory. See [trial_design()]
#' for the full account and for how the `dropout` vector itself is read.
#'
#' One field of the result changes meaning under dropout. No single \eqn{s^{*2}}
#' applies across strata, so `var_tte` reports the effective value obtained by
#' inverting the sample-size formula, rather than [slope_var()] at the full visit
#' schedule.
#'
#' @references
#' Dawson, J. D., and S. W. Lagakos. 1991. Analyzing laboratory marker changes in
#' AIDS clinical trials. \emph{Journal of Acquired Immune Deficiency Syndromes}
#' 4: 667--676.
#'
#' Dawson, J. D., and S. W. Lagakos. 1993. Size and power of two-sample tests of
#' repeated measures data. \emph{Biometrics} 49: 1022--1032.
#' \doi{10.2307/2532244}
#'
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#'
#' @name stage_two
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# sample size
# ---------------------------------------------------------------------------

#' Sample size for a trial comparing slopes
#'
#' "How many participants do I need?" --- the second stage of the two-stage
#' approach of Nash et al. (2021), solved for `n`. Combines slope and variance
#' estimates with the design of a proposed two-arm parallel trial and returns the
#' total sample size that achieves `power`.
#'
#' Use [slope_power()] for the converse question, "what power does this many
#' participants buy?". The two take different inputs and answer differently
#' shaped questions, so they are separate functions; they share their algebra
#' internally and are exact inverses up to the `ceiling()` that rounds a
#' fractional participant up.
#'
#' To power a trial for an effect that is a fraction `p` of a previously observed
#' one, obtain the sample size with `target = "observed"` and multiply it by
#' `p^-2` (Nash et al. 2021, section 4.1.3, whose p.594 worked example scales
#' N = 318 to 1,272 for `p = 0.5`).
#'
#' That shortcut is very slightly conservative. The scaling law is exact on the
#' unrounded per-arm requirement, but `n` has already been rounded up to a whole
#' participant per arm, so multiplying magnifies the rounding. Rounding last
#' instead gives 1,266 on that example, and agrees with what this function
#' returns for parameters whose reference slope has been halved directly:
#'
#' ```
#' ss <- slope_sample_size(params, design, target = "observed")
#' 2 * ceiling((qnorm(1 - ss$alpha / 2) + qnorm(ss$power))^2 /
#'             (p * ss$effect_size)^2)
#' ```
#'
#' The gap is bounded by `p^-2` participants per arm. Either is defensible ---
#' the published figure recruits a few people more than it strictly needs.
#'
#' @inheritParams stage_two
#' @param power Desired power, in (0, 1). Defaults to 0.8.
#'
#' @return An object of class `slope_sample_size`, a list with elements `n`,
#'   `n_per_arm`, `power`, `alpha`, `effectiveness`, `target`, `tte`, `var_tte`,
#'   `effect_size`, `slope_difference`, `reference_slope`, `params` and `design`.
#'   See [as.data.frame.slope_result()] for a tabular form whose column names are
#'   stable across designs and across both entry points.
#'
#' @inheritSection stage_two The reference slope
#' @inheritSection stage_two Dropout
#' @inherit stage_two references
#'
#' @examples
#' # No comparator: measured toward a slope of zero, fitted to all two
#' # hundred participants of `slpower1`. Reproduces the paper's p.588 result.
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#' slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33, power = 0.9)
#'
#' # Case/healthy-control comparator: measured toward the healthy controls'
#' # slope. Two cases and two controls, a subset of `slpower2`, whose visits
#' # are calendar dates and so are converted to years in the formula.
#' df2 <- slpower2[slpower2$id %in% c(1, 2, 251, 252), ]
#' pars2 <- slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, data = df2, healthy = case)
#' slope_sample_size(pars2, c(0, 1, 2), effectiveness = 0.33)
#'
#' # Randomised-trial comparator, target = "observed": size a repeat trial to
#' # detect the same effect the trial actually found, fitted to all one
#' # hundred and fifty participants of `slpower3`.
#' pars3 <- slope_params(sdmt ~ visit | id, data = slpower3, treated = treat)
#' slope_sample_size(pars3, c(0, 0.5, 2), target = "observed")
#'
#' @seealso [trial_design()] to build the `design` argument,
#'   [slope_power()] for the power of a given `n`,
#'   [slope_sample_size_grid()] to compare many designs at once,
#'   [slope_bootstrap()] for an interval around the result.
#' @export
slope_sample_size <- function(params, design,
                              effectiveness = 0.25,
                              target = c("effectiveness", "observed"),
                              power = 0.8, alpha = 0.05) {
  context <- "slope_sample_size()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  res <- solve_slope(params, design, effectiveness,
                     target = target, alpha = alpha,
                     n = NULL, power = power, context = context)
  structure(res, class = c("slope_sample_size", "slope_result"))
}

# ---------------------------------------------------------------------------
# power
# ---------------------------------------------------------------------------

#' Power of a trial comparing slopes
#'
#' "What power does this many participants buy?" --- the second stage of the
#' two-stage approach of Nash et al. (2021), solved for power. Combines slope and
#' variance estimates with the design of a proposed two-arm parallel trial and a
#' fixed total sample size.
#'
#' Use [slope_sample_size()] for the converse question, "how many participants do
#' I need?". The two take different inputs and answer differently shaped
#' questions, so they are separate functions; they share their algebra internally
#' and are exact inverses up to the `ceiling()` that rounds a fractional
#' participant up.
#'
#' @inheritParams stage_two
#' @param n Total number of participants across both arms. Required: it is the
#'   quantity whose power is being evaluated. Odd values are reduced by one so
#'   that the arms are equal, with the value as supplied kept in `n_requested`.
#'
#' @return An object of class `slope_power`, a list with elements `n`,
#'   `n_per_arm`, `n_requested`, `power`, `alpha`, `effectiveness`, `target`,
#'   `tte`, `var_tte`, `effect_size`, `slope_difference`, `reference_slope`,
#'   `params` and `design`. See [as.data.frame.slope_result()] for a tabular form
#'   whose column names are stable across designs and across both entry points.
#'
#' @inheritSection stage_two The reference slope
#' @inheritSection stage_two Dropout
#' @inherit stage_two references
#'
#' @examples
#' # No comparator: measured toward a slope of zero, fitted to all two
#' # hundred participants of `slpower1`. Reproduces the paper's p.588 result.
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_power(pars, c(0, 1, 2), n = 712, effectiveness = 0.33)
#'
#' # Case/healthy-control comparator: measured toward the healthy controls'
#' # slope. Two cases and two controls, a subset of `slpower2`, whose visits
#' # are calendar dates and so are converted to years in the formula.
#' df2 <- slpower2[slpower2$id %in% c(1, 2, 251, 252), ]
#' pars2 <- slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, data = df2, healthy = case)
#' slope_power(pars2, c(0, 1, 2), n = 40, effectiveness = 0.33)
#'
#' # Randomised-trial comparator, target = "observed": what power would a
#' # repeat trial have to detect the same effect the trial actually found,
#' # fitted to all one hundred and fifty participants of `slpower3`.
#' pars3 <- slope_params(sdmt ~ visit | id, data = slpower3, treated = treat)
#' slope_power(pars3, c(0, 0.5, 2), n = 396, target = "observed")
#'
#' @seealso [trial_design()] to build the `design` argument,
#'   [slope_sample_size()] for the `n` a target power needs,
#'   [slope_power_grid()] to compare many designs at once,
#'   [slope_bootstrap()] for an interval around the result.
#' @export
slope_power <- function(params, design, n,
                        effectiveness = 0.25,
                        target = c("effectiveness", "observed"),
                        alpha = 0.05) {
  context <- "slope_power()"
  # `is.null(n)` as well as `missing(n)`: solve_slope() picks its branch on
  # is.null(), so an explicit n = NULL -- the shape a programmatic caller gets
  # from do.call() with an unset element -- would otherwise slip past this guard
  # into the solve-for-n branch and fail complaining about `power`, an argument
  # this function does not have.
  if (missing(n) || is.null(n)) {
    stop(sprintf(paste0(
      "%s: `n` is required -- it is the sample size whose power is being\n",
      "  evaluated. To solve for the sample size that achieves a given power,\n",
      "  use slope_sample_size(params, design, power = 0.8)."),
      context), call. = FALSE)
  }
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  res <- solve_slope(params, design, effectiveness,
                     target = target, alpha = alpha,
                     n = n, power = NULL, context = context)
  # Kept separately from `n`, which is the even number actually used. Positioned
  # by name rather than by index: CONTRACT.md section 4.2 fixes it "after
  # n_per_arm", and solve_slope() assembles that list two hundred lines away.
  res <- append(res, list(n_requested = as.numeric(n)),
                after = match("n_per_arm", names(res)))
  structure(res, class = c("slope_power", "slope_result"))
}

# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

#' Render the follow-up schedule the way the Stata command does
#' @noRd
schedule_string <- function(design) {
  follow_up <- design$visits[-1L]
  if (!design$has_dropout) return(label_numeric(follow_up))
  paste(sprintf("%s (%s)", fmt_num(follow_up), fmt_num(design$dropout)),
        collapse = ", ")
}

#' The "Data characteristics" block, common to both print methods
#'
#' Layout follows the Stata command's output closely, which makes results
#' directly comparable with the worked examples in Nash et al. (2021).
#' @noRd
print_data_block <- function(x) {
  params <- x$params
  comparator <- params$comparator

  cat("\nData characteristics:\n")
  cat(fmt_line("number of observations in model", params$n_obs %||% NA, digits = 0L), "\n", sep = "")
  cat(fmt_line("number of participants in model", params$n_subjects %||% NA, digits = 0L), "\n", sep = "")

  # The comparator slope, and the difference from it, are shown exactly when the
  # target is measured against it: for healthy controls always, for a previous
  # trial's treated arm only under target = "observed" -- with
  # target = "effectiveness" that arm is ignored and the reference slope is zero
  # (CONTRACT.md section 5.3), so printing it would misdescribe the calculation.
  labels <- slope_labels(comparator)
  show_comparator <- comparator == "healthy" ||
    (comparator == "treated" && identical(x$target, "observed"))

  if (show_comparator) {
    cat(fmt_line("observed difference in slopes", x$slope_difference), "\n", sep = "")
  }
  cat(fmt_line(labels$own, params$slope), "\n", sep = "")
  if (show_comparator) {
    cat(fmt_line(labels$comparator, params$slope_comparator), "\n", sep = "")
  }
  invisible(x)
}

#' The tail of the "Parameters for planned study" block, common to both
#' @noRd
print_design_block <- function(x) {
  design <- x$design
  cat(fmt_line("effectiveness", x$effectiveness), "\n", sep = "")
  cat(fmt_line("target treatment difference in slopes", x$tte), "\n", sep = "")
  cat(fmt_line("number of follow-up visits", length(design$visits) - 1L, digits = 0L), "\n", sep = "")
  cat(fmt_line("schedule (and dropouts)", schedule_string(design)), "\n", sep = "")
  invisible(x)
}

#' Print a slope sample-size calculation
#'
#' @param x A `slope_sample_size` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#'
#' @export
print.slope_sample_size <- function(x, ...) {
  print_data_block(x)
  cat("\nParameters for planned study:\n")
  cat(fmt_line("alpha", x$alpha), "\n", sep = "")
  cat(fmt_line("power", x$power), "\n", sep = "")
  print_design_block(x)
  cat("\n  Estimated sample size:\n")
  cat(fmt_line("N", x$n, digits = 0L), "\n", sep = "")
  cat(fmt_line("N per arm", x$n_per_arm, digits = 0L), "\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Print a slope power calculation
#'
#' @param x A `slope_power` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_power(pars, c(0, 1, 2), n = 712, effectiveness = 0.33)
#'
#' @export
print.slope_power <- function(x, ...) {
  print_data_block(x)
  cat("\nParameters for planned study:\n")
  cat(fmt_line("alpha", x$alpha), "\n", sep = "")
  cat(fmt_line("specified N", x$n_requested, digits = 0L), "\n", sep = "")
  cat(fmt_line("actual N", x$n, digits = 0L), "\n", sep = "")
  cat(fmt_line("N per arm", x$n_per_arm, digits = 0L), "\n", sep = "")
  print_design_block(x)
  cat("\nEstimated power:\n")
  cat(fmt_line("power", x$power), "\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Coerce a stage-two result to a one-row data frame
#'
#' Column names are identical for every combination of data type, comparator and
#' target, and for both [slope_sample_size()] and [slope_power()], which makes
#' results from different designs and different questions safe to bind together.
#' This differs deliberately from the Stata command, whose returned matrix renames
#' its slope columns depending on which model was fitted.
#'
#' `solve_for` records which question produced the row --- `"n"` for
#' [slope_sample_size()], `"power"` for [slope_power()], `"n_floor"` for
#' [slope_sample_size_floor()]. It is derived from the object's class, and
#' exists so that a bound table stays interpretable; the functions themselves
#' have no such switch.
#'
#' A `"n_floor"` row has `n_follow_up = NA`, because the bound it reports holds
#' for every visit schedule and so is attached to none of them.
#'
#' @param x A `slope_sample_size`, `slope_power` or `slope_sample_size_floor`
#'   object.
#' @param row.names,optional Passed on for consistency with the generic; ignored.
#' @param ... Ignored.
#' @return A one-row data frame.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' res <- slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#' as.data.frame(res)
#'
#' # Stable column names make results from different designs and different
#' # questions safe to bind into one comparison table.
#' res2 <- slope_power(pars, c(0, 1, 2, 3), n = 712, effectiveness = 0.33)
#' rbind(as.data.frame(res), as.data.frame(res2))
#'
#' # The floor binds in too, as the row no schedule can beat.
#' rbind(as.data.frame(res),
#'       as.data.frame(slope_sample_size_floor(pars, effectiveness = 0.33)))
#'
#' @export
as.data.frame.slope_result <- function(x, row.names = NULL, optional = FALSE, ...) {
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
    # NA, not 0 or -1: a floor result carries no `design` at all, because the
    # bound holds whatever the schedule is (CONTRACT.md section 4.3).
    n_follow_up      = if (is.null(x$design)) NA_integer_ else length(x$design$visits) - 1L,
    n_obs            = as.numeric(x$params$n_obs %||% NA_real_),
    n_subjects       = as.numeric(x$params$n_subjects %||% NA_real_),
    solve_for        = if (inherits(x, "slope_sample_size_floor")) "n_floor"
                       else if (inherits(x, "slope_power")) "power"
                       else "n",
    row.names        = row.names,
    stringsAsFactors = FALSE
  )
}
