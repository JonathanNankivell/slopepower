# Methods for the accessor generics of nlme and stats.
#
# A `slope_params` object answers the questions an `lme` fit does -- what are
# the fixed effects, the random-effects covariance, the residual SD, the
# sampling covariance of the estimates -- through the same generics, so a
# reader used to nlme need not learn the object's field names. Each method
# answers in the package's own terms rather than the internal model's: the
# fixed effects are the two slopes, not the `sp_*` coefficients `$fit` carries,
# and the random effects are those of the untreated (or case) group, which is
# what stage two uses. `$fit` is still there for the model itself.
#
# Bootstrap results get `confint()`, the stats generic for an interval.

#' nlme and stats accessors for slope parameters
#'
#' The generics an [nlme::lme()] fit answers to, answered in the terms of a
#' `slope_params` object: the quantities stage two uses, named as the rest of
#' the package names them. They work on objects from [slope_params_manual()]
#' too, except `vcov()`, which needs a fitted model.
#'
#' * `fixef()`: the slope of the untreated (or case) group, `slope`, and when
#'   there is a comparator its slope too, `slope_comparator`.
#' * `getVarCov()`: the 2 x 2 covariance matrix of the random intercept and
#'   slope, from `sigma2_intercept`, `sigma2_slope` and `cov_intercept_slope`.
#'   Under `comparator = "healthy"` these are the cases' own.
#' * `sigma()`: the residual standard deviation, `sqrt(sigma2_residual)`. With
#'   visit-specific variances (`varIdent()`), it is the residual SD at the
#'   first visit time; see `$residual$sd_ratio` for the others.
#' * `vcov()`: the sampling covariance matrix of the slopes `fixef()` returns,
#'   from the fitted model. A matrix of `NA` for an object with no fitted
#'   model.
#'
#' The fitted model itself is `$fit`, for nlme's own methods on the internal
#' coefficients.
#'
#' @param object,obj A `slope_params` object.
#' @param ... Ignored.
#'
#' @return `fixef()` a named numeric vector; `getVarCov()` and `vcov()` named
#'   matrices; `sigma()` a single number.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' fixef(pars)
#' getVarCov(pars)
#' sigma(pars)
#' vcov(pars)
#' sqrt(vcov(pars)["slope", "slope"])  # slope_se(pars)
#'
#' @seealso [slope_se()], the standard error of the slope alone.
#' @name slope_params_methods
NULL

#' @rdname slope_params_methods
#' @export
fixef.slope_params <- function(object, ...) {
  check_params(object, "fixef()")
  if (identical(object$comparator, "none")) {
    c(slope = object$slope)
  } else {
    c(slope = object$slope, slope_comparator = object$slope_comparator)
  }
}

#' @rdname slope_params_methods
#' @export
getVarCov.slope_params <- function(obj, ...) {
  check_params(obj, "getVarCov()")
  re_cov_matrix(obj$sigma2_intercept, obj$sigma2_slope, obj$cov_intercept_slope)
}

#' @rdname slope_params_methods
#' @export
sigma.slope_params <- function(object, ...) {
  check_params(object, "sigma()")
  sqrt(object$sigma2_residual)
}

#' @rdname slope_params_methods
#' @export
vcov.slope_params <- function(object, ...) {
  check_params(object, "vcov()")
  slope_vcov(object, "vcov()")
}

#' Confidence intervals from a bootstrap
#'
#' The intervals a bootstrap built, in the shape [stats::confint()] returns
#' them: a row for the bootstrapped statistic and a row for the slope the
#' resampling perturbed, each named for the quantity it is -- `n_per_arm` or
#' `n_total`, `power`, `tte`, `slope`. They are read from the object's
#' `<quantity>_lower` and `<quantity>_upper` fields, not recomputed, so they
#' are the intervals `print()` shows: BCa or percentile as the object records,
#' widened to even trial sizes for a sample size.
#'
#' @param object A bootstrap result, from [slope_sample_size_boot()],
#'   [slope_power_boot()], [slope_params_boot()] or the bound bootstraps.
#' @param parm Which rows, by quantity name -- any the object stores, so
#'   `"n_total"` and `"n_per_arm"` are both available for a sample size -- or
#'   by position among the default rows, as in [stats::confint()]. Defaults to
#'   the rows `print()` shows: the statistic on the display basis, and
#'   `"slope"`.
#' @param level The interval's confidence level. It must be the one the
#'   bootstrap was run at, which is the default; for another, rerun the
#'   bootstrap with that `level` (and the same `seed`).
#' @param ... Ignored.
#' @param per_arm For a bootstrapped sample size, which row the default shows:
#'   `TRUE` for `n_per_arm`, `FALSE` for `n_total`. Defaults to `NULL`,
#'   meaning the basis the bootstrap was called with, as in `print()`.
#'
#' @return A matrix with one row per `parm` and columns for the lower and
#'   upper limits, labelled by percentage as [stats::confint()] labels them.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' \donttest{
#' b <- slope_sample_size_boot(pars, c(0, 1, 2), effectiveness = 0.33,
#'                             R = 50, seed = 42)
#' confint(b)
#' confint(b, "n_total")
#' }
#' @export
confint.slope_bootstrap <- function(object, parm, level = object$level, ..., per_arm = NULL) {
  context <- "confint()"
  shown <- boot_rows(object, display_basis(object, per_arm, context))
  rows <- boot_parm(parm, shown, boot_quantities(object), context)
  check_confint_level(level, object$level, context)
  out <- cbind(vapply(rows, function(q) object[[paste0(q, "_lower")]], numeric(1L)),
               vapply(rows, function(q) object[[paste0(q, "_upper")]], numeric(1L)))
  dimnames(out) <- list(rows, confint_labels(object$level))
  out
}

#' Resolve `confint()`'s `parm` to quantity names
#'
#' By name among every quantity stored, or by position among those shown by
#' default; missing means all of those shown.
#' @noRd
boot_parm <- function(parm, shown, available, context) {
  if (missing(parm)) return(shown)
  rows <- if (is.numeric(parm)) {
    if (all(parm %in% seq_along(shown))) shown[parm]
  } else if (is.character(parm) && all(parm %in% available)) {
    parm
  }
  if (is.null(rows)) {
    stop(sprintf("%s: `parm` must be among %s, or positions among %s.", context,
                 paste(sQuote(available), collapse = ", "),
                 paste(sQuote(shown), collapse = ", ")), call. = FALSE)
  }
  rows
}

#' Refuse a `level` other than the one a bootstrap's intervals were built at
#' @noRd
check_confint_level <- function(level, built, context) {
  if (!is.numeric(level) || length(level) != 1L || !isTRUE(all.equal(level, built))) {
    stop(sprintf(paste0(
      "%s: this bootstrap's intervals were built at level = %s, and are not\n",
      "  recomputed here. Rerun the bootstrap with the level wanted (and the same\n",
      "  `seed`, for the same replicates)."), context, format(built)),
      call. = FALSE)
  }
  invisible(level)
}

#' `confint()`'s column labels, "2.5 %" and "97.5 %"
#' @noRd
confint_labels <- function(level) {
  paste(format(100 * boot_probs(level), trim = TRUE, scientific = FALSE, digits = 3), "%")
}

# Re-exported so `fixef(pars)` works without attaching nlme, as lme4 does for
# the same generic.
#' @importFrom nlme fixef
#' @export
nlme::fixef

#' @importFrom nlme getVarCov
#' @export
nlme::getVarCov
