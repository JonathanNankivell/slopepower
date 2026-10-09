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
  nm <- c("(Intercept)", "time")
  matrix(c(obj$sigma2_intercept, obj$cov_intercept_slope,
           obj$cov_intercept_slope, obj$sigma2_slope),
         nrow = 2L, dimnames = list(nm, nm))
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
  context <- "vcov()"
  check_params(object, context)
  nm <- names(fixef(object))
  out <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
  fit <- object$fit
  if (is.null(fit)) return(out)
  b <- nlme::fixef(fit)
  terms <- slope_terms(object, b, context)
  if (is.null(terms)) return(out)
  # Row 1 sums every part, as slope_se() does; row 2, when there is a
  # comparator, is its own slope: the first part alone.
  A <- rbind(as.numeric(names(b) %in% terms),
             if (length(nm) > 1L) as.numeric(names(b) == terms[[1L]]))
  out[] <- A %*% as.matrix(stats::vcov(fit)) %*% t(A)
  out
}

#' Confidence intervals from a bootstrap
#'
#' The interval a bootstrap built, in the shape [stats::confint()] returns
#' one: a row for the bootstrapped statistic, named as `statistic` names it
#' (`"n"`, `"power"`, `"tte"` or `"slope"`), and a row for the slope the
#' resampling perturbed. Both are read from the object, not recomputed, so
#' they are the intervals `print()` shows -- BCa or percentile as the object
#' records, widened to even sample sizes for `n`, and for `n` on the same
#' basis `print()` uses: per arm unless the bootstrap was run with
#' `per_arm = FALSE`. The row is then named `n_per_arm`, so the basis is
#' never in doubt.
#'
#' @param object A bootstrap result, from [slope_sample_size_boot()],
#'   [slope_power_boot()], [slope_params_boot()] or the bound bootstraps.
#' @param parm Which rows, by name -- the statistic's and `"slope"` -- or by
#'   position, as in [stats::confint()]. Defaults to both (one row, for
#'   [slope_params_boot()], whose statistic is the slope).
#' @param level The interval's confidence level. It must be the one the
#'   bootstrap was run at, which is the default; for another, rerun the
#'   bootstrap with that `level` (and the same `seed`).
#' @param ... Ignored.
#' @param per_arm For a bootstrapped `n`: `TRUE` for participants per arm,
#'   `FALSE` for the trial total. Defaults to `NULL`, meaning the basis the
#'   bootstrap was called with, as in `print()`. `parm` may name the row
#'   `"n"` either way, or `"n_per_arm"` when that is what it prints as.
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
#' confint(b, "n")
#' }
#' @export
confint.slope_bootstrap <- function(object, parm, level = object$level, ..., per_arm = NULL) {
  context <- "confint()"
  # The basis print.slope_bootstrap() shows, by the same rule (boot_summary_frame()).
  halve <- isTRUE(object$lattice) && display_basis(object, per_arm, context)
  rows <- unique(c(object$statistic, "slope"))
  shown <- ifelse(halve & rows == "n", "n_per_arm", rows)
  if (missing(parm)) parm <- seq_along(rows)
  idx <- if (is.numeric(parm)) {
    if (all(parm %in% seq_along(rows))) parm
  } else if (is.character(parm)) {
    # The stored name or the printed one, so a row can be asked for as shown.
    ifelse(parm %in% rows, match(parm, rows), match(parm, shown))
  }
  if (is.null(idx) || anyNA(idx)) {
    stop(sprintf("%s: `parm` must be among %s, or their positions.", context,
                 paste(sQuote(unique(c(rows, shown))), collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !isTRUE(all.equal(level, object$level))) {
    stop(sprintf(paste0(
      "%s: this bootstrap's intervals were built at level = %s, and are not\n",
      "  recomputed here. Rerun the bootstrap with the level wanted (and the same\n",
      "  `seed`, for the same replicates)."), context, format(object$level)),
      call. = FALSE)
  }
  ci <- list(object$ci / if (halve) 2 else 1, object$slope_ci)[seq_along(rows)]
  a <- (1 - object$level) / 2
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
  out <- do.call(rbind, ci[idx])
  dimnames(out) <- list(shown[idx], pct)
  out
}

# Re-exported so `fixef(pars)` works without attaching nlme, as lme4 does for
# the same generic.
#' @importFrom nlme fixef
#' @export
nlme::fixef

#' @importFrom nlme getVarCov
#' @export
nlme::getVarCov
