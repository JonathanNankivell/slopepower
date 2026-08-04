# Layer 1: stage-one parameter estimation.
#
# Estimates the slope(s) and the between/within-subject variance components that
# stage two needs, either by fitting a linear mixed model to previously collected
# longitudinal data (`slope_params`) or by taking them directly from the
# literature (`slope_params_manual`).
#
# Internal column names -------------------------------------------------------
# All model fitting happens on a freshly assembled data frame whose columns have
# fixed internal names (`sp_y`, `sp_time`, `sp_subject`, ...). This is
# deliberate: variance components are extracted from the fitted object *by name*,
# and pinning the names here makes that extraction deterministic regardless of
# what the user called their variables. It is the mechanism by which this port
# avoids the positional `e(b)` indexing that makes the Stata implementation
# fragile.

# ---- group coercion ---------------------------------------------------------

#' Coerce a two-level group column to a 0/1 numeric indicator
#'
#' Accepts numeric 0/1, logical, two-level factor/character, or haven_labelled.
#' For factors and characters the first level maps to 0 and the second to 1.
#' @noRd
coerce_binary <- function(x, name, context) {
  if (inherits(x, "haven_labelled")) x <- as.numeric(x)
  if (is.logical(x)) return(as.numeric(x))
  if (is.factor(x) || is.character(x)) {
    f <- if (is.factor(x)) x else factor(x)
    if (nlevels(f) != 2L) {
      stop(sprintf("%s: `%s` must have exactly two levels; got %d.",
                   context, name, nlevels(f)), call. = FALSE)
    }
    return(as.numeric(f) - 1)
  }
  if (!is.numeric(x)) {
    stop(sprintf("%s: `%s` must be numeric, logical, factor or labelled.",
                 context, name), call. = FALSE)
  }
  u <- sort(unique(x[!is.na(x)]))
  if (length(u) != 2L) {
    stop(sprintf("%s: `%s` must have exactly two distinct values; got %d.",
                 context, name, length(u)), call. = FALSE)
  }
  if (!isTRUE(all.equal(u, c(0, 1)))) {
    stop(sprintf("%s: `%s` must be coded 0/1 (1 = %s); got values %s.",
                 context, name,
                 if (identical(name, "healthy")) "case" else "treated",
                 paste(format(u), collapse = "/")), call. = FALSE)
  }
  as.numeric(x)
}

#' Evaluate a bare column name (or expression) against `data`, then the caller
#' @noRd
eval_column <- function(expr, data, env, context, name) {
  if (is.null(expr)) return(NULL)
  out <- tryCatch(eval(expr, data, env), error = function(e) {
    stop(sprintf("%s: could not evaluate `%s`: %s", context, name,
                 conditionMessage(e)), call. = FALSE)
  })
  if (is.null(out)) {
    stop(sprintf("%s: `%s` evaluated to NULL.", context, name), call. = FALSE)
  }
  out
}

#' Coerce a time vector to numeric, unwrapping Date/POSIXct
#' @noRd
coerce_time <- function(x, context) {
  if (inherits(x, "Date")) return(as.numeric(x))
  if (inherits(x, "POSIXct")) return(as.numeric(x) / 86400)
  if (inherits(x, "haven_labelled")) x <- as.numeric(x)
  if (!is.numeric(x)) {
    stop(sprintf("%s: the time variable must be numeric (or a Date).", context),
         call. = FALSE)
  }
  as.numeric(x)
}

# ---- variance-component extraction -----------------------------------------

#' Pull the random-effects covariance for a named intercept/slope pair
#'
#' Extraction is by dimname, never by position.
#' @noRd
extract_re <- function(fit, int_name, slope_name, context) {
  G <- tryCatch(nlme::getVarCov(fit), error = function(e) {
    stop(sprintf("%s: could not extract the random-effects covariance: %s",
                 context, conditionMessage(e)), call. = FALSE)
  })
  nm <- dimnames(G)[[1L]]
  missing <- setdiff(c(int_name, slope_name), nm)
  if (length(missing)) {
    stop(sprintf("%s: expected random effect(s) %s in the fitted model; found %s.",
                 context, paste(sQuote(missing), collapse = ", "),
                 paste(sQuote(nm), collapse = ", ")), call. = FALSE)
  }
  list(sigma2_intercept = as.numeric(G[int_name, int_name]),
       sigma2_slope     = as.numeric(G[slope_name, slope_name]),
       sigma_cov        = as.numeric(G[int_name, slope_name]))
}

#' Residual variance for one level of a `varIdent` structure, by level name
#'
#' `nlme` parameterises heteroscedastic residuals as sigma * delta_g with
#' delta = 1 at the reference level, so the variance for group g is
#' (sigma * delta_g)^2. `allCoef = TRUE` returns every level including the
#' reference, named, which lets us look the level up rather than index it.
#' @noRd
extract_residual <- function(fit, level = NULL, context) {
  s <- stats::sigma(fit)
  vs <- fit$modelStruct$varStruct
  if (is.null(vs) || is.null(level)) return(s^2)
  dl <- stats::coef(vs, unconstrained = FALSE, allCoef = TRUE)
  if (!level %in% names(dl)) {
    stop(sprintf("%s: residual level %s not found; have %s.",
                 context, sQuote(level), paste(sQuote(names(dl)), collapse = ", ")),
         call. = FALSE)
  }
  (s * as.numeric(dl[[level]]))^2
}

# ---- model fitting ----------------------------------------------------------

#' lme control settings used for every fit
#' @noRd
slope_lme_control <- function() {
  nlme::lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50,
                   opt = "optim", tolerance = 1e-7, msTol = 1e-8,
                   returnObject = FALSE)
}

#' Fit while muffling the structurally inevitable singular-precision warning
#'
#' The two-block random-effects structure used for `comparator = "healthy"`
#' gives every subject a design matrix with two all-zero columns (a case loads
#' only on the case block, a control only on the control block). `nlme` notes
#' this as "Singular precision matrix in level -1, block 1". It is benign: the
#' unloaded random effects contribute nothing to the likelihood. This was
#' verified empirically -- the joint fit reproduces the two separate per-group
#' fits to six significant figures (see `common_variance` note in
#' [slope_params()]).
#' @noRd
fit_quietly <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Singular precision matrix", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' @noRd
fit_healthy_model <- function(dat, reduced, ctrl) {
  comparator_block <- if (reduced) {
    nlme::pdIdent(~ sp_control - 1)
  } else {
    nlme::pdSymm(~ sp_control + sp_control_time - 1)
  }
  rand <- list(sp_subject = nlme::pdBlocked(list(
    nlme::pdSymm(~ sp_case + sp_case_time - 1),
    comparator_block)))
  fit_quietly(
    nlme::lme(sp_y ~ sp_case * sp_time,
              random  = rand,
              weights = nlme::varIdent(form = ~ 1 | sp_grp),
              data    = dat,
              method  = "REML",
              control = ctrl)
  )
}

# ---- main entry point -------------------------------------------------------

#' Estimate slope and variance parameters from longitudinal data
#'
#' Stage one of the two-stage sample-size method of Nash et al. (2021). Fits a
#' linear mixed model to previously collected longitudinal data and extracts the
#' slope(s) and the between- and within-subject variance components that
#' [slope_power()] needs.
#'
#' @param formula A two-sided formula `outcome ~ time | subject`. The time term
#'   may be an expression, e.g. `sdmt ~ I(as.numeric(vdate) / 365) | id` to work
#'   in years when visits are recorded as dates. Choose units on the order of the
#'   study duration: fitting on a badly scaled axis (days over several years, say)
#'   leaves the random-slope variance near zero and the REML optimiser converges
#'   less precisely. There is no `scale()` argument as there is in Stata --
#'   express `visits` in [trial_design()] in whatever units are used here.
#' @param data A data frame in long format, one row per measurement.
#' @param healthy Optional bare column name identifying healthy controls, coded
#'   `1` for cases (subjects with the disease) and `0` for healthy controls.
#'   Mutually exclusive with `treated`.
#' @param treated Optional bare column name identifying the treated arm of a
#'   previously conducted trial, coded `1` for treated and `0` for the control
#'   arm. Mutually exclusive with `healthy`.
#' @param origin `"subject"` (default) shifts each subject's time so their first
#'   visit is time zero, reproducing the Stata command's behaviour and ensuring
#'   the random intercept is estimated at baseline. `"none"` leaves time as
#'   supplied.
#' @param common_variance Controls the random-effects structure for healthy
#'   controls. `NULL` (default) fits the full model and falls back automatically
#'   if it fails to converge; `TRUE` forces the reduced structure (a random
#'   intercept only, equivalent to the Stata `nocontvar` option); `FALSE` forces
#'   the full structure and errors on failure. Ignored unless `healthy` is given.
#' @param na.action Applied to the assembled model frame. Defaults to
#'   [stats::na.omit()].
#'
#' @details
#' Three scenarios are supported, matching paper section 2.3:
#'
#' * neither `healthy` nor `treated`: a single group of untreated subjects with
#'   the disease. The target treatment effect will be measured toward a slope of
#'   zero.
#' * `healthy`: observational data containing both cases and healthy controls.
#'   The target effect will be measured toward the healthy-control slope.
#' * `treated`: data from a previous trial. The observed treatment effect is
#'   available via `target = "observed"` in [slope_power()].
#'
#' The returned variance components are always those of the **untreated / case**
#' group. When `healthy` is supplied the controls contribute only their slope;
#' their variance components are estimated and discarded, per paper section 2.3.
#'
#' Note that for the `healthy` scenario the model factorises exactly into two
#' independent fits, one per group: the fixed effects `y ~ case * time` span the
#' same column space as separate per-group intercepts and slopes, the
#' random-effects blocks are independent, and the residual variances are
#' separate. Consequently `common_variance` cannot change the estimates returned
#' for the cases -- it only affects how many nuisance parameters are estimated
#' for the controls, and therefore whether the fit converges at all.
#'
#' @return An object of class `"slope_params"`.
#'
#' @references
#' Nash, S., Morgan, K. E., Frost, C. and Mulick, A. (2021). Power and
#' sample-size calculations for trials that compare slopes over time:
#' Introducing the slopepower command. \emph{The Stata Journal} 21(3): 575--601.
#'
#' @seealso [slope_params_manual()] to supply parameters directly,
#'   [slope_power()] for stage two.
#' @export
slope_params <- function(formula, data,
                         healthy = NULL, treated = NULL,
                         origin = c("subject", "none"),
                         common_variance = NULL,
                         na.action = stats::na.omit) {
  context <- "slope_params()"
  cl <- match.call()
  origin <- match.arg(origin)

  if (missing(data) || !is.data.frame(data)) {
    data <- tryCatch(as.data.frame(data), error = function(e) {
      stop(sprintf("%s: `data` must be a data frame.", context), call. = FALSE)
    })
  } else {
    data <- as.data.frame(data)
  }

  parts <- parse_slope_formula(formula, context)
  if (is.null(parts$subject)) {
    stop(sprintf("%s: `formula` must name the subject identifier, e.g. `%s ~ %s | id`.",
                 context, deparse(parts$outcome), deparse(parts$time)),
         call. = FALSE)
  }

  env <- environment(formula) %||% parent.frame()

  y       <- eval_column(parts$outcome, data, env, context, "outcome")
  tim     <- coerce_time(eval_column(parts$time, data, env, context, "time"), context)
  subject <- eval_column(parts$subject, data, env, context, "subject")

  if (!is.numeric(y)) {
    if (inherits(y, "haven_labelled")) y <- as.numeric(y) else
      stop(sprintf("%s: the outcome must be numeric.", context), call. = FALSE)
  }

  healthy_expr <- substitute(healthy)
  treated_expr <- substitute(treated)
  if (!is.null(healthy_expr) && !is.null(treated_expr)) {
    stop(sprintf("%s: supply only one of `healthy` and `treated`, not both.", context),
         call. = FALSE)
  }

  comparator <- if (!is.null(healthy_expr)) "healthy" else
    if (!is.null(treated_expr)) "treated" else "none"

  grp <- NULL
  if (comparator != "none") {
    gexpr <- if (comparator == "healthy") healthy_expr else treated_expr
    graw  <- eval_column(gexpr, data, env, context, comparator)
    grp   <- coerce_binary(graw, comparator, context)
    if (length(grp) != length(y)) {
      stop(sprintf("%s: `%s` has length %d but the data have %d rows.",
                   context, comparator, length(grp), length(y)), call. = FALSE)
    }
  }

  if (!is.null(common_variance) && comparator != "healthy") {
    warning(sprintf("%s: `common_variance` applies only when `healthy` is supplied; ignoring it.",
                    context), call. = FALSE)
    common_variance <- NULL
  }

  n <- length(y)
  if (length(tim) != n || length(subject) != n) {
    stop(sprintf("%s: outcome, time and subject must be the same length.", context),
         call. = FALSE)
  }

  dat <- data.frame(sp_y = y, sp_time = tim,
                    sp_subject = factor(as.character(subject)),
                    stringsAsFactors = FALSE)
  if (!is.null(grp)) dat$sp_case <- grp

  dat <- na.action(dat)
  if (nrow(dat) < 3L) {
    stop(sprintf("%s: fewer than 3 usable observations after removing missing values.",
                 context), call. = FALSE)
  }
  dat$sp_subject <- droplevels(dat$sp_subject)

  # per-subject time origin
  time_shifted <- FALSE
  if (origin == "subject") {
    first <- stats::ave(dat$sp_time, dat$sp_subject, FUN = min)
    if (any(abs(first) > 1e-12)) {
      time_shifted <- TRUE
      message(sprintf(paste0("%s: time did not start at zero for all subjects. ",
                             "Times have been shifted so each subject's first ",
                             "visit is time zero."), context))
    }
    dat$sp_time <- dat$sp_time - first
  }

  if (comparator != "none" && length(unique(dat$sp_case)) != 2L) {
    stop(sprintf("%s: `%s` must contain both groups after removing missing values.",
                 context, comparator), call. = FALSE)
  }

  ctrl <- slope_lme_control()
  reduced_used <- FALSE

  if (comparator == "none") {
    fit <- nlme::lme(sp_y ~ sp_time, random = ~ sp_time | sp_subject,
                     data = dat, method = "REML", control = ctrl)
    b <- nlme::fixef(fit)
    slope <- unname(b[["sp_time"]])
    slope_comparator <- NA_real_
    re <- extract_re(fit, "(Intercept)", "sp_time", context)
    s2r <- extract_residual(fit, NULL, context)

  } else if (comparator == "treated") {
    # Stata: mixed y time placebo#c.time || subject: time, cov(uns)
    # One common intercept (randomisation implies equal baselines), separate
    # slopes. A numeric placebo indicator keeps the coefficient mapping explicit.
    dat$sp_placebo_time <- (1 - dat$sp_case) * dat$sp_time
    fit <- nlme::lme(sp_y ~ sp_time + sp_placebo_time,
                     random = ~ sp_time | sp_subject,
                     data = dat, method = "REML", control = ctrl)
    b <- nlme::fixef(fit)
    slope            <- unname(b[["sp_time"]] + b[["sp_placebo_time"]])  # control arm
    slope_comparator <- unname(b[["sp_time"]])                           # treated arm
    re <- extract_re(fit, "(Intercept)", "sp_time", context)
    s2r <- extract_residual(fit, NULL, context)

  } else {
    dat$sp_control      <- 1 - dat$sp_case
    dat$sp_case_time    <- dat$sp_case * dat$sp_time
    dat$sp_control_time <- dat$sp_control * dat$sp_time
    dat$sp_grp <- factor(ifelse(dat$sp_case == 1, "case", "control"),
                         levels = c("control", "case"))

    fit <- NULL
    if (!isTRUE(common_variance)) {
      fit <- tryCatch(fit_healthy_model(dat, reduced = FALSE, ctrl = ctrl),
                      error = function(e) e)
      if (inherits(fit, "error")) {
        if (isFALSE(common_variance)) {
          stop(sprintf(paste0("%s: the full model did not converge and ",
                              "`common_variance = FALSE` forbids the reduced ",
                              "structure. Underlying error: %s"),
                       context, conditionMessage(fit)), call. = FALSE)
        }
        message(sprintf(paste0("%s: the full random-effects structure for healthy ",
                               "controls did not converge; falling back to a ",
                               "random intercept only for controls (equivalent to ",
                               "the Stata `nocontvar` option). This does not affect ",
                               "the estimates returned for cases."), context))
        fit <- NULL
      }
    }
    if (is.null(fit)) {
      reduced_used <- TRUE
      fit <- tryCatch(fit_healthy_model(dat, reduced = TRUE, ctrl = ctrl),
                      error = function(e) {
                        stop(sprintf("%s: the mixed model did not converge: %s",
                                     context, conditionMessage(e)), call. = FALSE)
                      })
    }

    b <- nlme::fixef(fit)
    slope            <- unname(b[["sp_time"]] + b[["sp_case:sp_time"]])  # cases
    slope_comparator <- unname(b[["sp_time"]])                            # controls
    re  <- extract_re(fit, "sp_case", "sp_case_time", context)
    s2r <- extract_residual(fit, "case", context)
  }

  new_slope_params(
    slope            = slope,
    slope_comparator = slope_comparator,
    comparator       = comparator,
    sigma2_intercept = re$sigma2_intercept,
    sigma2_slope     = re$sigma2_slope,
    sigma_cov        = re$sigma_cov,
    sigma2_residual  = s2r,
    n_obs            = nrow(dat),
    n_subjects       = nlevels(dat$sp_subject),
    common_variance  = reduced_used,
    time_shifted     = time_shifted,
    fit              = fit,
    call             = cl,
    context          = context
  )
}

# ---- direct construction ----------------------------------------------------

#' Construct slope parameters directly
#'
#' Builds a `"slope_params"` object from values supplied by hand, for planning a
#' trial from published estimates when no suitable dataset is available. Paper
#' section 2.6 anticipates exactly this situation.
#'
#' @param slope Slope of the untreated (or case) group, per unit time.
#' @param sigma2_intercept Between-subject variance of random intercepts.
#' @param sigma2_slope Between-subject variance of random slopes.
#' @param sigma_cov Covariance of random intercepts and slopes.
#' @param sigma2_residual Within-subject residual variance.
#' @param slope_comparator Slope of the healthy controls or the treated arm.
#'   Required unless `comparator = "none"`.
#' @param comparator One of `"none"`, `"healthy"` or `"treated"`.
#'
#' @return An object of class `"slope_params"`.
#' @seealso [slope_params()] to estimate these from data.
#' @export
slope_params_manual <- function(slope,
                                sigma2_intercept, sigma2_slope,
                                sigma_cov, sigma2_residual,
                                slope_comparator = NA_real_,
                                comparator = c("none", "healthy", "treated")) {
  context <- "slope_params_manual()"
  cl <- match.call()
  comparator <- match.arg(comparator)

  check_scalar(slope, "slope", context)
  check_variance(sigma2_intercept, "sigma2_intercept", context)
  check_variance(sigma2_slope, "sigma2_slope", context)
  check_scalar(sigma_cov, "sigma_cov", context)
  check_variance(sigma2_residual, "sigma2_residual", context)

  if (comparator == "none") {
    slope_comparator <- NA_real_
  } else {
    if (length(slope_comparator) != 1L || is.na(slope_comparator)) {
      stop(sprintf("%s: `slope_comparator` is required when `comparator` is %s.",
                   context, sQuote(comparator)), call. = FALSE)
    }
    check_scalar(slope_comparator, "slope_comparator", context)
  }

  new_slope_params(
    slope            = as.numeric(slope),
    slope_comparator = as.numeric(slope_comparator),
    comparator       = comparator,
    sigma2_intercept = as.numeric(sigma2_intercept),
    sigma2_slope     = as.numeric(sigma2_slope),
    sigma_cov        = as.numeric(sigma_cov),
    sigma2_residual  = as.numeric(sigma2_residual),
    n_obs            = NA_integer_,
    n_subjects       = NA_integer_,
    common_variance  = FALSE,
    time_shifted     = FALSE,
    fit              = NULL,
    call             = cl,
    context          = context
  )
}

#' Validate and build the object
#' @noRd
new_slope_params <- function(slope, slope_comparator, comparator,
                             sigma2_intercept, sigma2_slope, sigma_cov,
                             sigma2_residual, n_obs, n_subjects,
                             common_variance, time_shifted, fit, call,
                             context) {
  if (!is.finite(slope)) {
    stop(sprintf("%s: the estimated slope is not finite.", context), call. = FALSE)
  }
  for (nm in c("sigma2_intercept", "sigma2_slope", "sigma2_residual")) {
    v <- get(nm)
    if (!is.finite(v) || v <= 0) {
      stop(sprintf("%s: `%s` must be a positive number; got %s.",
                   context, nm, format(v)), call. = FALSE)
    }
  }
  if (!is.finite(sigma_cov)) {
    stop(sprintf("%s: `sigma_cov` is not finite.", context), call. = FALSE)
  }

  G <- matrix(c(sigma2_intercept, sigma_cov, sigma_cov, sigma2_slope), 2L, 2L)
  if (!is_positive_definite(G)) {
    stop(sprintf(paste0("%s: the implied random-effects covariance matrix is not ",
                        "positive definite (var_int = %g, var_slope = %g, cov = %g)."),
                 context, sigma2_intercept, sigma2_slope, sigma_cov), call. = FALSE)
  }

  structure(
    list(slope            = slope,
         slope_comparator = slope_comparator,
         comparator       = comparator,
         sigma2_intercept = sigma2_intercept,
         sigma2_slope     = sigma2_slope,
         sigma_cov        = sigma_cov,
         sigma2_residual  = sigma2_residual,
         n_obs            = n_obs,
         n_subjects       = n_subjects,
         common_variance  = common_variance,
         time_shifted     = time_shifted,
         fit              = fit,
         call             = call),
    class = "slope_params")
}

# ---- printing ---------------------------------------------------------------

#' Print stage-one slope parameters
#'
#' @param x A `"slope_params"` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.slope_params <- function(x, ...) {
  lab <- switch(x$comparator,
                none    = "single group (target: no change over time)",
                healthy = "cases and healthy controls",
                treated = "previous randomised trial")
  cat("Slope parameters (", lab, ")\n\n", sep = "")

  if (!is.na(x$n_obs)) {
    cat(fmt_line("number of observations in model", x$n_obs, digits = 0L), "\n")
    cat(fmt_line("number of participants in model", x$n_subjects, digits = 0L), "\n")
  } else {
    cat(fmt_line("source", "supplied directly"), "\n")
  }

  slope_label <- switch(x$comparator,
                        none    = "slope of cases",
                        healthy = "slope of cases",
                        treated = "slope of control arm")
  cat(fmt_line(slope_label, x$slope), "\n")
  if (!is.na(x$slope_comparator)) {
    comp_label <- if (x$comparator == "healthy")
      "slope of healthy controls" else "slope of experimental arm"
    cat(fmt_line(comp_label, x$slope_comparator), "\n")
    cat(fmt_line("observed difference in slopes",
                 x$slope - x$slope_comparator), "\n")
  }

  cat("\n")
  cat(fmt_line("variance of random intercepts", x$sigma2_intercept), "\n")
  cat(fmt_line("variance of random slopes", x$sigma2_slope), "\n")
  cat(fmt_line("covariance of intercept and slope", x$sigma_cov), "\n")
  cat(fmt_line("residual variance", x$sigma2_residual), "\n")

  if (isTRUE(x$common_variance)) {
    cat("\nNote: reduced random-effects structure used for healthy controls\n")
    cat("      (Stata `nocontvar`). Case estimates are unaffected.\n")
  }
  if (isTRUE(x$time_shifted)) {
    cat("\nNote: subject times were shifted so each first visit is time zero.\n")
  }
  invisible(x)
}
