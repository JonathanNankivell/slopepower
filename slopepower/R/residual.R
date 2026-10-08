# Layer 1 (continued) -- the within-participant residual structure.
#
# The model of Nash et al. (2021) gives each participant independent residuals
# with one variance, sigma2_residual. This file lets the residuals of one
# participant have a structured covariance instead, written the way `nlme`
# writes it: a `corStruct` for the correlation (`correlation = corAR1()`) and a
# `varIdent` for visit-specific variances (`weights = varIdent(form = ~ 1 |
# visit)`), so that `correlation = corSymm()` together with that `varIdent` is
# the fully unstructured residual covariance.
#
# What a `slope_params` object stores is plain data, not `nlme` objects -- the
# `residual` field of CONTRACT.md section 2 -- because stage two has to evaluate
# the structure at visit times that need not have occurred in the data, and
# because `slope_params_manual()` has to be able to state one without a fit.
# The correlation functions below restate `nlme`'s definitions; the test suite
# checks them against `nlme::getVarCov(type = "marginal")` on fitted models, so
# the two cannot quietly disagree.

#' The correlation structures a `slope_params` object can carry
#'
#' `nlme` class names, so that what the user typed, what the fitted object
#' stores and what the printed parameters say are one vocabulary. `"none"`
#' means independent residuals -- possibly with visit-specific variances.
#' @noRd
RESIDUAL_CORRELATIONS <- c("none", "corAR1", "corCAR1", "corExp", "corGaus", "corSymm")

#' Human-readable name of a correlation structure, for printing and messages
#' @noRd
residual_label <- function(correlation) {
  switch(correlation,
         none    = "independent",
         corAR1  = "AR(1) (corAR1)",
         corCAR1 = "continuous-time AR(1) (corCAR1)",
         corExp  = "exponential (corExp)",
         corGaus = "Gaussian (corGaus)",
         corSymm = "unstructured (corSymm)")
}

#' Does this residual structure exist only at a fixed set of visit times?
#'
#' `corSymm()` has one correlation per *pair of visits* and `varIdent()` one
#' variance per visit, so neither says anything about a time that was not one
#' of those visits. Every other supported structure is a function of time and
#' can be evaluated anywhere. Stage two, the floor and the printing all branch
#' on this one fact, so it is asked one way.
#' @noRd
residual_on_grid <- function(residual) {
  !is.null(residual) && !is.null(residual$times)
}

# ---- user arguments to slope_params() ---------------------------------------

#' Validate `correlation` and `weights` as given to [slope_params()]
#'
#' Both are `nlme` constructors, used the way `nlme` users already write them.
#' The one thing this port decides for them is the `form`: residuals are always
#' correlated *within a participant* and *as a function of the time variable*
#' of `formula`, because that is the only reading stage two can evaluate at a
#' planned schedule. So `corAR1()` with `nlme`'s default `form = ~ 1` -- which
#' in `nlme` itself means "order of observation", and miscounts the lag across
#' a missed visit -- is taken as `~ time | subject`. A `form` naming those two
#' variables explicitly is accepted, since that is the correct `nlme` spelling
#' of the same thing; any other `form` would describe a different model, and is
#' refused rather than silently replaced.
#'
#' @return `NULL` when neither is supplied, which is the original model exactly.
#'   Otherwise a list: `correlation` (an element of [RESIDUAL_CORRELATIONS]),
#'   `struct` (the user's `corStruct`, or `NULL`) and `by_visit` (`TRUE` when
#'   `weights` gives each visit its own variance).
#' @noRd
residual_spec <- function(correlation, weights, parts, context) {
  if (is.null(correlation) && is.null(weights)) return(NULL)
  cls <- "none"
  if (!is.null(correlation)) {
    cls <- check_correlation_class(correlation, context)
    check_struct_form(correlation, "correlation", parts, context)
  }
  by_visit <- FALSE
  if (!is.null(weights)) {
    check_weights_class(weights, context)
    time_text <- paste(deparse(parts$time), collapse = " ")
    cov <- nlme::getCovariateFormula(weights)
    grp <- nlme::getGroupsFormula(weights)
    if (!identical(cov[[2L]], 1) || is.null(grp) || !identical(grp[[2L]], parts$time)) {
      stop(sprintf(paste0(
        "%s: `weights` must give each visit time its own residual variance, written as ",
        "`varIdent(form = ~ 1 | %s)`, with the time variable of `formula` as the stratum.\n",
        "  `varIdent()` with no stratum is a constant variance in nlme, i.e. no `weights` ",
        "at all; other strata would make the variance depend on something other than time."),
        context, time_text), call. = FALSE)
    }
    if (length(unclass(weights)) || isTRUE(attr(weights, "fixed"))) {
      stop(sprintf(paste0(
        "%s: `weights` must not carry initial or fixed values here; slope_params() ",
        "estimates the variance ratios. To state them, use slope_params_manual()."),
        context), call. = FALSE)
    }
    by_visit <- TRUE
  }
  list(correlation = cls, struct = correlation, by_visit = by_visit)
}

#' The class of a user-supplied correlation structure, or a reasoned refusal
#'
#' Shared by both routes into the class. `corCompSymm()` is refused by name:
#' an exchangeable residual correlation adds the same amount to every pair of
#' a participant's visits, which is exactly what the random intercept already
#' does, so the two cannot be told apart and the fit is not identified.
#' @noRd
check_correlation_class <- function(correlation, context) {
  if (!inherits(correlation, "corStruct")) {
    stop(sprintf(paste0(
      "%s: `correlation` must be an nlme correlation structure, e.g. `nlme::corAR1()` ",
      "or `nlme::corSymm()`."), context), call. = FALSE)
  }
  cls <- class(correlation)[1L]
  if (identical(cls, "corCompSymm")) {
    stop(sprintf(paste0(
      "%s: corCompSymm() cannot be combined with the random intercept: an exchangeable ",
      "residual correlation adds the same covariance to every pair of a participant's ",
      "visits, which is what the random intercept already models, so the two are not ",
      "identified. Leave `correlation` unset -- the random intercept is that structure."),
      context), call. = FALSE)
  }
  if (!cls %in% RESIDUAL_CORRELATIONS[-1L]) {
    stop(sprintf("%s: correlation structure %s is not supported; use one of %s.",
                 context, sQuote(cls),
                 paste0(RESIDUAL_CORRELATIONS[-1L], "()", collapse = ", ")),
         call. = FALSE)
  }
  cls
}

#' Refuse a variance function other than `varIdent()`
#' @noRd
check_weights_class <- function(weights, context) {
  if (!inherits(weights, "varIdent")) {
    stop(sprintf(paste0(
      "%s: `weights` must be an nlme `varIdent()` stratified by visit time; ",
      "other variance functions are not supported."), context), call. = FALSE)
  }
  invisible(weights)
}

#' Accept a correlation `form` only if it means "within participant, by time"
#' @noRd
check_struct_form <- function(struct, name, parts, context) {
  cov <- nlme::getCovariateFormula(struct)
  grp <- nlme::getGroupsFormula(struct)
  cov_ok <- identical(cov[[2L]], 1) || identical(cov[[2L]], parts$time)
  grp_ok <- is.null(grp) || identical(grp[[2L]], parts$subject)
  if (!cov_ok || !grp_ok) {
    stop(sprintf(paste0(
      "%s: the `form` of `%s` must be omitted or be `~ %s | %s`. Residuals are ",
      "correlated within a participant as a function of the time variable of ",
      "`formula`; slope_params() sets that form itself."),
      context, name, paste(deparse(parts$time), collapse = " "),
      paste(deparse(parts$subject), collapse = " ")), call. = FALSE)
  }
  invisible(struct)
}

#' Check that the data can carry the residual structure asked for
#'
#' Run on the frame actually fitted -- after `na.action` and the time origin --
#' so that what is checked is what `nlme` will see. Three checks, each one an
#' error `nlme` would otherwise raise in its own words, or not raise at all:
#'
#' * at most one measurement per participant per time, without which a
#'   correlation *as a function of time* is undefined;
#' * integer times for `corAR1()`, which `nlme` defines on integer lags only;
#' * for `corSymm()` and `varIdent()`, a visit schedule shared across
#'   participants -- one parameter per distinct time means nothing when every
#'   participant was seen on different days.
#'
#' @return The sorted distinct times (the visit grid) when the structure needs
#'   one, else `NULL`.
#' @noRd
check_residual_data <- function(spec, time, subject, context) {
  if (anyDuplicated(data.frame(subject, time_key(time)))) {
    stop(sprintf(paste0(
      "%s: a residual correlation needs at most one measurement per participant at each ",
      "time, but some participants have two or more at the same time."),
      context), call. = FALSE)
  }
  if (identical(spec$correlation, "corAR1") && !all(is_whole(time))) {
    stop(sprintf(paste0(
      "%s: corAR1() is a discrete-time process, defined by nlme on integer times only, ",
      "and the time variable takes non-integer values. Use corCAR1(), its continuous-time ",
      "counterpart, which gives the same correlation Phi^|t - s| at any times."),
      context), call. = FALSE)
  }
  if (!identical(spec$correlation, "corSymm") && !spec$by_visit) return(NULL)

  key  <- time_key(time)
  grid <- sort(unique(key))
  seen <- tapply(subject, factor(key, levels = grid), function(s) length(unique(s)))
  if (length(grid) < 2L || any(seen < 2L)) {
    stop(sprintf(paste0(
      "%s: %s has a parameter for every distinct visit time, so it needs a visit schedule ",
      "shared across participants. Here %d of the %d distinct times are attended by fewer ",
      "than two participants. Record time as the scheduled visit rather than the date it ",
      "took place, or use a structure that is a function of time, such as corCAR1() or ",
      "corExp()."),
      context,
      if (spec$by_visit) "`weights = varIdent()`" else "corSymm()",
      sum(seen < 2L), length(grid)), call. = FALSE)
  }
  grid
}

#' One key per distinct time, robust to representation noise in the data
#' @noRd
time_key <- function(t) signif(t, 10L)

#' Is each element a whole number, to within representation noise?
#' @noRd
is_whole <- function(t) abs(t - round(t)) <= 1e-8 * pmax(1, abs(t))

#' Labels for the visit grid, as `varIdent()` levels and in printed output
#' @noRd
time_labels <- function(times) as.character(signif(times, 12L))

#' Add the columns the residual structure is fitted against to `dat`
#'
#' `sp_pos` indexes the visit grid for `corSymm()`, which `nlme` requires as
#' integer positions. `sp_visit` is the `varIdent()` stratum, labelled by time
#' so that `print(params$fit)` names the visits as the user knows them; under
#' `healthy` it is crossed with the group, because each group already has its
#' own residual variance and keeps it.
#' @noRd
add_residual_columns <- function(dat, spec, grid, comparator) {
  if (is.null(grid)) return(dat)
  pos <- match(time_key(dat$sp_time), grid)
  if (identical(spec$correlation, "corSymm")) dat$sp_pos <- pos
  if (spec$by_visit) {
    lab <- time_labels(grid)
    dat$sp_visit <- if (identical(comparator, "healthy")) {
      factor(paste(ifelse(dat$sp_case == 1, "case", "control"), lab[pos], sep = ":"),
             levels = c(paste("control", lab, sep = ":"), paste("case", lab, sep = ":")))
    } else {
      factor(lab[pos], levels = lab)
    }
  }
  dat
}

#' A call constructing a correlation structure on the internal columns
#'
#' Built as a *call*, not an object, so that `params$fit$call` reads as the
#' model that was fitted -- `correlation = nlme::corCAR1(form = ~sp_time |
#' sp_subject)` -- rather than as a deparsed object. Used both for the fit
#' (carrying any initial or fixed value the user gave) and by the bootstrap's
#' refitter (carrying the fitted value, when it was fixed).
#' @noRd
correlation_call <- function(correlation, value = NULL, nugget = FALSE, fixed = FALSE,
                             time = quote(sp_time), subject = quote(sp_subject),
                             position = quote(sp_pos)) {
  covariate <- if (identical(correlation, "corSymm")) position else time
  args <- list(form = call("~", call("|", covariate, subject)))
  if (length(value)) args$value <- unname(value)
  if (isTRUE(nugget)) args$nugget <- TRUE
  if (isTRUE(fixed)) args$fixed <- TRUE
  as.call(c(list(call("::", quote(nlme), as.name(correlation))), args))
}

#' The `correlation` and `weights` calls for a fit, by comparator
#'
#' Without a residual structure this is the original model exactly: no
#' correlation, and the per-group `varIdent()` under `healthy` only.
#' @noRd
residual_calls <- function(spec, comparator) {
  healthy <- identical(comparator, "healthy")
  weights <- if (!is.null(spec) && spec$by_visit) {
    quote(nlme::varIdent(form = ~ 1 | sp_visit))
  } else if (healthy) {
    quote(nlme::varIdent(form = ~ 1 | sp_grp))
  }
  correlation <- if (!is.null(spec) && !identical(spec$correlation, "none")) {
    s <- spec$struct
    correlation_call(spec$correlation, value = struct_value(s, spec$correlation),
                     nugget = isTRUE(attr(s, "nugget")), fixed = isTRUE(attr(s, "fixed")))
  }
  list(correlation = correlation, weights = weights)
}

#' The starting (or fixed) value a user's correlation structure carries
#'
#' On the natural scale, so that it can be written back into a constructor
#' call. `corExp()`, `corGaus()` and `corSymm()` store what was typed until
#' they are initialised; `corAR1()` and `corCAR1()` store it transformed, and
#' are initialised on a two-visit dummy participant to read it back -- the
#' value does not depend on the data for either.
#' @noRd
struct_value <- function(s, correlation) {
  v <- as.numeric(unclass(s))
  if (!length(v) || !correlation %in% c("corAR1", "corCAR1")) return(v)
  attr(s, "formula") <- stats::as.formula("~ sp_time | sp_subject", env = baseenv())
  s <- nlme::Initialize(s, data.frame(sp_time = c(0, 1), sp_subject = 1))
  as.numeric(stats::coef(s, unconstrained = FALSE))
}

# ---- extraction from a fit --------------------------------------------------

#' The residual variance and structure of a fitted model
#'
#' Replaces the bare [extract_residual()] call when a structure was fitted.
#' `sigma2_residual` keeps its meaning throughout -- the variance of one
#' residual of the untreated / case group -- and, with visit-specific
#' variances, is that variance at the *first* visit time, the reference level
#' of `varIdent()`; the others are stored as `nlme` reports them, as ratios of
#' standard deviations to it.
#' @noRd
residual_components <- function(fit, spec, grid, comparator, context) {
  level <- if (identical(comparator, "healthy")) "case" else NULL
  if (is.null(spec)) {
    return(list(sigma2_residual = extract_residual(fit, level, context), residual = NULL))
  }
  cs <- fit$modelStruct$corStruct
  coefs <- if (is.null(cs)) numeric() else
    correlation_coef(stats::coef(cs, unconstrained = FALSE), spec$correlation, length(grid))

  sd_ratio <- NULL
  if (spec$by_visit) {
    lab <- time_labels(grid)
    if (!is.null(level)) lab <- paste(level, lab, sep = ":")
    dl <- stats::coef(fit$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)
    if (!all(lab %in% names(dl))) {
      stop(sprintf("%s: residual levels %s not found; have %s.", context,
                   paste(sQuote(setdiff(lab, names(dl))), collapse = ", "),
                   paste(sQuote(names(dl)), collapse = ", ")), call. = FALSE)
    }
    delta <- as.numeric(dl[lab])
    s2r <- (stats::sigma(fit) * delta[1L])^2
    sd_ratio <- delta / delta[1L]
  } else {
    s2r <- extract_residual(fit, level, context)
  }
  list(sigma2_residual = s2r,
       residual = list(correlation = spec$correlation,
                       coef        = coefs,
                       fixed       = isTRUE(attr(spec$struct, "fixed")),
                       times       = if (!is.null(grid)) as.numeric(grid),
                       sd_ratio    = sd_ratio))
}

#' Name a correlation structure's natural-scale coefficients consistently
#'
#' `nlme` names them inconsistently between classes and even between
#' initialisations (`Phi` from a fit, `Phi1` from `Initialize()` on a time
#' covariate), so they are renamed here by position, once, to the names
#' [check_residual()] and [residual_cor()] read.
#' @noRd
correlation_coef <- function(x, correlation, n_times) {
  x <- as.numeric(x)
  switch(correlation,
    corAR1  = ,
    corCAR1 = c(Phi = x[[1L]]),
    corExp  = ,
    corGaus = if (length(x) == 2L) c(range = x[[1L]], nugget = x[[2L]]) else c(range = x[[1L]]),
    corSymm = x)
}

# ---- slope_params_manual() --------------------------------------------------

#' Read a stated residual structure from `nlme` constructors
#'
#' The same constructors as [slope_params()] accepts, now carrying values:
#' `corCAR1(0.6)`, `corExp(c(2, 0.1), nugget = TRUE)`, `corSymm(c(...))`, and
#' `varIdent(c("1" = 1.2, "2" = 1.5), form = ~ 1 | visit)`. The values are read
#' by initialising each object on the stated visit times -- the one way to get
#' them back on their natural scale, since an uninitialised `nlme` object stores
#' some on the natural scale and some transformed, class by class -- so `nlme`
#' applies its own range checks too.
#' @noRd
residual_manual <- function(correlation, weights, times, context) {
  if (is.null(correlation) && is.null(weights)) {
    if (!is.null(times)) {
      stop(sprintf(paste0("%s: `times` applies only with `correlation = corSymm()` or ",
                          "`weights = varIdent()`."), context), call. = FALSE)
    }
    return(NULL)
  }
  cls <- if (is.null(correlation)) "none" else check_correlation_class(correlation, context)
  on_grid <- identical(cls, "corSymm") || !is.null(weights)
  if (on_grid) {
    if (is.null(times)) {
      stop(sprintf(paste0("%s: `times` is required with %s: it gives the visit times the ",
                          "stated parameters belong to, e.g. `times = c(0, 1, 2, 3)`."),
                   context, if (identical(cls, "corSymm")) "corSymm()" else "varIdent()"),
           call. = FALSE)
    }
    if (!is.numeric(times) || length(times) < 2L || any(!is.finite(times)) ||
        is.unsorted(times, strictly = TRUE)) {
      stop(sprintf("%s: `times` must be a strictly increasing numeric vector of length >= 2.",
                   context), call. = FALSE)
    }
  } else if (!is.null(times)) {
    stop(sprintf(paste0("%s: `times` applies only with `correlation = corSymm()` or ",
                        "`weights = varIdent()`; %s() is a function of time."),
                 context, cls), call. = FALSE)
  }
  times <- if (on_grid) as.numeric(times)

  coefs <- numeric()
  if (!identical(cls, "none")) {
    # corAR1() and corCAR1() store a default (0 and 0.2) that cannot be told
    # from a stated value, so only the structures without one can be caught.
    if (!length(unclass(correlation))) {
      stop(sprintf(paste0("%s: state the value of the correlation parameter(s), e.g. ",
                          "`%s`; without one nlme would choose a starting value from data ",
                          "that do not exist here."),
                   context, switch(cls, corExp = "corExp(2)", corGaus = "corGaus(2)",
                                   corSymm = "corSymm(c(0.5, 0.3, 0.5))")), call. = FALSE)
    }
    probe <- times %||% c(0, 1)
    dummy <- data.frame(sp_time = probe, sp_pos = seq_along(probe), sp_subject = 1)
    cs <- correlation
    attr(cs, "formula") <- eval(correlation_call(cls)[["form"]])
    cs <- tryCatch(nlme::Initialize(cs, dummy), error = function(e) {
      stop(sprintf("%s: invalid `correlation`: %s", context, conditionMessage(e)), call. = FALSE)
    })
    coefs <- correlation_coef(stats::coef(cs, unconstrained = FALSE), cls, length(probe))
  }

  sd_ratio <- NULL
  if (!is.null(weights)) {
    check_weights_class(weights, context)
    sd_ratio <- manual_sd_ratio(weights, times, context)
  }
  list(correlation = cls, coef = coefs, fixed = FALSE, times = times, sd_ratio = sd_ratio)
}

#' Visit-specific SD ratios from a stated `varIdent()`, keyed by time
#'
#' `nlme` keys a `varIdent()`'s values by level name and ignores them (with a
#' warning) unless `form` names a stratum. Here the levels are the visit times,
#' so names are matched to `times` *as numbers* -- `"0.50"` is the visit at
#' 0.5 -- and the first time is the reference, with ratio 1.
#' @noRd
manual_sd_ratio <- function(weights, times, context) {
  v <- unclass(weights)
  attributes(v) <- list(names = names(v))
  if (is.null(nlme::getGroupsFormula(weights)) || !length(v) || is.null(names(v))) {
    stop(sprintf(paste0(
      "%s: state the visit-specific SD ratios by time, as nlme does by level, e.g. ",
      "`varIdent(c(\"%s\" = 1.2), form = ~ 1 | visit)`; values without a `form` stratum ",
      "are ignored by nlme."), context, time_labels(times[2L])), call. = FALSE)
  }
  # nlme stores a varIdent's values on the log scale until it is initialised.
  ratio <- exp(as.numeric(v))
  at <- vapply(suppressWarnings(as.numeric(names(v))), function(nm) {
    hit <- which(abs(times - nm) <= 1e-8 * max(1, abs(nm)))
    if (length(hit)) hit[1L] else NA_integer_
  }, integer(1L))
  if (anyNA(at)) {
    stop(sprintf("%s: `weights` names level(s) %s, which are not among `times` (%s).",
                 context, paste(sQuote(names(v)[is.na(at)]), collapse = ", "),
                 label_numeric(times)), call. = FALSE)
  }
  if (anyDuplicated(at)) {
    stop(sprintf("%s: `weights` names the same time more than once.", context), call. = FALSE)
  }
  if (any(at == 1L) && abs(ratio[at == 1L] - 1) > 1e-8) {
    stop(sprintf(paste0("%s: the first time, %s, is the reference level of `varIdent()`, ",
                        "so its SD ratio is 1 by definition; `sigma2_residual` is the ",
                        "variance there."), context, time_labels(times[1L])), call. = FALSE)
  }
  out <- rep(1, length(times))
  out[at] <- ratio
  out
}

# ---- validation -------------------------------------------------------------

#' Validate a `residual` field
#'
#' The one statement of its invariants, run by [new_slope_params()] at
#' construction and by [check_params()] at use, exactly as the variance
#' components themselves are. `NULL` -- independent residuals of constant
#' variance -- is always valid, and is also what an object built before this
#' field existed is read as.
#' @noRd
check_residual <- function(residual, context) {
  if (is.null(residual)) return(invisible(NULL))
  bad <- function(what) {
    stop(sprintf("%s: `params$residual` is not a valid residual structure: %s.",
                 context, what), call. = FALSE)
  }
  # `times` and `sd_ratio` are NULL for most structures, and assigning NULL to
  # a list element removes it, so their absence is read as NULL.
  if (!is.list(residual) || !all(c("correlation", "coef") %in% names(residual))) {
    bad("it must be a list with at least `correlation` and `coef`")
  }
  cls <- residual$correlation
  if (!is.character(cls) || length(cls) != 1L || !cls %in% RESIDUAL_CORRELATIONS) {
    bad(sprintf("`correlation` must be one of %s",
                paste(sQuote(RESIDUAL_CORRELATIONS), collapse = ", ")))
  }
  times <- residual$times
  if (!is.null(times) && (!is.numeric(times) || length(times) < 2L ||
                          any(!is.finite(times)) || is.unsorted(times, strictly = TRUE))) {
    bad("`times` must be strictly increasing finite numbers")
  }
  needs_times <- identical(cls, "corSymm") || !is.null(residual$sd_ratio)
  if (needs_times != !is.null(times)) {
    bad("`times` must be given exactly when the structure is corSymm or has `sd_ratio`")
  }
  r <- residual$sd_ratio
  if (!is.null(r) && (!is.numeric(r) || length(r) != length(times) ||
                      any(!is.finite(r)) || any(r <= 0) || abs(r[1L] - 1) > 1e-8)) {
    bad("`sd_ratio` must be positive, one per time, and 1 at the first time")
  }
  x <- residual$coef
  if (!is.numeric(x) || any(!is.finite(x))) bad("`coef` must be finite numbers")
  ok <- switch(cls,
    none    = length(x) == 0L,
    corAR1  = length(x) == 1L && abs(x[[1L]]) < 1,
    corCAR1 = length(x) == 1L && x[[1L]] > 0 && x[[1L]] < 1,
    corExp  = ,
    corGaus = length(x) %in% 1:2 && x[[1L]] > 0 &&
      (length(x) == 1L || (x[[2L]] >= 0 && x[[2L]] < 1)),
    corSymm = length(x) == length(times) * (length(times) - 1L) / 2 &&
      all(abs(x) < 1) && is_positive_definite(symm_matrix(x, length(times))))
  if (!ok) {
    bad(switch(cls,
      none    = "`coef` must be empty without a correlation",
      corAR1  = "Phi must lie in (-1, 1)",
      corCAR1 = "Phi must lie in (0, 1)",
      corExp  = ,
      corGaus = "range must be positive and any nugget in [0, 1)",
      corSymm = "the correlations must form a positive definite matrix over `times`"))
  }
  invisible(residual)
}

#' Fill a correlation matrix from `corSymm()`'s coefficient vector
#'
#' `nlme` lists the correlations row by row along the upper triangle --
#' (1,2), (1,3), ..., (2,3), ... -- which is the column-major order of the
#' lower triangle, so `lower.tri()` assignment places them directly.
#' @noRd
symm_matrix <- function(x, k) {
  m <- diag(k)
  m[lower.tri(m)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

# ---- evaluation at a schedule -----------------------------------------------

#' Residual covariance of one participant at visit times `t`
#'
#' \deqn{R_{ij} = \sigma^2_\epsilon \delta_i \delta_j \rho(t_i, t_j)}
#'
#' with \eqn{\delta} the visit-specific SD ratios (1 without `varIdent()`) and
#' \eqn{\rho} the correlation function. Each entry depends only on its own
#' pair of times, as the random-effects part of the covariance does, which is
#' what lets [effect_components()] slice a dropout stratum's covariance out of
#' the full schedule's rather than rebuild it.
#' @noRd
residual_cov <- function(params, t, context) {
  residual <- params$residual
  s2 <- params$sigma2_residual
  if (is.null(residual)) return(diag(s2, length(t)))
  pos <- if (residual_on_grid(residual)) grid_positions(t, residual$times, context)
  d <- if (is.null(residual$sd_ratio)) rep(1, length(t)) else residual$sd_ratio[pos]
  s2 * outer(d, d) * residual_cor(residual, t, pos, context)
}

#' The residual correlation matrix at visit times `t`
#'
#' Restates `nlme`'s definitions: `corAR1()` and `corCAR1()` are both
#' \eqn{\phi^{|t_i - t_j|}}, the former on integer times only; `corExp()` is
#' \eqn{(1 - n) e^{-d / r}} and `corGaus()` \eqn{(1 - n) e^{-(d / r)^2}} for
#' \eqn{d > 0}, with nugget \eqn{n} (0 without one) and range \eqn{r};
#' `corSymm()` is looked up by visit.
#' @noRd
residual_cor <- function(residual, t, pos, context) {
  x <- residual$coef
  d <- abs(outer(t, t, "-"))
  switch(residual$correlation,
    none    = diag(length(t)),
    corAR1  = {
      if (!all(is_whole(t))) {
        stop(sprintf(paste0(
          "%s: the residual correlation is corAR1(), which is defined on integer times only, ",
          "but the visit times include %s. Refit stage one with corCAR1(), its continuous-",
          "time counterpart, to plan visits between whole time units."),
          context, label_numeric(t[!is_whole(t)])), call. = FALSE)
      }
      x[["Phi"]]^round(d)
    },
    corCAR1 = x[["Phi"]]^d,
    corExp  = ,
    corGaus = {
      nugget <- if ("nugget" %in% names(x)) x[["nugget"]] else 0
      decay <- if (residual$correlation == "corExp") exp(-d / x[["range"]])
               else exp(-(d / x[["range"]])^2)
      out <- (1 - nugget) * decay
      diag(out) <- 1
      out
    },
    corSymm = symm_matrix(x, length(residual$times))[pos, pos, drop = FALSE])
}

#' Positions of planned visits in a structure's visit grid, or a refusal
#' @noRd
grid_positions <- function(t, times, context) {
  pos <- vapply(t, function(v) {
    hit <- which(abs(times - v) <= 1e-8 * max(1, abs(v)))
    if (length(hit)) hit[1L] else NA_integer_
  }, integer(1L))
  if (anyNA(pos)) {
    stop(sprintf(paste0(
      "%s: the residual structure has a parameter for each of the visit times %s, and ",
      "says nothing about time(s) %s. Plan visits among those times, or refit stage one ",
      "with a structure that is a function of time, such as corCAR1() or corExp()."),
      context, label_numeric(times), label_numeric(t[is.na(pos)])), call. = FALSE)
  }
  pos
}

# ---- bootstrap ---------------------------------------------------------------

#' `correlation` and `weights` arguments that refit the same structure
#'
#' For [make_refitter()], whose replicate frames carry the columns `time` and
#' `subject`. A correlation the user fixed stays fixed at the same value;
#' an estimated one is re-estimated on every replicate.
#' @noRd
residual_refit_args <- function(residual) {
  if (is.null(residual)) return(list())
  out <- list()
  if (!identical(residual$correlation, "none")) {
    out$correlation <- correlation_call(
      residual$correlation,
      value  = if (isTRUE(residual$fixed)) residual$coef,
      nugget = "nugget" %in% names(residual$coef),
      fixed  = isTRUE(residual$fixed),
      time = quote(time), subject = quote(subject), position = quote(time))
  }
  if (!is.null(residual$sd_ratio)) out$weights <- quote(nlme::varIdent(form = ~ 1 | time))
  out
}

# ---- printing ------------------------------------------------------------------

#' The residual lines of [print.slope_params()]
#' @noRd
print_residual <- function(x) {
  r <- x$residual
  cat_line(if (is.null(r$sd_ratio)) "residual variance" else
             sprintf("residual variance at time %s", time_labels(r$times[1L])),
           x$sigma2_residual)
  if (is.null(r)) return(invisible(x))
  cat_line("residual correlation", residual_label(r$correlation))
  if (r$correlation %in% c("corAR1", "corCAR1", "corExp", "corGaus")) {
    for (nm in names(r$coef)) {
      cat_line(paste0("correlation ", nm, if (isTRUE(r$fixed)) " (fixed)" else ""),
               r$coef[[nm]])
    }
  }
  if (identical(r$correlation, "corSymm")) {
    m <- symm_matrix(r$coef, length(r$times))
    lab <- time_labels(r$times)
    for (i in seq_along(lab)[-1L]) for (j in seq_len(i - 1L)) {
      cat_line(sprintf("correlation, times %s and %s", lab[j], lab[i]), m[i, j])
    }
  }
  if (!is.null(r$sd_ratio)) {
    lab <- time_labels(r$times)
    for (i in seq_along(lab)[-1L]) {
      cat_line(sprintf("residual SD ratio, time %s", lab[i]), r$sd_ratio[i])
    }
  }
  invisible(x)
}

#' The note every print method adds under a residual structure
#'
#' Shown with the parameters and with every stage-two result, as the covariate
#' note is and for the same reason: the treatment-effect variance is computed
#' for an analysis that models this residual structure, and nothing in the
#' printed sample size says so.
#' @noRd
residual_note <- function(residual) {
  if (is.null(residual)) return(invisible())
  text <- sprintf(paste0("residuals are %s%s%s. Plan for an analysis that models the same ",
                         "residual structure%s."),
                  if (identical(residual$correlation, "none")) "independent" else
                    paste("correlated,", residual_label(residual$correlation)),
                  if (is.null(residual$sd_ratio)) "" else ", with a variance for each visit",
                  if (!residual_on_grid(residual)) "" else
                    sprintf(", defined at times %s only", label_numeric(residual$times)),
                  if (!residual_on_grid(residual)) "" else
                    "; planned visits must be among those times")
  cat("\n", paste(strwrap(text, width = 72L, initial = "Note: ", prefix = "      "),
                  collapse = "\n"), "\n", sep = "")
  invisible()
}
