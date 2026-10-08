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
#' Accepts numeric 0/1, logical, haven_labelled, or a factor/character whose
#' levels are literally "0" and "1". Labelled factors such as
#' c("case", "control") are rejected rather than mapped by level order: which
#' level means "case" cannot be inferred, and guessing silently swaps the two
#' groups, so the fitted slope and every variance component would be taken from
#' the wrong one. The numeric path has always been strict about this; the
#' factor path used to trust alphabetical order, which made
#' `group = <"case"/"control" column>` return the healthy controls' slope
#' labelled as the cases'.
#'
#' @param meaning What "1" means for this column (e.g. `"case"` or
#'   `"treated"`), for the error messages below. The caller knows this; the
#'   coercion itself is agnostic to which of `slope_params()`'s arguments
#'   supplied `x`.
#' @noRd
coerce_binary <- function(x, name, context, meaning) {
  x <- unlabel(x)
  if (is.logical(x)) return(as.numeric(x))
  if (is.factor(x) || is.character(x)) {
    f <- if (is.factor(x)) x else factor(x)
    if (nlevels(f) != 2L) {
      stop(sprintf("%s: `%s` must have exactly two levels; got %d.",
                   context, name, nlevels(f)), call. = FALSE)
    }
    if (!identical(levels(f), c("0", "1"))) {
      stop(sprintf(paste0(
        "%s: `%s` is a %s with levels %s, so which level means \"%s\" cannot be ",
        "determined.\n  Recode it explicitly as 0/1, e.g.",
        "\n    data$%s <- as.integer(data$%s == \"%s\")\n  where 1 marks the %s group."),
        context, name, if (is.factor(x)) "factor" else "character vector",
        paste(sprintf("\"%s\"", levels(f)), collapse = " and "), meaning,
        name, name, levels(f)[2L], meaning), call. = FALSE)
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
                 context, name, meaning,
                 paste(format(u), collapse = "/")), call. = FALSE)
  }
  # all.equal() admits codes a hair off 0/1; snap them so the `== 1` tests
  # downstream see exactly what this check accepted.
  as.numeric(x == u[2L])
}

#' Look up a fixed effect by name, or stop
#'
#' The lookup itself -- including why an interaction has to be resolved rather
#' than spelled -- is [resolve_fixef_name()] in utils.R, shared with
#' [slope_se()]. This adds only the error for a name that is not there at all,
#' which the other caller answers with `NA` instead.
#' @noRd
fixef_term <- function(b, parts, context) {
  hit <- resolve_fixef_name(b, parts)
  if (is.na(hit)) {
    stop(sprintf("%s: fixed effect `%s` not found in the fitted model; have %s.",
                 context, paste(parts, collapse = ":"), paste(names(b), collapse = ", ")),
         call. = FALSE)
  }
  unname(b[[hit]])
}

#' The fixed-effect terms that sum to the slope, by comparator
#'
#' One list per comparator, each element a term spec as [resolve_fixef_name()]
#' expects: a single name, or the two parts of an interaction in model-formula
#' order. [slope_params()] sums each term's *value* -- via [fixef_term()] -- to
#' build the point estimate for `slope`; [slope_se()] sums the same terms'
#' variances and covariance to build its standard error. One mapping means the
#' two can never name a different set of coefficients for the same comparator,
#' which used to happen silently -- see the note in [slope_se()].
#' @noRd
slope_fixef_parts <- function(comparator) {
  switch(comparator,
    none    = list("sp_time"),
    treated = list("sp_time", "sp_placebo_time"),
    healthy = list("sp_time", c("sp_case", "sp_time"))
  )
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

#' Strip a `haven_labelled` wrapper, leaving anything else alone
#'
#' `haven::read_dta()` returns labelled columns, and this port's whole reason
#' for existing is that its data arrives from Stata. Every column the package
#' coerces has to shed the wrapper before its type can be tested, so the rule
#' is written once here rather than inline at each of the three columns
#' ([coerce_binary()], [coerce_time()] and the outcome in `slope_params()`).
#' Deliberately narrower than [coerce_time()], which also unwraps `Date` and
#' `POSIXct`: those mean something as a time and nothing as an outcome.
#' @noRd
unlabel <- function(x) {
  if (inherits(x, "haven_labelled")) as.numeric(x) else x
}

#' Coerce a time vector to numeric, unwrapping Date/POSIXct
#' @noRd
coerce_time <- function(x, context) {
  if (inherits(x, "Date")) return(as.numeric(x))
  if (inherits(x, "POSIXct")) return(as.numeric(x) / 86400)
  x <- unlabel(x)
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
extract_residual <- function(fit, level, context) {
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

#' The `nlme::lme()` control settings used by every fit
#'
#' Returns the [nlme::lmeControl()] settings that [slope_params()] passes to
#' every call to [nlme::lme()]: more iterations than `nlme`'s own defaults,
#' the `"optim"` optimiser, and tighter convergence tolerances. These were
#' chosen because the untightened defaults converged less precisely,
#' particularly for the two-block random-effects structure fitted when
#' `healthy` is supplied (see the `common_variance` note in
#' [slope_params()]).
#'
#' There is deliberately no argument to [slope_params()] for supplying a
#' different control object -- see "What these models do and do not include"
#' in `?slope_params` for why the model this package fits is fixed rather
#' than user-tunable. This function exists so the settings behind every fit
#' are inspectable and reproducible outside the package, not so they can be
#' overridden inside it.
#'
#' @return A list of control settings, as returned by [nlme::lmeControl()].
#'
#' @examples
#' slope_lme_control()
#'
#' @seealso [slope_params()], which uses this for every mixed-model fit.
#' @export
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

#' Fit the single-group or trial model
#'
#' A small helper for a reason that is not style. A formula written in a
#' function body captures that function's evaluation frame as its environment,
#' and the fitted object keeps it for life. Called directly from
#' [slope_params()], `lme()` would therefore pin `slope_params()`'s frame --
#' including the user's entire `data` argument, every column of it, used or not
#' -- inside `params$fit`, which is a contract field retained in every
#' `slope_sample_size` and `slope_power` result. Fitting from a small helper
#' frame instead drops that reference: on a 379 kB input frame it took a
#' serialized `slope_params` object from 523 kB to under 200 kB, and the saving
#' grows with the caller's data. The formula is therefore written in
#' [fixed_formula()], another small helper, and only passed through
#' `slope_params()`; it never closes over that frame. `fit_healthy_model()`
#' below has always had this property by accident of being a helper; this one
#' has it on purpose.
#' @noRd
fit_common_model <- function(dat, ctrl, fixed, resid) {
  eval(bquote(nlme::lme(.(fixed), random = ~ sp_time | sp_subject,
                        correlation = .(resid$correlation),
                        weights = .(resid$weights),
                        data = dat, method = "REML", control = ctrl)))
}

#' Append covariate terms to a fixed-effects formula
#'
#' With no covariates the formula is returned untouched, so an unadjusted fit
#' is exactly the model it always was. The new formula keeps the original's
#' environment -- the small helper frame -- for the reason given above.
#' `bquote()` at the call sites inlines the formula itself into `fit$call`
#' rather than the symbol `fixed`, so the call still reads as the model fitted.
#' @noRd
with_covariates <- function(f, cov_terms) {
  if (!length(cov_terms)) return(f)
  stats::as.formula(paste(deparse1(f), "+",
                          paste(cov_terms, collapse = " + ")),
                    env = environment(f))
}

#' @param resid The `correlation` and `weights` calls, from [residual_calls()].
#'   Under `healthy` the `weights` are always there: each group's own residual
#'   variance, crossed with the visit when visit-specific variances were asked
#'   for.
#' @noRd
fit_healthy_model <- function(dat, reduced, ctrl, fixed, resid) {
  comparator_block <- if (reduced) {
    nlme::pdIdent(~ sp_control - 1)
  } else {
    nlme::pdSymm(~ sp_control + sp_control_time - 1)
  }
  rand <- list(sp_subject = nlme::pdBlocked(list(
    nlme::pdSymm(~ sp_case + sp_case_time - 1),
    comparator_block)))
  fit_quietly(eval(bquote(
    nlme::lme(.(fixed),
              random      = rand,
              correlation = .(resid$correlation),
              weights     = .(resid$weights),
              data        = dat,
              method  = "REML",
              control = ctrl)
  )))
}

# ---- covariates -------------------------------------------------------------
#
# Covariates are handled in two passes, either side of `na.action`.
# [covariate_variables()] reads the *raw* variables the formula names, one value
# per visit, so that missing values can be filled within a participant and the
# rows still missing after that removed by `na.action` with everything else.
# [covariate_basis()] then expands the formula -- factors, interactions, bases
# such as poly() and splines::ns() -- on one row per *participant* of the data
# that survived. Expanding per visit instead, as this used to, had two faults:
# a data-dependent basis refused the `NA`s on follow-up rows that the
# baseline-only layout leaves (poly() errors on any `NA`), and its knots,
# orthogonalisation or scaling were weighted by how many visits each participant
# happened to have, so the same participants stored two ways gave two different
# adjustments.

#' Validate a covariate formula and read the raw variables it names
#'
#' Factors are dummy-coded later by [covariate_basis()]; here every variable is
#' returned as stored, one row per visit, by [stats::get_all_vars()] -- which
#' looks a name up in `data` first and then in the formula's environment, just
#' as [stats::model.frame()] would.
#' @return `NULL` without covariates, else a data frame of the raw variables.
#' @noRd
covariate_variables <- function(covariates, data, context) {
  if (is.null(covariates)) return(NULL)
  if (!inherits(covariates, "formula") || length(covariates) != 2L) {
    stop(sprintf("%s: `covariates` must be a one-sided formula, e.g. `~ age + sex`.",
                 context), call. = FALSE)
  }
  if ("." %in% all.vars(covariates)) {
    stop(sprintf(paste0("%s: `covariates` cannot use `.`; it would pull in every column, ",
                        "including the outcome, time and identifier. Name the covariates."),
                 context), call. = FALSE)
  }
  if (!is.null(attr(stats::terms(covariates), "offset"))) {
    stop(sprintf(paste0("%s: `covariates` cannot contain offset() terms; an offset is ",
                        "not an adjustment and would be dropped."), context), call. = FALSE)
  }
  raw <- tryCatch(stats::get_all_vars(covariates, data),
                  error = function(e) {
                    stop(sprintf("%s: could not evaluate `covariates`: %s",
                                 context, conditionMessage(e)), call. = FALSE)
                  })
  if (!ncol(raw)) {
    stop(sprintf("%s: `covariates` names no variables.", context), call. = FALSE)
  }
  raw
}

#' Carry each participant's recorded covariate value to their missing visits
#'
#' Long-format data often record a baseline covariate on the first visit only.
#' Left as `NA`, `na.action` would remove every follow-up visit -- leaving too
#' little repeated data to fit, or, if some participants had the value on every
#' row, fitting silently on those alone. A covariate has to be constant within a
#' participant anyway, so a value recorded on any visit is that participant's
#' value. A participant with no value recorded at all is still removed by
#' `na.action`, and one whose recorded values differ is left for
#' [check_baseline_covariates()] to refuse. Filled by position rather than by
#' value so factor, character and logical columns keep their type.
#' @noRd
fill_baseline_covariates <- function(raw, subject) {
  for (cc in names(raw)) {
    miss <- is.na(raw[[cc]])
    if (!any(miss) || all(miss)) next
    donor <- which(!miss)[match(subject, subject[!miss])]
    raw[[cc]][miss] <- raw[[cc]][donor[miss]]
  }
  raw
}

#' Refuse covariates that change during a participant's follow-up
#'
#' Time-varying covariates are refused for the same reason group membership is:
#' every column is read row by row, and a covariate that changes during follow-up
#' turns the time coefficient into something other than a slope. Checked on the
#' raw variables, so the message names the variable the user wrote; any
#' function of variables that are constant within a participant is constant too.
#' @noRd
check_baseline_covariates <- function(raw, subject, context) {
  for (cc in names(raw)) {
    x <- raw[[cc]]
    same <- if (is.numeric(x)) {
      tol <- covariate_tol(x)
      function(v) all(abs(v - v[1L]) <= tol)
    } else {
      x <- as.character(x)
      function(v) all(v == v[1L])
    }
    varies <- !tapply(x, subject, same)
    if (any(varies)) {
      stop(sprintf(paste0(
        "%s: covariates must be constant within a participant (baseline values), ",
        "but `%s` changes during follow-up for %d participant(s)."),
        context, cc, sum(varies)), call. = FALSE)
    }
  }
  invisible(raw)
}

#' Tolerance for comparing covariate values
#' @noRd
covariate_tol <- function(x) {
  sqrt(.Machine$double.eps) * max(1, abs(x), na.rm = TRUE)
}

#' Expand the covariate formula per participant, centre it, and spread it back
#' over visits as numeric columns `sp_cov_1`, ...
#'
#' Factors are dummy-coded by [stats::model.matrix()] with the usual treatment
#' contrasts and the intercept column dropped. The internal names keep the
#' by-name extraction of the slope terms untouched: no covariate can collide
#' with `sp_time` or `sp_case`.
#'
#' A column that does not vary *between* participants -- the dummy for a factor
#' level absent from the data used (an unused level, or one whose rows were all
#' removed as missing), or one a bootstrap or jackknife resample happened to
#' leave out -- carries nothing the intercept does not, and would make the
#' fixed-effects design singular. It is dropped, which is the fit `droplevels()`
#' would have given.
#'
#' Each kept column is centred at its mean over the participants flagged by
#' `reference`, so the fitted slope is the slope of a participant with those
#' participants' average covariate values, whatever the visit pattern. Without
#' centring it would be the slope at covariate = 0 (age zero, say) whenever
#' covariate-by-time terms are included. [slope_params()] passes the cases under
#' `healthy` -- the population the planned trial will enrol -- and everyone
#' otherwise.
#'
#' @param raw The raw variables, one row per visit of the data used, filled and
#'   checked by the two helpers above.
#' @param reference Logical, one per visit: whose mean to centre on.
#' @return A list: `X`, the kept columns, one row per visit; and `labels`, their
#'   `model.matrix()` names, keyed by internal name.
#' @noRd
covariate_basis <- function(covariates, raw, subject, reference, context) {
  first <- !duplicated(subject)
  # The model always has an intercept, so the covariate columns must be coded
  # against one: `~ 0 + sex` would otherwise give a full set of dummies,
  # collinear with the model's own intercept.
  tt <- stats::terms(covariates)
  attr(tt, "intercept") <- 1L
  P <- tryCatch(
    stats::model.matrix(tt, stats::model.frame(tt, raw[first, , drop = FALSE])),
    error = function(e) {
      stop(sprintf("%s: could not evaluate `covariates`: %s",
                   context, conditionMessage(e)), call. = FALSE)
    })
  P <- P[, colnames(P) != "(Intercept)", drop = FALSE]
  if (!ncol(P)) {
    stop(sprintf("%s: `covariates` produced no columns.", context), call. = FALSE)
  }
  labels <- stats::setNames(colnames(P), paste0("sp_cov_", seq_len(ncol(P))))
  colnames(P) <- names(labels)
  varies <- apply(P, 2L, function(x) any(abs(x - x[1L]) > covariate_tol(x)))
  P <- P[, varies, drop = FALSE]
  ref <- reference[first]
  P <- sweep(P, 2L, colMeans(P[ref, , drop = FALSE]))
  list(X = P[match(subject, subject[first]), , drop = FALSE],
       labels = labels[varies])
}

#' Refuse a fixed-effects design whose covariate terms cannot be separated
#'
#' A covariate that is the group indicator under another name (controls
#' recruited from a site of their own, adjusted for site), or a function of
#' other covariates, makes the fixed-effects design rank deficient. `nlme`
#' then fails with "Singularity in backsolve", and under `healthy` that
#' failure used to be taken for non-convergence: the fit fell back to the
#' reduced structure, failed the same way, and was reported as a model that
#' "did not converge". The problem is the design, not the optimiser, so it is
#' named here before anything is fitted. Pivoted QR places the aliased columns
#' last, and the covariate terms follow the comparator's own, so it is a
#' covariate term that is named.
#' @noRd
check_covariate_rank <- function(fixed, dat, labels, comparator, context) {
  M <- stats::model.matrix(fixed, dat)
  q <- qr(M)
  if (q$rank == ncol(M)) return(invisible())
  aliased <- colnames(M)[q$pivot[-seq_len(q$rank)]]
  readable <- vapply(strsplit(aliased, ":", fixed = TRUE), function(p) {
    p[p %in% names(labels)] <- labels[p[p %in% names(labels)]]
    p[p == "sp_time"] <- "time"
    p[p == "sp_case"] <- comparator
    paste(p, collapse = ":")
  }, character(1L))
  stop(sprintf(paste0(
    "%s: the covariates cannot be separated from the rest of the model: %s %s ",
    "an exact linear combination of the other terms%s. Remove the covariate%s ",
    "concerned, or any covariate that duplicates another."),
    context, paste(sprintf("`%s`", readable), collapse = ", "),
    if (length(readable) == 1L) "is" else "are",
    if (comparator == "none") "" else
      sprintf(" -- a covariate that is fixed by `%s` cannot be adjusted for", comparator),
    if (length(readable) == 1L) "" else "s"), call. = FALSE)
}

#' The fixed-effects formula for a comparator, with any covariate terms
#'
#' Built here rather than in [slope_params()] for the reason given at
#' [fit_common_model()]: the formula's environment is this small frame, not one
#' holding the user's data.
#' @noRd
fixed_formula <- function(comparator, cov_terms = character()) {
  base <- switch(comparator,
    none    = sp_y ~ sp_time,
    treated = sp_y ~ sp_time + sp_placebo_time,
    healthy = sp_y ~ sp_case * sp_time)
  with_covariates(base, cov_terms)
}

#' How a fit was adjusted, for print methods
#'
#' Shown with the parameters, and again with every stage-two result, because
#' the variance components are only right for a trial analysed with the same
#' adjustment, and nothing else in the printed numbers says they are adjusted.
#' @noRd
covariate_note <- function(covariates) {
  if (is.null(covariates)) return(invisible())
  text <- sprintf("adjusted for baseline covariates %s%s. Plan for an analysis with the same adjustment.",
                  paste(covariates$columns, collapse = ", "),
                  if (!isTRUE(covariates$time)) "" else
                    if (length(covariates$columns) == 1L) " and its interaction with time"
                    else ", and their interactions with time")
  cat("\n", paste(strwrap(text, width = 72L, initial = "Note: ", prefix = "      "),
                  collapse = "\n"), "\n", sep = "")
  invisible()
}

# ---- main entry point -------------------------------------------------------

#' Estimate slope and variance parameters from longitudinal data
#'
#' Stage one of the two-stage sample-size method of Nash et al. (2021). Fits a
#' linear mixed model to previously collected longitudinal data and extracts the
#' slope(s) and the between- and within-subject variance components that
#' [slope_sample_size()] and [slope_power()] need.
#'
#' @param formula A two-sided formula `outcome ~ time | subject`. The time term
#'   may be an expression, e.g. `sdmt ~ I(as.numeric(vdate) / 365) | id` to work
#'   in years when visits are recorded as dates. Choose units on the order of the
#'   study duration: fitting on a badly scaled axis (days over several years, say)
#'   leaves the random-slope variance near zero and the REML optimiser converges
#'   less precisely. There is no `scale()` argument as there is in Stata --
#'   express `visits` in stage two in whatever units are used here.
#' @param data A data frame in long format, one row per measurement.
#' @param comparator What, if anything, the untreated group in `data` is
#'   compared with: `"none"` (the default), `"healthy"` for observational data
#'   with healthy controls, or `"treated"` for data from a previous trial with
#'   a treated arm. See Details.
#' @param group Bare column name identifying the two groups, required unless
#'   `comparator = "none"` and refused with it. Under `"healthy"`, coded `1`
#'   for cases (subjects with the disease) and `0` for healthy controls; under
#'   `"treated"`, `1` for the treated arm and `0` for the control arm. Looked
#'   up in `data` first, then in the calling environment.
#' @param origin `"subject"` (default) shifts each subject's time so their first
#'   visit is time zero, reproducing the Stata command's behaviour and ensuring
#'   the random intercept is estimated at baseline. `"none"` leaves time as
#'   supplied.
#' @param common_variance Controls the random-effects structure for healthy
#'   controls. `NULL` (default) fits the full model and falls back automatically
#'   if it fails to converge; `TRUE` forces the reduced structure (a random
#'   intercept only, equivalent to the Stata `nocontvar` option); `FALSE` forces
#'   the full structure and errors on failure. Ignored unless
#'   `comparator = "healthy"`.
#' @param na.action Applied to the assembled model frame. Defaults to
#'   [stats::na.omit()].
#' @param covariates Optional one-sided formula of baseline covariates to adjust
#'   for, e.g. `~ age + sex`. Anything [stats::model.matrix()] accepts is
#'   allowed -- factors, interactions (`age * sex`), transformations
#'   (`log(age)`, `I(age^2)`) and bases (`poly(age, 2)`, `splines::ns(age, 3)`)
#'   -- as long as every resulting column is constant within a participant.
#'   Factors and character columns use treatment contrasts, and are always
#'   coded against an intercept, so `~ 0 + sex` is the same as `~ sex`.
#'   `offset()` terms and `.` are refused. A value recorded on any of a
#'   participant's visits is used for all of them, so a covariate entered on
#'   the baseline row only is enough; participants with no value at all are
#'   removed by `na.action`. Factors, interactions and bases are expanded on
#'   one row per participant, so a basis's knots or scaling do not depend on
#'   how many visits each participant has. A column that does not vary between
#'   participants (e.g. an unused factor level) is dropped. Every column is
#'   centred at its mean over participants -- over the cases only when
#'   `comparator = "healthy"` -- so `slope` is the slope of a participant with those
#'   average covariate values. The returned variance components are then the
#'   *adjusted* ones, and the object records the adjustment in `$covariates`.
#'   A covariate that is an exact linear combination of the group indicator or
#'   of other covariates is refused. See "Covariate adjustment" below.
#' @param covariate_time If `TRUE` (default), also include each covariate's
#'   interaction with time, so covariates can explain differences in slope as
#'   well as in baseline. This is what reduces `sigma2_slope`, and so the
#'   sample size. `FALSE` adjusts the intercept only. Ignored without
#'   `covariates`.
#' @param correlation Optional residual correlation structure, as an `nlme`
#'   `corStruct` exactly as it would be passed to [nlme::lme()]:
#'   [nlme::corAR1()], [nlme::corCAR1()], [nlme::corExp()],
#'   [nlme::corGaus()] or [nlme::corSymm()]. `NULL` (default) keeps the
#'   independent residuals of Nash et al. (2021). The correlation is always
#'   within a participant and a function of the time term of `formula`, so the
#'   `form` is set for you: leave it out, or write the equivalent
#'   `form = ~ time | subject` with the names used in `formula`. A starting
#'   value, a nugget (`corExp(nugget = TRUE)`) and `fixed = TRUE` are honoured
#'   as in `nlme`. See "Residual structures" below.
#' @param weights Optional visit-specific residual variances, written
#'   `nlme::varIdent(form = ~ 1 | time)` with the time term of `formula` as the
#'   stratum. Together with `correlation = nlme::corSymm()` this is the fully
#'   unstructured residual covariance. No other variance function is
#'   supported.
#'
#' @details
#' Three scenarios are supported, matching paper section 2.3:
#'
#' * `comparator = "none"`: a single group of untreated subjects with
#'   the disease. The target treatment effect will be measured toward a slope of
#'   zero.
#' * `comparator = "healthy"`: observational data containing both cases and
#'   healthy controls.
#'   The target effect will be measured toward the healthy-control slope.
#' * `comparator = "treated"`: data from a previous trial. The observed
#'   treatment effect is
#'   available via `target = "observed"` in [slope_sample_size()] and
#'   [slope_power()].
#'
#' The returned variance components are always those of the **untreated / case**
#' group. Under `comparator = "healthy"` the controls contribute only their slope;
#' their variance components are estimated and discarded, per paper section 2.3.
#'
#' Note that for the `"healthy"` scenario without covariates the model factorises
#' exactly into two independent fits, one per group: the fixed effects
#' `y ~ case * time` span the same column space as separate per-group intercepts
#' and slopes, the random-effects blocks are independent, and the residual
#' variances are separate. Consequently `common_variance` cannot change the
#' estimates returned for the cases -- it only affects how many nuisance
#' parameters are estimated for the controls, and therefore whether the fit
#' converges at all. With `covariates` this no longer holds; see "Covariate
#' adjustment".
#'
#' @section The models fitted:
#'
#' Write \eqn{y_{ij}}{y[ij]} for the outcome of participant \eqn{i}{i} at time
#' \eqn{t_{ij}}{t[ij]}, measured from that participant's own first visit under
#' the default `origin = "subject"`. All three scenarios share one
#' participant-level structure --- a random intercept and a random slope, with
#' independent residuals:
#'
#' \deqn{y_{ij} = \mu(t_{ij}) + a_i + b_i t_{ij} + \epsilon_{ij},}{
#'       y[ij] = mu(t[ij]) + a[i] + b[i] * t[ij] + e[ij],}
#' \deqn{\left(a_i, b_i\right) \sim N(0, G), \qquad
#'       \epsilon_{ij} \sim N(0, \sigma^2_\epsilon),}{
#'       (a[i], b[i]) ~ N(0, G),   e[ij] ~ N(0, sigma2_residual),}
#'
#' where \eqn{G}{G} is an unstructured 2 by 2 matrix with diagonal
#' \eqn{\sigma^2_a}{sigma2_intercept}, \eqn{\sigma^2_b}{sigma2_slope} and
#' off-diagonal \eqn{\sigma_{ab}}{sigma_cov}. By default residuals are
#' independent across visits and across participants; `correlation` and
#' `weights` replace that with a structured residual covariance within each
#' participant (see "Residual structures"). Only the mean \eqn{\mu}{mu} and the
#' number of variance parameters differ between the scenarios.
#'
#' \describe{
#'   \item{`comparator = "none"`}{
#'     \deqn{\mu(t) = \beta_0 + \beta_1 t}{mu(t) = b0 + b1 * t}
#'     with one \eqn{G}{G} and one \eqn{\sigma^2_\epsilon}{sigma2_residual}. The
#'     returned slope is \eqn{\beta_1}{b1}. Equivalent to
#'     `nlme::lme(y ~ t, random = ~ t | id)`.}
#'   \item{`comparator = "healthy", group = g`, with \eqn{g_i = 1}{g[i] = 1} for cases}{
#'     \deqn{\mu(t) = \beta_0 + \beta_g g_i + (\beta_1 + \beta_{1g} g_i) t}{
#'           mu(t) = b0 + bg * g[i] + (b1 + b1g * g[i]) * t}
#'     and, in addition, a **separate** \eqn{G}{G} and a **separate**
#'     \eqn{\sigma^2_\epsilon}{sigma2_residual} for each group. The slope of the
#'     cases is \eqn{\beta_1 + \beta_{1g}}{b1 + b1g} and that of the controls
#'     \eqn{\beta_1}{b1}; the variance components returned are the cases'.
#'     Without covariates the model factorises into two independent per-group
#'     fits (see above).}
#'   \item{`comparator = "treated", group = z`, with \eqn{z_i = 1}{z[i] = 1} for the treated arm}{
#'     \deqn{\mu(t) = \beta_0 + \beta_1 t + \beta_p (1 - z_i) t}{
#'           mu(t) = b0 + b1 * t + bp * (1 - z[i]) * t}
#'     with one \eqn{G}{G} and one
#'     \eqn{\sigma^2_\epsilon}{sigma2_residual} shared by both arms. Note the
#'     single intercept and the absence of a main effect of \eqn{z}{z}:
#'     randomisation makes the expected baseline equal in the two arms, so
#'     baseline is modelled as a correlated outcome rather than adjusted for.
#'     The slope of the treated arm is \eqn{\beta_1}{b1} and that of the control
#'     arm \eqn{\beta_1 + \beta_p}{b1 + bp}.}
#' }
#'
#' Stage two then assumes the planned trial will be analysed with the matching
#' model, \eqn{\mu(t) = \beta_0 + \beta_1 t + \beta_2 g_i t}{mu(t) = b0 + b1 * t
#' + b2 * g[i] * t}, in which the treatment effect \eqn{\beta_2}{b2} is the
#' difference in slopes and the arms again share an intercept.
#'
#' @section What these models do and do not include:
#'
#' The design matrices above are fixed by the method and are not what you would
#' write by hand in `lme4`. Both comparator models depart from the obvious
#' `y ~ group * time + (time | id)`: the `"treated"` model drops the group main
#' effect, and the `"healthy"` model gives each group its own residual variance,
#' which `lme4` cannot fit at all. Included, therefore:
#'
#' * one continuous, approximately Gaussian outcome;
#' * a mean that is linear in time, plus optional baseline covariates (see
#'   "Covariate adjustment");
#' * exactly one grouping level, the participant, whose random intercept and
#'   random slope have an unstructured covariance;
#' * residuals independent between participants, and within a participant
#'   either independent with a variance that is constant within a group (the
#'   default) or structured through `correlation` and `weights`;
#' * at most two groups, distinguished only by their slope (and, for `"healthy"`,
#'   by their variance components).
#'
#' Not included, and not obtainable by any argument to this function:
#'
#' * **time-varying covariates.** Baseline covariates can be entered through
#'   `covariates` (see "Covariate adjustment"); covariates that change during
#'   follow-up cannot. `formula` itself still takes a single time term.
#' * **baseline as a covariate.** Baseline is part of the outcome vector, under
#'   a common intercept for both arms (paper section 2.1). This is not the
#'   ANCOVA-style adjustment used by many trial analyses, and the two give
#'   different standard errors.
#' * **further levels of clustering.** Visits within participants within sites,
#'   clinics, families or therapists cannot be represented: the random effects
#'   have one grouping factor, and stage two treats participants as independent.
#'   For a design with meaningful site-level variation the sample size returned
#'   here will be too small.
#' * **non-linear trajectories** --- quadratic time, splines, change points ---
#'   and any estimand that is not a difference in slopes.
#' * **non-Gaussian outcomes**: binary, ordinal, count or time-to-event.
#' * **residual structures other than those listed under "Residual
#'   structures"**: compound symmetry (`corCompSymm()`, which the random
#'   intercept already is), ARMA, spatial structures other than the exponential
#'   and Gaussian, and variance functions other than `varIdent()` by visit.
#' * **more than two arms, unequal allocation, or cluster-randomised, crossover
#'   and stepped-wedge designs.** Stage two assumes two equal parallel arms.
#'
#' @section Covariate adjustment:
#'
#' With `covariates = ~ x1 + x2`, the mean in each scenario above gains
#' \eqn{\gamma_1 x_{1i} + \gamma_2 x_{2i}}{g1 * x1[i] + g2 * x2[i]} and, when
#' `covariate_time = TRUE`, \eqn{(\delta_1 x_{1i} + \delta_2 x_{2i}) t}{(d1 *
#' x1[i] + d2 * x2[i]) * t}, with the covariates centred over participants. The
#' coefficients are shared by both groups, the random effects and residuals are
#' unchanged, and the slope and variance components are extracted exactly as
#' before -- they are now conditional on the covariates.
#'
#' Under `comparator = "healthy"` the covariates are centred at the mean of the **cases**, not
#' of everyone, because the planned trial enrols cases: `slope` is the slope of
#' a case with the cases' average covariate values, and `slope_comparator` that
#' of a healthy control with the same values. Their difference does not depend
#' on where the covariates are centred. Under `comparator = "treated"` both arms come from the
#' trial population, so everyone is used.
#'
#' Two consequences. First, stage two then assumes the planned trial will be
#' analysed with the same adjustment: using adjusted variance components to
#' plan an unadjusted analysis overstates power. The adjustment is recorded in
#' the returned object's `$covariates` -- a list of the `columns` adjusted for
#' and whether their interactions with `time` were included, or `NULL` -- and
#' is printed with the parameters and with every stage-two result. Second,
#' because the covariate coefficients are shared across groups, the `"healthy"`
#' model no longer factorises into two independent fits, so `common_variance`
#' can now move the cases' estimates slightly.
#'
#' Baseline covariates are constant within a participant, so they reduce the
#' between-participant variances (`sigma2_intercept`, and `sigma2_slope` via the
#' time interactions) rather than `sigma2_residual`.
#'
#' @section Residual structures:
#'
#' With `correlation` and `weights` the residual term above becomes
#' \eqn{\epsilon_i \sim N(0, R_i)}{e[i] ~ N(0, R[i])} with
#' \deqn{(R_i)_{jk} = \sigma^2_\epsilon \, \delta_j \delta_k \,
#'       \rho(t_{ij}, t_{ik}),}{
#'       R[i][j, k] = sigma2_residual * delta[j] * delta[k] * rho(t[ij], t[ik]),}
#' the notation of `nlme`: \eqn{\rho}{rho} is the correlation and
#' \eqn{\delta}{delta} the residual SD ratios of `varIdent()`, 1 without it.
#' The structures, and the arguments that give them:
#'
#' | Structure | Argument | \eqn{\rho(s, t)}{rho(s, t)}, \eqn{d = |s - t|}{d = |s - t|} | Stata `mixed` |
#' |---|---|---|---|
#' | independent (default) | `correlation = NULL` | 0 | `residuals(independent)` |
#' | AR(1) | `nlme::corAR1()` | \eqn{\phi^d}{Phi^d}, integer times only | `residuals(ar 1, t())` |
#' | continuous-time AR(1) | `nlme::corCAR1()` | \eqn{\phi^d}{Phi^d} | `residuals(exponential, t())` |
#' | exponential | `nlme::corExp()` | \eqn{(1 - n) e^{-d/r}}{(1 - n) exp(-d / r)} | |
#' | Gaussian | `nlme::corGaus()` | \eqn{(1 - n) e^{-(d/r)^2}}{(1 - n) exp(-(d / r)^2)} | |
#' | unstructured correlation | `nlme::corSymm()` | one per pair of visit times | |
#' | unstructured covariance | `corSymm()` plus `weights` | one per pair, and a variance per time | `residuals(unstructured, t())` |
#'
#' Here \eqn{r}{r} is the range and \eqn{n}{n} the nugget, 0 unless
#' `nugget = TRUE`; a nugget splits the residual into measurement error and a
#' serially correlated part. `corAR1()` and `corCAR1()` agree wherever both
#' are defined; the first is `nlme`'s discrete-time process and so needs
#' integer times, both in the data and in any schedule planned from it. The
#' unstructured covariance adds `weights = nlme::varIdent(form = ~ 1 | time)`,
#' a variance for each visit time, to `corSymm()`.
#'
#' The fitted values are returned in `$residual` (see the "Value" section) and
#' used by every stage-two function, which then plans for a trial analysed
#' with the same residual structure. Three consequences:
#'
#' * `corSymm()` and `varIdent()` have a parameter for every distinct visit
#'   time, so they need the participants to share a visit schedule -- record
#'   time as the scheduled visit, not the date it took place -- and a planned
#'   schedule can only use those times. A structure that is a function of
#'   time (`corCAR1()`, `corExp()`, `corGaus()`) can price any schedule.
#' * Under `corSymm()` the random-effects variances and the residual variance
#'   are not separately identified: an unstructured correlation can absorb any
#'   covariance the random effects imply. The covariance they imply together
#'   at the visit times is identified, and that is the only thing stage two
#'   uses, but the individual components printed are one point on a ridge of
#'   equally good fits, and [slope_var_floor()] is refused.
#' * Under `comparator = "healthy"` the correlation parameters are shared by cases and
#'   controls -- `nlme` estimates one correlation structure per model -- while
#'   the residual variances, and with `varIdent()` the per-visit variances, stay
#'   separate per group. So the model no longer factorises into two per-group
#'   fits, and `common_variance` can move the cases' estimates slightly, as it
#'   can with covariates.
#'
#' @return An object of class `"slope_params"`. Its `$residual` component is
#'   `NULL` for independent residuals, or else a list: `correlation`, the
#'   `nlme` class name (`"none"` with `varIdent()` alone); `coef`, the fitted
#'   correlation parameters on their natural scale (`Phi`; `range` and any
#'   `nugget`; or `corSymm()`'s correlations in `nlme`'s order); `fixed`;
#'   `times`, the visit times a `corSymm()` or `varIdent()` structure is
#'   defined at, else `NULL`; and `sd_ratio`, the residual SD at each of those
#'   times relative to the first, else `NULL`. When `sd_ratio` is set,
#'   `sigma2_residual` is the residual variance at the first time.
#'
#'   Its `$fit` component is the
#'   fitted `"lme"` object itself, useful for `nlme`'s diagnostic plots (e.g.
#'   `plot(fit)`, `qqnorm(fit)`) -- note that these must reference the
#'   internal column names (`sp_y`, `sp_time`, `sp_subject`, ...) described
#'   above, not the originals from `data`.
#'
#' @examples
#' # No comparator: a single group of untreated subjects.
#' # Four of the two hundred participants of `slpower1`, kept small so the
#' # example runs quickly -- see `slpower1` for the paper's fit on the full data.
#' df <- slpower1[slpower1$id %in% 1:4, ]
#' slope_params(sdmt ~ visit | id, data = df)
#'
#' # comparator = "healthy": two cases and two healthy controls, a subset of
#' # `slpower2`.
#' # Visits are recorded as calendar dates there, so the time term converts
#' # them to years.
#' df2 <- slpower2[slpower2$id %in% c(1, 2, 251, 252), ]
#' slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, data = df2,
#'              comparator = "healthy", group = case)
#'
#' # comparator = "treated": data from a completed trial, a subset of
#' # `slpower3`. Fitting the
#' # random-effects structure shared by both arms needs more than a couple of
#' # subjects per arm to converge, so this excerpt keeps six per arm.
#' df3 <- slpower3[slpower3$id %in% c(1:6, 76:81), ]
#' slope_params(sdmt ~ visit | id, data = df3, comparator = "treated", group = treat)
#'
#' # Adjusting for a baseline covariate (simulated here for illustration).
#' df4 <- slpower1[slpower1$id %in% 1:20, ]
#' set.seed(1)
#' df4$age <- ave(df4$id, df4$id, FUN = function(i) round(runif(1, 30, 60)))
#' slope_params(sdmt ~ visit | id, data = df4, covariates = ~ age)
#'
#' # Serially correlated residuals: a continuous-time AR(1), on the first 50
#' # participants of `slpower1`.
#' df5 <- slpower1[slpower1$id %in% 1:50, ]
#' slope_params(sdmt ~ visit | id, data = df5, correlation = nlme::corCAR1())
#'
#' # A fully unstructured residual covariance, as `residuals(unstructured)`
#' # in Stata's `mixed`.
#' slope_params(sdmt ~ visit | id, data = df5,
#'              correlation = nlme::corSymm(),
#'              weights = nlme::varIdent(form = ~ 1 | visit))
#'
#' @references
#' Nash, S., Morgan, K. E., Frost, C. and Mulick, A. (2021). Power and
#' sample-size calculations for trials that compare slopes over time:
#' Introducing the slopepower command. \emph{The Stata Journal} 21(3): 575--601.
#' \doi{10.1177/1536867X211045512}
#'
#' @seealso [slope_params_manual()] to supply parameters directly,
#'   [slope_sample_size()] and [slope_power()] for stage two,
#'   [slope_params_boot()] for an interval around the fitted slope.
#' @export
slope_params <- function(formula, data,
                         comparator = c("none", "healthy", "treated"),
                         group = NULL,
                         origin = c("subject", "none"),
                         common_variance = NULL,
                         na.action = stats::na.omit,
                         covariates = NULL, covariate_time = TRUE,
                         correlation = NULL, weights = NULL) {
  context <- "slope_params()"
  cl <- match.call()
  # `group` is an argument of this call, so a symbol in it that is
  # not a column belongs to whoever wrote the call -- not to whoever built the
  # formula, which can be a different frame entirely (a wrapper passing its
  # own group vector alongside a formula made at top level).
  caller <- parent.frame()
  comparator <- match.arg(comparator)
  origin <- match.arg(origin)

  data <- tryCatch(as.data.frame(data), error = function(e) {
    stop(sprintf("%s: `data` must be a data frame.", context), call. = FALSE)
  })

  parts <- parse_slope_formula(formula, context)
  if (is.null(parts$subject)) {
    stop(sprintf("%s: `formula` must name the subject identifier, e.g. `%s ~ %s | id`.",
                 context, deparse(parts$outcome), deparse(parts$time)),
         call. = FALSE)
  }

  # Checked before anything is evaluated or fitted: these are about what was
  # typed, not about the data.
  spec <- residual_spec(correlation, weights, parts, context)

  env <- environment(formula) %||% parent.frame()

  y       <- eval_column(parts$outcome, data, env, context, "outcome")
  tim     <- coerce_time(eval_column(parts$time, data, env, context, "time"), context)
  subject <- eval_column(parts$subject, data, env, context, "subject")

  y <- unlabel(y)
  if (!is.numeric(y)) {
    stop(sprintf("%s: the outcome must be numeric.", context), call. = FALSE)
  }

  group_expr <- substitute(group)
  if (identical(comparator, "none") && !is.null(group_expr)) {
    stop(sprintf(paste0(
      "%s: `group` is supplied but comparator = \"none\".\n  Set comparator = ",
      "\"healthy\" (1 = case, 0 = healthy control) or \"treated\" (1 = treated, ",
      "0 = control) to say what it identifies."), context), call. = FALSE)
  }
  if (!identical(comparator, "none") && is.null(group_expr)) {
    stop(sprintf("%s: comparator = \"%s\" needs `group`, the column identifying %s.",
                 context, comparator,
                 if (comparator == "healthy") "cases (1) and healthy controls (0)"
                 else "the treated (1) and control (0) arms"), call. = FALSE)
  }

  grp <- NULL
  if (comparator != "none") {
    graw  <- eval_column(group_expr, data, caller, context, "group")
    grp   <- coerce_binary(graw, "group", context,
                           meaning = if (comparator == "healthy") "case" else "treated")
    if (length(grp) != length(y)) {
      stop(sprintf("%s: `group` has length %d but the data have %d rows.",
                   context, length(grp), length(y)), call. = FALSE)
    }
  }

  common_variance <- warn_unused_arg(
    common_variance, !is.null(common_variance) && comparator != "healthy", NULL,
    "%s: `common_variance` applies only when comparator = \"healthy\"; ignoring it.",
    context)

  n <- length(y)
  if (length(tim) != n || length(subject) != n) {
    stop(sprintf("%s: outcome, time and subject must be the same length.", context),
         call. = FALSE)
  }

  dat <- data.frame(sp_y = y, sp_time = tim,
                    sp_subject = factor(as.character(subject)),
                    stringsAsFactors = FALSE)
  if (!is.null(grp)) dat$sp_case <- grp
  # The raw covariates stay out of `dat`, which holds only the internal columns;
  # `sp_row` finds a surviving row's covariates again after `na.action`, and
  # `sp_cov_ok` is NA where any is still missing after filling, so `na.action`
  # removes those rows with the rest.
  raw <- covariate_variables(covariates, data, context)
  if (!is.null(raw)) {
    if (nrow(raw) != n) {
      stop(sprintf("%s: `covariates` gave %d rows but the data have %d.",
                   context, nrow(raw), n), call. = FALSE)
    }
    raw <- fill_baseline_covariates(raw, dat$sp_subject)
    dat$sp_row    <- seq_len(n)
    dat$sp_cov_ok <- ifelse(stats::complete.cases(raw), 1, NA_real_)
  }

  dat <- na.action(dat)
  if (nrow(dat) < 3L) {
    stop(sprintf("%s: fewer than 3 usable observations after removing missing values.",
                 context), call. = FALSE)
  }
  dat$sp_subject <- droplevels(dat$sp_subject)

  # A random-slope model needs subjects with repeat visits to identify the
  # slope variance at all; a subject seen once contributes nothing to it. Row
  # count alone does not catch this -- 3 rows can be 2 subjects, one of them
  # seen only once -- so check participants directly. Two is the bare
  # mathematical minimum (with one, there is no between-subject variance to
  # estimate); it is not a claim that two is *enough* for a trustworthy fit.
  n_repeat <- sum(table(dat$sp_subject) >= 2L)
  if (n_repeat < 2L) {
    stop(sprintf(paste0(
      "%s: too little repeated-measures data to fit a random-slope model: only %d ",
      "of %d participant(s) have more than one visit. At least 2 participants with ",
      "repeat visits are needed to identify the slope variance."),
      context, n_repeat, nlevels(dat$sp_subject)), call. = FALSE)
  }

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
  # On the times as fitted, after the origin shift: a visit grid is a set of
  # times since each participant's first visit.
  grid <- if (!is.null(spec)) check_residual_data(spec, dat$sp_time, dat$sp_subject, context)

  if (comparator != "none") {
    if (length(unique(dat$sp_case)) != 2L) {
      stop(sprintf("%s: `%s` must contain both groups after removing missing values.",
                   context, comparator), call. = FALSE)
    }
    # Group membership also has to be a property of the participant rather than
    # of the visit. Every model below reads `sp_case` row by row, so a
    # participant whose indicator changes
    # part-way through their follow-up is fitted as a case for some visits and a
    # control for others -- loading on both random-effects blocks at once for
    # `healthy`, and switching arms mid-trial for `treated`. Nothing about the
    # fit looks wrong afterwards; it just answers a different question. The
    # participant is also counted in both groups by the size check below, so a
    # dataset made entirely of such rows could pass that too.
    split_subject <- tapply(dat$sp_case, dat$sp_subject,
                            function(g) any(g != g[1L]))
    if (any(split_subject)) {
      offenders <- names(split_subject)[which(split_subject)]
      stop(sprintf(paste0(
        "%s: `%s` must be constant within a participant, but it changes during ",
        "follow-up for %d of %d participant(s) (%s%s). Group membership is a ",
        "property of the participant, not of the visit."),
        context, comparator, length(offenders), nlevels(dat$sp_subject),
        paste(sQuote(offenders[seq_len(min(5L, length(offenders)))]), collapse = ", "),
        if (length(offenders) > 5L) ", ..." else ""), call. = FALSE)
    }
    # Presence of both groups is not enough: for `healthy`/`treated` each group
    # gets its own variance components (the `healthy` model factorises into two
    # fully independent per-group fits -- see the note above), so a group too
    # small to have a between-subject variance of its own returns a fit that
    # looks exactly like a well-supported one instead of failing. Two per group
    # is the same bare mathematical minimum as the repeat-visit check above,
    # applied per group instead of overall; below it, `nlme` still "converges"
    # and returns an unstable, misleadingly precise-looking number.
    group_n <- tapply(dat$sp_subject, dat$sp_case, function(s) length(unique(s)))
    if (any(group_n < 2L)) {
      stop(sprintf(paste0(
        "%s: each level of `%s` needs at least 2 participants to identify its own ",
        "variance components; got %d (coded 0) and %d (coded 1)."),
        context, comparator, group_n[["0"]], group_n[["1"]]), call. = FALSE)
    }
    if (min(group_n) < 5L || max(group_n) / min(group_n) >= 5) {
      warning(sprintf(paste0(
        "%s: the two `%s` groups are small or unbalanced (%d coded 0, %d coded 1). ",
        "The variance components estimated for the smaller group -- and any sample ",
        "size or power computed from them -- may be unstable."),
        context, comparator, group_n[["0"]], group_n[["1"]]), call. = FALSE)
    }
  }

  # Covariates are expanded only now, on the participants and visits the model
  # will actually use (see the note on the covariates section).
  cov_terms <- character()
  adjusted  <- NULL
  if (!is.null(raw)) {
    if (anyNA(dat$sp_cov_ok)) {
      stop(sprintf(paste0(
        "%s: covariate values are missing for %d participant(s) after `na.action`, ",
        "which kept them. A participant with no recorded covariate value cannot be ",
        "adjusted; use `na.action = na.omit` to remove them."),
        context, length(unique(dat$sp_subject[is.na(dat$sp_cov_ok)]))), call. = FALSE)
    }
    raw <- raw[dat$sp_row, , drop = FALSE]
    dat$sp_row <- dat$sp_cov_ok <- NULL
    check_baseline_covariates(raw, dat$sp_subject, context)
    # Under `healthy` the planned trial enrols cases, so the case slope is
    # reported at the cases' own average covariate values. The pooled mean would
    # give the slope of a population that is part healthy -- shifted, when
    # covariates interact with time, by delta times the gap between the two
    # means. The difference in slopes is the same either way, because the
    # covariate coefficients are shared by both groups.
    reference <- if (comparator == "healthy") dat$sp_case == 1 else rep(TRUE, nrow(dat))
    basis <- covariate_basis(covariates, raw, dat$sp_subject, reference, context)
    if (!length(basis$labels)) {
      warning(sprintf(paste0(
        "%s: every `covariates` column takes the same value for all participants ",
        "in the data used, so no adjustment was made."), context), call. = FALSE)
    } else {
      dat <- cbind(dat, basis$X)
      cov_cols  <- names(basis$labels)
      cov_terms <- c(cov_cols,
                     if (isTRUE(covariate_time)) paste0(cov_cols, ":sp_time"))
      adjusted  <- list(columns = unname(basis$labels), time = isTRUE(covariate_time))
    }
  }

  if (comparator == "treated") {
    # Stata: mixed y time placebo#c.time || subject: time, cov(uns)
    # One common intercept (randomisation implies equal baselines), separate
    # slopes. A numeric placebo indicator keeps the coefficient mapping explicit.
    dat$sp_placebo_time <- (1 - dat$sp_case) * dat$sp_time
  } else if (comparator == "healthy") {
    dat$sp_control      <- 1 - dat$sp_case
    dat$sp_case_time    <- dat$sp_case * dat$sp_time
    dat$sp_control_time <- dat$sp_control * dat$sp_time
    dat$sp_grp <- factor(ifelse(dat$sp_case == 1, "case", "control"),
                         levels = c("control", "case"))
  }
  dat <- add_residual_columns(dat, spec, grid, comparator)
  resid <- residual_calls(spec, comparator)
  fixed <- fixed_formula(comparator, cov_terms)
  if (length(cov_terms)) {
    check_covariate_rank(fixed, dat, basis$labels, comparator, context)
  }

  ctrl <- slope_lme_control()
  reduced_used <- FALSE

  if (comparator != "healthy") {
    fit <- tryCatch(fit_common_model(dat, ctrl, fixed, resid), error = function(e) {
      if (is.null(spec)) stop(e)
      stop(sprintf(paste0("%s: the mixed model with %s residuals did not converge: %s"),
                   context, residual_label(spec$correlation), conditionMessage(e)),
           call. = FALSE)
    })

  } else {
    # Both outcomes are decided inside the handler, so `fit` only ever holds a
    # model or NULL. Capturing the condition into `fit` and testing its class
    # afterwards made "did the full model fail?" a fact about `fit`'s type, and
    # left a window in which the variable held a condition object rather than a
    # fit.
    fit <- NULL
    if (!isTRUE(common_variance)) {
      fit <- tryCatch(
        fit_healthy_model(dat, reduced = FALSE, ctrl = ctrl, fixed = fixed, resid = resid),
        error = function(e) {
          if (isFALSE(common_variance)) {
            stop(sprintf(paste0("%s: the full model did not converge and ",
                                "`common_variance = FALSE` forbids the reduced ",
                                "structure. Underlying error: %s"),
                         context, conditionMessage(e)), call. = FALSE)
          }
          # Without covariates the model factorises per group, so the reduced
          # structure cannot move the cases' estimates; the shared covariate
          # coefficients couple the groups, and then it can.
          message(sprintf(paste0("%s: the full random-effects structure for healthy ",
                                 "controls did not converge; falling back to a ",
                                 "random intercept only for controls (equivalent to ",
                                 "the Stata `nocontvar` option). %s"), context,
                          if (length(cov_terms)) {
                            paste0("Because the covariate coefficients are shared by ",
                                   "both groups, this can shift the estimates ",
                                   "returned for cases slightly.")
                          } else {
                            "This does not affect the estimates returned for cases."
                          }))
          NULL
        })
    }
    if (is.null(fit)) {
      reduced_used <- TRUE
      fit <- tryCatch(fit_healthy_model(dat, reduced = TRUE, ctrl = ctrl, fixed = fixed,
                                        resid = resid),
                      error = function(e) {
                        stop(sprintf("%s: the mixed model did not converge: %s",
                                     context, conditionMessage(e)), call. = FALSE)
                      })
    }
  }

  # One rule for all three models: the slope is the sum of the comparator's
  # fixed-effect parts, and the comparator's own slope is the first of them
  # (the bare time term) -- the treated arm under `treated`, the healthy
  # controls under `healthy`, and nothing at all under `none`, which has a
  # single part. Written per branch, the `treated` and `healthy` blocks were
  # identical apart from a comparator string already in scope, so the guarantee
  # slope_fixef_parts() exists to give -- that the values summed here and the
  # variances summed in slope_se() can never name different coefficients -- held
  # only by inspection. This is the same vapply() over the same mapping that
  # slope_se() (bootstrap.R) already uses.
  b <- nlme::fixef(fit)
  fixef_values <- vapply(slope_fixef_parts(comparator),
                         function(p) fixef_term(b, p, context), numeric(1L))
  slope <- sum(fixef_values)
  slope_comparator <- if (length(fixef_values) > 1L) fixef_values[[1L]] else NA_real_

  # Random-effects covariance and residual variance are extracted by the same
  # names either way: the case-specific random-slope block and residual level
  # for `healthy`, the shared intercept/slope block and homoscedastic residual
  # for the other two -- `treated`'s coefficient mapping differs from `none`'s,
  # but its variance components come from the same random-effects structure.
  re <- if (comparator == "healthy") {
    extract_re(fit, "sp_case", "sp_case_time", context)
  } else {
    extract_re(fit, "(Intercept)", "sp_time", context)
  }
  res <- residual_components(fit, spec, grid, comparator, context)

  new_slope_params(
    slope            = slope,
    slope_comparator = slope_comparator,
    comparator       = comparator,
    sigma2_intercept = re$sigma2_intercept,
    sigma2_slope     = re$sigma2_slope,
    sigma_cov        = re$sigma_cov,
    sigma2_residual  = res$sigma2_residual,
    residual         = res$residual,
    n_obs            = nrow(dat),
    n_subjects       = nlevels(dat$sp_subject),
    common_variance  = reduced_used,
    time_shifted     = time_shifted,
    covariates       = adjusted,
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
#' @param correlation Optional residual correlation, as the `nlme` constructor
#'   carrying the value to assume: `nlme::corAR1(0.6)`, `nlme::corCAR1(0.6)`,
#'   `nlme::corExp(2)`, `nlme::corExp(c(2, 0.1), nugget = TRUE)`,
#'   `nlme::corGaus(2)`, or `nlme::corSymm(c(...))` with the correlations in
#'   `nlme`'s order, (1,2), (1,3), ..., (2,3), .... See "Residual structures"
#'   in [slope_params()] for what each means. `corAR1()` and `corCAR1()` fall
#'   back to `nlme`'s default value when none is given, so state it.
#' @param weights Optional visit-specific residual SD ratios, as
#'   `nlme::varIdent(c("1" = 1.1, "2" = 1.3), form = ~ 1 | visit)`: named by
#'   visit time, relative to the first time in `times`, whose ratio is 1. The
#'   stratum in `form` is required by `nlme` and otherwise ignored.
#'   `sigma2_residual` is then the residual variance at the first time.
#' @param times The visit times that the parameters of `corSymm()` and
#'   `varIdent()` belong to, strictly increasing. Required with either, and
#'   refused otherwise. A planned schedule must then be among these times.
#'
#' @return An object of class `"slope_params"`.
#'
#' @examples
#' # A single group, powered toward a slope of zero. Figures from Table 1,
#' # p.595 of Nash et al. (2021).
#' slope_params_manual(
#'   slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
#'   sigma_cov = 5, sigma2_residual = 10
#' )
#'
#' # Case/healthy-control parameters taken from a published paper, with no
#' # dataset of individual participants available to fit slope_params() to.
#' slope_params_manual(
#'   slope = -1.672, slope_comparator = -0.5,
#'   sigma2_intercept = 100, sigma2_slope = 2,
#'   sigma_cov = 5, sigma2_residual = 10,
#'   comparator = "healthy"
#' )
#'
#' # Serially correlated residuals, with correlation 0.5 between visits a
#' # year apart.
#' slope_params_manual(
#'   slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
#'   sigma_cov = 5, sigma2_residual = 10,
#'   correlation = nlme::corCAR1(0.5)
#' )
#'
#' @seealso [slope_params()] to estimate these from data.
#' @export
slope_params_manual <- function(slope,
                                sigma2_intercept, sigma2_slope,
                                sigma_cov, sigma2_residual,
                                slope_comparator = NA_real_,
                                comparator = c("none", "healthy", "treated"),
                                correlation = NULL, weights = NULL, times = NULL) {
  context <- "slope_params_manual()"
  cl <- match.call()
  comparator <- match.arg(comparator)
  residual <- residual_manual(correlation, weights, times, context)

  # The five components are validated and coerced by new_slope_params() below,
  # which is the single validation point for both routes into the class.
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
    slope            = slope,
    slope_comparator = as.numeric(slope_comparator),
    comparator       = comparator,
    sigma2_intercept = sigma2_intercept,
    sigma2_slope     = sigma2_slope,
    sigma_cov        = sigma_cov,
    sigma2_residual  = sigma2_residual,
    residual         = residual,
    n_obs            = NA_integer_,
    n_subjects       = NA_integer_,
    common_variance  = FALSE,
    time_shifted     = FALSE,
    covariates       = NULL,
    fit              = NULL,
    call             = cl,
    context          = context
  )
}

#' Validate and build the object
#'
#' The single validation point for both routes into the class: the fitted one
#' through [slope_params()] and the direct one through [slope_params_manual()].
#' The checks also coerce -- `check_scalar()` returns its argument as a double --
#' so components reach the object in the type CONTRACT.md section 2 specifies
#' without a separate `as.numeric()` pass that would turn a non-numeric argument
#' into `NA` and a coercion warning before it could be reported properly.
#' @noRd
new_slope_params <- function(slope, slope_comparator, comparator,
                             sigma2_intercept, sigma2_slope, sigma_cov,
                             sigma2_residual, residual, n_obs, n_subjects,
                             common_variance, time_shifted, covariates, fit,
                             call, context) {
  v <- check_param_values(list(slope = slope, sigma2_intercept = sigma2_intercept,
                                sigma2_slope = sigma2_slope, sigma2_residual = sigma2_residual,
                                sigma_cov = sigma_cov), context)
  check_residual(residual, context)

  structure(
    list(slope            = v$slope,
         slope_comparator = slope_comparator,
         comparator       = comparator,
         sigma2_intercept = v$sigma2_intercept,
         sigma2_slope     = v$sigma2_slope,
         sigma_cov        = v$sigma_cov,
         sigma2_residual  = v$sigma2_residual,
         residual         = residual,
         n_obs            = n_obs,
         n_subjects       = n_subjects,
         common_variance  = common_variance,
         time_shifted     = time_shifted,
         covariates       = covariates,
         fit              = fit,
         call             = call),
    class = "slope_params")
}

# ---- printing ---------------------------------------------------------------

#' What the two slopes of a `slope_params` object are called in printed output
#'
#' Naming the parts of the object belongs with the class, not with each thing
#' that renders it: [print.slope_params()] and `print_data_block()` in power.R
#' show the same two quantities and must not disagree about what they are, and
#' the paper-parity tests pin these strings exactly. Whether a comparator slope
#' is shown at all is a separate, per-caller decision -- `print_data_block()`
#' hides it unless `target = "observed"` -- and stays where it is made.
#'
#' @return A list with `own` (the untreated / case / control-arm slope),
#'   `comparator`, and `difference` -- the label for the gap between the two,
#'   which does not vary with the comparator but belongs with its siblings
#'   rather than being spelled once here and once in [print_data_block()].
#' @noRd
slope_labels <- function(comparator) {
  labs <- if (identical(comparator, "treated")) {
    list(own = "slope of control arm", comparator = "slope of experimental arm")
  } else {
    list(own = "slope of cases", comparator = "slope of healthy controls")
  }
  c(labs, list(difference = "observed difference in slopes"))
}

#' Print stage-one slope parameters
#'
#' @param x A `"slope_params"` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#'
#' @examples
#' slope_params(sdmt ~ visit | id, data = slpower1)
#'
#' @export
print.slope_params <- function(x, ...) {
  lab <- switch(x$comparator,
                none    = "single group (target: no change over time)",
                healthy = "cases and healthy controls",
                treated = "previous randomised trial")
  cat("Slope parameters (", lab, ")\n\n", sep = "")

  if (!is.na(x$n_obs)) {
    cat_line("number of observations in model", x$n_obs, digits = 0L)
    cat_line("number of participants in model", x$n_subjects, digits = 0L)
  } else {
    cat_line("source", "supplied directly")
  }

  labels <- slope_labels(x$comparator)
  cat_line(labels$own, x$slope)
  if (!is.na(x$slope_comparator)) {
    cat_line(labels$comparator, x$slope_comparator)
    cat_line(labels$difference, x$slope - x$slope_comparator)
  }

  cat("\n")
  cat_line("variance of random intercepts", x$sigma2_intercept)
  cat_line("variance of random slopes", x$sigma2_slope)
  cat_line("covariance of intercept and slope", x$sigma_cov)
  print_residual(x)

  covariate_note(x$covariates)
  residual_note(x$residual)
  # A statement about the fit, so not shown for components supplied directly.
  if (identical(x$residual$correlation, "corSymm") && !is.null(x$fit)) {
    cat(paste0("\nNote: with an unstructured residual correlation the variance components\n",
               "      above are not separately identified -- only the covariance they imply\n",
               "      at the visit times is, and that is all stage two uses.\n"))
  }
  if (isTRUE(x$common_variance)) {
    cat("\nNote: reduced random-effects structure used for healthy controls\n")
    # The shared covariate coefficients couple the two groups; see the
    # fallback message in slope_params().
    cat(if (is.null(x$covariates)) "      (Stata `nocontvar`). Case estimates are unaffected.\n"
        else "      (Stata `nocontvar`). With covariates, this can shift the case\n      estimates slightly.\n")
  }
  if (isTRUE(x$time_shifted)) {
    cat("\nNote: subject times were shifted so each first visit is time zero.\n")
  }
  invisible(x)
}
