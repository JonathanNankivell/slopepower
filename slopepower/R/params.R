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
#' `healthy = <"case"/"control" column>` return the healthy controls' slope
#' labelled as the cases'.
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
    if (!identical(levels(f), c("0", "1"))) {
      meaning <- if (identical(name, "healthy")) "case" else "treated"
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
                 context, name,
                 if (identical(name, "healthy")) "case" else "treated",
                 paste(format(u), collapse = "/")), call. = FALSE)
  }
  as.numeric(x)
}

#' Look up a fixed effect by name, resolving interactions whichever way round
#' R happened to spell them
#'
#' `lme` names an interaction after the order its components appear in the model
#' formula, so `sp_case * sp_time` yields `sp_case:sp_time` while
#' `sp_time * sp_case` yields `sp_time:sp_case`. The two fits are identical --
#' only the label differs -- so extraction must not depend on which spelling
#' arose. Indexing a fixed spelling would break on an innocuous reordering of
#' the internal formula, loudly rather than silently, but break nonetheless.
#' @noRd
fixef_term <- function(b, parts, context) {
  cand <- unique(c(paste(parts, collapse = ":"), paste(rev(parts), collapse = ":")))
  hit <- cand[cand %in% names(b)]
  if (length(hit) == 0L) {
    stop(sprintf("%s: fixed effect `%s` not found in the fitted model; have %s.",
                 context, cand[1L], paste(names(b), collapse = ", ")), call. = FALSE)
  }
  unname(b[[hit[1L]]])
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
#' [slope_sample_size()] and [slope_power()] need.
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
#'   available via `target = "observed"` in [slope_sample_size()] and
#'   [slope_power()].
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
#' off-diagonal \eqn{\sigma_{ab}}{sigma_cov}. Residuals are independent across
#' visits and across participants; there is no serial correlation term. Only the
#' mean \eqn{\mu}{mu} and the number of variance parameters differ between the
#' scenarios.
#'
#' \describe{
#'   \item{Neither `healthy` nor `treated`}{
#'     \deqn{\mu(t) = \beta_0 + \beta_1 t}{mu(t) = b0 + b1 * t}
#'     with one \eqn{G}{G} and one \eqn{\sigma^2_\epsilon}{sigma2_residual}. The
#'     returned slope is \eqn{\beta_1}{b1}. Equivalent to
#'     `nlme::lme(y ~ t, random = ~ t | id)`.}
#'   \item{`healthy = g`, with \eqn{g_i = 1}{g[i] = 1} for cases}{
#'     \deqn{\mu(t) = \beta_0 + \beta_g g_i + (\beta_1 + \beta_{1g} g_i) t}{
#'           mu(t) = b0 + bg * g[i] + (b1 + b1g * g[i]) * t}
#'     and, in addition, a **separate** \eqn{G}{G} and a **separate**
#'     \eqn{\sigma^2_\epsilon}{sigma2_residual} for each group. The slope of the
#'     cases is \eqn{\beta_1 + \beta_{1g}}{b1 + b1g} and that of the controls
#'     \eqn{\beta_1}{b1}; the variance components returned are the cases'. The
#'     model factorises into two independent per-group fits (see above).}
#'   \item{`treated = z`, with \eqn{z_i = 1}{z[i] = 1} for the treated arm}{
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
#' `y ~ group * time + (time | id)`: the `treated` model drops the group main
#' effect, and the `healthy` model gives each group its own residual variance,
#' which `lme4` cannot fit at all. Included, therefore:
#'
#' * one continuous, approximately Gaussian outcome;
#' * a mean that is linear in time, and nothing else --- no other covariates;
#' * exactly one grouping level, the participant, whose random intercept and
#'   random slope have an unstructured covariance;
#' * independent residuals with a variance that is constant within a group;
#' * at most two groups, distinguished only by their slope (and, for `healthy`,
#'   by their variance components).
#'
#' Not included, and not obtainable by any argument to this function:
#'
#' * **covariate adjustment.** Age, sex, disease duration, centre and the rest
#'   cannot be entered. `formula` takes a single time term and rejects anything
#'   more. If the trial you are planning will adjust for prognostic covariates,
#'   its residual variance will be smaller than the one estimated here and the
#'   resulting sample size is conservative.
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
#' * **structured residuals**, such as AR(1) or other serial correlation within
#'   a participant.
#' * **more than two arms, unequal allocation, or cluster-randomised, crossover
#'   and stepped-wedge designs.** Stage two assumes two equal parallel arms.
#'
#' @return An object of class `"slope_params"`.
#'
#' @examples
#' # Neither `healthy` nor `treated`: a single group of untreated subjects.
#' # Four of the two hundred participants of `slpower1`, kept small so the
#' # example runs quickly -- see `slpower1` for the paper's fit on the full data.
#' df <- data.frame(
#'   id    = c(1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4),
#'   visit = c(0,1,2,3, 0,1,2,3, 0,1,2,3, 0,1,2,3),
#'   sdmt  = c(41,33,25,30, 16,14,13,6, 35,27,31,23, 38,33,33,20)
#' )
#' slope_params(sdmt ~ visit | id, data = df)
#'
#' # `healthy`: two cases and two healthy controls, a subset of `slpower2`.
#' # Visits are recorded as calendar dates there, so the time term converts
#' # them to years.
#' df2 <- data.frame(
#'   id    = rep(c(1, 2, 251, 252), each = 4),
#'   case  = rep(c(0, 0, 1, 1), each = 4),
#'   vdate = as.Date(c("2009-07-11", "2010-06-21", "2011-07-06", "2012-05-06",
#'                     "2009-07-24", "2010-03-22", "2011-04-14", "2012-06-17",
#'                     "2009-06-06", "2010-08-13", "2011-08-31", "2012-09-05",
#'                     "2009-06-27", "2010-04-11", "2011-06-22", "2012-10-13")),
#'   sdmt  = c(40,46,41,45, 43,42,43,39, 35,34,36,39, 25,16,18,12)
#' )
#' slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, data = df2, healthy = case)
#'
#' # `treated`: data from a completed trial, a subset of `slpower3`. Fitting the
#' # random-effects structure shared by both arms needs more than a couple of
#' # subjects per arm to converge, so this excerpt keeps six per arm.
#' df3 <- data.frame(
#'   id    = rep(c(1:6, 76:81), each = 3),
#'   treat = rep(c(0,0,0,0,0,0, 1,1,1,1,1,1), each = 3),
#'   visit = rep(c(0, 0.5, 2), times = 12),
#'   sdmt  = c(35,36,34, 32,34,29, 27,29,22, 15,22,13, 16,17,14, 29,27,28,
#'             29,33,35, 19,19,17, 19,14,7, 37,38,32, 19,18,23, 43,44,43)
#' )
#' slope_params(sdmt ~ visit | id, data = df3, treated = treat)
#'
#' @references
#' Nash, S., Morgan, K. E., Frost, C. and Mulick, A. (2021). Power and
#' sample-size calculations for trials that compare slopes over time:
#' Introducing the slopepower command. \emph{The Stata Journal} 21(3): 575--601.
#'
#' @seealso [slope_params_manual()] to supply parameters directly,
#'   [slope_sample_size()] and [slope_power()] for stage two.
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
    slope <- fixef_term(b, "sp_time", context)
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
    slope            <- fixef_term(b, "sp_time", context) +
                        fixef_term(b, "sp_placebo_time", context)   # control arm
    slope_comparator <- fixef_term(b, "sp_time", context)                 # treated arm
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
    slope            <- fixef_term(b, "sp_time", context) +
                        fixef_term(b, c("sp_case", "sp_time"), context)  # cases
    slope_comparator <- fixef_term(b, "sp_time", context)                  # controls
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
