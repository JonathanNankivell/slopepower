# Layer 4 -- uncertainty in the stage-one estimates.
#
# Section 2.6 of Nash et al. (2021) recommends bootstrapping the estimated sample
# size, because the variance components that drive it are themselves estimated and
# can be imprecise. The Stata recipe requires the user to assemble four separate
# options correctly -- cluster(), idcluster(), strata() and jack() -- and warns
# that it silently assumes no observations were excluded. Here the resampling
# scheme is fixed by construction: subjects are the sampling unit, replicates get
# fresh identifiers so a subject drawn twice counts as two people, and groups are
# stratified automatically.

#' Recover the prepared modelling frame from a fitted `slope_params` object
#' @noRd
boot_frame <- function(params, context) {
  if (is.null(params$fit)) {
    stop(sprintf(paste0("%s: `params` has no fitted model, so there is nothing to resample. ",
                        "Objects from slope_params_manual() carry parameter values only; ",
                        "bootstrap the study they came from instead."), context), call. = FALSE)
  }
  g <- tryCatch(nlme::getData(params$fit), error = function(e) NULL)
  if (!is.data.frame(g) || !all(c("sp_y", "sp_time", "sp_subject") %in% names(g))) {
    stop(sprintf(paste0("%s: could not recover the modelling data from `params$fit`. ",
                        "Refit with slope_params() in this session and try again."),
                 context), call. = FALSE)
  }
  out <- data.frame(y = as.numeric(g$sp_y),
                    time = as.numeric(g$sp_time),
                    subject = g$sp_subject,
                    stringsAsFactors = FALSE)
  if (!identical(params$comparator, "none")) {
    if (!"sp_case" %in% names(g)) {
      stop(sprintf("%s: the fitted model has no group indicator to stratify on.", context),
           call. = FALSE)
    }
    out$group <- as.integer(g$sp_case)
  }
  out
}

#' Build a closure that refits stage one on a resampled frame
#'
#' The call to `slope_params()` -- formula, `origin = "none"` (the times in the
#' recovered frame have already been re-origined per subject, so this avoids
#' repeating that work and the message that goes with it), and which comparator
#' argument to pass -- is identical for every replicate of a given bootstrap;
#' only the resampled frame differs. Assembling it once here, rather than once
#' per replicate, means [run_bootstrap()]'s loop -- run several hundred times,
#' plus once per subject for the BCa jackknife -- reconstructs only what
#' actually changes between replicates.
#' @noRd
make_refitter <- function(comparator) {
  args <- list(formula = y ~ time | subject, data = quote(frame), origin = "none")
  cl <- as.call(c(list(quote(slope_params)), args,
                  if (identical(comparator, "healthy")) list(healthy = quote(group))
                  else if (identical(comparator, "treated")) list(treated = quote(group))
                  else NULL))
  # Evaluated with the closure's own frame as the enclosure, so the free symbol
  # `slope_params` resolves through this function's environment into the package
  # namespace. Rooting it in `parent.frame()` instead looked the name up in
  # whichever frame happened to call `make_refitter()` -- found only when the
  # package is attached, and shadowed by a user object of the same name. See the
  # same reasoning, at length, beside the call `slopepower()` builds in compat.R.
  function(frame) {
    suppressMessages(suppressWarnings(eval(cl, list(frame = frame))))
  }
}

#' Build one cluster-resampled frame, stratified by group where present
#'
#' `groups` holds each stratum's members as *positions* in `subject_index`
#' rather than as subject-label strings, so that looking a pick up is a
#' constant-time integer index into the list rather than a linear scan of its
#' names -- the same random draws either way, since `sample()`'s draws depend
#' only on a vector's length, not its values.
#'
#' Drawing is done via `sample.int(length(pos), ...)`, indexed back into `pos`,
#' rather than `sample(pos, ...)` directly. `sample(x, size)` treats a
#' length-one `x` as a range to draw *from* (`1:x`) rather than a single value
#' to draw *with replacement*, so a stratum of exactly one subject -- one case
#' in a `healthy` comparison, say -- would silently have its "resample" drawn
#' from the wrong distribution instead of always returning that subject.
#' @noRd
resample_frame <- function(frame, subject_index, groups) {
  picks <- unlist(lapply(groups, function(pos) {
    pos[sample.int(length(pos), length(pos), replace = TRUE)]
  }), use.names = FALSE)
  picked <- subject_index[picks]
  out <- frame[unlist(picked, use.names = FALSE), , drop = FALSE]
  out$subject <- rep(seq_along(picks), lengths(picked))  # fresh identifiers:
  out                              # a subject drawn twice is two people
}

#' Standard error of the estimated slope
#'
#' The standard error of the untreated (or case) slope, taken from the fitted
#' mixed model's fixed-effects covariance matrix: for `comparator = "healthy"`
#' or `"treated"` this combines the variances of, and covariance between, the
#' two fixed-effect terms whose sum is the slope, not just the variance of a
#' single coefficient. [slope_bootstrap()] compares this to the slope itself
#' as the check recommended in section 2.6 of Nash et al. (2021): a slope
#' less than 2.5 times its standard error means bootstrap replicates can
#' straddle zero, at which point the resulting interval stops meaning
#' anything.
#'
#' @param params A `slope_params` object produced by [slope_params()]. Objects
#'   from [slope_params_manual()] carry no fitted model, so `NA_real_` is
#'   returned for them.
#'
#' @return A single non-negative number, or `NA_real_` if `params` has no
#'   fitted model, or if the slope terms could not be identified in it (with
#'   a warning).
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_se(pars)
#'
#' # slope_params_manual() objects carry no fitted model, so there is no
#' # standard error to report.
#' slope_se(slope_params_manual(
#'   slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
#'   sigma_cov = 5, sigma2_residual = 10
#' ))
#'
#' @seealso [slope_sigma()] and [slope_var()], the other quantities computed
#'   from a `slope_params` object; [slope_bootstrap()], which uses this to
#'   flag an unreliable interval.
#' @export
slope_se <- function(params) {
  context <- "slope_se()"
  check_params(params, context)
  fit <- params$fit
  if (is.null(fit)) return(NA_real_)
  b <- tryCatch(nlme::fixef(fit), error = function(e) NULL)
  V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(b) || is.null(V)) return(NA_real_)
  # Which terms sum to the slope, by comparator, is params.R's
  # slope_fixef_parts() -- the same mapping slope_params() itself sums the
  # *values* of, via fixef_term(), so the two can never name a different set
  # of coefficients. Every part -- not just an interaction -- is resolved via
  # resolve_fixef_name() rather than assuming a spelling, exactly as
  # fixef_term() does: for a single name it degenerates to a plain lookup, so
  # one call handles both shapes and this can never trust an unresolved name
  # the way indexing `p` directly would. Getting this wrong used to return NA
  # silently, which switched off the section 2.6 warning below that is the
  # entire reason for computing the standard error.
  terms <- vapply(slope_fixef_parts(params$comparator),
                  function(p) resolve_fixef_name(b, p),
                  character(1L))
  # `anyNA()` alone: resolve_fixef_name() returns either a member of names(b)
  # or NA, so a non-NA element being absent from names(b) is not a second
  # failure mode to test for.
  if (anyNA(terms)) {
    warning(sprintf(paste0(
      "%s: could not identify the slope terms in the fitted model (have %s); ",
      "returning NA. If this call came from slope_bootstrap(), its section 2.6 ",
      "check on the slope-to-standard-error ratio was skipped."),
      context, paste(names(b), collapse = ", ")), call. = FALSE)
    return(NA_real_)
  }
  k <- as.numeric(names(b) %in% terms)
  sqrt(drop(k %*% as.matrix(V) %*% k))
}

#' Reject anything left in `...`
#'
#' The methods take no pass-through arguments: everything the calculation needs
#' is already in the object being bootstrapped. Silently ignoring a stray
#' `design =` or `n =` would be the worst outcome, because the result would look
#' like a successful bootstrap of something else entirely -- and that is exactly
#' the shape of a call written against the pre-generic interface, where the
#' calculation was re-specified here rather than dispatched on.
#'
#' `context` rather than a hard-coded `"slope_bootstrap()"`: the same hazard,
#' and the same fix, applies to every generic in the package whose methods
#' carry `...` only because the generic does. [slope_sample_size_floor()] is
#' the second.
#' @noRd
reject_dots <- function(dots, advice, context) {
  if (length(dots) == 0L) return(invisible(NULL))
  nms <- names(dots) %||% rep("", length(dots))
  shown <- ifelse(nzchar(nms), nms, "<unnamed>")
  stop(sprintf("%s: unused argument%s (%s).\n  %s", context,
               if (length(dots) > 1L) "s" else "",
               paste(shown, collapse = ", "), advice), call. = FALSE)
}

#' `match.arg()` for `statistic`, with a message that names the object
#'
#' Which statistics are on offer is now a property of the object being
#' bootstrapped, so a rejected one almost always means the wrong object was
#' handed over -- `statistic = "power"` on a result that solved for the sample
#' size, say. Bare `match.arg()` reports only `'arg' should be one of ...`, which
#' names neither the argument nor the object and leaves the caller to work out
#' that the fix is upstream, in the call that built `x`. The matching itself
#' (including the "whole default vector passed through" case) is still
#' `match.arg()`'s; only the error message is replaced.
#' @noRd
match_statistic <- function(statistic, choices, advice) {
  matched <- tryCatch(match.arg(statistic, choices), error = function(e) NULL)
  if (is.null(matched)) {
    stop(sprintf("slope_bootstrap(): `statistic` must be %s, not %s.\n  %s",
                 paste(sQuote(choices), collapse = " or "),
                 sQuote(paste(as.character(statistic), collapse = ", ")), advice),
         call. = FALSE)
  }
  matched
}

#' The advice half of that message, for the two stage-two methods
#' @noRd
dots_advice_result <- function(entry) {
  paste0("The calculation comes from the object being bootstrapped, so `design`,\n",
         "  `effectiveness`, `target`, `alpha` and the sample size or power belong\n",
         "  to the ", entry, " call that produced it, not here.")
}

#' Re-solve a stage-two result against resampled parameters
#'
#' Rebuilds the call that produced `result`, substituting the parameters fitted
#' to one bootstrap replicate. Everything needed is stored on the result:
#' `design` is already a `trial_design` object, and `alpha`, `target` and the
#' solved-for input are carried alongside it.
#'
#' `effectiveness` is omitted under `target = "observed"` via
#' `maybe_add_effectiveness()`, shared with the grid functions' base args.
#' @noRd
resolve_args <- function(p, result) {
  args <- list(params = p, design = result$design,
               target = result$target, alpha = result$alpha)
  maybe_add_effectiveness(args, result$effectiveness, result$target)
}

#' Bootstrap a scalar computed from resampled stage-one parameters
#'
#' The resampling scheme, the section 2.6 check and the interval construction are
#' the same whatever is being bootstrapped; only `compute` differs, and each
#' method supplies one closure that reads its statistic off a refit.
#' @noRd
run_bootstrap <- function(params, compute, observed, statistic, R, type, level,
                          seed, progress) {
  context <- "slope_bootstrap()"
  type <- match.arg(type, c("bca", "percentile"))
  check_whole_number(R, "R", "replicates", context, lower = 1)
  check_probability(level, "level", context)
  if (!is.null(seed)) set.seed(seed)

  # `observed` is read off the object rather than recomputed. It is the same
  # number either way -- `compute` on the original parameters reproduces the
  # result it was handed -- but recomputing would re-emit any warning the
  # stage-two call made, and under this interface the caller has already run
  # that exact call themselves to produce the object.

  # Section 2.6: if the slope is not large relative to its standard error the
  # replicates can straddle zero, and an interval for the sample size stops
  # meaning anything.
  se <- slope_se(params)
  if (is.finite(se) && se > 0 && abs(params$slope) / se < 2.5) {
    warning(sprintf(paste0("%s: the estimated slope (%.4g) is only %.2f times its standard error ",
                           "(%.4g). Nash et al. (2021, section 2.6) suggest a ratio above 2.5; ",
                           "below it, bootstrap replicates may include slopes of both signs and ",
                           "the resulting interval can be meaningless."),
                    context, params$slope, abs(params$slope) / se, se), call. = FALSE)
  }

  frame <- boot_frame(params, context)
  subject_index <- split(seq_len(nrow(frame)), frame$subject)
  positions <- seq_along(subject_index)
  groups <- if (is.null(frame$group)) {
    list(positions)
  } else {
    sub_group <- vapply(subject_index, function(r) frame$group[r[1L]], numeric(1L))
    unname(split(positions, sub_group))
  }
  refitter <- make_refitter(params$comparator)

  # Warnings are suppressed for every replicate. Anything worth saying about the
  # calculation -- a target effect that makes the slope more extreme, say -- is a
  # property of the data, and the caller has already heard it once from the
  # stage-two call that built the object; repeating it several hundred times
  # here, attributed to a function they did not call, tells them nothing new.
  # Suppressed once around the whole loop rather than per replicate: messages
  # (the progress ticks below) are a different condition class and pass through
  # regardless.
  replicates <- rep(NA_real_, R)
  tick <- max(1L, floor(R / 10))
  suppressWarnings(for (b in seq_len(R)) {
    replicates[b] <- tryCatch(
      compute(refitter(resample_frame(frame, subject_index, groups))),
      error = function(e) NA_real_)
    if (isTRUE(progress) && b %% tick == 0L) {
      message(sprintf("%s: %d of %d replicates", context, b, R))
    }
  })

  failed <- is.na(replicates)
  n_failed <- sum(failed)
  good <- replicates[!failed]
  if (length(good) < 2L) {
    stop(sprintf("%s: %d of %d replicates failed; not enough succeeded to form an interval.",
                 context, n_failed, R), call. = FALSE)
  }
  if (n_failed > 0L) {
    warning(sprintf("%s: %d of %d replicates failed to converge and were discarded.",
                    context, n_failed, R), call. = FALSE)
  }

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  ci <- stats::quantile(good, probs, names = FALSE, type = 7)
  used_type <- "percentile"

  if (identical(type, "bca")) {
    bca <- bca_interval(good, observed, frame, subject_index, refitter,
                        compute, probs)
    if (is.null(bca)) {
      warning(sprintf(paste0("%s: the bias-correction could not be computed (every replicate ",
                             "falls on one side of the observed value); reporting a percentile ",
                             "interval instead."), context), call. = FALSE)
    } else {
      ci <- bca
      used_type <- "bca"
    }
  }

  structure(list(observed = observed, replicates = good, ci = ci,
                 type = used_type, statistic = statistic, R = R,
                 n_failed = n_failed, level = level, se = se),
            class = "slope_bootstrap")
}

#' Bootstrap the stage-one estimates
#'
#' Resamples subjects with replacement, refits the stage-one mixed model on each
#' replicate, and recomputes a quantity of interest. This propagates the
#' estimation uncertainty in the slope and variance components through to the
#' sample size or power, as recommended in section 2.6 of Nash et al. (2021).
#'
#' Hand it the result you want an interval around, not a fresh specification of
#' the calculation. A `slope_sample_size` object knows the design, effectiveness,
#' target power and significance level it was solved with, so the bootstrap
#' re-solves exactly that calculation on each replicate:
#'
#' ```
#' ss <- slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#' slope_bootstrap(ss, R = 999, type = "bca")
#' ```
#'
#' Dispatching on the result rather than on a `statistic` argument also makes it
#' impossible to bootstrap one of the calculation's own inputs -- an interval
#' around the `n` you supplied yourself is zero-width, and used to be an easy
#' call to write. Each method offers only quantities its object solved for or
#' derived.
#'
#' Subjects, not observations, are the sampling unit, and each drawn subject is
#' given a fresh identifier so that a subject selected twice is treated as two
#' people rather than as one person with twice as much data. When the parameters
#' came from two-group data, resampling is stratified so that each replicate has
#' the same number of cases and controls (or treated and control subjects) as the
#' original.
#'
#' Refitting a mixed model several hundred times is slow, and `type = "bca"` adds
#' a leave-one-subject-out jackknife on top, costing one further fit per subject.
#' Start with a small `R` to gauge the cost.
#'
#' @param x What to bootstrap: a `slope_sample_size` object from
#'   [slope_sample_size()], a `slope_power` object from [slope_power()], or a
#'   `slope_params` object from [slope_params()] for the fitted slope itself.
#'   The underlying parameters must come from [slope_params()] in either case;
#'   [slope_params_manual()] objects carry no data to resample. For the `print()`
#'   method, the `slope_bootstrap` object to show.
#' @param R Number of bootstrap replicates.
#' @param type `"bca"` (the default) for bias-corrected and accelerated
#'   intervals, or `"percentile"`. The paper recommends BCa because the
#'   distribution of estimated sample sizes is typically skewed.
#' @param ... Not used. The calculation is taken from `x`, so any argument here
#'   is an error rather than something silently ignored.
#' @param level Confidence level for the interval.
#' @param seed Optional integer seed, for reproducibility.
#' @param progress Report progress while resampling.
#'
#' @return An object of class `slope_bootstrap`, with elements `observed`, `replicates`,
#'   `ci`, `type`, `statistic`, `R`, `n_failed`, `level` and `se`.
#'
#' @references
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#'
#' @examples
#' # The parameters must come from slope_params(): bootstrapping resamples the
#' # underlying subjects, which slope_params_manual() objects do not carry.
#'
#' # No comparator: all two hundred participants of `slpower1`.
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' ss <- slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#'
#' # One mixed-model fit per replicate, so a real run wants a much larger R.
#' \donttest{
#' slope_bootstrap(ss, R = 100, seed = 42)
#'
#' # The same result also carries the target treatment effect.
#' slope_bootstrap(ss, R = 100, statistic = "tte", seed = 42)
#'
#' # An interval around the slope needs no trial design at all.
#' slope_bootstrap(pars, R = 100, seed = 42)
#' }
#'
#' # Case/healthy-control comparator: forty cases and forty healthy controls,
#' # simulated (rather than drawn from `slpower2`, whose 500 participants
#' # would make every replicate refit the slowest model in the package) and
#' # built subject by subject so each keeps one intercept and one slope
#' # across its four visits.
#' set.seed(2)
#' subj2 <- data.frame(id = 1:80, case = rep(c(1, 0), each = 40))
#' subj2$intercept <- rnorm(80, 50, 10)
#' subj2$slope     <- rnorm(80, ifelse(subj2$case == 1, -1.7, -0.3), 1.4)
#' sim2 <- merge(subj2, data.frame(visit = 0:3))
#' sim2$sdmt <- sim2$intercept + sim2$slope * sim2$visit + rnorm(nrow(sim2), 0, 3)
#' pars2 <- slope_params(sdmt ~ visit | id, data = sim2, healthy = case)
#'
#' \donttest{
#' # Bootstrapping the power of a fixed sample size, rather than the size itself.
#' pw2 <- slope_power(pars2, c(0, 1, 2), n = 400, effectiveness = 0.33)
#' slope_bootstrap(pw2, R = 100, seed = 42)
#' }
#'
#' # Randomised-trial comparator, target = "observed": all one hundred and
#' # fifty participants of `slpower3`, bootstrapping the sample size needed
#' # to detect the effect the trial actually found.
#' pars3 <- slope_params(sdmt ~ visit | id, data = slpower3, treated = treat)
#'
#' \donttest{
#' ss3 <- slope_sample_size(pars3, c(0, 0.5, 2), target = "observed")
#' slope_bootstrap(ss3, R = 100, seed = 42)
#' }
#'
#' @seealso [slope_sample_size()], [slope_power()], [slope_params()],
#'   [slope_se()] for the standard error behind the section 2.6 check
#' @export
slope_bootstrap <- function(x, R = 199, type = c("bca", "percentile"), ...,
                            level = 0.95, seed = NULL, progress = FALSE) {
  UseMethod("slope_bootstrap")
}

#' Shared body of the two stage-two `slope_bootstrap()` methods
#'
#' `slope_bootstrap.slope_sample_size()` and `slope_bootstrap.slope_power()`
#' differ only in which stage-two function is re-solved on each replicate,
#' which of the object's own inputs is held fixed while doing so, and which
#' statistics are on offer; everything else -- matching `statistic`, rejecting
#' `...`, and the call to `run_bootstrap()` -- is identical.
#' @noRd
bootstrap_stage_two <- function(x, fn, fixed_name, choices, advice, label,
                                R, type, statistic, level, seed, progress, dots) {
  statistic <- match_statistic(statistic, choices, advice)
  reject_dots(dots, dots_advice_result(label), "slope_bootstrap()")
  fixed <- stats::setNames(list(x[[fixed_name]]), fixed_name)
  # `compute` closes over `slim` -- just the four inputs resolve_args() reads --
  # rather than `x` itself, so it does not keep x$params$fit (the original fit
  # and its model frame) reachable through every one of several hundred replicates.
  slim <- x[c("design", "target", "alpha", "effectiveness")]
  compute <- function(p) do.call(fn, c(resolve_args(p, slim), fixed))[[statistic]]
  run_bootstrap(x$params, compute, x[[statistic]], statistic, R, type, level,
               seed, progress)
}

#' @describeIn slope_bootstrap Bootstrap the required sample size (the default)
#'   or the target treatment effect behind it.
#' @param statistic Which of the object's quantities to bootstrap. Each method
#'   offers only what its object solved for or derived, and defaults to the
#'   quantity the object exists to report.
#' @export
slope_bootstrap.slope_sample_size <- function(x, R = 199,
                                              type = c("bca", "percentile"),
                                              statistic = c("n", "tte"), ...,
                                              level = 0.95, seed = NULL,
                                              progress = FALSE) {
  # `power` is the target this result was solved to, and so is an input that is
  # held fixed across replicates -- it is what makes `n` vary.
  bootstrap_stage_two(x, slope_sample_size, "power", c("n", "tte"), paste0(
    "This result solved for the sample size, so `n` and the target treatment\n",
    "  effect `tte` behind it are what it can offer. For the power a fixed\n",
    "  sample size achieves, bootstrap a slope_power() result instead."),
    "slope_sample_size()", R, type, statistic, level, seed, progress, list(...))
}

#' @describeIn slope_bootstrap Bootstrap the power achieved (the default) or the
#'   target treatment effect behind it.
#' @export
slope_bootstrap.slope_power <- function(x, R = 199,
                                        type = c("bca", "percentile"),
                                        statistic = c("power", "tte"), ...,
                                        level = 0.95, seed = NULL,
                                        progress = FALSE) {
  # `x$n` rather than `x$n_requested`: the even number actually used, so the
  # replicates answer the question the observed value answered.
  bootstrap_stage_two(x, slope_power, "n", c("power", "tte"), paste0(
    "This result solved for the power a fixed sample size achieves, so `power`\n",
    "  and the target treatment effect `tte` behind it are what it can offer. For\n",
    "  the sample size a target power needs, bootstrap a slope_sample_size()\n",
    "  result instead."),
    "slope_power()", R, type, statistic, level, seed, progress, list(...))
}

#' @describeIn slope_bootstrap Bootstrap the fitted slope itself, which needs no
#'   trial design.
#' @export
slope_bootstrap.slope_params <- function(x, R = 199,
                                         type = c("bca", "percentile"), ...,
                                         level = 0.95, seed = NULL,
                                         progress = FALSE) {
  reject_dots(list(...), paste0(
    "Bootstrapping a `slope_params` object gives an interval for the fitted\n",
    "  slope, which needs no trial design. For an interval around a sample size\n",
    "  or a power, bootstrap the result instead:\n",
    "    slope_bootstrap(slope_sample_size(params, design, ...), R = 999)"),
    "slope_bootstrap()")
  run_bootstrap(x, function(p) p$slope, x$slope, "slope", R, type, level, seed,
                progress)
}

#' @describeIn slope_bootstrap Reject anything else, with a pointer to what is
#'   accepted.
#' @export
slope_bootstrap.default <- function(x, R = 199, type = c("bca", "percentile"),
                                    ..., level = 0.95, seed = NULL,
                                    progress = FALSE) {
  stop(sprintf(paste0(
    "slope_bootstrap(): cannot bootstrap an object of class %s. Pass the result\n",
    "  you want an interval around -- a slope_sample_size object from\n",
    "  slope_sample_size(), a slope_power object from slope_power(), or a\n",
    "  slope_params object from slope_params() for the slope itself.\n",
    "  The grid functions return one row per design; bootstrap the design you\n",
    "  settle on rather than the grid."),
    paste(sQuote(class(x)), collapse = "/")), call. = FALSE)
}

#' Bias-corrected and accelerated interval
#'
#' The acceleration comes from a leave-one-subject-out jackknife, matching the
#' clustering used for the bootstrap itself.
#' @noRd
bca_interval <- function(theta, observed, frame, subject_index, refitter,
                         compute, probs) {
  prop <- mean(theta < observed)
  if (prop <= 0 || prop >= 1) return(NULL)
  z0 <- stats::qnorm(prop)

  jack <- rep(NA_real_, length(subject_index))
  suppressWarnings(for (i in seq_along(subject_index)) {
    drop_rows <- subject_index[[i]]
    jack[i] <- tryCatch(compute(refitter(frame[-drop_rows, , drop = FALSE])),
                        error = function(e) NA_real_)
  })
  jack <- jack[!is.na(jack)]
  if (length(jack) < 3L) return(NULL)

  d <- mean(jack) - jack
  denom <- 6 * (sum(d^2))^1.5
  a <- if (denom == 0) 0 else sum(d^3) / denom

  z <- stats::qnorm(probs)
  adj <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  if (any(!is.finite(adj))) return(NULL)
  stats::quantile(theta, adj, names = FALSE, type = 7)
}

#' @describeIn slope_bootstrap Print a bootstrap result.
# No `@param x` here: this block and the generic share one help topic, and the
# later of the two wins. Documenting `x` as a `slope_bootstrap` object silently
# replaced the generic's description with the one class slope_bootstrap.default()
# refuses, so the help page told the reader to pass exactly the wrong thing.
#' @export
print.slope_bootstrap <- function(x, ...) {
  cat("<slope_bootstrap>\n")
  cat(sprintf("  Statistic:   %s\n", x$statistic))
  cat(sprintf("  Observed:    %s\n", format(x$observed, digits = 6)))
  cat(sprintf("  Replicates:  %d used%s\n", length(x$replicates),
              if (x$n_failed > 0L) sprintf(", %d failed", x$n_failed) else ""))
  cat(sprintf("  %.0f%% %s CI: [%s, %s]\n", 100 * x$level,
              if (identical(x$type, "bca")) "BCa" else "percentile",
              format(x$ci[1L], digits = 6), format(x$ci[2L], digits = 6)))
  invisible(x)
}
