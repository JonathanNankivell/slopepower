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
  # Covariate columns travel under their internal names, already numeric, so
  # make_refitter() can hand them back to slope_params() unchanged. Each
  # replicate recentres them on its own participants, and drops any that a
  # resample leaves constant (a rare factor level not drawn) rather than
  # failing on a singular fit.
  for (cc in grep("^sp_cov_[0-9]+$", names(g), value = TRUE)) out[[cc]] <- g[[cc]]
  out
}

#' Build a closure that refits stage one on a resampled frame
#'
#' The call to `slope_params()` -- formula, `origin = "none"` (the times in the
#' recovered frame have already been re-origined per subject, so this avoids
#' repeating that work and the message that goes with it), which comparator
#' argument to pass, and the random-effects structure to hold fixed -- is
#' identical for every replicate of a given bootstrap; only the resampled frame
#' differs. Assembling it once here, rather than once per replicate, means
#' [run_bootstrap()]'s loop -- run several hundred times, plus once per subject
#' for the BCa jackknife -- reconstructs only what actually changes between
#' replicates.
#'
#' `common_variance` is pinned to `params$common_variance`, which records the
#' structure the *observed* fit ended up with -- whether the caller asked for it
#' or `slope_params()` fell back to it -- so every replicate fits the model the
#' point estimate came from. Left unset, as it was, the argument defaulted to
#' `NULL` on every replicate, and all four combinations went wrong:
#'
#' * `TRUE`: the observed fit used the reduced structure; every replicate tried
#'   the full one. Different model, same interval.
#' * `FALSE`: the caller forbade the reduced structure, and `slope_params()`
#'   errors on it -- but a replicate passing `NULL` fell back to it silently and
#'   was counted as a success, so the interval contained fits the caller had
#'   ruled out.
#' * `NULL` either way: replicates that could not fit the full structure fell
#'   back one at a time, mixing two models within a single interval instead of
#'   being reported as failures.
#'
#' It matters only under `healthy`. Without covariates it matters only through
#' `slope_comparator`: the model factorises per group (see the `common_variance`
#' note in [slope_params()]), so the case estimates are invariant. With
#' covariates the shared coefficients couple the groups and the case estimates
#' can move too. Either way the *controls'* slope is not invariant, and that is
#' what `slope_difference` -- and so every stage-two answer -- is measured
#' against. Balanced complete data hides this entirely, because
#' GLS then coincides with OLS whatever the covariance structure; ragged
#' follow-up does not, and resampling preserves each subject's own visit
#' pattern, so an unbalanced study stays unbalanced in every replicate.
#' @noRd
make_refitter <- function(params) {
  comparator <- params$comparator
  args <- list(formula = y ~ time | subject, data = quote(frame), origin = "none")
  # Refit with the same adjustment the point estimate used; dropping it would
  # bootstrap the unadjusted model and report its interval as the adjusted one.
  b <- names(nlme::fixef(params$fit))
  cov_cols <- unique(regmatches(b, regexpr("sp_cov_[0-9]+", b)))
  if (length(cov_cols)) {
    args$covariates <- stats::as.formula(paste("~", paste(cov_cols, collapse = " + ")),
                                         env = baseenv())
    args$covariate_time <- any(grepl("sp_cov_[0-9]+:sp_time|sp_time:sp_cov_[0-9]+", b))
  }
  # Likewise the residual structure: a replicate fitted with independent
  # residuals would put the wrong model's spread around the point estimate.
  args <- c(args, residual_refit_args(params$residual))
  # And the optimiser settings. An object built before `control` was recorded
  # has none, and was fitted under the default, which is what NULL leaves.
  if (!is.null(params$control)) args$control <- params$control
  cl <- as.call(c(list(quote(slope_params)), args,
                  if (identical(comparator, "healthy")) {
                    # Only under `healthy`: for the other two slope_params()
                    # ignores `common_variance` with a warning, so passing it
                    # there would only ask for a warning -- one per refit,
                    # muffled below -- about an argument that changes nothing.
                    list(comparator = "healthy", group = quote(group),
                         common_variance = isTRUE(params$common_variance))
                  } else if (identical(comparator, "treated")) {
                    list(comparator = "treated", group = quote(group))
                  } else NULL))
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
#' single coefficient. The bootstrap functions compare this to the slope itself
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
#'   cov_intercept_slope = 5, sigma2_residual = 10
#' ))
#'
#' @seealso [slope_sigma()] and [slope_var()], the other quantities computed
#'   from a `slope_params` object; [slope_params_boot()] and its siblings,
#'   which use this to flag an unreliable interval.
#' @export
slope_se <- function(params) {
  context <- "slope_se()"
  check_params(params, context)
  sqrt(slope_vcov(params, context)[1L, 1L])
}

#' Sampling covariance of the slopes, from the fitted model
#'
#' The body of [vcov.slope_params()], and of [slope_se()] -- the square root
#' of its first element -- and [boot_check()], so the slope's standard error
#' and the slope difference's are read off one contrast of the fit's own
#' covariance matrix. Row 1 sums every term of the slope; row 2, when there is
#' a comparator, is its own slope: the first term alone. `NA` throughout when
#' there is no fitted model, or its coefficients cannot be read or named.
#' `params` arrives validated.
#' @noRd
slope_vcov <- function(params, context) {
  nm <- if (identical(params$comparator, "none")) "slope" else c("slope", "slope_comparator")
  out <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
  fit <- params$fit
  if (is.null(fit)) return(out)
  b <- tryCatch(nlme::fixef(fit), error = function(e) NULL)
  V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(b) || is.null(V)) return(out)
  terms <- slope_terms(params, b, context)
  if (is.null(terms)) return(out)
  A <- rbind(as.numeric(names(b) %in% terms),
             if (length(nm) > 1L) as.numeric(names(b) == terms[[1L]]))
  out[] <- A %*% as.matrix(V) %*% t(A)
  out
}

#' The fitted model's coefficient names that sum to the slope, or `NULL`
#'
#' Which terms sum to the slope, by comparator, is params.R's
#' slope_fixef_parts() -- the same mapping slope_params() itself sums the
#' *values* of, via fixef_term(), so the standard error and the slope can never
#' name a different set of coefficients. Every part -- not just an interaction
#' -- is resolved via resolve_fixef_name() rather than assuming a spelling,
#' exactly as fixef_term() does: for a single name it degenerates to a plain
#' lookup. Getting this wrong used to return NA silently, which switched off
#' the section 2.6 warning that is the entire reason [slope_se()] computes the
#' standard error; so a failure warns, and returns `NULL`. In comparator order:
#' the first term, when there are two, is the comparator's own slope.
#' @noRd
slope_terms <- function(params, b, context) {
  terms <- vapply(slope_fixef_parts(params$comparator),
                  function(p) resolve_fixef_name(b, p),
                  character(1L))
  # `anyNA()` alone: resolve_fixef_name() returns either a member of names(b)
  # or NA, so a non-NA element being absent from names(b) is not a second
  # failure mode to test for.
  if (anyNA(terms)) {
    warning(sprintf(paste0(
      "%s: could not identify the slope terms in the fitted model (have %s); ",
      "returning NA. If this call came from a bootstrap, its section 2.6 ",
      "check on the slope-to-standard-error ratio was skipped."),
      context, paste(names(b), collapse = ", ")), call. = FALSE)
    return(NULL)
  }
  terms
}

#' `match.arg()` for `statistic`, with a message that says where the others are
#'
#' Bare `match.arg()` reports only `'arg' should be one of ...`, which names
#' neither the argument nor the function, and a rejected statistic almost
#' always means the other bootstrap was wanted -- `statistic = "power"` asked of
#' the one that solves for the sample size, say. `advice` says which. The
#' matching itself (including the "whole default vector passed through" case)
#' is still `match.arg()`'s; only the error message is replaced.
#' @noRd
match_statistic <- function(statistic, choices, advice, context) {
  matched <- tryCatch(match.arg(statistic, choices), error = function(e) NULL)
  if (is.null(matched)) {
    stop(sprintf("%s: `statistic` must be %s, not %s.\n  %s", context,
                 paste(sQuote(choices), collapse = " or "),
                 sQuote(paste(as.character(statistic), collapse = ", ")), advice),
         call. = FALSE)
  }
  matched
}

#' Read and restore the caller's random number stream
#'
#' `set.seed()` writes to `.Random.seed` in the global environment, so a bare
#' call to it inside a function reseeds the *session*: every subsequent draw the
#' caller makes -- a simulation, a permutation test, a second bootstrap left
#' unseeded on purpose -- is silently determined by whatever seed they passed to
#' make this one result reproducible. Asking for a reproducible answer should not
#' be the thing that makes the rest of a script reproducible too, and the effect
#' leaves no trace: nothing about the returned object records that the stream
#' moved.
#'
#' `current_seed()` returns `NULL` when the caller has not drawn anything yet, in
#' which case there is no `.Random.seed` and R will initialise one from the clock
#' and process id on first use. `restore_seed()` puts that state back by removing
#' the variable rather than writing a placeholder, so the caller's next draw is as
#' unpredictable as it would have been. The existence check inside it is not
#' redundant with the one in `current_seed()`: the two run either side of the
#' bootstrap, and `.Random.seed` is only absent at the second if nothing in
#' between drew -- which, for a seeded call, cannot happen, but the guard costs
#' nothing and makes `restore_seed()` total.
#' @noRd
current_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' @rdname current_seed
#' @noRd
restore_seed <- function(old) {
  g <- globalenv()
  if (is.null(old)) {
    if (exists(".Random.seed", envir = g, inherits = FALSE)) rm(".Random.seed", envir = g)
  } else {
    assign(".Random.seed", old, envir = g)
  }
  invisible(NULL)
}

#' Widen an interval onto the lattice its statistic lives on
#'
#' `n` is the one statistic the bootstrap reports that is not continuous. Every
#' replicate of it is `2 * ceiling(...)`, so the whole distribution sits on the
#' even integers -- but `quantile()` interpolates between order statistics, and
#' a 95% interval for `slpower1` came back as `[519.45, 963.3]`. Neither
#' endpoint is a trial anyone can run.
#'
#' So the endpoints are moved out to the nearest sizes that are: the lower down,
#' the upper up. Outward rather than to the nearest, because rounding must not
#' be able to narrow a confidence interval -- the reported range always contains
#' the interpolated one, and every value in it is a trial that could actually be
#' fielded. Both are already on the lattice when the interpolation happened to
#' land there, so nothing moves in that case.
#'
#' The alternative -- reading the endpoints off with `quantile(type = 1)`, the
#' inverse ECDF, so they are order statistics and therefore lattice points by
#' construction -- is the textbook estimator for a discrete distribution, and it
#' is not usable here. Type 1 is a step function of `p`, and `p` is
#' `(1 - level) / 2`, which for `level = 0.95` is one ulp *above* 0.025 in
#' binary. When `p * B` is a whole number that ulp decides which side of the
#' step the answer falls on: at `R = 40` it moved the lower endpoint from the
#' first order statistic to the second, 452 to 518. `p * B` is a whole number
#' at exactly the replicate counts people choose -- 200, 1000 -- so the failure
#' is not a corner case. Interpolating first and widening afterwards moves the
#' knife edge onto the lattice, where a one-ulp perturbation maps a point to
#' itself.
#'
#' Deliberately keyed on the statistic, through [on_lattice()], rather than on
#' the values: `n`'s discreteness is a property of how [solve_slope()] builds
#' it, not something to be rediscovered per bootstrap from replicates that
#' happen to look integral.
#' @noRd
widen_to_lattice <- function(ci, statistic) {
  if (!on_lattice(statistic)) return(ci)
  c(2 * floor(ci[[1L]] / 2), 2 * ceiling(ci[[2L]] / 2))
}

#' Is this statistic reported on the even-participant lattice?
#'
#' The one place the rule is written: [widen_to_lattice()] widens by it, and
#' [run_bootstrap()] records it as the result's `lattice` flag for the print
#' methods. Two copies could disagree, leaving an interval widened but printed
#' without its per-arm divisor or basis note.
#' @noRd
on_lattice <- function(statistic) identical(statistic, "n")

#' The three arguments every bootstrap driver takes, checked once
#'
#' [run_bootstrap()] and [slope_sample_size_grid_boot()] are two drivers of one
#' resampling scheme and validate its arguments identically. Written twice, the
#' two could accept different things -- a `ci_method` spelling admitted by one and
#' not the other, say -- while claiming in their shared documentation to take
#' the same arguments.
#'
#' Returns the matched `ci_method` rather than validating in place, since
#' [match.arg()] is the one of the three that has a value to give back.
#' @noRd
check_boot_args <- function(R, ci_method, level, context) {
  check_whole_number(R, "R", "replicates", context, lower = 1)
  check_probability(level, "level", context)
  match.arg(ci_method, c("bca", "percentile"))
}

#' The interval's tail probabilities
#'
#' One spelling, because the exact floating-point value matters: see
#' [widen_to_lattice()] on how one ulp in `(1 - level) / 2` decides which side
#' of a `quantile()` step an endpoint falls on. Two drivers computing it in two
#' places is two things to keep in step with that reasoning.
#' @noRd
boot_probs <- function(level) c((1 - level) / 2, 1 - (1 - level) / 2)

#' Seed a bootstrap, and say what the caller must put back
#'
#' The seeding policy both drivers follow, in one place: a seeded call reads
#' the caller's stream, seeds, and hands back what to restore; `seed = NULL`
#' touches nothing and returns `NULL`, so an unseeded bootstrap draws from and
#' advances the stream as any other RNG-using function does -- back-to-back
#' unseeded calls must give different answers.
#'
#' The `on.exit()` that consumes the return value stays at the call site rather
#' than being registered from in here, because it is the *caller's* frame it has
#' to fire on: a driver can leave by several routes -- a "not enough replicates
#' succeeded" stop, an error no `tryCatch()` covers, a user interrupt part-way
#' through several hundred fits -- and all of them must leave the stream as they
#' found it. Registering it in a parent frame from here would work, but only by
#' reaching into a frame this function does not own; two readable lines at each
#' driver are worth more than one clever one here.
#' Note that `NULL` comes back in two different situations, and the caller must
#' not confuse them: `seed = NULL`, where there is nothing to put back because
#' nothing was disturbed, and a seeded call made in a session that had no
#' `.Random.seed` yet, where `NULL` is precisely what has to be restored --
#' [restore_seed()] reads it as "remove the variable `set.seed()` just created".
#' So the `on.exit()` at each call site is gated on `seed`, never on this
#' return value.
#' @return The stream to restore; see the note above on its `NULL`.
#' @noRd
seed_bootstrap <- function(seed) {
  if (is.null(seed)) return(NULL)
  old_seed <- current_seed()
  set.seed(seed)
  old_seed
}

#' One jackknife pass, taken on first use and kept
#'
#' Both drivers want the same laziness -- a bootstrap whose bias correction is
#' degenerate refits nothing, and one that needs several columns pays for the
#' refits once -- and both append the replicate slope as the *last* accessor so
#' that it can be read alongside the statistics proper. That convention used to
#' be re-expressed as a literal index at each read site, where getting it wrong
#' yields a plausible acceleration rather than an error; here `slope_col()`
#' knows the position because it is the function that chose it.
#' @param setup As returned by [boot_setup()].
#' @param computes The statistics' accessors, without the slope.
#' @return `list(col = function(k), slope_col = function())`, both reading the
#'   one memoised matrix.
#' @noRd
lazy_jackknife <- function(setup, computes) {
  jack <- NULL
  matrix_of <- function() {
    if (is.null(jack)) {
      jack <<- jackknife_values(setup$frame, setup$subject_index, setup$refitter,
                                c(computes, list(function(p) p$slope)))
    }
    jack
  }
  list(col = function(k) matrix_of()[, k],
       slope_col = function() matrix_of()[, length(computes) + 1L])
}

#' The summary of the refitted slopes every bootstrap result carries
#'
#' The slope is what the resampling actually perturbs, and every interval a
#' bootstrap reports is a function of it, so both drivers carry the same six
#' fields -- `run_bootstrap()` as list elements, the grid as attributes. Two
#' hand-maintained field lists is the drift [stage_two_result()] (power.R) was
#' extracted to remove for the result classes; this is the same extraction for
#' the slope block.
#'
#' `slope_int` is passed in rather than built here: under
#' `statistic = "slope"` [run_bootstrap()] has already built exactly this
#' interval as its main one and must reuse it, since a second [boot_interval()]
#' pass would warn twice about one failure.
#'
#' On `straddle`: this is section 2.6's sign-straddling hazard measured on the
#' replicates actually drawn, rather than the analytic standard error that
#' stands proxy for it. Measured on the replicate *slopes*, against the fitted
#' slope -- never on a statistic against its own observed value. What the paper
#' (p.583) warns about is bootstrap samples yielding "estimates of the mean
#' slope that are both negative and positive", so that some replicates describe
#' a trial reducing a rising slope and others one reducing a falling slope.
#' Comparing a statistic's own sign answers a different question, and for `n`
#' and `power` it is not a question at all: a sample size is a positive integer
#' and a power lies in [0, 1], so that comparison was identically 0 and the
#' check silently never fired.
#' @noRd
slope_replicate_summary <- function(observed_slope, good_slopes, slope_int, check,
                                    good_checks) {
  list(straddle = mean(sign(good_checks) != sign(check$observed)),
       straddle_of = check$label,
       slope_observed = observed_slope,
       slope_replicates = good_slopes,
       slope_mean = mean(good_slopes),
       slope_sd = stats::sd(good_slopes),
       slope_ci = slope_int$ci,
       slope_ci_method = slope_int$ci_method)
}

#' Resample setup shared by every replicate of a bootstrap
#'
#' The section 2.6 standard-error check, the recovered modelling frame, the
#' per-subject row index, the stratification groups and the refitting closure
#' all depend only on `params` -- never on what is being computed from a
#' replicate, or on how many things are. Lifted out of [run_bootstrap()] so
#' that [slope_sample_size_grid_boot()] can build this once for a whole grid
#' rather than once per cell: the resampling scheme does not know a grid cell
#' exists.
#' @noRd
boot_setup <- function(params, context, target = NULL) {
  # Section 2.6: if the slope is not large relative to its standard error the
  # replicates can straddle zero, and an interval for the sample size stops
  # meaning anything. Measured on whatever stage two divides by -- the slope
  # difference, when the target is a comparator's slope; see boot_check().
  check <- boot_check(params, target, context)
  check$observed <- observed <- check$of(params)
  se <- check$se
  if (is.finite(se) && se > 0 && abs(observed) / se < 2.5) {
    warning(sprintf(paste0("%s: the estimated %s (%.4g) is only %.2f times its standard error ",
                           "(%.4g). Nash et al. (2021, section 2.6) suggest a ratio above 2.5; ",
                           "below it, bootstrap replicates may include %ss of both signs and ",
                           "the resulting interval can be meaningless."),
                    context, check$label, observed, abs(observed) / se, se, check$label),
            call. = FALSE)
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

  list(frame = frame, subject_index = subject_index, groups = groups,
      refitter = make_refitter(params), check = check)
}

#' What became of a replicate that refitted but could not be solved, in the
#' run-time warning and the printed note alike
#' @noRd
BOOT_SOLVE_FAILED <- paste("refitted, but the calculation could not be solved on them",
                           "(for instance, a slope difference of zero)")

#' The warning text for discarded replicates, split by cause
#'
#' `n_failed` counts every replicate that yielded no statistic; `n_refit` those
#' among them whose stage-one refit failed, so the rest refitted and failed in
#' the stage-two calculation itself.
#' @noRd
boot_failure_text <- function(n_refit, n_failed, R) {
  parts <- c(
    if (n_refit > 0L) sprintf("%d of %d replicates failed to converge", n_refit, R),
    if (n_failed > n_refit) sprintf("%d of %d %s", n_failed - n_refit, R, BOOT_SOLVE_FAILED))
  paste0(paste(parts, collapse = "; "), "; all were discarded.")
}

#' What the section 2.6 check is measured on
#'
#' The paper's hazard is replicates whose slope changes sign, because a trial
#' powered to reduce a falling slope and one powered to reduce a rising slope
#' are different trials. What stage two actually divides by, though, is the
#' slope *difference*: the slope less the reference it is moved toward. That is
#' the slope itself when the reference is zero, but under
#' `comparator = "healthy"` or `target = "observed"` it is the slope less the
#' comparator's (see [reference_is_comparator()]). There a case slope of -1.0
#' against a control slope of -0.95 is far from zero while the difference
#' that sizes the trial is not, and checking the slope alone would pass a
#' meaningless interval. So the check follows the difference wherever the
#' difference is what is used.
#'
#' `target = NULL` is [slope_params_boot()], which bootstraps the slope itself
#' and is checked on it whatever the comparator.
#'
#' @return A list: `of`, reading the checked quantity off a `slope_params`
#'   object (the fit or a refit); `se`, its analytic standard error from the
#'   fitted model, `NA` without one; and `label`, its name for messages.
#'   [boot_setup()] adds `observed`, the quantity on the fit itself.
#' @noRd
boot_check <- function(params, target, context) {
  V <- slope_vcov(params, context)
  if (is.null(target) || !reference_is_comparator(params$comparator, target)) {
    return(list(of = function(p) p$slope, se = sqrt(V[1L, 1L]), label = "slope"))
  }
  list(of = function(p) p$slope - p$slope_comparator,
       se = sqrt(V[1L, 1L] + V[2L, 2L] - 2 * V[1L, 2L]), label = "slope difference")
}

#' Refit `R` resampled replicates, reading one or more statistics off each
#'
#' The loop [run_bootstrap()] used to run for its one `compute`, generalised to
#' a list of them: `jackknife_values()` already reads several accessors off one
#' refit for the price of one jackknife pass, and a bootstrap grid needs the
#' same thing here, for the same reason -- one refit priced out at every design
#' rather than one refit per cell. A single-element `computes` reproduces
#' [run_bootstrap()]'s own loop exactly, column for column.
#'
#' Warnings are suppressed for every replicate. Anything worth saying about the
#' calculation -- a target effect that makes the slope more extreme, say -- is a
#' property of the data, and the caller has already heard it once from the
#' stage-two call that built the object; repeating it several hundred times
#' here, attributed to a function they did not call, tells them nothing new.
#' Suppressed once around the whole loop rather than per replicate: messages
#' (the progress ticks) are a different condition class and pass through
#' regardless.
#'
#' The refit is caught separately from each `compute` so that the replicate's
#' own slope can be read off it: that -- not any one statistic -- is what a
#' sign-straddling check has to be measured on. A refit failure and a
#' `compute` failure both land the whole row (or that one column) on `NA`, so
#' which of the two failed does not change the accounting.
#'
#' @param setup As returned by [boot_setup()].
#' @param computes A list of `function(p) numeric(1)` closures, each reading
#'   one statistic off a refit. Failing independently: one column's `NA` does
#'   not cost the others their value for that replicate.
#' @return A list with `replicates`, an `R` x `length(computes)` matrix, and
#'   `slopes`, the length-`R` vector of each replicate's own refitted slope.
#' @noRd
boot_replicate_matrix <- function(setup, computes, R, progress, context) {
  frame <- setup$frame
  subject_index <- setup$subject_index
  groups <- setup$groups
  refitter <- setup$refitter

  replicates <- matrix(NA_real_, nrow = R, ncol = length(computes))
  slopes <- rep(NA_real_, R)
  checks <- rep(NA_real_, R)
  tick <- max(1L, floor(R / 10))
  suppressWarnings(for (b in seq_len(R)) {
    p <- tryCatch(refitter(resample_frame(frame, subject_index, groups)),
                  error = function(e) NULL)
    if (!is.null(p)) {
      slopes[b] <- p$slope
      checks[b] <- setup$check$of(p)
      replicates[b, ] <- eval_computes(computes, p)
    }
    if (isTRUE(progress) && b %% tick == 0L) {
      message(sprintf("%s: %d of %d replicates", context, b, R))
    }
  })

  list(replicates = replicates, slopes = slopes, checks = checks)
}

#' One statistic's interval, BCa where possible and percentile as its fallback
#'
#' The body of the `interval()` closure [run_bootstrap()] used to build inline,
#' lifted out so a bootstrap grid can build one interval per cell without
#' repeating the fallback logic. `jack_col` is a zero-argument function rather
#' than the column itself, so the jackknife it reads from stays lazy: under
#' `ci_method = "percentile"` it is never called, and the jackknife is never taken.
#'
#' The percentile quantile is computed only where it is actually used: as the
#' answer itself under `ci_method = "percentile"`, or as the fallback when a BCa
#' interval could not be built. A successful BCa interval never touches it, so
#' it is no longer computed -- and immediately discarded -- on every such call.
#' @noRd
boot_interval <- function(theta, observed, jack_col, ci_method, probs, context, what) {
  percentile <- function(theta) {
    list(ci = stats::quantile(theta, probs, names = FALSE, type = 7),
        ci_method = "percentile")
  }
  if (!identical(ci_method, "bca")) return(percentile(theta))
  bca <- bca_from_jack(theta, observed, jack_col(), probs)
  # A character return is the reason it could not be built. Previously either
  # failure was reported as "every replicate falls on one side of the observed
  # value", including the ones where the bias correction was fine and it was
  # the jackknife that could not be had -- so the warning named a cause the
  # caller could not act on and, in that case, was simply untrue.
  if (is.character(bca)) {
    warning(sprintf("%s: %s%s; reporting a percentile interval instead.",
                    context, bca, what), call. = FALSE)
    return(percentile(theta))
  }
  list(ci = bca, ci_method = "bca")
}

#' Bootstrap a scalar computed from resampled stage-one parameters
#'
#' The resampling scheme, the section 2.6 check and the interval construction are
#' the same whatever is being bootstrapped; only `compute` differs, and each
#' exported driver supplies one closure that reads its statistic off a refit.
#' The arguments arrive validated -- `ci_method` matched, `per_arm` checked --
#' by the driver, which has to check them before it computes `observed` anyway.
#' `cls` is the driver's own class, put in front of the shared
#' `"slope_bootstrap"` that the print method dispatches on.
#' @noRd
run_bootstrap <- function(params, compute, observed, statistic, R, ci_method, level,
                          seed, progress, per_arm, context, cls, target = NULL) {
  old_seed <- seed_bootstrap(seed)
  if (!is.null(seed)) on.exit(restore_seed(old_seed), add = TRUE)

  # `observed` is passed in rather than recomputed here. It is the same number
  # either way -- `compute` on the original parameters reproduces it -- but the
  # driver has already solved that exact calculation once, warnings and all,
  # and solving it again would repeat them.
  setup <- boot_setup(params, context, target)
  se <- setup$check$se

  mat <- boot_replicate_matrix(setup, list(compute), R, progress, context)
  replicates <- mat$replicates[, 1L]
  slopes <- mat$slopes

  failed <- is.na(replicates)
  n_failed <- sum(failed)
  good <- replicates[!failed]
  # No NAs to guard against: `slopes[b]` is filled whenever the refit that
  # `replicates[b]` also depends on succeeded, and new_slope_params() admits
  # only a finite slope.
  #
  # Under `statistic == "slope"`, `compute` above is `function(p) p$slope` --
  # the same read as `slopes[b] <- p$slope` -- so `good` and `slopes[!failed]`
  # are the same values under two names. Aliased here, and in the summary
  # statistics built from it below, rather than recomputed independently, so
  # the two can never drift apart the way two separately-taken means could.
  good_slopes <- if (identical(statistic, "slope")) good else slopes[!failed]
  if (length(good) < 2L) {
    stop(sprintf("%s: %d of %d replicates failed; not enough succeeded to form an interval.",
                 context, n_failed, R), call. = FALSE)
  }
  # Two different failures, reported apart: a refit that did not converge is
  # the optimiser's (and slope_lme_control()'s) business, while a refit whose
  # stage-two calculation failed -- a zero slope difference, say -- is the
  # data's. One count for both pointed the reader at the wrong one.
  n_refit_failed <- sum(is.na(slopes))
  if (n_failed > 0L) {
    warning(sprintf("%s: %s", context, boot_failure_text(n_refit_failed, n_failed, R)),
            call. = FALSE)
  }

  probs <- boot_probs(level)

  # One jackknife pass serves both intervals; lazy_jackknife() owns both the
  # memoisation and the convention that the slope accessor is appended last.
  jack <- lazy_jackknife(setup, list(compute))

  main <- boot_interval(good, observed, function() jack$col(1L), ci_method, probs, context,
                        "")
  # The replicate slopes get an interval of their own, so that the printed table
  # can show what the resampling did to the quantity every other row is derived
  # from -- a sample size two-thirds wider than its point estimate means one
  # thing if the slope barely moved and another if the slope moved with it.
  # Under `statistic = "slope"` the two are the same computation on the same
  # numbers, so it is reused rather than repeated: a second pass would warn
  # twice about one failure.
  slope_int <- if (identical(statistic, "slope")) {
    main
  } else {
    boot_interval(good_slopes, params$slope, jack$slope_col, ci_method, probs, context,
                 " for the replicate slopes")
  }

  # After the choice, so a percentile interval and a BCa one are reported on the
  # same lattice; before the object is built, so nothing downstream has to know
  # which statistic needs it. The slope is continuous and needs none of this.
  ci <- widen_to_lattice(main$ci, statistic)
  used_method <- main$ci_method
  # Decided once, here, and carried on the object as `lattice` rather than
  # re-tested against `statistic` wherever the choice matters again -- print
  # methods included. A discrete statistic added later, or a print-side
  # refactor that misses one of several such tests, cannot silently disagree
  # with the widening this function already did -- on_lattice() being the one
  # rule both read.
  lattice <- on_lattice(statistic)

  # The slope block -- `straddle` included -- is built by the same helper the
  # bootstrapped grid uses, so the two results describe the resampling in the
  # same terms. Under `statistic == "slope"` its `slope_mean`/`slope_sd` are
  # `boot_mean`/`boot_sd` recomputed on the identical vector: `compute` is then
  # `function(p) p$slope`, the same read as `slopes[b] <- p$slope`, so `good`
  # and `good_slopes` are one vector under two names and the two means cannot
  # differ. Recomputed rather than aliased, so the helper owes its caller
  # nothing about which statistic was asked for.
  slope_block <- slope_replicate_summary(params$slope, good_slopes, slope_int, setup$check,
                                         mat$checks[!failed])

  structure(c(list(observed = observed, replicates = good, ci = ci,
                   ci_method = used_method, statistic = statistic, R = R,
                   n_failed = n_failed, n_refit_failed = n_refit_failed, level = level,
                   se = se,
                   boot_mean = mean(good), boot_sd = stats::sd(good),
                   straddle = slope_block$straddle, straddle_of = slope_block$straddle_of,
                   lattice = lattice),
              slope_block[setdiff(names(slope_block), c("straddle", "straddle_of"))]),
            class = c(cls, "slope_bootstrap"), per_arm = per_arm)
}

#' Bootstrap a sample size or power calculation
#'
#' Resamples subjects with replacement, refits the stage-one mixed model on each
#' replicate, and re-solves the stage-two calculation against each refit. This
#' propagates the estimation uncertainty in the slope and variance components
#' through to the sample size or power, as recommended in section 2.6 of Nash
#' et al. (2021).
#'
#' Each function takes exactly the arguments of the calculation it bootstraps
#' --- [slope_sample_size()] or [slope_power()] --- followed by the bootstrap's
#' own, so a calculation is turned into its interval by changing the function
#' name:
#'
#' ```
#' slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33)
#' slope_sample_size_boot(pars, c(0, 1, 2), effectiveness = 0.33, R = 999)
#' ```
#'
#' `statistic` offers only the quantities the calculation solves for or
#' derives --- never one of its own inputs, whose interval would be
#' zero-width. To bootstrap the fitted slope alone, with no trial design, use
#' [slope_params_boot()].
#'
#' Subjects, not observations, are the sampling unit, and each drawn subject is
#' given a fresh identifier so that a subject selected twice is treated as two
#' people rather than as one person with twice as much data. When the parameters
#' came from two-group data, resampling is stratified so that each replicate has
#' the same number of cases and controls (or treated and control subjects) as the
#' original.
#'
#' Every replicate refits the random-effects structure the original fit ended up
#' with, so a `healthy` study whose controls were fitted with a random intercept
#' only --- because `common_variance = TRUE` asked for it, or because the full
#' structure failed to converge --- is bootstrapped under that same structure
#' throughout. A replicate that cannot fit it is discarded and counted in
#' `n_failed` rather than quietly refitted under the other one, so an interval
#' never mixes the two.
#'
#' A bootstrapped sample size is reported as sample sizes. Every replicate of
#' `n` is an even integer, so the interval's endpoints are widened out to the
#' nearest sizes a trial could actually be fielded at rather than left between
#' them, and the bias correction behind a BCa interval counts replicates tied
#' with the observed size as half below rather than wholly below. The other
#' statistics are continuous and neither applies to them.
#'
#' `R` defaults to 999. A confidence interval needs an order of magnitude more
#' replicates than a standard error does, because its endpoints are read from
#' the tails: at `R = 199` a 2.5% tail holds about five replicates, and the
#' interval wobbles with the seed rather than with the data. Bootstrapping
#' `slpower1`'s sample size six times over, changing only the seed, the reported
#' bounds spanned 504--532 and 1021--1061 at `R = 199`, against 504--512 and
#' 1027--1052 at `R = 999` --- tens of participants of pure Monte Carlo noise on
#' a number the caller is meant to plan a trial around. The paper's own worked
#' bootstrap uses 2000. The odd number is the usual convention, `R + 1` round
#' rather than `R`.
#'
#' That costs real time, because every replicate is a mixed-model fit and
#' `ci_method = "bca"` adds a leave-one-subject-out jackknife on top --- one
#' further fit per subject. On the paper's `slpower1`, 200 participants, a
#' default BCa bootstrap is roughly a minute and a half; on `slpower3`, whose
#' model is slower to fit, nearer two and a half. Pass a small `R` while setting
#' a calculation up, and leave the default for the answer you intend to report.
#'
#' @inheritParams slope_sample_size
#' @param params A `slope_params` object from [slope_params()]. Objects from
#'   [slope_params_manual()] carry no data to resample and are refused.
#' @param statistic Which quantity to bootstrap: the one the calculation solves
#'   for (the default) or the target treatment effect `tte` behind it.
#' @param R Number of bootstrap replicates. The default is sized for the
#'   interval this returns rather than for speed; see above before lowering it.
#' @param ci_method How the confidence interval is built: `"bca"` (the default)
#'   for bias-corrected and accelerated, or `"percentile"`. The paper
#'   recommends BCa because the distribution of estimated sample sizes is
#'   typically skewed.
#' @param level Confidence level for the interval.
#' @param seed Optional integer seed, for reproducibility. The caller's random
#'   number stream is restored afterwards, so a seeded call reproduces its own
#'   result without reseeding the session; leaving it `NULL` draws from, and
#'   advances, the stream as usual.
#' @param progress Report progress while resampling.
#' @param per_arm For `statistic = "n"`, which basis the printed result is on:
#'   `TRUE` (the default) for participants per arm, `FALSE` for the trial
#'   total. Display only -- `replicates`, `ci`, `boot_mean` and `boot_sd` are
#'   always the totals `slope_sample_size()` itself returns; `per_arm` is
#'   recorded as an attribute (`attr(x, "per_arm")`) and consulted by
#'   `print()`, which halves them for display and can be overridden
#'   afterwards with `print(x, per_arm = FALSE)`. Ignored, but still accepted,
#'   for any other `statistic`, none of which have arms.
#'
#' @return An object of class `c("slope_sample_size_boot", "slope_bootstrap")`
#'   (or `"slope_power_boot"`, or `"slope_params_boot"`), with elements
#'   `observed`, `replicates`, `ci`, `ci_method` (the method actually used,
#'   which falls back to `"percentile"` where BCa could not be built),
#'   `statistic`, `R`, `n_failed` (replicates that yielded no statistic, for
#'   any reason), `n_refit_failed` (those among them whose model refit failed
#'   to converge; the rest refitted, but the calculation could not be solved
#'   on them), `level`, `se`, `boot_mean` and `boot_sd`
#'   (the mean and SD of `replicates`), `straddle` (the proportion of retained
#'   replicates whose *refitted slope* has a different sign from the fitted
#'   slope, whatever `statistic` was asked for --- the section 2.6 hazard
#'   itself, rather than the analytic proxy for it in `se`), `straddle_of`
#'   (`"slope"`, or `"slope difference"` when the target treatment effect is
#'   measured toward a comparator's slope --- under `comparator = "healthy"`,
#'   or `target = "observed"` --- since that difference, not the slope, is then
#'   what the sample size depends on; `straddle` and `se` are measured on it),
#'   and `lattice`
#'   (whether `statistic` is `"n"`, decided once here and read by the print
#'   method rather than re-tested there). The `per_arm` argument is recorded
#'   as an attribute, not an element, since it changes nothing about these
#'   values.
#'
#'   The refitted slopes are summarised alongside, whatever `statistic` was
#'   asked for, in `slope_observed` (the fitted slope), `slope_replicates`,
#'   `slope_mean`, `slope_sd`, `slope_ci` and `slope_ci_method`. They are what
#'   the resampling actually perturbs --- every replicate of every other
#'   statistic is a function of them --- so an interval on the statistic is
#'   only as meaningful as the one on the slope beneath it. `slope_ci` is built
#'   the same way as `ci` and from the same jackknife, so it costs no extra
#'   model fits, but the two can fall back independently: compare
#'   `slope_ci_method` with `ci_method`. The print method does not show them
#'   --- it reports the statistic asked for, and `straddle` beneath it --- so
#'   these six fields are where a reader who wants the slope's own interval
#'   finds it.
#'
#'   Every summary is over the retained replicates --- the ones whose refit and
#'   whose statistic both succeeded; `n_failed` counts the rest, which are
#'   discarded rather than imputed.
#'
#'   For `statistic = "n"` each replicate is a whole trial: [slope_sample_size()]
#'   is re-solved on that replicate's refitted parameters, and it rounds *up* to
#'   a whole participant per arm before doubling for the total, so every value in
#'   `replicates` is an even integer. `boot_mean` and `boot_sd` therefore
#'   summarise sizes that were each rounded up individually; they do not round an
#'   average, and `boot_mean` is not itself a size a trial could be run at. The
#'   slope carries no such rounding.
#'
#' @references
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#' \doi{10.1177/1536867X211045512}
#'
#' @examples
#' # The parameters must come from slope_params(): bootstrapping resamples the
#' # underlying subjects, which slope_params_manual() objects do not carry.
#'
#' # No comparator: all two hundred participants of `slpower1`.
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#'
#' # One mixed-model fit per replicate, so a real run wants a much larger R.
#' \donttest{
#' slope_sample_size_boot(pars, c(0, 1, 2), effectiveness = 0.33, R = 100, seed = 42)
#'
#' # The same calculation also carries the target treatment effect.
#' slope_sample_size_boot(pars, c(0, 1, 2), effectiveness = 0.33,
#'                        statistic = "tte", R = 100, seed = 42)
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
#' pars2 <- slope_params(sdmt ~ visit | id, data = sim2,
#'                       comparator = "healthy", group = case)
#'
#' \donttest{
#' # Bootstrapping the power of a fixed sample size, rather than the size itself.
#' slope_power_boot(pars2, c(0, 1, 2), n = 400, effectiveness = 0.33,
#'                  R = 100, seed = 42)
#' }
#'
#' # Randomised-trial comparator, target = "observed": all one hundred and
#' # fifty participants of `slpower3`, bootstrapping the sample size needed
#' # to detect the effect the trial actually found.
#' pars3 <- slope_params(sdmt ~ visit | id, data = slpower3,
#'                       comparator = "treated", group = treat)
#'
#' \donttest{
#' slope_sample_size_boot(pars3, c(0, 0.5, 2), target = "observed",
#'                        R = 100, seed = 42)
#' }
#'
#' @seealso [slope_sample_size()] and [slope_power()], the calculations
#'   bootstrapped; [slope_params_boot()] for the fitted slope alone;
#'   [slope_se()] for the standard error behind the section 2.6 check;
#'   [slope_sample_size_floor_boot()] and [slope_power_ceiling_boot()] for the
#'   bounds over all designs;
#'   [slope_sample_size_grid_boot()] and [slope_power_grid_boot()] to bootstrap
#'   every cell of a grid at once, sharing one set of replicates
#' @export
slope_sample_size_boot <- function(params, visits, dropout = NULL,
                                   dropout_scale = c("incremental", "cumulative"),
                                   power = 0.8, effectiveness = 0.25,
                                   target = c("effectiveness", "observed"),
                                   alpha = 0.05, per_arm = TRUE,
                                   statistic = c("n", "tte"),
                                   R = 999, ci_method = c("bca", "percentile"),
                                   level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_sample_size_boot()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  statistic <- match_statistic(statistic, c("n", "tte"), paste0(
    "This bootstrap solves for the sample size, so `n` and the target treatment\n",
    "  effect `tte` behind it are what it can offer. For the power a fixed\n",
    "  sample size achieves, use slope_power_boot() instead."), context)
  ci_method <- check_boot_args(R, ci_method, level, context)
  per_arm <- check_per_arm(per_arm, context)
  design <- stage_two_design(params, visits, dropout, match.arg(dropout_scale), context)
  # `power` is the target the calculation is solved to, and so is an input held
  # fixed across replicates -- it is what makes `n` vary.
  x <- sample_size_result(params, design, power, effectiveness, target, alpha,
                          per_arm, context)
  bootstrap_stage_two(x, sample_size_result, "power", statistic, R, ci_method, level,
                      seed, progress, per_arm, context, "slope_sample_size_boot")
}

#' @rdname slope_sample_size_boot
#' @param n Total number of participants across both arms, as in
#'   [slope_power()]. Required.
#' @export
slope_power_boot <- function(params, visits, dropout = NULL,
                             dropout_scale = c("incremental", "cumulative"),
                             n, effectiveness = 0.25,
                             target = c("effectiveness", "observed"),
                             alpha = 0.05, per_arm = TRUE,
                             statistic = c("power", "tte"),
                             R = 999, ci_method = c("bca", "percentile"),
                             level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_power_boot()"
  # `is.null(n)` too; see the note on the same guard in slope_power().
  if (missing(n) || is.null(n)) {
    stop(sprintf(paste0(
      "%s: `n` is required -- it is the sample size whose power is being\n",
      "  bootstrapped. For an interval around the sample size a target power\n",
      "  needs, use slope_sample_size_boot()."), context), call. = FALSE)
  }
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  statistic <- match_statistic(statistic, c("power", "tte"), paste0(
    "This bootstrap solves for the power a fixed sample size achieves, so\n",
    "  `power` and the target treatment effect `tte` behind it are what it can\n",
    "  offer. For the sample size a target power needs, use\n",
    "  slope_sample_size_boot() instead."), context)
  ci_method <- check_boot_args(R, ci_method, level, context)
  per_arm <- check_per_arm(per_arm, context)
  design <- stage_two_design(params, visits, dropout, match.arg(dropout_scale), context)
  x <- power_result(params, design, n, effectiveness, target, alpha, per_arm, context)
  # `x$n` rather than the `n` supplied: the even number actually used, so the
  # replicates answer the question the observed value answered.
  bootstrap_stage_two(x, power_result, "n", statistic, R, ci_method, level,
                      seed, progress, per_arm, context, "slope_power_boot")
}

#' Bootstrap the sample-size floor or the power ceiling
#'
#' Intervals around the bounds of [slope_sample_size_floor()] and
#' [slope_power_ceiling()], by the subject-level resampling
#' [slope_sample_size_boot()] uses: each replicate refits stage one and
#' re-solves the bound against the refitted variance components. The bounds
#' depend only on `sigma2_slope`, `cov_intercept_slope` and `sigma2_intercept` (through
#' [slope_var_floor()]) and on the slope difference, so their uncertainty is
#' the uncertainty in those --- often a wider interval, relative to the point
#' estimate, than a stated design's, because the floor rests entirely on the
#' between-participant slope variance that a design's residual noise would
#' otherwise dilute.
#'
#' Each function takes exactly the arguments of the bound it bootstraps,
#' followed by the bootstrap's own, as [slope_sample_size_boot()] and
#' [slope_power_boot()] do for theirs; see those for the resampling scheme, the
#' choice of `R`, and how a bootstrapped `n` is reported.
#'
#' @inheritParams slope_sample_size_boot
#' @inheritParams slope_sample_size_floor
#' @param n For `slope_power_ceiling_boot()`: total number of participants, as
#'   in [slope_power_ceiling()]. Required.
#' @param statistic Which quantity to bootstrap: the bound (the default --
#'   `"n"` for the floor, `"power"` for the ceiling) or the target treatment
#'   effect `tte` behind it.
#'
#' @return An object of class `c("slope_sample_size_floor_boot",
#'   "slope_bootstrap")` or `c("slope_power_ceiling_boot", "slope_bootstrap")`,
#'   with the elements described under \sQuote{Value} in
#'   [slope_sample_size_boot()].
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' \donttest{
#' slope_sample_size_floor_boot(pars, effectiveness = 0.33, R = 100, seed = 42)
#' slope_power_ceiling_boot(pars, n = 200, effectiveness = 0.33, R = 100, seed = 42)
#' }
#'
#' @seealso [slope_sample_size_floor()] and [slope_power_ceiling()], the bounds
#'   bootstrapped; [slope_sample_size_boot()] and [slope_power_boot()] for a
#'   stated design.
#' @export
slope_sample_size_floor_boot <- function(params, power = 0.8, effectiveness = 0.25,
                                         target = c("effectiveness", "observed"),
                                         alpha = 0.05, per_arm = TRUE,
                                         statistic = c("n", "tte"),
                                         R = 999, ci_method = c("bca", "percentile"),
                                         level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_sample_size_floor_boot()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  statistic <- match_statistic(statistic, c("n", "tte"), paste0(
    "This bootstrap solves for the smallest sample size any design could need,\n",
    "  so `n` and the target treatment effect `tte` behind it are what it can\n",
    "  offer. For the highest power a fixed sample size could reach, use\n",
    "  slope_power_ceiling_boot() instead."), context)
  ci_method <- check_boot_args(R, ci_method, level, context)
  per_arm <- check_per_arm(per_arm, context)
  x <- floor_result(params, effectiveness, target, alpha, per_arm, context, power = power)
  bootstrap_stage_two(x, floor_result, "power", statistic, R, ci_method, level, seed,
                      progress, per_arm, context, "slope_sample_size_floor_boot")
}

#' @rdname slope_sample_size_floor_boot
#' @export
slope_power_ceiling_boot <- function(params, n, effectiveness = 0.25,
                                     target = c("effectiveness", "observed"),
                                     alpha = 0.05, per_arm = TRUE,
                                     statistic = c("power", "tte"),
                                     R = 999, ci_method = c("bca", "percentile"),
                                     level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_power_ceiling_boot()"
  # `is.null(n)` too; see the note on the same guard in slope_power().
  if (missing(n) || is.null(n)) {
    stop(sprintf(paste0(
      "%s: `n` is required -- it is the sample size whose highest achievable\n",
      "  power is being bootstrapped. For an interval around the smallest sample\n",
      "  size any design could need, use slope_sample_size_floor_boot()."), context),
      call. = FALSE)
  }
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  statistic <- match_statistic(statistic, c("power", "tte"), paste0(
    "This bootstrap solves for the highest power a fixed sample size could\n",
    "  reach, so `power` and the target treatment effect `tte` behind it are\n",
    "  what it can offer. For the smallest sample size any design could need,\n",
    "  use slope_sample_size_floor_boot() instead."), context)
  ci_method <- check_boot_args(R, ci_method, level, context)
  per_arm <- check_per_arm(per_arm, context)
  x <- floor_result(params, effectiveness, target, alpha, per_arm, context, n = n)
  # `x$n` rather than the `n` supplied: the even number actually used, as in
  # slope_power_boot().
  bootstrap_stage_two(x, floor_result, "n", statistic, R, ci_method, level, seed,
                      progress, per_arm, context, "slope_power_ceiling_boot")
}

#' Bootstrap the fitted slope
#'
#' An interval for the slope [slope_params()] estimated, by the same
#' subject-level resampling [slope_sample_size_boot()] uses --- stratified by
#' group, refitting the same model on every replicate --- but with no trial
#' design to re-solve against each refit. The interval for the slope is also
#' reported, alongside the main statistic, by every other bootstrap function;
#' this one is for when the slope is all that is wanted.
#'
#' @inheritParams slope_sample_size_boot
#'
#' @return An object of class `c("slope_params_boot", "slope_bootstrap")`, with
#'   the elements described under \sQuote{Value} in [slope_sample_size_boot()]
#'   and `statistic = "slope"`. The slope fields and the main ones then describe
#'   one calculation, and agree by construction.
#'
#' @examples
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' \donttest{
#' slope_params_boot(pars, R = 100, seed = 42)
#' }
#'
#' @seealso [slope_se()], the analytic standard error of the same slope;
#'   [slope_sample_size_boot()] and [slope_power_boot()] for intervals around a
#'   trial's sample size or power.
#' @export
slope_params_boot <- function(params, R = 999, ci_method = c("bca", "percentile"),
                              level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_params_boot()"
  check_params(params, context)
  ci_method <- check_boot_args(R, ci_method, level, context)
  # A slope has no arms, so the `per_arm` every result carries is fixed: it has
  # no effect, since `lattice` is FALSE for `statistic = "slope"` regardless.
  run_bootstrap(params, function(p) p$slope, params$slope, "slope", R, ci_method,
                level, seed, progress, per_arm = TRUE, context, "slope_params_boot")
}

#' Build one replicate-statistic closure in an environment of its own
#'
#' The things `compute` reads, and nothing else. See the note at its call site
#' on why an inline `function(p)` would retain the whole stage-two result --
#' fitted model included -- no matter what was copied out of it first.
#' @noRd
boot_stage_two_compute <- function(fn, slim, fixed, statistic, context) {
  function(p) {
    do.call(fn, c(list(params = p), slim,
                  list(per_arm = TRUE, context = context), fixed))[[statistic]]
  }
}

#' Shared body of the stage-two and bound bootstrap functions
#'
#' [slope_sample_size_boot()], [slope_power_boot()] and the two bound
#' bootstraps differ only in which calculation is re-solved on each replicate
#' (`fn`: a stage-two solver, or [floor_result()] for the bounds) and which of
#' its inputs is held fixed while doing so -- `fixed_name`, `"power"` or `"n"`;
#' everything after the observed result `x` has been solved is identical. A
#' bound has no design, so `x` carries no `design` field to pass on.
#'
#' `effectiveness` is carried as `x` reports it -- `NA` under
#' `target = "observed"`, where [target_components()] ignores it.
#' @noRd
bootstrap_stage_two <- function(x, fn, fixed_name, statistic, R, ci_method, level,
                                seed, progress, per_arm, context, cls) {
  fixed <- stats::setNames(list(x[[fixed_name]]), fixed_name)
  # `compute` is built by a factory rather than inline, so that its enclosing
  # environment holds only the arguments passed to it. Written inline, its
  # environment would be *this* frame, which holds `x` -- and so x$params$fit,
  # the original fit and its model frame -- keeping all of it reachable through
  # the closure for as long as the closure lives, however slim the `slim` copy
  # is. R closures capture environments, not the variables named in them.
  slim <- x[intersect(c("design", "target", "alpha", "effectiveness"), names(x))]
  compute <- boot_stage_two_compute(fn, slim, fixed, statistic, context)
  run_bootstrap(x$params, compute, x[[statistic]], statistic, R, ci_method, level,
                seed, progress, per_arm, context, cls, target = x$target)
}

#' Read every statistic off one refit, a failing statistic becoming `NA`
#'
#' Shared by the replicate and jackknife loops, which must agree on how a
#' statistic that cannot be computed off an otherwise good fit is recorded.
#' @noRd
eval_computes <- function(computes, p) {
  vapply(computes, function(f) tryCatch(f(p), error = function(e) NA_real_), numeric(1L))
}

#' Leave-one-subject-out fits, read for several quantities at once
#'
#' The refit is the whole cost of a jackknife; reading a second number off one
#' is free. Taking a list of accessors rather than a single `compute` is what
#' lets the sample size and the replicate slopes have BCa intervals apiece for
#' the price of one pass, which is what the printed table needs.
#'
#' A subject whose refit fails leaves a row of `NA`, so each column can be
#' cleaned independently -- one quantity failing to be computable off a fit does
#' not cost the others that subject.
#' @noRd
jackknife_values <- function(frame, subject_index, refitter, computes) {
  jack <- matrix(NA_real_, nrow = length(subject_index), ncol = length(computes))
  suppressWarnings(for (i in seq_along(subject_index)) {
    p <- tryCatch(refitter(frame[-subject_index[[i]], , drop = FALSE]),
                  error = function(e) NULL)
    if (!is.null(p)) {
      jack[i, ] <- eval_computes(computes, p)
    }
  })
  jack
}

#' Bias-corrected and accelerated interval, given a jackknife already taken
#'
#' The acceleration comes from a leave-one-subject-out jackknife, matching the
#' clustering used for the bootstrap itself; `jack` is that jackknife's values
#' for the one quantity this interval is for, taken by [jackknife_values()].
#'
#' Returns the interval, or a single string saying why it could not be built --
#' the two reasons being different enough that the caller's warning has to tell
#' them apart.
#'
#' The bias correction counts a replicate tied with the observed value as half
#' below rather than not below at all:
#' \eqn{z_0 = \Phi^{-1}((\#\{\theta^* < \hat\theta\} +
#' \#\{\theta^* = \hat\theta\}/2) / B)}. The textbook form counts only
#' those strictly below, which assumes a continuous statistic, where exact ties
#' have probability zero and the two agree. `n` is not continuous -- it lands on
#' the even integers -- and ties are ordinary: at an observed 96, 11 of 200
#' replicates came back exactly 96. Counting every one of them as "below" pushes
#' `z0` down and carries the whole interval with it, by more the smaller the
#' sample size, since smaller sizes collide on the lattice more often.
#'
#' Two things fall out of the half-count. A bootstrap distribution symmetric
#' about the observed value gives `z0 = 0` whether or not it has an atom there,
#' which is what "no bias to correct" ought to mean. And the degenerate branch
#' is no longer reachable by a pile of ties: with every replicate equal to the
#' observed value the proportion is 1/2, not 0, so BCa returns an interval
#' instead of falling back and blaming a one-sided distribution that was in fact
#' centred. What remains in that branch really is one-sided, which is what its
#' message now says.
#' @noRd
bca_from_jack <- function(theta, observed, jack, probs) {
  prop <- (sum(theta < observed) + sum(theta == observed) / 2) / length(theta)
  if (prop <= 0 || prop >= 1) {
    return(paste("the bias correction could not be computed (every replicate",
                 "falls strictly on one side of the observed value)"))
  }
  z0 <- stats::qnorm(prop)

  jack <- jack[!is.na(jack)]
  if (length(jack) < 3L) {
    return(paste("the acceleration could not be computed (under three",
                 "leave-one-subject-out fits succeeded)"))
  }

  d <- mean(jack) - jack
  denom <- 6 * (sum(d^2))^1.5
  a <- if (denom == 0) 0 else sum(d^3) / denom

  z <- stats::qnorm(probs)
  adj <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  if (any(!is.finite(adj))) {
    return("the bias-corrected endpoints were not finite")
  }
  stats::quantile(theta, adj, names = FALSE, type = 7)
}

#' Format a bootstrapped figure the way the printed columns are
#'
#' `format()` rather than [fmt_num()]: a per-arm figure and its total, or a
#' set of interval endpoints, must be padded to the *same* number of decimal
#' places so the column lines up digit-for-digit, which is what `format()`'s
#' vector-wide behaviour gives and what formatting each element on its own
#' takes away. `drop0trailing` keeps [fmt_num()]'s own convention, so a whole
#' number -- a sample size, say -- prints as "753" and not "753.0".
#' @noRd
boot_fmt <- function(v) format(v, digits = 6, trim = TRUE, drop0trailing = TRUE)

#' An interval as one character column, digits stacked over digits
#'
#' Both endpoints of every row are formatted in one [boot_fmt()] call, so the
#' whole column shares a decimal count, and then each side is padded against
#' the other rows before the "to" -- which is what makes "269 to  562" and
#' "538 to 1124" line up as a column rather than as two ragged numbers either
#' side of a fixed word.
#'
#' The result is one character column of a printed data frame, and
#' `print.data.frame()` right-justifies character columns, so the padding done
#' here is what keeps the interval aligned within its cells while R keeps the
#' cells themselves aligned.
#' @noRd
boot_interval_col <- function(lower, upper) {
  ends <- boot_fmt(c(lower, upper))
  k <- length(lower)
  los <- ends[seq_len(k)]
  his <- ends[k + seq_len(k)]
  paste(formatC(los, width = max(nchar(los))), "to",
        formatC(his, width = max(nchar(his))))
}

#' A labelled note, wrapped with a hanging indent
#'
#' The straddle note ran to 105 characters on one line and wrapped raggedly in
#' an 80-column terminal, which is where these are read. Wrapping here keeps the
#' label column of the notes lined up with each other and with the block above
#' them, whichever layout printed it.
#' @noRd
boot_note <- function(label, text, width = 78L) {
  lead <- formatC(paste0("  ", label, ":"), width = -15L)
  strwrap(text, width = width, initial = lead, prefix = strrep(" ", 15L))
}

#' The method, replicate count and interval level, as the first note
#'
#' These three used to ride in the span header of the table both print methods
#' drew. A data frame printed by R has no span header to put them in, and they
#' are not per-row facts, so they belong with the notes rather than in a column
#' repeating itself down the page. First of the notes because it says what the
#' three bootstrap columns above it are.
#' @noRd
boot_method_note <- function(ci_method, R, level) {
  boot_note("Bootstrap", sprintf("R = %d replicates; %.0f%% %s intervals.",
                                 R, 100 * level,
                                 if (identical(ci_method, "bca")) "BCa" else "percentile"))
}

#' Which basis every count on the page is on
#'
#' Printed first among the notes of [print.slope_bootstrap()] and
#' [print.slope_sample_size_grid_boot()], ahead of [boot_method_note()],
#' because it describes the whole frame -- every count column -- rather than
#' only the bootstrap columns that note explains. One sentence shared by both
#' print methods so a copy-edit to one cannot leave them describing the same
#' choice in different words, the reason [boot_method_note()] itself was
#' extracted.
#' @noRd
basis_note <- function(per_arm) {
  boot_note("Counts", sprintf(
    "`n` and the columns beside it are %s.",
    if (per_arm) "per arm" else "totals for the trial"))
}

#' How many replicates refit a slope on the wrong side of zero
#'
#' Printed unconditionally by both bootstrap print methods, including the 0/R
#' case. This is the measured form of the section 2.6 hazard, and a reader
#' checking whether the interval means anything needs to see that it was
#' checked -- an absent line is indistinguishable from a version of the package
#' that never looked. Suppressing it at zero also made the one number worth
#' reporting visible only when it was already bad news.
#'
#' The count is shown beside the percentage because the percentage alone hides
#' its own precision: "2.5%" off 40 replicates is one replicate, and a reader
#' deciding whether to worry needs to know it was one. Recovered by rounding
#' rather than stored, since `straddle` is a count over the replicates kept and
#' so the product is that count to within floating point.
#'
#' One sentence rather than two identical ones, because the two printers report
#' the same statistic and a copy-edit to one would otherwise leave them
#' describing it in different voices -- the drift [boot_method_note()] was
#' extracted to prevent, one note further on.
#' @noRd
boot_straddle_note <- function(straddle, n_used, of) {
  boot_note("Note", sprintf(paste0("%d/%d (%.1f%%) of replicates refit a %s on ",
                                   "the opposite side of zero from the fitted one."),
                            round(straddle * n_used), n_used, 100 * straddle,
                            of %||% "slope"))
}

#' The printed frame: one row for the statistic, on the chosen basis
#'
#' A sample size is the one bootstrapped statistic that comes in two units --
#' participants in total and participants per arm -- and `per_arm` picks
#' which one this frame shows, which is what `lattice` already records
#' (`run_bootstrap()` decides it once, so the print side never re-tests it
#' against `statistic`) together with the `per_arm` a caller requested or
#' [print.slope_bootstrap()] was asked to show instead.
#'
#' The per-arm entry is exactly half the total, in the arithmetic and not only
#' in the display. `solve_slope()` builds `n` as `2 * n_per_arm` and
#' `widen_to_lattice()` moves both interval endpoints out to even sizes, so
#' the half here is never rounded; the mean halves because it is linear.
#' @noRd
boot_summary_frame <- function(x, per_arm) {
  divisor <- boot_divisor(x, per_arm)
  data.frame(
    statistic  = x$statistic,
    calculated = x$observed / divisor,
    mean       = x$boot_mean / divisor,
    sd         = x$boot_sd / divisor,
    ci         = boot_interval_col(x$ci[1L] / divisor, x$ci[2L] / divisor),
    stringsAsFactors = FALSE
  )
}

#' What a bootstrap result's statistic is divided by on the display basis
#'
#' The one statement of the rule [boot_summary_frame()] and
#' [confint.slope_bootstrap()] both apply: a sample size, on the even lattice,
#' is halved to show it per arm; nothing else has arms. `per_arm` arrives
#' resolved by [display_basis()]. `isTRUE()`, so an object saved before
#' `lattice` existed is shown as stored rather than failing.
#' @noRd
boot_divisor <- function(x, per_arm) if (isTRUE(x$lattice) && per_arm) 2 else 1

#' Print a bootstrap result
#'
#' Printed as a data frame, by R itself: one row for the statistic, on the
#' chosen basis when it is a sample size, and columns for the calculated
#' value, the bootstrap mean and SD, and the interval. The method, replicate
#' count, interval level, and the failure and sign-straddling counts are
#' reported in the notes beneath, which print on every call whether or not
#' anything went wrong; for a sample size, a "Counts:" note ahead of them
#' says which basis -- per arm or total -- the row is on.
#'
#' The resampled slope is not shown. It is still on the object, in
#' `slope_observed` and the five fields beside it (see \sQuote{Value} in
#' [slope_sample_size_boot()]), and the straddle note still reports the one
#' thing about it that bears on whether the interval means anything.
#
# The result is shown as a data frame, printed by R itself, rather than as a
# hand-drawn table: every other tabular result in the package --
# slope_sample_size_grid() and slope_power_grid() -- is a plain data frame, and
# a reader who can read one of those can read this without learning a second
# layout. What the hand-drawn table had that a data frame does not -- a span
# header naming the method, the replicate count and the interval level -- moves
# into the notes below it, via boot_method_note().
#' @param x A result of [slope_sample_size_boot()], [slope_power_boot()],
#'   [slope_params_boot()], [slope_sample_size_floor_boot()] or
#'   [slope_power_ceiling_boot()].
#' @param ... Ignored.
#' @param per_arm Which basis to print a bootstrapped `n` on: `TRUE` for
#'   participants per arm, `FALSE` for the trial total. Defaults to `NULL`,
#'   meaning "whatever the bootstrap was called with" -- read from `x`'s
#'   `per_arm` attribute, or per arm if that is absent. Ignored, silently, for
#'   any statistic other than a sample size.
#' @return `x`, invisibly.
#' @export
print.slope_bootstrap <- function(x, ..., per_arm = NULL) {
  per_arm <- display_basis(x, per_arm, "print.slope_bootstrap()")

  cat(sprintf("<%s>\n\n", class(x)[1L]))
  print.data.frame(boot_summary_frame(x, per_arm))
  cat("\n")

  if (x$lattice) cat(basis_note(per_arm), sep = "\n")
  cat(boot_method_note(x$ci_method, x$R, x$level), sep = "\n")

  # Printed unconditionally, including the 0/R case, for the same reason as the
  # straddle note below: a clean run and a version of the package that never
  # checked must not look identical on the page. `n_failed` is already counted
  # out of `R` -- not out of the replicates that succeeded -- so the denominator
  # here is the number requested, matching what a reader asked for and what the
  # warning above (when there were any failures) already reported.
  # An object from before `n_refit_failed` existed counted every failure as a
  # convergence failure, and is printed as it always was.
  n_refit <- x$n_refit_failed %||% x$n_failed
  cat(boot_note("Note", sprintf(
    "%d/%d (%.1f%%) bootstrap samples failed to converge.",
    n_refit, x$R, 100 * n_refit / x$R)),
      sep = "\n")
  if (x$n_failed > n_refit) {
    cat(boot_note("Note", sprintf("%d/%d (%.1f%%) more %s.",
      x$n_failed - n_refit, x$R, 100 * (x$n_failed - n_refit) / x$R, BOOT_SOLVE_FAILED)),
      sep = "\n")
  }

  cat(boot_straddle_note(x$straddle, length(x$replicates), x$straddle_of), sep = "\n")

  # Where the rounding happens, and what follows from it. "mean 357.855" invites
  # the reading that 357.855 is the size the bootstrap recommends and that a
  # reader may round it back down to 356; it is an average of several hundred
  # trials each already rounded up, and rounding it again would be a third
  # rounding of one number. Two lines because this prints on every call --- the
  # full account, including which replicates are summarised and why the slope is
  # exempt, is in the `@return` section of ?slope_sample_size_boot.
  if (x$lattice) {
    cat(boot_note("Mean, SD", paste0(
      "each replicate is rounded up to a whole participant per arm before ",
      "averaging, so the mean is not a runnable ",
      if (per_arm) "arm" else "trial", " size.")),
        sep = "\n")
  }

  invisible(x)
}
