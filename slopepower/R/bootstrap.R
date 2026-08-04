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

#' Refit stage one on a resampled frame
#'
#' The times in the recovered frame have already been re-origined per subject, so
#' `origin = "none"` avoids repeating that work and the message that goes with it.
#' @noRd
refit_frame <- function(frame, comparator) {
  args <- list(formula = y ~ time | subject, data = frame, origin = "none")
  cl <- as.call(c(list(quote(slope_params)), args,
                  if (identical(comparator, "healthy")) list(healthy = quote(group))
                  else if (identical(comparator, "treated")) list(treated = quote(group))
                  else NULL))
  suppressMessages(suppressWarnings(eval(cl, list(frame = frame), parent.frame())))
}

#' Build one cluster-resampled frame, stratified by group where present
#' @noRd
resample_frame <- function(frame, subject_index, groups) {
  picks <- unlist(lapply(groups, function(ids) sample(ids, length(ids), replace = TRUE)),
                  use.names = FALSE)
  parts <- vector("list", length(picks))
  for (i in seq_along(picks)) {
    rows <- subject_index[[picks[i]]]
    part <- frame[rows, , drop = FALSE]
    part$subject <- i          # fresh identifier: a subject drawn twice is two people
    parts[[i]] <- part
  }
  do.call(rbind, parts)
}

#' Standard error of the reported slope, from the fitted model
#' @noRd
slope_se <- function(params) {
  fit <- params$fit
  if (is.null(fit)) return(NA_real_)
  b <- tryCatch(nlme::fixef(fit), error = function(e) return(NULL))
  V <- tryCatch(stats::vcov(fit), error = function(e) return(NULL))
  if (is.null(b) || is.null(V)) return(NA_real_)
  terms <- switch(params$comparator,
                  none    = "sp_time",
                  healthy = c("sp_time", "sp_case:sp_time"),
                  treated = c("sp_time", "sp_placebo_time"))
  if (!all(terms %in% names(b))) return(NA_real_)
  k <- as.numeric(names(b) %in% terms)
  sqrt(drop(k %*% as.matrix(V) %*% k))
}

#' Bootstrap the stage-one estimates
#'
#' Resamples subjects with replacement, refits the stage-one mixed model on each
#' replicate, and recomputes a statistic of interest. This propagates the
#' estimation uncertainty in the slope and variance components through to the
#' sample size or power, as recommended in section 2.6 of Nash et al. (2021).
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
#' @param params A `slope_params` object produced by [slope_params()]. Objects
#'   from [slope_params_manual()] cannot be bootstrapped: they carry no data.
#' @param R Number of bootstrap replicates.
#' @param type `"percentile"` (the default) or `"bca"` for bias-corrected and
#'   accelerated intervals. The paper recommends BCa because the distribution of
#'   estimated sample sizes is typically skewed.
#' @param statistic Which quantity to bootstrap. `"n"`, `"power"` and `"tte"`
#'   require the design arguments of [slope_power()] to be supplied in `...`.
#' @param ... Passed to [slope_power()]; typically `design`, `effectiveness`, and
#'   one of `n` or `power`.
#' @param level Confidence level for the interval.
#' @param seed Optional integer seed, for reproducibility.
#' @param progress Report progress while resampling.
#'
#' @return An object of class `slope_bootstrap`, with elements `observed`, `replicates`,
#'   `ci`, `type`, `statistic`, `R`, `n_failed` and `se`.
#'
#' @references
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#'
#' @examples
#' \dontrun{
#' pars <- slope_params(sdmt ~ visit | id, data = d1)
#' slope_bootstrap(pars, R = 200, design = c(0, 1, 2), effectiveness = 0.33)
#' }
#'
#' @seealso [slope_power()], [slope_params()]
#' @export
slope_bootstrap <- function(params, R = 199,
                            type = c("percentile", "bca"),
                            statistic = c("n", "power", "tte", "slope"),
                            ..., level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_bootstrap()"
  type <- match.arg(type)
  statistic <- match.arg(statistic)
  check_scalar(R, "R", context, lower = 1, upper = Inf, lower_open = FALSE)
  check_probability(level, "level", context)
  if (!is.null(seed)) set.seed(seed)

  dots <- list(...)
  needs_design <- !identical(statistic, "slope")
  if (needs_design && is.null(dots$design)) {
    stop(sprintf('%s: statistic = "%s" needs a trial design; pass `design` (and any of ',
                 context, statistic),
         "`effectiveness`, `n`, `power`, `alpha`) through `...`.", call. = FALSE)
  }

  compute <- function(p) {
    if (identical(statistic, "slope")) return(p$slope)
    res <- do.call(slope_power, c(list(params = p), dots))
    switch(statistic, n = res$n, power = res$power, tte = res$tte)
  }

  observed <- compute(params)

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
  ids <- names(subject_index)
  groups <- if (is.null(frame$group)) {
    list(ids)
  } else {
    sub_group <- vapply(subject_index, function(r) frame$group[r[1L]], numeric(1L))
    unname(split(ids, sub_group))
  }

  replicates <- rep(NA_real_, R)
  tick <- max(1L, floor(R / 10))
  for (b in seq_len(R)) {
    replicates[b] <- tryCatch(compute(refit_frame(resample_frame(frame, subject_index, groups),
                                                  params$comparator)),
                              error = function(e) NA_real_)
    if (isTRUE(progress) && b %% tick == 0L) {
      message(sprintf("%s: %d of %d replicates", context, b, R))
    }
  }

  n_failed <- sum(is.na(replicates))
  good <- replicates[!is.na(replicates)]
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
    bca <- bca_interval(good, observed, frame, subject_index, params$comparator,
                        compute, probs, context)
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

#' Bias-corrected and accelerated interval
#'
#' The acceleration comes from a leave-one-subject-out jackknife, matching the
#' clustering used for the bootstrap itself.
#' @noRd
bca_interval <- function(theta, observed, frame, subject_index, comparator,
                         compute, probs, context) {
  prop <- mean(theta < observed)
  if (prop <= 0 || prop >= 1) return(NULL)
  z0 <- stats::qnorm(prop)

  ids <- names(subject_index)
  jack <- rep(NA_real_, length(ids))
  for (i in seq_along(ids)) {
    drop_rows <- subject_index[[i]]
    jack[i] <- tryCatch(compute(refit_frame(frame[-drop_rows, , drop = FALSE], comparator)),
                        error = function(e) NA_real_)
  }
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
#' @param x A `slope_bootstrap` object.
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
