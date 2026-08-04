# Layer 2 --- the proposed trial design: visit schedule and dropout pattern.

#' Format numbers compactly for messages and printing
#' @noRd
fmt_num <- function(x) {
  vapply(x, function(v) format(v, trim = TRUE, drop0trailing = TRUE), character(1L))
}

#' Render a numeric vector as the R call that would create it
#' @noRd
fmt_call_vec <- function(x) {
  if (length(x) == 1L) fmt_num(x) else paste0("c(", paste(fmt_num(x), collapse = ", "), ")")
}

#' Describe a proposed trial's visit schedule and dropout pattern
#'
#' `trial_design()` bundles the design of a *future* trial --- when participants
#' are seen, and what proportion are expected to withdraw at each point --- into
#' a validated object that can be reused across many stage-two calls. It
#' carries no information about the outcome or its variability; that lives in the
#' [slope_params] object.
#'
#' @param visits Numeric vector of visit times, **including the baseline visit at
#'   time 0**, strictly increasing, of length at least 2. Times may be any real
#'   values --- `c(0, 0.5, 1, 1.5, 2)` is valid --- expressed in the same units as
#'   the `time` variable used to estimate the slope parameters.
#' @param dropout Optional numeric vector of expected dropout proportions, of
#'   length `length(visits) - 1`. `NULL` (the default) means no dropout.
#' @param dropout_type How `dropout` is expressed. `"incremental"` (the default)
#'   means element `j` is the proportion of participants whose **last attended
#'   visit is `visits[j]`**. `"cumulative"` means element `j` is the proportion
#'   who have withdrawn by `visits[j + 1]`, i.e. who fail to attend it. Whichever
#'   is supplied, the stored `dropout` field is always incremental.
#'
#' @details
#' # Baseline is explicit
#'
#' Stata's `schedule()` option lists follow-up visits only and assumes an implicit
#' baseline at time 0. This port requires baseline to be given explicitly, because
#' the implicit convention is a common source of off-by-one design errors. A
#' `visits` vector that does not start at 0 is rejected with a suggested correction
#' rather than silently repaired.
#'
#' # Real-valued visit times
#'
#' Stata restricts `schedule()` to ascending integers of at least 1 and provides a
#' `scale()` option to compensate, because it builds the covariance matrix on a
#' unit-integer grid. This port builds it directly at the requested times, so
#' `scale()` is unnecessary and has no equivalent: express `visits` in whatever
#' units the fitted slope uses.
#'
#' # Incremental versus cumulative dropout
#'
#' The two conventions describe the same strata from opposite ends, and confusing
#' them changes the trial design materially. A participant counted in
#' `dropout[j]` (incremental) attends `visits[1:j]` and misses everything after,
#' so their last attended visit is `visits[j]` and the first they miss is
#' `visits[j + 1]`. The print method shows both columns so the reading is never
#' ambiguous.
#'
#' Participants whose last attended visit is baseline contribute no follow-up
#' measurement and therefore no information about the slope. They are excluded
#' from the calculation entirely, and `trial_design()` warns when `dropout[1]` is
#' non-zero so that this is visible rather than silent.
#'
#' @return An object of class `trial_design`: a list with elements `visits`,
#'   `dropout` (always incremental), `has_dropout` and `dropout_type`.
#'
#' @examples
#' trial_design(c(0, 1, 2))
#' trial_design(c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
#' trial_design(seq(0, 3, by = 0.5), dropout = rep(0.025, 6))
#' trial_design(c(0, 1, 2, 3), dropout = c(0.05, 0.10, 0.15),
#'              dropout_type = "cumulative")
#'
#' @seealso [slope_params()], [slope_sample_size()], [slope_power()]
#' @export
trial_design <- function(visits,
                         dropout = NULL,
                         dropout_type = c("incremental", "cumulative")) {
  ctx <- "trial_design()"
  dropout_type <- match.arg(dropout_type)

  visits <- validate_visits(visits, ctx)
  n_visits <- length(visits)
  n_intervals <- n_visits - 1L

  dropout <- validate_dropout(dropout, n_intervals, dropout_type, visits, ctx)

  if (dropout[1L] > 0) {
    warning(sprintf(
      paste0("%s: dropout[1] = %s applies to participants whose last attended visit is ",
             "baseline (t = %s). With no follow-up measurement they contribute nothing to ",
             "the comparison of slopes and are excluded from the calculation. The Stata ",
             "original skips this stratum silently."),
      ctx, fmt_num(dropout[1L]), fmt_num(visits[1L])), call. = FALSE)
  }

  structure(
    list(
      visits       = visits,
      dropout      = dropout,
      has_dropout  = any(dropout > 0),
      dropout_type = dropout_type
    ),
    class = "trial_design"
  )
}

#' Validate the visit schedule
#' @noRd
validate_visits <- function(visits, ctx) {
  if (!is.numeric(visits) || length(visits) < 2L) {
    stop(sprintf(paste0("%s: `visits` must be a numeric vector of at least 2 visit times ",
                        "(baseline plus at least one follow-up); got %s of length %d."),
                 ctx, class(visits)[1L], length(visits)), call. = FALSE)
  }
  if (any(!is.finite(visits))) {
    stop(sprintf("%s: `visits` must be finite; element(s) %s are not.",
                 ctx, paste(which(!is.finite(visits)), collapse = ", ")), call. = FALSE)
  }
  visits <- as.numeric(visits)

  dups <- unique(visits[duplicated(visits)])
  if (length(dups) > 0L) {
    stop(sprintf(paste0("%s: `visits` must not contain repeated times; %s appear%s more than ",
                        "once. Repeated visit times make the covariance matrix singular."),
                 ctx, paste(fmt_num(dups), collapse = ", "),
                 if (length(dups) == 1L) "s" else ""), call. = FALSE)
  }
  if (is.unsorted(visits)) {
    stop(sprintf("%s: `visits` must be in increasing order; got %s.",
                 ctx, fmt_call_vec(visits)), call. = FALSE)
  }

  if (visits[1L] != 0) {
    msg <- sprintf(paste0("%s: `visits` must begin with the baseline visit at time 0; ",
                          "got visits[1] = %s."),
                   ctx, fmt_num(visits[1L]))
    if (visits[1L] > 0) {
      msg <- paste0(
        msg,
        sprintf(paste0("\n  Stata's schedule() lists follow-up visits only and assumes an ",
                       "implicit baseline at time 0; this port requires it explicitly.",
                       "\n  Did you mean: visits = %s"),
                fmt_call_vec(c(0, visits))))
    } else {
      msg <- paste0(msg, "\n  Visit times are measured from baseline, so none may be negative.")
    }
    stop(msg, call. = FALSE)
  }

  visits
}

#' Validate and normalise the dropout vector to incremental proportions
#' @noRd
validate_dropout <- function(dropout, n_intervals, dropout_type, visits, ctx) {
  if (is.null(dropout)) {
    return(rep(0, n_intervals))
  }

  if (!is.numeric(dropout)) {
    stop(sprintf("%s: `dropout` must be numeric or NULL; got %s.",
                 ctx, class(dropout)[1L]), call. = FALSE)
  }
  if (any(!is.finite(dropout))) {
    stop(sprintf("%s: `dropout` must be finite; element(s) %s are not.",
                 ctx, paste(which(!is.finite(dropout)), collapse = ", ")), call. = FALSE)
  }
  dropout <- as.numeric(dropout)

  if (length(dropout) != n_intervals) {
    stop(sprintf(paste0("%s: `dropout` must have one element per follow-up visit: ",
                        "length(visits) - 1 = %d, but length(dropout) = %d.",
                        "\n  visits = %s covers %d follow-up visit%s after baseline."),
                 ctx, n_intervals, length(dropout),
                 fmt_call_vec(visits), n_intervals,
                 if (n_intervals == 1L) "" else "s"), call. = FALSE)
  }

  if (any(dropout < 0)) {
    stop(sprintf("%s: `dropout` proportions must be non-negative; element(s) %s are not.",
                 ctx, paste(which(dropout < 0), collapse = ", ")), call. = FALSE)
  }

  if (identical(dropout_type, "cumulative")) {
    steps <- diff(dropout)
    bad <- which(steps < -DROPOUT_TOL)
    if (length(bad) > 0L) {
      j <- bad[1L] + 1L
      stop(sprintf(paste0("%s: cumulative `dropout` must be non-decreasing; element %d (%s) ",
                          "is smaller than element %d (%s). Once a participant has withdrawn ",
                          "they cannot return."),
                   ctx, j, fmt_num(dropout[j]), j - 1L, fmt_num(dropout[j - 1L])),
           call. = FALSE)
    }
    if (dropout[n_intervals] > 1 + DROPOUT_TOL) {
      stop(sprintf(paste0("%s: cumulative `dropout` cannot exceed 1; the final element is %s. ",
                          "Proportions, not percentages, are expected."),
                   ctx, fmt_num(dropout[n_intervals])), call. = FALSE)
    }
    dropout <- pmax(diff(c(0, dropout)), 0)
  } else {
    total <- sum(dropout)
    if (total > 1 + DROPOUT_TOL) {
      stop(sprintf(paste0("%s: incremental `dropout` proportions sum to %s, which exceeds 1. ",
                          "Each element is the proportion whose last attended visit is that ",
                          "visit, so they partition the randomised sample and cannot total ",
                          "more than 1.\n  If you meant cumulative proportions, pass ",
                          "dropout_type = \"cumulative\"."),
                   ctx, fmt_num(total)), call. = FALSE)
    }
  }

  dropout
}

#' @describeIn trial_design Print a trial design.
#' @param x A `trial_design` object.
#' @param ... Ignored.
#' @export
print.trial_design <- function(x, ...) {
  n_visits <- length(x$visits)
  k <- n_visits - 1L

  cat("<trial_design>\n")
  cat(sprintf("  Visits (%d):  %s\n", n_visits, paste(fmt_num(x$visits), collapse = ", ")))
  cat(sprintf("  Follow-up:   %d visit%s, last at t = %s\n",
              k, if (k == 1L) "" else "s", fmt_num(x$visits[n_visits])))

  if (!x$has_dropout) {
    cat("  Dropout:     none; all participants attend every visit\n")
    return(invisible(x))
  }

  cat(sprintf("  Dropout:     supplied as %s\n\n", x$dropout_type))
  cat("    last visit   first missed   proportion   cumulative\n")
  cum <- cumsum(x$dropout)
  for (j in seq_len(k)) {
    cat(sprintf("    %10s   %12s   %10s   %10s\n",
                fmt_num(x$visits[j]),
                fmt_num(x$visits[j + 1L]),
                formatC(x$dropout[j], format = "f", digits = 3),
                formatC(cum[j], format = "f", digits = 3)))
  }
  cat(sprintf("\n  Completers:  %s attend all %d visits\n",
              formatC(1 - sum(x$dropout), format = "f", digits = 3), n_visits))

  invisible(x)
}
