# Layer 2 --- the proposed trial design: visit schedule and dropout pattern.

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
#' # How dropout enters the calculation
#'
#' Dropout is accommodated by the pattern-mixture approach of Dawson and Lagakos
#' (1991, 1993), which is what section 2.5 of Nash et al. (2021) adopts and what
#' the Stata command's `dropouts()` option does. Participants are divided into
#' strata by the visits they attend: stratum `j` attends `visits[1:j]` and
#' nothing after, and the completers --- a proportion `1 - sum(dropout)` ---
#' attend everything. Each stratum is sized as though the entire trial followed
#' that one visit pattern, and the strata are combined as the reciprocal of the
#' weighted mean of the reciprocals of those stratum-specific sample sizes, the
#' weights being `dropout`. Equivalently, and this is how it is actually
#' computed, the squared standardised effect sizes are averaged with those
#' weights; see [slope_effect_size()].
#'
#' Withdrawers therefore still contribute the visits they did attend. That makes
#' the adjustment less conservative than inflating a completers-only sample size
#' by `1 / (1 - sum(dropout))`, and it is the appropriate one when the trial will
#' be analysed with a mixed model fitted to all observed measurements, as in
#' [slope_params()]. If the planned analysis instead discards partial records,
#' this will understate the sample size needed.
#'
#' Two assumptions come with it. Dropout is monotone and truncates a schedule
#' common to every participant: withdrawal is permanent, and intermittent
#' missingness --- a visit missed and follow-up resumed --- has no representation
#' here, as in the Stata original. And every stratum is given the same slope
#' difference and the same variance components, so withdrawal is assumed
#' unrelated to a participant's own trajectory; dropout driven by how fast
#' someone is declining is outside the model.
#'
#' @return An object of class `trial_design`: a list with elements `visits`,
#'   `dropout` (always incremental), `has_dropout` and `dropout_type`.
#'
#' @references
#' Dawson, J. D., and S. W. Lagakos. 1991. Analyzing laboratory marker changes in
#' AIDS clinical trials. \emph{Journal of Acquired Immune Deficiency Syndromes}
#' 4: 667--676.
#'
#' Dawson, J. D., and S. W. Lagakos. 1993. Size and power of two-sample tests of
#' repeated measures data. \emph{Biometrics} 49: 1022--1032.
#' \doi{10.2307/2532244}
#'
#' Frost, C., M. G. Kenward, and N. C. Fox. 2008. Optimizing the design of
#' clinical trials where the outcome is a rate. Can estimating a baseline rate in
#' a run-in period increase efficiency? \emph{Statistics in Medicine} 27:
#' 3717--3731. \doi{10.1002/sim.3280}
#'
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
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
  warn_baseline_dropout(dropout, visits, ctx)

  out <- structure(
    list(
      visits       = visits,
      dropout      = dropout,
      has_dropout  = any(dropout > 0),
      dropout_type = dropout_type
    ),
    class = "trial_design"
  )
  # Records exactly what was checked, so as_trial_design() can tell a design
  # that reaches it unchanged from the constructor -- already warned about,
  # right here -- from one that was hand-built or had `$dropout` edited
  # afterwards, which was never checked at all. See as_trial_design().
  attr(out, "slopepower_checked_dropout") <- dropout
  out
}

#' Warn when the first dropout stratum attends only the baseline visit
#'
#' Split out of [trial_design()] so [as_trial_design()] can apply exactly the
#' same check, in the same words, to a hand-built or subsequently edited
#' `trial_design` object -- one that never went through this constructor and
#' so was never warned about at all.
#' @noRd
warn_baseline_dropout <- function(dropout, visits, ctx) {
  if (dropout[1L] > 0) {
    # Classed rather than left as a plain warning: slope_power_grid() and
    # slope_sample_size_grid() collect this one specifically, by class, to
    # report it once per grid rather than once per cell.
    warning(warningCondition(
      sprintf(
        paste0("%s: dropout[1] = %s applies to participants whose last attended visit is ",
               "baseline (t = %s). With no follow-up measurement they contribute nothing to ",
               "the comparison of slopes and are excluded from the calculation. The Stata ",
               "original skips this stratum silently."),
        ctx, fmt_num(dropout[1L]), fmt_num(visits[1L])),
      class = "slopepower_baseline_dropout"))
  }
  invisible(dropout)
}

#' Validate the visit schedule
#' @noRd
validate_visits <- function(visits, ctx) {
  if (!is.numeric(visits) || length(visits) < 2L) {
    stop(sprintf(paste0("%s: `visits` must be a numeric vector of at least 2 visit times ",
                        "(baseline plus at least one follow-up); got %s of length %d."),
                 ctx, class(visits)[1L], length(visits)), call. = FALSE)
  }
  check_finite_vector(visits, "visits", ctx)
  visits <- as.numeric(visits)

  dups <- unique(visits[duplicated(visits)])
  if (length(dups) > 0L) {
    stop(sprintf(paste0("%s: `visits` must not contain repeated times; %s appear%s more than ",
                        "once. Repeated visit times make the covariance matrix singular."),
                 ctx, label_numeric(dups),
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
  dropout <- as.numeric(dropout)

  check_dropout_length(dropout, visits, "dropout", ctx)

  check_dropout_values(dropout, "dropout", ctx)

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
    check_dropout_total(dropout, "dropout", ctx)
  }

  dropout
}

#' The value rules an incremental dropout vector obeys, whatever built it
#'
#' Split out so that a hand-assembled `trial_design` reaching stage two through
#' [as_trial_design()] is held to the same rules, *and told about a breach in
#' the same words*, as one that came from the constructor. Wording that differs
#' by route means the same mistake gets a better diagnosis on one path than the
#' other, which is precisely backwards: the hand-built object is the one whose
#' author had no constructor to guide them.
#'
#' Split in two because `validate_dropout()` applies the total only to a vector
#' the user supplied as incremental -- a cumulative one is bounded by its own
#' final element instead, before conversion.
#' @noRd
check_dropout_values <- function(dropout, name, ctx) {
  check_finite_vector(dropout, name, ctx)
  if (any(dropout < 0)) {
    stop(sprintf("%s: `%s` proportions must be non-negative; element(s) %s are not.",
                 ctx, name, paste(which(dropout < 0), collapse = ", ")), call. = FALSE)
  }
  invisible(dropout)
}

#' @rdname check_dropout_values
#'
#' The length rule, shared for the same reason and by the same argument:
#' [as_trial_design()] used to write its own terser version of this message, so
#' a hand-built design -- the one route whose author had no constructor to
#' guide them -- got the *worse* of the two diagnoses for the identical
#' mistake. `name` carries the caller's spelling of the vector so the message
#' names what the user actually typed.
#' @noRd
check_dropout_length <- function(dropout, visits, name, ctx) {
  n_intervals <- length(visits) - 1L
  if (length(dropout) != n_intervals) {
    stop(sprintf(paste0("%s: `%s` must have one element per follow-up visit: ",
                        "length(visits) - 1 = %d, but length(%s) = %d.",
                        "\n  visits = %s covers %d follow-up visit%s after baseline."),
                 ctx, name, n_intervals, name, length(dropout),
                 fmt_call_vec(visits), n_intervals,
                 if (n_intervals == 1L) "" else "s"), call. = FALSE)
  }
  invisible(dropout)
}

#' @rdname check_dropout_values
#' @noRd
check_dropout_total <- function(dropout, name, ctx) {
  total <- sum(dropout)
  if (total > 1 + DROPOUT_TOL) {
    stop(sprintf(paste0("%s: incremental `%s` proportions sum to %s, which exceeds 1. ",
                        "Each element is the proportion whose last attended visit is that ",
                        "visit, so they partition the randomised sample and cannot total ",
                        "more than 1.\n  If you meant cumulative proportions, pass ",
                        "dropout_type = \"cumulative\"."),
                 ctx, name, fmt_num(total)), call. = FALSE)
  }
  invisible(dropout)
}

#' Coerce and validate the `design` argument of a stage-two call
#'
#' Accepts a `trial_design` object, or a bare numeric vector of visit times
#' which is passed to [trial_design()]. Lives here, beside the class it
#' validates, rather than with the calculations that call it.
#'
#' The rules for `visits` and `dropout` come from the constructor's own
#' validators, so that a hand-built `trial_design` is held to exactly the same
#' standard as one from [trial_design()]. Normalisation is not repeated here:
#' `dropout` is already incremental by the time a `trial_design` object exists.
#'
#' The baseline-only warning *is* repeated, but only when needed: a design
#' just returned by [trial_design()], unmodified, was already warned about
#' there, and repeating it here would warn twice for that, the ordinary,
#' path. A hand-built object, or one whose `$dropout` was edited after
#' construction, carries no record of ever having been warned about, and gets
#' checked now instead -- see the `slopepower_checked_dropout` attribute set
#' by [trial_design()].
#' @noRd
as_trial_design <- function(design, context) {
  if (is.numeric(design)) design <- trial_design(visits = design)
  if (!inherits(design, "trial_design")) {
    stop(sprintf("%s: `design` must be a `trial_design` object or a numeric vector of visit times.",
                 context), call. = FALSE)
  }
  visits <- validate_visits(design$visits, context)
  dropout <- design$dropout
  if (!is.numeric(dropout)) {
    stop(sprintf("%s: `design$dropout` must be numeric; got %s.",
                 context, class(dropout)[1L]), call. = FALSE)
  }
  check_dropout_length(dropout, visits, "design$dropout", context)
  check_dropout_values(dropout, "design$dropout", context)
  check_dropout_total(dropout, "design$dropout", context)

  if (!identical(attr(design, "slopepower_checked_dropout"), dropout)) {
    warn_baseline_dropout(dropout, visits, context)
  }

  if (!identical(design$dropout_type, "incremental") &&
      !identical(design$dropout_type, "cumulative")) {
    stop(sprintf("%s: `design$dropout_type` must be \"incremental\" or \"cumulative\".",
                 context), call. = FALSE)
  }
  # `dropout_type` records only how the user supplied the vector; `dropout`
  # itself is always incremental, converted once by the constructor. So a value
  # of "cumulative" is provenance, not an instruction, and must not be acted on
  # here -- doing so rejected every design built with dropout_type =
  # "cumulative", including the one in trial_design()'s own examples, with an
  # error telling the user to do what they had already done. A hand-built list
  # that puts cumulative values in `dropout` cannot be detected anyway: the
  # values are indistinguishable from valid incremental ones.
  #
  # `has_dropout` is what the calculation branches on, so it is derived here
  # rather than trusted. A hand-built design that omits it would otherwise fail
  # with "argument is of length zero", and one that sets it FALSE alongside a
  # non-zero `dropout` would silently report the unweighted s*^2.
  design$has_dropout <- any(dropout > 0)
  design$visits <- visits
  invisible(design)
}

#' @describeIn trial_design Print a trial design.
#' @param x A `trial_design` object.
#' @param ... Ignored.
#' @export
print.trial_design <- function(x, ...) {
  n_visits <- length(x$visits)
  k <- n_visits - 1L

  cat("<trial_design>\n")
  cat(sprintf("  Visits (%d):  %s\n", n_visits, label_numeric(x$visits)))
  cat(sprintf("  Follow-up:   %d visit%s, last at t = %s\n",
              k, if (k == 1L) "" else "s", fmt_num(x$visits[n_visits])))

  if (!x$has_dropout) {
    cat("  Dropout:     none; all participants attend every visit\n")
    return(invisible(x))
  }

  cat(sprintf("  Dropout:     supplied as %s\n\n", x$dropout_type))
  cat("    last visit   first missed   proportion   cumulative\n")
  # `dropout` and `cum` are already length k; only the visit columns need
  # offsetting against each other.
  cum <- cumsum(x$dropout)
  cat(sprintf("    %10s   %12s   %10s   %10s\n",
              fmt_num(x$visits[-n_visits]), fmt_num(x$visits[-1L]),
              formatC(x$dropout, format = "f", digits = 3),
              formatC(cum, format = "f", digits = 3)),
      sep = "")
  cat(sprintf("\n  Completers:  %s attend all %d visits\n",
              formatC(1 - sum(x$dropout), format = "f", digits = 3), n_visits))

  invisible(x)
}
