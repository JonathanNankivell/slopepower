# Layer 2 --- the proposed trial design: visit schedule and dropout pattern.

#' Render a numeric vector as the R call that would create it
#' @noRd
fmt_call_vec <- function(x) {
  if (length(x) == 1L) fmt_num(x) else paste0("c(", paste(fmt_num(x), collapse = ", "), ")")
}

# ---------------------------------------------------------------------------
# dropout rates
#
# `dropout_rate()` lives beside `trial_design()` rather than with the grid
# functions that were once its only consumer, because it is an alternative way
# of saying what `trial_design()`'s `dropout` argument means, not a feature of
# the grid. It was originally defined in grid.R, and the consequence was that
# the object could only be expanded by the grid: `trial_design(v,
# dropout_rate(0.05))` -- the call its own documentation pointed at -- was
# rejected as a non-numeric dropout.
# ---------------------------------------------------------------------------

#' A constant dropout rate per unit of time
#'
#' Dropout proportions are supplied to [trial_design()] per visit, so the same
#' underlying withdrawal rate needs a different vector for every candidate visit
#' schedule: "5% per year" over three years is `0.15` for a trial with a single
#' final visit, `rep(0.05, 3)` with annual visits, and `rep(0.025, 6)` with
#' six-monthly visits. `dropout_rate()` expresses the rate once and lets
#' [trial_design()] expand it correctly for whichever schedule it is given ---
#' which is what makes it useful to the grid functions, where one object drives
#' every row of a table of competing schedules.
#'
#' `type` chooses which of two withdrawal patterns `rate` describes:
#'
#' * `"linear"` (the default) applies `rate` to the *original randomised
#'   cohort*, matching the worked example in section 4.2 of Nash et al. (2021):
#'   the proportion whose last attended visit is `visits[j]` is
#'   `rate * (visits[j + 1] - visits[j]) / per`. The same number withdraw in
#'   every equal-length interval regardless of how many have already left, and
#'   the increments sum to `rate * total_duration / per`.
#' * `"cumulative"` applies `rate` as the proportion of *whoever is still in
#'   follow-up* who withdraws per `per` units of time, so the same proportion
#'   of the remaining participants drops out at every visit rather than the
#'   same proportion of the original cohort. The fraction still being followed
#'   at time `t` is `(1 - rate) ^ (t / per)`, decaying geometrically, and the
#'   proportion whose last attended visit is `visits[j]` is the drop in that
#'   fraction over the interval:
#'   `(1 - rate) ^ (visits[j] / per) - (1 - rate) ^ (visits[j + 1] / per)`.
#'   Because the survival fraction approaches but never passes zero, the total
#'   never exceeds 1 however long the trial runs, and `rate` --- itself a
#'   proportion of the remaining sample --- is restricted to `[0, 1]`.
#'
#' Because both expansions produce incremental proportions by construction, a
#' `dropout_rate` cannot be combined with `dropout_type = "cumulative"`; see
#' [trial_design()]. That argument is a different, unrelated choice --- how the
#' *vector itself* was written down --- from this object's `type`, which
#' describes how withdrawal behaves over time; the shared word is coincidence,
#' not the same idea twice.
#'
#' This object only produces the per-visit proportions. What the calculation then
#' does with them --- the Dawson and Lagakos (1991, 1993) pattern mixture, and
#' what it assumes about why people withdraw --- is described in
#' [trial_design()].
#'
#' @param rate Expected proportion withdrawing per `per` units of time. For
#'   `type = "linear"`, a proportion of the original randomised sample and so
#'   non-negative with no upper bound of its own (the total across all visits
#'   is what is capped at 1). For `type = "cumulative"`, a proportion of
#'   whoever remains and so restricted to `[0, 1]`.
#' @param per Length of time `rate` refers to, in the units of the `time` variable
#'   used to estimate the slope parameters. Defaults to 1, i.e. `rate` is a
#'   per-unit-time rate.
#' @param type Which withdrawal pattern `rate` describes: `"linear"` (the
#'   default) or `"cumulative"`. See Details.
#'
#' @return An object of class `dropout_rate`.
#'
#' @examples
#' dropout_rate(0.05)              # 5% of the original cohort per unit time
#' dropout_rate(0.10, per = 12)    # 10% per 12 months, if time is in months
#' dropout_rate(0.05, type = "cumulative")  # 5% of those remaining, per unit time
#'
#' # The same rate, expanded for two different schedules
#' trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.05))$dropout
#' trial_design(c(0, 1.5, 3), dropout = dropout_rate(0.05))$dropout
#'
#' @seealso [trial_design()], [slope_power_grid()], [slope_sample_size_grid()]
#' @export
dropout_rate <- function(rate, per = 1, type = c("linear", "cumulative")) {
  context <- "dropout_rate()"
  type <- match.arg(type)
  if (identical(type, "cumulative")) {
    # A cumulative rate is a proportion of the remaining sample and is the base
    # of a power in expand_dropout_rate(); anything above 1 would describe more
    # than everyone remaining leaving, and a non-integer exponent of a negative
    # base is complex, not a dropout proportion.
    check_scalar(rate, "rate", context, lower = 0, upper = 1, lower_open = FALSE,
                 upper_open = FALSE)
  } else {
    check_scalar(rate, "rate", context, lower = 0, upper = Inf, lower_open = FALSE)
  }
  check_scalar(per, "per", context, lower = 0, upper = Inf, lower_open = TRUE)
  structure(list(rate = as.numeric(rate), per = as.numeric(per), type = type),
            class = "dropout_rate")
}

#' @describeIn dropout_rate Print a dropout rate.
#' @param x A `dropout_rate` object.
#' @param ... Ignored.
#' @export
print.dropout_rate <- function(x, ...) {
  cat(sprintf("<dropout_rate> %s per %s unit%s of time (%s)\n",
              fmt_num(x$rate), fmt_num(x$per), if (x$per == 1) "" else "s", x$type))
  invisible(x)
}

#' Expand a `dropout_rate` into incremental proportions for one schedule
#'
#' The single place a rate becomes numbers, called by [validate_dropout()] on
#' the ordinary path and by `expand_dropout()` (grid.R) on the grid path, so the
#' two cannot expand the same object differently. `where` carries the grid's
#' cell label into the message; it is empty for a direct [trial_design()] call,
#' which has no cell to name.
#'
#' The two `type`s (see [dropout_rate()]) compute the incremental proportions
#' differently: `"linear"` applies `rate` to the original cohort, so the total
#' across the whole trial can exceed 1 and is checked for that here, with a
#' diagnosis in terms of `rate` and `per` rather than the expanded vector.
#' `"cumulative"` applies `rate` to whoever remains, compounding geometrically,
#' so the total is mathematically bounded below 1 and no equivalent check is
#' needed. Both callers run `increments` through `check_dropout_total()`
#' afterwards regardless, so a bug in this earlier, friendlier check could not
#' let an invalid total through uncaught.
#' @noRd
expand_dropout_rate <- function(spec, visits, ctx, where = "") {
  if (identical(spec$type, "cumulative")) {
    survival <- (1 - spec$rate) ^ (visits / spec$per)
    return(-diff(survival))
  }

  increments <- (spec$rate / spec$per) * diff(visits)
  total <- sum(increments)
  # The same bound `check_dropout_total()` enforces on `dropout` generally --
  # reusing its shared DROPOUT_TOL so the threshold itself cannot drift between
  # the two -- but diagnosed here in terms of `rate` and `per` rather than the
  # expanded vector, because the general check would otherwise report a
  # `dropout_rate()` mistake by naming a vector the caller never wrote.
  if (total > 1 + DROPOUT_TOL) {
    stop(sprintf(paste0("%s%s: a rate of %s per %s unit(s) of time over a trial lasting %s ",
                        "implies total dropout of %s, which exceeds 1."),
                 ctx, where, fmt_num(spec$rate), fmt_num(spec$per),
                 fmt_num(diff(range(visits))), fmt_num(total)), call. = FALSE)
  }
  increments
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
#' @param dropout Optional expected dropout. Either a numeric vector of
#'   proportions, of length `length(visits) - 1`, or a [dropout_rate()] object
#'   giving a constant rate per unit of time, which is expanded to one
#'   proportion per interval of `visits`. `NULL` (the default) means no dropout.
#' @param dropout_type How `dropout` is expressed. `"incremental"` (the default)
#'   means element `j` is the proportion of participants whose **last attended
#'   visit is `visits[j]`**. `"cumulative"` means element `j` is the proportion
#'   who have withdrawn by `visits[j + 1]`, i.e. who fail to attend it. Whichever
#'   is supplied, the stored `dropout` field is always incremental. A
#'   [dropout_rate()] produces incremental proportions by construction, so it
#'   cannot be combined with `"cumulative"`.
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
#' \doi{10.1177/1536867X211045512}
#'
#' @examples
#' trial_design(c(0, 1, 2))
#' trial_design(c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
#' trial_design(seq(0, 3, by = 0.5), dropout = rep(0.025, 6))
#' trial_design(c(0, 1, 2, 3), dropout = c(0.05, 0.10, 0.15),
#'              dropout_type = "cumulative")
#'
#' # The same 2.5% per unit time, stated once rather than per schedule
#' trial_design(seq(0, 3, by = 0.5), dropout = dropout_rate(0.025))
#'
#' @seealso [slope_params()], [slope_sample_size()], [slope_power()],
#'   [dropout_rate()]
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

  if (inherits(dropout, "dropout_rate")) {
    # Rejected rather than ignored. A rate expands to incremental proportions
    # by construction, so there is no sense in which it could have been
    # supplied cumulatively; accepting the pair would store incremental values
    # in a design whose `dropout_type` -- and whose print method -- announce
    # them as cumulative.
    if (identical(dropout_type, "cumulative")) {
      stop(sprintf(paste0("%s: a dropout_rate() cannot be combined with dropout_type = ",
                          "\"cumulative\"; it expands to the proportion withdrawing within ",
                          "each interval, which is already incremental. Drop the ",
                          "`dropout_type` argument, or supply the cumulative proportions ",
                          "as a numeric vector."),
                   ctx), call. = FALSE)
    }
    dropout <- expand_dropout_rate(dropout, visits, ctx)
  } else if (!is.numeric(dropout)) {
    stop(sprintf("%s: `dropout` must be numeric, a dropout_rate() object, or NULL; got %s.",
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

#' The length rule, shared for the same reason and by the same argument:
#' [as_trial_design()] used to write its own terser version of this message, so
#' a hand-built design -- the one route whose author had no constructor to
#' guide them -- got the *worse* of the two diagnoses for the identical
#' mistake. `name` carries the caller's spelling of the vector so the message
#' names what the user actually typed.
#' @rdname check_dropout_values
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
    # A `dropout_rate` reaching here means the object was assembled by hand or
    # had `$dropout` overwritten: `trial_design()` expands a rate on
    # construction, so one that survives into the stored field never went
    # through it. The field is incremental proportions by contract (CONTRACT.md
    # section 3) and expanding it now would leave the rest of the object --
    # `has_dropout` in particular -- describing something else, so this says
    # where the expansion belongs instead of doing it here.
    hint <- if (inherits(design$dropout, "dropout_rate")) {
      sprintf("\n  Pass the rate to trial_design(visits = %s, dropout = ...), which expands it.",
              fmt_call_vec(visits))
    } else {
      ""
    }
    stop(sprintf("%s: `design$dropout` must be numeric; got %s.%s",
                 context, class(dropout)[1L], hint), call. = FALSE)
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
