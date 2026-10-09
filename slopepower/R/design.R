# Layer 2 --- the proposed trial design: visit schedule and dropout pattern.

#' Render a numeric vector as the R call that would create it
#' @noRd
fmt_call_vec <- function(x) {
  if (length(x) == 1L) fmt_num(x) else paste0("c(", paste(fmt_num(x), collapse = ", "), ")")
}

# ---------------------------------------------------------------------------
# dropout rates
#
# `dropout_rate()` lives beside the design validators rather than with the grid
# functions that were once its only consumer, because it is an alternative way
# of saying what every stage-two function's `dropout` argument means, not a
# feature of the grid. It was originally defined in grid.R, and the consequence
# was that the object could only be expanded by the grid: passing one to a
# single-design calculation -- the call its own documentation pointed at -- was
# rejected as a non-numeric dropout.
# ---------------------------------------------------------------------------

#' A constant dropout rate per unit of time
#'
#' Dropout proportions are otherwise supplied per visit, so the same underlying
#' withdrawal rate needs a different vector for every candidate visit
#' schedule: "5% per year" over three years is `0.15` for a trial with a single
#' final visit, `rep(0.05, 3)` with annual visits, and `rep(0.025, 6)` with
#' six-monthly visits. `dropout_rate()` expresses the rate once and is expanded
#' correctly for whichever schedule it is paired with --- which is what makes it
#' useful to the grid functions, where one object drives every row of a table
#' of competing schedules. It is accepted as the `dropout` argument of every
#' stage-two function.
#'
#' `pattern` chooses which of two withdrawal patterns `rate` describes:
#'
#' * `"linear"` (the default) applies `rate` to the *original randomised
#'   cohort*, matching the worked example in section 4.2 of Nash et al. (2021):
#'   the proportion whose last attended visit is `visits[j]` is
#'   `rate * (visits[j + 1] - visits[j]) / per`. The same number withdraw in
#'   every equal-length interval regardless of how many have already left, so
#'   the proportion still in follow-up falls linearly, and the increments sum
#'   to `rate * total_duration / per`.
#' * `"geometric"` applies `rate` as the proportion of *whoever is still in
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
#' Both patterns expand to incremental proportions by construction, so a
#' `dropout_rate` cannot be combined with `dropout_scale = "cumulative"`, which
#' describes how a numeric `dropout` *vector* was written down. The two are
#' named apart deliberately: `pattern` is how withdrawal behaves over time,
#' `dropout_scale` how a vector of proportions is to be read.
#'
#' A rate always has some participants withdrawing before the first follow-up
#' visit, and so attending baseline only. They carry no slope information and
#' are left out of the calculation, as for any dropout; but since that share is
#' an inevitable consequence of the rate rather than a number anyone wrote, it
#' is not warned about, as a `dropout` vector with a non-zero first element is.
#'
#' This object only produces the per-visit proportions. What the calculation then
#' does with them --- the Dawson and Lagakos (1991, 1993) pattern mixture, and
#' what it assumes about why people withdraw --- is described in the "Dropout"
#' section of [slope_sample_size()].
#'
#' @param rate Expected proportion withdrawing per `per` units of time. For
#'   `pattern = "linear"`, a proportion of the original randomised sample and so
#'   non-negative with no upper bound of its own (the total across all visits
#'   is what is capped at 1). For `pattern = "geometric"`, a proportion of
#'   whoever remains and so restricted to `[0, 1]`.
#' @param per Length of time `rate` refers to, in the units of the `time` variable
#'   used to estimate the slope parameters. Defaults to 1, i.e. `rate` is a
#'   per-unit-time rate.
#' @param pattern Which withdrawal pattern `rate` describes: `"linear"` (the
#'   default) or `"geometric"`. See Details.
#'
#' @return An object of class `dropout_rate`.
#'
#' @examples
#' dropout_rate(0.05)              # 5% of the original cohort per unit time
#' dropout_rate(0.10, per = 12)    # 10% per 12 months, if time is in months
#' dropout_rate(0.05, pattern = "geometric")  # 5% of those remaining, per unit time
#'
#' # The same rate, priced against two different schedules
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sample_size(pars, c(0, 1, 2, 3), dropout = dropout_rate(0.05))
#' slope_sample_size(pars, c(0, 1.5, 3), dropout = dropout_rate(0.05))
#'
#' @seealso [slope_sample_size()], [slope_power()], [slope_power_grid()],
#'   [slope_sample_size_grid()]
#' @export
dropout_rate <- function(rate, per = 1, pattern = c("linear", "geometric")) {
  context <- "dropout_rate()"
  pattern <- match.arg(pattern)
  if (identical(pattern, "geometric")) {
    # A geometric rate is a proportion of the remaining sample and is the base
    # of a power in expand_dropout_rate(); anything above 1 would describe more
    # than everyone remaining leaving, and a non-integer exponent of a negative
    # base is complex, not a dropout proportion.
    check_scalar(rate, "rate", context, lower = 0, upper = 1, lower_open = FALSE,
                 upper_open = FALSE)
  } else {
    check_scalar(rate, "rate", context, lower = 0, upper = Inf, lower_open = FALSE)
  }
  check_scalar(per, "per", context, lower = 0, upper = Inf, lower_open = TRUE)
  structure(list(rate = as.numeric(rate), per = as.numeric(per), pattern = pattern),
            class = "dropout_rate")
}

#' @describeIn dropout_rate Print a dropout rate.
#' @param x A `dropout_rate` object.
#' @param ... Ignored.
#' @export
print.dropout_rate <- function(x, ...) {
  cat(sprintf("<dropout_rate> %s per %s unit%s of time (%s)\n",
              fmt_num(x$rate), fmt_num(x$per), if (x$per == 1) "" else "s", x$pattern))
  invisible(x)
}

#' Expand a `dropout_rate` into incremental proportions for one schedule
#'
#' The single place a rate becomes numbers, called by [validate_dropout()] for
#' every design, single or in a grid. A grid names the failing cell around the
#' whole error, in `grid_cell_error()`.
#'
#' The two `pattern`s (see [dropout_rate()]) compute the incremental proportions
#' differently: `"linear"` applies `rate` to the original cohort, so the total
#' across the whole trial can exceed 1 and is checked for that here, with a
#' diagnosis in terms of `rate` and `per` rather than the expanded vector.
#' `"geometric"` applies `rate` to whoever remains, so the total is
#' mathematically bounded below 1 and no equivalent check is needed. Both
#' callers run `increments` through `check_dropout_total()` afterwards
#' regardless, so a bug in this earlier, friendlier check could not let an
#' invalid total through uncaught.
#' @noRd
expand_dropout_rate <- function(spec, visits, ctx) {
  if (identical(spec$pattern, "geometric")) {
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
    stop(sprintf(paste0("%s: a rate of %s per %s unit(s) of time over a trial lasting %s ",
                        "implies total dropout of %s, which exceeds 1."),
                 ctx, fmt_num(spec$rate), fmt_num(spec$per),
                 fmt_num(diff(range(visits))), fmt_num(total)), call. = FALSE)
  }
  increments
}

#' Validate a proposed trial's visit schedule and dropout pattern
#'
#' Every stage-two entry point takes `visits`, `dropout` and `dropout_scale`
#' directly and passes them here, so the design of the *future* trial is
#' checked by one set of rules whichever function the caller reached it
#' through. The result is an internal `trial_design` object -- `visits`,
#' `dropout` (always incremental), `has_dropout` and `dropout_scale`, see
#' CONTRACT.md section 3 -- which the calculation reads and which every result
#' carries as its `design` field. `ctx` names the exported function the caller
#' typed, for its errors.
#' @noRd
build_trial_design <- function(visits, dropout, dropout_scale, ctx) {
  visits <- validate_visits(visits, ctx)
  n_intervals <- length(visits) - 1L

  # A rate puts someone in the baseline-only stratum whenever it is non-zero,
  # so warning would only say that a rate is a rate; see ?dropout_rate. The
  # warning is for a first element the caller wrote.
  from_rate <- inherits(dropout, "dropout_rate")
  dropout <- validate_dropout(dropout, n_intervals, dropout_scale, visits, ctx)
  if (!from_rate) warn_baseline_dropout(dropout, visits, ctx)

  structure(
    list(
      visits        = visits,
      dropout       = dropout,
      has_dropout   = any(dropout > 0),
      dropout_scale = dropout_scale
    ),
    class = "trial_design"
  )
}

#' Warn when the first dropout stratum attends only the baseline visit
#'
#' Kept apart from [build_trial_design()] so the condition class, which the
#' grid functions catch by name, is defined in one obvious place.
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
validate_dropout <- function(dropout, n_intervals, dropout_scale, visits, ctx) {
  if (is.null(dropout)) {
    return(rep(0, n_intervals))
  }

  if (inherits(dropout, "dropout_rate")) {
    # Rejected rather than ignored. A rate expands to incremental proportions
    # by construction, so there is no sense in which it could have been
    # supplied cumulatively; accepting the pair would store incremental values
    # in a design whose `dropout_scale` -- and whose print method -- announce
    # them as cumulative.
    if (identical(dropout_scale, "cumulative")) {
      stop(sprintf(paste0("%s: a dropout_rate() cannot be combined with dropout_scale = ",
                          "\"cumulative\"; it expands to the proportion withdrawing within ",
                          "each interval, which is already incremental. Drop the ",
                          "`dropout_scale` argument, or supply the cumulative proportions ",
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

  if (identical(dropout_scale, "cumulative")) {
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

#' The length rule. `name` carries the caller's spelling of the vector so the
#' message names what the user actually typed.
#' @rdname check_dropout_values
#' @noRd
check_dropout_length <- function(dropout, visits, name, ctx) {
  n_intervals <- length(visits) - 1L
  if (length(dropout) != n_intervals) {
    stop(sprintf(paste0("%s: `%s` must have one element per follow-up visit: ",
                        "length(visits) - 1 = %d, but length(%s) = %d.",
                        "\n  visits = %s covers %d follow-up visit%s after baseline.",
                        "\n  Use dropout_rate() to express a rate that applies across ",
                        "schedules of different lengths."),
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
                        "dropout_scale = \"cumulative\"."),
                 ctx, name, fmt_num(total)), call. = FALSE)
  }
  invisible(dropout)
}

#' Print the design a stage-two result carries
#'
#' Not a constructor anyone calls -- the design is built from `visits` and
#' `dropout` by the stage-two functions -- but every result keeps it as
#' `$design`, and printing that should show both readings of the dropout
#' vector rather than a bare list.
#' @param x A `trial_design` object.
#' @param ... Ignored.
#' @noRd
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

  cat(sprintf("  Dropout:     supplied as %s\n\n", x$dropout_scale))
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
