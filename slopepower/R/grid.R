# Layer 4 -- exploring a grid of candidate trial designs.
#
# Section 4.2 of Nash et al. (2021) compares nine combinations of visit schedule
# and dropout rate. Because this port separates parameter estimation from the
# sample-size calculation, the mixed model is fitted once and every design is
# evaluated against the same estimates -- the Stata original refits the model for
# each row of that table.

# ---------------------------------------------------------------------------
# dropout specifications
#
# `dropout_rate()` itself lives in design.R, beside the `trial_design()`
# argument it is an alternative spelling of. Only the grid's own handling of a
# specification -- naming the failing cell, and the length rule that a grid
# makes newly relevant -- is here.
# ---------------------------------------------------------------------------

#' Expand a dropout specification for one visit schedule
#'
#' The grid's wrapper around `expand_dropout_rate()` (design.R): it adds the
#' cell label to any message, and applies the length rule to a bare numeric
#' vector *before* `trial_design()` sees it, because a length mismatch means
#' something specific in a grid -- a fixed vector paired with a schedule it was
#' not written for -- that a single design cannot exhibit.
#'
#' @param spec `NULL`, a numeric vector of incremental proportions, or a
#'   `dropout_rate` object.
#' @param visits The visit times of the design being evaluated.
#' @param label The grid label of the dropout specification, named in every
#'   error message so the failing cell can be identified.
#' @noRd
expand_dropout <- function(spec, visits, context, label) {
  where <- sprintf(" (dropout = \"%s\")", label)

  if (is.null(spec)) return(NULL)

  if (inherits(spec, "dropout_rate")) {
    return(expand_dropout_rate(spec, visits, context, where))
  }

  if (is.numeric(spec)) {
    if (length(spec) != length(visits) - 1L) {
      stop(sprintf(paste0("%s%s: dropout vector has %d element(s) but the design %s has %d ",
                          "follow-up visit(s).\n  Use dropout_rate() to express a rate that ",
                          "applies across designs with different visit schedules."),
                   context, where, length(spec), fmt_call_vec(visits),
                   length(visits) - 1L), call. = FALSE)
    }
    return(spec)
  }

  stop(sprintf("%s%s: dropout must be NULL, a numeric vector, or a dropout_rate() object; got %s.",
               context, where, class(spec)[1L]), call. = FALSE)
}

# ---------------------------------------------------------------------------
# argument normalisation
# ---------------------------------------------------------------------------

#' Fill in any blank names of a list via `label`, then dedupe
#'
#' The common tail of `as_visits_list()` and `as_dropout_list()`, once each has
#' normalised its argument to a (possibly partly-named) list.
#' @noRd
fill_list_names <- function(x, label) {
  nms <- names(x) %||% rep("", length(x))
  blank <- !nzchar(nms)
  nms[blank] <- vapply(x[blank], label, character(1L))
  names(x) <- make.unique(nms, sep = "_")
  x
}

#' Normalise `visits` to a named list of numeric vectors
#' @noRd
as_visits_list <- function(visits, context) {
  if (is.numeric(visits)) visits <- stats::setNames(list(visits), label_numeric(visits))
  if (!is.list(visits) || length(visits) == 0L) {
    stop(sprintf(paste0("%s: `visits` must be a numeric vector of visit times, or a named list ",
                        "of such vectors."), context), call. = FALSE)
  }
  fill_list_names(visits, label_visits)
}

#' Normalise `dropout` to a named list of specifications
#' @noRd
as_dropout_list <- function(dropout, context) {
  if (is.null(dropout)) return(stats::setNames(list(NULL), "none"))
  if (is.numeric(dropout) || inherits(dropout, "dropout_rate")) {
    dropout <- stats::setNames(list(dropout), label_dropout(dropout))
  }
  if (!is.list(dropout) || length(dropout) == 0L) {
    stop(sprintf(paste0("%s: `dropout` must be NULL, a numeric vector, a dropout_rate() object, ",
                        "or a named list of these."), context), call. = FALSE)
  }
  fill_list_names(dropout, label_dropout)
}

#' Name an unnamed grid element after its own contents
#'
#' One labeller per axis, so that `as_visits_list()` and `as_dropout_list()`
#' read the same way. Both tolerate an element they cannot describe: a
#' non-numeric `visits` entry is rejected by `grid_impl()` with a message that
#' names the element, which needs the element to have a name first.
#' @noRd
label_visits <- function(x) if (is.numeric(x)) label_numeric(x) else "design"

#' @rdname label_visits
#' @noRd
label_dropout <- function(x) {
  if (is.null(x)) return("none")
  if (inherits(x, "dropout_rate")) {
    return(sprintf("%s per %s", fmt_num(x$rate), fmt_num(x$per)))
  }
  if (is.numeric(x)) return(label_numeric(x))
  "dropout"
}

# ---------------------------------------------------------------------------
# the grid
# ---------------------------------------------------------------------------

#' Wrap an error from one grid cell with the cell that produced it
#'
#' Shared by the `trial_design()` call and the `evaluate()` call in
#' `grid_impl()`'s loop below, so a design/dropout combination that fails
#' either step is reported the same way: named, rather than surfacing whatever
#' unqualified message the failing call happens to raise.
#' @noRd
grid_cell_error <- function(context, vname, dname) {
  function(e) {
    stop(sprintf("%s: design \"%s\" with dropout \"%s\" failed.\n  %s",
                 context, vname, dname, conditionMessage(e)), call. = FALSE)
  }
}

#' Walk the visits x dropout cross product, evaluating one design per cell
#'
#' `evaluate` takes a `trial_design` and returns a `slope_sample_size` or
#' `slope_power` object. Everything except that call is shared between the two
#' grids, so the two tables are guaranteed to have the same shape and the same
#' dropout expansion.
#' @noRd
grid_impl <- function(visits, dropout, evaluate, context) {
  visit_list <- as_visits_list(visits, context)
  drop_list <- as_dropout_list(dropout, context)

  n_cells <- length(visit_list) * length(drop_list)
  # Twelve preallocated columns, filled in place and made a data frame once at
  # the end. Assembling a one-row data frame per cell and rbind()ing the list
  # cost more than the calculation the grid exists to perform -- about 2 ms a
  # cell against 0.6 ms for the sample size itself, so roughly 60% of a grid's
  # total runtime went on building the table rather than filling it.
  out <- list(
    design        = character(n_cells),
    dropout       = character(n_cells),
    n_visits      = integer(n_cells),
    n_follow_up   = integer(n_cells),
    last_visit    = numeric(n_cells),
    dropout_total = numeric(n_cells),
    n             = numeric(n_cells),
    n_per_arm     = numeric(n_cells),
    power         = numeric(n_cells),
    tte           = numeric(n_cells),
    var_tte       = numeric(n_cells),
    effect_size   = numeric(n_cells)
  )
  baseline_only <- character(0L)
  tte_direction <- character(0L)
  k <- 0L

  for (di in seq_along(visit_list)) {
    v <- visit_list[[di]]
    vname <- names(visit_list)[di]
    if (!is.numeric(v)) {
      stop(sprintf("%s: element \"%s\" of `visits` is not numeric.",
                   context, vname), call. = FALSE)
    }
    for (dj in seq_along(drop_list)) {
      k <- k + 1L
      dname <- names(drop_list)[dj]

      inc <- expand_dropout(drop_list[[dj]], v, context, dname)

      # trial_design() warns when the first stratum attends baseline only. That is
      # correct and expected here -- it fires for most non-zero rates -- so it is
      # collected and reported once rather than once per cell. Matched by
      # condition class, not message text, so a copy-edit of the warning's
      # wording in design.R cannot silently break this. An invalid combination
      # (e.g. a `visits` element not starting at 0) is caught here too, and
      # named the same way a failure from `evaluate()` below is.
      des <- tryCatch(
        withCallingHandlers(
          trial_design(v, inc),
          slopepower_baseline_dropout = function(w) {
            baseline_only <<- c(baseline_only, sprintf("%s / %s", vname, dname))
            invokeRestart("muffleWarning")
          }
        ),
        error = grid_cell_error(context, vname, dname)
      )

      # effect_components() warns, by the same class mechanism, when the
      # target treatment effect makes the slope more extreme rather than less.
      # Collected and reported once here too, rather than once per cell.
      res <- tryCatch(
        withCallingHandlers(
          evaluate(des),
          slopepower_tte_direction = function(w) {
            tte_direction <<- c(tte_direction, sprintf("%s / %s", vname, dname))
            invokeRestart("muffleWarning")
          }
        ),
        error = grid_cell_error(context, vname, dname)
      )

      out$design[k]        <- vname
      out$dropout[k]       <- dname
      out$n_visits[k]      <- length(v)
      out$n_follow_up[k]   <- length(v) - 1L
      out$last_visit[k]    <- v[length(v)]
      out$dropout_total[k] <- sum(des$dropout)

      # The six result columns are read straight off `res`. That is the same
      # thing as.data.frame.slope_result() reports -- it copies each of these
      # six from `x` unchanged -- and routing through it instead was the single
      # most expensive line in the loop, at 1.3 ms a cell. The guarantee that
      # the two cannot drift is kept by test-grid.R, which pins a grid row
      # against as.data.frame() of the same object, rather than by paying for
      # the other twelve columns and discarding them.
      out$n[k]           <- res$n
      out$n_per_arm[k]   <- res$n_per_arm
      out$power[k]       <- res$power
      out$tte[k]         <- res$tte
      out$var_tte[k]     <- res$var_tte
      out$effect_size[k] <- res$effect_size
    }
  }

  report_collected(context, baseline_only, k,
                   paste0("a non-zero proportion is expected to attend the baseline visit ",
                          "only (%s). Those participants provide no follow-up measurement ",
                          "and so contribute nothing to the comparison of slopes."))
  report_collected(context, tte_direction, k,
                   paste0("the target treatment effect makes the slope more extreme (%s). ",
                          "The comparator slope is further from zero than the group being ",
                          "treated; check that `slope_comparator` is the intended target."))

  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Report one class of collected per-cell warning, once for the whole grid
#'
#' Both collectors in [grid_impl()] owe the user the same sentence -- how many
#' of how many combinations, and which -- and differ only in what they then say
#' about it. Writing the preamble twice let the two drift; `tail` supplies just
#' the clause that differs, with a single `%s` where the cell list goes.
#' @noRd
report_collected <- function(context, cells, k, tail) {
  if (length(cells) == 0L) return(invisible(NULL))
  warning(sprintf("%s: in %d of %d combinations %s", context, length(cells), k,
                  sprintf(tail, paste(cells, collapse = "; "))), call. = FALSE)
}

#' Compare the power of many candidate trial designs
#'
#' Evaluates the cross product of a set of visit schedules and a set of dropout
#' assumptions against a single set of stage-one parameter estimates, at one fixed
#' sample size, returning one row per combination. This is the calculation behind
#' table 1 of Nash et al. (2021), which explores how adding interim visits changes
#' the power of a trial under different rates of loss to follow-up.
#'
#' Because dropout proportions are supplied per visit, a single numeric vector
#' cannot describe the same withdrawal behaviour across designs with different
#' numbers of visits. Use [dropout_rate()] to state the rate once and have it
#' expanded correctly for each schedule. However the proportions are arrived at,
#' each cell of the grid handles them exactly as [slope_power()] does, by the
#' Dawson and Lagakos (1991, 1993) pattern mixture described in [trial_design()].
#'
#' Use [slope_sample_size_grid()] for the converse table: the sample size each
#' design needs to reach a target power.
#'
#' @param params A `slope_params` object.
#' @param visits A numeric vector of visit times, or a named list of such vectors.
#'   Each must begin at 0.
#' @param dropout `NULL`, a numeric vector of incremental proportions, a
#'   [dropout_rate()] object, or a named list mixing any of these.
#' @param n Total number of participants, held constant across the grid.
#'   Required.
#' @param effectiveness,target,alpha Passed to [slope_power()] and held constant
#'   across the grid.
#'
#' @return A data frame with one row per design and dropout combination, with
#'   columns `design`, `dropout`, `n_visits`, `n_follow_up`, `last_visit`,
#'   `dropout_total`, `n`, `n_per_arm`, `power`, `tte`, `var_tte` and
#'   `effect_size`. `n` is constant; `power` is what varies.
#'
#' @examples
#' # Table 1, p.595 of Nash et al. (2021): explores how adding interim visits
#' # changes power under three dropout rates, holding the sample size at the
#' # trial's planned N = 450. No comparator: measured toward a slope of zero.
#' # The stage-one estimates are the exact fit reported for `slpower1` (see
#' # `slope_params(sdmt ~ visit | id, data = slpower1)`), passed here via
#' # slope_params_manual() so the example does not have to refit the mixed
#' # model on all two hundred participants.
#' pars <- slope_params_manual(
#'   slope = -1.6725, sigma2_intercept = 111.786636, sigma2_slope = 2.350021,
#'   sigma_cov = 2.810881, sigma2_residual = 9.159780
#' )
#' slope_power_grid(
#'   pars, n = 450, effectiveness = 0.33,
#'   visits  = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05), `10pc` = dropout_rate(0.10))
#' )
#'
#' # Case/healthy-control comparator: two cases and two healthy controls, a
#' # subset of `slpower2` whose visits are calendar dates converted to years.
#' df2 <- slpower2[slpower2$id %in% c(1, 2, 251, 252), ]
#' pars2 <- slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, data = df2, healthy = case)
#' slope_power_grid(
#'   pars2, n = 40, effectiveness = 0.33,
#'   visits  = list(annual = c(0, 1, 2), six_month = seq(0, 2, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05))
#' )
#'
#' # Randomised-trial comparator, target = "observed": compare a design with
#' # only the trial's original two follow-up visits against a denser one, at
#' # the effect size the trial actually found, fitted to all one hundred and
#' # fifty participants of `slpower3`.
#' pars3 <- slope_params(sdmt ~ visit | id, data = slpower3, treated = treat)
#' slope_power_grid(
#'   pars3, n = 396, target = "observed",
#'   visits  = list(as_planned = c(0, 0.5, 2), denser = c(0, 0.5, 1, 2)),
#'   dropout = NULL
#' )
#'
#' @seealso [slope_power()], [slope_sample_size_grid()], [dropout_rate()]
#' @export
slope_power_grid <- function(params, visits, dropout = NULL, n,
                             effectiveness = 0.25,
                             target = c("effectiveness", "observed"),
                             alpha = 0.05) {
  context <- "slope_power_grid()"
  # `is.null(n)` too; see the note on the same guard in slope_power(). Without
  # it the failure surfaces from inside the loop, wrapped in a per-cell message,
  # and complains about `power` rather than the missing `n`.
  if (missing(n) || is.null(n)) {
    stop(sprintf(paste0(
      "%s: `n` is required -- this grid holds the sample size fixed and reports\n",
      "  the power each design achieves. For the sample size each design needs,\n",
      "  use slope_sample_size_grid()."), context), call. = FALSE)
  }
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)

  args <- maybe_add_effectiveness(list(params = params, n = n, target = target, alpha = alpha),
                                  effectiveness, target)

  grid_impl(visits, dropout,
            function(des) do.call(slope_power, c(args, list(design = des))),
            context)
}

#' Compare the sample size many candidate trial designs need
#'
#' The companion to [slope_power_grid()]: holds the target power fixed and reports
#' the sample size each combination of visit schedule and dropout assumption
#' requires. Section 4.2 of Nash et al. (2021) is the power version of this table;
#' this is the same exploration read the other way round.
#'
#' Dropout is handled cell by cell exactly as in [slope_sample_size()], by the
#' Dawson and Lagakos (1991, 1993) pattern mixture described in [trial_design()].
#' Use [dropout_rate()] rather than a fixed vector when the schedules being
#' compared have different numbers of visits.
#'
#' @inheritParams slope_power_grid
#' @param power Target power, held constant across the grid. Defaults to 0.8.
#' @param effectiveness,target,alpha Passed to [slope_sample_size()] and held
#'   constant across the grid.
#'
#' @return A data frame with the same columns as [slope_power_grid()]. `power` is
#'   constant; `n` and `n_per_arm` are what vary.
#'
#' @examples
#' # No comparator: measured toward a slope of zero, fitted to all two
#' # hundred participants of `slpower1`.
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#' slope_sample_size_grid(
#'   pars, power = 0.8, effectiveness = 0.33,
#'   visits  = list(annual = c(0, 1, 2, 3), six_month = seq(0, 3, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05))
#' )
#'
#' # Case/healthy-control comparator: two cases and two healthy controls, a
#' # subset of `slpower2` whose visits are calendar dates converted to years.
#' df2 <- slpower2[slpower2$id %in% c(1, 2, 251, 252), ]
#' pars2 <- slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, data = df2, healthy = case)
#' slope_sample_size_grid(
#'   pars2, power = 0.8, effectiveness = 0.33,
#'   visits  = list(annual = c(0, 1, 2), six_month = seq(0, 2, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05))
#' )
#'
#' # Randomised-trial comparator, target = "observed": the sample size needed
#' # to detect the same effect the trial actually found, comparing its
#' # original visit schedule against a denser one, fitted to all one hundred
#' # and fifty participants of `slpower3`.
#' pars3 <- slope_params(sdmt ~ visit | id, data = slpower3, treated = treat)
#' slope_sample_size_grid(
#'   pars3, power = 0.8, target = "observed",
#'   visits  = list(as_planned = c(0, 0.5, 2), denser = c(0, 0.5, 1, 2)),
#'   dropout = NULL
#' )
#'
#' @seealso [slope_sample_size()], [slope_power_grid()], [dropout_rate()]
#' @export
slope_sample_size_grid <- function(params, visits, dropout = NULL, power = 0.8,
                                   effectiveness = 0.25,
                                   target = c("effectiveness", "observed"),
                                   alpha = 0.05) {
  context <- "slope_sample_size_grid()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)

  args <- maybe_add_effectiveness(
    list(params = params, power = power, target = target, alpha = alpha),
    effectiveness, target)

  grid_impl(visits, dropout,
            function(des) do.call(slope_sample_size, c(args, list(design = des))),
            context)
}
