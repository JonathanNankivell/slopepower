# Layer 4 -- exploring a grid of candidate trial designs.
#
# Section 4.2 of Nash et al. (2021) compares nine combinations of visit schedule
# and dropout rate. Because this port separates parameter estimation from the
# sample-size calculation, the mixed model is fitted once and every design is
# evaluated against the same estimates -- the Stata original refits the model for
# each row of that table.

# ---------------------------------------------------------------------------
# dropout rates
# ---------------------------------------------------------------------------

#' A constant dropout rate per unit of time
#'
#' Dropout proportions are supplied to [trial_design()] per visit, so the same
#' underlying withdrawal rate needs a different vector for every candidate visit
#' schedule: "5% per year" over three years is `0.15` for a trial with a single
#' final visit, `rep(0.05, 3)` with annual visits, and `rep(0.025, 6)` with
#' six-monthly visits. `dropout_rate()` expresses the rate once and lets
#' the grid functions expand it correctly for each design.
#'
#' The expansion is linear in elapsed time, matching the way the worked example in
#' section 4.2 of Nash et al. (2021) is set up: the proportion whose last attended
#' visit is `visits[j]` is `rate * (visits[j + 1] - visits[j])`. Every dropout is
#' expressed as a fraction of the randomised cohort, not of those still in
#' follow-up, so the increments sum to `rate * total_duration`.
#'
#' @param rate Expected proportion of the randomised sample withdrawing per `per`
#'   units of time. Must be non-negative.
#' @param per Length of time `rate` refers to, in the units of the `time` variable
#'   used to estimate the slope parameters. Defaults to 1, i.e. `rate` is a
#'   per-unit-time rate.
#'
#' @return An object of class `dropout_rate`.
#'
#' @examples
#' dropout_rate(0.05)              # 5% per unit time
#' dropout_rate(0.10, per = 12)    # 10% per 12 months, if time is in months
#'
#' @seealso [slope_power_grid()], [slope_sample_size_grid()], [trial_design()]
#' @export
dropout_rate <- function(rate, per = 1) {
  context <- "dropout_rate()"
  check_scalar(rate, "rate", context, lower = 0, upper = Inf, lower_open = FALSE)
  check_scalar(per, "per", context, lower = 0, upper = Inf, lower_open = TRUE)
  structure(list(rate = as.numeric(rate), per = as.numeric(per)),
            class = "dropout_rate")
}

#' @describeIn dropout_rate Print a dropout rate.
#' @param x A `dropout_rate` object.
#' @param ... Ignored.
#' @export
print.dropout_rate <- function(x, ...) {
  cat(sprintf("<dropout_rate> %s per %s unit%s of time\n",
              format(x$rate), format(x$per), if (x$per == 1) "" else "s"))
  invisible(x)
}

#' Expand a dropout specification for one visit schedule
#'
#' @param spec `NULL`, a numeric vector of incremental proportions, or a
#'   `dropout_rate` object.
#' @param visits The visit times of the design being evaluated.
#' @noRd
expand_dropout <- function(spec, visits, context, label = NULL) {
  where <- if (is.null(label)) "" else sprintf(" (dropout = \"%s\")", label)

  if (is.null(spec)) return(NULL)

  if (inherits(spec, "dropout_rate")) {
    increments <- (spec$rate / spec$per) * diff(visits)
    total <- sum(increments)
    if (total > 1 + DROPOUT_TOL) {
      stop(sprintf(paste0("%s%s: a rate of %s per %s unit(s) of time over a trial lasting %s ",
                          "implies total dropout of %s, which exceeds 1."),
                   context, where, format(spec$rate), format(spec$per),
                   format(diff(range(visits))), format(total)), call. = FALSE)
    }
    return(increments)
  }

  if (is.numeric(spec)) {
    if (length(spec) != length(visits) - 1L) {
      stop(sprintf(paste0("%s%s: dropout vector has %d element(s) but the design %s has %d ",
                          "follow-up visit(s).\n  Use dropout_rate() to express a rate that ",
                          "applies across designs with different visit schedules."),
                   context, where, length(spec),
                   paste0("c(", paste(format(visits, trim = TRUE), collapse = ", "), ")"),
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

#' Render a numeric vector as a short label
#' @noRd
label_numeric <- function(x) paste(format(x, trim = TRUE, drop0trailing = TRUE), collapse = ", ")

#' Normalise `visits` to a named list of numeric vectors
#' @noRd
as_visits_list <- function(visits, context) {
  if (is.numeric(visits)) visits <- stats::setNames(list(visits), label_numeric(visits))
  if (!is.list(visits) || length(visits) == 0L) {
    stop(sprintf(paste0("%s: `visits` must be a numeric vector of visit times, or a named list ",
                        "of such vectors."), context), call. = FALSE)
  }
  nms <- names(visits)
  if (is.null(nms)) nms <- rep("", length(visits))
  blank <- !nzchar(nms)
  nms[blank] <- vapply(visits[blank], function(v) {
    if (is.numeric(v)) label_numeric(v) else "design"
  }, character(1L))
  names(visits) <- make.unique(nms, sep = "_")
  visits
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
  nms <- names(dropout)
  if (is.null(nms)) nms <- rep("", length(dropout))
  blank <- !nzchar(nms)
  nms[blank] <- vapply(dropout[blank], label_dropout, character(1L))
  names(dropout) <- make.unique(nms, sep = "_")
  dropout
}

#' @noRd
label_dropout <- function(x) {
  if (is.null(x)) return("none")
  if (inherits(x, "dropout_rate")) {
    return(sprintf("%s per %s", format(x$rate), format(x$per)))
  }
  if (is.numeric(x)) return(label_numeric(x))
  "dropout"
}

# ---------------------------------------------------------------------------
# the grid
# ---------------------------------------------------------------------------

#' Reject `effectiveness` alongside target = "observed"
#'
#' Raised in the grid wrappers rather than left to the stage-two functions,
#' because the loop below omits `effectiveness` from the call when the target is
#' the previously observed effect. Without this check the omission would bypass
#' the guard and silently return a whole table computed at effectiveness = 1.
#' @noRd
check_target_effectiveness <- function(target, supplied, context) {
  if (identical(target, "observed") && isTRUE(supplied)) {
    stop(sprintf(paste0(
      "%s: supply only one of `effectiveness` and target = \"observed\".\n",
      "  target = \"observed\" reuses the treatment effect observed in the ",
      "previous trial, which fixes effectiveness at 1."), context), call. = FALSE)
  }
  invisible(NULL)
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

  rows <- vector("list", length(visit_list) * length(drop_list))
  baseline_only <- character(0L)
  k <- 0L

  for (di in seq_along(visit_list)) {
    v <- visit_list[[di]]
    if (!is.numeric(v)) {
      stop(sprintf("%s: element \"%s\" of `visits` is not numeric.",
                   context, names(visit_list)[di]), call. = FALSE)
    }
    for (dj in seq_along(drop_list)) {
      k <- k + 1L
      dname <- names(drop_list)[dj]
      vname <- names(visit_list)[di]

      inc <- expand_dropout(drop_list[[dj]], v, context, dname)

      # trial_design() warns when the first stratum attends baseline only. That is
      # correct and expected here -- it fires for most non-zero rates -- so it is
      # collected and reported once rather than once per cell.
      des <- withCallingHandlers(
        trial_design(v, inc),
        warning = function(w) {
          if (grepl("last attended visit is", conditionMessage(w), fixed = TRUE)) {
            baseline_only <<- c(baseline_only, sprintf("%s / %s", vname, dname))
            invokeRestart("muffleWarning")
          }
        }
      )

      res <- tryCatch(evaluate(des), error = function(e) {
        stop(sprintf("%s: design \"%s\" with dropout \"%s\" failed.\n  %s",
                     context, vname, dname, conditionMessage(e)), call. = FALSE)
      })

      rows[[k]] <- data.frame(
        design        = vname,
        dropout       = dname,
        n_visits      = length(v),
        n_follow_up   = length(v) - 1L,
        last_visit    = v[length(v)],
        dropout_total = sum(des$dropout),
        n             = res$n,
        n_per_arm     = res$n_per_arm,
        power         = res$power,
        tte           = res$tte,
        var_tte       = res$var_tte,
        effect_size   = res$effect_size,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(baseline_only) > 0L) {
    warning(sprintf(paste0("%s: in %d of %d combinations a non-zero proportion is expected to ",
                           "attend the baseline visit only (%s). Those participants provide no ",
                           "follow-up measurement and so contribute nothing to the comparison of ",
                           "slopes."),
                    context, length(baseline_only), k,
                    paste(baseline_only, collapse = "; ")), call. = FALSE)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
#' expanded correctly for each schedule.
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
#' pars <- slope_params_manual(
#'   slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
#'   sigma_cov = 5, sigma2_residual = 10
#' )
#' slope_power_grid(
#'   pars, n = 450, effectiveness = 0.33,
#'   visits  = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05), `10pc` = dropout_rate(0.10))
#' )
#'
#' @seealso [slope_power()], [slope_sample_size_grid()], [dropout_rate()]
#' @export
slope_power_grid <- function(params, visits, dropout = NULL, n,
                             effectiveness = 0.25,
                             target = c("effectiveness", "observed"),
                             alpha = 0.05) {
  context <- "slope_power_grid()"
  if (missing(n)) {
    stop(sprintf(paste0(
      "%s: `n` is required -- this grid holds the sample size fixed and reports\n",
      "  the power each design achieves. For the sample size each design needs,\n",
      "  use slope_sample_size_grid()."), context), call. = FALSE)
  }
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)

  args <- list(params = params, n = n, target = target, alpha = alpha)
  if (!identical(target, "observed")) args$effectiveness <- effectiveness

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
#' @inheritParams slope_power_grid
#' @param power Target power, held constant across the grid. Defaults to 0.8.
#' @param effectiveness,target,alpha Passed to [slope_sample_size()] and held
#'   constant across the grid.
#'
#' @return A data frame with the same columns as [slope_power_grid()]. `power` is
#'   constant; `n` and `n_per_arm` are what vary.
#'
#' @examples
#' pars <- slope_params_manual(
#'   slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
#'   sigma_cov = 5, sigma2_residual = 10
#' )
#' slope_sample_size_grid(
#'   pars, power = 0.8, effectiveness = 0.33,
#'   visits  = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05))
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

  args <- list(params = params, power = power, target = target, alpha = alpha)
  if (!identical(target, "observed")) args$effectiveness <- effectiveness

  grid_impl(visits, dropout,
            function(des) do.call(slope_sample_size, c(args, list(design = des))),
            context)
}
