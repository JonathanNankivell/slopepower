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

#' Normalise a scalar grid argument to a named list of single values
#'
#' The third kind of axis. `visits` and `dropout` describe the trial itself and
#' each needs a reader of its own; every other argument the grid used to hold
#' constant -- `effectiveness`, `alpha`, and whichever of `n` and `power` the
#' grid is not solving for -- is one number per cell, so a single reader serves
#' them all. A bare vector is as good a spelling of an axis as a list is, since
#' none of these arguments accepts a vector of its own.
#'
#' The elements themselves are left to the stage-two call to validate:
#' `check_scalar()` there already states each argument's own bounds, and
#' [grid_cell_error()] names the cell whose value failed them.
#' @noRd
as_scalar_list <- function(x, arg, context) {
  if (is.numeric(x)) x <- as.list(x)
  if (!is.list(x) || length(x) == 0L) {
    stop(sprintf(paste0("%s: `%s` must be a number, or a numeric vector or named list of ",
                        "numbers -- one per level of the sensitivity analysis."),
                 context, arg), call. = FALSE)
  }
  fill_list_names(x, label_scalar)
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
    # `type` is carried into the label, like `per`, only when it is not the
    # default: "linear" is what every existing rate label already means, so
    # leaving it off keeps those labels unchanged, while a cumulative rate --
    # a different assumption from a linear one at the same rate and per --
    # would otherwise key the same row.
    suffix <- if (identical(x$type, "cumulative")) ", cumulative" else ""
    return(sprintf("%s per %s%s", fmt_num(x$rate), fmt_num(x$per), suffix))
  }
  if (is.numeric(x)) return(label_numeric(x))
  "dropout"
}

#' @rdname label_visits
#' @noRd
label_scalar <- function(x) if (is.numeric(x) && length(x) == 1L) fmt_num(x) else "value"

# ---------------------------------------------------------------------------
# the grid
# ---------------------------------------------------------------------------

#' Wrap an error from one grid cell with the cell that produced it
#'
#' Shared by the `trial_design()` call and the `evaluate()` call in
#' `grid_impl()`'s loops below, so a cell that fails either step is reported the
#' same way: named, rather than surfacing whatever unqualified message the
#' failing call happens to raise. The cell arrives as its coordinate --- one
#' named element per axis --- so an axis added to the grid names itself here
#' without this function having to learn about it.
#' @noRd
grid_cell_error <- function(context, labels) {
  function(e) {
    stop(sprintf("%s: %s failed.\n  %s", context, cell_label(labels),
                 conditionMessage(e)), call. = FALSE)
  }
}

#' Name one cell of the grid by its position on every axis
#' @noRd
cell_label <- function(labels) {
  paste(sprintf('%s "%s"', names(labels), labels), collapse = ", ")
}

#' The cross product of the axes, in nested-loop order
#'
#' One row per cell, one column per axis, holding an index into that axis's
#' list. `expand.grid()` varies its *first* argument fastest, which is the
#' opposite of how a reader expects a table to nest, so the axes go in reversed
#' and the columns come back out reversed: the first axis varies slowest.
#' @noRd
grid_cells <- function(axes) {
  rev(expand.grid(lapply(rev(axes), seq_along), KEEP.OUT.ATTRS = FALSE))
}

#' Normalise a grid's axes and build the columns that describe a cell
#'
#' The half of what used to be `grid_impl()` that has nothing to do with
#' solving a cell: normalising the three kinds of axis, building one
#' `trial_design` per visits/dropout pair (collecting the baseline-dropout
#' warning as it goes, since that is a property of the design, not of what is
#' later computed from it), the cross product of every axis, and the six
#' columns -- `design` through `dropout_total` -- that a cell has independently
#' of any calculation.
#'
#' Split out from the solving loop, now [grid_evaluate()], so that
#' [slope_sample_size_grid_boot()] can build this once and read the *design*
#' every cell needs -- `g$designs[[g$design_of[k]]]` -- to build its own
#' per-cell closures, without [grid_evaluate()]'s plain stage-two solve running
#' first and being thrown away. `grid_impl()` below still calls both, in the
#' same order, for the same two grids this always served.
#'
#' `scalars` is a named list of the remaining axes -- `effectiveness`, `alpha`,
#' and whichever of `n` and `power` the grid is not solving for -- each a single
#' value or several. They nest inside `visits` and `dropout`, which vary
#' slowest, so a grid that varies nothing else lists its cells in exactly the
#' order the visits x dropout loop always did.
#' @noRd
grid_axes <- function(visits, dropout, scalars, context) {
  visit_list <- as_visits_list(visits, context)
  drop_list  <- as_dropout_list(dropout, context)
  scalar_lists <- stats::setNames(
    lapply(names(scalars), function(nm) as_scalar_list(scalars[[nm]], nm, context)),
    names(scalars))

  n_designs <- length(visit_list) * length(drop_list)
  baseline_only <- character(0L)

  # The designs are built once per visits/dropout pair rather than once per
  # cell. A sensitivity axis re-prices the same trial, so rebuilding it at every
  # level of one would both cost more and report the baseline-dropout warning
  # below once per level of an axis that has nothing to do with it.
  designs <- vector("list", n_designs)
  for (di in seq_along(visit_list)) {
    v <- visit_list[[di]]
    vname <- names(visit_list)[di]
    if (!is.numeric(v)) {
      stop(sprintf("%s: element \"%s\" of `visits` is not numeric.",
                   context, vname), call. = FALSE)
    }
    for (dj in seq_along(drop_list)) {
      dname <- names(drop_list)[dj]

      inc <- expand_dropout(drop_list[[dj]], v, context, dname)

      # trial_design() warns when the first stratum attends baseline only. That is
      # correct and expected here -- it fires for most non-zero rates -- so it is
      # collected and reported once rather than once per design. Matched by
      # condition class, not message text, so a copy-edit of the warning's
      # wording in design.R cannot silently break this. An invalid combination
      # (e.g. a `visits` element not starting at 0) is caught here too, and
      # named the same way a failure from `evaluate()` in grid_evaluate() is.
      designs[[(di - 1L) * length(drop_list) + dj]] <- tryCatch(
        withCallingHandlers(
          trial_design(v, inc),
          slopepower_baseline_dropout = function(w) {
            baseline_only <<- c(baseline_only, sprintf("%s / %s", vname, dname))
            invokeRestart("muffleWarning")
          }
        ),
        error = grid_cell_error(context, c(design = vname, dropout = dname))
      )
    }
  }
  report_collected(context, baseline_only, n_designs,
                   paste0("a non-zero proportion is expected to attend the baseline visit ",
                          "only (%s). Those participants provide no follow-up measurement ",
                          "and so contribute nothing to the comparison of slopes."))

  axes <- c(list(design = visit_list, dropout = drop_list), scalar_lists)
  cells <- grid_cells(axes)
  n_cells <- nrow(cells)

  # Which axes a message has to name to identify a cell. The two design axes
  # always do -- they are what the grid is fundamentally a table of, and are
  # named even when there is only one of each, as they always were. A scalar
  # axis earns its place in the message only by having more than one level:
  # otherwise it is the same constant in every cell, and naming it would push
  # the message that actually matters further along the line.
  named <- lengths(axes) > 1L
  named[c("design", "dropout")] <- TRUE

  # Called from inside grid_evaluate()'s condition handlers, and nowhere else.
  # Lazy argument evaluation is what keeps it that way: `labels` is a promise
  # inside grid_cell_error()'s closure, forced only if the cell actually fails,
  # so offering a message that no cell will need costs nothing per cell.
  labels_at <- function(k) {
    stats::setNames(
      vapply(seq_along(axes), function(i) names(axes[[i]])[cells[[i]][k]], character(1L)),
      names(axes))[named]
  }

  # Every column below is a property of the visits/dropout pair the cell sits
  # on rather than of the calculation, so each is filled in one indexing
  # operation off the cell's position on those two axes -- not a row at a time
  # in grid_evaluate()'s loop, where a sensitivity axis would have it looked
  # up again at every level.
  di <- cells$design
  dj <- cells$dropout
  design_of <- (di - 1L) * length(drop_list) + dj
  n_visits <- unname(lengths(visit_list))

  out <- list(
    design           = names(visit_list)[di],
    dropout          = names(drop_list)[dj],
    scheduled_visits = n_visits[di],
    last_visit       = vapply(visit_list, function(v) v[length(v)], numeric(1L),
                              USE.NAMES = FALSE)[di],
    dropout_total    = vapply(designs, function(d) sum(d$dropout), numeric(1L),
                              USE.NAMES = FALSE)[design_of]
  )

  # The scalar axes' values for one cell, as an index into each axis rather
  # than the value itself, so grid_evaluate() can refill one shared args list
  # by position instead of reassembling and renaming it n_cells times.
  scalar_idx <- lapply(names(scalar_lists), function(nm) cells[[nm]])

  list(visit_list = visit_list, drop_list = drop_list, scalar_lists = scalar_lists,
      designs = designs, design_of = design_of, n_cells = n_cells,
      axes = axes, named = named, labels_at = labels_at, scalar_idx = scalar_idx,
      out = out)
}

#' The expected number of visits one participant attends under a design
#'
#' Stratum `j` of the Dawson-Lagakos pattern mixture (see [trial_design()])
#' attends `visits[1:j]` and nothing after -- `j` visits -- and the
#' completers, a proportion `1 - sum(dropout)`, attend all `length(visits)`.
#' The expectation over strata is therefore a dropout-weighted mean of visit
#' counts, which collapses to `length(visits)` itself when `dropout` is all
#' zero -- every participant is a completer.
#' @param design A `trial_design`.
#' @return A single number, generally not a whole one: an expectation, not a
#'   count any one participant can attend.
#' @noRd
expected_visits <- function(design) {
  sum(design$dropout * seq_along(design$dropout)) +
    (1 - sum(design$dropout)) * length(design$visits)
}

#' Walk the cross product of a normalised grid, evaluating one design per cell
#'
#' `evaluate(design, args)` takes a `trial_design` and the cell's values for the
#' scalar axes, and returns a `slope_sample_size` or `slope_power` object.
#' Everything except that call is shared between the two grids, so the two
#' tables are guaranteed to have the same shape and the same dropout expansion.
#'
#' `g` is what [grid_axes()] built. Kept as a second function rather than
#' folded back into one, so that a caller wanting only the axes -- a bootstrap
#' grid solving each cell several hundred times over its own resampled
#' parameters, rather than once -- can build `g` and read `g$designs` directly
#' instead of paying for this plain solve first and discarding it.
#' @noRd
grid_evaluate <- function(g, evaluate, context) {
  scalar_lists <- g$scalar_lists
  n_cells <- g$n_cells
  tte_direction <- character(0L)

  # Eight columns, filled in place and returned as a list, made a data frame
  # once by the caller. Assembling a one-row data frame per cell and
  # rbind()ing the list cost more than the calculation the grid exists to
  # perform -- about 2 ms a cell against 0.6 ms for the sample size itself, so
  # roughly 60% of a grid's total runtime went on building the table rather
  # than filling it.
  res_cols <- list(
    n             = numeric(n_cells),
    n_per_arm     = numeric(n_cells),
    power         = numeric(n_cells),
    alpha         = numeric(n_cells),
    effectiveness = numeric(n_cells),
    tte           = numeric(n_cells),
    var_tte       = numeric(n_cells),
    effect_size   = numeric(n_cells)
  )

  # The scalar axes' values for one cell, spliced into the stage-two call. The
  # names are the same in every cell, so the list is built once and refilled by
  # position rather than reassembled and renamed n_cells times.
  args <- lapply(scalar_lists, function(levels) levels[[1L]])

  for (k in seq_len(n_cells)) {
    for (j in seq_along(args)) args[[j]] <- scalar_lists[[j]][[g$scalar_idx[[j]][k]]]

    # effect_components() warns, by the same class mechanism the baseline
    # warning in grid_axes() uses, when the target treatment effect makes the
    # slope more extreme rather than less. Collected and reported once here
    # too, rather than once per cell.
    res <- tryCatch(
      withCallingHandlers(
        evaluate(g$designs[[g$design_of[k]]], args),
        slopepower_tte_direction = function(w) {
          tte_direction <<- c(tte_direction, paste(g$labels_at(k), collapse = " / "))
          invokeRestart("muffleWarning")
        }
      ),
      error = grid_cell_error(context, g$labels_at(k))
    )

    # The eight result columns are read straight off `res`. That is the same
    # thing as.data.frame.slope_result() reports -- it copies each of these
    # eight from `x` unchanged -- and routing through it instead was the single
    # most expensive line in the loop, at 1.3 ms a cell. The guarantee that
    # the two cannot drift is kept by test-grid.R, which pins a grid row
    # against as.data.frame() of the same object, rather than by paying for
    # the other ten columns and discarding them.
    res_cols$n[k]             <- res$n
    res_cols$n_per_arm[k]     <- res$n_per_arm
    res_cols$power[k]         <- res$power
    res_cols$alpha[k]         <- res$alpha
    res_cols$effectiveness[k] <- res$effectiveness
    res_cols$tte[k]           <- res$tte
    res_cols$var_tte[k]       <- res$var_tte
    res_cols$effect_size[k]   <- res$effect_size
  }

  report_collected(context, tte_direction, n_cells,
                   paste0("the target treatment effect makes the slope more extreme (%s). ",
                          "The comparator slope is further from zero than the group being ",
                          "treated; check that `slope_comparator` is the intended target."))

  # The anticipated visit burden: `n`/`n_per_arm` times the design's own
  # expected_visits(), one vectorised pass off each cell's design rather than
  # a row at a time in the loop above -- the same reason dropout_total is
  # computed that way in grid_axes().
  visits_pp <- vapply(g$designs, expected_visits, numeric(1L))[g$design_of]
  res_cols$visits_total   <- res_cols$n * visits_pp
  res_cols$visits_per_arm <- res_cols$n_per_arm * visits_pp

  res_cols
}

#' Build and solve a grid: [grid_axes()] then [grid_evaluate()]
#' @noRd
grid_impl <- function(visits, dropout, scalars, evaluate, context) {
  g <- grid_axes(visits, dropout, scalars, context)
  res <- grid_evaluate(g, evaluate, context)
  as.data.frame(c(g$out, res), stringsAsFactors = FALSE)
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
#' @param n Total number of participants. Required. A single value is held
#'   constant across the grid; several make it another axis.
#' @param effectiveness,alpha Passed to [slope_power()]. A single value is held
#'   constant across the grid; a numeric vector or named list of values makes
#'   that argument another axis of the grid, so that every design is priced at
#'   every level of it. Names are used to identify the cell in any message; the
#'   value itself is what the table reports, in the column of the same name.
#' @param target Passed to [slope_power()] and held constant across the grid.
#'
#' @return A data frame with one row per combination of the axes supplied, with
#'   columns `design`, `dropout`, `scheduled_visits`, `last_visit`,
#'   `dropout_total`, `n`, `n_per_arm`, `power`, `alpha`, `effectiveness`,
#'   `tte`, `var_tte`, `effect_size`, `visits_total` and `visits_per_arm`.
#'   `power` is what varies; `n` is constant unless it was itself given several
#'   values. `visits_total`/`visits_per_arm` are `n`/`n_per_arm` times the
#'   expected number of visits one participant attends -- `scheduled_visits`
#'   itself when `dropout_total` is zero, less otherwise, since a participant
#'   who withdraws still contributes the visits attended before doing so (see
#'   [trial_design()]'s "How dropout enters the calculation"). Cells are
#'   listed with `design` varying slowest, then `dropout`, then any further
#'   axes in the order the arguments are declared above.
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
#' # Sensitivity analysis: any of `n`, `effectiveness` and `alpha` may take
#' # several values, and each becomes another axis of the same table. Here the
#' # annual schedule is priced at three assumed effectivenesses and two sample
#' # sizes, for six rows in all.
#' slope_power_grid(
#'   pars, visits = 0:3, dropout = dropout_rate(0.05),
#'   n = c(`450` = 450, `600` = 600),
#'   effectiveness = list(`20pc` = 0.2, `25pc` = 0.25, `33pc` = 0.33)
#' )
#'
#' # A bare vector is as good a spelling of an axis as a named list. With one
#' # design and two axes the long table reads better reshaped, which is what
#' # the axes reporting their own value in a column of their own name is for:
#' # here the assumed effectiveness goes down the side and the significance
#' # level across the top.
#' sens <- slope_power_grid(
#'   pars, visits = 0:3, n = 450,
#'   effectiveness = c(0.20, 0.25, 0.33), alpha = c(0.05, 0.01)
#' )
#' round(xtabs(power ~ effectiveness + alpha, sens), 3)
#'
#' # Designs and assumptions vary together in one call: three schedules at two
#' # assumed effectivenesses is six rows, and the schedules can still be
#' # compared within each level because every one of them was priced at it.
#' slope_power_grid(
#'   pars, n = 450, effectiveness = c(0.25, 0.33),
#'   visits  = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
#'   dropout = dropout_rate(0.05)
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
#' @seealso [slope_power()], [slope_sample_size_grid()], [dropout_rate()],
#'   [slope_sample_size_grid_boot()] for a bootstrapped interval around the
#'   converse table -- the sample size each design needs
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

  grid_stage_two(params, visits, dropout, "n", n, effectiveness, target, alpha,
                 slope_power, context)
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
#' @param power Target power. Defaults to 0.8. A single value is held constant
#'   across the grid; several make it another axis.
#' @param effectiveness,alpha Passed to [slope_sample_size()]. A single value is
#'   held constant across the grid; a numeric vector or named list of values
#'   makes that argument another axis of the grid, so that every design is sized
#'   at every level of it. Names are used to identify the cell in any message;
#'   the value itself is what the table reports, in the column of the same name.
#' @param target Passed to [slope_sample_size()] and held constant across the
#'   grid.
#'
#' @return A data frame with the same columns as [slope_power_grid()]. `n` and
#'   `n_per_arm` are what vary; `power` is constant unless it was itself given
#'   several values.
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
#' # Sensitivity analysis: `power`, `effectiveness` and `alpha` may each take
#' # several values, and each becomes another axis of the same table. Here two
#' # dropout assumptions are sized at two target powers and two assumed
#' # effectivenesses, for eight rows in all.
#' slope_sample_size_grid(
#'   pars, visits = 0:3,
#'   dropout = list(`0pc` = dropout_rate(0), `5pc` = dropout_rate(0.05)),
#'   power = c(`80pc` = 0.8, `90pc` = 0.9),
#'   effectiveness = list(`20pc` = 0.2, `30pc` = 0.3)
#' )
#'
#' # A bare vector is as good a spelling of an axis as a named list, and with a
#' # single design the long table reshapes into the shape such a sweep is
#' # usually reported in -- the assumption being varied down the side, the
#' # target power across the top. The sample size scales with the inverse
#' # square of the effect targeted, so the effectiveness axis is by far the
#' # more expensive of the two to be wrong about.
#' sens <- slope_sample_size_grid(
#'   pars, visits = 0:3, effectiveness = c(0.20, 0.25, 0.33),
#'   power = c(0.8, 0.9)
#' )
#' xtabs(n ~ effectiveness + power, sens)
#'
#' # Varying alpha as well: a sweep over all three assumptions at once, for one
#' # design, is twelve rows and one call.
#' slope_sample_size_grid(
#'   pars, visits = 0:3, effectiveness = c(0.25, 0.33),
#'   power = c(0.8, 0.9), alpha = c(0.05, 0.01, 0.001)
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
#' @seealso [slope_sample_size()], [slope_power_grid()], [dropout_rate()],
#'   [slope_sample_size_grid_boot()] for a bootstrapped confidence interval
#'   around every cell of this table, sharing one set of replicates across
#'   the whole grid
#' @export
slope_sample_size_grid <- function(params, visits, dropout = NULL, power = 0.8,
                                   effectiveness = 0.25,
                                   target = c("effectiveness", "observed"),
                                   alpha = 0.05) {
  context <- "slope_sample_size_grid()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)

  grid_stage_two(params, visits, dropout, "power", power, effectiveness, target, alpha,
                 slope_sample_size, context)
}

#' Shared body of the two grid functions
#'
#' [slope_power_grid()] and [slope_sample_size_grid()] differ only in which
#' argument holds the value they hold fixed across the grid (`n` vs `power`)
#' and which stage-two function is re-solved per cell; everything else --
#' deciding whether `effectiveness` is an axis, and the call to `grid_impl()`
#' -- is identical, so both call through here. The same shape as
#' `bootstrap_stage_two()` in bootstrap.R, for the same reason.
#'
#' `fixed_name`/`fixed_value` rather than a pre-built one-element list: the
#' name comes from a literal at each call site, so `stats::setNames()` here
#' keeps that pairing in one place rather than repeating `list(n = n)` and
#' `list(power = power)` beside two otherwise-identical blocks.
#' @noRd
grid_stage_two <- function(params, visits, dropout, fixed_name, fixed_value,
                           effectiveness, target, alpha, fn, context) {
  # maybe_add_effectiveness() decides whether `effectiveness` belongs in the
  # call at all -- it must not be passed under target = "observed" -- which here
  # is the same question as whether it is an axis of the grid.
  scalars <- maybe_add_effectiveness(stats::setNames(list(fixed_value), fixed_name),
                                     effectiveness, target)
  scalars$alpha <- alpha

  grid_impl(visits, dropout, scalars,
            function(des, args) do.call(fn,
                                        c(list(params = params, design = des, target = target),
                                          args)),
            context)
}
