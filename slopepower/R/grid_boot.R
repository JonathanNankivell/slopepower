# Layer 4 -- a bootstrap interval around every cell of a sample-size or power
# grid.
#
# slope_sample_size_grid() and slope_power_grid() (grid.R) price many candidate
# designs against one stage-one fit; slope_sample_size_boot() and
# slope_power_boot() (bootstrap.R) put a confidence interval around one such
# price by resampling subjects and refitting. Section 2.6 of Nash et al. (2021)
# recommends exactly this interval when a design is being chosen -- which is
# exactly when a grid is in use -- so the two belong together.
#
# Calling a single-design bootstrap once per cell would refit the stage-one model
# n_cells times as often as necessary: the resampling scheme -- boot_setup(),
# boot_replicate_matrix(), jackknife_values(), all in bootstrap.R -- depends
# only on `params`, never on the design being priced. So one set of
# replicates, refit once, prices every cell in the table; only the 0.6ms
# stage-two solve differs between cells. Sharing the replicates also makes
# the comparison between cells paired rather than independently noisy: two
# designs bootstrapped from the same draws differ only because of the design,
# not because of two different resamples.

# ---------------------------------------------------------------------------
# per-cell closures
# ---------------------------------------------------------------------------

#' The scalar axes' values for every cell of a grid, one list per cell
#'
#' grid_evaluate() (grid.R) refills one shared `args` list by position as it
#' walks the cells -- the cheaper way to fill n_cells rows of a data frame.
#' This grid instead needs each cell's own values kept around after the walk,
#' to close over inside a per-cell replicate function. `n_cells` short lists
#' cost nothing next to the model refits such a grid spends its time on.
#' @noRd
grid_cell_args <- function(g) {
  lapply(seq_len(g$n_cells), function(k) {
    stats::setNames(
      lapply(seq_along(g$scalar_lists), function(j) g$scalar_lists[[j]][[g$scalar_idx[[j]][k]]]),
      names(g$scalar_lists))
  })
}

#' The replicate statistics a grid needs: one solved-for quantity per cell, one
#' target treatment effect per distinct `effectiveness`
#'
#' Every cell of the grid gets its own `main` closure, a `function(p)
#' numeric(1)` giving the sample size [slope_sample_size()] would report
#' (`statistic = "n"`), or the power [slope_power()] would (`statistic =
#' "power"`), against one resampled replicate's refitted parameters `p` and the
#' cell's own design and scalar values. It is those functions' own arithmetic
#' -- [effect_components()] then [scale_effect()] and either [size_per_arm()]
#' or the normal-tail power, as [solve_slope()] (power.R) runs them -- but the
#' effect size is shared, through [memo_effect_size()], by
#' every cell with the same design: re-solving it once per `power` (or `n`), `alpha`
#' and `effectiveness` level as well made a D x P x A x E grid do D x P x A x E
#' full solves per replicate, and per jackknife refit, where D would do.
#'
#' `tte` gets one closure per *distinct* `effectiveness` level instead of one
#' per cell. The target treatment effect is
#' `-effectiveness * (slope - reference_slope)`: it depends on the replicate
#' and on `effectiveness`, and on nothing else this grid varies -- not the
#' visit schedule, the dropout pattern, `power` or `alpha`. One closure per
#' cell would make a D x P x E grid recompute the same E values D x P times
#' over on every replicate, and again on every jackknife refit, and then
#' store D x P identical copies of each in the replicate and jackknife
#' matrices.
#'
#' It is read from [target_components()] -- the function [slope_sample_size()]
#' itself gets it from, so the two cannot report different effects -- rather
#' than off a second [slope_sample_size()] solve, which would price a whole
#' design to reach a number that does not depend on one. That also keeps a
#' replicate's `tte` free of a cell's design: a resample whose stage-two solve
#' fails for one schedule still has the same target effect every other
#' schedule has, which is what this grid's `@return` block promises of the
#' `tte_*` columns.
#'
#' Both kinds keep the flat `function(p) numeric(1)` shape
#' [boot_replicate_matrix()] and [jackknife_values()] (both bootstrap.R) take,
#' so neither has to learn that a column may serve more than one cell; the
#' `tte_of` map below is what remembers which.
#'
#' Each closure closes over its cell's design and scalar values, or over one
#' `effectiveness` level, alone -- not over `params` or the grid itself -- so a
#' several-hundred-replicate run does not keep the original fit, and its model
#' frame, reachable through every closure. `lapply()` rather than a loop, so
#' each closure captures its own values in a fresh call frame instead of the
#' last value a shared loop variable happened to hold.
#'
#' @return `list(main = <one closure per cell>, tte = <one closure per distinct
#'   `effectiveness` level>, tte_of = <for each cell, the index of its `tte`
#'   closure>)`.
#' @noRd
grid_boot_computes <- function(g, target, statistic, context) {
  cell_args <- grid_cell_args(g)

  # NA_real_ stands for the absent level of a grid solved for
  # `target = "observed"`, where `effectiveness` is not an axis and
  # `cell_args` carries none; match() pairs NA with NA, so every cell of such
  # a grid collapses onto the single closure it needs.
  eff <- vapply(cell_args, function(a) a$effectiveness %||% NA_real_, numeric(1L))
  eff_levels <- unique(eff)
  eff_of <- match(eff, eff_levels)

  # The unscaled effect size depends on the design alone -- `effectiveness`
  # scales it afterwards, and enters effect_components() only through `tte` --
  # so every cell of a design shares one solve per replicate, whatever its
  # `effectiveness`, `power` or `alpha`. Each design's solve is given the
  # `effectiveness` of its first cell only because effect_components() needs
  # one to validate.
  effect <- lapply(seq_along(g$designs), function(d) {
    e <- eff[match(d, g$design_of)]
    memo_effect_size(g$designs[[d]], if (is.na(e)) NULL else e, target, context)
  })

  # The effectiveness each level scales by. NA is the `target = "observed"`
  # grid, where target_components() fixes it at 1.
  eff_scale <- ifelse(is.na(eff_levels), 1, eff_levels)

  # alpha, power and n were validated by the point-estimate grid already solved
  # from these same cells, so z_alpha() is taken once here, not per replicate.
  main <- lapply(seq_len(g$n_cells), function(k) {
    effect_k <- effect[[g$design_of[k]]]
    scale_k <- eff_scale[eff_of[k]]
    z_a <- z_alpha(cell_args[[k]]$alpha, context)
    if (identical(statistic, "n")) {
      power_k <- cell_args[[k]]$power
      function(p) 2 * size_per_arm(scale_effect(effect_k(p), scale_k), z_a, power_k)$n_per_arm
    } else {
      n_k <- cell_args[[k]]$n
      function(p) power_at_n(scale_effect(effect_k(p), scale_k), z_a, n_k)$power
    }
  })

  tte <- lapply(eff_levels, function(e) {
    eff_e <- if (is.na(e)) NULL else e
    function(p) target_components(p, target, eff_e, context)$tte
  })

  list(main = main, tte = tte, tte_of = eff_of)
}

#' One design's effect size against a replicate, remembered for the replicate
#'
#' [boot_replicate_matrix()] and [jackknife_values()] call every column's
#' closure on the same `p` in turn, so the cells sharing a design ask for the
#' same effect size one after another. This solves it for the first and hands
#' it back to the rest: [effect_components()] -- the parameter check, the
#' design's positive-definiteness check and every dropout stratum's variance --
#' is the whole cost of a stage-two solve, and [scale_effect()] and
#' [size_per_arm()] after it are a line of arithmetic each.
#'
#' `identical()` is the key, so a hit is exact rather than heuristic; for the
#' same object it returns at once on the pointer. Only a success is
#' remembered: a replicate whose solve fails fails again for every cell
#' asking, each one recorded as its own `NA`, exactly as before. The single
#' replicate held between calls is the one in hand anyway.
#'
#' The value is `effect_size` unscaled; each cell applies its own
#' `effectiveness` with [scale_effect()], as [solve_slope()] (power.R) does.
#' @noRd
memo_effect_size <- function(design, effectiveness, target, context) {
  last_p <- NULL
  last <- NULL
  function(p) {
    if (!identical(p, last_p)) {
      last <<- effect_components(p, design, target, effectiveness, context)$effect_size
      last_p <<- p
    }
    last
  }
}

#' Flatten the compute list into the shape [boot_replicate_matrix()] and
#' [jackknife_values()] both take
#'
#' Every cell's solved-for quantity first, then one `tte` per distinct
#' `effectiveness` level: column `k` is cell `k`'s sample size or power, and
#' column `n_cells + cc$tte_of[k]` its target treatment effect -- the one
#' indexing rule the rest of [grid_boot_impl()] relies on.
#' @noRd
grid_boot_flatten <- function(cc) c(cc$main, cc$tte)

# ---------------------------------------------------------------------------
# one cell's interval
# ---------------------------------------------------------------------------

#' One cell's mean, SD, interval and failure count for one statistic
#'
#' `col` is that statistic's column of the shared replicate matrix; `jack_col`
#' is a zero-argument accessor into the shared jackknife, in the shape
#' [boot_interval()] expects. `statistic` names what is being bootstrapped and
#' is passed straight to [widen_to_lattice()], which is keyed on that name --
#' so which statistics live on the even lattice stays written in that one
#' function, rather than being re-decided as a boolean at each call site here.
#'
#' Fewer than two surviving replicates is not fatal here the way it is in
#' [run_bootstrap()] -- a grid that cost several minutes to resample must not
#' be discarded over one bad cell -- so a starved cell reports `NA` and is
#' named in the warning [grid_boot_impl()] raises once for the whole table,
#' via [report_collected()] (grid.R).
#' @noRd
grid_boot_cell_stat <- function(col, jack_col, observed, ci_method, probs, context, what,
                                statistic) {
  bad <- is.na(col)
  n_failed <- sum(bad)
  good <- col[!bad]
  if (length(good) < 2L) {
    return(list(mean = NA_real_, sd = NA_real_, ci = c(NA_real_, NA_real_),
               ci_method = NA_character_, n_failed = n_failed, starved = TRUE))
  }
  iv <- boot_interval(good, observed, jack_col, ci_method, probs, context, what)
  list(mean = mean(good), sd = stats::sd(good), ci = widen_to_lattice(iv$ci, statistic),
      ci_method = iv$ci_method,
      n_failed = n_failed, starved = FALSE)
}

# ---------------------------------------------------------------------------
# the driver
# ---------------------------------------------------------------------------

#' Bootstrap every cell of a sample-size or power grid, in one resampling pass
#'
#' The bootstrap counterparts of [slope_sample_size_grid()] and
#' [slope_power_grid()]: they price the same cross product of visit schedules
#' and dropout assumptions, and put a confidence interval around what each one
#' solves for --- the sample size a design needs, or the power it achieves ---
#' and around the target treatment effect behind it. Nash et al. (2021,
#' section 2.6) recommends the interval precisely because the variance
#' components behind a sample size are themselves estimated; that matters most
#' exactly when a design is being chosen from several, which is what a grid is
#' for.
#'
#' Every design shares one set of resampled replicates rather than being
#' bootstrapped independently: the resampling scheme -- which subjects are
#' drawn, and the refitted stage-one parameters that come back -- depends only
#' on `params`, never on the design being priced (see
#' [slope_sample_size_boot()]'s own resampling scheme, which this reuses
#' unchanged). So `R` replicates, refit once, price every cell in the table,
#' rather than `R` refits *per cell*: a nine-cell grid at the default
#' `R = 999` costs about what one [slope_sample_size_boot()] call does, not
#' nine times it. Sharing the replicates also makes the comparison between
#' cells paired: two designs' intervals move together with whatever the
#' resampling happened to do to the slope, so a real difference between two
#' designs is not confounded with two independently-drawn samples' worth of
#' Monte Carlo noise -- precisely the comparison a grid exists to support.
#'
#' A cell with fewer than two surviving replicates -- the solved-for quantity
#' and the target treatment effect can fail independently, since a resampled
#' slope difference can occasionally make the effect size non-finite -- is not
#' fatal to the whole grid the way it is to a single-design bootstrap. Losing
#' one cell of several does not justify discarding a run that may have taken
#' several minutes to resample; that cell's interval columns are `NA` instead,
#' and every such cell is named once in a warning. The grid only stops if not a
#' single cell could be bootstrapped at all.
#'
#' Unlike the single-design bootstraps, these take no `statistic`: every cell
#' gets an interval for both the quantity solved for and the target treatment
#' effect, the two a `statistic` would choose between. The `tte` interval
#' costs next to nothing beside the refits, since it is one column per
#' `effectiveness` level, not one per cell.
#'
#' @inheritParams slope_sample_size_grid
#' @inheritParams slope_sample_size_boot
#' @param n For `slope_power_grid_boot()`: total number of participants, as in
#'   [slope_power_grid()]. Required. A single value is held constant across the
#'   grid; several make it another axis.
#' @param per_arm Which basis the printed table's counts are on: `TRUE` (the
#'   default) for participants per arm, `FALSE` for the trial total. Unlike
#'   the plain grids' `per_arm` (which reshapes the returned columns), these
#'   grids always carry both bases -- see the `@return` block below -- and
#'   `per_arm` only sets which one the print method shows by default; the
#'   stored value can still be reached through the other columns, and a single
#'   call can be printed either way with `print(x, per_arm = ...)`.
#'
#' @return A data frame of class `c("slope_sample_size_grid_boot",
#'   "data.frame")` (or `"slope_power_grid_boot"`), one row per cell, with the
#'   columns [slope_sample_size_grid()] reports -- except that both bases are
#'   kept, whatever `per_arm` says: `n` is always the trial total, beside
#'   `n_per_arm`, and `visits` is split into `visits_total` and
#'   `visits_per_arm`. So `n` here is *not* the `n` of a plain grid called with
#'   the default `per_arm = TRUE`; compare `n_per_arm` with that. The table
#'   also has:
#'   \describe{
#'     \item{`n_mean`, `n_sd`, `n_lower`, `n_upper`}{For
#'       `slope_sample_size_grid_boot()`: the bootstrap mean, SD and confidence
#'       interval of that cell's sample size, on the trial-total basis of `n`.
#'       The interval is widened to the nearest even sizes a trial could
#'       actually be run at, as [slope_sample_size_boot()] does.}
#'     \item{`power_mean`, `power_sd`, `power_lower`, `power_upper`}{For
#'       `slope_power_grid_boot()`: the same, for the power that cell's fixed
#'       `n` achieves.}
#'     \item{`tte_mean`, `tte_sd`, `tte_lower`, `tte_upper`}{The same, for the
#'       target treatment effect behind that cell. Constant down the `design`
#'       and `dropout` axes -- it does not depend on the design, only on
#'       `effectiveness` -- so it repeats unless `effectiveness` is itself an
#'       axis of the grid, and cells sharing an `effectiveness` share one
#'       interval exactly rather than merely agreeing to rounding.}
#'     \item{`ci_method`, `tte_ci_method`}{`"bca"` or `"percentile"`: the method
#'       actually used for that cell's main and `tte` interval, which can fall
#'       back to percentile cell by cell -- and, since the two intervals have
#'       their own bias corrections and their own jackknife columns, one
#'       without the other -- even when `ci_method = "bca"` was requested. `NA`
#'       for an interval that could not be built at all.}
#'     \item{`n_failed`}{How many of the `R` replicates failed to yield the
#'       solved-for quantity for that cell -- refits that failed to converge,
#'       plus any resample whose stage-two solve itself failed -- out of `R`.}
#'   }
#'
#'   The table also carries, as attributes rather than columns because they
#'   are shared by every cell: `R`, `ci_method` (as requested), `statistic`
#'   (`"n"` or `"power"`), `level`, `se`, `n_refit_failed` (replicates whose
#'   stage-one refit failed, before any per-cell solve is attempted),
#'   `straddle`, `straddle_of` (as for [slope_sample_size_boot()]), `per_arm`
#'   (as requested, or overridden at print time), and
#'   `slope_observed`, `slope_mean`, `slope_sd`, `slope_ci`,
#'   `slope_ci_method`, `slope_replicates` -- the same summary of the refitted
#'   slopes a single-design bootstrap result carries, since the slope is what
#'   the resampling actually perturbs and every cell's interval is a function
#'   of it. The printed table does not show them --- it reports the cells and,
#'   beneath them, `straddle` --- so these six attributes are where a reader
#'   who wants the slope's own interval finds it. They describe the whole
#'   table, so a subset such as `x[1:3, ]` drops them *and* the class,
#'   returning a plain data frame of the surviving cells, printed as one.
#'
#' @examples
#' # No comparator: fitted to all two hundred participants of `slpower1`.
#' # A real run wants the default R = 999; this uses far fewer so the example
#' # runs quickly, the same trade-off slope_sample_size_boot()'s own examples
#' # make.
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#'
#' \donttest{
#' slope_sample_size_grid_boot(
#'   pars, power = 0.8, effectiveness = 0.33,
#'   visits  = list(annual = c(0, 1, 2, 3), six_month = seq(0, 3, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05)),
#'   R = 50, seed = 42)
#'
#' # The converse table: the power a fixed 450 participants achieves.
#' slope_power_grid_boot(
#'   pars, n = 450, effectiveness = 0.33,
#'   visits  = list(annual = c(0, 1, 2, 3), six_month = seq(0, 3, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05)),
#'   R = 50, seed = 42)
#' }
#'
#' @seealso [slope_sample_size_grid()] and [slope_power_grid()] for the point
#'   estimates alone, [slope_sample_size_boot()] and [slope_power_boot()] to
#'   bootstrap a single design, [dropout_rate()]
#' @export
slope_sample_size_grid_boot <- function(params, visits, dropout = NULL,
                                        dropout_scale = c("incremental", "cumulative"),
                                        power = 0.8, effectiveness = 0.25,
                                        target = c("effectiveness", "observed"),
                                        alpha = 0.05, per_arm = TRUE,
                                        R = 999, ci_method = c("bca", "percentile"),
                                        level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_sample_size_grid_boot()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  grid_boot_impl(params, visits, dropout, match.arg(dropout_scale), "power", power,
                 effectiveness, target, alpha, R, ci_method, level, seed, progress,
                 per_arm, context, "slope_sample_size_grid_boot")
}

#' @rdname slope_sample_size_grid_boot
#' @export
slope_power_grid_boot <- function(params, visits, dropout = NULL,
                                  dropout_scale = c("incremental", "cumulative"),
                                  n, effectiveness = 0.25,
                                  target = c("effectiveness", "observed"),
                                  alpha = 0.05, per_arm = TRUE,
                                  R = 999, ci_method = c("bca", "percentile"),
                                  level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_power_grid_boot()"
  # `is.null(n)` too; see the note on the same guard in slope_power().
  if (missing(n) || is.null(n)) {
    stop(sprintf(paste0(
      "%s: `n` is required -- this grid holds the sample size fixed and bootstraps\n",
      "  the power each design achieves. For the sample size each design needs,\n",
      "  use slope_sample_size_grid_boot()."), context), call. = FALSE)
  }
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  grid_boot_impl(params, visits, dropout, match.arg(dropout_scale), "n", n,
                 effectiveness, target, alpha, R, ci_method, level, seed, progress,
                 per_arm, context, "slope_power_grid_boot")
}

#' Shared body of the two bootstrapped grids
#'
#' `fixed_name` is the input held fixed -- `"power"` for the sample-size grid,
#' `"n"` for the power grid -- and decides the solved-for `statistic`, the
#' other of the two. Everything else, from the axes to the attribute list, is
#' the same for both, so that the two tables describe the resampling in the
#' same terms.
#' @noRd
grid_boot_impl <- function(params, visits, dropout, dropout_scale, fixed_name,
                           fixed_value, effectiveness, target, alpha, R, ci_method,
                           level, seed, progress, per_arm, context, cls) {
  statistic <- if (identical(fixed_name, "power")) "n" else "power"
  ci_method <- check_boot_args(R, ci_method, level, context)
  per_arm <- check_per_arm(per_arm, context)
  # Reproducible without reseeding the session, exactly as run_bootstrap()
  # (bootstrap.R) arranges for a single result; see seed_bootstrap() on why the
  # on.exit() stays here rather than moving into it.
  old_seed <- seed_bootstrap(seed)
  if (!is.null(seed)) on.exit(restore_seed(old_seed), add = TRUE)

  # The point estimates: built through grid_stage_two_spec() (grid.R), the same
  # axis set and the same per-cell closure the plain grid is itself built from,
  # so this table's `n`/`power`/`tte`/... columns cannot drift from that
  # function's -- and an axis added there reaches this grid too. Only the two
  # halves underneath grid_impl() are called separately, since `g` is needed on
  # its own to price each cell several hundred times over.
  spec <- grid_stage_two_spec(params, fixed_name, fixed_value, effectiveness, target, alpha,
                              context)
  g <- grid_axes(visits, dropout, dropout_scale, spec$scalars, context)
  pts <- grid_evaluate(g, spec$evaluate, context)

  cc <- grid_boot_computes(g, target, statistic, context)
  computes <- grid_boot_flatten(cc)

  setup <- boot_setup(params, context, target)
  mat <- boot_replicate_matrix(setup, computes, R, progress, context)

  failed_refit <- is.na(mat$slopes)
  n_refit_failed <- sum(failed_refit)
  good_slopes <- mat$slopes[!failed_refit]
  if (length(good_slopes) < 2L) {
    stop(sprintf(paste0("%s: %d of %d replicates failed to refit; not enough succeeded to ",
                        "bootstrap any cell."), context, n_refit_failed, R), call. = FALSE)
  }

  probs <- boot_probs(level)

  # One jackknife pass, covering every cell's main statistic and `tte` plus the
  # slope, taken on first use and kept. lazy_jackknife() (bootstrap.R) is the
  # same memoisation run_bootstrap() uses for its own single jackknife, and owns
  # the convention that the slope accessor is appended last -- so this driver
  # reads `slope_col()` rather than counting columns to find it.
  jack <- lazy_jackknife(setup, computes)

  slope_int <- boot_interval(good_slopes, params$slope, jack$slope_col, ci_method, probs,
                             context, " for the replicate slopes")

  main_res <- vector("list", g$n_cells)
  starved <- character(0L)
  for (k in seq_len(g$n_cells)) {
    label <- cell_values(g$labels_at(k))
    main_res[[k]] <- grid_boot_cell_stat(mat$replicates[, k], function() jack$col(k),
                                         pts[[statistic]][k], ci_method, probs, context,
                                         sprintf(" for cell %s", label),
                                         statistic = statistic)
    if (isTRUE(main_res[[k]]$starved)) starved <- c(starved, label)
  }
  if (length(starved) == g$n_cells) {
    stop(sprintf(paste0("%s: every cell had fewer than 2 surviving replicates for `%s`; no ",
                        "interval could be built."), context, statistic), call. = FALSE)
  }

  # One interval per distinct `effectiveness` level rather than one per cell,
  # matching the one replicate column per level grid_boot_computes() filled:
  # every cell at a level shares that column, so a per-cell pass would rebuild
  # the identical interval -- and, at ci_method = "bca", read the identical
  # jackknife column -- once for every design and dropout in the grid. Spread
  # back over the cells by `tte_of` when the columns are assembled below, so
  # the table still reports every cell's own row.
  tte_res <- lapply(seq_along(cc$tte), function(j) {
    ti <- g$n_cells + j
    # The first cell at this level, for a label naming a real row of the grid.
    k <- match(j, cc$tte_of)
    grid_boot_cell_stat(mat$replicates[, ti], function() jack$col(ti),
                        pts$tte[k], ci_method, probs, context,
                        sprintf(" for cell %s (tte)", cell_values(g$labels_at(k))),
                        statistic = "tte")
  })[cc$tte_of]

  report_collected(context, starved, g$n_cells,
                   sprintf(paste0("fewer than two replicates succeeded for `%s` (%%s), so no ",
                                  "interval could be built for it there. Its `%s_*` columns are ",
                                  "NA for those rows; its `tte_*` columns are built separately ",
                                  "and are unaffected."), statistic, statistic))

  extract <- function(res, field, template) vapply(res, function(r) r[[field]], template)

  # `ci` is always a pair, so extract()ing it gives a 2-by-cells matrix whose
  # rows are the two endpoint columns -- the same helper as every other column
  # rather than four more hand-written closures beside it.
  main_ci <- extract(main_res, "ci", numeric(2L))
  tte_ci <- extract(tte_res, "ci", numeric(2L))

  added <- stats::setNames(
    list(extract(main_res, "mean", numeric(1L)),
         extract(main_res, "sd", numeric(1L)),
         main_ci[1L, ],
         main_ci[2L, ]),
    paste0(statistic, c("_mean", "_sd", "_lower", "_upper")))
  added <- c(added, list(
    tte_mean = extract(tte_res, "mean", numeric(1L)),
    tte_sd = extract(tte_res, "sd", numeric(1L)),
    tte_lower = tte_ci[1L, ],
    tte_upper = tte_ci[2L, ],
    tte_ci_method = extract(tte_res, "ci_method", character(1L)),
    ci_method = extract(main_res, "ci_method", character(1L)),
    n_failed = extract(main_res, "n_failed", integer(1L))
  ))

  df <- as.data.frame(c(g$out, pts, added), stringsAsFactors = FALSE)

  # `added_cols` rather than a constant listing these names a second time: the
  # printed cells frame is defined by subtracting them from the table, so a
  # column added above and forgotten in a hand-written list would leak straight
  # into that frame. Carried as an attribute, and so stripped by `[` along with
  # the rest of the grid-wide summary.
  do.call(structure,
          c(list(df, class = c(cls, "data.frame"),
                 R = R, ci_method = ci_method, statistic = statistic, level = level,
                 se = setup$check$se, n_refit_failed = n_refit_failed,
                 added_cols = names(added),
                 named = g$named, per_arm = per_arm),
            slope_replicate_summary(params$slope, good_slopes, slope_int, setup$check,
                                    mat$checks[!failed_refit])))
}

#' Subsetting drops the grid-wide summary, not just the class
#'
#' Every attribute [slope_sample_size_grid_boot()] and [slope_power_grid_boot()]
#' add -- `R`, the slope block, `named` and the rest -- describes the *whole
#' table*: `R` replicates were drawn once for every cell, `slope_replicates`
#' holds the ones that priced all `nrow(x)` of them, and `named` says which axes
#' vary across the lot. A `[` subset can change which cells are present, or how
#' many, or leave only one level of an axis standing, without changing any of
#' that -- and base `[.data.frame` does not reliably strip attributes it does
#' not recognise, so, empirically, those become stale rather than absent:
#' `slope_replicates` still summarised the resampling behind cells a row subset
#' had dropped. So this method strips them itself, on every subset, along with
#' the class -- explicitly, rather than relying on incidental attribute-dropping
#' the print method could not safely assume.
#'
#' Every attribute *but* the base data-frame ones (`names`, `row.names`,
#' `class`) is stripped, rather than naming each of the grid's own attributes
#' here too -- the second list that could drift from the first if one gained
#' an attribute the other forgot, the same hazard `stage_two_result()`
#' (power.R) exists to avoid for the result classes' field lists.
#'
#' A subsetted table is then unambiguously a plain data frame: still every
#' cell that survived the subset, just no longer carrying a summary that no
#' longer matches it.
#' @param x A `slope_sample_size_grid_boot` or `slope_power_grid_boot` object to
#'   subset.
#' @param ... Passed on to `[.data.frame`.
#' @export
`[.slope_sample_size_grid_boot` <- function(x, ...) {
  out <- NextMethod()
  if (inherits(out, "data.frame")) {
    class(out) <- setdiff(class(out), c("slope_sample_size_grid_boot",
                                        "slope_power_grid_boot"))
    extra <- setdiff(names(attributes(out)), c("names", "row.names", "class"))
    for (a in extra) attr(out, a) <- NULL
  }
  out
}

#' @rdname sub-.slope_sample_size_grid_boot
#' @export
`[.slope_power_grid_boot` <- `[.slope_sample_size_grid_boot`

# ---------------------------------------------------------------------------
# printing
# ---------------------------------------------------------------------------

#' One statistic's interval as a single printed column
#'
#' `lower`/`upper` joined by `boot_interval_col()` (bootstrap.R) into the
#' character column a data frame can hold, with the two markers the hand-drawn
#' table used to carry in its own cells:
#'
#' * `"--"` where the cell had fewer than two surviving replicates, so there is
#'   no interval to show. Not `NA`: the note below the table names `"--"`, and
#'   a reader should find in the table the thing the note pointed at.
#' * `" *"` where the interval that could be built is not the one asked for.
#'   `used` is the per-cell `ci_method` -- the main statistic's or `tte`'s,
#'   since either can fall back without the other; `asked` the grid-wide
#'   `ci_method`.
#'
#' Both required. They were optional while only the `n` column carried the
#' marker, and a `tte` column that forgot them silently claimed a method it had
#' not used -- so there is no longer a caller that should be allowed to omit them.
#' @noRd
grid_boot_ci_col <- function(lower, upper, used, asked) {
  out <- boot_interval_col(lower, upper)
  mixed <- !is.na(used) & used != asked
  out[mixed] <- paste0(out[mixed], " *")
  ifelse(is.na(lower), "--", out)
}

#' The second printed frame: the target treatment effect, one row per
#' `effectiveness` level
#'
#' `tte` depends on `effectiveness` and on nothing else this grid varies --
#' see [slope_sample_size_grid_boot()]'s `@return`, and grid_boot_computes(),
#' which computes one of them per level for exactly that reason. So this
#' frame is keyed by `effectiveness` alone: keying it by the `design` and
#' `dropout` columns the first frame leads with would print every interval
#' once per design, which is the repetition the frame is gated on
#' `effectiveness` varying to avoid in the first place. A grid over two
#' designs, two dropout patterns and two effectiveness levels wants two rows
#' here, not eight rows holding two distinct intervals.
#'
#' Built rather than printed, so [print_grid_boot()] can gate the notes beneath
#' both frames on what they actually show.
#' @noRd
grid_boot_tte_frame <- function(x, ci_method) {
  first <- !duplicated(x$effectiveness)
  data.frame(
    effectiveness = x$effectiveness[first],
    tte = x$tte[first],
    tte_mean = x$tte_mean[first],
    tte_sd = x$tte_sd[first],
    # `used`/`asked` as the first frame passes them: `tte`'s interval has its
    # own bias correction and its own jackknife column, so it can fall back to
    # percentile where the main statistic's did not, and the `*` has to say so
    # here too or the header note claims a method this column did not use.
    tte_ci = grid_boot_ci_col(x$tte_lower[first], x$tte_upper[first],
                              x$tte_ci_method[first], ci_method),
    stringsAsFactors = FALSE)
}

#' Print a bootstrapped grid
#'
#' Printed as data frames, by R itself, rather than as a hand-drawn table: a
#' grid *is* a data frame, the plain grids print as one, and this shows exactly
#' that table plus three columns for the solved-for quantity --- `n_mean`,
#' `n_sd`, and the interval as one `n_ci` column reading "924 to 1704", or
#' their `power_` counterparts. A second frame follows it for the target
#' treatment effect when `effectiveness` varies, since that is the only axis it
#' depends on. Everything shared by every cell --- the method, the replicate
#' count, the interval level, and the failure and sign-straddling counts --- is
#' reported in the notes beneath, which print on every call whether or not
#' anything went wrong.
#'
#' `n` and `visits` --- and for a sample-size grid `n_mean`, `n_sd` and `n_ci`
#' --- are shown on one basis, participants per arm by default, rather than as
#' the `n`/`n_per_arm` and `visits_total`/`visits_per_arm` pairs the underlying
#' table carries; a "Counts:" note above the other notes says which basis is on
#' the page. This is display only: the object itself is unaffected, and
#' `print(x, per_arm = !attr(x, "per_arm"))` shows the other basis of the same
#' table without rebuilding it.
#'
#' The resampled slope is not shown. It is still on the object, in the
#' `slope_observed` attribute and the five beside it, and the straddle note
#' still reports the one thing about it that bears on whether these intervals
#' mean anything.
#'
#' Falls back to `print.data.frame()` when the bootstrap attributes are
#' missing. A row-subsetted table is still every bit a data frame; it just no
#' longer carries the summary this method exists to show.
#'
#' @param x A `slope_sample_size_grid_boot` or `slope_power_grid_boot` object.
#' @param per_arm Which basis to print the table's counts on: `TRUE` for
#'   participants per arm, `FALSE` for the trial total. Defaults to `NULL`,
#'   meaning "whatever the call that built `x` requested" -- read from its
#'   `per_arm` attribute, or per arm if that is absent. Passing it explicitly
#'   prints the other basis of the same object without rebuilding it.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.slope_sample_size_grid_boot <- function(x, ..., per_arm = NULL) {
  print_grid_boot(x, per_arm, "print.slope_sample_size_grid_boot()", ...)
}

#' @rdname print.slope_sample_size_grid_boot
#' @export
print.slope_power_grid_boot <- function(x, ..., per_arm = NULL) {
  print_grid_boot(x, per_arm, "print.slope_power_grid_boot()", ...)
}

# The grid columns are recovered by dropping the `added_cols` attribute -- which
# grid_boot_impl() fills from the names of the block it actually built -- rather
# than by listing them here, through this class's own `[` method, which already
# returns a plain data frame stripped of the grid-wide summary.
#' @noRd
print_grid_boot <- function(x, per_arm, context, ...) {
  if (is.null(attr(x, "R"))) {
    print.data.frame(x, ...)
    return(invisible(x))
  }

  ci_method <- attr(x, "ci_method")
  statistic <- attr(x, "statistic") %||% "n"
  level <- attr(x, "level")
  R <- attr(x, "R")
  n_refit_failed <- attr(x, "n_refit_failed")
  named <- attr(x, "named")
  per_arm <- display_basis(x, per_arm, context)

  # `[` on these classes strips the class and every grid-wide attribute, so this
  # is already the plain data frame print.data.frame() should be handed --
  # no unclass() here, and no second place that knows which attributes exist.
  # basis_columns() (grid.R) then reduces its n/n_per_arm and
  # visits_total/visits_per_arm pairs to the one basis being printed -- the
  # same reduction the plain grids apply to their returned data, here applied
  # only to what is shown.
  cells <- basis_columns(x[, setdiff(names(x), attr(x, "added_cols")), drop = FALSE], per_arm)

  # The bootstrap summary of `n` is total-only on the object (CONTRACT.md
  # section 4.4); halving it here to match `cells$n` is exact, not
  # approximate -- widen_to_lattice() already moved n_lower/n_upper out to even
  # sizes before they were stored, and the mean and SD are linear. A power has
  # no arms and is shown as stored.
  divisor <- if (on_lattice(statistic) && per_arm) 2 else 1
  col <- function(suffix) x[[paste0(statistic, suffix)]] / divisor
  cells[[paste0(statistic, "_mean")]] <- col("_mean")
  cells[[paste0(statistic, "_sd")]] <- col("_sd")
  cells[[paste0(statistic, "_ci")]] <- grid_boot_ci_col(col("_lower"), col("_upper"),
                                                       x$ci_method, ci_method)

  # The target treatment effect does not depend on the design, only on
  # `effectiveness` -- see slope_sample_size_grid_boot()'s @return -- so a
  # second frame for it earns its place only when that axis actually varies;
  # otherwise every row would repeat the same interval.
  #
  # Single-bracket indexing, not `[[`: `named` carries no `effectiveness`
  # element at all for a grid solved with `target = "observed"`, where it is
  # not an axis, and `[[` on an absent name errors where `[` returns NA.
  # Built before anything is printed so that the notes below can be gated on
  # what the two frames actually show.
  tte <- if (isTRUE(unname(named["effectiveness"]))) grid_boot_tte_frame(x, ci_method)

  cat(sprintf("<%s>\n\n", class(x)[1L]))
  print.data.frame(cells)
  cat("\n")
  if (!is.null(tte)) {
    print.data.frame(tte)
    cat("\n")
  }

  # Every interval column on the page, so that the two markers' notes below
  # follow the markers wherever they appear rather than only where the main
  # statistic put them. `tte` is NULL when its frame was not shown, and c()
  # drops it.
  shown_ci <- c(cells[[paste0(statistic, "_ci")]], tte$tte_ci)

  cat(basis_note(per_arm), sep = "\n")
  cat(boot_method_note(ci_method, R, level), sep = "\n")

  cat(boot_note("Note", sprintf(paste0(
    "%d/%d (%.1f%%) bootstrap replicates failed to refit the stage-one model, and were ",
    "discarded from every cell."), n_refit_failed, R, 100 * n_refit_failed / R)), sep = "\n")

  # `x$n_failed` counts every replicate that yielded no solved-for value for
  # the cell, refit failures included; the refits are reported on their own
  # above, so what is left to report here is the difference -- the losses that
  # were the stage-two solve's own.
  extra_failed <- x$n_failed - n_refit_failed
  if (any(extra_failed > 0L)) {
    cat(boot_note("Note", sprintf(paste0(
      "cells also lost replicates solving for `%s` itself, beyond the refits above ",
      "(%d-%d of %d failed per cell)."), statistic, min(extra_failed), max(extra_failed),
      R)), sep = "\n")
  }
  if (any(shown_ci == "--", na.rm = TRUE)) {
    cat(boot_note("Note", paste(
      "cells marked \"--\" had fewer than two surviving replicates for that quantity;",
      "no interval could be built for them.")), sep = "\n")
  }

  cat(boot_straddle_note(attr(x, "straddle"), length(attr(x, "slope_replicates")),
                         attr(x, "straddle_of")),
      sep = "\n")

  if (on_lattice(statistic)) {
    cat(boot_note("Mean, SD", paste0(
      "each replicate of `n` is rounded up to a whole participant per arm before averaging, ",
      "so its mean is not a runnable ", if (per_arm) "arm" else "trial", " size.")), sep = "\n")
  }

  if (any(grepl("*", shown_ci, fixed = TRUE), na.rm = TRUE)) {
    cat("  * percentile interval; BCa could not be built there.\n")
  }

  invisible(x)
}
