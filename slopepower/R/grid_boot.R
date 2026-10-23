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
      lapply(seq_along(g$scalar_lists), function(j) g$scalar_at(j, k)),
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
#'   closure>)`. Flattened as `c(main, tte)` for [boot_replicate_matrix()] and
#'   [jackknife_values()]: column `k` is cell `k`'s sample size or power, and
#'   column `n_cells + tte_of[k]` its target treatment effect -- the one
#'   indexing rule the rest of [grid_boot_impl()] relies on.
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
  # Each level as effect_components() and target_components() take it: NULL
  # for the `target = "observed"` grid's absent level.
  eff_arg <- lapply(eff_levels, function(e) if (is.na(e)) NULL else e)

  # The unscaled effect size depends on the design alone -- `effectiveness`
  # scales it afterwards, and enters effect_components() only through `tte` --
  # so every cell of a design shares one solve per replicate, whatever its
  # `effectiveness`, `power` or `alpha`. Each design's solve is given the
  # `effectiveness` of its first cell only because effect_components() needs
  # one to validate.
  effect <- lapply(seq_along(g$designs), function(d) {
    memo_effect_size(g$designs[[d]], eff_arg[[eff_of[match(d, g$design_of)]]], target,
                     context)
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
    if (on_lattice(statistic)) {
      power_k <- cell_args[[k]]$power
      function(p) 2 * size_per_arm(scale_effect(effect_k(p), scale_k), z_a, power_k)$n_per_arm
    } else {
      n_k <- cell_args[[k]]$n_total
      function(p) power_at_n(scale_effect(effect_k(p), scale_k), z_a, n_k)$power
    }
  })

  tte <- lapply(eff_arg, function(eff_e) {
    force(eff_e)
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
#' `identical()` on the parameters is the key, so a hit is exact rather than
#' heuristic. Only a success is
#' remembered: a replicate whose solve fails fails again for every cell
#' asking, each one recorded as its own `NA`, exactly as before. The single
#' replicate held between calls is the one in hand anyway.
#'
#' The value is `effect_size` unscaled; each cell applies its own
#' `effectiveness` with [scale_effect()], as [solve_slope()] (power.R) does.
#' @noRd
memo_effect_size <- function(design, effectiveness, target, context) {
  last_key <- NULL
  last <- NULL
  function(p) {
    # Keyed on the parameters without the fitted model, so the closure does not
    # hold a replicate's whole `lme` fit -- several hundred kB per design --
    # between calls. Everything effect_components() reads is still in the key.
    key <- p
    key$fit <- key$call <- NULL
    if (!identical(key, last_key)) {
      last <<- effect_components(p, design, target, effectiveness, context)$effect_size
      last_key <<- key
    }
    last
  }
}

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
#' @param n_total,n_per_arm For `slope_power_grid_boot()`: the sample size held
#'   fixed, as in [slope_power_grid()]. Supply exactly one. A single value is
#'   held constant across the grid; several make it another axis.
#' @param per_arm Which basis the printed table's counts are on: `TRUE` (the
#'   default) for participants per arm, `FALSE` for the trial total. Display
#'   only, as for the plain grids: both bases are always columns, and a grid
#'   can be printed either way with `print(x, per_arm = ...)`.
#'
#' @return A data frame of class `c("slope_sample_size_grid_boot",
#'   "slope_sample_size_grid", "data.frame")` (or the power equivalents), one
#'   row per cell, with every column [slope_sample_size_grid()] reports, under
#'   the same names. A bootstrapped quantity `q` adds `q_mean`, `q_sd`,
#'   `q_lower` and `q_upper` -- the same names a single-design bootstrap
#'   result uses for its fields (see [slope_sample_size_boot()]):
#'   \describe{
#'     \item{`n_total_*`, `n_per_arm_*`}{For `slope_sample_size_grid_boot()`:
#'       the bootstrap mean, SD and confidence interval of that cell's sample
#'       size, on both bases, the per-arm columns exactly half the totals. The
#'       interval is widened to the nearest even trial sizes, as
#'       [slope_sample_size_boot()] does.}
#'     \item{`power_*`}{For `slope_power_grid_boot()`: the same, for the power
#'       that cell's fixed sample size achieves.}
#'     \item{`tte_*`}{The same, for the target treatment effect behind that
#'       cell. Constant down the `design` and `dropout` axes -- it does not
#'       depend on the design, only on `effectiveness` -- so it repeats unless
#'       `effectiveness` is itself an axis of the grid, and cells sharing an
#'       `effectiveness` share one interval exactly.}
#'     \item{`n_ci_method` (or `power_ci_method`), `tte_ci_method`}{`"bca"` or
#'       `"percentile"`: the method actually used for that cell's interval,
#'       which can fall back to percentile cell by cell -- and, since the two
#'       intervals have their own bias corrections and their own jackknife
#'       columns, one without the other -- even when `ci_method = "bca"` was
#'       requested. `NA` for an interval that could not be built at all.}
#'     \item{`replicates_failed`}{How many of the `R` replicates failed to
#'       yield the solved-for quantity for that cell -- refits that failed, plus
#'       any resample whose stage-two solve itself failed.}
#'   }
#'
#'   The table also carries, as attributes rather than columns because they
#'   are shared by every cell: `R`, `requested_ci_method`, `statistic` (`"n"`
#'   or `"power"`), `level`, `slope_se`, `refits_failed` (replicates whose
#'   stage-one refit failed, before any per-cell solve is attempted),
#'   `straddle`, `straddle_of` (as for [slope_sample_size_boot()]), `named`
#'   (which axes vary), `per_arm`, and the slope block a single-design
#'   bootstrap carries -- `slope`, `slope_replicates`, `slope_mean`,
#'   `slope_sd`, `slope_lower`, `slope_upper`, `slope_ci_method` -- since the
#'   slope is what the resampling perturbs and every cell's interval is a
#'   function of it. They describe the whole table, so a subset such as
#'   `x[1:3, ]` drops them *and* the class, returning a plain data frame of
#'   the surviving cells. [confint()][confint.slope_sample_size_grid_boot]
#'   returns one quantity's interval in every cell, as a matrix.
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
#'   pars, n_total = 450, effectiveness = 0.33,
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
                                  n_total = NULL, n_per_arm = NULL,
                                  effectiveness = 0.25,
                                  target = c("effectiveness", "observed"),
                                  alpha = 0.05, per_arm = TRUE,
                                  R = 999, ci_method = c("bca", "percentile"),
                                  level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_power_grid_boot()"
  n_total <- resolve_n(n_total, n_per_arm, paste0(
    "this grid holds the sample size fixed and bootstraps\n",
    "  the power each design achieves. For the sample size each design needs,\n",
    "  use slope_sample_size_grid_boot()."), context, axis = TRUE)
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  grid_boot_impl(params, visits, dropout, match.arg(dropout_scale), "n_total", n_total,
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
  # function's -- and an axis added there reaches this grid too. grid_axes()
  # and grid_evaluate() are called here directly, since `g` is needed on its
  # own to price each cell several hundred times over.
  spec <- grid_stage_two_spec(params, fixed_name, fixed_value, effectiveness, target, alpha,
                              context)
  g <- grid_axes(visits, dropout, dropout_scale, spec$scalars, context)
  pts <- grid_evaluate(g, spec$evaluate, grid_result_cols(fixed_name), context)

  cc <- grid_boot_computes(g, target, statistic, context)
  computes <- c(cc$main, cc$tte)

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
                                         pts[[stat_field(statistic)]][k], ci_method, probs,
                                         context,
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
  # rows are the two endpoint columns.
  main_ci <- extract(main_res, "ci", numeric(2L))
  tte_ci <- extract(tte_res, "ci", numeric(2L))
  main_mean <- extract(main_res, "mean", numeric(1L))
  main_sd <- extract(main_res, "sd", numeric(1L))

  # The same field names a single-design bootstrap uses (quantity_block(),
  # bootstrap.R), as columns: a sample size on both bases, the per-arm block
  # exactly half the total -- the interval endpoints are already even totals,
  # and a mean and SD are linear.
  added <- if (on_lattice(statistic)) {
    c(cell_block("n_total", main_mean, main_sd, main_ci[1L, ], main_ci[2L, ]),
      cell_block("n_per_arm", main_mean / 2, main_sd / 2, main_ci[1L, ] / 2, main_ci[2L, ] / 2))
  } else {
    cell_block(statistic, main_mean, main_sd, main_ci[1L, ], main_ci[2L, ])
  }
  added[[paste0(statistic, "_ci_method")]] <- extract(main_res, "ci_method", character(1L))
  added <- c(added,
             cell_block("tte", extract(tte_res, "mean", numeric(1L)),
                        extract(tte_res, "sd", numeric(1L)), tte_ci[1L, ], tte_ci[2L, ]),
             list(tte_ci_method = extract(tte_res, "ci_method", character(1L)),
                  replicates_failed = extract(main_res, "n_failed", integer(1L))))

  df <- as.data.frame(c(g$out, pts, added), stringsAsFactors = FALSE)

  # What every cell shares is carried as attributes, not columns (CONTRACT.md
  # section 4.4), and stripped by `[` with the class.
  do.call(structure,
          c(list(df, class = c(cls, grid_class(fixed_name), "data.frame"),
                 R = R, requested_ci_method = ci_method, statistic = statistic,
                 level = level, slope_se = setup$check$se, refits_failed = n_refit_failed,
                 named = g$named, per_arm = per_arm),
            slope_replicate_summary(params$slope, good_slopes, slope_int, setup$check,
                                    mat$checks[!failed_refit])))
}

#' One bootstrapped quantity's four summary columns, named after it
#'
#' The grid's counterpart of [quantity_block()] (bootstrap.R): the estimate is
#' already a column of the plain grid, and the replicates are not kept per
#' cell.
#' @noRd
cell_block <- function(q, mean, sd, lower, upper) {
  stats::setNames(list(mean, sd, lower, upper),
                  paste0(q, c("_mean", "_sd", "_lower", "_upper")))
}

# ---------------------------------------------------------------------------
# printing
# ---------------------------------------------------------------------------

#' The target-effect frame: one row per `effectiveness` level
#'
#' `tte` depends on `effectiveness` alone, so its interval is the same down
#' every design; printed apart, once per level, when that axis varies.
#' @noRd
grid_boot_tte_frame <- function(x, show_method) {
  first <- !duplicated(x$effectiveness)
  cols <- c("effectiveness", "tte", "tte_mean", "tte_sd", "tte_lower", "tte_upper",
            if (show_method) "tte_ci_method")
  out <- x[first, cols, drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @rdname print.slope_sample_size_grid
#' @export
print.slope_sample_size_grid_boot <- function(x, ..., per_arm = NULL) {
  print_grid_boot(x, per_arm, "Bootstrapped sample size grid",
                  "print.slope_sample_size_grid_boot()", ...)
}

#' @rdname print.slope_sample_size_grid
#' @export
print.slope_power_grid_boot <- function(x, ..., per_arm = NULL) {
  print_grid_boot(x, per_arm, "Bootstrapped power grid", "print.slope_power_grid_boot()", ...)
}

# The table is the grid on one display basis (basis_columns(), grid.R), every
# column under its stored name. The target-effect columns move to a frame of
# their own, and a `<stat>_ci_method` column is shown only where some cell fell
# back from the method asked for -- otherwise it would repeat the method the
# notes already state.
print_grid_boot <- function(x, per_arm, kind, context, ...) {
  requested <- attr(x, "requested_ci_method")
  statistic <- attr(x, "statistic")
  R <- attr(x, "R")
  refits_failed <- attr(x, "refits_failed")
  named <- attr(x, "named")
  per_arm <- display_basis(x, per_arm, context)

  plain <- strip_grid(x)
  method_col <- paste0(statistic, "_ci_method")
  mixed <- function(m) any(!is.na(m) & m != requested)
  hide <- c(grep("^tte_", names(plain), value = TRUE), "replicates_failed",
            if (!mixed(plain[[method_col]])) method_col)
  cells <- basis_columns(plain[setdiff(names(plain), hide)], per_arm)

  # Single-bracket indexing, not `[[`: `named` carries no `effectiveness`
  # element at all for a grid solved with `target = "observed"`, where it is
  # not an axis, and `[[` on an absent name errors where `[` returns NA.
  tte <- if (isTRUE(unname(named["effectiveness"]))) {
    grid_boot_tte_frame(plain, mixed(plain$tte_ci_method))
  }

  cat_title(grid_title(x, kind))
  print.data.frame(cells, ...)
  cat("\n")
  if (!is.null(tte)) {
    print.data.frame(tte, ...)
    cat("\n")
  }

  cat(boot_method_note(requested, R, level = attr(x, "level")), sep = "\n")
  cat(boot_refit_note(refits_failed, R, ", and were discarded from every cell"), sep = "\n")

  # `replicates_failed` counts every replicate that yielded no solved-for value
  # for the cell, refit failures included; the refits are reported on their
  # own above, so what is left to report here is the difference -- the losses
  # that were the stage-two solve's own.
  extra_failed <- x$replicates_failed - refits_failed
  if (any(extra_failed > 0L)) {
    cat(boot_note("Note", sprintf(paste0(
      "cells also lost replicates solving for `%s` itself, beyond the refits above ",
      "(%d-%d of %d failed per cell)."), statistic, min(extra_failed), max(extra_failed),
      R)), sep = "\n")
  }
  lowers <- c(cells[[grep("_lower$", names(cells), value = TRUE)[1L]]], tte$tte_lower)
  if (anyNA(lowers)) {
    cat(boot_note("Note", paste(
      "an interval shown as NA had fewer than two surviving replicates in that cell;",
      "no interval could be built there.")), sep = "\n")
  }
  if (method_col %in% names(cells) || "tte_ci_method" %in% names(tte)) {
    cat(boot_note("Note", paste(
      "the `_ci_method` column marks the cells where BCa could not be built and a",
      "percentile interval was used.")), sep = "\n")
  }

  cat(boot_straddle_note(attr(x, "straddle"), length(attr(x, "slope_replicates")),
                         attr(x, "straddle_of")),
      sep = "\n")

  if (on_lattice(statistic)) cat(boot_mean_note(per_arm), sep = "\n")

  invisible(x)
}

#' Confidence intervals from a bootstrapped grid
#'
#' One quantity's interval in every cell, as a matrix in the shape
#' [stats::confint()] returns: a row per cell, named by the values of the axes
#' that vary, and a column for each limit. Read from the grid's
#' `<quantity>_lower` and `<quantity>_upper` columns, not recomputed.
#'
#' @param object A grid from [slope_sample_size_grid_boot()] or
#'   [slope_power_grid_boot()].
#' @param parm The one quantity whose intervals to return: `"n_per_arm"`,
#'   `"n_total"` or `"power"` (whichever the grid bootstrapped), or `"tte"`.
#'   Defaults to the bootstrapped statistic, a sample size on the display
#'   basis.
#' @param level The level the bootstrap was run at, which is the default.
#' @param ... Ignored.
#' @param per_arm For a sample-size grid, which basis the default `parm` is on.
#'   Defaults to the grid's own.
#' @return A matrix with one row per cell and columns for the lower and upper
#'   limits, labelled by percentage.
#' @export
confint.slope_sample_size_grid_boot <- function(object, parm, level = attr(object, "level"),
                                                ..., per_arm = NULL) {
  context <- "confint()"
  statistic <- attr(object, "statistic")
  if (missing(parm)) {
    parm <- main_quantity(statistic, on_lattice(statistic) && display_basis(object, per_arm, context))
  }
  available <- c(stat_quantities(statistic), "tte")
  if (!is.character(parm) || length(parm) != 1L || !parm %in% available) {
    stop(sprintf("%s: `parm` must be one of %s.", context,
                 paste(sQuote(available), collapse = ", ")), call. = FALSE)
  }
  check_confint_level(level, attr(object, "level"), context)
  out <- cbind(object[[paste0(parm, "_lower")]], object[[paste0(parm, "_upper")]])
  dimnames(out) <- list(grid_row_labels(object), confint_labels(attr(object, "level")))
  out
}

#' @rdname confint.slope_sample_size_grid_boot
#' @export
confint.slope_power_grid_boot <- confint.slope_sample_size_grid_boot

#' A label per grid row, from the values of the axes that vary
#'
#' The design and dropout names always, then each scalar axis with more than
#' one level -- a fixed sample size read from its requested total.
#' @noRd
grid_row_labels <- function(x) {
  named <- attr(x, "named")
  axes <- names(named)[named]
  cols <- lapply(axes, function(nm) {
    v <- x[[if (identical(nm, "n_total")) "n_requested_total" else nm]]
    if (is.numeric(v)) fmt_num(v) else as.character(v)
  })
  do.call(paste, c(cols, sep = " / "))
}
