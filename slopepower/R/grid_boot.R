# Layer 4 -- a bootstrap interval around every cell of a sample-size grid.
#
# slope_sample_size_grid() (grid.R) prices many candidate designs against one
# stage-one fit; slope_bootstrap() (bootstrap.R) puts a confidence interval
# around one such price by resampling subjects and refitting. Section 2.6 of
# Nash et al. (2021) recommends exactly this interval when a design is being
# chosen -- which is exactly when a grid is in use -- so the two belong
# together.
#
# Calling slope_bootstrap() once per cell would refit the stage-one model
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

#' The replicate statistics a grid needs: one sample size per cell, one target
#' treatment effect per distinct `effectiveness`
#'
#' Every cell of the grid gets its own `n` closure, a `function(p) numeric(1)`
#' re-solving [slope_sample_size()] against one resampled replicate's refitted
#' parameters `p` and the cell's own design and scalar values. `power` is
#' spliced in as the fixed argument the way `bootstrap_stage_two()`
#' (bootstrap.R) does for a single result; [resolve_args()] -- the same
#' helper -- supplies `design`, `target`, `alpha` and `effectiveness`.
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
#' Each closure closes over `slim`/`power_k`, or over one `effectiveness`
#' level, alone -- not over `params` or the grid itself -- so a
#' several-hundred-replicate run does not keep the original fit, and its model
#' frame, reachable through every closure. `lapply()` rather than a loop, so
#' each closure captures its own values in a fresh call frame instead of the
#' last value a shared loop variable happened to hold.
#'
#' @return `list(n = <one closure per cell>, tte = <one closure per distinct
#'   `effectiveness` level>, tte_of = <for each cell, the index of its `tte`
#'   closure>)`.
#' @noRd
grid_boot_computes <- function(g, target, context) {
  cell_args <- grid_cell_args(g)

  n <- lapply(seq_len(g$n_cells), function(k) {
    slim <- list(design = g$designs[[g$design_of[k]]], target = target,
                alpha = cell_args[[k]]$alpha, effectiveness = cell_args[[k]]$effectiveness)
    power_k <- cell_args[[k]]$power
    function(p) do.call(slope_sample_size, c(resolve_args(p, slim), list(power = power_k)))$n
  })

  # NA_real_ stands for the absent level of a grid solved for
  # `target = "observed"`, where `effectiveness` is not an axis and
  # `cell_args` carries none; match() pairs NA with NA, so every cell of such
  # a grid collapses onto the single closure it needs.
  eff <- vapply(cell_args,
                function(a) if (is.null(a$effectiveness)) NA_real_ else a$effectiveness,
                numeric(1L))
  eff_levels <- unique(eff)
  tte <- lapply(eff_levels, function(e) {
    eff_e <- if (is.na(e)) NULL else e
    function(p) target_components(p, target, eff_e, context)$tte
  })

  list(n = n, tte = tte, tte_of = match(eff, eff_levels))
}

#' Flatten the compute list into the shape [boot_replicate_matrix()] and
#' [jackknife_values()] both take
#'
#' Every cell's `n` first, then one `tte` per distinct `effectiveness` level:
#' column `k` is cell `k`'s sample size, and column `n_cells + cc$tte_of[k]`
#' its target treatment effect -- the one indexing rule the rest of
#' [slope_sample_size_grid_boot()] relies on. Two blocks rather than the `n`
#' and `tte` pairs this used to interleave, since there is no longer one `tte`
#' column per cell to pair each `n` with.
#' @noRd
grid_boot_flatten <- function(cc) c(cc$n, cc$tte)

#' A compact row label for one cell of a bootstrap grid's printed table
#'
#' `g$labels_at(k)` (grid.R) names only the axes that vary -- `design` and
#' `dropout` always, any scalar axis with more than one level besides --
#' exactly the set [grid_cell_error()]'s messages already use. Joined by " / "
#' rather than [cell_label()]'s `axis "value"` form, which is built for an
#' error message naming its axes; here the axis names are already the column
#' headers a data-frame row sits under, so only the values are needed.
#' @noRd
grid_boot_row_label <- function(labels) paste(labels, collapse = " / ")

# ---------------------------------------------------------------------------
# one cell's interval
# ---------------------------------------------------------------------------

#' One cell's mean, SD, interval and failure count for one statistic
#'
#' `col` is that statistic's column of the shared replicate matrix; `jack_col`
#' is a zero-argument accessor into the shared jackknife, in the shape
#' [boot_interval()] expects.
#'
#' Fewer than two surviving replicates is not fatal here the way it is in
#' [run_bootstrap()] -- a grid that cost several minutes to resample must not
#' be discarded over one bad cell -- so a starved cell reports `NA` and is
#' named in the warning [slope_sample_size_grid_boot()] raises once for the
#' whole table, via [report_collected()] (grid.R).
#' @noRd
grid_boot_cell_stat <- function(col, jack_col, observed, type, probs, context, what, lattice) {
  bad <- is.na(col)
  n_failed <- sum(bad)
  good <- col[!bad]
  if (length(good) < 2L) {
    return(list(mean = NA_real_, sd = NA_real_, ci = c(NA_real_, NA_real_),
               type = NA_character_, n_failed = n_failed, starved = TRUE))
  }
  iv <- boot_interval(good, observed, jack_col, type, probs, context, what)
  ci <- if (lattice) widen_to_lattice(iv$ci, "n") else iv$ci
  list(mean = mean(good), sd = stats::sd(good), ci = ci, type = iv$type,
      n_failed = n_failed, starved = FALSE)
}

# ---------------------------------------------------------------------------
# the driver
# ---------------------------------------------------------------------------

#' Bootstrap the sample size every cell of a grid needs, in one resampling pass
#'
#' The bootstrap counterpart of [slope_sample_size_grid()]: prices the same
#' cross product of visit schedules and dropout assumptions, and puts a
#' confidence interval around the sample size -- and the target treatment
#' effect behind it -- that each one needs. Nash et al. (2021, section 2.6)
#' recommends the interval precisely because the variance components behind a
#' sample size are themselves estimated; that matters most exactly when a
#' design is being chosen from several, which is what a grid is for.
#'
#' Every design shares one set of resampled replicates rather than being
#' bootstrapped independently: the resampling scheme -- which subjects are
#' drawn, and the refitted stage-one parameters that come back -- depends only
#' on `params`, never on the design being priced (see [slope_bootstrap()]'s
#' own resampling scheme, which this reuses unchanged). So `R` replicates,
#' refit once, price every cell in the table, rather than `R` refits *per
#' cell*: a nine-cell grid at the default `R = 999` costs about what one
#' [slope_bootstrap()] call does, not nine times it. Sharing the replicates
#' also makes the comparison between cells paired: two designs' intervals move
#' together with whatever the resampling happened to do to the slope, so a
#' real difference between two designs is not confounded with two
#' independently-drawn samples' worth of Monte Carlo noise -- precisely the
#' comparison a grid exists to support.
#'
#' A cell with fewer than two surviving replicates -- both the sample size and
#' the target treatment effect can fail independently, since a resampled
#' slope difference can occasionally make the effect size non-finite -- is not
#' fatal to the whole grid the way it is to a single [slope_bootstrap()] call.
#' Losing one cell of several does not justify discarding a run that may have
#' taken several minutes to resample; that cell's interval columns are `NA`
#' instead, and every such cell is named once in a warning. The grid only
#' stops if not a single cell could be bootstrapped at all.
#'
#' @inheritParams slope_sample_size_grid
#' @inheritParams slope_bootstrap
#'
#' @return A data frame of class `c("slope_sample_size_grid_boot",
#'   "data.frame")`, one row per cell, with the fifteen columns
#'   [slope_sample_size_grid()] reports plus:
#'   \describe{
#'     \item{`n_mean`, `n_sd`, `n_lower`, `n_upper`}{The bootstrap mean, SD
#'       and confidence interval of that cell's sample size. The interval is
#'       widened to the nearest even sizes a trial could actually be run at,
#'       as a single [slope_bootstrap()] call does for `statistic = "n"`.}
#'     \item{`tte_mean`, `tte_sd`, `tte_lower`, `tte_upper`}{The same, for the
#'       target treatment effect behind that sample size. Constant down the
#'       `design` and `dropout` axes -- it does not depend on the design,
#'       only on `effectiveness` -- so it repeats unless `effectiveness` is
#'       itself an axis of the grid, and cells sharing an `effectiveness`
#'       share one interval exactly rather than merely agreeing to rounding.}
#'     \item{`ci_type`, `tte_ci_type`}{`"bca"` or `"percentile"`: the method
#'       actually used for that cell's `n` and `tte` interval, which can fall
#'       back to percentile cell by cell -- and, since the two intervals have
#'       their own bias corrections and their own jackknife columns, one
#'       without the other -- even when `type = "bca"` was requested. `NA` for
#'       an interval that could not be built at all.}
#'     \item{`n_failed`}{How many of the `R` replicates failed to yield a
#'       sample size for that cell -- refits that failed to converge, plus
#'       any resample whose stage-two solve itself failed -- out of `R`.}
#'   }
#'
#'   The table also carries, as attributes rather than columns because they
#'   are shared by every cell: `R`, `type` (as requested), `level`, `se`,
#'   `n_refit_failed` (replicates whose stage-one refit failed, before any
#'   per-cell solve is attempted), `straddle`, and `slope_observed`,
#'   `slope_mean`, `slope_sd`, `slope_ci`, `slope_type`, `slope_replicates` --
#'   the same summary of the refitted slopes a single [slope_bootstrap()]
#'   result carries, since the slope is what the resampling actually
#'   perturbs and every cell's interval is a function of it. The printed
#'   table does not show them --- it reports the cells and, beneath them,
#'   `straddle` --- so these six attributes are where a reader who wants the
#'   slope's own interval finds it. Base
#'   `[.data.frame` drops attributes it does not know, so `x[1:3, ]` keeps
#'   the class but not these; [print.slope_sample_size_grid_boot()] falls
#'   back to a plain data-frame print when they are absent.
#'
#' @examples
#' # No comparator: fitted to all two hundred participants of `slpower1`.
#' # A real run wants the default R = 999; this uses far fewer so the example
#' # runs quickly, the same trade-off slope_bootstrap()'s own examples make.
#' pars <- slope_params(sdmt ~ visit | id, data = slpower1)
#'
#' \donttest{
#' slope_sample_size_grid_boot(
#'   pars, power = 0.8, effectiveness = 0.33,
#'   visits  = list(annual = c(0, 1, 2, 3), six_month = seq(0, 3, 0.5)),
#'   dropout = list(none = NULL, `5pc` = dropout_rate(0.05)),
#'   R = 50, seed = 42)
#' }
#'
#' @seealso [slope_sample_size_grid()] for the point estimates alone,
#'   [slope_bootstrap()] to bootstrap a single design, [dropout_rate()]
#' @export
slope_sample_size_grid_boot <- function(params, visits, dropout = NULL, power = 0.8,
                                        effectiveness = 0.25,
                                        target = c("effectiveness", "observed"),
                                        alpha = 0.05,
                                        R = 999, type = c("bca", "percentile"),
                                        level = 0.95, seed = NULL, progress = FALSE) {
  context <- "slope_sample_size_grid_boot()"
  target <- match.arg(target)
  check_target_effectiveness(target, !missing(effectiveness), context)
  type <- check_boot_args(R, type, level, context)
  # Reproducible without reseeding the session, exactly as run_bootstrap()
  # (bootstrap.R) arranges for a single result; see seed_bootstrap() on why the
  # on.exit() stays here rather than moving into it.
  old_seed <- seed_bootstrap(seed)
  if (!is.null(seed)) on.exit(restore_seed(old_seed), add = TRUE)

  # The point estimates: built through grid_stage_two_spec() (grid.R), the same
  # axis set and the same per-cell closure slope_sample_size_grid() itself is
  # built from, so this table's `n`/`tte`/... columns cannot drift from that
  # function's -- and an axis added there reaches this grid too. Only the two
  # halves underneath grid_impl() are called separately, since `g` is needed on
  # its own to price each cell several hundred times over.
  spec <- grid_stage_two_spec(params, "power", power, effectiveness, target, alpha,
                              slope_sample_size)
  g <- grid_axes(visits, dropout, spec$scalars, context)
  pts <- grid_evaluate(g, spec$evaluate, context)

  cc <- grid_boot_computes(g, target, context)
  computes <- grid_boot_flatten(cc)

  setup <- boot_setup(params, context)
  mat <- boot_replicate_matrix(setup, computes, R, progress, context)

  failed_refit <- is.na(mat$slopes)
  n_refit_failed <- sum(failed_refit)
  good_slopes <- mat$slopes[!failed_refit]
  if (length(good_slopes) < 2L) {
    stop(sprintf(paste0("%s: %d of %d replicates failed to refit; not enough succeeded to ",
                        "bootstrap any cell."), context, n_refit_failed, R), call. = FALSE)
  }

  probs <- boot_probs(level)

  # One jackknife pass, covering every cell's `n` and `tte` plus the slope,
  # taken on first use and kept. lazy_jackknife() (bootstrap.R) is the same
  # memoisation run_bootstrap() uses for its own single jackknife, and owns the
  # convention that the slope accessor is appended last -- so this driver reads
  # `slope_col()` rather than counting columns to find it.
  jack <- lazy_jackknife(setup$frame, setup$subject_index, setup$refitter, computes)

  slope_int <- boot_interval(good_slopes, params$slope, jack$slope_col, type, probs, context,
                             " for the replicate slopes")

  n_res <- vector("list", g$n_cells)
  starved <- character(0L)
  for (k in seq_len(g$n_cells)) {
    label <- grid_boot_row_label(g$labels_at(k))
    n_res[[k]] <- grid_boot_cell_stat(mat$replicates[, k], function() jack$col(k),
                                      pts$n[k], type, probs, context,
                                      sprintf(" for cell %s", label), lattice = TRUE)
    if (isTRUE(n_res[[k]]$starved)) starved <- c(starved, label)
  }
  if (length(starved) == g$n_cells) {
    stop(sprintf(paste0("%s: every cell had fewer than 2 surviving replicates for `n`; no ",
                        "interval could be built."), context), call. = FALSE)
  }

  # One interval per distinct `effectiveness` level rather than one per cell,
  # matching the one replicate column per level grid_boot_computes() filled:
  # every cell at a level shares that column, so a per-cell pass would rebuild
  # the identical interval -- and, at type = "bca", read the identical
  # jackknife column -- once for every design and dropout in the grid. Spread
  # back over the cells by `tte_of` when the columns are assembled below, so
  # the table still reports every cell's own row.
  tte_res <- lapply(seq_along(cc$tte), function(j) {
    ti <- g$n_cells + j
    # The first cell at this level, for a label naming a real row of the grid.
    k <- match(j, cc$tte_of)
    grid_boot_cell_stat(mat$replicates[, ti], function() jack$col(ti),
                        pts$tte[k], type, probs, context,
                        sprintf(" for cell %s (tte)", grid_boot_row_label(g$labels_at(k))),
                        lattice = FALSE)
  })[cc$tte_of]

  report_collected(context, starved, g$n_cells,
                   paste0("fewer than two replicates succeeded for `n` (%s), so no interval ",
                          "could be built for it there. Its `n_*`/`tte_*` columns are NA for ",
                          "those rows."))

  extract <- function(res, field, template) vapply(res, function(r) r[[field]], template)

  added <- list(
    n_mean = extract(n_res, "mean", numeric(1L)),
    n_sd = extract(n_res, "sd", numeric(1L)),
    n_lower = vapply(n_res, function(r) r$ci[1L], numeric(1L)),
    n_upper = vapply(n_res, function(r) r$ci[2L], numeric(1L)),
    tte_mean = extract(tte_res, "mean", numeric(1L)),
    tte_sd = extract(tte_res, "sd", numeric(1L)),
    tte_lower = vapply(tte_res, function(r) r$ci[1L], numeric(1L)),
    tte_upper = vapply(tte_res, function(r) r$ci[2L], numeric(1L)),
    tte_ci_type = extract(tte_res, "type", character(1L)),
    ci_type = extract(n_res, "type", character(1L)),
    n_failed = extract(n_res, "n_failed", integer(1L))
  )

  df <- as.data.frame(c(g$out, pts, added), stringsAsFactors = FALSE)

  # `added_cols` rather than a constant listing these names a second time: the
  # printed cells frame is defined by subtracting them from the table, so a
  # column added above and forgotten in a hand-written list would leak straight
  # into that frame. Carried as an attribute, and so stripped by `[` along with
  # the rest of the grid-wide summary.
  do.call(structure,
          c(list(df, class = c("slope_sample_size_grid_boot", "data.frame"),
                 R = R, type = type, level = level, se = setup$se,
                 n_refit_failed = n_refit_failed,
                 added_cols = names(added),
                 named = g$named),
            slope_replicate_summary(params$slope, good_slopes, slope_int)))
}

#' Subsetting drops the grid-wide summary, not just the class
#'
#' Every attribute [slope_sample_size_grid_boot()] adds -- `R`, the slope
#' block, `named` and the rest -- describes the *whole table*: `R` replicates
#' were drawn once for every cell, `slope_replicates` holds the ones that
#' priced all `nrow(x)` of them, and `named` says which axes vary across the
#' lot. A `[` subset can change which cells are present, or how many, or
#' leave only one level of an axis standing, without changing any of that --
#' and base `[.data.frame` does not reliably strip attributes it does not
#' recognise, so, empirically, those become stale rather than absent:
#' `slope_replicates` still summarised the resampling behind cells a row
#' subset had dropped. So this method strips them itself, on every subset,
#' along with the class -- explicitly, rather than relying on incidental
#' attribute-dropping [print.slope_sample_size_grid_boot()] could not safely
#' assume.
#'
#' Every attribute *but* the base data-frame ones (`names`, `row.names`,
#' `class`) is stripped, rather than naming each of [slope_sample_size_grid_boot()]'s
#' own attributes here too -- the second list that could drift from the
#' first if one gained an attribute the other forgot, the same hazard
#' `stage_two_result()` (power.R) exists to avoid for the result classes'
#' field lists.
#'
#' A subsetted table is then unambiguously a plain data frame: still every
#' cell that survived the subset, just no longer carrying a summary that no
#' longer matches it.
#' @param x A `slope_sample_size_grid_boot` object to subset.
#' @param ... Passed on to `[.data.frame`.
#' @export
`[.slope_sample_size_grid_boot` <- function(x, ...) {
  out <- NextMethod()
  if (inherits(out, "data.frame")) {
    class(out) <- setdiff(class(out), "slope_sample_size_grid_boot")
    extra <- setdiff(names(attributes(out)), c("names", "row.names", "class"))
    for (a in extra) attr(out, a) <- NULL
  }
  out
}

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
#'   `used` is the per-cell `ci_type` -- `n`'s or `tte`'s, since either can fall
#'   back without the other; `asked` the grid-wide `type`.
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
#' Built rather than printed, so [print.slope_sample_size_grid_boot()] can
#' gate the notes beneath both frames on what they actually show.
#' @noRd
grid_boot_tte_frame <- function(x, type) {
  first <- !duplicated(x$effectiveness)
  data.frame(
    effectiveness = x$effectiveness[first],
    tte = x$tte[first],
    tte_mean = x$tte_mean[first],
    tte_sd = x$tte_sd[first],
    # `used`/`asked` as the first frame passes them: `tte`'s interval has its
    # own bias correction and its own jackknife column, so it can fall back to
    # percentile where `n`'s did not, and the `*` has to say so here too or
    # the header note claims a method this column did not use.
    tte_ci = grid_boot_ci_col(x$tte_lower[first], x$tte_upper[first],
                              x$tte_ci_type[first], type),
    stringsAsFactors = FALSE)
}

#' @describeIn slope_sample_size_grid_boot Print a bootstrapped sample-size grid.
#'
#' Printed as data frames, by R itself, rather than as a hand-drawn table: a
#' grid *is* a data frame, [slope_sample_size_grid()] prints as one, and this
#' shows exactly that table plus three columns --- `n_mean`, `n_sd`, and the
#' interval as one `n_ci` column reading "924 to 1704". A second frame follows
#' it for the target treatment effect when `effectiveness` varies, since that
#' is the only axis it depends on. Everything shared by every cell --- the
#' method, the replicate count, the interval level, and the failure and
#' sign-straddling counts --- is reported in the notes beneath, which print on
#' every call whether or not anything went wrong.
#'
#' The resampled slope is not shown. It is still on the object, in the
#' `slope_observed` attribute and the five beside it, and the straddle note
#' still reports the one thing about it that bears on whether these intervals
#' mean anything.
#'
#' Falls back to `print.data.frame()` when the bootstrap attributes are
#' missing -- as they are after e.g. `x[1:3, ]`, since base `[.data.frame`
#' drops attributes it does not know about. A row-subsetted table is still
#' every bit a data frame; it just no longer carries the summary this method
#' exists to show.
#'
#' @param x The `slope_sample_size_grid_boot` object to print. Written out
#'   explicitly, unlike `print.slope_bootstrap()`'s own `x` -- this topic's
#'   `@inheritParams slope_bootstrap` would otherwise leak *that* generic's
#'   unrelated `x` (an object to bootstrap, not one to print) in here, since
#'   both share this help topic and both happen to have a parameter of the
#'   same name.
#' @param ... Not used.
# The grid columns are recovered by dropping the `added_cols` attribute -- which
# slope_sample_size_grid_boot() fills from the names of the block it actually
# built -- rather than by listing them here, through this class's own `[`
# method, which already returns a plain data frame stripped of the grid-wide
# summary.
#' @export
print.slope_sample_size_grid_boot <- function(x, ...) {
  if (is.null(attr(x, "R"))) {
    print.data.frame(x, ...)
    return(invisible(x))
  }

  type <- attr(x, "type")
  level <- attr(x, "level")
  R <- attr(x, "R")
  n_refit_failed <- attr(x, "n_refit_failed")
  named <- attr(x, "named")

  # `[` on this class strips the class and every grid-wide attribute, so this
  # is already the plain data frame print.data.frame() should be handed --
  # no unclass() here, and no second place that knows which attributes exist.
  cells <- x[, setdiff(names(x), attr(x, "added_cols")), drop = FALSE]
  cells$n_mean <- x$n_mean
  cells$n_sd <- x$n_sd
  cells$n_ci <- grid_boot_ci_col(x$n_lower, x$n_upper, x$ci_type, type)

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
  tte <- if (isTRUE(unname(named["effectiveness"]))) grid_boot_tte_frame(x, type)

  cat("<slope_sample_size_grid_boot>\n\n")
  print.data.frame(cells)
  cat("\n")
  if (!is.null(tte)) {
    print.data.frame(tte)
    cat("\n")
  }

  # Every interval column on the page, so that the two markers' notes below
  # follow the markers wherever they appear rather than only where `n` put
  # them. `tte` is NULL when its frame was not shown, and c() drops it.
  shown_ci <- c(cells$n_ci, tte$tte_ci)

  cat(boot_method_note(type, R, level), sep = "\n")

  cat(boot_note("Note", sprintf(paste0(
    "%d/%d (%.1f%%) bootstrap replicates failed to refit the stage-one model, and were ",
    "discarded from every cell."), n_refit_failed, R, 100 * n_refit_failed / R)), sep = "\n")

  # `x$n_failed` counts every replicate that yielded no `n` for the cell,
  # refit failures included; the refits are reported on their own above, so
  # what is left to report here is the difference -- the losses that were the
  # stage-two solve's own.
  extra_failed <- x$n_failed - n_refit_failed
  if (any(extra_failed > 0L)) {
    cat(boot_note("Note", sprintf(paste0(
      "cells also lost replicates solving for `n` itself, beyond the refits above ",
      "(%d-%d of %d failed per cell)."), min(extra_failed), max(extra_failed), R)), sep = "\n")
  }
  if (any(shown_ci == "--", na.rm = TRUE)) {
    cat(boot_note("Note", paste(
      "cells marked \"--\" had fewer than two surviving replicates for that quantity;",
      "no interval could be built for them.")), sep = "\n")
  }

  straddle <- attr(x, "straddle")
  n_used <- length(attr(x, "slope_replicates"))
  cat(boot_note("Note", sprintf(paste0(
    "%d/%d (%.1f%%) of replicates refit a slope on the opposite side of zero from the ",
    "fitted one."), round(straddle * n_used), n_used, 100 * straddle)), sep = "\n")

  cat(boot_note("Mean, SD", paste(
    "each replicate of `n` is rounded up to a whole participant per arm before averaging,",
    "so its mean is not a runnable trial size.")), sep = "\n")

  if (any(grepl("*", shown_ci, fixed = TRUE), na.rm = TRUE)) {
    cat("  * percentile interval; BCa could not be built there.\n")
  }

  invisible(x)
}
