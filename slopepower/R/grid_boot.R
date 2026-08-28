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

#' Two replicate-statistic closures per cell: the sample size and its target
#' treatment effect
#'
#' Every cell of the grid gets its own pair of `function(p) numeric(1)`
#' closures, each re-solving [slope_sample_size()] against one resampled
#' replicate's refitted parameters `p` and the cell's own design and scalar
#' values. `power` is spliced in as the fixed argument the way
#' `bootstrap_stage_two()` (bootstrap.R) does for a single result;
#' [resolve_args()] -- the same helper -- supplies `design`, `target`,
#' `alpha` and `effectiveness`.
#'
#' `tte` is re-solved by its own call rather than read off an `n`-solve
#' already made for the replicate: [boot_replicate_matrix()] and
#' [jackknife_values()] (both bootstrap.R) take a flat list of
#' `function(p) numeric(1)` accessors, and keeping every accessor to that
#' shape means neither has to learn that two of them happen to come from the
#' same underlying solve. The price is one extra call to
#' [slope_sample_size()] per cell per replicate -- arithmetic, not a model
#' fit, and small next to what a replicate already costs.
#'
#' Each closure closes over `slim` and `power_k` alone, not over `params` or
#' the grid itself, so a several-hundred-replicate run does not keep the
#' original fit -- and its model frame -- reachable through every closure.
#' `lapply()` rather than a loop, so each cell's closures capture their own
#' `slim`/`power_k` in a fresh call frame instead of the last value a shared
#' loop variable happened to hold.
#'
#' @return A list of `g$n_cells` elements, each `list(n = <closure>,
#'   tte = <closure>)`.
#' @noRd
grid_boot_computes <- function(g, target) {
  cell_args <- grid_cell_args(g)
  lapply(seq_len(g$n_cells), function(k) {
    slim <- list(design = g$designs[[g$design_of[k]]], target = target,
                alpha = cell_args[[k]]$alpha, effectiveness = cell_args[[k]]$effectiveness)
    power_k <- cell_args[[k]]$power
    list(
      n = function(p) do.call(slope_sample_size, c(resolve_args(p, slim), list(power = power_k)))$n,
      tte = function(p) do.call(slope_sample_size,
                                c(resolve_args(p, slim), list(power = power_k)))$tte
    )
  })
}

#' Flatten the per-cell compute list into the shape [boot_replicate_matrix()]
#' and [jackknife_values()] both take
#'
#' `n` and `tte` interleaved, cell by cell, so column `2 * k - 1` is always
#' cell `k`'s sample size and `2 * k` its target treatment effect -- the one
#' indexing rule the rest of [slope_sample_size_grid_boot()] relies on.
#' @noRd
grid_boot_flatten <- function(cell_computes) {
  unlist(lapply(cell_computes, function(cc) list(cc$n, cc$tte)), recursive = FALSE)
}

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
#'   "data.frame")`, one row per cell, with the fourteen columns
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
#'       itself an axis of the grid.}
#'     \item{`ci_type`}{`"bca"` or `"percentile"`: the method actually used
#'       for that cell's `n` interval, which can fall back to percentile cell
#'       by cell even when `type = "bca"` was requested. `NA` for a cell that
#'       could not be bootstrapped at all.}
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
#'   perturbs and every cell's interval is a function of it. Base
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
  type <- match.arg(type, c("bca", "percentile"))
  check_whole_number(R, "R", "replicates", context, lower = 1)
  check_probability(level, "level", context)

  # Seeded calls are reproducible without reseeding the session, exactly as
  # run_bootstrap() (bootstrap.R) arranges for a single result; see its own
  # note on why on.exit() rather than a restore at the end.
  if (!is.null(seed)) {
    old_seed <- current_seed()
    on.exit(restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  # The point estimates: built by exactly the code path slope_sample_size_grid()
  # itself uses (grid_stage_two()'s own evaluate() closure, reproduced here),
  # so this table's `n`/`tte`/... columns can never drift from that function's.
  scalars <- maybe_add_effectiveness(stats::setNames(list(power), "power"), effectiveness, target)
  scalars$alpha <- alpha
  g <- grid_axes(visits, dropout, scalars, context)
  pts <- grid_evaluate(g,
                       function(des, args) do.call(slope_sample_size,
                                                   c(list(params = params, design = des,
                                                          target = target), args)),
                       context)

  computes <- grid_boot_flatten(grid_boot_computes(g, target))
  n_slots <- 2L * g$n_cells

  setup <- boot_setup(params, context)
  mat <- boot_replicate_matrix(setup, computes, R, progress, context)

  failed_refit <- is.na(mat$slopes)
  n_refit_failed <- sum(failed_refit)
  good_slopes <- mat$slopes[!failed_refit]
  if (length(good_slopes) < 2L) {
    stop(sprintf(paste0("%s: %d of %d replicates failed to refit; not enough succeeded to ",
                        "bootstrap any cell."), context, n_refit_failed, R), call. = FALSE)
  }

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)

  # One jackknife pass, covering every cell's `n` and `tte` plus the slope,
  # taken on first use and kept -- the same laziness run_bootstrap() applies
  # to its own single jackknife, generalised from one column pair to
  # `n_slots + 1`.
  jack <- NULL
  jack_matrix <- function() {
    if (is.null(jack)) {
      jack <<- jackknife_values(setup$frame, setup$subject_index, setup$refitter,
                                c(computes, list(function(p) p$slope)))
    }
    jack
  }

  slope_int <- boot_interval(good_slopes, params$slope,
                             function() jack_matrix()[, n_slots + 1L], type, probs, context,
                             " for the replicate slopes")

  n_res <- vector("list", g$n_cells)
  tte_res <- vector("list", g$n_cells)
  starved <- character(0L)
  for (k in seq_len(g$n_cells)) {
    ni <- 2L * k - 1L
    ti <- 2L * k
    label <- grid_boot_row_label(g$labels_at(k))
    n_res[[k]] <- grid_boot_cell_stat(mat$replicates[, ni], function() jack_matrix()[, ni],
                                      pts$n[k], type, probs, context,
                                      sprintf(" for cell %s", label), lattice = TRUE)
    tte_res[[k]] <- grid_boot_cell_stat(mat$replicates[, ti], function() jack_matrix()[, ti],
                                        pts$tte[k], type, probs, context,
                                        sprintf(" for cell %s (tte)", label), lattice = FALSE)
    if (isTRUE(n_res[[k]]$starved)) starved <- c(starved, label)
  }
  if (length(starved) == g$n_cells) {
    stop(sprintf(paste0("%s: every cell had fewer than 2 surviving replicates for `n`; no ",
                        "interval could be built."), context), call. = FALSE)
  }
  report_collected(context, starved, g$n_cells,
                   paste0("fewer than two replicates succeeded for `n` (%s), so no interval ",
                          "could be built for it there. Its `n_*`/`tte_*` columns are NA for ",
                          "those rows."))

  extract <- function(res, field, template) vapply(res, function(r) r[[field]], template)

  out <- c(g$out, pts, list(
    n_mean = extract(n_res, "mean", numeric(1L)),
    n_sd = extract(n_res, "sd", numeric(1L)),
    n_lower = vapply(n_res, function(r) r$ci[1L], numeric(1L)),
    n_upper = vapply(n_res, function(r) r$ci[2L], numeric(1L)),
    tte_mean = extract(tte_res, "mean", numeric(1L)),
    tte_sd = extract(tte_res, "sd", numeric(1L)),
    tte_lower = vapply(tte_res, function(r) r$ci[1L], numeric(1L)),
    tte_upper = vapply(tte_res, function(r) r$ci[2L], numeric(1L)),
    ci_type = extract(n_res, "type", character(1L)),
    n_failed = extract(n_res, "n_failed", integer(1L))
  ))

  df <- as.data.frame(out, stringsAsFactors = FALSE)

  structure(df, class = c("slope_sample_size_grid_boot", "data.frame"),
           R = R, type = type, level = level, se = setup$se,
           n_refit_failed = n_refit_failed,
           straddle = mean(sign(good_slopes) != sign(params$slope)),
           slope_observed = params$slope,
           slope_mean = mean(good_slopes), slope_sd = stats::sd(good_slopes),
           slope_ci = slope_int$ci, slope_type = slope_int$type,
           slope_replicates = good_slopes,
           named = g$named,
           row_labels = vapply(seq_len(g$n_cells), function(k) grid_boot_row_label(g$labels_at(k)),
                               character(1L)))
}

#' Subsetting drops the grid-wide summary, not just the class
#'
#' Every attribute [slope_sample_size_grid_boot()] adds -- `R`, the slope
#' block, `row_labels` and the rest -- describes the *whole table*: `R`
#' replicates were drawn once for every cell, and `row_labels` names all
#' `nrow(x)` of them. A `[` subset can change which cells are present, or how
#' many, without changing any of that -- and base `[.data.frame` does not
#' reliably strip attributes it does not recognise, so, empirically, those
#' become stale rather than absent: `row_labels` still named the original
#' cells after a row subset dropped some of them, longer than the data it
#' would label. So this method strips them itself, on every subset, along
#' with the class -- explicitly, rather than relying on incidental
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

#' A width-driven table of one statistic across a labelled set of rows
#'
#' Shared by both tables [print.slope_sample_size_grid_boot()] prints -- the
#' sample size (with the shared slope row prepended) and, when it varies, the
#' target treatment effect -- so the column-width and span-header arithmetic
#' `boot_n_table()` (bootstrap.R) already worked out for a fixed three rows is
#' written once here for however many `labels` names, rather than twice.
#'
#' `sep_after`, if given, draws a second rule after that row -- used to set
#' the shared slope row apart from the per-cell rows beneath it, the same
#' break `boot_n_table()` draws between "Slope" and "Sample size".
#' @noRd
boot_grid_stat_table <- function(labels, calc, means, sds, cis, ci_head, span, sep_after = NULL) {
  w <- c(max(nchar(labels)),
        max(nchar(c("Calculated", calc))),
        max(nchar(c("Mean", means))),
        max(nchar(c("SD", sds))),
        max(nchar(c(ci_head, cis))))

  span_w <- sum(w[3:5]) + 6L
  if (nchar(span) > span_w) {
    w[5L] <- w[5L] + nchar(span) - span_w
    span_w <- nchar(span)
  }

  rule <- paste0("  |", paste(vapply(w + 2L, strrep, character(1L), x = "-"), collapse = "|"), "|")
  align <- c("l", "r", "r", "r", "r")
  row <- function(i) boot_table_row(c(labels[i], calc[i], means[i], sds[i], cis[i]), w, align)

  body <- character(0L)
  for (i in seq_along(labels)) {
    body <- c(body, row(i))
    if (identical(sep_after, i)) body <- c(body, rule)
  }

  c(boot_table_row(c("", "Calculated", span), c(w[1:2], span_w), c("l", "r", "l")),
    boot_table_row(c("", "", "Mean", "SD", ci_head), w, c("l", "l", "r", "r", "r")),
    rule, body)
}

#' @describeIn slope_sample_size_grid_boot Print a bootstrapped sample-size grid.
#'
#' Reuses `boot_table_row()`, `boot_ci_cells()` and `boot_note()`
#' (bootstrap.R) -- all three already width-driven and statistic-agnostic --
#' via `boot_grid_stat_table()`, so a table with several rows costs no new
#' formatting machinery beyond generalising `boot_n_table()`'s own fixed
#' three rows to however many the grid has.
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
  row_labels <- attr(x, "row_labels")
  named <- attr(x, "named")

  fmt <- function(v) format(v, digits = 6, trim = TRUE, drop0trailing = TRUE)
  dash_if_na <- function(v, s) ifelse(is.na(v), "--", s)

  ci_head <- sprintf("%.0f%% CI", 100 * level)
  method <- if (identical(type, "bca")) "BCa" else "percentile"
  span <- sprintf("Bootstrapped (%s, R = %d%s)", method, R,
                  if (n_refit_failed > 0L) sprintf(", %d refit failed", n_refit_failed) else "")

  slope_ci <- attr(x, "slope_ci")
  slope_mixed <- !identical(attr(x, "slope_type"), type)
  slope_ci_str <- boot_ci_cells(fmt(slope_ci[1L]), fmt(slope_ci[2L]))
  if (slope_mixed) slope_ci_str <- paste0(slope_ci_str, " *")

  n_mixed <- !is.na(x$ci_type) & x$ci_type != type
  n_ci_str <- boot_ci_cells(fmt(x$n_lower), fmt(x$n_upper))
  n_ci_str[n_mixed] <- paste0(n_ci_str[n_mixed], " *")
  n_ci_str <- dash_if_na(x$n_lower, n_ci_str)

  cat("<slope_sample_size_grid_boot>\n\n")
  cat(boot_grid_stat_table(
    labels = c("Slope", row_labels),
    calc = c(fmt(attr(x, "slope_observed")), fmt(x$n)),
    means = c(fmt(attr(x, "slope_mean")), dash_if_na(x$n_mean, fmt(x$n_mean))),
    sds = c(fmt(attr(x, "slope_sd")), dash_if_na(x$n_sd, fmt(x$n_sd))),
    cis = c(slope_ci_str, n_ci_str),
    ci_head = ci_head, span = span, sep_after = 1L), sep = "\n")
  cat("\n")

  # The target treatment effect does not depend on the design, only on
  # `effectiveness` -- see slope_sample_size_grid_boot()'s @return -- so a
  # second table for it earns its place only when that axis actually varies;
  # otherwise every row would repeat the same interval.
  show_tte <- isTRUE(named[["effectiveness"]])
  if (show_tte) {
    tte_ci_str <- dash_if_na(x$tte_lower, boot_ci_cells(fmt(x$tte_lower), fmt(x$tte_upper)))
    cat(boot_grid_stat_table(
      labels = row_labels, calc = fmt(x$tte),
      means = dash_if_na(x$tte_mean, fmt(x$tte_mean)),
      sds = dash_if_na(x$tte_sd, fmt(x$tte_sd)),
      cis = tte_ci_str, ci_head = ci_head,
      span = sprintf("Target treatment effect, bootstrapped (%s, R = %d)", method, R)),
        sep = "\n")
    cat("\n")
  }

  cat(boot_note("Note", sprintf(paste0(
    "%d/%d (%.1f%%) bootstrap replicates failed to refit the stage-one model, and were ",
    "discarded from every cell."), n_refit_failed, R, 100 * n_refit_failed / R)), sep = "\n")

  extra_failed <- x$n_failed - n_refit_failed
  if (any(extra_failed > 0L)) {
    cat(boot_note("Note", sprintf(paste0(
      "cells also lost replicates solving for `n` itself, beyond the refits above ",
      "(%d-%d of %d failed per cell)."), min(x$n_failed), max(x$n_failed), R)), sep = "\n")
  }
  if (any(is.na(x$ci_type))) {
    cat(boot_note("Note", paste(
      "cells marked \"--\" had fewer than two surviving replicates for `n`; no interval",
      "could be built for them.")), sep = "\n")
  }

  straddle <- attr(x, "straddle")
  n_used <- length(attr(x, "slope_replicates"))
  cat(boot_note("Note", sprintf(paste0(
    "%d/%d (%.1f%%) of replicates refit a slope on the opposite side of zero from the ",
    "fitted one."), round(straddle * n_used), n_used, 100 * straddle)), sep = "\n")

  cat(boot_note("Mean, SD", paste(
    "each replicate of `n` is rounded up to a whole participant per arm before averaging,",
    "so its mean is not a runnable trial size.")), sep = "\n")

  if (slope_mixed || any(n_mixed)) {
    cat("  * percentile interval; BCa could not be built there.\n")
  }

  invisible(x)
}
