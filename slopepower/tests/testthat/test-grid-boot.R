# Layer 4 --- slope_sample_size_grid_boot().
#
# The bootstrap counterpart of slope_sample_size_grid(): every design shares
# one set of resampled replicates rather than being bootstrapped
# independently, because the resampling scheme depends only on the stage-one
# `params`, never on the design being priced. These tests exist to pin that
# sharing -- a one-cell grid must reproduce slope_sample_size_boot() exactly, and a
# many-cell grid must draw the same replicates as a one-cell one at the same
# seed -- as much as to test the table itself.
#
# R is kept small throughout, as in test-bootstrap.R: these test the
# machinery, not the coverage.

# small_fit() is in helper-fits.R, shared with test-bootstrap.R.

# The same, with a treated arm: the only shape `target = "observed"` accepts,
# since it reuses a previous trial's own observed treatment effect.
small_treated_fit <- function() {
  subj <- local({
    set.seed(11)
    data.frame(id = 1:12, arm = rep(0:1, 6),
               a = rnorm(12, 50, 8), b = rnorm(12, -2, 0.6))
  })
  subj$b <- subj$b + 0.8 * subj$arm
  d <- merge(subj, data.frame(visit = 0:3))
  d$y <- d$a + d$b * d$visit + rnorm(nrow(d), 0, 2)
  slope_params(y ~ visit | id, d, comparator = "treated", group = arm)
}

# --- equivalence to the single-design bootstraps -----------------------------

test_that("a one-cell grid reproduces slope_sample_size_boot() exactly, at the same seed", {
  pars <- small_fit()
  design <- c(0, 1, 2, 3)

  one <- suppressWarnings(slope_sample_size_boot(
    pars, design, effectiveness = 0.33, R = 12, seed = 2))
  grid <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = design, effectiveness = 0.33, R = 12, seed = 2))

  # Same names, same numbers: a grid row's summary columns are the single
  # result's fields, for both bases of the sample size.
  for (q in c("n_total", "n_per_arm")) {
    for (suffix in c("", "_mean", "_sd", "_lower", "_upper")) {
      nm <- paste0(q, suffix)
      expect_equal(grid[[nm]], one[[nm]], info = nm)
    }
  }
  expect_equal(grid$replicates_failed, one$replicates_failed)
  expect_equal(grid$n_ci_method, one$n_ci_method)

  # The slope block -- built once for the whole grid -- must agree with the
  # single result's own slope summary too, since both are the same
  # computation on the same draws.
  expect_equal(attr(grid, "slope_lower"), one$slope_lower)
  expect_equal(attr(grid, "slope_upper"), one$slope_upper)
  expect_equal(attr(grid, "slope_mean"), one$slope_mean)
  expect_equal(attr(grid, "slope_sd"), one$slope_sd)
  expect_equal(attr(grid, "straddle"), one$straddle)
  expect_identical(attr(grid, "slope_replicates"), one$slope_replicates)
})

test_that("a one-cell grid's point estimates match slope_sample_size_grid()", {
  # Both carry both bases of every count (CONTRACT.md section 4.5), so every
  # column of the plain grid is a column of the bootstrapped one, by name.
  pars <- small_fit()
  plain <- slope_sample_size_grid(pars, visits = c(0, 1, 2, 3), effectiveness = 0.33)
  boot <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 1))

  expect_s3_class(boot, "slope_sample_size_grid")
  for (nm in names(plain)) expect_equal(boot[[nm]], plain[[nm]], info = nm)
})

# --- shared replicates -------------------------------------------------------

test_that("every cell of a grid draws the same replicates, at the same seed", {
  pars <- small_fit()

  one <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 10, seed = 7))
  many <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 10, seed = 7,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  # Can only agree if the same subjects were resampled and refit once, not
  # once per cell: four cells' worth of independent draws would not, in
  # general, reproduce the one-cell run's slopes.
  expect_identical(attr(one, "slope_replicates"), attr(many, "slope_replicates"))
})

# --- table shape --------------------------------------------------------------

test_that("slope_sample_size_grid_boot() returns one row per combination, with the expected columns", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 6, seed = 3,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  expect_s3_class(out, "slope_sample_size_grid_boot")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 4L)
  expect_true(all(c("n_total_mean", "n_total_sd", "n_total_lower", "n_total_upper",
                    "n_per_arm_mean", "n_per_arm_sd", "n_per_arm_lower", "n_per_arm_upper",
                    "tte_mean", "tte_sd", "tte_lower", "tte_upper",
                    "tte_ci_method", "n_ci_method", "replicates_failed") %in% names(out)))
  # Every per-arm summary is exactly half its total.
  for (suffix in c("_mean", "_sd", "_lower", "_upper")) {
    expect_equal(out[[paste0("n_total", suffix)]], 2 * out[[paste0("n_per_arm", suffix)]])
  }
})

test_that("a bootstrapped `n` interval is always widened to even sizes", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 8, seed = 5,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  expect_true(all(out$n_total_lower %% 2 == 0))
  expect_true(all(out$n_total_upper %% 2 == 0))
})

test_that("`ci_method = \"percentile\"` is honoured cell by cell, without a jackknife", {
  pars <- small_fit()
  out <- slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3), effectiveness = 0.33,
                                     R = 6, seed = 1, ci_method = "percentile")
  expect_true(all(out$n_ci_method == "percentile"))
})

test_that("the target treatment effect is constant across design and dropout, and varies with effectiveness", {
  pars <- small_fit()

  # Fixed effectiveness: tte (and its bootstrap summary) does not depend on
  # the design, so every cell reports the same value.
  same_eff <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 4,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))
  expect_equal(length(unique(same_eff$tte)), 1L)
  expect_equal(length(unique(same_eff$tte_mean)), 1L)

  # effectiveness itself an axis: tte varies with it, one value per level.
  vary_eff <- slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3),
                                          effectiveness = c(0.25, 0.33),
                                          R = 5, seed = 4)
  expect_equal(length(unique(vary_eff$tte)), 2L)
})

# --- randomness --------------------------------------------------------------

test_that("a seeded call restores the caller's random number stream", {
  pars <- small_fit()
  set.seed(99)
  before <- .Random.seed
  invisible(suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 123)))
  expect_identical(.Random.seed, before)
})

# --- failure accounting -------------------------------------------------------

test_that("grid_boot_cell_stat() reports a cell as starved rather than erroring", {
  grid_boot_cell_stat <- slopepower:::grid_boot_cell_stat
  # Only one non-NA value: too few to form any interval.
  res <- grid_boot_cell_stat(c(10, NA, NA, NA), function() c(9, 11, 10), observed = 10,
                             ci_method = "bca", probs = c(0.025, 0.975), context = "test",
                             what = "", statistic = "n")
  expect_true(res$starved)
  expect_true(is.na(res$mean))
  expect_true(all(is.na(res$ci)))
  expect_equal(res$n_failed, 3L)
})

test_that("too few refits succeeding aborts the whole grid, as it does for a single-design bootstrap", {
  # R = 1 means at most one replicate can possibly succeed -- fewer than the
  # two run_bootstrap() (bootstrap.R) itself requires for an interval -- so
  # this is deterministic however well-behaved the fit is, unlike forcing an
  # actual convergence failure.
  pars <- small_fit()
  expect_error(
    suppressWarnings(slope_sample_size_grid_boot(pars, visits = c(0, 1, 2), R = 1, seed = 1)),
    "not enough succeeded")
})

# --- printing ------------------------------------------------------------------

# The main frame and nothing else: everything between the title and its blank
# line, and the next blank line.
frame_block <- function(out) {
  blanks <- which(out == "")
  out[seq(blanks[1L] + 1L, blanks[2L] - 1L)]
}

# The column names of a printed frame. R wraps a wide frame into blocks, each
# headed by its own column names and followed by rows that begin with a row
# number, so the names are read off the lines that are not rows.
printed_columns <- function(frame) {
  headers <- frame[!grepl("^[0-9]", frame)]
  Filter(nzchar, unlist(strsplit(trimws(headers), " +")))
}

test_that("print.slope_sample_size_grid_boot() prints one basis under the stored names", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))
  expect_identical(attr(out, "per_arm"), TRUE)

  lines <- capture.output(print(out))
  expect_identical(lines[1L], "Bootstrapped sample size grid (4 rows)")

  printed <- printed_columns(frame_block(lines))
  expect_true(all(printed %in% names(out)))
  expect_true(all(c("n_per_arm", "visits_per_arm", "n_per_arm_mean", "n_per_arm_sd",
                    "n_per_arm_lower", "n_per_arm_upper") %in% printed))
  expect_false(any(grepl("_total(_|$)", printed)))
  expect_false(any(grepl("^tte_|ci_method|replicates_failed", printed)))

  # The printed values are the stored per-arm ones, not a re-derivation.
  expect_true(any(grepl(format(out$n_per_arm_mean[1L], digits = 7), lines, fixed = TRUE)))

  expect_true(any(grepl("Bootstrap:   R = 5 replicates; 95% BCa intervals.",
                        lines, fixed = TRUE)))
  # A clean run (no refit failures, no starved cells) still reports both
  # notes, per the convention set for print.slope_bootstrap().
  expect_true(any(grepl("bootstrap samples failed to refit", lines)))
  expect_true(any(grepl("of replicates refit a slope on the opposite side", lines)))
})

test_that("print.slope_sample_size_grid_boot() shows the totals under per_arm = FALSE", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  printed <- printed_columns(frame_block(capture.output(print(out, per_arm = FALSE))))
  expect_true(all(c("n_total", "visits_total", "n_total_mean", "n_total_lower") %in% printed))
  expect_false(any(grepl("_per_arm", printed)))

  # The object itself is untouched by which basis was printed.
  expect_identical(attr(out, "per_arm"), TRUE)
  expect_equal(out$n_total, 2L * out$n_per_arm)
  expect_equal(out$visits_total, 2L * out$visits_per_arm)
})

test_that("print(x, per_arm = ...) overrides the basis a grid was built with", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 6,
    per_arm = FALSE))
  expect_identical(attr(out, "per_arm"), FALSE)

  default_lines  <- capture.output(print(out))
  override_lines <- capture.output(print(out, per_arm = TRUE))
  expect_true("n_total" %in% printed_columns(frame_block(default_lines)))
  expect_true("n_per_arm" %in% printed_columns(frame_block(override_lines)))

  expect_error(print(out, per_arm = NA), "per_arm")
})

test_that("print.slope_sample_size_grid_boot() does not report the resampled slope", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits = list(a = c(0, 1, 2, 3), b = c(0, 1, 2))))

  frame <- frame_block(capture.output(print(out)))
  expect_false(any(grepl("slope", frame, ignore.case = TRUE)))

  # Still on the object, and the straddle note still reports the one thing
  # about it that bears on whether these intervals mean anything.
  expect_equal(attr(out, "slope"), pars$slope)
  expect_length(c(attr(out, "slope_lower"), attr(out, "slope_upper")), 2L)
})

test_that("print.slope_sample_size_grid_boot() shows a starved cell and a fallen-back interval", {
  # Forced onto a finished grid rather than provoked, so the test does not
  # depend on which resample happens to starve a cell.
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits = list(a = c(0, 1, 2, 3), b = c(0, 1, 2))))
  clean <- capture.output(print(out))
  expect_false("n_ci_method" %in% printed_columns(frame_block(clean)))
  expect_false(any(grepl("_ci_method` column", clean, fixed = TRUE)))

  out$n_per_arm_lower[1L] <- NA_real_
  out$n_per_arm_upper[1L] <- NA_real_
  out$n_ci_method[1L] <- NA_character_
  out$n_ci_method[2L] <- "percentile"

  lines <- capture.output(print(out))
  frame <- frame_block(lines)
  # The method column appears once a cell differs from the method asked for.
  expect_true("n_ci_method" %in% printed_columns(frame))
  expect_true(any(grepl("percentile", frame, fixed = TRUE)))
  expect_true(any(grepl("\\bNA\\b", frame)))

  # And each is explained beneath, in a note that appears only when the thing
  # it explains is on the page.
  expect_true(any(grepl("an interval shown as NA had fewer than two surviving", lines)))
  expect_true(any(grepl("_ci_method` column marks the cells", lines, fixed = TRUE)))
})

test_that("print.slope_sample_size_grid_boot() shows a fallen-back tte interval too", {
  # `tte`'s interval has its own bias correction and its own jackknife column,
  # so it can fall back to percentile where `n`'s did not.
  pars <- small_fit()
  out <- slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3),
                                     effectiveness = c(0.25, 0.33), R = 5, seed = 6)
  expect_true(all(out$n_ci_method == "bca"))
  out$tte_ci_method[1L] <- "percentile"
  out$tte_lower[2L] <- NA_real_
  out$tte_upper[2L] <- NA_real_
  out$tte_ci_method[2L] <- NA_character_

  lines <- capture.output(print(out))
  expect_true(any(grepl("tte_ci_method", lines, fixed = TRUE)))
  expect_true(any(grepl("_ci_method` column marks the cells", lines, fixed = TRUE)))
  expect_true(any(grepl("an interval shown as NA had fewer than two surviving", lines)))
})

test_that("the extra-failure note reports the losses beyond the refits, not the totals", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 10, seed = 6,
    visits = list(a = c(0, 1, 2, 3), b = c(0, 1, 2))))

  # Three replicates lost to the refit, shared by every cell; one cell lost
  # four more to its own stage-two solve, the other none. The note is about
  # the "more", so it must read 0-4 and not the 3-7 totals.
  attr(out, "refits_failed") <- 3L
  out$replicates_failed <- c(3L, 7L)
  # boot_note() wraps a note across lines, so the notes are read as one string
  # rather than line by line: the phrase this is about spans a break.
  notes <- function(x) gsub(" +", " ", paste(capture.output(print(x)), collapse = " "))
  expect_match(notes(out),
               "beyond the refits above (0-4 of 10 failed per cell).", fixed = TRUE)

  # And a grid that lost nothing beyond the refits does not raise the note.
  out$replicates_failed <- c(3L, 3L)
  expect_false(grepl("beyond the refits above", notes(out), fixed = TRUE))
})

test_that("print.slope_sample_size_grid_boot() adds a target-effect frame only when effectiveness varies", {
  pars <- small_fit()
  fixed <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 6))
  varying <- slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3),
                                         effectiveness = c(0.25, 0.33), R = 5, seed = 6)

  expect_false(any(grepl("tte_mean", capture.output(print(fixed)), fixed = TRUE)))

  lines <- capture.output(print(varying))
  for (col in c("tte_mean", "tte_sd", "tte_lower", "tte_upper")) {
    expect_true(any(grepl(col, lines, fixed = TRUE)))
  }
  # Keyed by `effectiveness` alone: it is the only axis tte depends on.
  expect_true(any(grepl("^ *effectiveness +tte +tte_mean", lines)))
})

test_that("the target-effect frame carries one row per effectiveness level, not per cell", {
  # tte does not depend on the design, so a frame keyed by design and dropout
  # would print the same two intervals four times over.
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = c(0.25, 0.33), R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))
  expect_equal(nrow(out), 8L)

  frame <- slopepower:::grid_boot_tte_frame(strip_grid(out), FALSE)
  expect_equal(nrow(frame), 2L)
  expect_equal(frame$effectiveness, c(0.25, 0.33))
  expect_false(any(c("design", "dropout") %in% names(frame)))

  # And the cells sharing a level share the interval exactly, not merely to
  # rounding: they were priced off one replicate column, not four identical ones.
  for (col in c("tte", "tte_mean", "tte_sd", "tte_lower", "tte_upper")) {
    expect_equal(length(unique(out[[col]])), 2L, info = col)
  }
})

test_that("a grid solved for target = \"observed\" prints", {
  # `effectiveness` is not an axis of such a grid -- it is not passed at all --
  # so the `named` vector the tte frame is gated on carries no element of that
  # name, and the guard must not subscript one out of bounds.
  pars <- small_treated_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    target = "observed", R = 6, seed = 8))

  lines <- capture.output(print(out))
  expect_identical(lines[1L], "Bootstrapped sample size grid (2 rows)")
  # No second frame: with one target effect for the whole grid there is
  # nothing for it to show that the notes do not.
  expect_false(any(grepl("tte_mean", lines, fixed = TRUE)))
  # The notes beneath still print, which they cannot if the frames aborted.
  expect_true(any(grepl("bootstrap samples failed to refit", lines)))
})

test_that("confint() on a bootstrapped grid reads each cell's stored interval", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  ci <- confint(out)
  expect_equal(dim(ci), c(4L, 2L))
  expect_identical(colnames(ci), c("2.5 %", "97.5 %"))
  expect_identical(rownames(ci), paste(out$design, out$dropout, sep = " / "))
  expect_equal(unname(ci[, 1L]), out$n_per_arm_lower)
  expect_equal(unname(ci[, 2L]), out$n_per_arm_upper)
  expect_equal(unname(confint(out, per_arm = FALSE)[, 1L]), out$n_total_lower)
  expect_equal(unname(confint(out, "tte")[, 2L]), out$tte_upper)
  expect_error(confint(out, "n"), "must be one of")
  expect_error(confint(out, level = 0.9), "level = 0.95")
})

# --- row/column subsetting ------------------------------------------------------

test_that("subsetting a bootstrapped grid drops the class and the grid-wide summary", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  row_sub <- out[1L, ]
  expect_false(inherits(row_sub, "slope_sample_size_grid_boot"))
  expect_null(attr(row_sub, "R"))
  expect_null(attr(row_sub, "slope_replicates"))
  expect_null(attr(row_sub, "named"))
  expect_null(attr(row_sub, "per_arm"))
  expect_equal(row_sub$n_total, out$n_total[1L])

  col_sub <- out[, c("design", "n_total")]
  expect_false(inherits(col_sub, "slope_sample_size_grid_boot"))

  # The original object is untouched by taking a subset of it.
  expect_s3_class(out, "slope_sample_size_grid_boot")
  expect_false(is.null(attr(out, "R")))
})

# --- reused validation --------------------------------------------------------

test_that("slope_sample_size_grid_boot() rejects effectiveness alongside target = \"observed\"", {
  pars <- small_fit()
  expect_error(
    slope_sample_size_grid_boot(pars, visits = c(0, 1, 2), effectiveness = 0.33,
                                target = "observed", R = 4),
    "target = \"observed\""
  )
})

test_that("slope_sample_size_grid_boot() collects the baseline-dropout warning once, as the plain grid does", {
  pars <- small_fit()
  expect_warning(
    slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3), effectiveness = 0.33,
                                dropout = c(0.3, 0.1, 0.1), R = 4, seed = 1),
    "baseline visit only"
  )
})

test_that("a one-cell power grid reproduces slope_power_boot() exactly, at the same seed", {
  pars <- small_fit()
  design <- c(0, 1, 2, 3)

  one <- suppressWarnings(slope_power_boot(
    pars, design, n_total = 60, effectiveness = 0.33, R = 12, seed = 2))
  grid <- suppressWarnings(slope_power_grid_boot(
    pars, visits = design, n_total = 60, effectiveness = 0.33, R = 12, seed = 2))

  expect_s3_class(grid, "slope_power_grid_boot")
  expect_identical(attr(grid, "statistic"), "power")
  for (suffix in c("", "_mean", "_sd", "_lower", "_upper")) {
    nm <- paste0("power", suffix)
    expect_equal(grid[[nm]], one[[nm]], info = nm)
  }
  expect_equal(grid$replicates_failed, one$replicates_failed)
  expect_equal(grid$power_ci_method, one$power_ci_method)
  expect_equal(attr(grid, "slope_lower"), one$slope_lower)
})

test_that("a power grid bootstrap prints its power intervals and no rounding note", {
  pars <- small_fit()
  out <- suppressWarnings(slope_power_grid_boot(
    pars, visits = list(a = c(0, 1, 2), b = c(0, 1, 2, 3)), n_total = c(60, 80),
    effectiveness = 0.33, R = 6, seed = 1))
  expect_equal(nrow(out), 4L)
  expect_true(all(c("power_mean", "power_sd", "power_lower", "power_upper",
                    "tte_ci_method", "power_ci_method", "replicates_failed") %in% names(out)))
  expect_false(any(grepl("^n_(total|per_arm)_(mean|lower)$", names(out))))
  printed <- capture.output(print(out))
  expect_identical(printed[1L], "Bootstrapped power grid (4 rows)")
  expect_true(any(grepl("power_lower", printed, fixed = TRUE)))
  expect_false(any(grepl("rounded up", printed, fixed = TRUE)))
  # Subsetting drops the class and the grid-wide summary, as for the
  # sample-size grid.
  expect_false(inherits(out[1:2, ], "slope_power_grid_boot"))
  expect_null(attr(out[1:2, ], "R"))
  expect_error(slope_power_grid_boot(pars, visits = c(0, 1, 2)),
               "one of `n_total` and `n_per_arm` is required")
})
