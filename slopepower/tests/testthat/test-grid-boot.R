# Layer 4 --- slope_sample_size_grid_boot().
#
# The bootstrap counterpart of slope_sample_size_grid(): every design shares
# one set of resampled replicates rather than being bootstrapped
# independently, because the resampling scheme depends only on the stage-one
# `params`, never on the design being priced. These tests exist to pin that
# sharing -- a one-cell grid must reproduce slope_bootstrap() exactly, and a
# many-cell grid must draw the same replicates as a one-cell one at the same
# seed -- as much as to test the table itself.
#
# R is kept small throughout, as in test-bootstrap.R: these test the
# machinery, not the coverage.

# Ten subjects, a small model that fits quickly and (at type = "bca") needs
# only ten leave-one-subject-out refits rather than slpower1's two hundred.
# Built on demand, not at load time, so it is not paid for by every other
# test file.
small_fit <- function() {
  subj <- local({
    set.seed(4)
    data.frame(id = 1:10, a = rnorm(10, 50, 8), b = rnorm(10, -2, 0.6))
  })
  d <- merge(subj, data.frame(visit = 0:3))
  d$y <- d$a + d$b * d$visit + rnorm(nrow(d), 0, 2)
  slope_params(y ~ visit | id, d)
}

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
  slope_params(y ~ visit | id, d, treated = arm)
}

# --- equivalence to slope_bootstrap() ---------------------------------------

test_that("a one-cell grid reproduces slope_bootstrap() exactly, at the same seed", {
  pars <- small_fit()
  design <- c(0, 1, 2, 3)

  one <- suppressWarnings(slope_bootstrap(
    slope_sample_size(pars, design, effectiveness = 0.33), R = 12, seed = 2))
  grid <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = design, effectiveness = 0.33, R = 12, seed = 2))

  expect_equal(grid$n_lower, one$ci[1L])
  expect_equal(grid$n_upper, one$ci[2L])
  expect_equal(grid$n_mean, one$boot_mean)
  expect_equal(grid$n_sd, one$boot_sd)
  expect_equal(grid$n_failed, one$n_failed)
  expect_equal(grid$ci_type, one$type)

  # The slope block -- built once for the whole grid -- must agree with the
  # single result's own slope summary too, since both are the same
  # computation on the same draws.
  expect_equal(attr(grid, "slope_ci"), one$slope_ci)
  expect_equal(attr(grid, "slope_mean"), one$slope_mean)
  expect_equal(attr(grid, "slope_sd"), one$slope_sd)
  expect_equal(attr(grid, "straddle"), one$straddle)
  expect_identical(attr(grid, "slope_replicates"), one$slope_replicates)
})

test_that("a one-cell grid's point estimates match slope_sample_size_grid()", {
  # The boot grid's own data frame always carries both bases unreduced
  # (CONTRACT.md section 4.4), while slope_sample_size_grid() itself returns
  # one basis at a time, under one name ("n"/"visits") for either -- so the
  # equivalence is checked once per basis for the two columns that differ,
  # and once, basis-independent, for the rest.
  pars <- small_fit()
  plain_arm   <- slope_sample_size_grid(pars, visits = c(0, 1, 2, 3), effectiveness = 0.33)
  plain_total <- slope_sample_size_grid(pars, visits = c(0, 1, 2, 3), effectiveness = 0.33,
                                        per_arm = FALSE)
  boot <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 1))

  expect_equal(boot$n, plain_total$n)
  expect_equal(boot$n_per_arm, plain_arm$n)
  expect_equal(boot$visits_total, plain_total$visits)
  expect_equal(boot$visits_per_arm, plain_arm$visits)

  shared <- setdiff(names(plain_total), c("n", "visits"))
  for (nm in shared) expect_equal(boot[[nm]], plain_total[[nm]], info = nm)
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
  expect_true(all(c("n_mean", "n_sd", "n_lower", "n_upper",
                    "tte_mean", "tte_sd", "tte_lower", "tte_upper",
                    "tte_ci_type", "ci_type", "n_failed") %in% names(out)))
})

test_that("a bootstrapped `n` interval is always widened to even sizes", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 8, seed = 5,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  expect_true(all(out$n_lower %% 2 == 0))
  expect_true(all(out$n_upper %% 2 == 0))
})

test_that("`type = \"percentile\"` is honoured cell by cell, without a jackknife", {
  pars <- small_fit()
  out <- slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3), effectiveness = 0.33,
                                     R = 6, seed = 1, type = "percentile")
  expect_true(all(out$ci_type == "percentile"))
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
                             type = "bca", probs = c(0.025, 0.975), context = "test",
                             what = "", lattice = TRUE)
  expect_true(res$starved)
  expect_true(is.na(res$mean))
  expect_true(all(is.na(res$ci)))
  expect_equal(res$n_failed, 3L)
})

test_that("too few refits succeeding aborts the whole grid, as it does for slope_bootstrap()", {
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

# The frames and nothing else: everything between the class header on the first
# line and the first note. "Counts:" -- which basis n/visits are on -- is now
# the first note printed, ahead of "Bootstrap:", so it is what this cuts at.
frame_block <- function(out) out[seq(2L, min(grep("Counts:", out, fixed = TRUE)) - 1L)]

test_that("print.slope_sample_size_grid_boot() prints the grid's own columns plus three, per arm by default", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))
  expect_identical(attr(out, "per_arm"), TRUE)

  lines <- capture.output(print(out))
  expect_true(any(grepl("<slope_sample_size_grid_boot>", lines, fixed = TRUE)))

  # Nothing left of the pipe table: this prints as the plain grid does, with
  # three columns added -- reduced, like the plain grid itself, to one basis
  # (basis_columns(), grid.R). R wraps a wide frame into blocks, each block
  # headed by its own column names and followed by rows that begin with a row
  # number, so the printed column set is read off the lines that are not rows.
  expect_false(any(startsWith(lines, "  |")))
  headers <- frame_block(lines)[!grepl("^[0-9]", frame_block(lines))]
  printed <- Filter(nzchar, unlist(strsplit(trimws(headers), " +")))
  added <- attr(out, "added_cols")
  raw_cells <- out[, setdiff(names(out), added), drop = FALSE]
  arm_names <- names(slopepower:::basis_columns(raw_cells, TRUE))
  expect_equal(printed, c(arm_names, "n_mean", "n_sd", "n_ci"))
  expect_false("n_per_arm" %in% printed)
  expect_false("visits_total" %in% printed)

  # The three added columns are the object's own values, halved to the
  # per-arm basis being shown here -- not a re-derivation.
  expect_true(any(grepl(format(out$n_mean[1L] / 2, digits = 7), lines, fixed = TRUE)))

  # "Counts:" states the basis, ahead of the bootstrap method note that used
  # to ride in the span header of the old hand-drawn table.
  expect_true(any(grepl("Counts:.*per arm", lines)))
  expect_true(any(grepl("Bootstrap:   R = 5 replicates; 95% BCa intervals.",
                        lines, fixed = TRUE)))

  # A clean run (no refit failures, no starved cells) still reports both
  # notes, per the convention set for slope_bootstrap()'s own print method.
  expect_true(any(grepl("bootstrap replicates failed to refit", lines)))
  expect_true(any(grepl("of replicates refit a slope on the opposite side", lines)))
})

test_that("print.slope_sample_size_grid_boot() shows the trial total under per_arm = FALSE, exactly double the per-arm print", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  lines <- capture.output(print(out, per_arm = FALSE))
  headers <- frame_block(lines)[!grepl("^[0-9]", frame_block(lines))]
  printed <- Filter(nzchar, unlist(strsplit(trimws(headers), " +")))
  added <- attr(out, "added_cols")
  raw_cells <- out[, setdiff(names(out), added), drop = FALSE]
  total_names <- names(slopepower:::basis_columns(raw_cells, FALSE))
  expect_equal(printed, c(total_names, "n_mean", "n_sd", "n_ci"))
  expect_true(any(grepl("Counts:.*total", lines)))

  # The object itself is untouched by which basis was printed, and the two
  # bases are exact halves of one another -- endpoint for endpoint, not only
  # in the point estimate.
  expect_identical(attr(out, "per_arm"), TRUE)
  expect_equal(out$n, 2L * out$n_per_arm)
  expect_equal(out$visits_total, 2L * out$visits_per_arm)
  expect_equal(out$n_mean, 2L * (out$n_mean / 2L))
})

test_that("print(x, per_arm = ...) overrides the basis a grid was built with", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 6,
    per_arm = FALSE))
  expect_identical(attr(out, "per_arm"), FALSE)

  default_lines  <- capture.output(print(out))
  override_lines <- capture.output(print(out, per_arm = TRUE))
  expect_true(any(grepl("Counts:.*total", default_lines)))
  expect_true(any(grepl("Counts:.*per arm", override_lines)))
  expect_false(identical(default_lines, override_lines))

  expect_error(print(out, per_arm = NA), "per_arm")
})

test_that("print.slope_sample_size_grid_boot() no longer reports the resampled slope", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits = list(a = c(0, 1, 2, 3), b = c(0, 1, 2))))

  frame <- frame_block(capture.output(print(out)))
  expect_false(any(grepl("slope", frame, ignore.case = TRUE)))
  expect_false(any(grepl(format(attr(out, "slope_mean"), digits = 6), frame, fixed = TRUE)))

  # Still on the object, and the straddle note still reports the one thing
  # about it that bears on whether these intervals mean anything.
  expect_equal(attr(out, "slope_observed"), pars$slope)
  expect_length(attr(out, "slope_ci"), 2L)
})

test_that("print.slope_sample_size_grid_boot() marks a starved cell and a fallen-back interval", {
  # Both markers used to live in the cells of the hand-drawn table; they still
  # do, in the one character column that holds the interval. Forced onto a
  # finished grid rather than provoked, so the test does not depend on which
  # resample happens to starve a cell.
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits = list(a = c(0, 1, 2, 3), b = c(0, 1, 2))))
  out$n_lower[1L] <- NA_real_
  out$n_upper[1L] <- NA_real_
  out$ci_type[1L] <- NA_character_
  out$ci_type[2L] <- "percentile"

  # R wraps a wide frame into blocks; the interval is the last column of the
  # last of them, so its rows are the two lines under the header that ends in
  # "n_ci", and the cell is whatever follows the last run of padding spaces.
  lines <- capture.output(print(out))
  frame <- frame_block(lines)
  header <- grep("n_ci$", frame)
  expect_length(header, 1L)
  ci <- sub(".*  +", "", frame[header + 1:2])
  expect_equal(ci[1L], "--")
  expect_match(ci[2L], "[*]$")

  # And each marker is explained beneath, in a note that appears only when the
  # thing it explains is on the page.
  expect_true(any(grepl("cells marked \"--\" had fewer than two surviving", lines)))
  expect_true(any(grepl("* percentile interval; BCa could not be built there.",
                        lines, fixed = TRUE)))
})

test_that("print.slope_sample_size_grid_boot() marks a fallen-back tte interval too", {
  # `tte`'s interval has its own bias correction and its own jackknife column,
  # so it can fall back to percentile where `n`'s did not. Forced, as the
  # markers above are, rather than provoked.
  pars <- small_fit()
  out <- slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3),
                                     effectiveness = c(0.25, 0.33), R = 5, seed = 6)
  expect_true(all(out$ci_type == "bca"))
  out$tte_ci_type[1L] <- "percentile"
  out$tte_lower[2L] <- NA_real_
  out$tte_upper[2L] <- NA_real_
  out$tte_ci_type[2L] <- NA_character_

  lines <- capture.output(print(out))
  frame <- frame_block(lines)
  header <- grep("tte_ci$", frame)
  expect_length(header, 1L)
  ci <- sub(".*  +", "", frame[header + 1:2])
  expect_match(ci[1L], "[*]$")
  expect_equal(ci[2L], "--")

  # Both markers are explained beneath even though it was `tte`, not `n`,
  # that earned them.
  expect_true(any(grepl("* percentile interval; BCa could not be built there.",
                        lines, fixed = TRUE)))
  expect_true(any(grepl("cells marked \"--\" had fewer than two surviving", lines)))
})

test_that("the extra-failure note reports the losses beyond the refits, not the totals", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 10, seed = 6,
    visits = list(a = c(0, 1, 2, 3), b = c(0, 1, 2))))

  # Three replicates lost to the refit, shared by every cell; one cell lost
  # four more to its own stage-two solve, the other none. The note is about
  # the "more", so it must read 0-4 and not the 3-7 totals.
  out <- structure(out, n_refit_failed = 3L)
  out$n_failed <- c(3L, 7L)
  # boot_note() wraps a note across lines, so the notes are read as one string
  # rather than line by line: the phrase this is about spans a break.
  notes <- function(x) gsub(" +", " ", paste(capture.output(print(x)), collapse = " "))
  expect_match(notes(out),
               "beyond the refits above (0-4 of 10 failed per cell).", fixed = TRUE)

  # And a grid that lost nothing beyond the refits does not raise the note.
  out$n_failed <- c(3L, 3L)
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
  for (col in c("tte_mean", "tte_sd", "tte_ci")) {
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

  frame <- slopepower:::grid_boot_tte_frame(out, "bca")
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
  expect_true(any(grepl("<slope_sample_size_grid_boot>", lines, fixed = TRUE)))
  # No second frame: with one target effect for the whole grid there is
  # nothing for it to show that the notes do not.
  expect_false(any(grepl("tte_mean", lines, fixed = TRUE)))
  # The notes beneath still print, which they cannot if the frames aborted.
  expect_true(any(grepl("bootstrap replicates failed to refit", lines)))
})

test_that("print.slope_sample_size_grid_boot() falls back to a plain print for a missing summary", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2))))

  lines <- capture.output(print(structure(out, R = NULL, class = class(out))))
  expect_false(any(grepl("<slope_sample_size_grid_boot>", lines, fixed = TRUE)))
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
  expect_equal(row_sub$n, out$n[1L])

  col_sub <- out[, c("design", "n")]
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
                                dropout = dropout_rate(0.3), R = 4, seed = 1),
    "baseline visit only"
  )
})

test_that("slope_bootstrap() points at slope_sample_size_grid_boot() for a grid result", {
  pars <- small_fit()
  g <- slope_sample_size_grid(pars, visits = c(0, 1, 2), effectiveness = 0.33)
  expect_error(slope_bootstrap(g), "slope_sample_size_grid_boot")
})
