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
  pars <- small_fit()
  plain <- slope_sample_size_grid(pars, visits = c(0, 1, 2, 3), effectiveness = 0.33)
  boot <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 1))

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
  expect_true(all(c("n_mean", "n_sd", "n_lower", "n_upper",
                    "tte_mean", "tte_sd", "tte_lower", "tte_upper",
                    "ci_type", "n_failed") %in% names(out)))
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

test_that("print.slope_sample_size_grid_boot() prints a width-consistent table", {
  pars <- small_fit()
  out <- suppressWarnings(slope_sample_size_grid_boot(
    pars, power = 0.8, effectiveness = 0.33, R = 5, seed = 6,
    visits  = list(a = c(0, 1, 2, 3), b = c(0, 1, 2)),
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05))))

  lines <- capture.output(print(out))
  expect_true(any(grepl("<slope_sample_size_grid_boot>", lines, fixed = TRUE)))
  expect_true(any(grepl("Slope", lines, fixed = TRUE)))
  expect_true(any(grepl("a / none", lines, fixed = TRUE)))

  # Every row of the table itself starts "  |" and, being built by the same
  # boot_table_row() every column in it, is one consistent width -- the same
  # invariant test-bootstrap.R pins for boot_n_table().
  rows <- lines[startsWith(lines, "  |")]
  expect_true(length(rows) > 1L)
  expect_equal(length(unique(nchar(rows))), 1L)

  # A clean run (no refit failures, no starved cells) still reports both
  # notes, per the convention set for slope_bootstrap()'s own print method.
  expect_true(any(grepl("bootstrap replicates failed to refit", lines)))
  expect_true(any(grepl("of replicates refit a slope on the opposite side", lines)))
})

test_that("print.slope_sample_size_grid_boot() adds a target-effect table only when effectiveness varies", {
  pars <- small_fit()
  fixed <- suppressWarnings(slope_sample_size_grid_boot(
    pars, visits = c(0, 1, 2, 3), effectiveness = 0.33, R = 5, seed = 6))
  varying <- slope_sample_size_grid_boot(pars, visits = c(0, 1, 2, 3),
                                         effectiveness = c(0.25, 0.33), R = 5, seed = 6)

  expect_false(any(grepl("Target treatment effect", capture.output(print(fixed)))))
  expect_true(any(grepl("Target treatment effect", capture.output(print(varying)))))
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
  expect_null(attr(row_sub, "row_labels"))
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
