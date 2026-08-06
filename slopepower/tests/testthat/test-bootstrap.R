# Layer 4 --- slope_bootstrap().
#
# Replaces the Stata incantation
#   bootstrap r(sampsize), cluster(id) idcluster(id2) strata(case) bca
#             jack(n(r(obs_in_model))): slopepower ...
# with resampling that is correct by construction: subjects, not rows, are the
# unit of resampling.
#
# R is kept small throughout; these test the machinery, not the coverage.

test_that("slope_bootstrap() returns replicates with a spread, for the slope", {
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 15, seed = 1))
  expect_true(is.list(b))
  expect_equal(b$statistic, "slope")
  reps <- b$replicates %||% b$t
  expect_true(!is.null(reps))
  expect_gt(stats::sd(reps, na.rm = TRUE), 0)
  # the replicate slopes should surround the observed one
  expect_lt(abs(mean(reps, na.rm = TRUE) - paper_fit("slpower1")$slope), 0.5)
})

test_that("slope_bootstrap() is reproducible under a fixed seed", {
  a <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 10, seed = 42))
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 10, seed = 42))
  expect_equal(a$replicates %||% a$t, b$replicates %||% b$t)
})

test_that("slope_bootstrap() dispatches on what it is handed", {
  p <- paper_fit("slpower1")
  ss <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  pw <- slope_power(p, c(0, 1, 2), n = 712, effectiveness = 0.33)

  expect_equal(suppressWarnings(slope_bootstrap(ss, R = 3, seed = 1))$statistic, "n")
  expect_equal(suppressWarnings(slope_bootstrap(pw, R = 3, seed = 1))$statistic, "power")
  expect_equal(suppressWarnings(slope_bootstrap(p, R = 3, seed = 1))$statistic, "slope")

  # anything else is refused with a pointer to what is accepted
  expect_error(slope_bootstrap(data.frame(x = 1)), "cannot bootstrap an object")
  expect_error(slope_bootstrap(42), "slope_sample_size")
})

test_that("slope_bootstrap() can bootstrap the sample size itself", {
  ss <- slope_sample_size(paper_fit("slpower1"), trial_design(c(0, 1, 2)),
                          effectiveness = 0.33)
  b <- suppressWarnings(slope_bootstrap(ss, R = 15, seed = 7))
  reps <- b$replicates %||% b$t
  expect_true(all(reps > 0, na.rm = TRUE))
  expect_gt(stats::sd(reps, na.rm = TRUE), 0)
  expect_equal(b$observed, ss$n)
})

test_that("slope_bootstrap() validates R and level", {
  p <- paper_fit("slpower1")
  expect_error(slope_bootstrap(p, R = 0), "R")
  expect_error(slope_bootstrap(p, R = 5, level = 1), "level")
  expect_error(slope_bootstrap(p, R = 5, level = 0), "level")
})

test_that("slope_bootstrap() warns when the slope is weak relative to its error", {
  # Paper section 2.6: if |slope| / se(slope) < 2.5 the replicate slopes can
  # straddle zero and the resulting interval is meaningless.
  d <- load_paper_data("slpower1")
  d$sdmt <- d$sdmt + 1.66 * d$time         # flatten the slope, keep the structure
  flat <- slope_params(sdmt ~ time | id, d)
  expect_warning(slope_bootstrap(flat, R = 5, seed = 1))
})

test_that("slope_bootstrap() refuses parameters that carry no data", {
  manual <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                                sigma2_slope = 2, sigma_cov = 5,
                                sigma2_residual = 10)
  expect_error(slope_bootstrap(manual, R = 5), "no fitted model")
  expect_error(
    slope_bootstrap(slope_sample_size(manual, c(0, 1, 2), effectiveness = 0.33), R = 5),
    "no fitted model")
})
