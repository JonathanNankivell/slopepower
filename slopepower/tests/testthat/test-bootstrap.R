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
  skip_without_paper_data()
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 15,
                                        statistic = "slope", seed = 1))
  expect_true(is.list(b))
  reps <- b$replicates %||% b$t
  expect_true(!is.null(reps))
  expect_gt(stats::sd(reps, na.rm = TRUE), 0)
  # the replicate slopes should surround the observed one
  expect_lt(abs(mean(reps, na.rm = TRUE) - paper_fit("slpower1")$slope), 0.5)
})

test_that("slope_bootstrap() is reproducible under a fixed seed", {
  skip_without_paper_data()
  a <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 10,
                                        statistic = "slope", seed = 42))
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 10,
                                        statistic = "slope", seed = 42))
  expect_equal(a$replicates %||% a$t, b$replicates %||% b$t)
})

test_that("a sample-size statistic requires a trial design", {
  skip_without_paper_data()
  expect_error(slope_bootstrap(paper_fit("slpower1"), R = 5, statistic = "n"),
               "design")
})

test_that("slope_bootstrap() can bootstrap the sample size itself", {
  skip_without_paper_data()
  b <- suppressWarnings(slope_bootstrap(
    paper_fit("slpower1"), R = 15, statistic = "n", seed = 7,
    design = trial_design(c(0, 1, 2)), effectiveness = 0.33))
  reps <- b$replicates %||% b$t
  expect_true(all(reps > 0, na.rm = TRUE))
  expect_gt(stats::sd(reps, na.rm = TRUE), 0)
})

test_that("slope_bootstrap() validates R and level", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")
  expect_error(slope_bootstrap(p, R = 0, statistic = "slope"), "R")
  expect_error(slope_bootstrap(p, R = 5, statistic = "slope", level = 1), "level")
  expect_error(slope_bootstrap(p, R = 5, statistic = "slope", level = 0), "level")
})

test_that("slope_bootstrap() warns when the slope is weak relative to its error", {
  # Paper section 2.6: if |slope| / se(slope) < 2.5 the replicate slopes can
  # straddle zero and the resulting interval is meaningless.
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  d$sdmt <- d$sdmt + 1.66 * d$time         # flatten the slope, keep the structure
  flat <- slope_params(sdmt ~ time | id, d)
  expect_warning(slope_bootstrap(flat, R = 5, statistic = "slope", seed = 1))
})
