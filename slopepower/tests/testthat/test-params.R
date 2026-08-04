# Layer 1 --- slope_params() and slope_params_manual(). See CONTRACT.md sec 2.

# Mixed-model fits are the slow part of this suite. Cache them in the global
# environment so that test-params.R and test-paper-parity.R together fit each
# dataset once. (Defined in both files rather than a shared helper because this
# agent owns only test-*.R.)
# --- slope_params_manual() --------------------------------------------------

test_that("slope_params_manual() returns exactly the contract fields", {
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                           sigma2_slope = 2, sigma_cov = 5,
                           sigma2_residual = 10)
  expect_s3_class(p, "slope_params")
  expect_setequal(names(p), c(
    "slope", "slope_comparator", "comparator", "sigma2_intercept",
    "sigma2_slope", "sigma_cov", "sigma2_residual", "n_obs", "n_subjects",
    "common_variance", "time_shifted", "fit", "call"))
  expect_length(names(p), 13L)
})

test_that("manually supplied parameters are stored unchanged", {
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                           sigma2_slope = 2, sigma_cov = 5,
                           sigma2_residual = 10)
  expect_equal(p$slope, -1.672)
  expect_equal(p$sigma2_intercept, 100)
  expect_equal(p$sigma2_slope, 2)
  expect_equal(p$sigma_cov, 5)
  expect_equal(p$sigma2_residual, 10)
  expect_identical(p$comparator, "none")
  expect_true(is.na(p$slope_comparator))
  expect_true(is.na(p$n_obs))
  expect_true(is.na(p$n_subjects))
  expect_null(p$fit)
})

test_that("slope_params_manual() rejects non-positive variance components", {
  base <- list(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
               sigma_cov = 5, sigma2_residual = 10)
  for (nm in c("sigma2_intercept", "sigma2_slope", "sigma2_residual")) {
    args <- base; args[[nm]] <- 0
    expect_error(do.call(slope_params_manual, args), nm)
    args[[nm]] <- -1
    expect_error(do.call(slope_params_manual, args), nm)
  }
})

test_that("slope_params_manual() rejects a non-positive-definite covariance", {
  # var_int * var_slope < cov^2
  expect_error(
    slope_params_manual(slope = -1, sigma2_intercept = 1, sigma2_slope = 1,
                        sigma_cov = 5, sigma2_residual = 10),
    "positive definite"
  )
})

test_that("a comparator slope is required unless comparator is 'none'", {
  for (cmp in c("healthy", "treated")) {
    expect_error(
      slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                          sigma_cov = 5, sigma2_residual = 10, comparator = cmp),
      "slope_comparator"
    )
  }
  expect_no_error(
    slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                        sigma_cov = 5, sigma2_residual = 10,
                        slope_comparator = 0.5, comparator = "healthy")
  )
})

test_that("a comparator slope is discarded when comparator is 'none'", {
  p <- slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                           sigma_cov = 5, sigma2_residual = 10,
                           slope_comparator = 0.5, comparator = "none")
  expect_true(is.na(p$slope_comparator))
})

test_that("slope_params_manual() rejects an unknown comparator", {
  expect_error(
    slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                        sigma_cov = 5, sigma2_residual = 10,
                        slope_comparator = 0.5, comparator = "nonsense")
  )
})

# --- slope_params(): scenario dispatch --------------------------------------

test_that("the scenario follows which group argument is supplied", {
  expect_identical(paper_fit("slpower1")$comparator, "none")
  expect_identical(paper_fit("slpower2")$comparator, "healthy")
  expect_identical(paper_fit("slpower3")$comparator, "treated")
})

test_that("healthy and treated are mutually exclusive", {
  d <- load_paper_data("slpower2")
  d$treat <- d$case
  expect_error(slope_params(sdmt ~ time | id, d, healthy = case, treated = treat),
               "only one of")
})

test_that("the formula must name a subject identifier", {
  d <- load_paper_data("slpower1")
  expect_error(slope_params(sdmt ~ time, d), "subject identifier")
})

test_that("a group variable must be binary and 0/1 coded", {
  d <- load_paper_data("slpower2")
  d$bad <- d$case + 1                      # coded 1/2
  expect_error(slope_params(sdmt ~ time | id, d, healthy = bad), "0/1")
  d$constant <- 1L
  expect_error(slope_params(sdmt ~ time | id, d, healthy = constant),
               "two distinct values")
})

test_that("slope_params() reports the number of observations and subjects", {
  expect_equal(paper_fit("slpower1")$n_obs, 800L)
  expect_equal(paper_fit("slpower1")$n_subjects, 200L)
  expect_equal(paper_fit("slpower2")$n_obs, 2000L)
  expect_equal(paper_fit("slpower2")$n_subjects, 500L)
  expect_equal(paper_fit("slpower3")$n_obs, 450L)
  expect_equal(paper_fit("slpower3")$n_subjects, 150L)
})

test_that("fitted variance components are positive and positive definite", {
  for (nm in c("slpower1", "slpower2", "slpower3")) {
    p <- paper_fit(nm)
    expect_gt(p$sigma2_intercept, 0)
    expect_gt(p$sigma2_slope, 0)
    expect_gt(p$sigma2_residual, 0)
    G <- matrix(c(p$sigma2_intercept, p$sigma_cov,
                  p$sigma_cov, p$sigma2_slope), 2, 2)
    expect_gt(det(G), 0)
  }
})

test_that("slpower1 variance components land near the simulation truth", {
  # Paper appendix: sigma2_a = 100, sigma2_b = 2, sigma_ab = 5, sigma2_e = 10.
  # Sampling noise over 200 subjects, so these are loose bounds.
  p <- paper_fit("slpower1")
  expect_gt(p$sigma2_intercept, 70);  expect_lt(p$sigma2_intercept, 140)
  expect_gt(p$sigma2_slope, 1);       expect_lt(p$sigma2_slope, 4)
  expect_gt(p$sigma2_residual, 6);    expect_lt(p$sigma2_residual, 14)
})

# --- per-subject time re-origining ------------------------------------------

test_that("time is shifted so each subject starts at zero", {
  # slpower2 records calendar visit dates, so every subject starts elsewhere.
  expect_true(paper_fit("slpower2")$time_shifted)
  # slpower1 and slpower3 already start at 0.
  expect_false(paper_fit("slpower1")$time_shifted)
  expect_false(paper_fit("slpower3")$time_shifted)
})

test_that("the re-origining emits a message, as Stata warns", {
  d <- load_paper_data("slpower2")
  expect_message(slope_params(sdmt ~ time | id, d, healthy = case),
                 "shifted")
})

test_that('origin = "none" leaves time alone and changes the fit', {
  d <- load_paper_data("slpower2")
  shifted   <- paper_fit("slpower2")
  unshifted <- suppressMessages(
    slope_params(sdmt ~ time | id, d, healthy = case, origin = "none"))
  expect_false(unshifted$time_shifted)
  # the slope is barely affected, but the intercept variance is not comparable
  expect_equal(unshifted$slope, shifted$slope, tolerance = 0.05)
  expect_false(isTRUE(all.equal(unshifted$sigma2_intercept,
                                shifted$sigma2_intercept)))
})

test_that("an explicit unit change rescales the slope proportionally", {
  d <- load_paper_data("slpower1")
  yearly <- paper_fit("slpower1")
  d$time <- d$time / 0.5                  # half-year axis
  halved <- slope_params(sdmt ~ time | id, d)
  expect_equal(halved$slope, yearly$slope / 2, tolerance = 1e-6)
})

# --- common_variance --------------------------------------------------------

test_that("common_variance is ignored with a warning outside the healthy case", {
  d <- load_paper_data("slpower1")
  expect_warning(slope_params(sdmt ~ time | id, d, common_variance = TRUE),
                 "only when `healthy`")
})

test_that("common_variance does not change the case estimates", {
  # The healthy model factorises into two independent per-group fits, so
  # reducing the controls' random-effects block cannot move the cases' numbers.
  d <- load_paper_data("slpower2")
  full    <- paper_fit("slpower2")
  reduced <- suppressMessages(
    slope_params(sdmt ~ time | id, d, healthy = case, common_variance = TRUE))
  expect_true(reduced$common_variance)
  expect_false(full$common_variance)
  expect_equal(reduced$slope, full$slope, tolerance = 1e-5)
  expect_equal(reduced$sigma2_intercept, full$sigma2_intercept, tolerance = 1e-4)
  expect_equal(reduced$sigma2_slope, full$sigma2_slope, tolerance = 1e-4)
  expect_equal(reduced$sigma2_residual, full$sigma2_residual, tolerance = 1e-4)
})

# --- printing ---------------------------------------------------------------

test_that("print.slope_params() runs for each scenario", {
  expect_output(print(paper_fit("slpower1")), "single group")
  expect_output(print(paper_fit("slpower2")), "healthy controls")
  expect_output(print(paper_fit("slpower3")), "previous randomised trial")
  out <- capture.output(print(paper_fit("slpower2")))
  expect_true(any(grepl("observed difference in slopes", out)))
})

test_that("print.slope_params() marks manually supplied parameters", {
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                           sigma2_slope = 2, sigma_cov = 5,
                           sigma2_residual = 10)
  expect_output(print(p), "supplied directly")
})
