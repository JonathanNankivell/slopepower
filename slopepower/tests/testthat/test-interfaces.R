# The interface pieces that make the functions agree with each other and with
# nlme: slope_var()'s dropout, slope_params()'s `subset` and `control`, the
# accessor methods, and confint() on a bootstrap.

# --- slope_var() with dropout --------------------------------------------------

test_that("slope_var() with dropout is the var_tte slope_power() reports", {
  p <- paper_fit("slpower1")
  v <- c(0, 1, 2, 3)
  d <- c(0, 0.1, 0.1)
  expect_equal(slope_var(p, v, dropout = d),
               slope_power(p, v, dropout = d, n = 300)$var_tte)
  expect_equal(slope_var(p, v, dropout = cumsum(d), dropout_scale = "cumulative"),
               slope_var(p, v, dropout = d))
  expect_gt(slope_var(p, v, dropout = d), slope_var(p, v))
})

test_that("slope_var() with no dropout is the same however it is written", {
  p <- paper_fit("slpower1")
  expect_equal(slope_var(p, c(0, 1, 2), dropout = c(0, 0)), slope_var(p, c(0, 1, 2)))
  # Without dropout any increasing times still do, as they always have.
  expect_true(is.finite(slope_var(p, c(1, 2, 3))))
})

test_that("slope_var() with dropout validates the design as slope_power() does", {
  p <- paper_fit("slpower1")
  expect_error(slope_var(p, c(1, 2, 3), dropout = c(0, 0.1)), "baseline visit at time 0")
  expect_error(suppressWarnings(slope_var(p, c(0, 1, 2), dropout = c(1, 0))),
               "nothing to estimate the slope")
})

# --- slope_params(): subset and control ----------------------------------------

test_that("`subset` fits the same model as filtering the data first", {
  a <- slope_params(sdmt ~ visit | id, data = slpower1, subset = id <= 100)
  b <- slope_params(sdmt ~ visit | id, data = slpower1[slpower1$id <= 100, ])
  expect_equal(a$slope, b$slope)
  expect_equal(a$sigma2_slope, b$sigma2_slope)
  expect_equal(a$n_subjects, 100L)
  # Row numbers select the same rows; NA in a logical subset counts as FALSE.
  rows <- which(slpower1$id <= 100)
  expect_equal(slope_params(sdmt ~ visit | id, data = slpower1, subset = rows)$slope, a$slope)
  keep <- ifelse(slpower1$id <= 100, TRUE, NA)
  expect_equal(slope_params(sdmt ~ visit | id, data = slpower1, subset = keep)$slope, a$slope)
})

test_that("a malformed `subset` is refused", {
  expect_error(slope_params(sdmt ~ visit | id, data = slpower1, subset = c(TRUE, FALSE)),
               "one value per row")
  expect_error(slope_params(sdmt ~ visit | id, data = slpower1, subset = "a"),
               "distinct row numbers")
  for (bad in list(c(-1, 2), 0, c(1, 1, 2)))
    expect_error(slope_params(sdmt ~ visit | id, data = slpower1, subset = bad),
                 "distinct row numbers")
})

test_that("`control` defaults to slope_lme_control() and is recorded", {
  p <- paper_fit("slpower1")
  expect_identical(p$control, slope_lme_control())
  expect_null(slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                                  cov_intercept_slope = 5, sigma2_residual = 10)$control)
})

test_that("slope_lme_control() overrides some settings and keeps the rest", {
  ctrl <- slope_lme_control(maxIter = 500)
  expect_equal(ctrl$maxIter, 500)
  expect_equal(ctrl[setdiff(names(ctrl), "maxIter")],
               slope_lme_control()[setdiff(names(ctrl), "maxIter")])
})

test_that("slope_lme_control() refuses settings it cannot pass on", {
  expect_error(slope_lme_control(500), "must be named")
  expect_error(slope_lme_control(maxiter = 500), "maxiter.? is not a setting")
})

test_that("`control` settings that change the model rather than the effort are refused", {
  expect_error(slope_params(sdmt ~ visit | id, data = slpower1,
                            control = slope_lme_control(returnObject = TRUE)),
               "returnObject` must be FALSE")
  expect_error(slope_params(sdmt ~ visit | id, data = slpower1,
                            control = slope_lme_control(sigma = 1)),
               "sigma` must be 0")
})

test_that("a partial `control` list is refused rather than filled from nlme's defaults", {
  expect_error(slope_params(sdmt ~ visit | id, data = slpower1, control = list(maxIter = 500)),
               "slope_lme_control\\(maxIter = 500\\)")
})

test_that("a bootstrap refits every replicate under the fit's own `control`", {
  p <- slope_params(sdmt ~ visit | id, data = slpower1,
                    control = slope_lme_control(maxIter = 300))
  cl <- environment(make_refitter(p))$cl
  expect_equal(cl$control$maxIter, 300)
  # An object from before `control` was recorded refits under the default.
  p$control <- NULL
  expect_null(environment(make_refitter(p))$cl$control)
})

# --- nlme and stats accessors ----------------------------------------------------

test_that("fixef(), getVarCov() and sigma() read the package's own components", {
  m <- slope_params_manual(slope = -1.5, sigma2_intercept = 100, sigma2_slope = 2,
                           cov_intercept_slope = 5, sigma2_residual = 10,
                           slope_comparator = 0.5, comparator = "healthy")
  expect_equal(fixef(m), c(slope = -1.5, slope_comparator = 0.5))
  expect_equal(getVarCov(m), matrix(c(100, 5, 5, 2), 2,
                                    dimnames = rep(list(c("(Intercept)", "time")), 2)))
  expect_equal(sigma(m), sqrt(10))
  expect_equal(fixef(paper_fit("slpower1")), c(slope = paper_fit("slpower1")$slope))
})

test_that("vcov() agrees with slope_se() under every comparator", {
  for (f in c("slpower1", "slpower2", "slpower3")) {
    p <- paper_fit(f)
    V <- vcov(p)
    expect_identical(dimnames(V), rep(list(names(fixef(p))), 2))
    expect_equal(sqrt(V["slope", "slope"]), slope_se(p))
    expect_true(isSymmetric(V))
  }
})

test_that("an object saved before the rename says how to update it", {
  p <- unclass(paper_fit("slpower1"))
  p$sigma_cov <- p$cov_intercept_slope
  p$cov_intercept_slope <- NULL
  class(p) <- "slope_params"
  expect_error(slope_var(p, c(0, 1, 2)), "params\\$cov_intercept_slope <- params\\$sigma_cov")
})

test_that("vcov() is NA without a fitted model", {
  m <- slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                           cov_intercept_slope = 5, sigma2_residual = 10)
  expect_true(is.na(vcov(m)))
})

# --- confint() on a bootstrap ------------------------------------------------------

test_that("confint() returns the intervals the bootstrap stored", {
  p <- paper_fit("slpower1")
  b <- suppressWarnings(slope_sample_size_boot(p, c(0, 1, 2), R = 20, seed = 1,
                                               ci_method = "percentile", per_arm = FALSE))
  ci <- confint(b)
  expect_identical(dimnames(ci), list(c("n", "slope"), c("2.5 %", "97.5 %")))
  expect_equal(unname(ci["n", ]), b$ci)
  expect_equal(unname(ci["slope", ]), b$slope_ci)
  expect_identical(rownames(confint(b, "slope")), "slope")
  expect_error(confint(b, "power"), "`parm` must be among")
  expect_error(confint(b, level = 0.9), "built at level = 0.95")

  # On the basis print() shows: per arm by default, the total on request.
  pa <- suppressWarnings(slope_sample_size_boot(p, c(0, 1, 2), R = 20, seed = 1,
                                                ci_method = "percentile", per_arm = TRUE))
  expect_identical(rownames(confint(pa)), c("n_per_arm", "slope"))
  expect_equal(unname(confint(pa, "n")[1, ]), pa$ci / 2)
  expect_equal(unname(confint(pa, "n", per_arm = FALSE)[1, ]), pa$ci)

  bp <- suppressWarnings(slope_params_boot(p, R = 10, seed = 1, ci_method = "percentile"))
  expect_identical(rownames(confint(bp)), "slope")
})
