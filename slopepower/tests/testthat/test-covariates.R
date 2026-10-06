# Covariate adjustment in slope_params().

cov_data <- function(seed = 42) {
  set.seed(seed)
  d <- slpower1
  ids <- unique(d$id)
  age <- stats::setNames(stats::rnorm(length(ids), 45, 10), ids)
  sex <- stats::setNames(factor(sample(c("F", "M"), length(ids), TRUE)), ids)
  d$age <- age[as.character(d$id)]
  d$sex <- sex[as.character(d$id)]
  # age is prognostic for slope, so adjusting for it should shrink sigma2_slope
  d$sdmt <- d$sdmt - 0.08 * (d$age - 45) * d$visit
  d
}

test_that("adjusted fit matches a hand-built lme with centred covariates", {
  d <- cov_data()
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + sex))
  d$t <- d$visit - stats::ave(d$visit, d$id, FUN = min)
  first <- !duplicated(d$id)
  d$age_c <- d$age - mean(d$age[first])
  d$m     <- (d$sex == "M") - mean(d$sex[first] == "M")
  h <- nlme::lme(sdmt ~ t + age_c + m + age_c:t + m:t, random = ~ t | id,
                 data = d, method = "REML", control = slope_lme_control())
  G <- nlme::getVarCov(h)
  expect_equal(p$slope, nlme::fixef(h)[["t"]], tolerance = 1e-6)
  expect_equal(p$sigma2_slope, G["t", "t"], tolerance = 1e-4)
  expect_equal(p$sigma2_intercept, G["(Intercept)", "(Intercept)"], tolerance = 1e-4)
  expect_equal(p$sigma2_residual, stats::sigma(h)^2, tolerance = 1e-4)
})

test_that("a prognostic covariate-by-time term shrinks sigma2_slope; intercept-only does not", {
  d <- cov_data()
  u  <- suppressMessages(slope_params(sdmt ~ visit | id, d))
  a  <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age))
  a0 <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age,
                                      covariate_time = FALSE))
  expect_lt(a$sigma2_slope, 0.8 * u$sigma2_slope)
  expect_equal(a0$sigma2_slope, u$sigma2_slope, tolerance = 1e-3)
})

test_that("no covariates leaves the fit unchanged", {
  d <- cov_data()
  u <- suppressMessages(slope_params(sdmt ~ visit | id, d))
  n <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = NULL))
  expect_identical(u[c("slope", "sigma2_slope", "sigma2_residual")],
                   n[c("slope", "sigma2_slope", "sigma2_residual")])
  expect_identical(names(nlme::fixef(u$fit)), c("(Intercept)", "sp_time"))
})

test_that("covariates work with treated and healthy comparators", {
  set.seed(3)
  t3 <- slpower3
  t3$x <- stats::ave(t3$id, t3$id, FUN = function(i) stats::rnorm(1))
  p <- slope_params(sdmt ~ visit | id, t3, treated = treat, covariates = ~ x)
  expect_true(all(c("sp_cov_1", "sp_time:sp_cov_1") %in% names(nlme::fixef(p$fit))))
  expect_false(is.na(p$slope_comparator))

  h2 <- slpower2
  h2$x <- stats::ave(h2$id, h2$id, FUN = function(i) stats::rnorm(1))
  q <- suppressMessages(slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, h2,
                                     healthy = case, covariates = ~ x))
  expect_equal(q$comparator, "healthy")
  expect_false(is.na(q$slope_comparator))
})

test_that("bootstrap refits keep the covariate adjustment", {
  d <- cov_data()
  a <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + sex))
  fr <- boot_frame(a, "test")
  expect_true(all(c("sp_cov_1", "sp_cov_2") %in% names(fr)))
  refit <- make_refitter(a)(fr)
  expect_equal(refit$sigma2_slope, a$sigma2_slope, tolerance = 1e-6)
  expect_equal(names(nlme::fixef(refit$fit)), names(nlme::fixef(a$fit)))
})

test_that("bad covariate specifications are refused", {
  d <- cov_data()
  d$tv <- d$visit
  expect_error(suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ tv)),
               "constant within a participant")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = "age"),
               "one-sided formula")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = ~ nope),
               "could not evaluate `covariates`")
})

test_that("missing covariate values drop rows through na.action", {
  d <- cov_data()
  d$age[1] <- NA
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age))
  expect_equal(p$n_obs, nrow(d) - 1L)
})

test_that("model.matrix-style covariate formulas are expanded as documented", {
  d <- cov_data()
  labs <- function(f) unname(attr(covariate_matrix(f, d, "test"), "labels"))
  expect_equal(labs(~ age * sex), c("age", "sexM", "age:sexM"))
  expect_equal(labs(~ 0 + sex), "sexM")      # always coded against an intercept
  expect_equal(labs(~ sex - 1), "sexM")
  expect_length(labs(~ splines::ns(age, 3)), 3L)
  # poly() is computed through a QR, so equal ages differ in the last bit;
  # that must not count as varying within a participant
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ poly(age, 2)))
  expect_s3_class(p, "slope_params")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = ~ .), "cannot use `.`")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = ~ age + offset(age)),
               "offset")
  expect_error(suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + visit)),
               "`visit` changes during follow-up")
})
