# Layer 3 (continued) --- the bound over all designs. See CONTRACT.md section
# 4.3 for the result object and 5.7 for the mathematics.

# --- slope_var_floor() ------------------------------------------------------

test_that("slope_var_floor() is twice the conditional variance of the random slope", {
  p <- ref_params()
  expect_equal(slope_var_floor(p),
               2 * (p$sigma2_slope - p$sigma_cov^2 / p$sigma2_intercept))

  # Equivalently, twice the Schur complement of the random-effects covariance.
  G <- matrix(c(p$sigma2_intercept, p$sigma_cov, p$sigma_cov, p$sigma2_slope), 2L, 2L)
  expect_equal(slope_var_floor(p), 2 * (G[2, 2] - G[2, 1] * G[1, 2] / G[1, 1]))
})

test_that("slope_var_floor() is positive whenever params are valid", {
  # Positivity is the determinant of the random-effects covariance over
  # sigma2_intercept, which check_re_covariance() already forces to be positive.
  # Near-singular is the interesting case: the floor shrinks toward zero but
  # never reaches it.
  expect_gt(slope_var_floor(ref_params()), 0)
  near <- slope_params_manual(slope = -1, sigma2_intercept = 100,
                              sigma2_slope = 2, sigma_cov = 14.1,
                              sigma2_residual = 10)
  expect_gt(slope_var_floor(near), 0)
  expect_lt(slope_var_floor(near), 0.05)
})

test_that("slope_var() is strictly above the floor for every schedule", {
  p <- ref_params()
  floor_value <- slope_var_floor(p)

  set.seed(11)
  gaps <- vapply(seq_len(500), function(i) {
    v <- sort(c(0, runif(sample.int(6L, 1L), 0.01, 30)))
    slope_var(p, v) - floor_value
  }, numeric(1))
  expect_true(all(gaps > 0))

  # ... and approached as the schedule becomes long and dense.
  expect_lt(slope_var(p, seq(0, 50, length.out = 2001)) - floor_value, 1e-3)
})

test_that("lengthening a two-visit schedule stops short of the floor", {
  # The limit of a two-visit design as the gap grows is the same expression
  # with the residual variance added to the baseline variance: strictly larger
  # than the floor, which is why reaching it needs density and not just
  # duration. Documented in slope_var_floor()'s details.
  p <- ref_params()
  two_visit_limit <- 2 * (p$sigma2_slope -
                            p$sigma_cov^2 / (p$sigma2_intercept + p$sigma2_residual))

  expect_gt(two_visit_limit, slope_var_floor(p))
  expect_equal(slope_var(p, c(0, 1e5)), two_visit_limit, tolerance = 1e-4)
  expect_gt(slope_var(p, c(0, 1e5)), slope_var_floor(p))
})

test_that("slope_var_floor() rejects a non-params object", {
  expect_error(slope_var_floor(c(0, 1, 2)), "slope_var_floor\\(\\).*slope_params")
})

# --- slope_sample_size_floor(), params method -------------------------------

test_that("the floor bounds every design, and no design beats it", {
  p <- ref_params()
  flr <- slope_sample_size_floor(p, effectiveness = 0.33)

  set.seed(12)
  for (i in seq_len(100)) {
    v <- sort(c(0, runif(sample.int(5L, 1L), 0.05, 20)))
    expect_gte(slope_sample_size(p, v, effectiveness = 0.33)$n, flr$n)
  }

  # Dropout only ever raises the requirement, so the bound survives it.
  d <- trial_design(c(0, 1, 2, 3), dropout = c(0, 0.1, 0.2))
  expect_gte(slope_sample_size(p, d, effectiveness = 0.33)$n, flr$n)
})

test_that("the floor is the limit a long dense schedule actually reaches", {
  p <- ref_params()
  flr <- slope_sample_size_floor(p, effectiveness = 0.33)
  expect_identical(slope_sample_size(p, seq(0, 100, length.out = 1001),
                                     effectiveness = 0.33)$n,
                   flr$n)
})

test_that("the floor uses the same equation (6) as slope_sample_size()", {
  p <- ref_params()
  flr <- slope_sample_size_floor(p, effectiveness = 0.33, power = 0.9, alpha = 0.01)

  z_sum_sq <- (qnorm(1 - 0.01 / 2) + qnorm(0.9))^2
  expect_identical(flr$n_per_arm,
                   ceiling(z_sum_sq / (flr$effect_size * flr$effectiveness)^2))
  expect_identical(flr$n, 2 * flr$n_per_arm)
  expect_equal(flr$var_tte, slope_var_floor(p))
  expect_equal(flr$effect_size, flr$slope_difference / sqrt(flr$var_tte))
})

test_that("the floor scales as effectiveness^-2", {
  p <- ref_params()
  a <- slope_sample_size_floor(p, effectiveness = 0.5)
  b <- slope_sample_size_floor(p, effectiveness = 0.25)
  # Exact on the unrounded per-arm requirement; `ceiling()` costs the last
  # percent, as it does for slope_sample_size() itself.
  expect_equal(b$n / a$n, 4, tolerance = 0.02)
})

test_that("the floor resolves the reference slope the way stage two does", {
  healthy <- ref_params("healthy", slope = -1.715, slope_comparator = 0.975)
  flr <- slope_sample_size_floor(healthy, effectiveness = 0.33)
  expect_equal(flr$reference_slope, 0.975)
  expect_equal(flr$slope_difference, -1.715 - 0.975)
  expect_equal(flr$tte, -0.33 * (-1.715 - 0.975))

  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)
  obs <- slope_sample_size_floor(treated, target = "observed")
  expect_equal(obs$reference_slope, -1.104)
  expect_identical(obs$target, "observed")
  expect_true(is.na(obs$effectiveness))

  # target = "effectiveness" on trial data ignores the treated arm (5.3).
  eff <- slope_sample_size_floor(treated, effectiveness = 0.33)
  expect_equal(eff$reference_slope, 0)
})

test_that("the floor enforces the same guards as the stage-two entry points", {
  p <- ref_params()
  expect_error(slope_sample_size_floor(p, effectiveness = 0.33, power = 1.2), "power")
  expect_error(slope_sample_size_floor(p, effectiveness = 0.33, alpha = 0), "alpha")
  expect_error(slope_sample_size_floor(p, effectiveness = 0), "effectiveness")
  expect_error(slope_sample_size_floor(p, effectiveness = 0.33, target = "observed"),
               "only one of")
  expect_error(slope_sample_size_floor(ref_params("none"), target = "observed"),
               "comparator")
  expect_error(slope_sample_size_floor(ref_params("healthy", slope = 0.5,
                                                  slope_comparator = 0.5)),
               "slope difference is exactly zero")
})

test_that("the floor rejects arguments it has no use for", {
  p <- ref_params()
  # `design` is the one a reader of slope_sample_size() would reach for, and
  # silently ignoring it would make the result look design-specific.
  expect_error(slope_sample_size_floor(p, design = c(0, 1, 2)), "unused argument")
  expect_error(slope_sample_size_floor(p, n = 400), "unused argument")
})

test_that("slope_sample_size_floor() rejects an object it cannot use", {
  expect_error(slope_sample_size_floor(c(0, 1, 2)), "cannot compute a floor")
  expect_error(slope_sample_size_floor(trial_design(c(0, 1, 2))), "cannot compute a floor")
})

# --- slope_sample_size_floor(), result method -------------------------------

test_that("the result method reuses the settings of the call that produced it", {
  p <- ref_params()
  ss <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33,
                          power = 0.9, alpha = 0.01)
  from_result <- slope_sample_size_floor(ss)
  from_params <- slope_sample_size_floor(p, effectiveness = 0.33,
                                         power = 0.9, alpha = 0.01)
  expect_equal(from_result[setdiff(names(from_result), "params")],
               from_params[setdiff(names(from_params), "params")])
  expect_lte(from_result$n, ss$n)
})

test_that("the result method carries target = \"observed\" through", {
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)
  ss <- slope_sample_size(treated, c(0, 2, 3), target = "observed")
  flr <- slope_sample_size_floor(ss)
  expect_identical(flr$target, "observed")
  expect_equal(flr$tte, ss$tte)
  expect_lte(flr$n, ss$n)
})

test_that("a slope_power result contributes the power it achieves", {
  p <- ref_params()
  pw <- slope_power(p, c(0, 1, 2), n = 450, effectiveness = 0.33)
  flr <- slope_sample_size_floor(pw)
  expect_equal(flr$power, pw$power)
  # The smallest n that could reach that power, so no larger than the n that did.
  expect_lte(flr$n, pw$n)
})

test_that("a saturated power is refused rather than reported as infinite", {
  p <- ref_params()
  pw <- slope_power(p, c(0, 1, 2), n = 2e5, effectiveness = 0.33)
  expect_identical(pw$power, 1)
  expect_error(slope_sample_size_floor(pw), "power is 1 to within double precision")
})

test_that("the result method rejects settings that belong to the original call", {
  ss <- slope_sample_size(ref_params(), c(0, 1, 2), effectiveness = 0.33)
  expect_error(slope_sample_size_floor(ss, power = 0.9), "unused argument")
})

# --- the result object ------------------------------------------------------

test_that("the floor result has exactly the documented fields", {
  flr <- slope_sample_size_floor(ref_params(), effectiveness = 0.33)
  expect_identical(names(flr),
                   c("n", "n_per_arm", "power", "alpha", "effectiveness", "target",
                     "tte", "var_tte", "effect_size", "slope_difference",
                     "reference_slope", "params"))
  expect_identical(class(flr), c("slope_sample_size_floor", "slope_result"))
  # No `design`: that is the whole point of the class (CONTRACT.md 4.3).
  expect_null(flr$design)
})

test_that("a floor row binds together with the other two entry points", {
  p <- ref_params()
  rows <- rbind(
    as.data.frame(slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)),
    as.data.frame(slope_power(p, c(0, 1, 2), n = 450, effectiveness = 0.33)),
    as.data.frame(slope_sample_size_floor(p, effectiveness = 0.33))
  )
  expect_identical(nrow(rows), 3L)
  expect_identical(rows$solve_for, c("n", "power", "n_floor"))
  expect_identical(rows$n_follow_up, c(2L, 2L, NA_integer_))
  expect_false(anyNA(rows$var_tte))
})

test_that("printing a floor names the schedule it does not depend on", {
  out <- capture.output(print(slope_sample_size_floor(ref_params(),
                                                     effectiveness = 0.33)))
  expect_true(any(grepl("Lower bound on sample size", out)))
  expect_true(any(grepl("visit schedule.*any", out)))
  expect_true(any(grepl("limiting s\\*\\^2", out)))
})

test_that("the floor's printed N follows per_arm, per arm by default", {
  flr <- slope_sample_size_floor(ref_params(), effectiveness = 0.33)
  expect_identical(attr(flr, "per_arm"), TRUE)

  out_arm   <- capture.output(print(flr))
  out_total <- capture.output(print(flr, per_arm = FALSE))
  expect_true(any(grepl(sprintf("N per arm = %d$", flr$n_per_arm), out_arm)))
  expect_false(any(grepl("^ *N = ", out_arm)))
  expect_true(any(grepl(sprintf("^ *N = %d$", flr$n), out_total)))
  expect_false(any(grepl("N per arm", out_total)))

  expect_error(print(flr, per_arm = "no"), "per_arm")
})

# --- against the paper ------------------------------------------------------

test_that("the floor reproduces the vignette's slpower1 figure", {
  skip_if_not_installed("nlme")
  p1 <- slope_params(sdmt ~ visit | id, data = load_paper_data("slpower1"))
  flr <- slope_sample_size_floor(p1, effectiveness = 0.33)

  # Derived by hand in the "What s* is" vignette, section 6.
  z <- qnorm(0.975) + qnorm(0.8)
  by_hand <- 2 * ceiling(z^2 * 2 *
    (p1$sigma2_slope - p1$sigma_cov^2 / p1$sigma2_intercept) /
    (0.33 * abs(p1$slope))^2)
  expect_identical(flr$n, by_hand)

  # A fifty-year trial with visits every five weeks reaches it exactly, and
  # the paper's own p.588 design (N = 712) needs three times as many.
  expect_identical(slope_sample_size(p1, seq(0, 50, length.out = 501),
                                     effectiveness = 0.33)$n, flr$n)
  expect_gt(slope_sample_size(p1, c(0, 1, 2), effectiveness = 0.33)$n, 2 * flr$n)
})
