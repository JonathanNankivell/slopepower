# Layer 3 --- the mathematical core. See CONTRACT.md section 5.

# A reference parameter set: the simulation truth used to generate slpower1
# (paper appendix), with the slope the paper reports for it.
ref_params <- function(comparator = "none", slope = -1.672,
                       slope_comparator = NA_real_) {
  slope_params_manual(
    slope            = slope,
    sigma2_intercept = 100,
    sigma2_slope     = 2,
    sigma_cov        = 5,
    sigma2_residual  = 10,
    slope_comparator = slope_comparator,
    comparator       = comparator
  )
}

# --- slope_sigma() ----------------------------------------------------------

test_that("slope_sigma() reproduces the closed form printed on paper p.579", {
  p <- ref_params()
  t <- c(0, 1.3, 2.7)
  s2a <- p$sigma2_intercept; s2b <- p$sigma2_slope
  sab <- p$sigma_cov;        s2e <- p$sigma2_residual

  # Sigma* as printed in the paper, for baseline plus two follow-ups.
  expected <- matrix(c(
    s2a + s2e,                      s2a + t[2] * sab,                          s2a + t[3] * sab,
    s2a + t[2] * sab,               s2a + 2 * t[2] * sab + t[2]^2 * s2b + s2e, s2a + (t[2] + t[3]) * sab + t[2] * t[3] * s2b,
    s2a + t[3] * sab,               s2a + (t[2] + t[3]) * sab + t[2] * t[3] * s2b, s2a + 2 * t[3] * sab + t[3]^2 * s2b + s2e
  ), nrow = 3, byrow = TRUE)

  expect_equal(unname(slope_sigma(p, t)), expected)
})

test_that("slope_sigma() is symmetric, positive definite and correctly sized", {
  p <- ref_params()
  for (t in list(c(0, 1, 2), c(0, 0.5, 1, 1.5, 2), c(0, 3), seq(0, 3, 0.25))) {
    s <- slope_sigma(p, t)
    expect_equal(dim(s), c(length(t), length(t)))
    expect_equal(unname(s), unname(t(s)))
    expect_true(all(eigen(s, symmetric = TRUE, only.values = TRUE)$values > 0))
  }
})

test_that("the residual variance appears only on the diagonal", {
  p <- ref_params()
  t <- c(0, 1, 2)
  s <- slope_sigma(p, t)
  # off-diagonal [1,2] has no residual term
  expect_equal(unname(s[1, 2]), p$sigma2_intercept + t[2] * p$sigma_cov)
  expect_equal(unname(s[1, 1]), p$sigma2_intercept + p$sigma2_residual)
})

# --- slope_var() ------------------------------------------------------------

test_that("slope_var() is positive and decreases as visits are added", {
  p <- ref_params()
  v2 <- slope_var(p, c(0, 2))
  v5 <- slope_var(p, c(0, 0.5, 1, 1.5, 2))
  expect_gt(v2, 0)
  expect_lt(v5, v2)   # more spread-out information -> smaller variance
  expect_lt(slope_var(p, 0:3), slope_var(p, c(0, 3)))
})

test_that("adding a visit never increases the treatment-effect variance", {
  # A GLS estimator cannot be made worse by observing more, so s*^2 must be
  # monotone non-increasing under insertion of any extra visit. This is the
  # property that makes the paper's Table 1 conclusion hold.
  p <- ref_params()
  bases <- list(c(0, 2), c(0, 3), c(0, 1, 3), c(0, 0.5, 2))
  for (base in bases) {
    for (add in c(0.25, 1.1, 2.5)) {
      if (add %in% base) next
      expect_lte(slope_var(p, sort(c(base, add))), slope_var(p, base) + 1e-12)
    }
  }
})

test_that("an interior visit can be worth nothing for estimating a slope", {
  # A visit that adds no spread about the centroid contributes no information
  # about the slope: here c(0, 2) and c(0, 1, 2) give exactly the same s*^2.
  # (Whether the gain is exactly zero depends on the variance components, so
  # this pins the specific case rather than a general law.)
  p <- ref_params()
  expect_equal(slope_var(p, c(0, 1, 2)), slope_var(p, c(0, 2)))
  # while an off-centre visit does help
  expect_lt(slope_var(p, c(0, 0.5, 2)), slope_var(p, c(0, 2)))
})

test_that("slope_var() falls as follow-up lengthens", {
  p <- ref_params()
  expect_lt(slope_var(p, c(0, 1, 3)), slope_var(p, c(0, 1, 2)))
})

test_that("slope_var() rejects degenerate visit vectors", {
  p <- ref_params()
  expect_error(slope_var(p, 0), "at least two finite times")
  expect_error(slope_var(p, c(0, 0)), "strictly increasing")
  expect_error(slope_var(p, c(2, 1)), "strictly increasing")
})

# --- the effectiveness scale conversion ------------------------------------

test_that("the effectiveness factor converts scales rather than double counting", {
  # CONTRACT.md 5.5: effectiveness appears inside tte and again in the N
  # formula. N should reduce exactly to eq. (6), N = {(z + z) s* / beta2}^2.
  p <- ref_params()
  e <- 0.33
  visits <- c(0, 1, 2)
  r <- slope_power(p, visits, effectiveness = e, power = 0.8, alpha = 0.05)

  s_star <- sqrt(slope_var(p, visits))
  beta2  <- -e * p$slope
  z <- stats::qnorm(0.975) + stats::qnorm(0.8)
  expect_equal(r$n_per_arm, ceiling((z * s_star / beta2)^2))
  expect_equal(r$tte, beta2)
})

test_that("halving effectiveness quadruples the sample size", {
  p <- ref_params()
  n_full <- slope_power(p, c(0, 1, 2), effectiveness = 0.4)$n_per_arm
  n_half <- slope_power(p, c(0, 1, 2), effectiveness = 0.2)$n_per_arm
  # exact up to the two ceiling() roundings
  expect_lt(abs(n_half - 4 * n_full), 5)
})

# --- reference slope dispatch ----------------------------------------------

test_that("the reference slope follows the comparator and target", {
  none    <- ref_params("none")
  healthy <- ref_params("healthy", slope = -1.715, slope_comparator = 0.975)
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)

  expect_equal(slope_power(none, c(0, 1, 2))$reference_slope, 0)
  expect_equal(slope_power(healthy, c(0, 1, 2))$reference_slope, 0.975)
  # target = "effectiveness" on trial data ignores the treated arm (Stata's
  # model-3 default) -- the reference is zero, not the treated slope.
  expect_equal(slope_power(treated, c(0, 1, 2))$reference_slope, 0)
  expect_equal(slope_power(treated, c(0, 1, 2), target = "observed")$reference_slope,
               -1.104)
})

test_that("slope_difference and tte follow the contract formulas", {
  healthy <- ref_params("healthy", slope = -1.715, slope_comparator = 0.975)
  r <- slope_power(healthy, c(0, 1, 2), effectiveness = 0.33)
  expect_equal(r$slope_difference, -1.715 - 0.975)
  expect_equal(r$tte, -0.33 * (-1.715 - 0.975))
})

test_that('target = "observed" forces effectiveness to 1 and reports it as NA', {
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)
  r <- slope_power(treated, c(0, 2, 3), target = "observed")
  expect_true(is.na(r$effectiveness))
  expect_identical(r$target, "observed")
  expect_equal(r$tte, -(-1.852 - -1.104))
  expect_equal(r$slope_difference, -1.852 - -1.104)
})

test_that('target = "observed" requires trial data and rejects effectiveness', {
  none    <- ref_params("none")
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)
  expect_error(slope_power(none, c(0, 1, 2), target = "observed"),
               'comparator = "treated"')
  expect_error(
    slope_power(treated, c(0, 2, 3), target = "observed", effectiveness = 0.5),
    "only one of"
  )
})

# --- solving in both directions --------------------------------------------

test_that("N and power round-trip consistently", {
  p <- ref_params()
  for (target_power in c(0.7, 0.8, 0.9, 0.95)) {
    rn <- slope_power(p, c(0, 1, 2), effectiveness = 0.33, power = target_power)
    rp <- slope_power(p, c(0, 1, 2), effectiveness = 0.33, n = rn$n)
    expect_gte(rp$power, target_power)
    expect_lt(rp$power - target_power, 0.005)   # only the ceiling() gap
  }
})

test_that("solve_for and n_requested record which direction was taken", {
  p <- ref_params()
  rn <- slope_power(p, c(0, 1, 2), effectiveness = 0.33)
  expect_identical(rn$solve_for, "n")
  expect_true(is.na(rn$n_requested))
  expect_equal(rn$power, 0.8)

  rp <- slope_power(p, c(0, 1, 2), effectiveness = 0.33, n = 450)
  expect_identical(rp$solve_for, "power")
  expect_equal(rp$n_requested, 450)
})

test_that("an odd n is reduced by one to keep the arms equal", {
  p <- ref_params()
  r <- slope_power(p, c(0, 1, 2), effectiveness = 0.33, n = 451)
  expect_equal(r$n, 450)
  expect_equal(r$n_per_arm, 225)
  expect_equal(r$n_requested, 451)
})

test_that("n is always twice n_per_arm", {
  p <- ref_params()
  for (pw in c(0.8, 0.9)) {
    r <- slope_power(p, c(0, 1, 2), effectiveness = 0.33, power = pw)
    expect_equal(r$n, 2 * r$n_per_arm)
  }
})

# --- dropout ----------------------------------------------------------------

test_that("dropout increases the required sample size", {
  p <- ref_params()
  no_drop <- slope_power(p, c(0, 1, 2), effectiveness = 0.33)$n
  drop    <- slope_power(p, trial_design(c(0, 1, 2), c(0, 0.2)),
                         effectiveness = 0.33)$n
  expect_gt(drop, no_drop)
})

test_that("a zero dropout vector matches no dropout at all", {
  p <- ref_params()
  a <- slope_power(p, c(0, 1, 2), effectiveness = 0.33)
  b <- slope_power(p, trial_design(c(0, 1, 2), c(0, 0)), effectiveness = 0.33)
  expect_equal(a$n, b$n)
  expect_equal(a$effect_size, b$effect_size)
})

test_that("baseline-only dropouts contribute nothing to the effect size", {
  # Stratum j = 1 is skipped: a single measurement carries no slope information.
  p <- ref_params()
  d <- suppressWarnings(trial_design(c(0, 1, 2), c(0.2, 0)))
  r <- suppressWarnings(slope_power(p, d, effectiveness = 0.33))
  full <- slope_power(p, c(0, 1, 2), effectiveness = 0.33)
  # 20% contribute zero, so the squared effect size is 80% of the complete one
  expect_equal(r$effect_size^2, 0.8 * full$effect_size^2)
})

test_that("var_tte equals slope_var when there is no dropout", {
  p <- ref_params()
  r <- slope_power(p, c(0, 1, 2), effectiveness = 0.33)
  expect_equal(r$var_tte, slope_var(p, c(0, 1, 2)))
})

test_that("var_tte under dropout back-solves the effective variance", {
  p <- ref_params()
  d <- trial_design(c(0, 1, 2), c(0, 0.1))
  r <- slope_power(p, d, effectiveness = 0.33)
  z <- stats::qnorm(1 - r$alpha / 2) + stats::qnorm(r$power)
  expect_equal(r$var_tte, r$n_per_arm * r$tte^2 / z^2)
  # and it sits between the complete-case and worst-stratum variances
  expect_gt(r$var_tte, slope_var(p, c(0, 1, 2)))
})

test_that("total dropout errors rather than returning a missing sample size", {
  p <- ref_params()
  d <- suppressWarnings(trial_design(c(0, 1, 2), c(1, 0)))
  expect_error(suppressWarnings(slope_power(p, d, effectiveness = 0.33)),
               "effect size is zero")
})

# --- guards -----------------------------------------------------------------

test_that("slope_power() enforces the contract guards", {
  p <- ref_params()
  expect_error(slope_power(p, c(0, 1, 2), n = 100, power = 0.8), "only one of")
  expect_error(slope_power(p, c(0, 1, 2), alpha = 0), "alpha")
  expect_error(slope_power(p, c(0, 1, 2), alpha = 1), "alpha")
  expect_error(slope_power(p, c(0, 1, 2), power = 0), "power")
  expect_error(slope_power(p, c(0, 1, 2), power = 1), "power")
  expect_error(slope_power(p, c(0, 1, 2), effectiveness = 0), "effectiveness")
  expect_error(slope_power(p, c(0, 1, 2), effectiveness = 1.5), "effectiveness")
  expect_error(slope_power(p, c(0, 1, 2), n = 1), "n")
  expect_error(slope_power(p, c(0, 1, 2), n = 100.5), "whole number")
})

test_that("effectiveness of exactly 1 is allowed", {
  p <- ref_params()
  expect_no_error(slope_power(p, c(0, 1, 2), effectiveness = 1))
})

test_that("a zero slope difference errors instead of returning NA", {
  # Stata returns N = . silently here (ado:636 divides by zero).
  p <- ref_params("healthy", slope = 1, slope_comparator = 1)
  expect_error(slope_power(p, c(0, 1, 2)), "slope difference is exactly zero")
})

test_that("slope_power() rejects objects that are not slope_params", {
  expect_error(slope_power(list(slope = 1), c(0, 1, 2)), "slope_params")
})

test_that("a design must begin at baseline", {
  p <- ref_params()
  expect_error(slope_power(p, c(1, 2)), "baseline")
})

# --- result object ----------------------------------------------------------

test_that("slope_power() returns exactly the contract fields", {
  p <- ref_params()
  r <- slope_power(p, c(0, 1, 2), effectiveness = 0.33)
  expect_s3_class(r, "slope_power")
  expect_setequal(names(r), c(
    "n", "n_per_arm", "n_requested", "power", "alpha", "effectiveness",
    "target", "tte", "var_tte", "effect_size", "slope_difference",
    "reference_slope", "solve_for", "params", "design"))
  expect_length(names(r), 15L)
})

test_that("as.data.frame() column names are stable across every scenario", {
  # Fixes the Stata r(table) defect, where columns silently renamed between
  # models (slope_cases/slope_controls vs slope_untreated/slope_treated).
  scenarios <- list(
    slope_power(ref_params("none"), c(0, 1, 2), effectiveness = 0.33),
    slope_power(ref_params("healthy", -1.715, 0.975), c(0, 1, 2), effectiveness = 0.33),
    slope_power(ref_params("treated", -1.852, -1.104), c(0, 1, 2), effectiveness = 0.33),
    slope_power(ref_params("treated", -1.852, -1.104), c(0, 2, 3), target = "observed"),
    slope_power(ref_params("none"), c(0, 1, 2), effectiveness = 0.33, n = 450),
    slope_power(ref_params("none"), trial_design(c(0, 1, 2), c(0, 0.1)), effectiveness = 0.33)
  )
  frames <- lapply(scenarios, as.data.frame)
  reference <- names(frames[[1]])
  for (f in frames) expect_identical(names(f), reference)
  for (f in frames) expect_equal(nrow(f), 1L)
  expect_no_error(do.call(rbind, frames))
})

test_that("slope_effect_size() agrees with the value inside slope_power()", {
  p <- ref_params()
  d <- trial_design(c(0, 1, 2), c(0, 0.1))
  expect_equal(slope_effect_size(p, d),
               slope_power(p, d, effectiveness = 0.33)$effect_size)
})

test_that("the effect size takes the sign of the slope difference", {
  falling <- ref_params("none", slope = -1.672)
  rising  <- ref_params("none", slope =  1.672)
  expect_lt(slope_effect_size(falling, c(0, 1, 2)), 0)
  expect_gt(slope_effect_size(rising,  c(0, 1, 2)), 0)
  # only the magnitude drives the sample size
  expect_equal(slope_power(falling, c(0, 1, 2), effectiveness = 0.33)$n,
               slope_power(rising,  c(0, 1, 2), effectiveness = 0.33)$n)
})

# --- printing ---------------------------------------------------------------

test_that("print.slope_power() runs for every scenario, including NA fields", {
  cases <- list(
    slope_power(ref_params("none"), c(0, 1, 2), effectiveness = 0.33),
    slope_power(ref_params("healthy", -1.715, 0.975), c(0, 1, 2), effectiveness = 0.33),
    slope_power(ref_params("treated", -1.852, -1.104), c(0, 2, 3), target = "observed"),
    slope_power(ref_params("none"), c(0, 1, 2), effectiveness = 0.33, n = 450)
  )
  for (r in cases) expect_output(print(r), "Parameters for planned study")
  # manual params have no n_obs; Stata prints "." for missing values
  out <- capture.output(print(cases[[1]]))
  expect_true(any(grepl("number of observations in model = \\.", out)))
  # effectiveness is missing under target = "observed"
  out_obs <- capture.output(print(cases[[3]]))
  expect_true(any(grepl("effectiveness = \\.", out_obs)))
})
