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
  r <- slope_sample_size(p, visits, effectiveness = e, power = 0.8, alpha = 0.05)

  s_star <- sqrt(slope_var(p, visits))
  beta2  <- -e * p$slope
  z <- stats::qnorm(0.975) + stats::qnorm(0.8)
  expect_equal(r$n_per_arm, ceiling((z * s_star / beta2)^2))
  expect_equal(r$tte, beta2)
})

test_that("halving effectiveness quadruples the sample size", {
  p <- ref_params()
  n_full <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.4)$n_per_arm
  n_half <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.2)$n_per_arm
  # exact up to the two ceiling() roundings
  expect_lt(abs(n_half - 4 * n_full), 5)
})

# --- reference slope dispatch ----------------------------------------------

test_that("the reference slope follows the comparator and target", {
  none    <- ref_params("none")
  healthy <- ref_params("healthy", slope = -1.715, slope_comparator = 0.975)
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)

  expect_equal(slope_sample_size(none, c(0, 1, 2))$reference_slope, 0)
  expect_equal(slope_sample_size(healthy, c(0, 1, 2))$reference_slope, 0.975)
  # target = "effectiveness" on trial data ignores the treated arm (Stata's
  # model-3 default) -- the reference is zero, not the treated slope.
  expect_equal(slope_sample_size(treated, c(0, 1, 2))$reference_slope, 0)
  expect_equal(slope_sample_size(treated, c(0, 1, 2), target = "observed")$reference_slope,
               -1.104)
})

test_that("the reference slope is resolved identically by both entry points", {
  # The dispatch lives in shared machinery, so asking for a sample size and
  # asking for power must never disagree about what is being targeted.
  healthy <- ref_params("healthy", slope = -1.715, slope_comparator = 0.975)
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)
  for (p in list(ref_params("none"), healthy, treated)) {
    a <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
    b <- slope_power(p, c(0, 1, 2), n = 400, effectiveness = 0.33)
    expect_equal(a$reference_slope, b$reference_slope)
    expect_equal(a$slope_difference, b$slope_difference)
    expect_equal(a$tte, b$tte)
    expect_equal(a$effect_size, b$effect_size)
  }
  a <- slope_sample_size(treated, c(0, 2, 3), target = "observed")
  b <- slope_power(treated, c(0, 2, 3), n = 400, target = "observed")
  expect_equal(a$reference_slope, b$reference_slope)
  expect_equal(a$tte, b$tte)
})

test_that("slope_difference and tte follow the contract formulas", {
  healthy <- ref_params("healthy", slope = -1.715, slope_comparator = 0.975)
  r <- slope_sample_size(healthy, c(0, 1, 2), effectiveness = 0.33)
  expect_equal(r$slope_difference, -1.715 - 0.975)
  expect_equal(r$tte, -0.33 * (-1.715 - 0.975))
})

test_that('target = "observed" forces effectiveness to 1 and reports it as NA', {
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)
  r <- slope_sample_size(treated, c(0, 2, 3), target = "observed")
  expect_true(is.na(r$effectiveness))
  expect_identical(r$target, "observed")
  expect_equal(r$tte, -(-1.852 - -1.104))
  expect_equal(r$slope_difference, -1.852 - -1.104)
})

test_that('target = "observed" requires trial data and rejects effectiveness', {
  none    <- ref_params("none")
  treated <- ref_params("treated", slope = -1.852, slope_comparator = -1.104)
  expect_error(slope_sample_size(none, c(0, 1, 2), target = "observed"),
               'comparator = "treated"')
  expect_error(
    slope_sample_size(treated, c(0, 2, 3), target = "observed", effectiveness = 0.5),
    "only one of"
  )
  # and the same on the power side
  expect_error(slope_power(none, c(0, 1, 2), n = 400, target = "observed"),
               'comparator = "treated"')
  expect_error(
    slope_power(treated, c(0, 2, 3), n = 400, target = "observed", effectiveness = 0.5),
    "only one of"
  )
})

# --- the two entry points are inverses --------------------------------------

test_that("slope_sample_size() and slope_power() invert each other", {
  p <- ref_params()
  for (target_power in c(0.7, 0.8, 0.9, 0.95)) {
    rn <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33, power = target_power)
    rp <- slope_power(p, c(0, 1, 2), n = rn$n, effectiveness = 0.33)
    expect_gte(rp$power, target_power)
    expect_lt(rp$power - target_power, 0.005)   # only the ceiling() gap
  }
})

test_that("the round trip holds under dropout and a non-default alpha too", {
  p <- ref_params()
  d <- trial_design(c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
  for (a in c(0.01, 0.05, 0.10)) {
    rn <- slope_sample_size(p, d, effectiveness = 0.33, power = 0.9, alpha = a)
    rp <- slope_power(p, d, n = rn$n, effectiveness = 0.33, alpha = a)
    expect_gte(rp$power, 0.9)
    expect_lt(rp$power - 0.9, 0.005)

    # var_tte is back-solved separately in the two branches and they do NOT
    # agree exactly. In the power branch the n cancels and the result is the
    # true dropout-weighted s*^2; in the sample-size branch n has been rounded
    # up, so the reported value is inflated by exactly that rounding. The gap is
    # therefore bounded by one participant per arm, and never negative.
    expect_gte(rn$var_tte, rp$var_tte)
    expect_lt(rn$var_tte / rp$var_tte, rn$n_per_arm / (rn$n_per_arm - 1))
  }
})

test_that("var_tte from slope_power() is exact, and free of n", {
  # pnorm/qnorm invert, so (z_a + qnorm(power))^2 = scaled_effect^2 * n_per_arm
  # and the n_per_arm in the numerator cancels: the same design must report the
  # same effective variance at every sample size.
  p <- ref_params()
  d <- trial_design(c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
  vals <- vapply(c(100, 328, 1000, 5000),
                 function(nn) slope_power(p, d, n = nn, effectiveness = 0.33)$var_tte,
                 numeric(1L))
  expect_equal(vals, rep(vals[1L], length(vals)), tolerance = 1e-10)
})

test_that("the class records which question was asked", {
  p <- ref_params()
  rn <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  expect_s3_class(rn, "slope_sample_size")
  expect_s3_class(rn, "slope_result")
  expect_equal(rn$power, 0.8)

  rp <- slope_power(p, c(0, 1, 2), n = 450, effectiveness = 0.33)
  expect_s3_class(rp, "slope_power")
  expect_s3_class(rp, "slope_result")
  expect_equal(rp$n_requested, 450)

  # as.data.frame() still records it, so rows from both bind meaningfully
  expect_identical(as.data.frame(rn)$solve_for, "n")
  expect_identical(as.data.frame(rp)$solve_for, "power")
})

test_that("an odd n is reduced by one to keep the arms equal", {
  p <- ref_params()
  r <- slope_power(p, c(0, 1, 2), n = 451, effectiveness = 0.33)
  expect_equal(r$n, 450)
  expect_equal(r$n_per_arm, 225)
  expect_equal(r$n_requested, 451)
})

test_that("n is always twice n_per_arm", {
  p <- ref_params()
  for (pw in c(0.8, 0.9)) {
    r <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33, power = pw)
    expect_equal(r$n, 2 * r$n_per_arm)
  }
  for (nn in c(450, 451)) {
    r <- slope_power(p, c(0, 1, 2), n = nn, effectiveness = 0.33)
    expect_equal(r$n, 2 * r$n_per_arm)
  }
})

# --- dropout ----------------------------------------------------------------

test_that("dropout increases the required sample size", {
  p <- ref_params()
  no_drop <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)$n
  drop    <- slope_sample_size(p, trial_design(c(0, 1, 2), c(0, 0.2)),
                               effectiveness = 0.33)$n
  expect_gt(drop, no_drop)
})

test_that("a zero dropout vector matches no dropout at all", {
  p <- ref_params()
  a <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  b <- slope_sample_size(p, trial_design(c(0, 1, 2), c(0, 0)), effectiveness = 0.33)
  expect_equal(a$n, b$n)
  expect_equal(a$effect_size, b$effect_size)
})

test_that("baseline-only dropouts contribute nothing to the effect size", {
  # Stratum j = 1 is skipped: a single measurement carries no slope information.
  p <- ref_params()
  d <- suppressWarnings(trial_design(c(0, 1, 2), c(0.2, 0)))
  r <- suppressWarnings(slope_sample_size(p, d, effectiveness = 0.33))
  full <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  # 20% contribute zero, so the squared effect size is 80% of the complete one
  expect_equal(r$effect_size^2, 0.8 * full$effect_size^2)
})

test_that("var_tte equals slope_var when there is no dropout", {
  p <- ref_params()
  expect_equal(slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)$var_tte,
               slope_var(p, c(0, 1, 2)))
  expect_equal(slope_power(p, c(0, 1, 2), n = 450, effectiveness = 0.33)$var_tte,
               slope_var(p, c(0, 1, 2)))
})

test_that("var_tte under dropout back-solves the effective variance", {
  p <- ref_params()
  d <- trial_design(c(0, 1, 2), c(0, 0.1))
  r <- slope_sample_size(p, d, effectiveness = 0.33)
  z <- stats::qnorm(1 - r$alpha / 2) + stats::qnorm(r$power)
  expect_equal(r$var_tte, r$n_per_arm * r$tte^2 / z^2)
  # and it sits between the complete-case and worst-stratum variances
  expect_gt(r$var_tte, slope_var(p, c(0, 1, 2)))
})

test_that("total dropout errors rather than returning a missing sample size", {
  p <- ref_params()
  d <- suppressWarnings(trial_design(c(0, 1, 2), c(1, 0)))
  expect_error(suppressWarnings(slope_sample_size(p, d, effectiveness = 0.33)),
               "effect size is zero")
  expect_error(suppressWarnings(slope_power(p, d, n = 400, effectiveness = 0.33)),
               "effect size is zero")
})

# --- guards -----------------------------------------------------------------

test_that("slope_sample_size() enforces the contract guards", {
  p <- ref_params()
  expect_error(slope_sample_size(p, c(0, 1, 2), alpha = 0), "alpha")
  expect_error(slope_sample_size(p, c(0, 1, 2), alpha = 1), "alpha")
  expect_error(slope_sample_size(p, c(0, 1, 2), power = 0), "power")
  expect_error(slope_sample_size(p, c(0, 1, 2), power = 1), "power")
  expect_error(slope_sample_size(p, c(0, 1, 2), effectiveness = 0), "effectiveness")
  expect_error(slope_sample_size(p, c(0, 1, 2), effectiveness = 1.5), "effectiveness")
})

test_that("slope_power() enforces the contract guards", {
  p <- ref_params()
  expect_error(slope_power(p, c(0, 1, 2), n = 1), "`n`")
  expect_error(slope_power(p, c(0, 1, 2), n = 100.5), "whole number")
  expect_error(slope_power(p, c(0, 1, 2), n = 100, alpha = 0), "alpha")
  expect_error(slope_power(p, c(0, 1, 2), n = 100, alpha = 1), "alpha")
  expect_error(slope_power(p, c(0, 1, 2), n = 100, effectiveness = 0), "effectiveness")
  expect_error(slope_power(p, c(0, 1, 2), n = 100, effectiveness = 1.5), "effectiveness")
})

test_that("neither entry point accepts the other's question as an argument", {
  # This is the point of the split: the ambiguity that used to be resolved by
  # which argument was left NULL is now impossible to express. Supplying `n` to
  # the sample-size calculation, or `power` to the power calculation, is an
  # error from R's own argument matching rather than a silent mode switch.
  p <- ref_params()
  expect_error(slope_sample_size(p, c(0, 1, 2), n = 450), "unused argument")
  expect_error(slope_power(p, c(0, 1, 2), n = 450, power = 0.9), "unused argument")
  # and `n` is not optional: there is no default sample size to fall back on
  expect_error(slope_power(p, c(0, 1, 2)), "`n` is required")
  expect_match(tryCatch(slope_power(p, c(0, 1, 2)), error = conditionMessage),
               "slope_sample_size")
})

test_that("effectiveness of exactly 1 is allowed", {
  p <- ref_params()
  expect_no_error(slope_sample_size(p, c(0, 1, 2), effectiveness = 1))
  expect_no_error(slope_power(p, c(0, 1, 2), n = 400, effectiveness = 1))
})

test_that("a zero slope difference errors instead of returning NA", {
  # Stata returns N = . silently here (ado:636 divides by zero).
  p <- ref_params("healthy", slope = 1, slope_comparator = 1)
  expect_error(slope_sample_size(p, c(0, 1, 2)), "slope difference is exactly zero")
  expect_error(slope_power(p, c(0, 1, 2), n = 400), "slope difference is exactly zero")
})

test_that("both entry points reject objects that are not slope_params", {
  expect_error(slope_sample_size(list(slope = 1), c(0, 1, 2)), "slope_params")
  expect_error(slope_power(list(slope = 1), c(0, 1, 2), n = 400), "slope_params")
})

test_that("a design must begin at baseline", {
  p <- ref_params()
  expect_error(slope_sample_size(p, c(1, 2)), "baseline")
  expect_error(slope_power(p, c(1, 2), n = 400), "baseline")
})

# --- result objects ---------------------------------------------------------

test_that("slope_sample_size() returns exactly the contract fields", {
  r <- slope_sample_size(ref_params(), c(0, 1, 2), effectiveness = 0.33)
  expect_setequal(names(r), c(
    "n", "n_per_arm", "power", "alpha", "effectiveness", "target", "tte",
    "var_tte", "effect_size", "slope_difference", "reference_slope",
    "params", "design"))
  expect_length(names(r), 13L)
  # n_requested belongs to the power branch: nothing was requested here
  expect_false("n_requested" %in% names(r))
  # solve_for is gone from the objects; the class carries it
  expect_false("solve_for" %in% names(r))
})

test_that("slope_power() returns exactly the contract fields", {
  r <- slope_power(ref_params(), c(0, 1, 2), n = 450, effectiveness = 0.33)
  expect_setequal(names(r), c(
    "n", "n_per_arm", "n_requested", "power", "alpha", "effectiveness",
    "target", "tte", "var_tte", "effect_size", "slope_difference",
    "reference_slope", "params", "design"))
  expect_length(names(r), 14L)
  expect_false("solve_for" %in% names(r))
})

test_that("as.data.frame() column names are stable across every scenario", {
  # Fixes the Stata r(table) defect, where columns silently renamed between
  # models (slope_cases/slope_controls vs slope_untreated/slope_treated). The
  # split adds a second axis to hold stable: sample-size rows and power rows
  # must still bind together.
  scenarios <- list(
    slope_sample_size(ref_params("none"), c(0, 1, 2), effectiveness = 0.33),
    slope_sample_size(ref_params("healthy", -1.715, 0.975), c(0, 1, 2), effectiveness = 0.33),
    slope_sample_size(ref_params("treated", -1.852, -1.104), c(0, 1, 2), effectiveness = 0.33),
    slope_sample_size(ref_params("treated", -1.852, -1.104), c(0, 2, 3), target = "observed"),
    slope_power(ref_params("none"), c(0, 1, 2), n = 450, effectiveness = 0.33),
    slope_power(ref_params("treated", -1.852, -1.104), c(0, 2, 3), n = 450, target = "observed"),
    slope_sample_size(ref_params("none"), trial_design(c(0, 1, 2), c(0, 0.1)),
                      effectiveness = 0.33)
  )
  frames <- lapply(scenarios, as.data.frame)
  reference <- names(frames[[1]])
  for (f in frames) expect_identical(names(f), reference)
  for (f in frames) expect_equal(nrow(f), 1L)
  bound <- do.call(rbind, frames)
  expect_equal(nrow(bound), length(scenarios))
  expect_setequal(unique(bound$solve_for), c("n", "power"))
})

test_that("slope_effect_size() agrees with the value inside both entry points", {
  p <- ref_params()
  d <- trial_design(c(0, 1, 2), c(0, 0.1))
  es <- slope_effect_size(p, d)
  expect_equal(es, slope_sample_size(p, d, effectiveness = 0.33)$effect_size)
  expect_equal(es, slope_power(p, d, n = 400, effectiveness = 0.33)$effect_size)
})

test_that("the effect size takes the sign of the slope difference", {
  falling <- ref_params("none", slope = -1.672)
  rising  <- ref_params("none", slope =  1.672)
  expect_lt(slope_effect_size(falling, c(0, 1, 2)), 0)
  expect_gt(slope_effect_size(rising,  c(0, 1, 2)), 0)
  # only the magnitude drives the sample size
  expect_equal(slope_sample_size(falling, c(0, 1, 2), effectiveness = 0.33)$n,
               slope_sample_size(rising,  c(0, 1, 2), effectiveness = 0.33)$n)
})

# --- printing ---------------------------------------------------------------

test_that("print.slope_sample_size() runs for every scenario, including NA fields", {
  cases <- list(
    slope_sample_size(ref_params("none"), c(0, 1, 2), effectiveness = 0.33),
    slope_sample_size(ref_params("healthy", -1.715, 0.975), c(0, 1, 2), effectiveness = 0.33),
    slope_sample_size(ref_params("treated", -1.852, -1.104), c(0, 2, 3), target = "observed")
  )
  for (r in cases) expect_output(print(r), "Parameters for planned study")
  for (r in cases) expect_output(print(r), "Estimated sample size")
  # manual params have no n_obs; Stata prints "." for missing values
  out <- capture.output(print(cases[[1]]))
  expect_true(any(grepl("number of observations in model = \\.", out)))
  # the target power is an input here, so it is shown and no N is specified
  expect_true(any(grepl("^ +power = 0\\.800", out)))
  expect_false(any(grepl("specified N", out)))
  # effectiveness is missing under target = "observed"
  expect_true(any(grepl("effectiveness = \\.", capture.output(print(cases[[3]])))))
})

test_that("print.slope_power() shows the requested N and the achieved power", {
  r <- slope_power(ref_params("none"), c(0, 1, 2), n = 451, effectiveness = 0.33)
  out <- capture.output(print(r))
  expect_true(any(grepl("Parameters for planned study", out)))
  expect_true(any(grepl("specified N = 451", out)))
  expect_true(any(grepl("actual N = 450", out)))
  expect_true(any(grepl("Estimated power", out)))
  expect_false(any(grepl("Estimated sample size", out)))
})

test_that("the printed schedule renders the dropout at each follow-up visit", {
  # The "schedule (and dropouts)" line is built by slopepower.ado:285-292 and
  # printed at :744; reproducing it is the whole job of schedule_string(),
  # because that line is how a reader checks a result against the worked
  # examples in Nash et al. (2021). Two renderings, chosen by has_dropout, and
  # both are pinned here: the file's other print tests exercise the no-dropout
  # branch without asserting what it produces, and nothing reached the dropout
  # branch at all.
  schedule_of <- function(x) {
    line <- grep("schedule (and dropouts)", capture.output(print(x)),
                 fixed = TRUE, value = TRUE)
    sub("^[^=]*= ", "", line)
  }
  p <- ref_params("healthy", -1.715, 0.975)

  # The design of the paper's p.593 power example: 5% of the randomised cohort
  # withdraws after each of two annual visits. trial_design() warns because
  # dropout[1] covers baseline-only attenders (CONTRACT.md 5.4) -- expected for
  # a flat rate, and not what is under test here.
  d593 <- suppressWarnings(trial_design(c(0, 1, 2), c(0.05, 0.05)))
  expect_identical(schedule_of(slope_power(p, d593, n = 200, effectiveness = 0.33)),
                   "1 (0.05), 2 (0.05)")

  # A zero keeps its "(0)" rather than the visit being dropped from the list --
  # ado:288 is an explicit else branch for exactly that. This is the p.588
  # extended design, where only the final visit loses anyone.
  expect_identical(
    schedule_of(slope_sample_size(p, trial_design(c(0, 1, 2, 5), c(0, 0, 0.1)),
                                  effectiveness = 0.33)),
    "1 (0), 2 (0), 5 (0.1)")

  # Stata's string is built from the integer schedule positions because Sigma is
  # assembled on a unit grid; this port carries real visit times (CONTRACT.md
  # 5.1), so the six-monthly row of Table 1, p.595 prints half-year times rather
  # than 1..6.
  d_half <- suppressWarnings(trial_design(seq(0, 3, 0.5), rep(0.025, 6)))
  expect_identical(
    schedule_of(suppressWarnings(slope_power(p, d_half, n = 450, effectiveness = 0.33))),
    "0.5 (0.025), 1 (0.025), 1.5 (0.025), 2 (0.025), 2.5 (0.025), 3 (0.025)")

  # Without dropout the parenthesised rates vanish entirely -- ado:248-252 uses
  # a separate loop for that case -- on both print methods.
  expect_identical(schedule_of(slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)),
                   "1, 2")
  expect_identical(schedule_of(slope_power(p, seq(0, 3, 0.5), n = 450, effectiveness = 0.33)),
                   "0.5, 1, 1.5, 2, 2.5, 3")
})

test_that("the two print methods share the data-characteristics block verbatim", {
  # Split methods must not drift: the block above "Parameters for planned study"
  # is produced by shared code and depends only on params, not on the question.
  p <- ref_params("healthy", -1.715, 0.975)
  head_of <- function(x) {
    out <- capture.output(print(x))
    out[seq_len(which(out == "Parameters for planned study:") - 1L)]
  }
  a <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  b <- slope_power(p, c(0, 1, 2), n = 450, effectiveness = 0.33)
  expect_identical(head_of(a), head_of(b))
})
