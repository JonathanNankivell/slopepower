# Regressions for defects found in code review of commit 82a4a56.
# Each of these silently returned a wrong number rather than erroring.

test_that("a labelled factor or character group column is rejected, not guessed", {
  skip_without_paper_data()
  d <- load_paper_data("slpower2")
  d$g_num <- d$case
  d$g_chr <- ifelse(d$case == 1, "case", "control")
  d$g_fct <- factor(d$g_chr)

  # Alphabetical level order puts "case" first, so mapping by level order would
  # assign case -> 0 and control -> 1: the exact inverse of the documented
  # coding. That returned the healthy controls' slope labelled as the cases',
  # and with it the controls' residual variance (10.699 rather than 10.354) --
  # the very error the by-name variance extraction exists to prevent.
  expect_error(
    suppressMessages(slope_params(sdmt ~ time | id, d, healthy = g_chr)),
    "cannot be determined"
  )
  expect_error(
    suppressMessages(slope_params(sdmt ~ time | id, d, healthy = g_fct)),
    "cannot be determined"
  )

  # The numeric path is unaffected and remains the reference.
  p <- suppressMessages(slope_params(sdmt ~ time | id, d, healthy = g_num))
  expect_equal(p$slope, -1.7153, tolerance = 1e-3)
  expect_equal(p$sigma2_residual, 10.354, tolerance = 1e-2)

  # A factor coded literally "0"/"1" is unambiguous and must still work.
  d$g_01 <- factor(as.character(d$case), levels = c("0", "1"))
  p01 <- suppressMessages(slope_params(sdmt ~ time | id, d, healthy = g_01))
  expect_equal(p01$slope, p$slope, tolerance = 1e-8)
  expect_equal(p01$sigma2_residual, p$sigma2_residual, tolerance = 1e-8)
})

test_that("hand-built trial_design objects have has_dropout derived, not trusted", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")

  # Omitting the field used to fail with "argument is of length zero".
  bare <- structure(list(visits = c(0, 1, 2), dropout = c(0, 0),
                         dropout_type = "incremental"), class = "trial_design")
  expect_silent(r <- slope_power(p, bare, effectiveness = 0.33))
  expect_equal(r$n, slope_power(p, trial_design(c(0, 1, 2)), effectiveness = 0.33)$n)

  # has_dropout = FALSE alongside a non-zero dropout vector used to report the
  # unweighted s*^2 (5.96990) in place of the dropout-weighted value (6.36497),
  # while still returning the correct N -- so the error was invisible in N alone.
  lying <- structure(list(visits = c(0, 1, 2, 5), dropout = c(0, 0, 0.1),
                          has_dropout = FALSE, dropout_type = "incremental"),
                     class = "trial_design")
  correct <- suppressWarnings(
    slope_power(p, trial_design(c(0, 1, 2, 5), c(0, 0, 0.1)), effectiveness = 0.33))
  fixed <- suppressWarnings(slope_power(p, lying, effectiveness = 0.33))
  expect_equal(fixed$var_tte, correct$var_tte, tolerance = 1e-10)
  expect_equal(fixed$n, correct$n)

  # A design built with dropout_type = "cumulative" must work. `dropout_type`
  # records only how the user supplied the vector -- `dropout` is always stored
  # incrementally -- so acting on the label rejected legitimate constructor
  # output, including trial_design()'s own documented example, with an error
  # telling the user to do exactly what they had done.
  cumul <- suppressWarnings(
    trial_design(c(0, 1, 2, 3), dropout = c(0.05, 0.10, 0.15), dropout_type = "cumulative"))
  incr <- suppressWarnings(
    trial_design(c(0, 1, 2, 3), dropout = c(0.05, 0.05, 0.05)))
  from_cumul <- suppressWarnings(slope_power(p, cumul, effectiveness = 0.33))
  from_incr <- suppressWarnings(slope_power(p, incr, effectiveness = 0.33))
  expect_equal(from_cumul$n, from_incr$n)
  expect_equal(from_cumul$var_tte, from_incr$var_tte, tolerance = 1e-10)
})

test_that("the covariate guard catches term removal and offsets, not just addition", {
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  d$age <- 50

  # terms() applies formula semantics -- `-` removes a term and offset() terms
  # are excluded from term.labels -- so a count of term labels alone let these
  # through to be evaluated arithmetically and fitted as the time variable.
  # On slpower1 that returned slopes of -0.043 and -0.020 for the true -1.6725.
  expect_error(slope_params(sdmt ~ time - age | id, d), "single time term")
  expect_error(slope_params(sdmt ~ time + offset(age) | id, d), "single time term")
  expect_error(slope_params(sdmt ~ offset(age) | id, d), "single time term")
  expect_error(slope_params(sdmt ~ time * age | id, d), "single time term")
  expect_error(slope_params(sdmt ~ time:age | id, d), "single time term")

  # Transformations remain legal: they are ordinary calls, not combining
  # operators.
  ref <- slope_params(sdmt ~ time | id, d)
  d$days <- d$time * 365
  expect_equal(slope_params(sdmt ~ I(days / 365) | id, d)$slope, ref$slope,
               tolerance = 1e-6)
})

test_that("slope_bootstrap() rejects a statistic that would bootstrap its own input", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")

  # slope_power() solves for whichever of n and power is absent, so these
  # combinations returned the user's own input in every replicate -- a
  # zero-width interval that looked like a successful bootstrap.
  expect_error(
    slope_bootstrap(p, R = 5, statistic = "power", design = c(0, 1, 2),
                    effectiveness = 0.33),
    'requires `n`')
  expect_error(
    slope_bootstrap(p, R = 5, statistic = "n", design = c(0, 1, 2), n = 400,
                    effectiveness = 0.33),
    'must not be supplied')
})

test_that("slope_effect_size() takes no effectiveness argument", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")

  # It used to accept one and ignore it, returning the same number for every
  # value. The result is on the slope-difference scale, so a caller
  # reconstructing N from it must apply the effectiveness factor themselves.
  expect_error(slope_effect_size(p, c(0, 1, 2), effectiveness = 0.33), "unused argument")

  es <- slope_effect_size(p, c(0, 1, 2))
  for (e in c(0.1, 0.33, 0.9)) {
    n_manual <- 2 * ceiling((stats::qnorm(0.975) + stats::qnorm(0.8))^2 / (es * e)^2)
    expect_equal(slope_power(p, c(0, 1, 2), effectiveness = e)$n, n_manual)
  }
})

test_that("slope_power_grid() enforces the effectiveness/observed guard", {
  skip_without_paper_data()
  p3 <- paper_fit("slpower3")

  # slope_power() errors on this combination; the grid omitted the argument
  # instead of checking it, and returned a full table computed at
  # effectiveness = 1 with no warning.
  expect_error(
    slope_power(p3, trial_design(c(0, 2, 3)), target = "observed", effectiveness = 0.9),
    "only one of"
  )
  expect_error(
    slope_power_grid(p3, visits = list(a = c(0, 2, 3)),
                     target = "observed", effectiveness = 0.9),
    "only one of"
  )

  # Not supplying it remains valid and reproduces the published N = 318.
  g <- suppressWarnings(
    slope_power_grid(p3, visits = list(a = c(0, 2, 3)),
                     dropout = list(d = c(0.2, 0.1)), target = "observed"))
  expect_equal(g$n[1], 318)
})
