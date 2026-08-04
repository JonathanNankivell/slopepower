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

  # A design claiming its stored dropout is cumulative is rejected: the stored
  # vector is always incremental, so honouring the claim would double-convert.
  cumul <- structure(list(visits = c(0, 1, 2), dropout = c(0.05, 0.05),
                          has_dropout = TRUE, dropout_type = "cumulative"),
                     class = "trial_design")
  expect_error(slope_power(p, cumul, effectiveness = 0.33), "already be incremental")
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
