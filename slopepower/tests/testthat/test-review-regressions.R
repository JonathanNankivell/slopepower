# Regressions for defects found in code review of commit 82a4a56.
# Each of these silently returned a wrong number rather than erroring.

test_that("a labelled factor or character group column is rejected, not guessed", {
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
  p <- paper_fit("slpower1")

  # Omitting the field used to fail with "argument is of length zero".
  bare <- structure(list(visits = c(0, 1, 2), dropout = c(0, 0),
                         dropout_type = "incremental"), class = "trial_design")
  expect_silent(r <- slope_sample_size(p, bare, effectiveness = 0.33))
  expect_equal(r$n, slope_sample_size(p, trial_design(c(0, 1, 2)), effectiveness = 0.33)$n)

  # has_dropout = FALSE alongside a non-zero dropout vector used to report the
  # unweighted s*^2 (5.96990) in place of the dropout-weighted value (6.36497),
  # while still returning the correct N -- so the error was invisible in N alone.
  lying <- structure(list(visits = c(0, 1, 2, 5), dropout = c(0, 0, 0.1),
                          has_dropout = FALSE, dropout_type = "incremental"),
                     class = "trial_design")
  correct <- suppressWarnings(
    slope_sample_size(p, trial_design(c(0, 1, 2, 5), c(0, 0, 0.1)), effectiveness = 0.33))
  fixed <- suppressWarnings(slope_sample_size(p, lying, effectiveness = 0.33))
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
  from_cumul <- suppressWarnings(slope_sample_size(p, cumul, effectiveness = 0.33))
  from_incr <- suppressWarnings(slope_sample_size(p, incr, effectiveness = 0.33))
  expect_equal(from_cumul$n, from_incr$n)
  expect_equal(from_cumul$var_tte, from_incr$var_tte, tolerance = 1e-10)
})

test_that("the covariate guard catches term removal and offsets, not just addition", {
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

test_that("slope_bootstrap() cannot be asked to bootstrap its own input", {
  p <- paper_fit("slpower1")
  ss <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  pw <- slope_power(p, c(0, 1, 2), n = 712, effectiveness = 0.33)

  # slope_bootstrap() used to re-specify the calculation through `...`, so a
  # caller could name a statistic that was also one of the inputs: every
  # replicate returned the number they had passed in, giving a zero-width
  # interval that looked like a successful bootstrap. Four hand-written guards
  # existed only to catch those combinations. Dispatching on the result removes
  # the possibility rather than policing it -- a result object already knows
  # which of n and power it solved for -- so the regression is now pinned by
  # what each method will accept at all.
  # The refusal must still diagnose, not just reject: bare match.arg() reports
  # `'arg' should be one of ...`, naming neither the argument nor the object, and
  # the fix is upstream in the call that built the object rather than here.
  err <- tryCatch(slope_bootstrap(ss, R = 5, statistic = "power"),
                  error = conditionMessage)
  expect_match(err, "`statistic`")
  expect_match(err, "bootstrap a slope_power\\(\\) result")
  expect_false(grepl("'arg'", err, fixed = TRUE))

  expect_match(tryCatch(slope_bootstrap(pw, R = 5, statistic = "n"),
                        error = conditionMessage),
               "bootstrap a slope_sample_size\\(\\)\\s+result")

  # And the inputs cannot be smuggled in alongside: they are not silently
  # ignored, which would be the worst outcome of all.
  expect_error(slope_bootstrap(ss, R = 5, n = 400), "unused argument")
  expect_error(slope_bootstrap(pw, R = 5, power = 0.8), "unused argument")

  # A call written against the old interface fails loudly rather than quietly
  # bootstrapping the slope and discarding the design.
  err <- tryCatch(
    slope_bootstrap(p, R = 5, statistic = "n", design = c(0, 1, 2),
                    effectiveness = 0.33),
    error = conditionMessage)
  expect_match(err, "unused argument")
  expect_match(err, "slope_bootstrap\\(slope_sample_size\\(")

  # "tte" is reachable from either result: it depends on neither n nor power.
  b <- suppressWarnings(slope_bootstrap(ss, R = 3, statistic = "tte"))
  expect_equal(b$observed, ss$tte)
  expect_equal(suppressWarnings(slope_bootstrap(pw, R = 3, statistic = "tte"))$observed,
               ss$tte)
})

test_that("slope_bootstrap() re-solves the calculation the result was built with", {
  # The bootstrap must hold every input fixed and vary only the parameters. The
  # observed value therefore has to reproduce the result it was handed, for each
  # of the settings carried on the object -- design, effectiveness, alpha and
  # the target power or sample size.
  p <- paper_fit("slpower1")
  ss <- slope_sample_size(p, trial_design(c(0, 1, 2, 3), dropout = c(0, 0.1, 0.1)),
                          effectiveness = 0.4, power = 0.9, alpha = 0.01)
  expect_equal(suppressWarnings(slope_bootstrap(ss, R = 3))$observed, ss$n)

  pw <- slope_power(p, c(0, 1, 2), n = 401, effectiveness = 0.33, alpha = 0.1)
  # n = 401 is rounded down to 400 for equal arms; the replicates must answer
  # for the number actually used, not the one requested.
  expect_equal(suppressWarnings(slope_bootstrap(pw, R = 3))$observed, pw$power)

  # target = "observed" stores effectiveness as NA and rejects it as an input,
  # so it must be omitted when the call is rebuilt rather than passed along.
  p3 <- paper_fit("slpower3")
  ss3 <- slope_sample_size(p3, c(0, 0.5, 2), target = "observed")
  expect_equal(suppressWarnings(slope_bootstrap(ss3, R = 3))$observed, ss3$n)
})

test_that("slope_bootstrap() does not repeat the stage-two call's own warning", {
  # The caller has already run the stage-two call themselves to produce the
  # object being bootstrapped, and heard anything it had to say. Recomputing the
  # observed value here said it a second time, attributed to slope_bootstrap()
  # -- a function that had not made the choice being warned about.
  set.seed(2)
  s <- data.frame(id = 1:60, case = rep(c(1, 0), each = 30))
  s$intercept <- rnorm(60, 50, 10)
  # controls declining faster than cases, so the comparator slope is further
  # from zero and the target effect makes the treated slope more extreme
  s$slope <- rnorm(60, ifelse(s$case == 1, -0.3, -1.7), 1.2)
  d <- merge(s, data.frame(visit = 0:3))
  d$sdmt <- d$intercept + d$slope * d$visit + rnorm(nrow(d), 0, 3)
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, healthy = case))

  ss <- suppressWarnings(slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33))
  msgs <- character()
  withCallingHandlers(
    slope_bootstrap(ss, R = 2, seed = 1),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("makes the slope more extreme", msgs)))
  # the observed value is still exactly the one the result reported
  expect_equal(suppressWarnings(slope_bootstrap(ss, R = 2, seed = 1))$observed, ss$n)
})

test_that("an explicit n = NULL is caught by the guard, not by the solver", {
  # solve_slope() picks its branch on is.null(n) while the guards were checking
  # missing(n), so `n = NULL` -- what do.call() produces from an argument list
  # with an unset element -- slipped through into the solve-for-n branch and
  # failed with "`power` must be a single finite number", naming an argument
  # slope_power() does not have.
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                           sigma2_slope = 2, sigma_cov = 5, sigma2_residual = 10)

  for (call_it in list(
    function() slope_power(p, c(0, 1, 2), n = NULL),
    function() do.call(slope_power, list(params = p, design = c(0, 1, 2), n = NULL))
  )) {
    err <- tryCatch(call_it(), error = conditionMessage)
    expect_match(err, "`n` is required")
    expect_false(grepl("`power` must be", err))
  }

  err_grid <- tryCatch(slope_power_grid(p, visits = list(a = c(0, 1, 2)), n = NULL),
                       error = conditionMessage)
  expect_match(err_grid, "`n` is required")
  # the old failure surfaced from inside the loop, wrapped per cell
  expect_false(grepl("failed", err_grid))
})

test_that("slope_effect_size() takes no effectiveness argument", {
  p <- paper_fit("slpower1")

  # It used to accept one and ignore it, returning the same number for every
  # value. The result is on the slope-difference scale, so a caller
  # reconstructing N from it must apply the effectiveness factor themselves.
  expect_error(slope_effect_size(p, c(0, 1, 2), effectiveness = 0.33), "unused argument")

  es <- slope_effect_size(p, c(0, 1, 2))
  for (e in c(0.1, 0.33, 0.9)) {
    n_manual <- 2 * ceiling((stats::qnorm(0.975) + stats::qnorm(0.8))^2 / (es * e)^2)
    expect_equal(slope_sample_size(p, c(0, 1, 2), effectiveness = e)$n, n_manual)
  }
})

test_that("both grids enforce the effectiveness/observed guard", {
  p3 <- paper_fit("slpower3")

  # The stage-two functions error on this combination; the grid omitted the
  # argument instead of checking it, and returned a full table computed at
  # effectiveness = 1 with no warning. The omission is still there -- it has to
  # be -- so the check must be raised in each grid wrapper, both of them.
  expect_error(
    slope_sample_size(p3, trial_design(c(0, 2, 3)), target = "observed",
                      effectiveness = 0.9),
    "only one of"
  )
  expect_error(
    slope_power(p3, trial_design(c(0, 2, 3)), n = 400, target = "observed",
                effectiveness = 0.9),
    "only one of"
  )
  expect_error(
    slope_sample_size_grid(p3, visits = list(a = c(0, 2, 3)),
                           target = "observed", effectiveness = 0.9),
    "only one of"
  )
  expect_error(
    slope_power_grid(p3, visits = list(a = c(0, 2, 3)), n = 400,
                     target = "observed", effectiveness = 0.9),
    "only one of"
  )

  # Not supplying it remains valid and reproduces the published N = 318.
  g <- suppressWarnings(
    slope_sample_size_grid(p3, visits = list(a = c(0, 2, 3)),
                           dropout = list(d = c(0.2, 0.1)), target = "observed"))
  expect_equal(g$n[1], 318)
})
