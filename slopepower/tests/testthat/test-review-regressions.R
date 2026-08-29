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
  p <- ref_params()

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
                           dropout = list(d = c(0.2, 0.1)), target = "observed",
                           per_arm = FALSE))
  expect_equal(g$n[1], 318)
})

# --- regressions from the full-package review of 9540a0f --------------------
# Each of these silently returned a wrong number, crashed, or stayed silent
# where CONTRACT.md requires an error or a warning.

test_that("check_params() rejects a non-PD random-effects matrix even when the marginal covariance stays PD", {
  # new_slope_params() checks this at construction for both public routes into
  # the class, so it is unreachable through slope_params() or
  # slope_params_manual(); a hand-built object bypasses it entirely. A large
  # sigma2_residual keeps the *marginal* Sigma built by sigma_at() positive
  # definite even though the random-effects matrix G = [[1, 100], [100, 1]]
  # is not (eigenvalues 101, -99), so that check alone used to let this
  # through, and slope_sample_size() returned N = 57,694,340 instead of
  # erroring.
  bad <- structure(
    list(slope = -1, slope_comparator = NA_real_, comparator = "none",
        sigma2_intercept = 1, sigma2_slope = 1, sigma_cov = 100,
        sigma2_residual = 1e6, n_obs = NA_integer_, n_subjects = NA_integer_,
        common_variance = FALSE, time_shifted = FALSE, fit = NULL, call = NULL),
    class = "slope_params")
  expect_true(is_positive_definite(matrix(c(bad$sigma2_intercept, bad$sigma_cov,
                                            bad$sigma_cov, bad$sigma2_slope), 2L)) == FALSE)
  expect_error(slope_sample_size(bad, c(0, 1, 2), effectiveness = 0.33),
              "positive definite")
})

test_that("slope_bootstrap() rejects a non-integer R instead of crashing sprintf() on a failed replicate", {
  p <- paper_fit("slpower1")
  ss <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  # R = 10.7 used to reach sprintf('%d of %d replicates failed...', n_failed, R)
  # as soon as one replicate failed to converge, and fail there instead with
  # "invalid format '%d'; use format %f, %e, %g or %a for numeric objects" --
  # an unrelated crash rather than a clean validation error.
  expect_error(slope_bootstrap(ss, R = 10.7, seed = 1), "whole number")
})

test_that("as_trial_design() warns on dropout[1] > 0 for a design that was never validated by trial_design()", {
  p <- paper_fit("slpower1")

  # A `trial_design` object built by the constructor, then edited directly:
  # the object is genuine, but this exact dropout vector never went through
  # trial_design()'s own warning check.
  mutated <- trial_design(c(0, 1, 2))
  mutated$dropout <- c(0.4, 0)
  mutated$has_dropout <- TRUE
  expect_warning(slope_sample_size(p, mutated, effectiveness = 0.33),
                "contribute nothing")

  # A hand-built object never carries the constructor's
  # `slopepower_checked_dropout` attribute at all.
  bare <- structure(list(visits = c(0, 2, 3), dropout = c(0.2, 0.1),
                        has_dropout = TRUE, dropout_type = "incremental"),
                    class = "trial_design")
  expect_warning(slope_sample_size(p, bare, effectiveness = 0.33),
                "contribute nothing")

  # The ordinary path -- build, then use immediately, unmodified -- still
  # warns exactly once (at construction), not twice.
  n_warn <- 0
  withCallingHandlers(
    {
      d <- trial_design(c(0, 2, 3), dropout = c(0.2, 0.1))
      slope_sample_size(p, d, effectiveness = 0.33)
    },
    warning = function(w) { n_warn <<- n_warn + 1; invokeRestart("muffleWarning") }
  )
  expect_equal(n_warn, 1L)
})

test_that("an invalid design in one grid cell is named, like an evaluate() failure is", {
  p <- paper_fit("slpower1")
  # "b" doesn't start at 0; this used to raise trial_design()'s own message
  # verbatim, with no mention of slope_sample_size_grid() or which named
  # design failed -- the wrapping every other grid-cell failure gets.
  err <- expect_error(
    slope_sample_size_grid(p, visits = list(a = c(0, 1, 2), b = c(1, 2, 3)),
                           power = 0.8, effectiveness = 0.33),
    "design \"b\""
  )
  expect_match(conditionMessage(err), "slope_sample_size_grid\\(\\)")
})

test_that("the tte-direction warning is deduplicated once per grid, like the baseline-dropout warning", {
  # comparator = "treated" with a comparator slope more extreme than the
  # treated one: effect_components() warns on every cell under
  # target = "observed". Used to fire once per cell instead of once per grid.
  p3 <- slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                            sigma_cov = 5, sigma2_residual = 10,
                            slope_comparator = -3, comparator = "treated")
  n_warn <- 0
  withCallingHandlers(
    slope_sample_size_grid(p3, visits = list(a = c(0, 1, 2), b = c(0, 1, 3), c = c(0, 2, 3)),
                           target = "observed", power = 0.8),
    warning = function(w) { n_warn <<- n_warn + 1; invokeRestart("muffleWarning") }
  )
  expect_equal(n_warn, 1L)
})

test_that("slopepower() validates dropouts before fitting stage one, not after", {
  d <- load_paper_data("slpower1")
  # 3 dropouts for a 2-visit schedule. If this error came from inside
  # trial_design(), it can only have been raised before the REML fit that
  # slope_params() performs -- there is no other way to observe the ordering
  # from outside, since the failure message is identical either way.
  expect_error(
    slopepower(d, "sdmt", "id", "visit", schedule = c(1, 2),
              dropouts = c(0.1, 0.1, 0.1), obs = TRUE, nocontrols = TRUE,
              effectiveness = 0.33),
    "one element per follow-up visit"
  )
})

# --- regressions from the full-package review of eeb5f4b --------------------

test_that("a multi-level subject expression is refused, not evaluated arithmetically", {
  # The mirror of the covariate guard above, on the other side of the bar. The
  # subject term is *evaluated*, not expanded as a formula, so `site/id` became
  # the quotient and `id + site` the sum, and factor() of either invented
  # participant identifiers that pooled unrelated rows. The fit then succeeded:
  # on slpower1 `sdmt ~ visit | site/id` returned a slope of -0.753 and claimed
  # 534 participants, against the true -1.6725 over 200. Nested and crossed
  # groupings are how nlme and lme4 users spell a second clustering level, and
  # this package has exactly one (see ?slope_params, "further levels of
  # clustering"), so the expression can only ever be a mistake.
  d <- load_paper_data("slpower1")
  d$site <- rep(1:4, length.out = nrow(d))

  for (f in list(sdmt ~ time | site/id, sdmt ~ time | id + site,
                 sdmt ~ time | id:site)) {
    err <- tryCatch(slope_params(f, d), error = conditionMessage)
    expect_match(err, "subject identifier must be a single grouping term")
    expect_match(err, "one level of clustering")
  }

  # A transformation of the identifier is still a single term, exactly as it is
  # on the time side: `factor()` and friends are ordinary calls, not operators
  # that combine model terms.
  ref <- slope_params(sdmt ~ time | id, d)
  expect_equal(slope_params(sdmt ~ time | factor(id), d)$slope, ref$slope,
               tolerance = 1e-12)
  expect_equal(slope_params(sdmt ~ time | factor(id), d)$n_subjects, 200L)
})

test_that("a group indicator that changes within a participant is refused", {
  # `sp_case` is read row by row by every model, so a participant coded 1 at
  # some visits and 0 at others is fitted as a case for part of their follow-up
  # and a control for the rest -- loading on both random-effects blocks under
  # `healthy`, switching arms mid-trial under `treated`. The fit converges and
  # looks entirely ordinary. It also defeats the group-size check, which counts
  # such a participant once in each group.
  d <- load_paper_data("slpower2")
  flip <- d$id <= 5 & d$vdate > as.Date("2010-06-01")
  d$case[flip] <- 1 - d$case[flip]
  err <- tryCatch(suppressMessages(slope_params(sdmt ~ time | id, d, healthy = case)),
                  error = conditionMessage)
  expect_match(err, "`healthy` must be constant within a participant")
  expect_match(err, "5 of 500 participant")

  # And `treated` is held to the same rule, in its own name.
  d3 <- load_paper_data("slpower3")
  d3$treat[d3$id == 1 & d3$visit == 2] <- 1 - d3$treat[d3$id == 1 & d3$visit == 2]
  expect_error(slope_params(sdmt ~ time | id, d3, treated = treat),
               "`treated` must be constant within a participant")

  # The unmodified data are unaffected: the guard costs the paper's own fits
  # nothing.
  expect_equal(paper_fit("slpower3")$n_subjects, 150L)
})

test_that("an alpha too small for 1 - alpha/2 errors instead of returning N = Inf", {
  # z_alpha() keeps Stata's `invnormal(1 - 0.5 * alpha)` spelling, because
  # `qnorm(alpha/2, lower.tail = FALSE)` differs in the last ulp and ceiling()
  # can turn that into an off-by-one in a published N. The cost of parity is
  # that 1 - alpha/2 saturates at 1 below about 2.2e-16, qnorm(1) is Inf and
  # ceiling(Inf) is Inf: slope_sample_size() reported n = Inf, "no sample size
  # is enough", where the true requirement is 7,584 -- an order of magnitude
  # from the 712 the same call gives at alpha = 0.05. Stata reports a missing N
  # here; CONTRACT.md 6 says the port says so instead.
  p <- ref_params()
  for (a in c(1e-16, 1e-17, 1e-300)) {
    err <- tryCatch(slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33, alpha = a),
                    error = conditionMessage)
    expect_match(err, "too small to compute a critical value")
    expect_false(is.infinite(suppressWarnings(as.numeric(err))))
  }
  # The floor shares the arithmetic, so it shares the guard.
  expect_error(slope_sample_size_floor(p, effectiveness = 0.33, alpha = 1e-17),
               "too small to compute a critical value")

  # Everything an actual study would use is untouched, and still finite.
  for (a in c(0.05, 0.01, 1e-8, 1e-14)) {
    expect_true(is.finite(slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33,
                                            alpha = a)$n))
  }
  # And the default is bit-for-bit what it was: parity with the .ado, not a
  # numerically nicer formula, is what z_alpha() is written to preserve.
  expect_identical(slopepower:::z_alpha(0.05, "ctx"), stats::qnorm(1 - 0.05 / 2))
})

test_that("whole numbers print in full rather than in scientific notation", {
  # fmt_line()'s digits = 0 branch exists to print a whole number as a whole
  # number, and bare format() does not: it switches to scientific notation
  # whenever that is the shorter string. So an N of exactly 100000 printed as
  # "1e+05" in a block that is otherwise a transcription of Stata's %5.0f
  # output -- and so did the observation count of a 100000-row model.
  fmt_line <- slopepower:::fmt_line
  expect_match(fmt_line("N", 1e5, digits = 0L), "N = 100000$")
  expect_match(fmt_line("N", 1e8, digits = 0L), "N = 100000000$")
  expect_match(fmt_line("N", 712, digits = 0L), "N = 712$")

  p <- ref_params()
  out <- capture.output(print(slope_power(p, c(0, 1, 2), n = 100000,
                                          effectiveness = 0.33),
                              per_arm = FALSE))
  expect_true(any(grepl("specified N = 100000", out, fixed = TRUE)))
  expect_true(any(grepl("actual N = 100000", out, fixed = TRUE)))
  expect_false(any(grepl("1e+05", out, fixed = TRUE)))
})
