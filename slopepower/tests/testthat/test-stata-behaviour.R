# Behaviours of the Stata original that are now established fact rather than
# inference, pinned as R expectations.
#
# These are separate from test-stata-reference.R on purpose. That file replays
# a table and compares columns; this one states, in R, what we learned the
# original actually does -- so the knowledge survives even if the fixtures are
# ever regenerated, and so a reader can see the claims without parsing a CSV.
#
# Every number quoted here came out of `slopepower.ado` under Stata on
# 2026-08-04. Where the R port deviates, the deviation is asserted too, with
# the reason.

# ---------------------------------------------------------------------------
# The residual variance under residuals(independent, by(case))
# ---------------------------------------------------------------------------
#
# This was the port's most dangerous trap and the last open question about the
# .ado, and they turn out to be the same question seen from two sides.
#
# Stata stores the two residual parameters as `lnsig_e:_cons`, the base log-SD
# of the REFERENCE group, and `r_lns2ose:_cons`, a log-RATIO offset. So
# slopepower.ado:371's exp(b[13]+b[14])^2 is the cases' variance and is
# correct. Stata's own output block prints
#     0: var(e) = 10.6989      (controls, the reference level)
#     1: var(e) = 10.35425     (cases)
#
# nlme has the same reference-level convention with the opposite trap:
# varIdent puts the reference level in sigma(fit), and the reference here is
# `control`, so a naive sigma(fit)^2 returns 10.699 -- the controls' variance,
# a 3.3% error into every model-1 sample size. The port extracts by name.

test_that("model 1 takes the cases' residual variance, not the controls'", {
  skip_without_paper_data()
  p <- paper_fit("slpower2")

  # Stata: 1: var(e) = 10.35425 for the cases.
  expect_equal(p$sigma2_residual, 10.35425, tolerance = 1e-4)

  # And emphatically not the controls' 10.6989. The two differ by 3.3%, which
  # is small enough to look like noise in a slope and large enough to move N.
  expect_gt(abs(p$sigma2_residual - 10.6989), 0.3)

  # The reference level really is `control`, so this is not a lucky ordering:
  # sigma(fit)^2 on the raw nlme object is the number we must NOT be using.
  naive <- stats::sigma(p$fit)^2
  expect_equal(naive, 10.6989, tolerance = 1e-3)
  expect_false(isTRUE(all.equal(naive, p$sigma2_residual, tolerance = 1e-3)))
})

test_that("nocontvar leaves the cases' parameters alone", {
  skip_without_paper_data()
  d <- load_paper_data("slpower2")
  full <- paper_fit("slpower2")
  reduced <- suppressMessages(
    slope_params(sdmt ~ time | id, d, healthy = case, common_variance = TRUE))

  # Stata, scale(365), full vs nocontvar:
  #   var(tcase)  2.259091  vs 2.259092
  #   var(cse)  111.9636    vs 111.9638
  #   1: var(e)  10.354254  vs 10.354257
  # Collapsing the controls' random-effects block is meant to buy convergence
  # when the controls are sparse; it is not meant to move the cases, and in
  # Stata it does not, to six figures.
  expect_equal(reduced$sigma2_slope, full$sigma2_slope, tolerance = 1e-4)
  expect_equal(reduced$sigma2_intercept, full$sigma2_intercept, tolerance = 1e-4)
  expect_equal(reduced$sigma2_residual, full$sigma2_residual, tolerance = 1e-4)
  expect_equal(reduced$slope, full$slope, tolerance = 1e-4)
  expect_true(reduced$common_variance)
  expect_false(full$common_variance)

  # The controls' variance is what nocontvar changes, and it changes a lot:
  # Stata's 0: var(e) goes 10.6989 -> 11.36152. The port discards the controls'
  # components, so the visible consequence is only that the fits differ.
  expect_false(isTRUE(all.equal(stats::logLik(reduced$fit),
                                stats::logLik(full$fit))))
})

# ---------------------------------------------------------------------------
# The dropout-length guard is live in Stata, so this is parity, not a fix
# ---------------------------------------------------------------------------
#
# slopepower.ado:239 and :276 write `sched_length ' with a space inside the
# macro name, which looked like it would make the guard at :279 dead code.
# Stata trims the name (the probe `counter ' + 1 evaluates to 8), both counters
# are correct, and a wrong-length dropouts() list returns _rc = 198. The R
# port's length check therefore restates Stata rather than repairing it.

test_that("a dropout list of the wrong length is refused, as in Stata", {
  # Stata: schedule(1 2 3) dropouts(0.05 0.05) -> 198,
  #        "Dropout list must correspond with visit schedule"
  expect_error(trial_design(c(0, 1, 2, 3), c(0.05, 0.05)), "dropout")
  # Stata: schedule(1 2) dropouts(0.05 0.05 0.05) -> 198
  expect_error(trial_design(c(0, 1, 2), c(0.05, 0.05, 0.05)), "dropout")

  # The well-formed call. It warns, because dropout[1] > 0 means some
  # participants attend baseline only -- Stata skips that stratum silently and
  # the port says so -- but it must not error. Stata gives N = 484 here.
  expect_warning(des <- trial_design(c(0, 1, 2, 3), c(0.05, 0.05, 0.05)),
                 "baseline")
  expect_equal(length(des$dropout), 3)
})

test_that("a dropout list summing to exactly 1 in decimal is accepted", {
  # Stata accepts dropouts(0.3 0.3 0.4) and returns N = 1418: its running
  # subtraction does not go negative. The R port compares with a 1e-8
  # tolerance and reaches the same conclusion by a different route, so this is
  # agreement rather than divergence -- worth pinning because naive
  # accumulation makes 1 - 0.3 - 0.3 - 0.4 equal -5.6e-17.
  des <- suppressWarnings(trial_design(c(0, 1, 2, 3), c(0.3, 0.3, 0.4)))
  expect_equal(sum(des$dropout), 1)
  expect_equal(des$dropout, c(0.3, 0.3, 0.4))
  expect_no_error(suppressWarnings(trial_design(c(0, 1, 2, 3),
                                                c(0.25, 0.25, 0.5))))
  # Over one is refused by both: Stata's "Dropouts cannot exceed 100%", rc 198.
  expect_error(trial_design(c(0, 1, 2, 3), c(0.4, 0.4, 0.4)))
})

test_that("everyone dropping out at the first visit is an error, not a number", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")
  # Stata returns _rc = 0 with N = missing: the effect size is 0 because every
  # stratum but the completers is the baseline-only stratum, which carries no
  # slope information, and the completer stratum has weight 0. A silent missing
  # is a worse answer than an error, so the port errors. CONTRACT.md 6.
  des <- suppressWarnings(trial_design(c(0, 1, 2, 3), c(1, 0, 0)))
  expect_error(slope_sample_size(p, des, effectiveness = 0.33), "drop out")
})

# ---------------------------------------------------------------------------
# var_tte at power saturation
# ---------------------------------------------------------------------------

test_that("var_tte stays finite where Stata's back-solve goes missing", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")
  # The baseline-only warning is asserted where it belongs, above; here it is
  # incidental to the design under test.
  des <- suppressWarnings(trial_design(c(0, 1, 2, 3), c(0.05, 0.05, 0.05)))

  # Solving for power with dropout, Stata reports var_tte by inverting the
  # sample size formula through the power it has just computed. Once power
  # saturates at exactly 1, invnormal(1) is missing and var_tte goes with it.
  # Observed: finite at n = 2000 (9.391673916773652), missing from n = 10000.
  # The port uses the algebraically equivalent closed form, which has no
  # n_per_arm in it and so cannot degenerate.
  small <- slope_power(p, des, n = 2000, effectiveness = 0.33)
  expect_equal(small$var_tte, 9.391673916773652, tolerance = 1e-3)

  for (nn in c(10000, 40000, 200000, 1e6)) {
    r <- slope_power(p, des, n = nn, effectiveness = 0.33)
    expect_equal(r$power, 1)
    expect_true(is.finite(r$var_tte))
    expect_gt(r$var_tte, 0)
    # And it is the same number regardless of n, which is the whole point of
    # the closed form: the effective s*^2 does not depend on the sample size.
    expect_equal(r$var_tte, small$var_tte, tolerance = 1e-8)
  }
})

# ---------------------------------------------------------------------------
# The reference-slope rules on RCT data
# ---------------------------------------------------------------------------

test_that("rct data without usetrt ignores the observed treated slope", {
  skip_without_paper_data()
  p <- paper_fit("slpower3")

  # Both slopes are estimated and reported -- Stata prints them in the data
  # characteristics block -- but with effectiveness the reference is ZERO, not
  # the treated arm. Stata drops the treated slope from r(table) to make the
  # point. This is the default and it is deliberate.
  eff <- slope_sample_size(p, c(0, 1, 2), effectiveness = 1)
  expect_equal(eff$reference_slope, 0)
  expect_equal(eff$slope_difference, p$slope)
  expect_equal(eff$tte, -p$slope)

  # usetrt measures toward the observed treated slope instead, with
  # effectiveness pinned to 1. Same data, same design, different question.
  obs <- slope_sample_size(p, c(0, 1, 2), target = "observed")
  expect_equal(obs$reference_slope, p$slope_comparator)
  expect_equal(obs$slope_difference, p$slope - p$slope_comparator)

  # The pair differs only in the reference slope, so the sample sizes differ
  # too -- and by a lot, because the treated arm's slope is most of the way to
  # the control arm's.
  expect_gt(obs$n, eff$n)
})

test_that("target = observed reports effectiveness as NA, matching Stata", {
  skip_without_paper_data()
  p <- paper_fit("slpower3")
  r <- slope_sample_size(p, c(0, 1, 2), target = "observed")
  # Stata sets `local effectiveness = .` before building r(table) under usetrt,
  # so the reported effectiveness is missing rather than 1. It is 1 internally
  # in both implementations; what is being pinned is the reporting.
  expect_true(is.na(r$effectiveness))
  expect_equal(r$target, "observed")
  expect_equal(as.data.frame(r)$effectiveness, NA_real_)
})

# ---------------------------------------------------------------------------
# Schedule validity: Stata's numlist rejects these before slopepower sees them
# ---------------------------------------------------------------------------

test_that("visit schedules Stata's numlist refuses are refused here too", {
  # schedule(numlist ascending integer >=1) with baseline 0 implied, so on the
  # R side the equivalent statements are about `visits`. The rc values are the
  # ones Stata actually returned.
  expect_error(trial_design(c(0, 0, 1, 2)))    # rc 124, elements out of order
  expect_error(trial_design(c(0, 1, 2, 2)))    # rc 124, SCHED-repeat
  expect_error(trial_design(c(0, 3, 2, 1)))    # rc 124, SCHED-descending
  expect_error(trial_design(c(0)))             # no follow-up visit at all
  expect_error(trial_design(c(1, 2, 3)))       # must start at baseline 0

  # What Stata cannot express and the port can: a non-integer schedule. Stata
  # returns rc 126 ("noninteger elements") for schedule(1 1.5 2) and can only
  # reach such times through scale(). CONTRACT.md 5.1.
  expect_silent(trial_design(c(0, 1, 1.5, 2)))
})

test_that("scale() has an exact analogue in real-valued visit times", {
  skip_without_paper_data()
  # Stata scale(0.001) schedule(1 2) on slpower1 reports slope -0.0016725 --
  # the year-scale -1.6725 divided by 1000 -- and N = 453940728. The port has
  # no scale argument: it reads the schedule in the units of the fitted time
  # variable, so the same design is visits = c(0, 0.001, 0.002).
  p <- paper_fit("slpower1")
  a <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  b <- slope_sample_size(p, c(0, 0.001, 0.002), effectiveness = 0.33)

  # N is scale-invariant only if the visits move with the axis; here they do
  # not, so the tiny-spacing design needs vastly more people. That asymmetry is
  # the point of dropping scale(): the visits carry the units.
  expect_gt(b$n, a$n)
  expect_equal(b$n, 453940728, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# n handling
# ---------------------------------------------------------------------------

test_that("n is floored to an even number, as Stata does", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")
  # Stata: local actual_n = 2 * floor(`given_n' / 2). n(999) prints
  # "specified N = 999 / actual N = 998" and power 0.913.
  r <- slope_power(p, c(0, 1, 2), n = 999, effectiveness = 0.33)
  expect_equal(r$n_requested, 999)
  expect_equal(r$n, 998)
  expect_equal(r$n_per_arm, 499)
  expect_equal(r$power, 0.913, tolerance = 5e-4)

  # Stata refuses n below 2, and non-integer n ("n must be a whole number
  # greater than or equal to 2", rc 198).
  expect_error(slope_power(p, c(0, 1, 2), n = 1))
  expect_error(slope_power(p, c(0, 1, 2), n = 0))
  expect_error(slope_power(p, c(0, 1, 2), n = 450.5))
  # n = 2 is the smallest Stata accepts: one person per arm, the s* design.
  expect_silent(slope_power(p, c(0, 1, 2), n = 2, effectiveness = 0.33))
})

test_that("alpha, power and effectiveness boundaries match Stata's", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")
  # Stata's messages, all rc 198:
  #   "Alpha must be a number greater than 0 and less than 1"
  #   "Power must be a value greater than 0" / "strictly less than 1"
  #   "Effectiveness must be strictly greater than 0" / "less than or equal 1"
  expect_error(slope_sample_size(p, c(0, 1, 2), alpha = 0))
  expect_error(slope_sample_size(p, c(0, 1, 2), alpha = 1))
  expect_error(slope_sample_size(p, c(0, 1, 2), power = 0))
  expect_error(slope_sample_size(p, c(0, 1, 2), power = 1))
  expect_error(slope_sample_size(p, c(0, 1, 2), effectiveness = 0))
  expect_error(slope_sample_size(p, c(0, 1, 2), effectiveness = 1.0001))

  # effectiveness = 1 is allowed by both: "eliminate the whole slope". Stata
  # gives N = 78 on slpower1 with schedule(1 2).
  r <- slope_sample_size(p, c(0, 1, 2), effectiveness = 1)
  expect_equal(r$n, 78)

  # power just under 1 is allowed and prints as 1.000. Stata: N = 3514.
  r <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33, power = 0.99999)
  expect_equal(r$n, 3514)
})

test_that("the documented defaults are the ones Stata uses", {
  skip_without_paper_data()
  p <- paper_fit("slpower1")
  # Stata defaults effectiveness to 0.25 and, when neither power nor n is
  # given, power to 0.8. The port keeps the effectiveness default in the
  # ordinary functions but keeps the power fallback only in the compatibility
  # wrapper: slope_power() has no default sample size, and slope_sample_size()
  # states power = 0.8 in its own signature.
  expect_equal(formals(slope_sample_size)$effectiveness, 0.25)
  expect_equal(formals(slope_sample_size)$power, 0.8)
  expect_equal(formals(slope_power)$effectiveness, 0.25)
  expect_true(is.symbol(formals(slope_power)$n))   # required, no default
})
