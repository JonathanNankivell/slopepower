# Layer 2 --- the trial design every stage-two function builds from `visits`,
# `dropout` and `dropout_scale`. See CONTRACT.md sections 3 and 6.
#
# Three of these tests pin behaviour that differs deliberately from the Stata
# original, where the corresponding checks are either absent, vacuous or
# documented incorrectly.
#
# The design is built by the one internal validator every stage-two entry point
# calls, so these exercise it directly: every rule here holds identically for
# slope_sample_size(), slope_power(), slope_effect_size() and the bootstraps.
# A few tests at the end confirm the route through the public functions.

design_of <- function(visits, dropout = NULL,
                      dropout_scale = c("incremental", "cumulative")) {
  slopepower:::build_trial_design(visits, dropout, match.arg(dropout_scale),
                                  "slope_sample_size()")
}

test_that("the design returns exactly the contract fields", {
  d <- design_of(c(0, 1, 2))
  expect_s3_class(d, "trial_design")
  expect_setequal(names(d), c("visits", "dropout", "has_dropout", "dropout_scale"))
  expect_length(names(d), 4L)
})

test_that("a design without dropout has a zero dropout vector", {
  d <- design_of(c(0, 1, 2, 5))
  expect_equal(d$visits, c(0, 1, 2, 5))
  expect_equal(d$dropout, rep(0, 3))
  expect_false(d$has_dropout)
  expect_identical(d$dropout_scale, "incremental")
})

test_that("dropout is stored with one element per follow-up visit", {
  d <- design_of(c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
  expect_length(d$dropout, length(d$visits) - 1L)
  expect_equal(d$dropout, c(0, 0, 0.1))
  expect_true(d$has_dropout)
})

test_that("visit times may be any real values (no scale() equivalent needed)", {
  # Stata restricts schedule() to ascending integers >= 1 and supplies scale()
  # to compensate. This port builds the covariance at the requested times.
  d <- design_of(seq(0, 3, by = 0.5))
  expect_equal(d$visits, c(0, 0.5, 1, 1.5, 2, 2.5, 3))
  expect_length(d$dropout, 6L)
})

# --- cumulative / incremental equivalence ----------------------------------

test_that("cumulative dropout converts to the same design as incremental", {
  cumulative  <- suppressWarnings(
    design_of(c(0, 1, 2, 3), dropout = c(0.05, 0.10, 0.15),
                 dropout_scale = "cumulative"))
  incremental <- suppressWarnings(
    design_of(c(0, 1, 2, 3), dropout = c(0.05, 0.05, 0.05)))

  expect_equal(cumulative$dropout, incremental$dropout)
  expect_equal(cumulative$visits, incremental$visits)
  # only the recorded provenance differs
  expect_identical(cumulative$dropout_scale, "cumulative")
  expect_identical(incremental$dropout_scale, "incremental")
})

test_that("cumulative dropout must be non-decreasing", {
  expect_error(
    design_of(c(0, 1, 2), dropout = c(0.2, 0.1), dropout_scale = "cumulative"),
    "non-decreasing"
  )
})

test_that("cumulative dropout cannot exceed 1", {
  expect_error(
    design_of(c(0, 1, 2), dropout = c(0.5, 1.4), dropout_scale = "cumulative"),
    "cannot exceed 1"
  )
})

# --- dropout_rate() on the ordinary (non-grid) path -------------------------
#
# dropout_rate() was originally defined in grid.R and expanded only by
# grid_impl(), so a single design with dropout_rate(0.05) -- the call
# dropout_rate()'s own documentation pointed at -- failed with "`dropout` must
# be numeric or NULL; got dropout_rate". These pin the design validator as the
# place the expansion happens, which is what makes the grid's use of it a
# special case rather than the only case.

test_that("the design expands a dropout_rate() to per-interval proportions", {
  d <- suppressWarnings(design_of(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  expect_equal(d$dropout, rep(0.05, 3))
  expect_true(d$has_dropout)
  # Expansion does not change the provenance the object records: the vector it
  # produced is incremental, which is what "incremental" already means here.
  expect_identical(d$dropout_scale, "incremental")
})

test_that("the same rate expands differently for different schedules", {
  # The whole point of the object: 5% per unit time over 3 units is the same
  # total dropout however many visits it is spread across.
  annual <- suppressWarnings(design_of(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  halves <- suppressWarnings(design_of(seq(0, 3, by = 0.5), dropout = dropout_rate(0.05)))
  single <- suppressWarnings(design_of(c(0, 3), dropout = dropout_rate(0.05)))

  expect_equal(annual$dropout, rep(0.05, 3))
  expect_equal(halves$dropout, rep(0.025, 6))
  expect_equal(single$dropout, 0.15)
  expect_equal(sum(annual$dropout), sum(halves$dropout))
  expect_equal(sum(annual$dropout), sum(single$dropout))
})

test_that("`per` scales the rate to the units of the visit times", {
  d <- suppressWarnings(design_of(c(0, 6, 12, 24), dropout = dropout_rate(0.10, per = 12)))
  expect_equal(d$dropout, c(0.05, 0.05, 0.10))
})

test_that("a zero rate yields a design with no dropout", {
  d <- design_of(c(0, 1, 2), dropout = dropout_rate(0))
  expect_equal(d$dropout, c(0, 0))
  expect_false(d$has_dropout)
})

test_that("a rate implying total dropout above 1 errors in terms of the rate", {
  err <- expect_error(design_of(c(0, 1, 2, 3), dropout = dropout_rate(0.5)),
                      "exceeds 1")
  # Diagnosed by `rate` and `per`, not by the expanded vector: naming a vector
  # the caller never wrote would send them looking for the wrong mistake.
  expect_match(conditionMessage(err), "rate of 0.5 per 1")
  expect_match(conditionMessage(err), "lasting 3")
})

test_that("a dropout_rate() cannot be combined with dropout_scale = cumulative", {
  err <- expect_error(
    design_of(c(0, 1, 2, 3), dropout = dropout_rate(0.05),
                 dropout_scale = "cumulative"),
    "cannot be combined")
  expect_match(conditionMessage(err), "already incremental")
})

# --- dropout_rate(pattern = "geometric") ------------------------------------
#
# The alternative to the default "linear" model: `rate` is the proportion of
# whoever is still in follow-up who withdraws per unit time, compounding
# geometrically, rather than a proportion of the original cohort.

test_that("a geometric rate applies to whoever remains, not the original cohort", {
  d <- suppressWarnings(design_of(c(0, 1, 2, 3),
                                     dropout = dropout_rate(0.05, pattern = "geometric")))
  # 5% of 1, then 5% of the 0.95 remaining, then 5% of the 0.9025 remaining.
  expect_equal(d$dropout, c(0.05, 0.05 * 0.95, 0.05 * 0.95^2))
  # Each interval's survivors are 95% of the last, so the total is strictly
  # less than a linear rate at the same `rate` and `per` would give.
  linear <- suppressWarnings(design_of(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  expect_lt(sum(d$dropout), sum(linear$dropout))
})

test_that("a geometric rate of 1 drops everyone after the first interval", {
  d <- suppressWarnings(design_of(c(0, 1, 2, 3),
                                     dropout = dropout_rate(1, pattern = "geometric")))
  expect_equal(d$dropout, c(1, 0, 0))
})

test_that("a geometric rate above 1 is rejected as not a proportion of survivors", {
  err <- expect_error(dropout_rate(1.5, pattern = "geometric"), "must be in \\[0, 1\\]")
  expect_match(conditionMessage(err), "1.5")
})

test_that("a geometric rate's total never exceeds 1, however long the trial", {
  # Unlike the linear model, geometric decay can approach but never pass zero
  # survival, so there is no rate/duration combination for which this errors.
  # (At extreme rate/duration combinations, survival underflows to exactly 0
  # in double precision and the total lands on exactly 1 rather than below
  # it -- still not over, which is the bound that matters.)
  d <- suppressWarnings(design_of(c(0, 100), dropout = dropout_rate(0.9, pattern = "geometric")))
  expect_lte(sum(d$dropout), 1)
})

test_that("`per` scales a geometric rate the same way it scales a linear one", {
  d <- suppressWarnings(design_of(c(0, 6, 12, 24),
                                     dropout = dropout_rate(0.10, per = 12, pattern = "geometric")))
  survival <- (1 - 0.10)^(c(0, 6, 12, 24) / 12)
  expect_equal(d$dropout, -diff(survival))
})

test_that("print.dropout_rate() names the pattern", {
  expect_output(print(dropout_rate(0.05)), "linear")
  expect_output(print(dropout_rate(0.05, pattern = "geometric")), "geometric")
})

test_that("a rate and the vector it expands to give the same design", {
  rate <- suppressWarnings(design_of(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  hand <- suppressWarnings(design_of(c(0, 1, 2, 3), dropout = rep(0.05, 3)))
  expect_equal(rate, hand)
})

test_that("a rate does not warn about the baseline-only stratum it always has", {
  # Every non-zero rate puts someone in the baseline-only stratum, so the
  # warning would only say that a rate is a rate (see ?dropout_rate). The same
  # proportions written out as a vector still warn: that first element was
  # typed by someone.
  expect_silent(d <- design_of(c(0, 1, 2), dropout = dropout_rate(0.05)))
  expect_gt(d$dropout[1L], 0)
  expect_warning(design_of(c(0, 1, 2), dropout = d$dropout), "contribute nothing")
})

# --- the Stata-form checks this layer restates or repairs --------------------

test_that("dropout summing to exactly 1 is accepted", {
  # 1 - .3 - .3 - .4 is -5.551e-17, so a bare `< 0` guard would reject
  # c(0.3, 0.3, 0.4) -- legal, and exactly 1 in decimal. Summing once and
  # comparing against a tolerance is what this layer does instead.
  #
  # Stata reaches the same answer by a different route, which is why this is a
  # restatement and not a repair: slopepower.ado:265 does accumulate by
  # subtraction, but a Stata local round-trips through a decimal string (.7 -
  # .3 stores as ".4"), so its residue is exactly 0. It accepts the list and
  # returns N = 1418 -- see the end-to-end pin in test-stata-behaviour.R.
  d <- suppressWarnings(design_of(c(0, 1, 2, 3), dropout = c(0.3, 0.3, 0.4)))
  expect_equal(sum(d$dropout), 1)
  expect_equal(1 - sum(d$dropout), 0)
})

test_that("a dropout vector of the wrong length errors", {
  # slopepower.ado:239,276 build the length counters with a space inside the
  # macro name, which looked like it would make the guard at :279 dead code.
  # It does not: Stata trims the name and the guard fires on a list that is
  # short or long (stata-reference/00_open_questions.log, and the pin in
  # test-stata-behaviour.R). So this restates Stata's check rather than
  # repairing it -- but it is the check every stage-two path depends on, so it
  # is pinned here too.
  expect_error(
    design_of(c(0, 1, 2), dropout = c(0.05, 0.05, 0.05)),
    "one element per follow-up visit"
  )
  expect_error(
    design_of(c(0, 1, 2, 5), dropout = c(0.05)),
    "one element per follow-up visit"
  )
})

test_that("incremental dropout summing above 1 errors and suggests cumulative", {
  expect_error(design_of(c(0, 1, 2), dropout = c(0.5, 0.6)), "exceeds 1")
  expect_error(design_of(c(0, 1, 2), dropout = c(0.5, 0.6)), "cumulative")
})

# --- baseline is explicit ---------------------------------------------------

test_that("visits must begin at baseline 0, with a corrected call suggested", {
  err <- expect_error(design_of(c(1, 2)), "baseline visit at time 0")
  expect_match(conditionMessage(err), "c(0, 1, 2)", fixed = TRUE)
})

test_that("negative visit times error", {
  expect_error(design_of(c(-1, 0, 1)), "baseline visit at time 0")
})

# --- other guards -----------------------------------------------------------

test_that("duplicate, unsorted, short and non-finite visits error", {
  expect_error(design_of(c(0, 1, 1)), "repeated times")
  expect_error(design_of(c(0, 2, 1)), "increasing order")
  expect_error(design_of(0), "at least 2 visit times")
  expect_error(design_of(c(0, NA, 2)), "finite")
  expect_error(design_of(c(0, 1, Inf)), "finite")
})

test_that("negative and non-numeric dropout error", {
  expect_error(design_of(c(0, 1, 2), dropout = c(-0.1, 0.2)), "non-negative")
  expect_error(design_of(c(0, 1, 2), dropout = c("a", "b")), "numeric")
  expect_error(design_of(c(0, 1, 2), dropout = c(0.1, NA)), "finite")
})

test_that("dropout before the first follow-up visit warns", {
  # Those participants attend baseline only and carry no slope information.
  # Stata skips this stratum silently (ado:581, :613).
  expect_warning(design_of(c(0, 2, 3), dropout = c(0.2, 0.1)),
                 "contribute nothing")
  expect_silent(design_of(c(0, 2, 3), dropout = c(0, 0.1)))
})

test_that("a design prints with and without dropout", {
  expect_output(print(design_of(c(0, 1, 2))), "trial_design")
  expect_output(print(design_of(c(0, 1, 2))), "none")
  d <- suppressWarnings(design_of(c(0, 2, 3), dropout = c(0.2, 0.1)))
  out <- capture.output(print(d))
  expect_true(any(grepl("last visit", out)))
  expect_true(any(grepl("first missed", out)))   # both readings shown
  expect_true(any(grepl("Completers", out)))
})

# --- the route through the public functions ---------------------------------

test_that("every stage-two function reads `dropout` through `dropout_scale`", {
  # One cumulative vector and the incremental vector it means, priced by each
  # entry point: the answers must agree exactly, because the conversion happens
  # once, in the validator, before any calculation sees the vector.
  p <- ref_params()
  v <- c(0, 1, 2, 3)
  cum <- c(0, 0.05, 0.15)
  inc <- c(0, 0.05, 0.10)
  expect_identical(
    slope_sample_size(p, v, cum, dropout_scale = "cumulative")$n,
    slope_sample_size(p, v, inc)$n)
  expect_equal(
    slope_power(p, v, cum, dropout_scale = "cumulative", n = 400)$power,
    slope_power(p, v, inc, n = 400)$power)
  expect_equal(
    slope_effect_size(p, v, cum, dropout_scale = "cumulative"),
    slope_effect_size(p, v, inc))
  expect_identical(
    slope_sample_size(p, v, cum, dropout_scale = "cumulative")$design$dropout_scale,
    "cumulative")
})

test_that("a stage-two function names itself in a design error", {
  p <- ref_params()
  expect_error(slope_power(p, c(1, 2), n = 100), "^slope_power\\(\\): `visits` must begin")
  expect_error(slope_effect_size(p, c(0, 1, 2), dropout = 0.1),
               "^slope_effect_size\\(\\): `dropout` must have one element")
  expect_error(slope_sample_size(p, c(0, 1, 2), dropout = dropout_rate(0.05),
                                 dropout_scale = "cumulative"),
               "cannot be combined with dropout_scale")
})
