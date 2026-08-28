# Layer 2 --- trial_design(). See CONTRACT.md sections 3 and 6.
#
# Three of these tests pin behaviour that differs deliberately from the Stata
# original, where the corresponding checks are either absent, vacuous or
# documented incorrectly.

test_that("trial_design() returns exactly the contract fields", {
  d <- trial_design(c(0, 1, 2))
  expect_s3_class(d, "trial_design")
  expect_setequal(names(d), c("visits", "dropout", "has_dropout", "dropout_type"))
  expect_length(names(d), 4L)
})

test_that("a design without dropout has a zero dropout vector", {
  d <- trial_design(c(0, 1, 2, 5))
  expect_equal(d$visits, c(0, 1, 2, 5))
  expect_equal(d$dropout, rep(0, 3))
  expect_false(d$has_dropout)
  expect_identical(d$dropout_type, "incremental")
})

test_that("dropout is stored with one element per follow-up visit", {
  d <- trial_design(c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
  expect_length(d$dropout, length(d$visits) - 1L)
  expect_equal(d$dropout, c(0, 0, 0.1))
  expect_true(d$has_dropout)
})

test_that("visit times may be any real values (no scale() equivalent needed)", {
  # Stata restricts schedule() to ascending integers >= 1 and supplies scale()
  # to compensate. This port builds the covariance at the requested times.
  d <- trial_design(seq(0, 3, by = 0.5))
  expect_equal(d$visits, c(0, 0.5, 1, 1.5, 2, 2.5, 3))
  expect_length(d$dropout, 6L)
})

# --- cumulative / incremental equivalence ----------------------------------

test_that("cumulative dropout converts to the same design as incremental", {
  cumulative  <- suppressWarnings(
    trial_design(c(0, 1, 2, 3), dropout = c(0.05, 0.10, 0.15),
                 dropout_type = "cumulative"))
  incremental <- suppressWarnings(
    trial_design(c(0, 1, 2, 3), dropout = c(0.05, 0.05, 0.05)))

  expect_equal(cumulative$dropout, incremental$dropout)
  expect_equal(cumulative$visits, incremental$visits)
  # only the recorded provenance differs
  expect_identical(cumulative$dropout_type, "cumulative")
  expect_identical(incremental$dropout_type, "incremental")
})

test_that("cumulative dropout must be non-decreasing", {
  expect_error(
    trial_design(c(0, 1, 2), dropout = c(0.2, 0.1), dropout_type = "cumulative"),
    "non-decreasing"
  )
})

test_that("cumulative dropout cannot exceed 1", {
  expect_error(
    trial_design(c(0, 1, 2), dropout = c(0.5, 1.4), dropout_type = "cumulative"),
    "cannot exceed 1"
  )
})

# --- dropout_rate() on the ordinary (non-grid) path -------------------------
#
# dropout_rate() was originally defined in grid.R and expanded only by
# grid_impl(), so trial_design(v, dropout_rate(0.05)) -- the call dropout_rate()'s
# own documentation pointed at -- failed with "`dropout` must be numeric or
# NULL; got dropout_rate". These pin the constructor as the place the expansion
# happens, which is what makes the grid's use of it a special case rather than
# the only case.

test_that("trial_design() expands a dropout_rate() to per-interval proportions", {
  d <- suppressWarnings(trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  expect_equal(d$dropout, rep(0.05, 3))
  expect_true(d$has_dropout)
  # Expansion does not change the provenance the object records: the vector it
  # produced is incremental, which is what "incremental" already means here.
  expect_identical(d$dropout_type, "incremental")
})

test_that("the same rate expands differently for different schedules", {
  # The whole point of the object: 5% per unit time over 3 units is the same
  # total dropout however many visits it is spread across.
  annual <- suppressWarnings(trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  halves <- suppressWarnings(trial_design(seq(0, 3, by = 0.5), dropout = dropout_rate(0.05)))
  single <- suppressWarnings(trial_design(c(0, 3), dropout = dropout_rate(0.05)))

  expect_equal(annual$dropout, rep(0.05, 3))
  expect_equal(halves$dropout, rep(0.025, 6))
  expect_equal(single$dropout, 0.15)
  expect_equal(sum(annual$dropout), sum(halves$dropout))
  expect_equal(sum(annual$dropout), sum(single$dropout))
})

test_that("`per` scales the rate to the units of the visit times", {
  d <- suppressWarnings(trial_design(c(0, 6, 12, 24), dropout = dropout_rate(0.10, per = 12)))
  expect_equal(d$dropout, c(0.05, 0.05, 0.10))
})

test_that("a zero rate yields a design with no dropout", {
  d <- trial_design(c(0, 1, 2), dropout = dropout_rate(0))
  expect_equal(d$dropout, c(0, 0))
  expect_false(d$has_dropout)
})

test_that("a rate implying total dropout above 1 errors in terms of the rate", {
  err <- expect_error(trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.5)),
                      "exceeds 1")
  # Diagnosed by `rate` and `per`, not by the expanded vector: naming a vector
  # the caller never wrote would send them looking for the wrong mistake.
  expect_match(conditionMessage(err), "rate of 0.5 per 1")
  expect_match(conditionMessage(err), "lasting 3")
})

test_that("a dropout_rate() cannot be combined with dropout_type = cumulative", {
  err <- expect_error(
    trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.05),
                 dropout_type = "cumulative"),
    "cannot be combined")
  expect_match(conditionMessage(err), "already incremental")
})

# --- dropout_rate(type = "cumulative") ---------------------------------------
#
# The alternative to the default "linear" model: `rate` is the proportion of
# whoever is still in follow-up who withdraws per unit time, compounding
# geometrically, rather than a proportion of the original cohort.

test_that("a cumulative rate applies to whoever remains, not the original cohort", {
  d <- suppressWarnings(trial_design(c(0, 1, 2, 3),
                                     dropout = dropout_rate(0.05, type = "cumulative")))
  # 5% of 1, then 5% of the 0.95 remaining, then 5% of the 0.9025 remaining.
  expect_equal(d$dropout, c(0.05, 0.05 * 0.95, 0.05 * 0.95^2))
  # Each interval's survivors are 95% of the last, so the total is strictly
  # less than a linear rate at the same `rate` and `per` would give.
  linear <- suppressWarnings(trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  expect_lt(sum(d$dropout), sum(linear$dropout))
})

test_that("a cumulative rate of 1 drops everyone after the first interval", {
  d <- suppressWarnings(trial_design(c(0, 1, 2, 3),
                                     dropout = dropout_rate(1, type = "cumulative")))
  expect_equal(d$dropout, c(1, 0, 0))
})

test_that("a cumulative rate above 1 is rejected as not a proportion of survivors", {
  err <- expect_error(dropout_rate(1.5, type = "cumulative"), "must be in \\[0, 1\\]")
  expect_match(conditionMessage(err), "1.5")
})

test_that("a cumulative rate's total never exceeds 1, however long the trial", {
  # Unlike the linear model, geometric decay can approach but never pass zero
  # survival, so there is no rate/duration combination for which this errors.
  # (At extreme rate/duration combinations, survival underflows to exactly 0
  # in double precision and the total lands on exactly 1 rather than below
  # it -- still not over, which is the bound that matters.)
  d <- suppressWarnings(trial_design(c(0, 100), dropout = dropout_rate(0.9, type = "cumulative")))
  expect_lte(sum(d$dropout), 1)
})

test_that("`per` scales a cumulative rate the same way it scales a linear one", {
  d <- suppressWarnings(trial_design(c(0, 6, 12, 24),
                                     dropout = dropout_rate(0.10, per = 12, type = "cumulative")))
  survival <- (1 - 0.10)^(c(0, 6, 12, 24) / 12)
  expect_equal(d$dropout, -diff(survival))
})

test_that("print.dropout_rate() names the type", {
  expect_output(print(dropout_rate(0.05)), "linear")
  expect_output(print(dropout_rate(0.05, type = "cumulative")), "cumulative")
})

test_that("a rate and the vector it expands to give the same design", {
  rate <- suppressWarnings(trial_design(c(0, 1, 2, 3), dropout = dropout_rate(0.05)))
  hand <- suppressWarnings(trial_design(c(0, 1, 2, 3), dropout = rep(0.05, 3)))
  expect_equal(rate, hand)
})

test_that("a rate warns about the baseline-only stratum like any other dropout", {
  # Not suppressed here: every non-zero rate puts someone in the baseline-only
  # stratum, and on this path -- unlike the grid, which collects the warning and
  # reports it once per table -- the user is told each time.
  expect_warning(trial_design(c(0, 1, 2), dropout = dropout_rate(0.05)),
                 "contribute nothing")
})

test_that("a dropout_rate() in a hand-built design is refused with the fix", {
  # `trial_design()` expands on construction, so a rate surviving into the
  # stored field means the object never went through the constructor. Expanding
  # it at the stage-two boundary would leave `has_dropout` describing something
  # else, so as_trial_design() names the constructor instead.
  d <- structure(list(visits = c(0, 1, 2, 3), dropout = dropout_rate(0.05),
                      has_dropout = TRUE, dropout_type = "incremental"),
                 class = "trial_design")
  err <- expect_error(slope_sample_size(ref_params(), d), "must be numeric")
  expect_match(conditionMessage(err), "trial_design(visits = c(0, 1, 2, 3)", fixed = TRUE)
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
  d <- suppressWarnings(trial_design(c(0, 1, 2, 3), dropout = c(0.3, 0.3, 0.4)))
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
    trial_design(c(0, 1, 2), dropout = c(0.05, 0.05, 0.05)),
    "one element per follow-up visit"
  )
  expect_error(
    trial_design(c(0, 1, 2, 5), dropout = c(0.05)),
    "one element per follow-up visit"
  )
})

test_that("incremental dropout summing above 1 errors and suggests cumulative", {
  expect_error(trial_design(c(0, 1, 2), dropout = c(0.5, 0.6)), "exceeds 1")
  expect_error(trial_design(c(0, 1, 2), dropout = c(0.5, 0.6)), "cumulative")
})

# --- baseline is explicit ---------------------------------------------------

test_that("visits must begin at baseline 0, with a corrected call suggested", {
  err <- expect_error(trial_design(c(1, 2)), "baseline visit at time 0")
  expect_match(conditionMessage(err), "c(0, 1, 2)", fixed = TRUE)
})

test_that("negative visit times error", {
  expect_error(trial_design(c(-1, 0, 1)), "baseline visit at time 0")
})

# --- other guards -----------------------------------------------------------

test_that("duplicate, unsorted, short and non-finite visits error", {
  expect_error(trial_design(c(0, 1, 1)), "repeated times")
  expect_error(trial_design(c(0, 2, 1)), "increasing order")
  expect_error(trial_design(0), "at least 2 visit times")
  expect_error(trial_design(c(0, NA, 2)), "finite")
  expect_error(trial_design(c(0, 1, Inf)), "finite")
})

test_that("negative and non-numeric dropout error", {
  expect_error(trial_design(c(0, 1, 2), dropout = c(-0.1, 0.2)), "non-negative")
  expect_error(trial_design(c(0, 1, 2), dropout = c("a", "b")), "numeric")
  expect_error(trial_design(c(0, 1, 2), dropout = c(0.1, NA)), "finite")
})

test_that("dropout before the first follow-up visit warns", {
  # Those participants attend baseline only and carry no slope information.
  # Stata skips this stratum silently (ado:581, :613).
  expect_warning(trial_design(c(0, 2, 3), dropout = c(0.2, 0.1)),
                 "contribute nothing")
  expect_silent(trial_design(c(0, 2, 3), dropout = c(0, 0.1)))
})

test_that("print.trial_design() runs for designs with and without dropout", {
  expect_output(print(trial_design(c(0, 1, 2))), "trial_design")
  expect_output(print(trial_design(c(0, 1, 2))), "none")
  d <- suppressWarnings(trial_design(c(0, 2, 3), dropout = c(0.2, 0.1)))
  out <- capture.output(print(d))
  expect_true(any(grepl("last visit", out)))
  expect_true(any(grepl("first missed", out)))   # both readings shown
  expect_true(any(grepl("Completers", out)))
})
