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
