# Layer 4 --- slopepower(), the Stata-style compatibility wrapper.
#
# Its purpose is to make the paper's published commands directly translatable
# and to serve as a parity harness, so these tests mirror the Stata syntax of
# each worked example as closely as the R interface allows.

test_that("slopepower() reproduces the p.588 single-group example", {
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  # Stata: slopepower sdmt, schedule(1 2) subject(id) time(visit) obs
  #        nocontrols effectiveness(0.33)
  r <- suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = c(1, 2), obs = TRUE, nocontrols = TRUE, effectiveness = 0.33))
  expect_s3_class(r, "slope_power")
  expect_equal(r$n, 712)
  expect_equal(r$n_per_arm, 356)
})

test_that("schedule() lists follow-up visits only, baseline implied", {
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  # Stata: schedule(1 2 5) dropouts(0 0 0.1) -> N = 328
  r <- suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = c(1, 2, 5), dropouts = c(0, 0, 0.1),
    obs = TRUE, nocontrols = TRUE, effectiveness = 0.33))
  expect_equal(r$n, 328)
  # the design carries the implicit baseline
  expect_equal(r$design$visits, c(0, 1, 2, 5))
})

test_that("slopepower() supports the scale() option", {
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  # Stata: schedule(1 2 3 4) scale(0.5) -> N = 620, on a half-year axis
  r <- suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = 1:4, scale = 0.5, obs = TRUE, nocontrols = TRUE,
    effectiveness = 0.33))
  expect_equal(r$n, 620)
  expect_lt(abs(r$tte - 0.276), 5e-4)
})

test_that("slopepower() reproduces the p.590 case-control example", {
  skip_without_paper_data()
  d <- load_paper_data("slpower2")
  # Stata: slopepower sdmt, schedule(1 2) scale(365) subject(id) time(vdate)
  #        obs casecon(case) effectiveness(0.33)
  r <- suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = c(1, 2), obs = TRUE, casecon = "case", effectiveness = 0.33))
  expect_equal(r$n, 296)
  expect_equal(r$n_per_arm, 148)
})

test_that("slopepower() computes power when n is given (p.593)", {
  skip_without_paper_data()
  d <- load_paper_data("slpower2")
  # the baseline-only dropout warning is asserted in test-design.R; incidental here
  r <- suppressWarnings(suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = c(1, 2), dropouts = c(0.05, 0.05), obs = TRUE,
    casecon = "case", effectiveness = 0.33, n = 200)))
  expect_lt(abs(r$power - 0.597), 5e-4)
  expect_equal(r$n, 200)
})

test_that("slopepower() reproduces the p.594 usetrt example", {
  skip_without_paper_data()
  d <- load_paper_data("slpower3")
  # Stata: slopepower sdmt, schedule(2 3) subject(id) time(visit) rct
  #        treat(treat) usetrt dropout(0.2 0.1)
  r <- suppressWarnings(suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = c(2, 3), dropouts = c(0.2, 0.1), rct = TRUE,
    treat = "treat", usetrt = TRUE)))
  expect_equal(r$n, 318)
  expect_equal(r$n_per_arm, 159)
  expect_true(is.na(r$effectiveness))
})

# --- Stata's model-selection rules ------------------------------------------

test_that("slopepower() enforces the obs/rct model selection rules", {
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  base <- list(d, depvar = "sdmt", subject = "id", time = "time",
               schedule = c(1, 2))
  expect_error(do.call(slopepower, c(base, list(obs = TRUE, rct = TRUE))),
               "both")
  expect_error(do.call(slopepower, base), "no model specified")
  expect_error(
    do.call(slopepower, c(base, list(rct = TRUE, nocontrols = TRUE))),
    "nocontrols")
  # observational data with controls needs casecon
  expect_error(do.call(slopepower, c(base, list(obs = TRUE))), "casecon")
  # rct needs treat
  expect_error(do.call(slopepower, c(base, list(rct = TRUE))), "treat")
})

test_that("slopepower() warns when an option does not apply to the model", {
  skip_without_paper_data()
  d <- load_paper_data("slpower3")
  # usetrt outside an RCT is ignored with a warning, as in Stata
  expect_warning(suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
    obs = TRUE, nocontrols = TRUE, effectiveness = 0.33, usetrt = TRUE)),
    "usetrt")
})

test_that("slopepower() rejects effectiveness together with usetrt", {
  skip_without_paper_data()
  d <- load_paper_data("slpower3")
  expect_error(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time", schedule = c(2, 3),
    rct = TRUE, treat = "treat", usetrt = TRUE, effectiveness = 0.5),
    "only one of")
})

test_that("slopepower() validates column names", {
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  expect_error(slopepower(d, depvar = "nope", subject = "id", time = "time",
                          schedule = c(1, 2), obs = TRUE, nocontrols = TRUE),
               "not a column")
  expect_error(slopepower(d, depvar = c("a", "b"), subject = "id", time = "time",
                          schedule = c(1, 2), obs = TRUE, nocontrols = TRUE),
               "single column name")
})

test_that("slopepower() rejects a non-positive schedule, as Stata does", {
  skip_without_paper_data()
  d <- load_paper_data("slpower1")
  # Stata declares schedule(numlist ascending integer >=1); baseline is implicit
  # and must not be listed
  expect_error(slopepower(d, depvar = "sdmt", subject = "id", time = "time",
                          schedule = c(0, 1, 2), obs = TRUE, nocontrols = TRUE,
                          effectiveness = 0.33))
})
