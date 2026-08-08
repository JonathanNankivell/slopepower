# Layer 4 --- slopepower(), the Stata-style compatibility wrapper.
#
# Its purpose is to make the paper's published commands directly translatable
# and to serve as a parity harness, so these tests mirror the Stata syntax of
# each worked example as closely as the R interface allows.

test_that("slopepower() reproduces the p.588 single-group example", {
  d <- load_paper_data("slpower1")
  # Stata: slopepower sdmt, schedule(1 2) subject(id) time(visit) obs
  #        nocontrols effectiveness(0.33)
  r <- suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = c(1, 2), obs = TRUE, nocontrols = TRUE, effectiveness = 0.33))
  # no n given, so the Stata wrapper asks the sample-size question
  expect_s3_class(r, "slope_sample_size")
  expect_s3_class(r, "slope_result")
  expect_equal(r$n, 712)
  expect_equal(r$n_per_arm, 356)
})

test_that("schedule() lists follow-up visits only, baseline implied", {
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
  d <- load_paper_data("slpower2")
  # the baseline-only dropout warning is asserted in test-design.R; incidental here
  r <- suppressWarnings(suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = c(1, 2), dropouts = c(0.05, 0.05), obs = TRUE,
    casecon = "case", effectiveness = 0.33, n = 200)))
  # n given, so it asks the power question
  expect_s3_class(r, "slope_power")
  expect_lt(abs(r$power - 0.597), 5e-4)
  expect_equal(r$n, 200)
})

test_that("slopepower() keeps Stata's single bimodal interface, and guards it", {
  # The Stata command takes n() or power(), never both, and picks the
  # calculation from which was supplied. This wrapper mirrors that; the split
  # into slope_sample_size() and slope_power() is for new code, not for parity.
  d <- load_paper_data("slpower1")
  base <- list(d, depvar = "sdmt", subject = "id", time = "time",
               schedule = c(1, 2), obs = TRUE, nocontrols = TRUE,
               effectiveness = 0.33)
  expect_error(
    do.call(slopepower, c(base, list(n = 450, power = 0.9))),
    "only one of")
  # neither given defaults to power 0.8, as documented
  r <- suppressMessages(do.call(slopepower, base))
  expect_s3_class(r, "slope_sample_size")
  expect_equal(r$power, 0.8)
  expect_equal(r$n, 712)
  # and an explicit power reaches the same place
  r9 <- suppressMessages(do.call(slopepower, c(base, list(power = 0.9))))
  expect_s3_class(r9, "slope_sample_size")
  expect_gt(r9$n, r$n)
})

test_that("slopepower() reproduces the p.594 usetrt example", {
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
  d <- load_paper_data("slpower3")
  # usetrt outside an RCT is ignored with a warning, as in Stata
  expect_warning(suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
    obs = TRUE, nocontrols = TRUE, effectiveness = 0.33, usetrt = TRUE)),
    "usetrt")
})

test_that("slopepower() warns and then really ignores casecon off model 1", {
  # slopepower.ado:131-133 warns whenever casecon() is given with a model other
  # than 1 and then carries on without it -- the warning is not followed by an
  # `exit`, and nothing downstream reads the macro. The port must do the same,
  # so the assertion that matters is not the warning but the equality below: a
  # warning that left the option in force would be the real defect.
  d3 <- load_paper_data("slpower3")

  # Model 3. This is the reference grid's IGN-casecon-model3 row, casecon(treat)
  # and all, listed in KNOWN_DIVERGENCES only because Stata's r(table) has no
  # column saying "ignored" for the comparison to check. Stata returned N = 624.
  expect_warning(
    ign <- suppressMessages(slopepower(
      d3, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
      rct = TRUE, treat = "treat", casecon = "treat", effectiveness = 0.33)),
    "casecon")
  plain <- suppressMessages(slopepower(
    d3, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
    rct = TRUE, treat = "treat", effectiveness = 0.33))
  expect_equal(ign, plain)
  expect_equal(ign$n, 624)
  # The fit is still the RCT one: casecon did not quietly become the grouping.
  expect_equal(ign$params$comparator, "treated")

  # Model 2, the other half of `model != 1`. Declaring nocontrols on data that
  # does have controls is exactly the mistake the warning is for: the healthy
  # subjects stay in the model, pooled with the cases, which is why Stata
  # refuses to let casecon() pass silently.
  d2 <- load_paper_data("slpower2")
  expect_warning(
    ign2 <- suppressMessages(slopepower(
      d2, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
      obs = TRUE, nocontrols = TRUE, casecon = "case", effectiveness = 0.33)),
    "casecon")
  plain2 <- suppressMessages(slopepower(
    d2, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
    obs = TRUE, nocontrols = TRUE, effectiveness = 0.33))
  expect_equal(ign2, plain2)
  expect_equal(ign2$params$comparator, "none")
  # the controls are still in the model, pooled, rather than split off
  expect_equal(ign2$params$n_subjects, 500L)
})

test_that("slopepower() warns and then really ignores treat off an RCT", {
  # slopepower.ado:195-197 means to warn and ignore here as it does for casecon,
  # but the message is written "...RCT data`. Treatment variable..." -- a stray
  # backtick opens a macro reference that never closes, so Stata aborts instead
  # of warning. That is one of the three original defects the port declines to
  # reproduce (see the "Differences from the Stata command" section of
  # ?slopepower); the intended behaviour is the one implemented, so there is no
  # reference-grid row to compare against, only the .ado's evident intent.
  d <- load_paper_data("slpower3")
  expect_warning(
    ign <- suppressMessages(slopepower(
      d, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
      obs = TRUE, nocontrols = TRUE, treat = "treat", effectiveness = 0.33)),
    "treat")
  plain <- suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
    obs = TRUE, nocontrols = TRUE, effectiveness = 0.33))
  expect_equal(ign, plain)

  # Ignoring treat() is not cosmetic on these data: with the arms pooled the
  # single fitted slope is roughly the average of the two, well away from the
  # -1.852 the same data give when treat is honoured (p.594). If the option had
  # survived the warning, this is the number that would have moved.
  expect_equal(ign$params$comparator, "none")
  expect_true(is.na(ign$params$slope_comparator))
  expect_lt(abs(ign$params$slope - (-1.478)), 5e-4)
})

test_that("slopepower() rejects effectiveness together with usetrt", {
  d <- load_paper_data("slpower3")
  expect_error(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time", schedule = c(2, 3),
    rct = TRUE, treat = "treat", usetrt = TRUE, effectiveness = 0.5),
    "only one of")
})

test_that("slopepower() validates column names", {
  d <- load_paper_data("slpower1")
  expect_error(slopepower(d, depvar = "nope", subject = "id", time = "time",
                          schedule = c(1, 2), obs = TRUE, nocontrols = TRUE),
               "not a column")
  expect_error(slopepower(d, depvar = c("a", "b"), subject = "id", time = "time",
                          schedule = c(1, 2), obs = TRUE, nocontrols = TRUE),
               "single column name")
})

# --- the schedule: Stata's numlist, and the one restriction lifted -----------
#
# Stata declares schedule(numlist ascending integer >=1), so most of what is
# checked here is validation the port has to restate: numlist does it before
# slopepower.ado ever runs, and there is no numlist in R. The exception is the
# integer requirement, which the port drops (CONTRACT.md 5.1) -- covered by the
# last test in this section. test-stata-behaviour.R makes the same statements
# about `visits` on trial_design(); these are about the wrapper's own guards,
# which sit on the other side of the interface and are reached by a different
# argument in different units.

test_that("slopepower() rejects a non-positive schedule, as Stata does", {
  d <- load_paper_data("slpower1")
  # Stata declares schedule(numlist ascending integer >=1); baseline is implicit
  # and must not be listed
  expect_error(slopepower(d, depvar = "sdmt", subject = "id", time = "time",
                          schedule = c(0, 1, 2), obs = TRUE, nocontrols = TRUE,
                          effectiveness = 0.33))
})

test_that("slopepower() rejects a schedule that is not a list of visit times", {
  d <- load_paper_data("slpower1")
  sched <- function(s) {
    slopepower(d, depvar = "sdmt", subject = "id", time = "time", schedule = s,
               obs = TRUE, nocontrols = TRUE, effectiveness = 0.33)
  }
  # None of these are numlists at all, so Stata never reaches slopepower.ado
  # with them; the equivalent R inputs have to be turned away here instead. The
  # message is matched, not just the failure, because the positivity and
  # ordering guards are two and three lines further on and would otherwise be
  # indistinguishable from this one.
  msg <- "must be a numeric vector of at least one follow-up visit"
  expect_error(sched("1 2"), msg)          # the Stata syntax, written literally
  expect_error(sched(numeric(0)), msg)     # schedule() is not optional
  expect_error(sched(TRUE), msg)
  expect_error(sched(c(1, NA)), msg)       # numlist has no missing values
  expect_error(sched(c(1, Inf)), msg)
})

test_that("slopepower() rejects a schedule that is not strictly increasing", {
  d <- load_paper_data("slpower1")
  sched <- function(s) {
    slopepower(d, depvar = "sdmt", subject = "id", time = "time", schedule = s,
               obs = TRUE, nocontrols = TRUE, effectiveness = 0.33)
  }
  # `ascending` in numlist means strictly ascending: Stata returns rc 124 for
  # both schedule(3 2 1) and schedule(1 2 2), the SCHED-descending and
  # SCHED-repeat rows of the reference grid. A repeated visit is not merely
  # redundant here -- two rows of Sigma at the same time are collinear, so the
  # design would be singular rather than wrong by a little.
  expect_error(sched(c(2, 1)), "strictly increasing")
  expect_error(sched(c(3, 2, 1)), "strictly increasing")
  expect_error(sched(c(1, 1, 2)), "strictly increasing")
  expect_error(sched(c(1, 2, 2)), "strictly increasing")
})

test_that("slopepower() announces fractional visits, and gets them right", {
  d <- load_paper_data("slpower1")
  # The advertised divergence. Stata returns rc 126 ("noninteger elements") for
  # schedule(1 1.5 2) and can only reach six-monthly visits through scale(); the
  # port evaluates Sigma at the times asked for (CONTRACT.md 5.1), so the p.589
  # design can be written as it is meant. The message exists to tell a reader
  # transcribing a Stata command that they no longer need scale().
  expect_message(
    r <- slopepower(d, depvar = "sdmt", subject = "id", time = "time",
                    schedule = c(0.5, 1, 1.5, 2), obs = TRUE, nocontrols = TRUE,
                    effectiveness = 0.33),
    "non-integer")

  # A notice attached to a wrong number would be worse than no notice, so the
  # arithmetic is checked too: p.589 gives N = 620 for six-monthly visits over
  # two years, and p.588 the year-scale tte of 0.552.
  expect_equal(r$design$visits, c(0, 0.5, 1, 1.5, 2))
  expect_equal(r$n, 620)
  expect_equal(r$n_per_arm, 310)
  expect_lt(abs(r$tte - 0.552), 5e-4)

  # And it is the same trial as the scale(0.5) schedule(1 2 3 4) route tested
  # above, viewed on the other axis: N counts people and does not move, while
  # tte and the slope are per year here and per half-year there. CONTRACT.md's
  # parity table requires both spellings to give 620.
  s <- suppressMessages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time",
    schedule = 1:4, scale = 0.5, obs = TRUE, nocontrols = TRUE,
    effectiveness = 0.33))
  expect_equal(r$n, s$n)
  expect_equal(r$tte, 2 * s$tte)
  expect_equal(r$params$slope, 2 * s$params$slope)

  # The notice is specific to fractional times, and in particular is not the
  # scale() one: a schedule Stata could have written draws neither.
  msgs <- capture_messages(slopepower(
    d, depvar = "sdmt", subject = "id", time = "time", schedule = c(1, 2),
    obs = TRUE, nocontrols = TRUE, effectiveness = 0.33))
  expect_false(any(grepl("non-integer|`scale`", msgs)))
})
