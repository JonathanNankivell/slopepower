# Layer 4 --- slope_power_grid(), slope_sample_size_grid() and dropout_rate().
#
# The point of the grid is section 4.2 of the paper: fit stage one once, then
# price out many trial designs. In Stata each row of Table 1 refits the mixed
# model; here the fit is reused.

# --- dropout_rate() ---------------------------------------------------------

test_that("dropout_rate() expands to per-interval proportions", {
  # A rate of 5% per year gives 0.15 across a single three-year interval, and
  # 0.05 per year when there are annual visits: exactly the conversion the
  # paper's appendix code does by hand for each row of Table 1.
  expect_s3_class(dropout_rate(0.05), "dropout_rate")
  d_final  <- suppressWarnings(slope_power_grid(
    ref_params(),
    visits = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
    dropout = list(`5pc` = dropout_rate(0.05)), n = 450, effectiveness = 0.33))
  expect_equal(d_final$dropout_total, c(0.15, 0.15, 0.15))
})

test_that("dropout_rate() errors when the implied total exceeds 1", {
  p <- ref_params()
  expect_error(
    slope_power_grid(p, visits = list(long = c(0, 10)),
                     dropout = list(fast = dropout_rate(0.5)),
                     n = 450, effectiveness = 0.33),
    "exceeds 1"
  )
})

test_that("print.dropout_rate() runs", {
  expect_output(print(dropout_rate(0.05)), "dropout_rate")
})

# --- slope_power_grid() -----------------------------------------------------

test_that("slope_power_grid() returns one row per design/dropout combination", {
  out <- suppressWarnings(slope_power_grid(
    paper_fit("slpower1"),
    visits  = table1_visits,
    dropout = list(none = NULL, `5pc` = dropout_rate(0.05)),
    n = 450, effectiveness = 0.33
  ))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), length(table1_visits) * 2L)
  expect_true(all(c("design", "dropout", "n", "n_per_arm", "power",
                    "tte", "var_tte", "effect_size") %in% names(out)))
  expect_setequal(unique(out$design), names(table1_visits))
  expect_setequal(unique(out$dropout), c("none", "5pc"))
})

test_that("slope_power_grid() agrees with slope_power() cell by cell", {
  p <- paper_fit("slpower1")
  out <- suppressWarnings(slope_power_grid(
    p, visits = table1_visits, dropout = list(none = NULL),
    n = 450, effectiveness = 0.33))

  for (nm in names(table1_visits)) {
    direct <- slope_power(p, trial_design(table1_visits[[nm]]),
                          n = 450, effectiveness = 0.33)$power
    expect_equal(out$power[out$design == nm], direct, tolerance = 1e-12)
  }
})

test_that("a grid row carries exactly what as.data.frame() reports for that cell", {
  # grid_impl() reads its six result columns straight off the result object
  # rather than through as.data.frame.slope_result(), because routing every
  # cell through that method to keep 6 of its 18 columns was the single most
  # expensive line in the loop. This is the guarantee that bought back: the two
  # must report the same numbers for the same object, so that a grid row and a
  # bound as.data.frame() row can never disagree.
  shared <- c("n", "n_per_arm", "power", "tte", "var_tte", "effect_size")
  p <- paper_fit("slpower1")

  for (nm in names(table1_visits)) {
    v <- table1_visits[[nm]]
    for (drop in list(NULL, rep(0.05, length(v) - 1L))) {
      des <- suppressWarnings(trial_design(v, drop))
      out <- suppressWarnings(slope_power_grid(
        p, visits = list(only = v), dropout = list(only = drop),
        n = 450, effectiveness = 0.33))
      direct <- as.data.frame(suppressWarnings(
        slope_power(p, des, n = 450, effectiveness = 0.33)))
      expect_equal(unlist(out[1L, shared]), unlist(direct[shared]),
                   tolerance = 0, info = nm)
    }
  }

  # The same for the sample-size grid, whose rows come off a different class.
  out <- slope_sample_size_grid(p, visits = list(only = c(0, 1, 2)),
                                dropout = list(only = NULL), effectiveness = 0.33)
  direct <- as.data.frame(slope_sample_size(p, trial_design(c(0, 1, 2)),
                                            effectiveness = 0.33))
  expect_equal(unlist(out[1L, shared]), unlist(direct[shared]), tolerance = 0)
})

test_that("slope_power_grid() reproduces all nine Table 1 powers in one call", {
  out <- suppressWarnings(slope_power_grid(
    paper_fit("slpower1"),
    visits  = table1_visits,
    dropout = list(none = NULL,
                   `5pc`  = dropout_rate(0.05),
                   `10pc` = dropout_rate(0.10)),
    n = 450, effectiveness = 0.33
  ))
  expect_equal(nrow(out), 9L)

  for (i in seq_len(nrow(paper_table1))) {
    got <- out$power[out$design == paper_table1$design[i] &
                       out$dropout == paper_table1$dropout[i]]
    expect_length(got, 1L)
    expect_lt(abs(got - paper_table1$power[i]), 5e-4)
  }
})

test_that("slope_sample_size_grid() reports the N each design needs", {
  out <- suppressWarnings(slope_sample_size_grid(
    paper_fit("slpower1"),
    visits = list(annual = c(0, 1, 2)),
    dropout = list(none = NULL),
    power = 0.8, effectiveness = 0.33
  ))
  expect_equal(nrow(out), 1L)
  expect_equal(out$n, 712)
  expect_equal(out$n_per_arm, 356)
  # power is the held-constant input on this side of the split
  expect_equal(out$power, 0.8)
})

test_that("slope_sample_size_grid() accepts explicit dropout vectors", {
  out <- suppressWarnings(slope_sample_size_grid(
    paper_fit("slpower1"),
    visits = list(extended = c(0, 1, 2, 5)),
    dropout = list(tenpc_final = c(0, 0, 0.1)),
    power = 0.8, effectiveness = 0.33))
  expect_equal(out$n, 328)
})

test_that("slope_sample_size_grid() agrees with slope_sample_size() cell by cell", {
  p <- paper_fit("slpower1")
  out <- suppressWarnings(slope_sample_size_grid(
    p, visits = table1_visits, dropout = list(none = NULL),
    power = 0.9, effectiveness = 0.33))

  for (nm in names(table1_visits)) {
    direct <- slope_sample_size(p, trial_design(table1_visits[[nm]]),
                                power = 0.9, effectiveness = 0.33)$n
    expect_equal(out$n[out$design == nm], direct)
  }
})

test_that("the two grids are inverses, design by design", {
  # Same shared loop, same dropout expansion: feeding one grid's N back into the
  # other must recover the target power in every cell.
  p <- paper_fit("slpower1")
  sizes <- suppressWarnings(slope_sample_size_grid(
    p, visits = table1_visits, dropout = list(`5pc` = dropout_rate(0.05)),
    power = 0.8, effectiveness = 0.33))
  expect_identical(names(sizes),
                   names(suppressWarnings(slope_power_grid(
                     p, visits = table1_visits,
                     dropout = list(`5pc` = dropout_rate(0.05)),
                     n = 450, effectiveness = 0.33))))
  for (i in seq_len(nrow(sizes))) {
    des <- suppressWarnings(trial_design(
      table1_visits[[sizes$design[i]]],
      dropout = dropout_rate(0.05)$rate * diff(table1_visits[[sizes$design[i]]])))
    back <- slope_power(p, des, n = sizes$n[i], effectiveness = 0.33)
    expect_gte(back$power, 0.8)
    expect_lt(back$power - 0.8, 0.005)
  }
})

test_that("each grid requires the input its question needs", {
  p <- paper_fit("slpower1")
  # slope_power_grid() holds n fixed; without it there is nothing to evaluate
  expect_error(slope_power_grid(p, visits = list(a = c(0, 1, 2))), "`n` is required")
  expect_match(tryCatch(slope_power_grid(p, visits = list(a = c(0, 1, 2))),
                        error = conditionMessage),
               "slope_sample_size_grid")
  # and neither accepts the other's input
  expect_error(slope_power_grid(p, visits = list(a = c(0, 1, 2)), n = 450, power = 0.8),
               "unused argument")
  expect_error(slope_sample_size_grid(p, visits = list(a = c(0, 1, 2)), n = 450),
               "unused argument")
})

test_that("a dropout vector of the wrong length for a design errors helpfully", {
  err <- expect_error(
    slope_power_grid(paper_fit("slpower1"),
                     visits = list(annual = 0:3),
                     dropout = list(two = c(0.05, 0.05)),
                     n = 450, effectiveness = 0.33))
  expect_match(conditionMessage(err), "dropout_rate")
})

test_that("slope_power_grid() collects the baseline-only warning once", {
  # trial_design() warns per design; the grid should report once, not nine times
  w <- capture_warnings(slope_power_grid(
    paper_fit("slpower1"), visits = table1_visits,
    dropout = list(`10pc` = dropout_rate(0.10)),
    n = 450, effectiveness = 0.33))
  expect_length(w, 1L)
  expect_match(w, "baseline visit only")
})

# --- automatic row labels ---------------------------------------------------
#
# `design` and `dropout` are the row identifiers of the returned table: they are
# the columns a caller subsets on to pull out one result, as the Table 1 tests
# above do. What the grid puts in them when the inputs are *not* named is
# therefore part of the interface rather than cosmetic, and every test above
# supplies named lists. These supply the unnamed forms and assert the labels
# together with the numbers on the rows they identify.

test_that("a bare numeric `visits` vector is labelled by its visit times", {
  # The one-design shorthand -- `visits = c(0, 1, 2)` instead of a list of one --
  # and `dropout` left at its default, which stands for a single no-dropout row.
  p <- paper_fit("slpower1")
  out <- slope_sample_size_grid(p, visits = c(0, 1, 2),
                                power = 0.8, effectiveness = 0.33)
  expect_equal(nrow(out), 1L)
  expect_identical(out$design, "0, 1, 2")
  expect_identical(out$dropout, "none")
  expect_equal(out$dropout_total, 0)
  # CONTRACT.md 7, paper p.588: the shorthand must reach the same calculation
  # the list form does, so the headline N of the first worked example stands.
  expect_equal(out$n, 712)
  expect_equal(out$n_per_arm, 356)
})

test_that("a bare dropout vector or dropout_rate() is labelled by its contents", {
  p <- paper_fit("slpower1")
  vec <- suppressWarnings(slope_power_grid(
    p, visits = list(annual = 0:3), dropout = c(0.05, 0.05, 0.05),
    n = 450, effectiveness = 0.33))
  expect_identical(vec$dropout, "0.05, 0.05, 0.05")
  expect_equal(vec$dropout_total, 0.15)

  # A dropout_rate() is labelled by the rate it states, not by the vector it
  # expands to. That is the only stable choice: the same object deliberately
  # yields a different vector for every schedule (0.15 across one interval,
  # rep(0.05, 3) across three), so a vector label would differ from row to row
  # of a single column.
  rate <- suppressWarnings(slope_power_grid(
    p, visits = table1_visits, dropout = dropout_rate(0.05),
    n = 450, effectiveness = 0.33))
  expect_identical(unique(rate$dropout), "0.05 per 1")
  expect_equal(rate$dropout_total, rep(0.15, 3))
  # and the rows so labelled are the `5pc` column of Table 1, p.595
  for (i in seq_len(nrow(rate))) {
    expect_lt(abs(rate$power[i] - paper_table1$power[paper_table1$dropout == "5pc" &
                                                       paper_table1$design == rate$design[i]]),
              5e-4)
  }

  # `per` is carried into the label as well as into the expansion, because 10%
  # per twelve units of time and 10% per unit time are different assumptions
  # and would otherwise key the same row.
  per <- suppressWarnings(slope_power_grid(
    p, visits = list(annual = 0:3), dropout = dropout_rate(0.1, per = 12),
    n = 450, effectiveness = 0.33))
  expect_identical(per$dropout, "0.1 per 12")
  expect_equal(per$dropout_total, 0.025)
})

test_that("a wholly unnamed list is labelled element by element and deduped", {
  # names() is NULL rather than a vector of blanks when nothing at all is named,
  # which is a distinct case from the partly-named list below.
  p <- paper_fit("slpower1")
  out <- slope_power_grid(p, visits = list(c(0, 3), 0:3, c(0, 3)),
                          n = 450, effectiveness = 0.33)
  # The repeated schedule generates a repeated label, and make.unique() suffixes
  # it rather than emitting two rows the caller cannot tell apart.
  expect_identical(out$design, c("0, 3", "0, 1, 2, 3", "0, 3_1"))
  expect_equal(out$n_visits, c(2L, 4L, 2L))
  # the suffixed row is genuinely the same design, evaluated a second time
  expect_equal(out$power[3], out$power[1])
  # first two rows of the `none` column of Table 1, p.595
  expect_lt(abs(out$power[1] - 0.798), 5e-4)
  expect_lt(abs(out$power[2] - 0.817), 5e-4)
})

test_that("labels fill in only the blanks of a partly named list", {
  p <- paper_fit("slpower1")
  out <- suppressWarnings(slope_power_grid(
    p,
    visits  = list(final_only = c(0, 3), 0:3),
    dropout = list(NULL, `5pc` = dropout_rate(0.05)),
    n = 450, effectiveness = 0.33))
  expect_identical(out$design,
                   c("final_only", "final_only", "0, 1, 2, 3", "0, 1, 2, 3"))
  # A NULL element inside a list is labelled "none", the same string the
  # `dropout = NULL` shorthand produces, so the two routes to a no-dropout row
  # are keyed identically and results from the two calls line up.
  expect_identical(out$dropout, rep(c("none", "5pc"), 2))
  expect_equal(out$dropout_total, c(0, 0.15, 0, 0.15))
  # visits varies slowest, dropout fastest -- the order the Table 1 tests rely on
  expect_equal(out$n_visits, c(2L, 2L, 4L, 4L))
})

test_that("duplicate names supplied by the caller are deduped too", {
  # make.unique() runs over the whole name vector, not just the filled-in
  # blanks, so reusing a name still yields distinguishable rows.
  p <- paper_fit("slpower1")
  out <- slope_power_grid(p, visits = list(sparse = c(0, 3), sparse = 0:3),
                          n = 450, effectiveness = 0.33)
  expect_identical(out$design, c("sparse", "sparse_1"))
  expect_equal(out$n_visits, c(2L, 4L))
})

test_that("visits and dropout must be of a shape the grid can normalise", {
  p <- paper_fit("slpower1")
  # Not a numeric vector and not a non-empty list. An empty list is rejected
  # rather than returning a zero-row table, because a grid with no cells is a
  # mistake in the call rather than a result.
  expect_error(slope_power_grid(p, visits = "0 1 2", n = 450),
               "must be a numeric vector of visit times")
  expect_error(slope_power_grid(p, visits = list(), n = 450),
               "must be a numeric vector of visit times")
  expect_error(slope_sample_size_grid(p, visits = list(), power = 0.8),
               "must be a numeric vector of visit times")
  expect_error(slope_power_grid(p, visits = list(a = c(0, 1, 2)), dropout = "5%", n = 450),
               "must be NULL, a numeric vector, a dropout_rate")
  expect_error(slope_power_grid(p, visits = list(a = c(0, 1, 2)), dropout = list(), n = 450),
               "must be NULL, a numeric vector, a dropout_rate")
})

test_that("a bad list element is reported against the label of its own row", {
  # Normalisation only checks the shape of the container; a bad *element*
  # surfaces later, from inside the loop. With many cells the message has to say
  # which one, so it quotes the row label -- including the fallback label the
  # element gets when it cannot be labelled from its contents.
  p <- paper_fit("slpower1")
  expect_error(slope_power_grid(p, visits = list(bad = "0 1 2"), n = 450),
               'element "bad" of `visits` is not numeric')
  expect_error(slope_power_grid(p, visits = list("0 1 2"), n = 450),
               'element "design" of `visits` is not numeric')

  err <- expect_error(slope_power_grid(p, visits = list(a = c(0, 1, 2)),
                                       dropout = list("5%"), n = 450))
  expect_match(conditionMessage(err), '(dropout = "dropout")', fixed = TRUE)
  expect_match(conditionMessage(err), "got character", fixed = TRUE)
})

test_that("a cell that fails is reported by its labels, with the cause attached", {
  # A stage-two failure is not the grid's own error, so it is re-raised with the
  # two labels naming the offending cell and the original message kept beneath.
  # Here the labels are generated ones, which is the case that matters: the
  # caller has nothing else to identify the row by.
  p <- ref_params(comparator = "healthy", slope = 1, slope_comparator = 1)
  err <- expect_error(slope_power_grid(p, visits = c(0, 1, 2), n = 450,
                                       effectiveness = 0.33))
  expect_match(conditionMessage(err),
               'design "0, 1, 2" with dropout "none" failed', fixed = TRUE)
  # CONTRACT.md 6: a zero slope difference errors rather than yielding N = .
  expect_match(conditionMessage(err), "slope difference is exactly zero", fixed = TRUE)
})
