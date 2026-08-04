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
    slope_params_manual(slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
                        sigma_cov = 5, sigma2_residual = 10),
    visits = list(final_only = c(0, 3), annual = 0:3, six_month = seq(0, 3, 0.5)),
    dropout = list(`5pc` = dropout_rate(0.05)), n = 450, effectiveness = 0.33))
  expect_equal(d_final$dropout_total, c(0.15, 0.15, 0.15))
})

test_that("dropout_rate() errors when the implied total exceeds 1", {
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
                           sigma_cov = 5, sigma2_residual = 10)
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
  skip_without_paper_data()
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
  skip_without_paper_data()
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

test_that("slope_power_grid() reproduces all nine Table 1 powers in one call", {
  skip_without_paper_data()
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
  skip_without_paper_data()
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
  skip_without_paper_data()
  out <- suppressWarnings(slope_sample_size_grid(
    paper_fit("slpower1"),
    visits = list(extended = c(0, 1, 2, 5)),
    dropout = list(tenpc_final = c(0, 0, 0.1)),
    power = 0.8, effectiveness = 0.33))
  expect_equal(out$n, 328)
})

test_that("slope_sample_size_grid() agrees with slope_sample_size() cell by cell", {
  skip_without_paper_data()
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
  skip_without_paper_data()
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
  skip_without_paper_data()
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
  skip_without_paper_data()
  err <- expect_error(
    slope_power_grid(paper_fit("slpower1"),
                     visits = list(annual = 0:3),
                     dropout = list(two = c(0.05, 0.05)),
                     n = 450, effectiveness = 0.33))
  expect_match(conditionMessage(err), "dropout_rate")
})

test_that("slope_power_grid() collects the baseline-only warning once", {
  skip_without_paper_data()
  # trial_design() warns per design; the grid should report once, not nine times
  w <- capture_warnings(slope_power_grid(
    paper_fit("slpower1"), visits = table1_visits,
    dropout = list(`10pc` = dropout_rate(0.10)),
    n = 450, effectiveness = 0.33))
  expect_length(w, 1L)
  expect_match(w, "baseline visit only")
})
