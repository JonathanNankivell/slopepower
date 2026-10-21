# The naming rule for outputs (CONTRACT.md section 4), checked across every
# object the package returns:
#
#   * a count on a planned trial always comes as a `_total`/`_per_arm` pair,
#     the total twice the per-arm figure, and never bare -- no `n`, no
#     `n_requested`, no `visits`;
#   * a bootstrapped quantity `q` has the fields `q`, `q_replicates`,
#     `q_mean`, `q_sd`, `q_lower` and `q_upper` -- no bare `ci`, `observed`,
#     `boot_mean`, `boot_sd` or `se`;
#   * a scalar result's fields carry the same names as its grid's columns.

# Ten subjects: quick to fit and to refit, so the bootstraps below are cheap.
names_fit <- function() {
  subj <- local({
    set.seed(4)
    data.frame(id = 1:10, a = rnorm(10, 50, 8), b = rnorm(10, -2, 0.6))
  })
  d <- merge(subj, data.frame(visit = 0:3))
  d$y <- d$a + d$b * d$visit + rnorm(nrow(d), 0, 2)
  slope_params(y ~ visit | id, d)
}

every_output <- function() {
  p <- names_fit()
  v <- c(0, 1, 2, 3)
  boot <- function(expr) suppressWarnings(expr)
  list(
    sample_size = slope_sample_size(p, v, effectiveness = 0.33),
    power_total = slope_power(p, v, n_total = 61, effectiveness = 0.33),
    power_arm = slope_power(p, v, n_per_arm = 30, effectiveness = 0.33),
    floor = slope_sample_size_floor(p, effectiveness = 0.33),
    ceiling = slope_power_ceiling(p, n_total = 61, effectiveness = 0.33),
    df_sample_size = as.data.frame(slope_sample_size(p, v, effectiveness = 0.33)),
    df_power = as.data.frame(slope_power(p, v, n_total = 61, effectiveness = 0.33)),
    ss_grid = slope_sample_size_grid(p, list(a = v, b = c(0, 1, 2)), effectiveness = 0.33),
    power_grid = slope_power_grid(p, list(a = v, b = c(0, 1, 2)), n_per_arm = c(20, 30),
                                  effectiveness = 0.33),
    ss_boot = boot(slope_sample_size_boot(p, v, effectiveness = 0.33, R = 6, seed = 1)),
    tte_boot = boot(slope_sample_size_boot(p, v, effectiveness = 0.33, statistic = "tte",
                                           R = 6, seed = 1)),
    power_boot = boot(slope_power_boot(p, v, n_total = 61, effectiveness = 0.33, R = 6,
                                       seed = 1)),
    floor_boot = boot(slope_sample_size_floor_boot(p, effectiveness = 0.33, R = 6, seed = 1)),
    ceiling_boot = boot(slope_power_ceiling_boot(p, n_total = 61, effectiveness = 0.33,
                                                 R = 6, seed = 1)),
    params_boot = boot(slope_params_boot(p, R = 6, seed = 1)),
    ss_grid_boot = boot(slope_sample_size_grid_boot(p, list(a = v, b = c(0, 1, 2)),
                                                    effectiveness = 0.33, R = 6, seed = 1)),
    power_grid_boot = boot(slope_power_grid_boot(p, list(a = v, b = c(0, 1, 2)),
                                                 n_total = 61, effectiveness = 0.33,
                                                 R = 6, seed = 1))
  )
}

outputs <- NULL
get_outputs <- function() {
  if (is.null(outputs)) outputs <<- every_output()
  outputs
}

# Field or column names, and attribute names, of an output.
output_names <- function(x) {
  c(names(x), setdiff(names(attributes(x)), c("names", "row.names", "class")))
}

test_that("no output carries an ambiguous name", {
  banned <- "^(n|n_requested|visits|ci|observed|boot_mean|boot_sd|se|n_failed|n_refit_failed|dropout_total|ci_method|slope_ci|slope_observed)$"
  for (nm in names(get_outputs())) {
    x <- get_outputs()[[nm]]
    found <- grep(banned, output_names(x), value = TRUE)
    expect_identical(found, character(0), info = nm)
  }
})

test_that("every _total has a _per_arm twin of half its value", {
  for (nm in names(get_outputs())) {
    x <- get_outputs()[[nm]]
    fields <- names(x)
    totals <- grep("_total(_|$)", fields, value = TRUE)
    for (t in totals) {
      twin <- sub("_total(_|$)", "_per_arm\\1", t)
      expect_true(twin %in% fields, info = paste(nm, t))
      # A standard deviation halves like the rest: it scales with its count.
      expect_equal(x[[t]], 2 * x[[twin]], info = paste(nm, t))
    }
  }
})

test_that("a bootstrapped quantity always has the same six fields", {
  suffixes <- c("", "_replicates", "_mean", "_sd", "_lower", "_upper")
  quantities <- list(ss_boot = c("n_total", "n_per_arm", "slope"),
                     tte_boot = c("tte", "slope"),
                     power_boot = c("power", "slope"),
                     floor_boot = c("n_total", "n_per_arm", "slope"),
                     ceiling_boot = c("power", "slope"),
                     params_boot = "slope")
  for (nm in names(quantities)) {
    x <- get_outputs()[[nm]]
    for (q in quantities[[nm]]) {
      expect_true(all(paste0(q, suffixes) %in% names(x)), info = paste(nm, q))
    }
    expect_true(paste0(x$statistic, "_ci_method") %in% names(x), info = nm)
    expect_true("slope_ci_method" %in% names(x), info = nm)
  }
  # And a grid has the same summary columns, per cell, for its statistic.
  gb <- get_outputs()$ss_grid_boot
  for (q in c("n_total", "n_per_arm", "tte")) {
    expect_true(all(paste0(q, c("_mean", "_sd", "_lower", "_upper")) %in% names(gb)), info = q)
  }
})

test_that("a scalar result's fields are named as its grid's columns", {
  o <- get_outputs()
  grid_cols <- names(o$ss_grid)
  # `design` is a label in a grid and the design itself in a result.
  shared <- setdiff(intersect(grid_cols, names(o$sample_size)), "design")
  expect_true(all(c("n_total", "n_per_arm", "visits_total", "visits_per_arm", "power",
                    "tte", "var_tte", "effect_size") %in% shared))
  expect_equal(unlist(strip_grid(o$ss_grid)[1L, shared]), unlist(o$sample_size[shared]))

  # as.data.frame() carries every result column of a grid row too.
  result_cols <- setdiff(grid_cols, c("design", "dropout", "scheduled_visits", "last_visit",
                                      "dropout_cumulative"))
  expect_true(all(result_cols %in% names(o$df_sample_size)))
  expect_true(all(result_cols %in% names(o$df_power)))

  # A bootstrap's statistic fields are its grid's columns for the same cell.
  gb <- o$ss_grid_boot
  b <- o$ss_boot
  for (nm in c("n_total_mean", "n_per_arm_lower", "n_per_arm_upper")) {
    expect_equal(gb[[nm]][1L], b[[nm]], info = nm)
  }
})

test_that("confint() reads the stored _lower and _upper fields", {
  o <- get_outputs()
  for (nm in c("ss_boot", "tte_boot", "power_boot", "floor_boot", "ceiling_boot",
               "params_boot")) {
    b <- o[[nm]]
    ci <- confint(b)
    for (q in rownames(ci)) {
      expect_equal(unname(ci[q, ]), c(b[[paste0(q, "_lower")]], b[[paste0(q, "_upper")]]),
                   info = paste(nm, q))
    }
  }
})

test_that("a sample size given per arm and in total land in the same fields", {
  o <- get_outputs()
  expect_equal(o$power_arm$n_total, 60)
  expect_equal(o$power_arm$n_per_arm, 30)
  expect_equal(o$power_arm$n_requested_total, 60)
  expect_equal(o$power_arm$n_requested_per_arm, 30)
  # 61 is evened down to 60; the request is kept as asked.
  expect_equal(o$power_total$n_total, 60)
  expect_equal(o$power_total$n_requested_total, 61)
  expect_equal(o$power_total$n_requested_per_arm, 30.5)
  expect_equal(o$power_total$power, o$power_arm$power)
  # A per-arm axis of a grid is labelled as given.
  expect_equal(o$power_grid$n_requested_per_arm, rep(c(20, 30), 2L))
})

test_that("every print opens with a title line and a blank line", {
  titles <- c(sample_size = "Sample size", power_total = "Power",
              floor = "Sample size: lower bound over all visit schedules",
              ceiling = "Power: upper bound over all visit schedules",
              ss_grid = "Sample size grid (2 rows)", power_grid = "Power grid (4 rows)",
              ss_boot = "Bootstrapped sample size", power_boot = "Bootstrapped power",
              tte_boot = "Bootstrapped target treatment effect",
              params_boot = "Bootstrapped slope",
              ss_grid_boot = "Bootstrapped sample size grid (2 rows)",
              power_grid_boot = "Bootstrapped power grid (2 rows)")
  for (nm in names(titles)) {
    out <- capture.output(print(get_outputs()[[nm]]))
    expect_identical(out[1:2], c(titles[[nm]], ""), info = nm)
  }
  out <- capture.output(print(get_outputs()$sample_size$params))
  expect_identical(out[1L], "Slope parameters: single group, target no change over time")
  expect_match(capture.output(print(dropout_rate(0.05))), "^Dropout rate: ")
  expect_identical(capture.output(print(get_outputs()$sample_size$design))[1L], "Trial design")
})
