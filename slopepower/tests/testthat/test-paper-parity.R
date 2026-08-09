# Reproduction of every published result in Nash et al. (2021),
# Stata Journal 21(3): 575-601. See CONTRACT.md section 7.
#
# These are the tests that establish the port is faithful. Sample sizes are
# integers after ceiling() and are asserted exactly: an off-by-one here is a
# real finding about where the underlying real-valued N sits, not a tolerance to
# be widened.

# Shared fit cache; see the note in test-params.R.
# The paper prints slopes, tte and power to 3 d.p.; compare on that scale.
expect_printed <- function(object, expected, label, digits = 3) {
  tol <- 0.5 * 10^(-digits) + 1e-9
  testthat::expect_lt(abs(object - expected), tol,
                      label = sprintf("%s (got %.6f, paper %.3f)",
                                      label, object, expected))
}

# ---------------------------------------------------------------------------
# Fitted slopes, "Data characteristics" blocks
# ---------------------------------------------------------------------------

test_that("slpower1: the fitted slope matches the paper (p.588)", {
  p <- paper_fit("slpower1")
  # The fit gives -1.6724999999999988, which clears the .0005 rounding boundary
  # against -1.672 by 1.2e-15 -- about one ulp. Comparing at 3 d.p. would pass,
  # but only just: an nlme or BLAS change worth one ulp would flip it, for
  # reasons having nothing to do with the port. Pin at 4 d.p. instead.
  expect_equal(round(p$slope, 4), -1.6725)
})

test_that("slpower2: case and control slopes match the paper (p.590)", {
  p <- paper_fit("slpower2")
  expect_printed(p$slope, -1.715, "slope of cases")
  expect_printed(p$slope_comparator, 0.975, "slope of healthy controls")
  expect_printed(p$slope - p$slope_comparator, -2.690, "observed difference")
})

test_that("slpower3: control and treated arm slopes match the paper (p.594)", {
  p <- paper_fit("slpower3")
  expect_printed(p$slope, -1.852, "slope of control arm")
  expect_printed(p$slope_comparator, -1.104, "slope of experimental arm")
  expect_printed(p$slope - p$slope_comparator, -0.747, "observed difference")
})

# ---------------------------------------------------------------------------
# Section 4.1.1 --- single-group data (slpower1)
# ---------------------------------------------------------------------------

test_that("p.588: annual visits over two years gives N = 712", {
  r <- slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33)
  expect_printed(r$tte, 0.552, "target treatment difference")
  expect_equal(r$n, 712)
  expect_equal(r$n_per_arm, 356)
  expect_equal(r$power, 0.8)
  expect_equal(r$alpha, 0.05)
})

test_that("p.588: extended follow-up with 10% dropout gives N = 328", {
  d <- trial_design(c(0, 1, 2, 5), dropout = c(0, 0, 0.1))
  r <- slope_sample_size(paper_fit("slpower1"), d, effectiveness = 0.33)
  expect_printed(r$tte, 0.552, "target treatment difference")
  expect_equal(r$n, 328)
  expect_equal(r$n_per_arm, 164)
})

test_that("p.589: six-monthly visits over two years gives N = 620", {
  r <- slope_sample_size(paper_fit("slpower1"), c(0, 0.5, 1, 1.5, 2),
                         effectiveness = 0.33)
  expect_equal(r$n, 620)
  expect_equal(r$n_per_arm, 310)
})

test_that("p.589: the same design on a rescaled time axis also gives N = 620", {
  # Stata needs scale(0.5) with schedule(1 2 3 4) because it builds the
  # covariance on a unit-integer grid. This port builds it at the requested
  # times, so the year axis with half-year visits and the half-year axis with
  # integer visits are the same calculation. Agreement is what makes the
  # scale() option unnecessary rather than merely absent.
  d <- load_paper_data("slpower1")
  d$time <- d$time / 0.5
  p <- slope_params(sdmt ~ time | id, d)
  expect_printed(p$slope, -0.836, "slope on the half-year axis")

  r <- slope_sample_size(p, 0:4, effectiveness = 0.33)
  expect_printed(r$tte, 0.276, "tte on the half-year axis")
  expect_equal(r$n, 620)
  expect_equal(r$n_per_arm, 310)

  # identical to the year-axis result
  r_year <- slope_sample_size(paper_fit("slpower1"), c(0, 0.5, 1, 1.5, 2),
                              effectiveness = 0.33)
  expect_equal(r$n, r_year$n)
  expect_equal(r$effect_size, r_year$effect_size, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Section 4.1.2 --- cases and healthy controls (slpower2)
# ---------------------------------------------------------------------------

test_that("p.590: cases versus healthy controls gives N = 296", {
  p <- paper_fit("slpower2")
  r <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  expect_printed(r$tte, 0.888, "target treatment difference")
  expect_equal(r$n, 296)
  expect_equal(r$n_per_arm, 148)
  # the effect is measured toward the healthy-control slope, not toward zero
  expect_equal(r$reference_slope, p$slope_comparator)
})

test_that("p.593: N = 200 with 5% dropout per visit gives power 0.597", {
  # dropout[1] = 0.05 attends baseline only and contributes nothing; the
  # warning is expected here and is asserted rather than suppressed.
  d <- expect_warning(trial_design(c(0, 1, 2), dropout = c(0.05, 0.05)),
                      "contribute nothing")
  r <- slope_power(paper_fit("slpower2"), d, n = 200, effectiveness = 0.33)
  expect_printed(r$tte, 0.888, "target treatment difference")
  expect_printed(r$power, 0.597, "estimated power")
  expect_equal(r$n, 200)
  expect_equal(r$n_per_arm, 100)
})

# ---------------------------------------------------------------------------
# Section 4.1.3 --- a previous trial (slpower3)
# ---------------------------------------------------------------------------

test_that("p.594: targeting the previously observed effect gives N = 318", {
  # dropout[1] = 0.2 means 20% attend baseline only and contribute nothing.
  # trial_design() warns about that; Stata skips the stratum silently.
  d <- expect_warning(trial_design(c(0, 2, 3), dropout = c(0.2, 0.1)),
                      "contribute nothing")
  r <- slope_sample_size(paper_fit("slpower3"), d, target = "observed")
  expect_printed(r$tte, 0.747, "target treatment difference")
  expect_equal(r$n, 318)
  expect_equal(r$n_per_arm, 159)
  expect_true(is.na(r$effectiveness))
})

test_that("p.594: powering for a fraction p of the observed effect scales as p^-2", {
  # Paper section 4.1.3: "multiply the sample size above by 4" for p = 0.5,
  # "so we would need a sample size of 1,272".
  d <- suppressWarnings(trial_design(c(0, 2, 3), dropout = c(0.2, 0.1)))
  p3 <- paper_fit("slpower3")
  ss <- slope_sample_size(p3, d, target = "observed")
  expect_equal(ss$n, 318)
  expect_equal(ss$n * 4, 1272)

  # The published shortcut is slightly conservative, and the two figures it
  # sits between are what this pins. The scaling law is exact on the unrounded
  # per-arm requirement; `n` is already ceiling()-rounded (see solve_slope()),
  # so multiplying magnifies that rounding. Round last instead and the answer
  # is 1,266 -- six people fewer, and the number the rest of the package is
  # consistent with. Neither is wrong; a reader comparing against the paper
  # needs to know both exist, which an assertion that 318 * 4 is 1272 -- true
  # of arithmetic, not of this package -- never told them.
  exact <- 2 * ceiling((stats::qnorm(1 - ss$alpha / 2) + stats::qnorm(ss$power))^2 /
                       (0.5 * ss$effect_size)^2)
  expect_equal(exact, 1266)

  # And that it really is the same calculation, not a coincidence of rounding:
  # halving the previous trial's observed effect at source and re-solving must
  # land on the same 1,266.
  halved <- slope_params_manual(
    slope            = p3$slope,
    slope_comparator = p3$slope + (p3$slope_comparator - p3$slope) / 2,
    sigma2_intercept = p3$sigma2_intercept, sigma2_slope = p3$sigma2_slope,
    sigma_cov        = p3$sigma_cov,        sigma2_residual = p3$sigma2_residual,
    comparator       = "treated")
  expect_equal(slope_sample_size(halved, d, target = "observed")$n, 1266)

  # The gap is bounded by p^-2 participants per arm, which is what makes the
  # shortcut safe to use: it over-recruits, never under-recruits.
  expect_gte(ss$n * 4, exact)
  expect_lt((ss$n * 4 - exact) / 2, 0.5^-2)
})

# ---------------------------------------------------------------------------
# Table 1 (p.595) --- design exploration
# ---------------------------------------------------------------------------

test_that("Table 1: all nine published powers reproduce", {
  p <- paper_fit("slpower1")

  for (i in seq_len(nrow(paper_table1))) {
    design_name  <- paper_table1$design[i]
    dropout_name <- paper_table1$dropout[i]
    expected     <- paper_table1$power[i]

    visits  <- table1_visits[[design_name]]
    dropout <- table1_dropout[[dropout_name]][[design_name]]

    d <- suppressWarnings(trial_design(visits, dropout = dropout))
    r <- slope_power(p, d, n = 450, effectiveness = 0.33)

    expect_printed(r$power, expected,
                   sprintf("Table 1 [%s / %s]", design_name, dropout_name))
  }
})

test_that("Table 1: Stata's own output rounds to the published percentages", {
  # The nine cells above are pinned only to the paper's printed figures, which
  # are percentages to one decimal place: a +/- 5e-4 window on power, and the
  # loosest assertion anywhere in the parity suite. Grid 6 runs the appendix's
  # nine commands (p.601) through slopepower.ado and records the powers at 17
  # significant digits, which closes the gap two ways at once.
  #
  # This test is the first half: Stata's own answer must round to what the
  # journal printed. That is a check on the transcription -- the appendix
  # converts "5% per year" into dropouts(0.15) for the final-visit-only design
  # and dropouts(0.05 0.05 0.05) for annual visits, and if this port had read
  # that conversion wrongly the printed table is the only thing that would ever
  # have said so. The second half is in test-stata-reference.R, where grid 6
  # joins onto grid 5's fits and the powers are pinned to 1e-12.
  skip_without_stata_reference("06_table1")

  ref <- read_stata_reference("06_table1")
  cells <- ref[startsWith(ref$tag, "T1-"), , drop = FALSE]
  expect_equal(nrow(cells), 9)
  expect_true(all(cells$rc == 0))

  tag_of <- c(final_only = "final", annual = "annual", six_month = "sixmonth")

  for (i in seq_len(nrow(paper_table1))) {
    tag <- sprintf("T1-%s-%s", tag_of[[paper_table1$design[i]]],
                   paper_table1$dropout[i])
    row <- cells[cells$tag == tag, , drop = FALSE]
    expect_equal(nrow(row), 1, info = tag)
    expect_printed(row$o_power, paper_table1$power[i],
                   sprintf("Stata, Table 1 [%s]", tag))
  }
})

test_that("Table 1: extra visits buy more power as dropout worsens", {
  # The paper's substantive conclusion: interim visits are worth
  # disproportionately more when dropout is high, because a mixed model uses
  # the data collected before withdrawal.
  p <- paper_fit("slpower1")

  power_of <- function(design_name, dropout_name) {
    d <- suppressWarnings(trial_design(
      table1_visits[[design_name]],
      dropout = table1_dropout[[dropout_name]][[design_name]]))
    slope_power(p, d, n = 450, effectiveness = 0.33)$power
  }

  gain_none <- power_of("six_month", "none")  - power_of("final_only", "none")
  gain_10pc <- power_of("six_month", "10pc")  - power_of("final_only", "10pc")
  expect_gt(gain_10pc, gain_none)

  # and within each dropout level, more visits is never worse
  for (dn in c("none", "5pc", "10pc")) {
    expect_gt(power_of("annual", dn),    power_of("final_only", dn))
    expect_gt(power_of("six_month", dn), power_of("annual", dn))
  }
})

# ---------------------------------------------------------------------------
# Cross-checks
# ---------------------------------------------------------------------------

test_that("every published sample size round-trips to its stated power", {
  cases <- list(
    list(p = "slpower1", visits = c(0, 1, 2), dropout = NULL, n = 712),
    list(p = "slpower1", visits = c(0, 1, 2, 5), dropout = c(0, 0, 0.1), n = 328),
    list(p = "slpower1", visits = c(0, 0.5, 1, 1.5, 2), dropout = NULL, n = 620),
    list(p = "slpower2", visits = c(0, 1, 2), dropout = NULL, n = 296)
  )
  for (cs in cases) {
    d <- suppressWarnings(trial_design(cs$visits, dropout = cs$dropout))
    r <- slope_power(paper_fit(cs$p), d, n = cs$n, effectiveness = 0.33)
    expect_gte(r$power, 0.8)
    expect_lt(r$power, 0.806)   # the ceiling() overshoot only
  }
})

test_that("the published examples tabulate with identical columns", {
  rows <- list(
    as.data.frame(slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33)),
    as.data.frame(slope_sample_size(paper_fit("slpower2"), c(0, 1, 2), effectiveness = 0.33)),
    as.data.frame(slope_sample_size(paper_fit("slpower3"),
                                    suppressWarnings(trial_design(c(0, 2, 3), c(0.2, 0.1))),
                                    target = "observed"))
  )
  reference <- names(rows[[1]])
  for (r in rows) expect_identical(names(r), reference)
  combined <- do.call(rbind, rows)
  expect_equal(nrow(combined), 3L)
  expect_equal(combined$n, c(712, 296, 318))
})
