# Parity against the golden-reference tables in ../stata-reference/.
#
# test-paper-parity.R pins the 15 numbers Nash et al. printed. This file pins
# everything else: the option-space corners the paper had no reason to visit,
# where a transcription slip would never show up against the published values.
#
# Every test here skips unless the corresponding CSV exists, so the suite is
# unchanged on a machine without Stata. See stata-reference/README.md for how
# to produce them.
#
# Tolerances. Slopes and variance components come from two different REML
# implementations and agree to roughly 8 significant figures, not to the last
# bit; tte and var_tte inherit that. Sample sizes are integers after ceiling()
# and are asserted exactly --- an off-by-one means the real-valued N sits within
# rounding distance of an integer, which is worth knowing rather than papering
# over.

# Tolerances, calibrated against the first full run rather than guessed.
#
# The fixed effects agree to 1e-15 on slpower1 and the variance components to
# about 1e-5 relative. That asymmetry is not a bug and not a coincidence:
# slpower1 is balanced, so the GLS slope is the OLS slope whatever the
# covariance structure, and the slope is therefore insensitive to exactly the
# quantity the two REML optimisers disagree about. Tightening nlme's control
# moves var_tte from 8.36388645 to 8.36394022 against Stata's 8.36394172, which
# is what confirms the disagreement is convergence rather than algebra.
#
# So these tolerances measure REML agreement, not correctness of the port. The
# tolerance-free test of the arithmetic is the one in
# test_that("... Stata's own variance components ..."), which removes the fit
# from the comparison entirely; keep that one tight and these ones honest.
# Measured, not guessed. Across all 542 comparable rows the worst disagreement
# in every quantity is the same row --- U3-visit<2, the two-timepoint slpower3
# subset, where a random-slope model has one degree of freedom per subject and
# the likelihood is nearly flat, so the two optimisers stop in visibly
# different places. Everything else on slpower1 agrees in the slope to 1e-8.
#
# A loose maximum on its own would hide a systematic error, so each quantity is
# checked twice: a loose bound on the worst row, and a tight bound on the
# median. Algebra that had drifted would move the median; an ill-conditioned
# fit moves only its own row.
TOL_SLOPE     <- 2e-5   # relative, worst row 1.0e-5
TOL_DERIV     <- 1e-2   # relative, on tte and var_tte; worst row 5.6e-3
TOL_POWER     <- 1e-5   # absolute, worst row 4.1e-6
TOL_N         <- 1e-2   # relative, worst row 5.5e-3
TOL_MEDIAN    <- 1e-7   # relative, on the median row of slope, tte and N
TOL_MEDIAN_VAR <- 5e-5  # relative, on the median row of var_tte
#
# The two median tolerances differ by three orders of magnitude, and the gap is
# the whole story of this comparison. var_tte is a pure function of the four
# variance components, so it inherits the REML disagreement directly and its
# median sits at 6e-6. The slope does not: on balanced data the GLS estimate is
# the OLS estimate whatever the covariance structure, so its median sits at
# 1e-15. Same fits, same rows --- the difference is which quantity depends on
# the part the optimisers disagree about.

# Rows the R port refuses, or answers differently, on purpose. Every entry here
# was checked against the actual Stata output rather than predicted from the
# .ado: an unexplained entry in this vector is worse than an unexplained
# failure.
KNOWN_DIVERGENCES <- c(
  # dropouts(1 0 0): everyone drops out at the first follow-up, so every
  # stratum but the completers is the baseline-only stratum and the completer
  # stratum has weight 0. The effect size is 0 and N is undefined. Stata
  # returns rc = 0 with N = missing; the R port errors, per CONTRACT.md 6.
  "ZERO-allbaseline",
  # Solving for power with dropout, Stata back-solves var_tte through the power
  # it just computed. Once power saturates at exactly 1, invnormal(1) is
  # missing and the reported var_tte goes missing with it --- observed from
  # n = 10000 upward on this design, finite at n = 2000. The R port uses the
  # algebraically equivalent closed form, which stays finite.
  "SAT-drop-10000", "SAT-drop-40000",
  "SAT-drop-200000", "SAT-drop-1000000",
  # Stata warns and then ignores these; the R port has no way to express an
  # ignored argument, so the call shape simply does not exist.
  "IGN-nocontvar-model2", "IGN-usetrt-model2",
  "IGN-casecon-model3", "IGN-nocontvar-model3",
  # Stata falls back to power(0.8) when neither power nor n is given. In the R
  # port that default survives only in the slopepower() compatibility wrapper.
  "NEITHER-default", "DEFAULT-eff"
)

# Rows that test Stata's own syntax validation and have no R counterpart at
# all: the R API cannot express the malformed call, so there is nothing to
# compare. These are not divergences, they are untranslatable.
STATA_ONLY <- c(
  # One bimodal entry point in Stata, two functions in R: `power` and `n` are
  # arguments of different functions and cannot both be passed.
  "BOTH-power-n",
  # Likewise usetrt and effectiveness: `target = "observed"` is not an
  # alternative spelling of an effectiveness value.
  "CONF-usetrt-eff",
  # Missing or wrongly coded group variables. In R you do not name a treatment
  # variable and omit it; you either pass `treated =` or you do not.
  "MISS-treat", "MISS-casecon", "MISS-both-models", "MISS-no-model",
  "MISS-rct-nocont", "GRP-12", "GRP-3levels"
)

# Verified NOT to be divergences on the first full run, recorded so they are
# not "fixed" back into the list above by someone reading the .ado:
#
#   LEN-short, LEN-long    Both return rc = 198. The dropout-length guard at
#                          slopepower.ado:279 is live --- Stata trims the
#                          trailing space in `sched_length ', so the counters
#                          are correct and the guard fires. The R port's strict
#                          length check restates Stata rather than fixing it.
#   SUM1-exact             dropouts(0.3 0.3 0.4) returns rc = 0, N = 1418.
#                          Stata's accumulation does not go negative, so both
#                          implementations accept a dropout list summing to 1.

# ---------------------------------------------------------------------------
# The comparison, shared by every grid
# ---------------------------------------------------------------------------

compare_stata_grid <- function(stem) {
  cmp <- stata_replay_table(stem)
  cmp <- cmp[!cmp$tag %in% STATA_ONLY, , drop = FALSE]
  ok <- cmp[cmp$rc == 0 & is.na(cmp$r_error), , drop = FALSE]
  ok <- ok[!ok$tag %in% KNOWN_DIVERGENCES, , drop = FALSE]

  # 1. Both sides agreed the input was usable.
  disagreed <- cmp[xor(cmp$rc == 0, is.na(cmp$r_error)), , drop = FALSE]
  disagreed <- disagreed[!disagreed$tag %in% KNOWN_DIVERGENCES, , drop = FALSE]
  testthat::expect_equal(
    nrow(disagreed), 0,
    info = stata_mismatch_report(disagreed, c("rc", "r_error"),
                                 paste0(stem, ": accept/refuse"))
  )

  if (nrow(ok) == 0) return(invisible(cmp))

  # 2. Sample size.
  #
  # Not asserted exactly, and the reason matters. N is ceiling()ed, so it is a
  # step function of the effect size: wherever the real-valued N lands within
  # 1e-5 of an integer, the REML gap between the two fits is enough to step it.
  # Two-timepoint designs are the worst case --- the random-slope model is
  # barely identified there and the likelihood is nearly flat, which is why
  # U3-visit<2 moves by 18 in 3256 while everything on slpower1 moves by 2 in
  # 200000. Exactness is tested where it is meaningful: in the manual-params
  # test below, which feeds Stata's own variance components in and so removes
  # the fit from the comparison entirely.
  #
  # The rate is checked as well as the size, because a systematic break would
  # show up as many rows moving a little, not a few rows moving a lot.
  n_ok <- !is.na(ok$o_n) & !is.na(ok$r_n)
  bad <- ok[n_ok & rel_diff(ok$o_n, ok$r_n) > TOL_N, , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_n", "r_n"), paste0(stem, ": N"))
  )
  # How many rows are allowed to miss N by any amount. A rate alone is the
  # wrong instrument on a small grid: grid 4 has 31 comparable rows and two of
  # them are EFF-0.0001 and SCALE-0.001, deliberately extreme rows where N runs
  # to 7.7e9 and 4.5e8, and at that magnitude a 1e-5 relative gap in the
  # variance components moves N by thousands with certainty. So allow a
  # proportion or a small absolute count, whichever is more generous, and let
  # the exact test carry the real burden.
  inexact <- sum(ok$o_n[n_ok] != ok$r_n[n_ok])
  testthat::expect_lte(inexact, max(3, 0.05 * sum(n_ok)))

  # The median guard, applied to every quantity at once. See the note on the
  # tolerances: this is what would catch algebra that had drifted uniformly,
  # which a maximum loose enough to admit an ill-conditioned fit cannot.
  med <- function(a, b) stats::median(rel_diff(a, b), na.rm = TRUE)
  testthat::expect_lt(med(ok$o_slope, ok$r_slope), TOL_MEDIAN)
  testthat::expect_lt(med(ok$o_tte, ok$r_tte), TOL_MEDIAN)
  testthat::expect_lt(med(ok$o_var_tte, ok$r_var_tte), TOL_MEDIAN_VAR)
  testthat::expect_lt(med(ok$o_n, ok$r_n), TOL_MEDIAN)

  # 3. Power.
  bad <- ok[!is.na(ok$o_power) & !is.na(ok$r_power) &
              abs(ok$o_power - ok$r_power) > TOL_POWER, , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_power", "r_power"),
                                 paste0(stem, ": power"))
  )

  # 4. The fitted slope, and the comparator slope where Stata reports one.
  bad <- ok[!is.na(ok$o_slope) &
              rel_diff(ok$o_slope, ok$r_slope) > TOL_SLOPE, , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_slope", "r_slope"),
                                 paste0(stem, ": slope"))
  )

  has_comp <- ok[!is.na(ok$o_slope_comp), , drop = FALSE]
  bad <- has_comp[rel_diff(has_comp$o_slope_comp, has_comp$r_slope_comp) >
                    TOL_SLOPE, , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_slope_comp", "r_slope_comp"),
                                 paste0(stem, ": comparator slope"))
  )

  # 5. tte.
  bad <- ok[!is.na(ok$o_tte) &
              rel_diff(ok$o_tte, ok$r_tte) > TOL_DERIV, , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_tte", "r_tte"), paste0(stem, ": tte"))
  )

  # 6. var_tte, on the rows where the two are defined the same way; the
  #    predicate is CONTRACT.md 5.6 and lives in the helper.
  same_def <- ok[stata_var_tte_comparable(ok), , drop = FALSE]
  bad <- same_def[!is.na(same_def$o_var_tte) &
                    rel_diff(same_def$o_var_tte, same_def$r_var_tte) > TOL_DERIV,
                  , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_var_tte", "r_var_tte"),
                                 paste0(stem, ": var_tte"))
  )

  # 7. The fit used the same data. Guarded on both columns, not just o_obs: a
  #    missing o_subjects would otherwise leave the mask NA, select an all-NA
  #    row, and print the failure as NA rather than as the mismatch report.
  bad <- ok[!is.na(ok$o_obs) & !is.na(ok$o_subjects) &
              (ok$o_obs != ok$r_obs |
                 ok$o_subjects != ok$r_subjects), , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_obs", "r_obs", "o_subjects",
                                        "r_subjects"),
                                 paste0(stem, ": rows in model"))
  )

  # 8. The effectiveness both sides report. On these rows it is an echo of the
  #    input, or of the default when `effin` is blank, so this is not a test of
  #    the arithmetic --- it is a test that the two defaults agree and that the
  #    value survives the round trip. Asserted exactly, with no tolerance: both
  #    sides parse the same decimal string to the same double, and all 525
  #    comparable rows agree bit for bit. A tolerance here would only hide the
  #    one thing this can catch.
  #
  #    usetrt rows are excluded because the port deliberately reports NA there
  #    rather than Stata's back-computed value; that divergence is pinned in
  #    test-stata-behaviour.R instead.
  eff <- ok[ok$usetrt != 1 & !is.na(ok$o_effectiveness), , drop = FALSE]
  bad <- eff[eff$o_effectiveness != eff$r_effectiveness, , drop = FALSE]
  testthat::expect_equal(
    nrow(bad), 0,
    info = stata_mismatch_report(bad, c("o_effectiveness", "r_effectiveness"),
                                 paste0(stem, ": effectiveness"))
  )

  invisible(cmp)
}

# ---------------------------------------------------------------------------

test_that("the committed fixtures are the ones the manifest describes", {
  skip_without_stata_reference("01_grid_arithmetic")
  mf <- file.path(stata_reference_dir(), "MANIFEST")
  skip_if_not(file.exists(mf), "reading regenerated CSVs, not the fixtures")

  lines <- readLines(mf)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  m <- utils::read.csv(text = paste(lines, collapse = "\n"), header = FALSE,
                       col.names = c("stem", "rows", "cols"),
                       stringsAsFactors = FALSE)

  # A fixture regenerated without rerunning make-fixtures.R, or a .do file
  # edited without regenerating, both show up here rather than as a confusing
  # partial mismatch three tests later.
  for (i in seq_len(nrow(m))) {
    x <- read_stata_reference(m$stem[i])
    expect_equal(nrow(x), m$rows[i], info = m$stem[i])
    expect_equal(ncol(x), m$cols[i], info = m$stem[i])
  }
})

# One test_that per grid, so a failure names the grid rather than one of four
# indistinguishable iterations --- but driven off DESIGN_GRIDS so that the set
# of grids under test is stated once, here and in the exact test below.
for (stem in names(DESIGN_GRIDS)) {
  local({
    stem <- stem
    test_that(DESIGN_GRIDS[[stem]], {
      skip_without_stata_reference(stem)
      compare_stata_grid(stem)
    })
  })
}

# ---------------------------------------------------------------------------
# The tolerance-free comparison
# ---------------------------------------------------------------------------
#
# Every test above compares two REML fits and therefore cannot be tighter than
# the two optimisers agree. This one takes the fit out: it reads Stata's own
# variance components from grid 5, feeds them to slope_params_manual(), and
# requires the arithmetic to agree to machine precision. If this passes and the
# grids above only miss by 1e-5, the port is right and the optimisers differ.
# If this fails, the port is wrong, whatever the grids above say.

# How far apart the two REML implementations land on each fitted quantity, at
# nlme's stock settings. Measured over all 37 fits in grid 5; the maxima are
# the worst single fit with a little headroom, the medians are the typical row.
#
# Deliberately calibrated to nlme's DEFAULTS rather than to a tightened
# lmeControl. Tightening closes most of the gap (var_tte on slpower1 moves from
# 8.36388645 to 8.36394022 against Stata's 8.36394172) but it changes what the
# package does for every user to make a test look better, which is backwards.
# These numbers describe the port as shipped.
#
# The maxima are dominated by fits that barely identify a random slope:
# F3-id<=40|id>110 is 80 subjects at three visits, and F3-visit<2 and
# F1-visit==0|visit==3 have two timepoints per subject. The slope variance is
# the worst-determined component in exactly those fits, which is why its
# tolerance is the loosest one here.
TOL_FIT <- c(intercept = 1e-3, slope_var = 5e-2, cov = 5e-2,
             residual = 1e-2, slope = 1e-3)
TOL_FIT_MEDIAN <- c(intercept = 1e-5, slope_var = 1e-4, cov = 5e-4,
                    residual = 1e-4, slope = 1e-6)

test_that("the fitted variance components match Stata's", {
  skip_without_stata_reference("05_variance_components")

  fits <- read_stata_reference("05_variance_components")
  expect_equal(nrow(fits), 37)

  d <- do.call(rbind, lapply(seq_len(nrow(fits)), function(i) {
    row <- fits[i, , drop = FALSE]
    p <- stata_fit(row)
    # Stata fitted on time/scale; put the R components on that axis too.
    f <- stata_axis_factor(row)
    data.frame(
      vtag      = row$vtag,
      intercept = rel_diff(row$v_intercept, p$sigma2_intercept),
      slope_var = rel_diff(row$v_slope,     p$sigma2_slope * f^2),
      cov       = rel_diff(row$cov_islope,  p$sigma_cov * f),
      residual  = rel_diff(row$v_residual,  p$sigma2_residual),
      slope     = rel_diff(row$slope,       p$slope * f),
      stringsAsFactors = FALSE
    )
  }))

  for (q in names(TOL_FIT)) {
    bad <- d[d[[q]] > TOL_FIT[[q]], c("vtag", q), drop = FALSE]
    expect_equal(nrow(bad), 0,
                 info = stata_mismatch_report(bad, q, paste("component", q)))
    expect_lt(stats::median(d[[q]]), TOL_FIT_MEDIAN[[q]])
  }

  # The comparator slope where there is one. slpower2's model-1 fit is the
  # unbalanced one, so this is the strictest real test of the fixed effects.
  hc <- fits[!is.na(fits$slope_comp), , drop = FALSE]
  expect_gt(nrow(hc), 0)
  for (i in seq_len(nrow(hc))) {
    row <- hc[i, , drop = FALSE]
    p <- stata_fit(row)
    expect_lt(rel_diff(row$slope_comp, p$slope_comparator * stata_axis_factor(row)),
              TOL_FIT[["slope"]])
  }
})

test_that("with Stata's own variance components the arithmetic is exact", {
  skip_without_stata_reference("05_variance_components")

  fits <- read_stata_reference("05_variance_components")
  expect_true(all(fits$converged == 1))

  # Every design grid, grid 4 included: the loose N bound above is excused on
  # the grounds that this exact test carries the real burden, which is only
  # true if the grid it excuses reaches this loop.
  # Per grid, not for the test as a whole: skipping the whole test because one
  # grid is absent would silently drop the exact comparison for every grid that
  # is present, which is the opposite of what a missing CSV should cost.
  compared <- 0
  for (stem in names(DESIGN_GRIDS)) {
    if (!have_stata_reference(stem)) next
    compared <- compared + 1
    grid <- read_stata_reference(stem)
    grid <- grid[grid$rc == 0 & !grid$tag %in% c(STATA_ONLY, KNOWN_DIVERGENCES),
                 , drop = FALSE]
    j <- stata_join_fits(grid, fits)
    expect_gt(nrow(j$grid), 0)

    res <- lapply(seq_len(nrow(j$grid)), function(i) {
      row <- j$grid[i, , drop = FALSE]
      p <- stata_manual_params(j$fits[i, , drop = FALSE])
      # Stata's own axis: the schedule as integers, no rescaling.
      visits <- c(0, as.numeric(strsplit(row$sched, "[[:space:]]+")[[1]]))
      des <- suppressWarnings(trial_design(visits, stata_dropout(row)))
      args <- list(params = p, design = des, alpha = row$alpha)
      if (isTRUE(row$usetrt == 1)) args$target <- "observed"
      else if (nzchar(row$effin)) args$effectiveness <- as.numeric(row$effin)
      # No tryCatch. Every row reaching here survived the filter above, so both
      # implementations accept it and a refusal is a finding: letting the error
      # surface names the call and the reason, where catching it would drop the
      # row and leave the suite green over one row fewer. It also keeps `res`
      # row-for-row with j$grid, so the mask below can be built from the frame
      # it indexes rather than from a different object of assumed equal length.
      out <- suppressWarnings(
        if (row$mode == "n") do.call(slope_power, c(args, list(n = row$nin)))
        else do.call(slope_sample_size, c(args, list(power = row$pin)))
      )
      data.frame(tag = row$tag, drops = row$drops, mode = row$mode,
                 o_n = row$o_n, r_n = out$n,
                 o_power = row$o_power, r_power = out$power,
                 o_var_tte = row$o_var_tte, r_var_tte = out$var_tte,
                 stringsAsFactors = FALSE)
    })
    res <- do.call(rbind, res)

    bad <- res[!is.na(res$o_n) & res$o_n != res$r_n, , drop = FALSE]
    expect_equal(nrow(bad), 0,
                 info = stata_mismatch_report(bad, c("o_n", "r_n"),
                                              paste0(stem, ": exact N")))

    bad <- res[!is.na(res$o_power) &
                 abs(res$o_power - res$r_power) > 1e-12, , drop = FALSE]
    expect_equal(nrow(bad), 0,
                 info = stata_mismatch_report(bad, c("o_power", "r_power"),
                                              paste0(stem, ": exact power")))

    # var_tte only where the two are defined the same way, per CONTRACT.md 5.6.
    v <- res[stata_var_tte_comparable(res) & !is.na(res$o_var_tte), ,
             drop = FALSE]
    bad <- v[rel_diff(v$o_var_tte, v$r_var_tte) > 1e-12, , drop = FALSE]
    expect_equal(nrow(bad), 0,
                 info = stata_mismatch_report(bad, c("o_var_tte", "r_var_tte"),
                                              paste0(stem, ": exact var_tte")))
  }

  # Grid 5 was present, so at least one design grid must have been too.
  expect_gt(compared, 0)
})

# ---------------------------------------------------------------------------
# Checks that need no Stata run: internal invariants the grids make visible
# ---------------------------------------------------------------------------

test_that("scale() has an exact analogue in real-valued visit times", {
  d <- load_paper_data("slpower1")
  p_year <- slope_params(sdmt ~ time | id, d)

  d_half <- d
  d_half$time <- d_half$time / 0.5
  p_half <- slope_params(sdmt ~ time | id, d_half)

  # Stata's scale(0.5) schedule(1 2 3 4) against the R port's visits at
  # 0, 0.5, 1, 1.5, 2 on the unscaled axis. The two fits are algebraically the
  # same model; only the time unit differs.
  a <- slope_sample_size(p_year, c(0, 0.5, 1, 1.5, 2), effectiveness = 0.33)
  b <- slope_sample_size(p_half, 0:4, effectiveness = 0.33)

  expect_equal(a$n, b$n)
  expect_equal(a$effect_size, b$effect_size, tolerance = 1e-8)
  # The half-year axis halves the slope: the paper prints -1.672 on the year
  # axis and -0.836 on the rescaled one (p.588, p.589).
  expect_equal(a$slope_difference, b$slope_difference * 2, tolerance = 1e-8)
})
