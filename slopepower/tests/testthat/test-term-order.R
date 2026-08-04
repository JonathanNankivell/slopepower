# Invariance to the order of terms in the model formula.
#
# The public formula is `outcome ~ time | subject`, whose right-hand side is a
# single term, so there is nothing for a user to permute there. The concern is
# really about the fixed-effect formulae `params.R` builds internally --
# `sp_y ~ sp_time + sp_placebo_time` and `sp_y ~ sp_case * sp_time`, the
# analogues of `sdmt ~ visit + visit:treat + (visit | id)` -- where reordering
# changes the design matrix column order and, for interactions, the coefficient
# *name*. The fit is mathematically identical either way; extraction must be too.

# `fixef_term()` is internal -- it is the extraction mechanism under test, not
# part of the public surface.
fixef_term <- slopepower:::fixef_term

sp_prep2 <- function() {
  d <- load_paper_data("slpower2")
  d$sp_y <- d$sdmt
  d$sp_time <- ave(d$time, d$id, FUN = function(x) x - min(x))
  d$sp_subject <- factor(d$id)
  d$sp_case <- as.integer(d$case)
  d
}

sp_prep3 <- function() {
  d <- load_paper_data("slpower3")
  d$sp_y <- d$sdmt
  d$sp_time <- d$time
  d$sp_subject <- factor(d$id)
  d$sp_treat <- as.integer(d$treat)
  d$sp_placebo_time <- (1 - d$sp_treat) * d$sp_time
  d
}

test_that("reordering main effects leaves the fit and the extraction unchanged", {
  d <- sp_prep3()

  forms <- list(
    sp_y ~ sp_time + sp_placebo_time,
    sp_y ~ sp_placebo_time + sp_time
  )
  fits <- lapply(forms, function(f) {
    nlme::lme(f, random = ~ sp_time | sp_subject, data = d, method = "REML")
  })

  # Same model, so the same likelihood -- only the column order differs.
  expect_equal(as.numeric(stats::logLik(fits[[1]])),
               as.numeric(stats::logLik(fits[[2]])), tolerance = 1e-10)
  expect_setequal(names(nlme::fixef(fits[[1]])), names(nlme::fixef(fits[[2]])))

  # The extraction params.R performs must agree despite the reordering.
  got <- lapply(fits, function(f) {
    b <- nlme::fixef(f)
    c(slope      = fixef_term(b, "sp_time", "t") + fixef_term(b, "sp_placebo_time", "t"),
      comparator = fixef_term(b, "sp_time", "t"))
  })
  # Not bit-identical: REML is iterative, and reordering the design matrix
  # columns changes the optimiser's path. Agreement is limited by the
  # convergence tolerance (`lmeControl(tolerance = 1e-7)`), not by the algebra;
  # observed differences are around 4e-10.
  expect_equal(got[[1]], got[[2]], tolerance = 1e-7)
  expect_equal(unname(got[[1]]["slope"]), -1.8517, tolerance = 1e-3)
})

test_that("an interaction resolves whichever way round R spells it", {
  d <- sp_prep2()

  # `sp_case * sp_time` names the interaction sp_case:sp_time; `sp_time *
  # sp_case` names it sp_time:sp_case. Both orderings of the equivalent
  # expanded form appear too, since that is how a reordering would arrive in
  # practice.
  forms <- list(
    "case * time"             = sp_y ~ sp_case * sp_time,
    "time * case"             = sp_y ~ sp_time * sp_case,
    "time + case + case:time" = sp_y ~ sp_time + sp_case + sp_case:sp_time,
    "case:time + time + case" = sp_y ~ sp_case:sp_time + sp_time + sp_case
  )
  fits <- lapply(forms, function(f) {
    nlme::lme(f, random = ~ sp_time | sp_subject, data = d, method = "REML")
  })

  ll <- vapply(fits, function(f) as.numeric(stats::logLik(f)), numeric(1))
  expect_equal(max(ll) - min(ll), 0, tolerance = 1e-8)

  # Both spellings genuinely occur -- otherwise this test proves nothing.
  spellings <- vapply(fits, function(f) grep(":", names(nlme::fixef(f)), value = TRUE),
                      character(1))
  expect_true(all(c("sp_case:sp_time", "sp_time:sp_case") %in% spellings))

  got <- lapply(fits, function(f) {
    b <- nlme::fixef(f)
    c(cases    = fixef_term(b, "sp_time", "t") + fixef_term(b, c("sp_case", "sp_time"), "t"),
      controls = fixef_term(b, "sp_time", "t"))
  })
  # See the note above on why this is 1e-7 rather than exact.
  for (i in seq_along(got)[-1]) expect_equal(got[[i]], got[[1]], tolerance = 1e-7)
})

test_that("fixef_term is order-invariant where positional indexing is not", {
  b <- c("(Intercept)" = 1, "sp_time" = -1.7, "sp_case:sp_time" = 0.5)
  b_perm <- c("(Intercept)" = 1, "sp_time:sp_case" = 0.5, "sp_time" = -1.7)

  expect_equal(fixef_term(b, c("sp_case", "sp_time"), "t"),
               fixef_term(b_perm, c("sp_case", "sp_time"), "t"))
  expect_equal(fixef_term(b, "sp_time", "t"), fixef_term(b_perm, "sp_time", "t"))

  # The negative control: position 3 is the interaction in one and the slope in
  # the other. This is the failure mode the Stata original is exposed to via its
  # hard-coded e(b) indices, and the reason extraction here is by name.
  expect_false(isTRUE(all.equal(unname(b[3]), unname(b_perm[3]))))

  expect_error(fixef_term(b, "sp_nonesuch", "ctx"), "not found in the fitted model")
})

test_that("results do not depend on column or row order of the data", {
  d <- load_paper_data("slpower2")
  ref <- suppressMessages(slope_params(sdmt ~ time | id, d, healthy = case))

  rev_cols <- suppressMessages(
    slope_params(sdmt ~ time | id, d[, rev(names(d)), drop = FALSE], healthy = case))
  set.seed(42)
  shuffled <- suppressMessages(
    slope_params(sdmt ~ time | id, d[sample(nrow(d)), , drop = FALSE], healthy = case))

  for (got in list(rev_cols, shuffled)) {
    expect_equal(got$slope, ref$slope, tolerance = 1e-8)
    expect_equal(got$slope_comparator, ref$slope_comparator, tolerance = 1e-8)
    expect_equal(got$sigma2_intercept, ref$sigma2_intercept, tolerance = 1e-6)
    expect_equal(got$sigma2_slope, ref$sigma2_slope, tolerance = 1e-6)
    expect_equal(got$sigma_cov, ref$sigma_cov, tolerance = 1e-6)
    expect_equal(got$sigma2_residual, ref$sigma2_residual, tolerance = 1e-6)
  }
})

test_that("a transformation of time is one term, but a covariate is rejected", {
  d <- load_paper_data("slpower2")
  d$vdays <- d$time * 365
  d$age <- rep(50, nrow(d))

  # One term, however written: the day axis rescaled inline must reproduce the
  # year-axis fit. (Only to 3 d.p. -- rescaling is exact algebraically but not
  # numerically, since a day axis leaves sigma^2_b around 1e-5.)
  ref <- suppressMessages(slope_params(sdmt ~ time | id, d, healthy = case))
  inline <- suppressMessages(slope_params(sdmt ~ I(vdays / 365) | id, d, healthy = case))
  expect_equal(inline$slope, ref$slope, tolerance = 1e-3)

  # Two terms is always a mistake: `time + age` used to be evaluated as the
  # arithmetic sum and fitted as though it were the time variable, returning a
  # silently meaningless slope rather than an error.
  expect_error(slope_params(sdmt ~ time + age | id, d, healthy = case),
               "single time term")
  expect_error(slope_params(sdmt ~ time * age | id, d, healthy = case),
               "single time term")
})
