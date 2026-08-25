# Layer 4 --- slope_bootstrap().
#
# Replaces the Stata incantation
#   bootstrap r(sampsize), cluster(id) idcluster(id2) strata(case) bca
#             jack(n(r(obs_in_model))): slopepower ...
# with resampling that is correct by construction: subjects, not rows, are the
# unit of resampling.
#
# R is kept small throughout; these test the machinery, not the coverage.

# A slope flattened to just under the section 2.6 threshold, so that a handful
# of replicates genuinely refit a slope of the other sign. Three tests need it;
# built on demand rather than at load time so that a fit is not paid for by
# every other test in the file.
flat_fit <- function() {
  d <- load_paper_data("slpower1")
  d$sdmt <- d$sdmt + 1.66 * d$time
  slope_params(sdmt ~ time | id, d)
}

test_that("slope_bootstrap() returns replicates with a spread, for the slope", {
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 15, seed = 1))
  expect_true(is.list(b))
  expect_equal(b$statistic, "slope")
  reps <- b$replicates
  expect_true(!is.null(reps))
  expect_gt(stats::sd(reps, na.rm = TRUE), 0)
  # the replicate slopes should surround the observed one
  expect_lt(abs(mean(reps, na.rm = TRUE) - paper_fit("slpower1")$slope), 0.5)
})

test_that("slope_bootstrap() reports the mean, SD and straddle of its replicates", {
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 15, seed = 1))
  expect_equal(b$boot_mean, mean(b$replicates))
  expect_equal(b$boot_sd, stats::sd(b$replicates))
  # A well-estimated decline slope: no replicate should cross zero.
  expect_equal(b$straddle, 0)
})

test_that("slope_bootstrap() straddle counts replicates of the opposite sign", {
  # Flattened as in the section 2.6 warning test above, so some replicate
  # slopes land on the other side of zero from the (small negative) observed
  # slope.
  b <- suppressWarnings(slope_bootstrap(flat_fit(), R = 30, seed = 1))
  expect_gt(b$straddle, 0)
  # For `statistic = "slope"` the replicates *are* the refitted slopes, so the
  # definition below coincides with the general one.
  expect_equal(b$straddle,
              mean(sign(b$replicates) != sign(b$observed)))
})

test_that("slope_bootstrap() straddle is measured on the slope, not the statistic", {
  # The point of the check: a sample size is a positive integer, so its own
  # sign can never straddle zero. What straddles is the slope each replicate
  # refits, and the same resampling (same seed, same subjects drawn) must
  # therefore report the same proportion whichever statistic is asked for.
  flat <- flat_fit()
  ref <- suppressWarnings(slope_bootstrap(flat, R = 30, seed = 1))
  ss <- suppressWarnings(
    slope_sample_size(flat, c(0, 1, 2), effectiveness = 0.33))
  b <- suppressWarnings(
    slope_bootstrap(ss, R = 30, type = "percentile", seed = 1))

  expect_true(all(b$replicates > 0))          # nothing to straddle here
  expect_equal(mean(sign(b$replicates) != sign(b$observed)), 0)
  expect_gt(b$straddle, 0)
  expect_equal(b$straddle, ref$straddle)
})

test_that("slope_bootstrap() is reproducible under a fixed seed", {
  a <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 10, seed = 42))
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 10, seed = 42))
  expect_equal(a$replicates, b$replicates)
})

test_that("slope_bootstrap() dispatches on what it is handed", {
  p <- paper_fit("slpower1")
  ss <- slope_sample_size(p, c(0, 1, 2), effectiveness = 0.33)
  pw <- slope_power(p, c(0, 1, 2), n = 712, effectiveness = 0.33)

  expect_equal(suppressWarnings(slope_bootstrap(ss, R = 3, seed = 1))$statistic, "n")
  expect_equal(suppressWarnings(slope_bootstrap(pw, R = 3, seed = 1))$statistic, "power")
  expect_equal(suppressWarnings(slope_bootstrap(p, R = 3, seed = 1))$statistic, "slope")

  # anything else is refused with a pointer to what is accepted
  expect_error(slope_bootstrap(data.frame(x = 1)), "cannot bootstrap an object")
  expect_error(slope_bootstrap(42), "slope_sample_size")
})

test_that("slope_bootstrap() can bootstrap the sample size itself", {
  ss <- slope_sample_size(paper_fit("slpower1"), trial_design(c(0, 1, 2)),
                          effectiveness = 0.33)
  b <- suppressWarnings(slope_bootstrap(ss, R = 15, seed = 7))
  reps <- b$replicates
  expect_true(all(reps > 0, na.rm = TRUE))
  expect_gt(stats::sd(reps, na.rm = TRUE), 0)
  expect_equal(b$observed, ss$n)
})

test_that("slope_bootstrap() validates R and level", {
  p <- paper_fit("slpower1")
  expect_error(slope_bootstrap(p, R = 0), "R")
  expect_error(slope_bootstrap(p, R = 5, level = 1), "level")
  expect_error(slope_bootstrap(p, R = 5, level = 0), "level")
})

test_that("slope_bootstrap() warns when the slope is weak relative to its error", {
  # Paper section 2.6: if |slope| / se(slope) < 2.5 the replicate slopes can
  # straddle zero and the resulting interval is meaningless.
  d <- load_paper_data("slpower1")
  d$sdmt <- d$sdmt + 1.66 * d$time         # flatten the slope, keep the structure
  flat <- slope_params(sdmt ~ time | id, d)
  expect_warning(slope_bootstrap(flat, R = 5, seed = 1))
})

# --- resample_frame() ---------------------------------------------------
# `resample_frame` is internal -- these exercise the resampling primitive
# directly rather than through a full fit, since a single-subject stratum
# needs a stage-one fit that has since been rejected by slope_params() itself
# (see test-params.R's "comparator group of 1 subject" regression) for
# healthy/treated data, but resample_frame() itself must still be correct for
# any stratum it is handed.
resample_frame <- slopepower:::resample_frame

test_that("resample_frame() always returns the sole member of a size-1 stratum", {
  # sample(x, size) treats a length-one x as a range to draw *from* (1:x)
  # rather than a value to draw *with replacement*, so a stratum whose lone
  # *position* is not 1 used to have its "resample" drawn from 1:pos instead
  # of always returning pos. (A lone position of exactly 1 cannot expose this,
  # since 1:1 = 1 either way -- hence subject 3 below, not subject 1.)
  frame <- data.frame(marker = rep(c("A", "B", "C"), each = 2),
                      time = rep(c(0, 1), 3), subject = rep(1:3, each = 2))
  subject_index <- split(seq_len(nrow(frame)), frame$subject)   # positions 1,2,3
  one_stratum <- 3L   # subject 3 ("C") is the sole member of its stratum

  markers <- vapply(1:200, function(i) {
    unique(resample_frame(frame, subject_index, list(one_stratum))$marker)
  }, character(1L))
  expect_true(all(markers == "C"))
})

test_that("resample_frame() keeps each replicate's stratum sizes fixed", {
  frame <- data.frame(y = rnorm(6), time = rep(c(0, 1), 3), subject = rep(1:3, each = 2))
  subject_index <- split(seq_len(nrow(frame)), frame$subject)
  groups <- list(3L, 1:2)   # stratum sizes 1 (position 3) and 2 (positions 1:2)
  for (i in 1:20) {
    out <- resample_frame(frame, subject_index, groups)
    expect_equal(length(unique(out$subject)), 3L)
    expect_equal(nrow(out), nrow(frame))
  }
})

test_that("slope_bootstrap() refuses parameters that carry no data", {
  manual <- ref_params()
  expect_error(slope_bootstrap(manual, R = 5), "no fitted model")
  expect_error(
    slope_bootstrap(slope_sample_size(manual, c(0, 1, 2), effectiveness = 0.33), R = 5),
    "no fitted model")
})

# --- failed replicates, and the progress ticks --------------------------
# A resample can be too degenerate to refit. With only a handful of subjects a
# draw that picks the same person two or three times leaves those "people"
# identical, so there is no between-subject variance for `nlme` to estimate and
# the replicate is lost. run_bootstrap() drops those and reports how many,
# which is the caller's only signal that the interval rests on fewer replicates
# than they asked for -- `$replicates` holds the survivors, so `R` cannot be
# recovered from its length.
#
# Three subjects at four visits is one more than the smallest fit
# slope_params() will accept (it requires 2 participants with repeat visits),
# and small enough that degenerate draws are common: between a fifth and a
# third of resamples fail. Its slope is 3.09 times its standard error,
# comfortably over the section 2.6 threshold, so the expectations below see the
# discard warning and nothing else. Built once, at file scope, for the same
# reason as helper-fits.R's cached fits: the mixed-model refits are the slow
# part.
boot_tiny_fit <- local({
  set.seed(11)
  d <- data.frame(id = rep(1:3, each = 4), visit = rep(0:3, 3))
  d$y <- rep(c(50, 40, 45), each = 4) +
    rep(c(-2, -1.5, -1), each = 4) * d$visit + rnorm(12, 0, 2)
  slope_params(y ~ visit | id, d)
})

test_that("slope_bootstrap() discards replicates that fail to refit, and says so", {
  w <- capture_warnings(
    b <- slope_bootstrap(boot_tiny_fit, R = 10, seed = 1, type = "percentile"))
  expect_gt(b$n_failed, 0L)
  # The count in the message is the count on the object -- the number that
  # tells the caller how much of the requested R they actually got. `R` itself
  # records what was asked for, not what survived.
  expect_match(w, sprintf("%d of 10 replicates failed to converge and were discarded",
                          b$n_failed), fixed = TRUE)
  expect_equal(b$R, 10)
  expect_equal(length(b$replicates), 10L - b$n_failed)
  expect_true(all(is.finite(b$replicates)))
})

test_that("slope_bootstrap() refuses to form an interval from under two replicates", {
  # Two is the minimum quantile() can interpolate between, so R = 1 reaches the
  # guard with nothing having gone wrong at all: the message counts zero
  # failures out of one. Asserted on a well-behaved fit precisely because it
  # depends on no convergence behaviour whatsoever.
  expect_error(slope_bootstrap(paper_fit("slpower1"), R = 1, seed = 1),
               "0 of 1 replicates failed; not enough succeeded")
  # And the case the message is actually written for: a seed under which every
  # resample of the three-subject fit is degenerate. The count is left out of
  # the pattern -- which replicates fail is the optimiser's business, that too
  # few did is the contract.
  expect_error(slope_bootstrap(boot_tiny_fit, R = 3, seed = 7, type = "percentile"),
               "of 3 replicates failed; not enough succeeded")
})

test_that("slope_bootstrap(progress = TRUE) ticks every max(1, R/10) replicates", {
  # The max() is load-bearing: without it any R < 10 gives a tick of zero, and
  # `b %% 0` is NaN, so the `== 0L` test is NA and the loop errors rather than
  # merely printing nothing. R = 5 therefore ticks on every replicate.
  msgs <- capture_messages(suppressWarnings(
    slope_bootstrap(boot_tiny_fit, R = 5, seed = 1, type = "percentile",
                    progress = TRUE)))
  expect_length(msgs, 5L)
  expect_match(msgs[5L], "5 of 5 replicates")

  # R = 20 ticks every second replicate: ten messages either way, which is the
  # point of the R/10.
  msgs20 <- capture_messages(suppressWarnings(
    slope_bootstrap(boot_tiny_fit, R = 20, seed = 1, type = "percentile",
                    progress = TRUE)))
  expect_length(msgs20, 10L)
  expect_match(msgs20[1L], "2 of 20 replicates")
  expect_match(msgs20[10L], "20 of 20 replicates")

  # Off by default. Replicate failures are warnings, not messages, so
  # suppressing them leaves nothing behind -- the ticks above are the only
  # thing that passes through the suppressWarnings() around the resample loop.
  expect_no_message(suppressWarnings(
    slope_bootstrap(boot_tiny_fit, R = 5, seed = 1, type = "percentile")))
})

# --- bca_interval() and its fallbacks -----------------------------------
# The acceleration needs a leave-one-subject-out jackknife, matching the
# clustering of the bootstrap itself, and there are fits for which it cannot be
# had. bca_interval() answers NULL in each such case and run_bootstrap() keeps
# the percentile interval it has already computed, so the observable
# consequence of every fallback is the same: `type = "bca"` in, `$type ==
# "percentile"` out.
bca_interval <- slopepower:::bca_interval

# Two subjects is the fewest slope_params() will fit. Leaving one out leaves
# one, which it refuses outright -- 2 participants with repeat visits are
# needed to identify the slope variance -- so no jackknife value is obtainable
# at all. That is a deterministic structural guard, not an optimiser failure.
boot_pair_fit <- local({
  set.seed(5)
  d <- data.frame(id = rep(1:2, each = 6), visit = rep(0:5, 2))
  d$y <- rep(c(50, 42), each = 6) +
    rep(c(-2.4, -1.4), each = 6) * d$visit + rnorm(12, 0, 1.2)
  slope_params(y ~ visit | id, d)
})

test_that("bca_interval() falls back when under three jackknife values survive", {
  frame <- slopepower:::boot_frame(boot_pair_fit, "test")
  subject_index <- split(seq_len(nrow(frame)), frame$subject)
  expect_length(subject_index, 2L)

  # `theta` straddles `observed`, so mean(theta < observed) is 0.4 and the bias
  # correction z0 is perfectly computable: the NULL below is the jackknife's,
  # not the already-covered qnorm(0)/qnorm(1) guard above it.
  expect_null(bca_interval(theta = c(1, 2, 3, 4, 5), observed = 3,
                           frame = frame, subject_index = subject_index,
                           refitter = slopepower:::make_refitter(boot_pair_fit),
                           compute = function(p) p$slope,
                           probs = c(0.025, 0.975)))
})

test_that("slope_bootstrap() reports a percentile interval when BCa cannot be had", {
  expect_warning(
    b <- slope_bootstrap(boot_pair_fit, R = 6, seed = 2, type = "bca"),
    "bias-correction could not be computed")
  expect_equal(b$type, "percentile")
  # No replicate failed under this seed, so the fallback is not a shortage of
  # replicates; and they straddle the observed slope, so it is not the bias
  # correction either. What is missing is the jackknife.
  expect_equal(b$n_failed, 0)
  prop <- mean(b$replicates < boot_pair_fit$slope)
  expect_gt(prop, 0)
  expect_lt(prop, 1)
  # The interval reported is exactly the percentile one already in hand.
  expect_equal(b$ci, stats::quantile(b$replicates, c(0.025, 0.975),
                                     names = FALSE, type = 7))
})

# Ten subjects: enough for the jackknife to produce three or more usable
# values, so BCa succeeds. This is the only bootstrap in this file that reaches
# the accelerated branch, and at 12 replicates plus 10 leave-one-out refits it
# is the most expensive thing in it -- hence once, at file scope. Its
# discard warning is asserted on `boot_tiny_fit` above rather than here.
boot_bca_fit <- local({
  set.seed(4)
  subj <- data.frame(id = 1:10, a = rnorm(10, 50, 8), b = rnorm(10, -2, 0.6))
  d <- merge(subj, data.frame(visit = 0:3))
  d$y <- d$a + d$b * d$visit + rnorm(nrow(d), 0, 2)
  slope_params(y ~ visit | id, d)
})
boot_bca <- suppressWarnings(
  slope_bootstrap(boot_bca_fit, R = 12, seed = 2, type = "bca"))

test_that("slope_bootstrap() returns a genuine BCa interval when it can", {
  expect_equal(boot_bca$type, "bca")
  expect_length(boot_bca$ci, 2L)
  expect_lt(boot_bca$ci[1L], boot_bca$observed)
  expect_gt(boot_bca$ci[2L], boot_bca$observed)
  # Shifted relative to the percentile interval on the same replicates, which
  # is the whole reason the paper recommends BCa: the distribution of the
  # bootstrapped quantity is skewed and the plain quantiles are biased.
  pctl <- stats::quantile(boot_bca$replicates, c(0.025, 0.975),
                          names = FALSE, type = 7)
  expect_false(isTRUE(all.equal(boot_bca$ci, pctl)))
})

# --- printing ------------------------------------------------------------

test_that("print.slope_bootstrap() labels a BCa interval and counts the failures", {
  expect_output(print(boot_bca), "<slope_bootstrap>")
  out <- capture.output(print(boot_bca))
  expect_true(any(grepl("Statistic:   slope", out, fixed = TRUE)))
  expect_true(any(grepl("95% BCa CI:", out, fixed = TRUE)))
  # The failure count is appended to the replicate line only when there is one,
  # so that a clean run does not read as though something went wrong.
  expect_true(any(grepl(sprintf("%d used, %d failed", length(boot_bca$replicates),
                                boot_bca$n_failed), out, fixed = TRUE)))
  # Returns its argument invisibly, so `b <- print(b)` is safe and a bare
  # `print(b)` at the prompt does not print the object twice.
  capture.output(res <- withVisible(print(boot_bca)))
  expect_false(res$visible)
  expect_identical(res$value, boot_bca)
})

test_that("print.slope_bootstrap() labels a percentile interval and omits a zero count", {
  b <- slope_bootstrap(boot_pair_fit, R = 6, seed = 2, type = "percentile")
  expect_equal(b$n_failed, 0)
  out <- capture.output(print(b))
  expect_true(any(grepl("95% percentile CI:", out, fixed = TRUE)))
  expect_true(any(grepl("6 used", out, fixed = TRUE)))
  expect_false(any(grepl("failed", out, fixed = TRUE)))
})

test_that("print.slope_bootstrap() always reports the bootstrap mean and SD", {
  out <- capture.output(print(boot_bca))
  expect_true(any(grepl(sprintf("mean %s, SD %s",
                                format(boot_bca$boot_mean, digits = 6),
                                format(boot_bca$boot_sd, digits = 6)),
                        out, fixed = TRUE)))
})

test_that("print.slope_bootstrap() notes straddling only when it occurs", {
  expect_false(any(grepl("opposite side of zero", capture.output(print(boot_bca)),
                         fixed = TRUE)))

  b <- suppressWarnings(slope_bootstrap(flat_fit(), R = 30, seed = 1))
  out <- capture.output(print(b))
  expect_true(any(grepl(sprintf("%.1f%% of replicates", 100 * b$straddle),
                        out, fixed = TRUE)))
})

# --- the random-effects structure held fixed across replicates -----------
# `common_variance` matters only under `healthy`, and only through the
# controls' slope: the model factorises per group, so the case estimates are
# invariant either way. It is therefore `slope_comparator` --- and every
# stage-two answer measured against it --- that moves when a replicate fits a
# different structure from the one the point estimate came from.

# Ragged control follow-up, deliberately: with balanced complete data GLS
# coincides with OLS whatever the covariance structure, and the two fits agree
# to every digit, so a balanced fixture could not tell the structures apart.
ragged_healthy <- function() {
  set.seed(2)
  subj <- data.frame(id = 1:80, case = rep(c(1, 0), each = 40))
  subj$intercept <- stats::rnorm(80, 50, 10)
  subj$slope <- stats::rnorm(80, ifelse(subj$case == 1, -1.7, -0.3), 1.4)
  d <- merge(subj, data.frame(visit = 0:3))
  d$sdmt <- d$intercept + d$slope * d$visit + stats::rnorm(nrow(d), 0, 3)
  set.seed(11)
  d[!(d$case == 0 & d$visit > 0 & stats::runif(nrow(d)) < 0.35), ]
}

test_that("make_refitter() pins the structure the observed fit ended up with", {
  d <- ragged_healthy()
  full <- suppressWarnings(
    slope_params(sdmt ~ visit | id, d, healthy = case, common_variance = FALSE))
  red <- suppressWarnings(
    slope_params(sdmt ~ visit | id, d, healthy = case, common_variance = TRUE))
  expect_false(full$common_variance)
  expect_true(red$common_variance)

  expect_equal(
    slopepower:::make_refitter(full)(slopepower:::boot_frame(full, "test"))$common_variance,
    FALSE)
  expect_equal(
    slopepower:::make_refitter(red)(slopepower:::boot_frame(red, "test"))$common_variance,
    TRUE)

  # Not merely a label: the two structures disagree about the controls' slope,
  # which is what makes pinning it worth doing. The case slope is invariant.
  expect_equal(full$slope, red$slope)
  expect_false(isTRUE(all.equal(full$slope_comparator, red$slope_comparator)))
})

test_that("slope_bootstrap() replicates keep the reduced structure when it was used", {
  d <- ragged_healthy()
  red <- suppressWarnings(
    slope_params(sdmt ~ visit | id, d, healthy = case, common_variance = TRUE))

  seen <- c()
  refitter <- slopepower:::make_refitter(red)
  frame <- slopepower:::boot_frame(red, "test")
  subject_index <- split(seq_len(nrow(frame)), frame$subject)
  sub_group <- vapply(subject_index, function(r) frame$group[r[1L]], numeric(1L))
  groups <- unname(split(seq_along(subject_index), sub_group))

  set.seed(4)
  for (i in seq_len(5L)) {
    p <- refitter(slopepower:::resample_frame(frame, subject_index, groups))
    seen <- c(seen, p$common_variance)
  }
  expect_true(all(seen))
})

test_that("make_refitter() leaves common_variance off the non-healthy calls", {
  # slope_params() warns that `common_variance` is ignored without `healthy`;
  # passing it anyway would earn that warning once per replicate.
  for (fit in list(paper_fit("slpower1"), paper_fit("slpower3"))) {
    cl <- environment(slopepower:::make_refitter(fit))$cl
    expect_false("common_variance" %in% names(as.list(cl)))
  }
})

# --- slope_se() ----------------------------------------------------------
# The section 2.6 check inside run_bootstrap() is slope_se()'s only caller in
# the package, so these pin its contract directly rather than through it.

test_that("slope_se() is the standard error of the slope the object reports", {
  # comparator "none": the slope is a single fixed effect, so its standard
  # error is that coefficient's own.
  p1 <- paper_fit("slpower1")
  V1 <- as.matrix(stats::vcov(p1$fit))
  expect_equal(slope_se(p1), sqrt(V1["sp_time", "sp_time"]))

  # "healthy": params.R takes the case slope as b[sp_time] + b[sp_case:sp_time],
  # so its standard error needs both variances *and* twice the covariance. The
  # interaction's name is resolved rather than spelled out, for the reason
  # test-term-order.R exists. Neither coefficient's own standard error would
  # do -- here the sum is estimated more precisely than the interaction alone,
  # so taking a single diagonal element would overstate the uncertainty and
  # fire section 2.6 on a slope that does not deserve it.
  p2 <- paper_fit("slpower2")
  b2 <- nlme::fixef(p2$fit)
  V2 <- as.matrix(stats::vcov(p2$fit))
  int2 <- slopepower:::resolve_fixef_name(b2, c("sp_case", "sp_time"))
  k2 <- as.numeric(names(b2) %in% c("sp_time", int2))
  expect_equal(slope_se(p2), sqrt(drop(k2 %*% V2 %*% k2)))
  expect_lt(slope_se(p2), sqrt(V2[int2, int2]))

  # "treated": the reported slope is the control arm's,
  # b[sp_time] + b[sp_placebo_time].
  p3 <- paper_fit("slpower3")
  b3 <- nlme::fixef(p3$fit)
  V3 <- as.matrix(stats::vcov(p3$fit))
  k3 <- as.numeric(names(b3) %in% c("sp_time", "sp_placebo_time"))
  expect_equal(slope_se(p3), sqrt(drop(k3 %*% V3 %*% k3)))

  for (p in list(p1, p2, p3)) {
    expect_length(slope_se(p), 1L)
    expect_gt(slope_se(p), 0)
  }
})

test_that("slope_se() returns NA for parameters that carry no fitted model", {
  # slope_params_manual() objects hold parameter values and nothing else
  # (CONTRACT.md section 2: `fit` is NULL), so there is no covariance matrix to
  # read a standard error off. NA rather than an error, because the section 2.6
  # check has to be skippable: manual parameters are a legitimate input to the
  # rest of the package.
  manual <- ref_params()
  expect_identical(slope_se(manual), NA_real_)
  expect_error(slope_se(42), "must be a")
})

test_that("slope_se() is the denominator of the ratio slope_bootstrap() warns on", {
  # CONTRACT.md section 6 / paper section 2.6: the warning fires exactly when
  # abs(slope) / slope_se(params) < 2.5, and the ratio it quotes is that one.
  # Pinning the printed number to slope_se() is the point -- the comment in
  # bootstrap.R records that a silently NA standard error once switched the
  # check off altogether, which no test of the warning's mere presence sees.
  strong <- paper_fit("slpower1")
  expect_gt(abs(strong$slope) / slope_se(strong), 2.5)
  expect_no_warning(slope_bootstrap(strong, R = 3, seed = 1, type = "percentile"))

  d <- load_paper_data("slpower1")
  d$sdmt <- d$sdmt + 1.66 * d$time         # flatten the slope, keep the structure
  flat <- slope_params(sdmt ~ time | id, d)
  ratio <- abs(flat$slope) / slope_se(flat)
  expect_lt(ratio, 2.5)
  expect_warning(slope_bootstrap(flat, R = 3, seed = 1, type = "percentile"),
                 sprintf("only %.2f times its standard error", ratio))
})
