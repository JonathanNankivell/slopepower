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

# --- and does so without reseeding the caller ----------------------------
# `seed` exists to make one result reproducible. It must not also decide what
# every later draw in the caller's script does, which a bare set.seed() inside
# the function would: .Random.seed lives in the global environment.

test_that("a seeded slope_bootstrap() leaves the caller's stream untouched", {
  set.seed(99)
  before <- .Random.seed
  suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 4,
                                   type = "percentile", seed = 7))
  expect_identical(.Random.seed, before)

  # The state being equal is the mechanism; this is what it buys the caller.
  set.seed(99)
  want <- stats::runif(3)
  set.seed(99)
  suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 4,
                                   type = "percentile", seed = 7))
  expect_equal(stats::runif(3), want)
})

test_that("the stream is restored even when the bootstrap errors out", {
  # R = 1 cannot yield the two replicates an interval needs, so this leaves by
  # the stop() rather than by returning -- the reason the restore is an
  # on.exit() and not a line at the end of the function.
  set.seed(99)
  before <- .Random.seed
  expect_error(
    suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 1,
                                     type = "percentile", seed = 7)),
    "not enough succeeded")
  expect_identical(.Random.seed, before)
})

test_that("a seeded bootstrap leaves no .Random.seed where there was none", {
  # A session that has not drawn yet has no .Random.seed at all, and R seeds
  # itself from the clock on first use. Restoring that means removing the
  # variable, not writing a placeholder that would fix the caller's next draw.
  #
  # Emptying the global stream is this test's setup, so it puts back whatever
  # was there on the way out: a later test drawing without a seed of its own
  # must not depend on whether this one ran first.
  outer <- slopepower:::current_seed()
  on.exit(slopepower:::restore_seed(outer), add = TRUE)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 4,
                                   type = "percentile", seed = 7))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("an unseeded slope_bootstrap() still draws from, and advances, the stream", {
  # The other half of the contract: NULL touches nothing, so consecutive
  # unseeded calls must differ rather than repeat.
  set.seed(1)
  a <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 5,
                                        type = "percentile"))
  b <- suppressWarnings(slope_bootstrap(paper_fit("slpower1"), R = 5,
                                        type = "percentile"))
  expect_false(isTRUE(all.equal(a$replicates, b$replicates)))
})

test_that("restore_seed() is total on both states", {
  old <- slopepower:::current_seed()
  on.exit(slopepower:::restore_seed(old), add = TRUE)

  set.seed(5)
  saved <- slopepower:::current_seed()
  expect_type(saved, "integer")
  stats::runif(1)
  slopepower:::restore_seed(saved)
  expect_identical(.Random.seed, saved)

  rm(".Random.seed", envir = globalenv())
  expect_null(slopepower:::current_seed())
  slopepower:::restore_seed(NULL)          # nothing to remove; must not error
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(5)
  slopepower:::restore_seed(NULL)          # now there is; must remove it
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
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

test_that("the default R is one number, not five that agree by inspection", {
  # The generic and its four methods each carry `R` in their own signature, and
  # a caller reaches the default through whichever method dispatch picks. They
  # have to agree, and nothing but this checks that they do.
  defaults <- vapply(
    list(slope_bootstrap,
         slopepower:::slope_bootstrap.slope_sample_size,
         slopepower:::slope_bootstrap.slope_power,
         slopepower:::slope_bootstrap.slope_params,
         slopepower:::slope_bootstrap.default),
    function(f) eval(formals(f)$R), numeric(1L))
  expect_equal(defaults, rep(999, 5L))
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

# --- bca_from_jack() and its fallbacks -----------------------------------
# The acceleration needs a leave-one-subject-out jackknife, matching the
# clustering of the bootstrap itself, and there are fits for which it cannot be
# had. bca_from_jack() answers with a string saying which in each such case and
# run_bootstrap() keeps the percentile interval it has already computed, so the
# observable consequence of every fallback is the same -- `type = "bca"` in,
# `$type == "percentile"` out -- while the warning distinguishes them.
bca_from_jack <- slopepower:::bca_from_jack
jackknife_values <- slopepower:::jackknife_values

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

# The jackknife taken once here and reused across both tests below, exactly as
# run_bootstrap() takes it once per bootstrap and passes the same column to
# every interval() call that needs it: with only two subjects, one left out at
# a time, it can never reach the three values bca_from_jack() requires,
# whatever `theta` and `observed` are being tested.
boot_pair_frame <- slopepower:::boot_frame(boot_pair_fit, "test")
boot_pair_jack <- jackknife_values(
  boot_pair_frame, split(seq_len(nrow(boot_pair_frame)), boot_pair_frame$subject),
  slopepower:::make_refitter(boot_pair_fit), list(function(p) p$slope))[, 1L]

test_that("bca_from_jack() falls back when under three jackknife values survive", {
  expect_length(boot_pair_jack, 2L)

  # `theta` straddles `observed` -- two below, one tied, two above, so the
  # half-counted proportion is 0.5 and z0 is perfectly computable. What fails
  # below is the jackknife, and the returned reason has to say so rather than
  # blaming the bias correction.
  expect_match(bca_from_jack(theta = c(1, 2, 3, 4, 5), observed = 3,
                             jack = boot_pair_jack, probs = c(0.025, 0.975)),
               "acceleration could not be computed")
})

test_that("bca_from_jack() half-counts replicates tied with the observed value", {
  # Every replicate tied with the observed value. Counting ties as "below"
  # makes the proportion 0 and qnorm(0) = -Inf, so this used to be rejected as
  # one-sided; half-counted it is 0.5, meaning no bias to correct, and the run
  # gets as far as the jackknife.
  expect_match(bca_from_jack(theta = rep(3, 5), observed = 3,
                             jack = boot_pair_jack, probs = c(0.025, 0.975)),
               "acceleration could not be computed")

  # Genuinely one-sided is still rejected, and still named as such -- before the
  # jackknife (unused here, since this never reaches it) would even matter.
  expect_match(bca_from_jack(theta = c(4, 5, 6), observed = 3,
                             jack = boot_pair_jack, probs = c(0.025, 0.975)),
               "every replicate falls strictly on one side")
})

test_that("slope_bootstrap() reports a percentile interval when BCa cannot be had", {
  expect_warning(
    b <- slope_bootstrap(boot_pair_fit, R = 6, seed = 2, type = "bca"),
    "acceleration could not be computed")
  expect_equal(b$type, "percentile")
  # No replicate failed under this seed, so the fallback is not a shortage of
  # replicates; and they straddle the observed slope, so it is not the bias
  # correction either. What is missing is the jackknife.
  expect_equal(b$n_failed, 0)
  prop <- (sum(b$replicates < boot_pair_fit$slope) +
             sum(b$replicates == boot_pair_fit$slope) / 2) / length(b$replicates)
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

# --- a sample-size interval stays on the lattice n lives on --------------

test_that("a bootstrapped n is reported as sizes a trial could be run at", {
  pars <- paper_fit("slpower1")
  ss <- suppressWarnings(
    slope_sample_size(pars, c(0, 1, 2), effectiveness = 0.33))
  b <- suppressWarnings(
    slope_bootstrap(ss, R = 40, type = "percentile", seed = 3))

  # Every replicate is 2 * ceiling(...); so is every endpoint. Left as
  # quantile() returns them the interval ran to 963.3 participants, which is
  # not a trial anyone can field.
  expect_true(all(b$replicates %% 2 == 0))
  expect_true(all(b$ci %% 2 == 0))

  # Widened outward, so the reported range always contains the interpolated
  # one: rounding must not be able to narrow a confidence interval.
  interp <- stats::quantile(b$replicates, c(0.025, 0.975),
                            names = FALSE, type = 7)
  expect_lte(b$ci[1L], interp[1L])
  expect_gte(b$ci[2L], interp[2L])
  expect_equal(b$ci, slopepower:::widen_to_lattice(interp, "n"))
})

test_that("widen_to_lattice() widens only outward, and only for n", {
  expect_equal(widen <- slopepower:::widen_to_lattice(c(519.45, 963.3), "n"),
               c(518, 964))
  # Already achievable: nothing moves.
  expect_equal(slopepower:::widen_to_lattice(c(520, 964), "n"), c(520, 964))
  # The continuous statistics are left exactly alone.
  for (st in c("power", "tte", "slope")) {
    expect_equal(slopepower:::widen_to_lattice(c(-1.5, 2.5), st), c(-1.5, 2.5))
  }
})

test_that("the sample-size interval does not swing on floating-point in `probs`", {
  # `probs` is (1 - level) / 2, which for level = 0.95 is one ulp above 0.025.
  # Reading endpoints off with quantile(type = 1) -- order statistics, so
  # lattice points by construction -- makes that ulp decide which side of a
  # step the answer falls on whenever p * B is a whole number, which is exactly
  # what R = 200 or 1000 gives. Interpolating and widening afterwards moves the
  # knife edge onto the lattice, where a perturbation maps a point to itself.
  set.seed(1)
  reps <- 2 * sample(200:600, 1000, replace = TRUE)
  exact <- c(0.025, 0.975)
  fp <- c((1 - 0.95) / 2, 1 - (1 - 0.95) / 2)
  expect_false(identical(exact, fp))

  for (B in c(40L, 1000L)) {
    g <- reps[seq_len(B)]
    widened <- function(p) {
      slopepower:::widen_to_lattice(
        stats::quantile(g, p, names = FALSE, type = 7), "n")
    }
    expect_equal(widened(exact), widened(fp))
    # The construction that was rejected, shown swinging on the same data.
    expect_false(identical(stats::quantile(g, exact, names = FALSE, type = 1),
                           stats::quantile(g, fp, names = FALSE, type = 1)))
  }
})

test_that("the continuous statistics keep the interpolating quantile", {
  b <- suppressWarnings(
    slope_bootstrap(paper_fit("slpower1"), R = 15, type = "percentile", seed = 1))
  expect_equal(b$ci, stats::quantile(b$replicates, c(0.025, 0.975),
                                     names = FALSE, type = 7))
})

# --- printing ------------------------------------------------------------

# One body row of the printed data frame, split into tokens.
#
# `print.data.frame()` separates cells by runs of spaces and prefixes each row
# with its row name, so a row is: row number, statistic, calculated, mean, SD,
# and then the three tokens of the interval ("269", "to", "562") -- the one
# cell with a space inside it, and the reason this splits on whitespace rather
# than trying to recover cell boundaries.
frame_row <- function(out, statistic) {
  row <- out[grepl(paste0("^ *[0-9]+ +", statistic, " "), out)]
  expect_length(row, 1L)
  strsplit(trimws(row), " +")[[1L]]
}

# The frame and nothing else: everything between the class header on the first
# line and the first note. "Counts:" -- which basis a sample size is on --
# prints ahead of "Bootstrap:" for a lattice statistic, and not at all for any
# other, so whichever of the two appears first is where the frame ends.
frame_block <- function(out) {
  out[seq(2L, min(grep("Counts:|Bootstrap:", out)) - 1L)]
}

test_that("print.slope_bootstrap() prints a data frame, not a hand-drawn table", {
  expect_output(print(boot_bca), "<slope_bootstrap>")
  out <- capture.output(print(boot_bca))

  # Nothing left of the pipe table: the columns are R's own, named in a header
  # R spaces to its own widths, so the names are matched but not the gaps.
  expect_false(any(startsWith(out, "  |")))
  expect_true(any(grepl("statistic +calculated +mean +sd +ci$", out)))

  # Returns its argument invisibly, so `b <- print(b)` is safe and a bare
  # `print(b)` at the prompt does not print the object twice.
  capture.output(res <- withVisible(print(boot_bca)))
  expect_false(res$visible)
  expect_identical(res$value, boot_bca)
})

test_that("print.slope_bootstrap() names the method, replicate count and level once", {
  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  for (type in c("bca", "percentile")) {
    b <- suppressWarnings(slope_bootstrap(ss, R = 6, seed = 5, type = type))
    out <- capture.output(print(b))
    label <- if (identical(b$type, "bca")) "BCa" else "percentile"

    # A data frame has no span header, so the three facts that used to sit in
    # one -- method, replicate count, interval level -- lead the notes instead.
    expect_true(any(grepl(sprintf("Bootstrap:   R = 6 replicates; 95%% %s intervals.", label),
                          out, fixed = TRUE)))
    expect_false(any(grepl("Bootstrapped (", out, fixed = TRUE)))
    expect_false(any(grepl("Replicates:", out, fixed = TRUE)))

    # `type` on a single result is already the type actually used, so with the
    # slope gone there is nothing left for a footnote to except.
    expect_false(any(grepl("*", out, fixed = TRUE)))
  }
})

test_that("print.slope_bootstrap() shows a sample size on one basis at a time, per arm by default", {
  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  b <- suppressWarnings(slope_bootstrap(ss, R = 6, seed = 5))
  expect_identical(attr(b, "per_arm"), TRUE)

  out_arm <- capture.output(print(b))
  out_tot <- capture.output(print(b, per_arm = FALSE))
  expect_false(any(grepl("n_per_arm", out_arm, fixed = TRUE)))
  expect_false(any(grepl("n_per_arm", out_tot, fixed = TRUE)))

  per <- frame_row(out_arm, "n")
  tot <- frame_row(out_tot, "n")
  expect_length(tot, 8L)

  # The calculated column is the object's own two figures, not a re-derivation.
  expect_equal(as.numeric(tot[3L]), ss$n)
  expect_equal(as.numeric(per[3L]), ss$n_per_arm)
  expect_equal(2 * as.numeric(per[3L]), as.numeric(tot[3L]))

  # Both interval endpoints halve exactly: widen_to_lattice() has already moved
  # them out to even sizes, so a per-arm endpoint is a whole participant.
  ci_of <- function(cells) as.numeric(cells[c(6L, 8L)])
  expect_equal(ci_of(tot), b$ci)
  expect_equal(2 * ci_of(per), b$ci)
  expect_true(all(b$ci %% 2 == 0))

  # The mean and SD are the object's, halved for the per-arm row; the column is
  # rounded for printing, so the numbers are compared rather than the strings.
  expect_equal(as.numeric(tot[4L]), b$boot_mean, tolerance = 1e-6)
  expect_equal(as.numeric(tot[5L]), b$boot_sd, tolerance = 1e-6)
  expect_equal(2 * as.numeric(per[4L]), b$boot_mean, tolerance = 1e-6)
  expect_equal(2 * as.numeric(per[5L]), b$boot_sd, tolerance = 1e-6)

  # "Counts:" says which basis is on the page, and print(x, per_arm = ...)
  # overrides the object's own attribute without recomputing anything.
  expect_true(any(grepl("Counts:.*per arm", out_arm)))
  expect_true(any(grepl("Counts:.*total", out_tot)))
  expect_error(print(b, per_arm = "x"), "per_arm")
})

test_that("print.slope_bootstrap() gives a statistic with no arms a single row", {
  # The rows are units, and only a sample size comes in two of them. A slope, a
  # power or a target treatment effect has nothing to divide by.
  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  for (b in list(boot_bca, suppressWarnings(slope_bootstrap(ss, R = 6, statistic = "tte",
                                                            seed = 5)))) {
    out <- capture.output(print(b))
    row <- frame_row(out, b$statistic)
    expect_false(any(grepl("n_per_arm", out, fixed = TRUE)))
    expect_equal(as.numeric(row[3L]), b$observed, tolerance = 1e-6)
    expect_equal(as.numeric(row[4L]), b$boot_mean, tolerance = 1e-6)
    expect_equal(as.numeric(row[5L]), b$boot_sd, tolerance = 1e-6)
    expect_equal(as.numeric(row[c(6L, 8L)]), b$ci, tolerance = 1e-5)
  }
})

test_that("print.slope_bootstrap() no longer reports the resampled slope", {
  # The slope is still on the object -- and the straddle note below still
  # reports the one thing about it that bears on the interval -- but it is not
  # a row of the table any more.
  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  b <- suppressWarnings(slope_bootstrap(ss, R = 6, seed = 5))
  frame <- frame_block(capture.output(print(b)))

  expect_false(any(grepl("slope", frame, ignore.case = TRUE)))
  expect_false(any(grepl(format(b$slope_mean, digits = 6), frame, fixed = TRUE)))
  # Still there to be read off the object.
  expect_equal(b$slope_observed, ss$params$slope)
  expect_length(b$slope_ci, 2L)
})

test_that("the slope's summary is the replicates the straddle is measured on", {
  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  b <- suppressWarnings(slope_bootstrap(ss, R = 8, seed = 5))

  expect_length(b$slope_replicates, length(b$replicates))
  expect_equal(b$slope_observed, ss$params$slope)
  expect_equal(b$slope_mean, mean(b$slope_replicates))
  expect_equal(b$slope_sd, sd(b$slope_replicates))
  expect_equal(b$straddle,
               mean(sign(b$slope_replicates) != sign(b$slope_observed)))

  # For statistic = "slope" the two summaries are one computation, reused rather
  # than repeated -- a second pass would warn twice about a single failure.
  sb <- suppressWarnings(slope_bootstrap(ss$params, R = 8, seed = 5))
  expect_equal(sb$slope_ci, sb$ci)
  expect_equal(sb$slope_replicates, sb$replicates)
  expect_equal(sb$slope_type, sb$type)
})

test_that("the printed note describes what boot_mean is actually a mean of", {
  # The note claims each replicate is a whole trial, rounded up to a whole
  # participant per arm and doubled, and that the mean averages those already
  # rounded sizes. Each clause is checked here, so the wording cannot drift away
  # from solve_slope() without a test failing.
  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  b <- suppressWarnings(slope_bootstrap(ss, R = 20, seed = 5))

  # "rounds up to a whole participant per arm, the total being twice that":
  # every replicate is an even integer, which only 2 * ceiling(...) can be.
  expect_true(all(b$replicates == round(b$replicates)))
  expect_true(all(b$replicates %% 2 == 0))

  # "over the retained replicates" -- the failures are discarded, not imputed.
  expect_equal(length(b$replicates), b$R - b$n_failed)

  # "summarise sizes already rounded individually rather than rounding a
  # summary": the mean is of the rounded values, and is not itself rounded.
  expect_equal(b$boot_mean, mean(b$replicates))
  expect_equal(b$boot_sd, sd(b$replicates))

  # "not itself a size a trial could be run at": with a spread of replicates the
  # average lands off the lattice they live on.
  expect_gt(b$boot_sd, 0)
  expect_false(b$boot_mean %% 2 == 0)

  # And the slope, which the note exempts, carries no rounding.
  expect_false(all(b$slope_replicates == round(b$slope_replicates)))

  # The printed note is the short form; the full account lives in the help page,
  # so what is pinned here is that the printed one stays two lines and that each
  # of its two claims is the one checked above.
  # Default print is per arm, so the note's tail names the basis actually
  # shown -- "arm size" -- not the total the object's own fields carry above.
  out <- capture.output(print(b))
  note <- out[seq(grep("Mean, SD:", out, fixed = TRUE), length(out))]
  expect_length(note, 2L)
  # Matched within a line: the note is wrapped, so a phrase spanning the break
  # would be pinning where the wrap falls rather than what the note says.
  expect_true(any(grepl("rounded up to a whole participant per arm",
                        note, fixed = TRUE)))
  expect_true(any(grepl("the mean is not a runnable arm size", note, fixed = TRUE)))

  # Printed on the total basis instead, the same note names the total.
  note_tot <- capture.output(print(b, per_arm = FALSE))
  note_tot <- note_tot[seq(grep("Mean, SD:", note_tot, fixed = TRUE), length(note_tot))]
  expect_true(any(grepl("the mean is not a runnable trial size", note_tot, fixed = TRUE)))

  # Only the sample size rounds, so only it gets the note.
  expect_false(any(grepl("Mean, SD:", capture.output(print(boot_bca)), fixed = TRUE)))
})

test_that("print.slope_bootstrap() wraps its notes to a readable width", {
  # These are read in an 80-column terminal; the straddle note alone ran to 105
  # characters and wrapped raggedly.
  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  out <- capture.output(print(suppressWarnings(slope_bootstrap(ss, R = 6, seed = 5))))
  # Two notes now share the label; the wrapping check wants all of them, from
  # the first "Note:" line to the end.
  notes <- out[seq(min(grep("Note:", out, fixed = TRUE)), length(out))]
  expect_true(all(nchar(notes) <= 78L))
  # Continuation lines land in the value column, under the note's first word.
  continued <- notes[!grepl(":", notes, fixed = TRUE)]
  expect_gt(length(continued), 0L)
  expect_true(all(startsWith(continued, strrep(" ", 15L))))
  expect_false(any(startsWith(continued, strrep(" ", 16L))))
})

test_that("print.slope_bootstrap() always reports the straddle, zero included", {
  # The 0/R case is the point: a reader has to be able to tell "checked, and
  # nothing crossed zero" from a print method that never checked at all.
  expect_equal(boot_bca$straddle, 0)
  expect_true(any(grepl(sprintf("Note:        0/%d (0.0%%) of replicates refit a slope",
                                length(boot_bca$replicates)),
                        capture.output(print(boot_bca)), fixed = TRUE)))

  b <- suppressWarnings(slope_bootstrap(flat_fit(), R = 30, seed = 1))
  expect_gt(b$straddle, 0)
  out <- capture.output(print(b))
  # The count is printed beside the percentage because the percentage alone
  # hides how few replicates it can stand for.
  n_used <- length(b$replicates)
  expect_true(any(grepl(sprintf("%d/%d (%.1f%%) of replicates",
                                round(b$straddle * n_used), n_used, 100 * b$straddle),
                        out, fixed = TRUE)))
})

test_that("print.slope_bootstrap() always reports the convergence failures, zero included", {
  # Same rationale as the straddle note: "0 failed" and "never checked" must
  # not look the same on the page.
  expect_equal(boot_bca$n_failed, 0)
  expect_true(any(grepl(sprintf("Note:        0/%d (0.0%%) bootstrap samples failed to converge.",
                                boot_bca$R),
                        capture.output(print(boot_bca)), fixed = TRUE)))

  ss <- suppressWarnings(
    slope_sample_size(paper_fit("slpower1"), c(0, 1, 2), effectiveness = 0.33))
  b <- suppressWarnings(slope_bootstrap(ss, R = 6, seed = 5))
  b$n_failed <- 2L
  out <- capture.output(print(b))
  expect_true(any(grepl(sprintf("2/%d (%.1f%%) bootstrap samples failed to converge.",
                                b$R, 100 * 2 / b$R),
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
