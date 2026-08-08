# Layer 1 --- slope_params() and slope_params_manual(). See CONTRACT.md sec 2.

# Mixed-model fits are the slow part of this suite. Cache them in the global
# environment so that test-params.R and test-paper-parity.R together fit each
# dataset once. (Defined in both files rather than a shared helper because this
# agent owns only test-*.R.)
# --- slope_params_manual() --------------------------------------------------

test_that("slope_params_manual() returns exactly the contract fields", {
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                           sigma2_slope = 2, sigma_cov = 5,
                           sigma2_residual = 10)
  expect_s3_class(p, "slope_params")
  expect_setequal(names(p), c(
    "slope", "slope_comparator", "comparator", "sigma2_intercept",
    "sigma2_slope", "sigma_cov", "sigma2_residual", "n_obs", "n_subjects",
    "common_variance", "time_shifted", "fit", "call"))
  expect_length(names(p), 13L)
})

test_that("manually supplied parameters are stored unchanged", {
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                           sigma2_slope = 2, sigma_cov = 5,
                           sigma2_residual = 10)
  expect_equal(p$slope, -1.672)
  expect_equal(p$sigma2_intercept, 100)
  expect_equal(p$sigma2_slope, 2)
  expect_equal(p$sigma_cov, 5)
  expect_equal(p$sigma2_residual, 10)
  expect_identical(p$comparator, "none")
  expect_true(is.na(p$slope_comparator))
  expect_true(is.na(p$n_obs))
  expect_true(is.na(p$n_subjects))
  expect_null(p$fit)
})

test_that("slope_params_manual() rejects non-positive variance components", {
  base <- list(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
               sigma_cov = 5, sigma2_residual = 10)
  for (nm in c("sigma2_intercept", "sigma2_slope", "sigma2_residual")) {
    args <- base; args[[nm]] <- 0
    expect_error(do.call(slope_params_manual, args), nm)
    args[[nm]] <- -1
    expect_error(do.call(slope_params_manual, args), nm)
  }
})

test_that("slope_params_manual() rejects a non-positive-definite covariance", {
  # var_int * var_slope < cov^2
  expect_error(
    slope_params_manual(slope = -1, sigma2_intercept = 1, sigma2_slope = 1,
                        sigma_cov = 5, sigma2_residual = 10),
    "positive definite"
  )
})

test_that("a comparator slope is required unless comparator is 'none'", {
  for (cmp in c("healthy", "treated")) {
    expect_error(
      slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                          sigma_cov = 5, sigma2_residual = 10, comparator = cmp),
      "slope_comparator"
    )
  }
  expect_no_error(
    slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                        sigma_cov = 5, sigma2_residual = 10,
                        slope_comparator = 0.5, comparator = "healthy")
  )
})

test_that("a comparator slope is discarded when comparator is 'none'", {
  p <- slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                           sigma_cov = 5, sigma2_residual = 10,
                           slope_comparator = 0.5, comparator = "none")
  expect_true(is.na(p$slope_comparator))
})

test_that("slope_params_manual() rejects an unknown comparator", {
  expect_error(
    slope_params_manual(slope = -1, sigma2_intercept = 100, sigma2_slope = 2,
                        sigma_cov = 5, sigma2_residual = 10,
                        slope_comparator = 0.5, comparator = "nonsense")
  )
})

# --- slope_lme_control() ----------------------------------------------------

test_that("slope_lme_control() pins the settings behind every fit", {
  ctrl <- slope_lme_control()
  expect_equal(ctrl$maxIter, 200)
  expect_equal(ctrl$msMaxIter, 200)
  expect_equal(ctrl$niterEM, 50)
  expect_identical(ctrl$opt, "optim")
  expect_equal(ctrl$tolerance, 1e-7)
  expect_equal(ctrl$msTol, 1e-8)

  # These are not decoration. CONTRACT.md sec 7.1 calibrates the fitted-vs-Stata
  # tolerances to *these* numbers and says so explicitly: loosening them back
  # toward nlme's defaults widens the gap to `mixed`, and tightening them further
  # narrows it, so either change means retightening those tolerances too.
  defaults <- nlme::lmeControl()
  expect_gt(ctrl$maxIter, defaults$maxIter)
  expect_gt(ctrl$msMaxIter, defaults$msMaxIter)
  expect_gt(ctrl$niterEM, defaults$niterEM)
  expect_lt(ctrl$tolerance, defaults$tolerance)
  expect_lt(ctrl$msTol, defaults$msTol)
  expect_false(identical(ctrl$opt, defaults$opt))

  # returnObject = FALSE is what turns a failure to converge into an error
  # instead of a plausible-looking object at whatever point the optimiser
  # stopped. The healthy-model fallback tested below exists only because of it.
  expect_false(ctrl$returnObject)
})

test_that("slope_lme_control() is the control slope_params() actually uses", {
  # The help page promises the settings behind every fit are "inspectable and
  # reproducible outside the package". That is only true if refitting the model
  # by hand with the exported control lands on the packaged fit exactly, so
  # reproduce slpower1's single-group fit -- internal column names and all --
  # and demand agreement to the last bit rather than to a tolerance.
  d <- load_paper_data("slpower1")
  dat <- data.frame(sp_y = d$sdmt, sp_time = d$visit,
                    sp_subject = factor(as.character(d$id)))
  by_hand <- nlme::lme(sp_y ~ sp_time, random = ~ sp_time | sp_subject,
                       data = dat, method = "REML",
                       control = slope_lme_control())
  ref <- paper_fit("slpower1")
  expect_equal(unname(nlme::fixef(by_hand)[["sp_time"]]), ref$slope,
               tolerance = 1e-12)
  expect_equal(stats::sigma(by_hand)^2, ref$sigma2_residual, tolerance = 1e-12)
  G <- nlme::getVarCov(by_hand)
  expect_equal(as.numeric(G["(Intercept)", "(Intercept)"]), ref$sigma2_intercept,
               tolerance = 1e-12)
  expect_equal(as.numeric(G["sp_time", "sp_time"]), ref$sigma2_slope,
               tolerance = 1e-12)
})

test_that("slope_params() has no control argument to override the settings", {
  # Deliberate, and documented as such in ?slope_lme_control: the model this
  # package fits is fixed by the method, so the control object exists to be
  # inspected rather than replaced. A silently ignored `control =` would be the
  # worst of both worlds.
  d <- load_paper_data("slpower1")
  expect_error(slope_params(sdmt ~ visit | id, d, control = slope_lme_control()),
               "unused argument")
})

# --- slope_params(): scenario dispatch --------------------------------------

test_that("the scenario follows which group argument is supplied", {
  expect_identical(paper_fit("slpower1")$comparator, "none")
  expect_identical(paper_fit("slpower2")$comparator, "healthy")
  expect_identical(paper_fit("slpower3")$comparator, "treated")
})

test_that("healthy and treated are mutually exclusive", {
  d <- load_paper_data("slpower2")
  d$treat <- d$case
  expect_error(slope_params(sdmt ~ time | id, d, healthy = case, treated = treat),
               "only one of")
})

test_that("the formula must name a subject identifier", {
  d <- load_paper_data("slpower1")
  expect_error(slope_params(sdmt ~ time, d), "subject identifier")
})

test_that("a group variable must be binary and 0/1 coded", {
  d <- load_paper_data("slpower2")
  d$bad <- d$case + 1                      # coded 1/2
  expect_error(slope_params(sdmt ~ time | id, d, healthy = bad), "0/1")
  d$constant <- 1L
  expect_error(slope_params(sdmt ~ time | id, d, healthy = constant),
               "two distinct values")
})

test_that("slope_params() reports the number of observations and subjects", {
  expect_equal(paper_fit("slpower1")$n_obs, 800L)
  expect_equal(paper_fit("slpower1")$n_subjects, 200L)
  expect_equal(paper_fit("slpower2")$n_obs, 2000L)
  expect_equal(paper_fit("slpower2")$n_subjects, 500L)
  expect_equal(paper_fit("slpower3")$n_obs, 450L)
  expect_equal(paper_fit("slpower3")$n_subjects, 150L)
})

test_that("fitted variance components are positive and positive definite", {
  for (nm in c("slpower1", "slpower2", "slpower3")) {
    p <- paper_fit(nm)
    expect_gt(p$sigma2_intercept, 0)
    expect_gt(p$sigma2_slope, 0)
    expect_gt(p$sigma2_residual, 0)
    G <- matrix(c(p$sigma2_intercept, p$sigma_cov,
                  p$sigma_cov, p$sigma2_slope), 2, 2)
    expect_gt(det(G), 0)
  }
})

test_that("slpower1 variance components land near the simulation truth", {
  # Paper appendix: sigma2_a = 100, sigma2_b = 2, sigma_ab = 5, sigma2_e = 10.
  # Sampling noise over 200 subjects, so these are loose bounds.
  p <- paper_fit("slpower1")
  expect_gt(p$sigma2_intercept, 70);  expect_lt(p$sigma2_intercept, 140)
  expect_gt(p$sigma2_slope, 1);       expect_lt(p$sigma2_slope, 4)
  expect_gt(p$sigma2_residual, 6);    expect_lt(p$sigma2_residual, 14)
})

# --- input coercion ----------------------------------------------------------
# What `formula` and the group arguments will accept, and what they refuse. The
# level-naming trap for factor and character group columns is a separate story
# and lives in test-review-regressions.R; these are the surrounding branches.

test_that("coerce_time() unwraps Date and POSIXct onto an axis of days", {
  # Date support is on the help page, so it is a promise, not an accident:
  # as.numeric() of a Date is days since 1970-01-01.
  expect_equal(slopepower:::coerce_time(as.Date("1970-01-11"), "ctx"), 10)
  # POSIXct counts seconds, so it is divided by 86400 -- without which the two
  # classes would disagree by a factor of 86400 on the same instant, and a
  # `visits` schedule written for one would be nonsense for the other.
  expect_equal(
    slopepower:::coerce_time(as.POSIXct("1970-01-11 12:00:00", tz = "UTC"), "ctx"),
    10.5)

  # A labelled time variable is unwrapped to its codes and then taken as
  # numeric, which is what makes `visit` from a Stata .dta usable as it stands.
  # The class is set by hand because `haven` is only in Suggests.
  expect_equal(
    slopepower:::coerce_time(structure(c(0, 1, 2), class = "haven_labelled",
                                       labels = c(baseline = 0)), "ctx"),
    c(0, 1, 2))

  # Anything else is refused rather than coerced by as.numeric(): a factor of
  # visit numbers would otherwise silently become its integer codes, which are
  # the level *ranks* and not the times.
  expect_error(slopepower:::coerce_time(c("0", "1"), "ctx"),
               "the time variable must be numeric \\(or a Date\\)")
  expect_error(slopepower:::coerce_time(factor(c(0, 1)), "ctx"),
               "the time variable must be numeric \\(or a Date\\)")
})

test_that("a Date or POSIXct time column fits, on an axis of days", {
  d <- load_paper_data("slpower1")
  yearly <- paper_fit("slpower1")

  # slpower1's annual visits, rewritten as calendar dates. Fitting on days
  # rather than years is the badly scaled axis ?slope_params warns about, so
  # compare the slope after converting back rather than expecting the variance
  # components to match to many digits.
  d$vdate <- as.Date("2010-01-01") + d$visit * 365
  d$vtime <- as.POSIXct(paste(d$vdate, "18:00:00"), tz = "UTC")
  as_date  <- suppressMessages(slope_params(sdmt ~ vdate | id, d))
  as_posix <- suppressMessages(slope_params(sdmt ~ vtime | id, d))

  expect_equal(as_date$slope * 365, yearly$slope, tolerance = 1e-5)
  expect_true(as_date$time_shifted)
  # The two land on the same axis: /86400 puts POSIXct in days, and the constant
  # time of day is removed again by the per-subject re-origining.
  expect_equal(as_posix$slope, as_date$slope, tolerance = 1e-12)
  expect_equal(as_posix$sigma2_residual, as_date$sigma2_residual, tolerance = 1e-12)

  d$visit_chr <- as.character(d$visit)
  expect_error(slope_params(sdmt ~ visit_chr | id, d),
               "the time variable must be numeric \\(or a Date\\)")
})

test_that("a logical group column is taken as coded, not guessed at", {
  # TRUE/FALSE is the natural thing to hand `treated =` or `healthy =`, and
  # unlike the c("case", "control") factor of test-review-regressions.R there is
  # nothing to determine: TRUE is 1. So it needs no recoding by the user, and
  # must reproduce the numeric fit exactly rather than approximately.
  d <- load_paper_data("slpower3")
  numeric_fit <- paper_fit("slpower3")
  d$trt <- d$treat == 1
  logical_fit <- slope_params(sdmt ~ time | id, d, treated = trt)
  expect_equal(logical_fit$slope, numeric_fit$slope, tolerance = 1e-12)
  expect_equal(logical_fit$slope_comparator, numeric_fit$slope_comparator,
               tolerance = 1e-12)
  expect_equal(logical_fit$sigma2_slope, numeric_fit$sigma2_slope, tolerance = 1e-12)

  # And the direction is real: negating the column swaps which arm is which, so
  # the control-arm slope becomes the experimental one and back again. If TRUE
  # were mapped to 0 this would be the identity instead.
  d$trt_flipped <- d$treat == 0
  flipped <- slope_params(sdmt ~ time | id, d, treated = trt_flipped)
  expect_equal(flipped$slope, numeric_fit$slope_comparator, tolerance = 1e-12)
  expect_equal(flipped$slope_comparator, numeric_fit$slope, tolerance = 1e-12)
})

test_that("a haven_labelled group column is unwrapped to its numeric codes", {
  # Stata data read with haven arrives labelled, which is the whole reason this
  # branch exists -- the paper's own datasets are .dta files. `haven` is only in
  # Suggests, so the class is set by hand here: coerce_binary() dispatches on
  # inherits(), and which branch it takes must not depend on haven being
  # installed. The labels are what the numeric path would otherwise reject.
  lab <- structure(c(0, 1, 0, 1), class = "haven_labelled",
                   labels = c(control = 0, case = 1))
  expect_equal(slopepower:::coerce_binary(lab, "healthy", "slope_params()", "case"),
               c(0, 1, 0, 1))
})

test_that("a group column with more than two levels is rejected on the count", {
  # Three levels is not the ambiguity that test-review-regressions.R covers:
  # there, two levels with unusable names could at least be recoded 0/1 by the
  # user, and the error says how. Here there is no two-arm question to ask at
  # all, so the count is checked first and the message is a different one.
  expect_error(
    slopepower:::coerce_binary(factor(c("a", "b", "c")), "healthy",
                               "slope_params()", "case"),
    "`healthy` must have exactly two levels; got 3")
  expect_error(
    slopepower:::coerce_binary(c("a", "b", "c"), "healthy",
                               "slope_params()", "case"),
    "`healthy` must have exactly two levels; got 3")
  # The numeric path counts distinct values instead of levels and says so, which
  # is the difference between "this factor has an unused level" and "these data
  # have three groups in them".
  expect_error(
    slopepower:::coerce_binary(c(0, 1, 2), "treated", "slope_params()", "treated"),
    "`treated` must have exactly two distinct values; got 3")
})

test_that("a group column of an unsupported type is rejected on the type", {
  # A Date is a numeric vector underneath, but is.numeric() is FALSE for it, so
  # it reaches the type error rather than being read as day numbers and then
  # complained about for having the wrong values.
  expect_error(
    slopepower:::coerce_binary(as.Date("2020-01-01") + 0:1, "healthy",
                               "slope_params()", "case"),
    "`healthy` must be numeric, logical, factor or labelled")
  expect_error(
    slopepower:::coerce_binary(list(0, 1), "treated", "slope_params()", "treated"),
    "`treated` must be numeric, logical, factor or labelled")
})

test_that("a group argument that evaluates to NULL is reported, not ignored", {
  # `healthy = g` where g is NULL in the caller's environment. substitute() sees
  # the symbol and so selects the healthy scenario; the NULL only surfaces when
  # the column is evaluated. Quietly falling back to comparator = "none" there
  # would answer a different question from the one asked -- powering toward a
  # slope of zero rather than toward the controls' slope -- and would do it
  # without saying so.
  d <- load_paper_data("slpower1")
  g <- NULL
  expect_error(slope_params(sdmt ~ visit | id, d, healthy = g),
               "`healthy` evaluated to NULL")
})

# --- per-subject time re-origining ------------------------------------------

test_that("time is shifted so each subject starts at zero", {
  # slpower2 records calendar visit dates, so every subject starts elsewhere.
  expect_true(paper_fit("slpower2")$time_shifted)
  # slpower1 and slpower3 already start at 0.
  expect_false(paper_fit("slpower1")$time_shifted)
  expect_false(paper_fit("slpower3")$time_shifted)
})

test_that("the re-origining emits a message, as Stata warns", {
  d <- load_paper_data("slpower2")
  expect_message(slope_params(sdmt ~ time | id, d, healthy = case),
                 "shifted")
})

test_that('origin = "none" leaves time alone and changes the fit', {
  d <- load_paper_data("slpower2")
  shifted   <- paper_fit("slpower2")
  unshifted <- suppressMessages(
    slope_params(sdmt ~ time | id, d, healthy = case, origin = "none"))
  expect_false(unshifted$time_shifted)
  # the slope is barely affected, but the intercept variance is not comparable
  expect_equal(unshifted$slope, shifted$slope, tolerance = 0.05)
  expect_false(isTRUE(all.equal(unshifted$sigma2_intercept,
                                shifted$sigma2_intercept)))
})

test_that("an explicit unit change rescales the slope proportionally", {
  d <- load_paper_data("slpower1")
  yearly <- paper_fit("slpower1")
  d$time <- d$time / 0.5                  # half-year axis
  halved <- slope_params(sdmt ~ time | id, d)
  expect_equal(halved$slope, yearly$slope / 2, tolerance = 1e-6)
})

# --- common_variance --------------------------------------------------------

test_that("common_variance is ignored with a warning outside the healthy case", {
  d <- load_paper_data("slpower1")
  expect_warning(slope_params(sdmt ~ time | id, d, common_variance = TRUE),
                 "only when `healthy`")
})

test_that("common_variance does not change the case estimates", {
  # The healthy model factorises into two independent per-group fits, so
  # reducing the controls' random-effects block cannot move the cases' numbers.
  d <- load_paper_data("slpower2")
  full    <- paper_fit("slpower2")
  reduced <- suppressMessages(
    slope_params(sdmt ~ time | id, d, healthy = case, common_variance = TRUE))
  expect_true(reduced$common_variance)
  expect_false(full$common_variance)
  expect_equal(reduced$slope, full$slope, tolerance = 1e-5)
  expect_equal(reduced$sigma2_intercept, full$sigma2_intercept, tolerance = 1e-4)
  expect_equal(reduced$sigma2_slope, full$sigma2_slope, tolerance = 1e-4)
  expect_equal(reduced$sigma2_residual, full$sigma2_residual, tolerance = 1e-4)
})

# --- the healthy-model convergence fallback ---------------------------------
# In Stata the reduced random-effects structure for healthy controls is
# something the user asks for, with `nocontvar`. Here it is that *and* what the
# port retreats to on its own when the full two-block structure will not
# converge. `common_variance` chooses between the two behaviours: TRUE demands
# the reduced structure, FALSE forbids it, NULL (the default) falls back and
# says so. CONTRACT.md sec 2 records why the retreat is safe -- the healthy
# model factorises into two independent per-group fits, so dropping the
# controls' random slope cannot move the cases' estimates, which are the only
# variance components the object carries.

# Data on which the *full* structure genuinely has no maximum to find, while
# the reduced one does. Nothing here is random: the failure is a property of the
# data rather than of which seed happened to be drawn.
#
# The healthy controls' trajectories are exactly linear -- no within-subject
# noise at all -- so their residual variance is identically zero. `nlme` carries
# it as a log-scale ratio against the cases' residual (the varIdent structure),
# and zero has no representation on that scale, so the optimiser runs off toward
# -Inf and `returnObject = FALSE` in slope_lme_control() turns that into an
# error. The reduced structure has no control-specific random slope, so the
# between-control variation in slope that it cannot represent lands in a
# strictly positive residual instead, and it converges.
#
# A constant `control_slopes` makes the controls exactly parallel as well as
# exactly linear. Then the reduced structure has nothing to absorb either, its
# own residual is zero too, and both fits fail.
#
# The cases are ordinary: 30 subjects at four visits, with an intercept, a
# slope and a within-subject wobble, all written as deterministic functions of
# the subject index so the whole dataset is fixed.
noiseless_controls <- function(control_slopes) {
  visits <- 0:3
  cases <- do.call(rbind, lapply(seq_len(30), function(i) data.frame(
    id    = sprintf("case%02d", i),
    case  = 1,
    visit = visits,
    sdmt  = (50 + 12 * cos(1.9 * i)) +          # subject intercept
            (-1.7 + 1.2 * sin(3.3 * i)) * visits +  # subject slope
            3 * sin(2.4 * i + 1.7 * visits)         # within-subject residual
  )))
  intercepts <- seq(35, 65, length.out = length(control_slopes))
  controls <- do.call(rbind, lapply(seq_along(control_slopes), function(i)
    data.frame(id    = sprintf("ctrl%03d", i),
               case  = 0,
               visit = 0:2,
               sdmt  = intercepts[i] + control_slopes[i] * (0:2))))
  rbind(cases, controls)
}

# The spread of control slopes is what the reduced structure absorbs into its
# residual; a single value leaves it nothing to absorb.
varying_control_slopes <- seq(-3, 2, length.out = 100)
parallel_control_slopes <- rep(-0.5, 100)

test_that("the healthy fit falls back to a reduced structure and says so", {
  d <- noiseless_controls(varying_control_slopes)
  expect_message(p <- slope_params(sdmt ~ visit | id, d, healthy = case),
                 "falling back to a random intercept only for controls")
  expect_true(p$common_variance)

  # The fallback must land on exactly the fit `common_variance = TRUE` asks for
  # outright -- the same model, reached by a different route -- not merely on a
  # similar one.
  forced <- suppressMessages(
    slope_params(sdmt ~ visit | id, d, healthy = case, common_variance = TRUE))
  expect_true(forced$common_variance)
  for (nm in c("slope", "slope_comparator", "sigma2_intercept", "sigma2_slope",
               "sigma_cov", "sigma2_residual")) {
    expect_equal(p[[nm]], forced[[nm]], tolerance = 1e-12)
  }
})

test_that("common_variance = FALSE forbids the fallback and errors instead", {
  # Supplying FALSE is a statement that the reduced structure would answer a
  # different question, so the retreat is not the port's to make: it reports
  # what happened and stops. The nlme failure is carried through rather than
  # swallowed, because "did not converge" alone leaves nothing to act on.
  d <- noiseless_controls(varying_control_slopes)
  msgs <- character()
  err <- withCallingHandlers(
    tryCatch(slope_params(sdmt ~ visit | id, d, healthy = case,
                          common_variance = FALSE),
             error = conditionMessage),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_match(err, "the full model did not converge")
  expect_match(err, "`common_variance = FALSE` forbids the reduced structure")
  expect_match(err, "Underlying error: .")
  # and the fallback message is not emitted on the way past
  expect_length(msgs, 0L)
})

test_that("a reduced structure that also fails is an error, not a third try", {
  d <- noiseless_controls(parallel_control_slopes)
  # Under the default the fallback is announced first and then fails: the
  # message is not a promise that the refit worked.
  expect_error(suppressMessages(slope_params(sdmt ~ visit | id, d, healthy = case)),
               "the mixed model did not converge")

  # The same error is reached when the reduced structure was asked for outright,
  # in which case the full model was never attempted and there is nothing to
  # announce.
  msgs <- character()
  err <- withCallingHandlers(
    tryCatch(slope_params(sdmt ~ visit | id, d, healthy = case,
                          common_variance = TRUE),
             error = conditionMessage),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_match(err, "the mixed model did not converge")
  expect_length(msgs, 0L)
})

# --- data adequacy regressions -----------------------------------------------
# Code-review regressions: too little data to identify the random-slope model
# used to fit silently and return an unstable estimate rather than fail.

test_that("a subject seen only once does not count toward the repeat-visit floor", {
  # 2 subjects, 3 rows, one subject with a single visit: fewer than 3 usable
  # rows would already be rejected, but 3 rows is enough to slip past that
  # guard while still being unidentifiable (only 1 of 2 subjects contributes
  # any information about the random-slope variance).
  d <- data.frame(id = c(1, 1, 2), visit = c(0, 1, 0), sdmt = c(40, 38, 35))
  expect_error(slope_params(sdmt ~ visit | id, d), "repeat visits")
})

test_that("a comparator group of 1 subject is rejected, not silently fit", {
  set.seed(2)
  subj <- data.frame(id = 1:41, case = c(1, rep(0, 40)))
  subj$intercept <- rnorm(41, 50, 10)
  subj$slope <- rnorm(41, ifelse(subj$case == 1, -1.7, -0.3), 0.5)
  sim <- merge(subj, data.frame(visit = 0:3))
  sim$sdmt <- sim$intercept + sim$slope * sim$visit + rnorm(nrow(sim), 0, 3)
  expect_error(
    suppressMessages(slope_params(sdmt ~ visit | id, sim, healthy = case)),
    "at least 2 participants")
})

test_that("small or unbalanced comparator groups warn but still fit", {
  set.seed(2)
  subj <- data.frame(id = 1:22, case = c(rep(1, 2), rep(0, 20)))
  subj$intercept <- rnorm(22, 50, 10)
  subj$slope <- rnorm(22, ifelse(subj$case == 1, -1.7, -0.3), 1.4)
  sim <- merge(subj, data.frame(visit = 0:3))
  sim$sdmt <- sim$intercept + sim$slope * sim$visit + rnorm(nrow(sim), 0, 3)
  expect_warning(
    suppressMessages(slope_params(sdmt ~ visit | id, sim, healthy = case)),
    "small or unbalanced")
})

# --- printing ---------------------------------------------------------------

test_that("print.slope_params() runs for each scenario", {
  expect_output(print(paper_fit("slpower1")), "single group")
  expect_output(print(paper_fit("slpower2")), "healthy controls")
  expect_output(print(paper_fit("slpower3")), "previous randomised trial")
  out <- capture.output(print(paper_fit("slpower2")))
  expect_true(any(grepl("observed difference in slopes", out)))
})

test_that("print.slope_params() flags a reduced random-effects structure", {
  # Nothing else in the printout distinguishes a reduced fit from a full one --
  # the variance components shown are the cases' either way -- so the note is
  # the only place a reader learns that the controls' random slope was dropped.
  # It names Stata's `nocontvar` because that is the option a Stata user would
  # have had to set by hand to get here, and it says the case estimates are
  # unaffected because CONTRACT.md sec 2 is why the retreat was allowed at all.
  d <- noiseless_controls(varying_control_slopes)
  p <- suppressMessages(
    slope_params(sdmt ~ visit | id, d, healthy = case, common_variance = TRUE))
  out <- capture.output(print(p))
  expect_true(any(grepl("reduced random-effects structure used for healthy controls",
                        out, fixed = TRUE)))
  expect_true(any(grepl("Stata `nocontvar`", out, fixed = TRUE)))
  expect_true(any(grepl("Case estimates are unaffected", out, fixed = TRUE)))

  # And absent when the full structure was used, so its presence means something.
  expect_false(any(grepl("nocontvar", capture.output(print(paper_fit("slpower2"))),
                         fixed = TRUE)))
})

test_that("print.slope_params() marks manually supplied parameters", {
  p <- slope_params_manual(slope = -1.672, sigma2_intercept = 100,
                           sigma2_slope = 2, sigma_cov = 5,
                           sigma2_residual = 10)
  expect_output(print(p), "supplied directly")
})
