# Residual structures: `correlation` and `weights` in slope_params() and
# slope_params_manual(), the `residual` field, and its use in stage two.
# See CONTRACT.md sections 2, 5.1 and 5.7.
#
# The anchor for everything here is nlme itself: the covariance stage two builds
# at a participant's own visit times must be the marginal covariance nlme
# reports for that participant, for every structure. The correlation functions
# in residual.R restate nlme's definitions, and these tests are what stop the
# two from drifting apart.

library(nlme)

# Fits on slpower1 take well under a second each; cache them anyway, since
# several tests share each one.
residual_fit <- function(name) {
  cache <- sp_fit_cache()
  key <- paste0("residual_", name)
  if (!is.null(cache[[key]])) return(cache[[key]])
  d <- load_paper_data("slpower1")
  p <- switch(name,
    ar1   = slope_params(sdmt ~ time | id, d, correlation = corAR1()),
    car1  = slope_params(sdmt ~ time | id, d, correlation = corCAR1()),
    exp   = slope_params(sdmt ~ time | id, d, correlation = corExp(nugget = TRUE)),
    gaus  = slope_params(sdmt ~ time | id, d, correlation = corGaus()),
    symm  = slope_params(sdmt ~ time | id, d, correlation = corSymm()),
    vi    = slope_params(sdmt ~ time | id, d, weights = varIdent(form = ~ 1 | time)),
    un    = slope_params(sdmt ~ time | id, d, correlation = corSymm(),
                         weights = varIdent(form = ~ 1 | time)))
  cache[[key]] <- p
  p
}

# The marginal covariance nlme reports for one participant, against
# slope_sigma() at that participant's own (fitted) visit times.
expect_matches_nlme <- function(p, id) {
  g <- getData(p$fit)
  t <- g$sp_time[g$sp_subject == id]
  V <- unclass(getVarCov(p$fit, individual = id, type = "marginal")[[1L]])
  S <- unclass(slope_sigma(p, t))
  attributes(V) <- attributes(S) <- list(dim = dim(S))
  expect_equal(S, V, tolerance = 1e-8)
}

# --- agreement with nlme ------------------------------------------------------

test_that("slope_sigma() is nlme's marginal covariance for every structure", {
  for (nm in c("ar1", "car1", "exp", "gaus", "symm", "vi", "un")) {
    p <- residual_fit(nm)
    expect_false(is.null(p$residual), info = nm)
    expect_matches_nlme(p, "1")
  }
})

test_that("the structure is matched under `healthy`, with per-group variances kept", {
  d <- load_paper_data("slpower2")
  p <- suppressMessages(slope_params(sdmt ~ time | id, d, comparator = "healthy", group = case,
                                     correlation = corCAR1()))
  expect_matches_nlme(p, as.character(d$id[d$case == 1][1L]))

  # A shared visit grid, so visit-specific variances can be fitted too.
  d$visit <- stats::ave(d$time, d$id, FUN = seq_along) - 1
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, comparator = "healthy", group = case,
                                     correlation = corCAR1(),
                                     weights = varIdent(form = ~ 1 | visit)))
  expect_matches_nlme(p, as.character(d$id[d$case == 1][1L]))
  expect_equal(p$residual$times, 0:3)
  expect_equal(p$residual$sd_ratio[1L], 1)
})

test_that("the structure is matched under `treated`, on unequally spaced visits", {
  d <- load_paper_data("slpower3")
  p <- slope_params(sdmt ~ visit | id, d, comparator = "treated", group = treat, correlation = corSymm(),
                    weights = varIdent(form = ~ 1 | visit))
  expect_matches_nlme(p, "1")
  expect_equal(p$residual$times, c(0, 0.5, 2))
  expect_length(p$residual$coef, 3L)
})

test_that("corAR1() and corCAR1() agree on integer times", {
  a <- residual_fit("ar1")
  b <- residual_fit("car1")
  expect_equal(a$residual$coef[["Phi"]], b$residual$coef[["Phi"]], tolerance = 1e-2)
  expect_equal(slope_sample_size(a, 0:2, effectiveness = 0.33)$n,
               slope_sample_size(b, 0:2, effectiveness = 0.33)$n)
})

test_that("leaving `correlation` and `weights` unset is the original model exactly", {
  p <- paper_fit("slpower1")
  expect_null(p$residual)
  expect_null(p$fit$modelStruct$corStruct)
  expect_identical(slope_sigma(p, 0:3)[2, 2],
                   p$sigma2_intercept + p$sigma2_slope + 2 * p$sigma_cov + p$sigma2_residual)
})

test_that("a starting value and fixed = TRUE are honoured as in nlme", {
  d <- load_paper_data("slpower1")
  p <- slope_params(sdmt ~ time | id, d, correlation = corCAR1(0.4, fixed = TRUE))
  expect_equal(p$residual$coef[["Phi"]], 0.4)
  expect_true(p$residual$fixed)
  expect_match(paste(deparse(p$fit$call), collapse = ""), "fixed = TRUE", fixed = TRUE)
})

test_that("an explicit `form` naming the formula's own time and subject is accepted", {
  d <- load_paper_data("slpower1")
  p <- slope_params(sdmt ~ time | id, d, correlation = corCAR1(form = ~ time | id))
  expect_equal(p$residual$coef, residual_fit("car1")$residual$coef)
})

# --- refusals in slope_params() ----------------------------------------------

test_that("slope_params() refuses structures it cannot use, saying why", {
  d <- load_paper_data("slpower1")
  fit <- function(...) slope_params(sdmt ~ time | id, d, ...)
  expect_error(fit(correlation = corCompSymm()), "random intercept already models")
  expect_error(fit(correlation = corARMA(p = 2)), "not supported")
  expect_error(fit(correlation = "ar1"), "nlme correlation structure")
  expect_error(fit(correlation = corAR1(form = ~ time | site)), "must be omitted or be")
  expect_error(fit(correlation = corAR1(form = ~ visit | id)), "must be omitted or be")
  expect_error(fit(weights = varIdent()), "with no stratum is a constant variance")
  expect_error(fit(weights = varIdent(form = ~ 1 | id)), "varIdent\\(form = ~ 1 \\| time\\)")
  expect_error(fit(weights = varPower()), "must be an nlme `varIdent\\(\\)`")
  expect_error(fit(weights = varIdent(c("1" = 2), form = ~ 1 | time)), "initial or fixed")
})

test_that("corAR1() is refused on non-integer times, pointing to corCAR1()", {
  d <- load_paper_data("slpower3")
  expect_error(slope_params(sdmt ~ visit | id, d, comparator = "treated", group = treat, correlation = corAR1()),
               "corCAR1")
})

test_that("corSymm() and varIdent() need a shared visit schedule", {
  d <- load_paper_data("slpower2")   # visits recorded as dates
  expect_error(suppressMessages(slope_params(sdmt ~ time | id, d, comparator = "healthy", group = case,
                                             correlation = corSymm())),
               "shared across participants")
  expect_error(suppressMessages(slope_params(sdmt ~ time | id, d, comparator = "healthy", group = case,
                                             weights = varIdent(form = ~ 1 | time))),
               "shared across participants")
})

test_that("two measurements at one time within a participant are refused", {
  d <- load_paper_data("slpower1")
  d <- rbind(d, d[d$id == 1 & d$time == 1, ])
  expect_error(slope_params(sdmt ~ time | id, d, correlation = corCAR1()),
               "at most one measurement per participant")
})

# --- stage two -----------------------------------------------------------------

test_that("a per-visit structure prices schedules among its times only", {
  p <- residual_fit("un")
  expect_silent(slope_sample_size(p, c(0, 1, 3), effectiveness = 0.33))
  expect_error(slope_sample_size(p, c(0, 0.5, 2), effectiveness = 0.33),
               "says nothing about time\\(s\\) 0.5")
  expect_error(slope_var(residual_fit("vi"), c(0, 1, 4)), "time\\(s\\) 4")
})

test_that("corAR1() refuses a planned schedule off the integers", {
  expect_error(slope_var(residual_fit("ar1"), c(0, 0.5, 1)), "integer times only")
  expect_silent(slope_var(residual_fit("car1"), c(0, 0.5, 1)))
})

test_that("dropout strata slice the structured covariance correctly", {
  # Each stratum's s*^2 is computed on a leading submatrix of the full
  # schedule's covariance. Rebuild it from slope_var() on each prefix instead:
  # with a gap in the schedule (no visit 2), a corSymm() lookup that indexed by
  # position rather than by time would disagree.
  p <- residual_fit("un")
  visits <- c(0, 1, 3)
  drop <- c(0, 0.2)
  es <- slope_effect_size(p, visits, dropout = drop)
  d <- p$slope
  manual <- sqrt((1 - sum(drop)) * d^2 / slope_var(p, visits) +
                 drop[2] * d^2 / slope_var(p, visits[1:2]))
  expect_equal(abs(es), manual)
})

test_that("every stage-two print carries the residual note", {
  out <- capture.output(slope_sample_size(residual_fit("car1"), 0:2, effectiveness = 0.33))
  expect_true(any(grepl("continuous-time AR(1)", out, fixed = TRUE)))
  out <- capture.output(slope_power(residual_fit("un"), 0:2, n = 400, effectiveness = 0.33))
  out <- gsub("\\s+", " ", paste(out, collapse = " "))
  expect_match(out, "planned visits must be among those times", fixed = TRUE)
})

# --- the floor -----------------------------------------------------------------

test_that("the floor still bounds every schedule under a serial correlation", {
  for (nm in c("car1", "exp")) {
    p <- residual_fit(nm)
    fl <- slope_var_floor(p)
    set.seed(11)
    for (i in 1:20) {
      v <- c(0, sort(stats::runif(sample(2:8, 1L), 0.1, 20)))
      expect_gt(slope_var(p, v), fl)
    }
    expect_lt(slope_var(p, seq(0, 200, by = 0.5)) / fl, 1.01)
  }
})

test_that("the floor is refused for a per-visit structure", {
  expect_error(slope_var_floor(residual_fit("un")), "no floor over all schedules")
  expect_error(slope_var_floor(residual_fit("vi")), "slope_var\\(params, c\\(0, 1, 2, 3\\)\\)")
  expect_error(slope_sample_size_floor(residual_fit("symm"), effectiveness = 0.33),
               "not separately identified")
})

# --- slope_params_manual() -----------------------------------------------------

manual <- function(...) {
  slope_params_manual(slope = -1.672, sigma2_intercept = 100, sigma2_slope = 2,
                      sigma_cov = 5, sigma2_residual = 10, ...)
}

test_that("a stated serial correlation is used as written", {
  p <- manual(correlation = corCAR1(0.5))
  expect_equal(p$residual$coef, c(Phi = 0.5))
  t <- c(0, 1, 2.5)
  expect_equal(unname(slope_sigma(p, t) - slope_sigma(ref_params(), t)),
               10 * (0.5^abs(outer(t, t, "-")) - diag(3)))

  p <- manual(correlation = corExp(c(2, 0.1), nugget = TRUE))
  expect_equal(p$residual$coef, c(range = 2, nugget = 0.1))
  expect_equal(unname(slope_sigma(p, c(0, 1))[1, 2] - slope_sigma(ref_params(), c(0, 1))[1, 2]),
               10 * 0.9 * exp(-1 / 2))
})

test_that("a stated unstructured covariance is keyed by time", {
  p <- manual(correlation = corSymm(c(0.3, 0.2, 0.4)),
              weights = varIdent(c("1" = 1.2, "2.0" = 1.5), form = ~ 1 | visit),
              times = c(0, 1, 2))
  expect_equal(p$residual$sd_ratio, c(1, 1.2, 1.5))
  R <- slope_sigma(p, c(0, 2)) - slope_sigma(ref_params(), c(0, 2)) + diag(10, 2)
  expect_equal(unname(R), 10 * matrix(c(1, 0.2 * 1.5, 0.2 * 1.5, 1.5^2), 2))
})

test_that("independent residuals stated as a structure change nothing", {
  zero <- manual(correlation = corSymm(c(0, 0, 0)), times = 0:2)
  expect_equal(slope_var(zero, 0:2), slope_var(ref_params(), 0:2))
  expect_equal(slope_var(manual(correlation = corCAR1(1e-12)), c(0, 1, 3)),
               slope_var(ref_params(), c(0, 1, 3)))
})

test_that("slope_params_manual() refuses incomplete or invalid structures", {
  expect_error(manual(correlation = corExp()), "state the value")
  expect_error(manual(correlation = corSymm(c(0.3, 0.2, 0.4))), "`times` is required")
  expect_error(manual(correlation = corCAR1(0.5), times = 0:2), "function of time")
  expect_error(manual(times = 0:2), "applies only with")
  expect_error(manual(correlation = corSymm(c(0.99, 0.99, -0.99)), times = 0:2),
               "positive-definite")
  expect_error(manual(weights = varIdent(c("5" = 1.2), form = ~ 1 | v), times = 0:2),
               "not among `times`")
  expect_error(manual(weights = varIdent(c("0" = 1.2), form = ~ 1 | v), times = 0:2),
               "reference level")
  # nlme itself warns that it ignores values without a stratum.
  expect_error(manual(weights = suppressWarnings(varIdent(c("1" = 1.2))), times = 0:2),
               "by time")
  expect_error(manual(correlation = corCompSymm(0.3)), "random intercept")
})

# --- the field itself -----------------------------------------------------------

test_that("a hand-edited residual field is re-checked before use", {
  p <- manual(correlation = corCAR1(0.5))
  p$residual$coef[["Phi"]] <- 1.5
  expect_error(slope_var(p, 0:2), "Phi must lie in \\(0, 1\\)")
  p <- manual(correlation = corSymm(c(0.3, 0.2, 0.4)), times = 0:2)
  p$residual$times <- NULL
  expect_error(slope_var(p, 0:2), "`times` must be given exactly when")
})

test_that("an object without a residual field, from before it existed, still works", {
  p <- ref_params()
  p$residual <- NULL
  expect_false("residual" %in% names(p))
  expect_equal(slope_var(p, 0:2), slope_var(ref_params(), 0:2))
})

test_that("print.slope_params() shows the structure and its parameters", {
  out <- capture.output(print(residual_fit("un")))
  expect_true(any(grepl("residual variance at time 0", out, fixed = TRUE)))
  expect_true(any(grepl("correlation, times 0 and 3", out, fixed = TRUE)))
  expect_true(any(grepl("residual SD ratio, time 3", out, fixed = TRUE)))
  expect_true(any(grepl("not separately identified", out, fixed = TRUE)))
  out <- capture.output(print(manual(correlation = corSymm(c(0, 0, 0)), times = 0:2)))
  expect_false(any(grepl("not separately identified", out, fixed = TRUE)))
})

# --- bootstrap ------------------------------------------------------------------

test_that("bootstrap replicates refit the same residual structure", {
  p <- residual_fit("un")
  refit <- make_refitter(p)
  cl <- environment(refit)$cl
  expect_match(paste(deparse(cl$correlation), collapse = ""), "corSymm")
  expect_match(paste(deparse(cl$weights), collapse = ""), "varIdent")
  # Refitting the observed frame itself must give back the observed fit.
  q <- refit(boot_frame(p, "test"))
  expect_equal(q$residual, p$residual, tolerance = 1e-4)
  expect_equal(q$sigma2_residual, p$sigma2_residual, tolerance = 1e-4)

  fixed <- slope_params(sdmt ~ time | id, load_paper_data("slpower1"),
                        correlation = corCAR1(0.4, fixed = TRUE))
  cl <- environment(make_refitter(fixed))$cl
  expect_equal(eval(cl$correlation[["value"]]), 0.4)
  expect_true(cl$correlation[["fixed"]])
})
