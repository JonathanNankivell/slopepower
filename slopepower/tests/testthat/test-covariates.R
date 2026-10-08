# Covariate adjustment in slope_params().

cov_data <- function(seed = 42) {
  set.seed(seed)
  d <- slpower1
  ids <- unique(d$id)
  age <- stats::setNames(stats::rnorm(length(ids), 45, 10), ids)
  sex <- stats::setNames(factor(sample(c("F", "M"), length(ids), TRUE)), ids)
  d$age <- age[as.character(d$id)]
  d$sex <- sex[as.character(d$id)]
  # age is prognostic for slope, so adjusting for it should shrink sigma2_slope
  d$sdmt <- d$sdmt - 0.08 * (d$age - 45) * d$visit
  d
}

test_that("adjusted fit matches a hand-built lme with centred covariates", {
  d <- cov_data()
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + sex))
  d$t <- d$visit - stats::ave(d$visit, d$id, FUN = min)
  first <- !duplicated(d$id)
  d$age_c <- d$age - mean(d$age[first])
  d$m     <- (d$sex == "M") - mean(d$sex[first] == "M")
  h <- nlme::lme(sdmt ~ t + age_c + m + age_c:t + m:t, random = ~ t | id,
                 data = d, method = "REML", control = slope_lme_control())
  G <- nlme::getVarCov(h)
  expect_equal(p$slope, nlme::fixef(h)[["t"]], tolerance = 1e-6)
  expect_equal(p$sigma2_slope, G["t", "t"], tolerance = 1e-4)
  expect_equal(p$sigma2_intercept, G["(Intercept)", "(Intercept)"], tolerance = 1e-4)
  expect_equal(p$sigma2_residual, stats::sigma(h)^2, tolerance = 1e-4)
})

test_that("a prognostic covariate-by-time term shrinks sigma2_slope; intercept-only does not", {
  d <- cov_data()
  u  <- suppressMessages(slope_params(sdmt ~ visit | id, d))
  a  <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age))
  a0 <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age,
                                      covariate_time = FALSE))
  expect_lt(a$sigma2_slope, 0.8 * u$sigma2_slope)
  expect_equal(a0$sigma2_slope, u$sigma2_slope, tolerance = 1e-3)
})

test_that("no covariates leaves the fit unchanged", {
  d <- cov_data()
  u <- suppressMessages(slope_params(sdmt ~ visit | id, d))
  n <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = NULL))
  expect_identical(u[c("slope", "sigma2_slope", "sigma2_residual")],
                   n[c("slope", "sigma2_slope", "sigma2_residual")])
  expect_identical(names(nlme::fixef(u$fit)), c("(Intercept)", "sp_time"))
})

test_that("covariates work with treated and healthy comparators", {
  set.seed(3)
  t3 <- slpower3
  t3$x <- stats::ave(t3$id, t3$id, FUN = function(i) stats::rnorm(1))
  p <- slope_params(sdmt ~ visit | id, t3, comparator = "treated", group = treat, covariates = ~ x)
  expect_true(all(c("sp_cov_1", "sp_time:sp_cov_1") %in% names(nlme::fixef(p$fit))))
  expect_false(is.na(p$slope_comparator))

  h2 <- slpower2
  h2$x <- stats::ave(h2$id, h2$id, FUN = function(i) stats::rnorm(1))
  q <- suppressMessages(slope_params(sdmt ~ I(as.numeric(vdate) / 365) | id, h2,
                                     comparator = "healthy", group = case, covariates = ~ x))
  expect_equal(q$comparator, "healthy")
  expect_false(is.na(q$slope_comparator))
})

test_that("bootstrap refits keep the covariate adjustment", {
  d <- cov_data()
  a <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + sex))
  fr <- boot_frame(a, "test")
  expect_true(all(c("sp_cov_1", "sp_cov_2") %in% names(fr)))
  refit <- make_refitter(a)(fr)
  expect_equal(refit$sigma2_slope, a$sigma2_slope, tolerance = 1e-6)
  expect_equal(names(nlme::fixef(refit$fit)), names(nlme::fixef(a$fit)))
})

test_that("bad covariate specifications are refused", {
  d <- cov_data()
  d$tv <- d$visit
  expect_error(suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ tv)),
               "constant within a participant")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = "age"),
               "one-sided formula")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = ~ nope),
               "could not evaluate `covariates`")
})

test_that("a participant with no covariate value is dropped through na.action", {
  d <- cov_data()
  d$age[d$id == d$id[1]] <- NA
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age))
  expect_equal(p$n_obs, sum(d$id != d$id[1]))
  expect_equal(p$n_subjects, length(unique(d$id)) - 1L)
})

test_that("a covariate recorded on the baseline row only is carried to every visit", {
  d <- cov_data()
  full <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + sex))
  baseline <- !duplicated(d$id)
  d$age[!baseline] <- NA
  d$sex[!baseline] <- NA
  # some participants recorded on every visit, the rest at baseline only: the
  # fit must not quietly fall back to the fully recorded ones
  keep <- d$id %in% unique(d$id)[1:5]
  d$age[keep] <- cov_data()$age[keep]
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + sex))
  expect_equal(p$n_obs, nrow(d))
  expect_equal(p$slope, full$slope, tolerance = 1e-10)
  expect_equal(p$sigma2_slope, full$sigma2_slope, tolerance = 1e-10)
})

test_that("recorded covariate values that disagree within a participant are still refused", {
  d <- cov_data()
  rows <- which(d$id == d$id[1])
  d$age[rows[1]] <- 40
  d$age[rows[2]] <- 50
  d$age[rows[-(1:2)]] <- NA
  expect_error(suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age)),
               "`age` changes during follow-up for 1 participant")
})

test_that("an unused factor level does not make the fit singular", {
  d <- cov_data()
  d$site <- factor(ifelse(d$id %% 2 == 0, "A", "B"), levels = c("A", "B", "C"))
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ site))
  d$site <- droplevels(d$site)
  q <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ site))
  expect_equal(p$slope, q$slope, tolerance = 1e-10)
  expect_equal(p$sigma2_slope, q$sigma2_slope, tolerance = 1e-10)
  expect_identical(names(nlme::fixef(p$fit)),
                   c("(Intercept)", "sp_time", "sp_cov_1", "sp_time:sp_cov_1"))
})

test_that("a level whose rows are all removed as missing is dropped too", {
  d <- cov_data()
  first <- d$id == d$id[1]
  d$site <- factor(ifelse(first, "C", ifelse(d$id %% 2 == 0, "A", "B")))
  d$sdmt[first] <- NA
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ site))
  expect_length(grep("^sp_cov_[0-9]+$", names(nlme::fixef(p$fit))), 1L)
})

test_that("covariates constant across all participants warn and adjust nothing", {
  d <- cov_data()
  d$k <- 7
  u <- suppressMessages(slope_params(sdmt ~ visit | id, d))
  expect_warning(p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ k)),
                 "no adjustment was made")
  expect_identical(names(nlme::fixef(p$fit)), c("(Intercept)", "sp_time"))
  expect_equal(p$slope, u$slope)
})

test_that("bootstrap and jackknife refits survive a rare factor level going missing", {
  d <- cov_data()
  ids <- unique(d$id)
  d$grade <- factor(ifelse(d$id == ids[1], "rare", ifelse(d$id %% 2 == 0, "A", "B")))
  a <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ grade))
  fr <- boot_frame(a, "test")
  refit <- make_refitter(a)
  # leave-one-out without the only "rare" participant, as the BCa jackknife does
  p <- refit(fr[fr$subject != fr$subject[d$id == ids[1]][1], , drop = FALSE])
  expect_s3_class(p, "slope_params")
  expect_length(grep("^sp_cov_[0-9]+$", names(nlme::fixef(p$fit))), 1L)
  b <- slope_params_boot(a, R = 40, seed = 2)
  expect_equal(b$n_failed, 0)
})

test_that("model.matrix-style covariate formulas are expanded as documented", {
  d <- cov_data()
  labs <- function(f) {
    unname(covariate_basis(f, covariate_variables(f, d, "test"), d$id,
                           rep(TRUE, nrow(d)), "test")$labels)
  }
  expect_equal(labs(~ age * sex), c("age", "sexM", "age:sexM"))
  expect_equal(labs(~ 0 + sex), "sexM")      # always coded against an intercept
  expect_equal(labs(~ sex - 1), "sexM")
  expect_length(labs(~ splines::ns(age, 3)), 3L)
  # poly() is computed through a QR, so equal ages differ in the last bit;
  # that must not count as varying within a participant
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ poly(age, 2)))
  expect_s3_class(p, "slope_params")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = ~ .), "cannot use `.`")
  expect_error(slope_params(sdmt ~ visit | id, d, covariates = ~ age + offset(age)),
               "offset")
  expect_error(suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + visit)),
               "`visit` changes during follow-up")
})

# ---- review fixes: covariate adjustment ------------------------------------

healthy_cov_data <- function(seed = 7) {
  set.seed(seed)
  h <- slpower2
  ids <- unique(h$id)
  case <- tapply(h$case, h$id, `[`, 1L)[as.character(ids)]
  # cases older than controls, and age prognostic for slope
  age <- stats::setNames(stats::rnorm(length(ids), ifelse(case == 1, 55, 40), 8), ids)
  h$age <- age[as.character(h$id)]
  h$t <- as.numeric(h$vdate) / 365
  h$sdmt <- h$sdmt - 0.05 * (h$age - 47) * (h$t - stats::ave(h$t, h$id, FUN = min))
  h
}

test_that("the reduced-structure note does not claim case estimates are unaffected under covariates", {
  h <- healthy_cov_data()
  p <- suppressMessages(slope_params(sdmt ~ t | id, h, comparator = "healthy", group = case,
                                     covariates = ~ age, common_variance = TRUE))
  out <- capture.output(print(p))
  expect_false(any(grepl("Case estimates are unaffected", out, fixed = TRUE)))
  expect_true(any(grepl("can shift the case", out, fixed = TRUE)))
})

test_that("poly() works with a covariate recorded on the baseline row only", {
  d <- cov_data()
  full <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ poly(age, 2)))
  d$age[duplicated(d$id)] <- NA
  p <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ poly(age, 2)))
  expect_equal(p$slope, full$slope, tolerance = 1e-10)
  expect_equal(p$sigma2_slope, full$sigma2_slope, tolerance = 1e-10)
  # and a participant with no value at all is removed rather than failing poly()
  d$age[d$id == d$id[1]] <- NA
  q <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ poly(age, 2)))
  expect_equal(q$n_subjects, length(unique(d$id)) - 1L)
})

test_that("a data-dependent basis does not depend on how many visits participants have", {
  d <- cov_data()
  f <- ~ splines::ns(age, 3)
  per_participant <- function(dd) {
    b <- covariate_basis(f, covariate_variables(f, dd, "test"), dd$id,
                         rep(TRUE, nrow(dd)), "test")$X
    b[!duplicated(dd$id), , drop = FALSE]
  }
  # duplicate every visit of the first ten participants: their weight per visit
  # doubles, but they are still the same participants
  extra <- d[d$id %in% unique(d$id)[1:10], ]
  expect_equal(unname(per_participant(rbind(d, extra))), unname(per_participant(d)))
})

test_that("a covariate aliased with the group or another covariate is named, not fitted", {
  h <- healthy_cov_data()
  h$case_copy <- h$case
  expect_error(suppressWarnings(suppressMessages(
    slope_params(sdmt ~ t | id, h, comparator = "healthy", group = case, covariates = ~ case_copy))),
    "fixed by `healthy`")
  d <- cov_data()
  d$age2 <- 2 * d$age + 1
  expect_error(suppressMessages(
    slope_params(sdmt ~ visit | id, d, covariates = ~ age + age2)),
    "`age2`, `time:age2` are an exact linear combination")
})

test_that("under healthy, covariates are centred on the cases", {
  h <- healthy_cov_data()
  p <- suppressMessages(slope_params(sdmt ~ t | id, h, comparator = "healthy", group = case, covariates = ~ age))
  g <- nlme::getData(p$fit)
  first <- !duplicated(g$sp_subject)
  expect_equal(mean(g$sp_cov_1[first & g$sp_case == 1]), 0, tolerance = 1e-10)
  # the case slope is the slope at the cases' mean age, which a pooled
  # centring would have shifted by delta * (case mean - pooled mean)
  ages <- tapply(h$age, h$id, `[`, 1L)
  cases <- tapply(h$case, h$id, `[`, 1L) == 1
  delta <- nlme::fixef(p$fit)[["sp_time:sp_cov_1"]]
  b <- nlme::fixef(p$fit)
  expect_equal(p$slope, b[["sp_time"]] + b[["sp_case:sp_time"]])
  pooled_slope <- p$slope + delta * (mean(ages) - mean(ages[cases]))
  expect_gt(abs(p$slope - pooled_slope), 0.1)
})

test_that("the adjustment is recorded on the object and printed through stage two", {
  d <- cov_data()
  u <- suppressMessages(slope_params(sdmt ~ visit | id, d))
  expect_null(u$covariates)
  a <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age + sex))
  expect_identical(a$covariates, list(columns = c("age", "sexM"), time = TRUE))
  a0 <- suppressMessages(slope_params(sdmt ~ visit | id, d, covariates = ~ age,
                                      covariate_time = FALSE))
  expect_false(a0$covariates$time)
  expect_true(grepl("adjusted for baseline covariates age, sexM, and their interactions",
                    paste(trimws(capture.output(print(a))), collapse = " "), fixed = TRUE))
  n <- slope_sample_size(a, c(0, 1, 2), effectiveness = 0.5)
  expect_true(any(grepl("adjusted for baseline covariates",
                        capture.output(print(n)), fixed = TRUE)))
  expect_false(any(grepl("adjusted for",
                         capture.output(print(slope_sample_size(u, c(0, 1, 2),
                                                                effectiveness = 0.5))),
                         fixed = TRUE)))
})
