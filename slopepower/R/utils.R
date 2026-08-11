# Shared internal helpers. Read-only for the layer implementations.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Slack allowed when comparing dropout proportions against their bounds
#'
#' Three uses, all in design.R and grid.R: an incremental vector's total
#' against 1, a cumulative vector's final element against 1, and a cumulative
#' vector's step-to-step differences against 0 (monotonicity). One constant so
#' the three cannot drift apart.
#'
#' Why any slack at all: a dropout vector that partitions the sample exactly
#' --- `c(0.3, 0.3, 0.4)` --- is legal and must be accepted, and floating-point
#' arithmetic on decimal proportions need not land on 1 exactly. In practice
#' `sum()` is well behaved here: across every 2- and 3-element partition of 1
#' into hundredths, and 200,000 random partitions into 2-6 parts of hundredths
#' or thousandths, it never once exceeded 1. So this is defensive rather than
#' load-bearing, and it is deliberately far larger than the error it guards
#' against --- proportions are a user-facing quantity, and nobody means a
#' dropout total of 1 + 1e-9.
#'
#' It is *not* a repair of the Stata original, though it reads like one.
#' `slopepower.ado:265` accumulates by repeated subtraction, which in doubles
#' would take `1 - 0.3 - 0.3 - 0.4` to -5.551e-17 and reject the list; but a
#' Stata local round-trips through a decimal string, so its residue is exactly
#' 0 and it accepts the list too. Verified, not assumed --- see DIVERGENCES.md
#' "Claims checked and rejected".
#' @noRd
DROPOUT_TOL <- 1e-8

#' Validate a single finite numeric scalar
#' @noRd
check_scalar <- function(x, name, context,
                         lower = -Inf, upper = Inf,
                         lower_open = TRUE, upper_open = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s: `%s` must be a single finite number.", context, name), call. = FALSE)
  }
  lo_ok <- if (lower_open) x > lower else x >= lower
  hi_ok <- if (upper_open) x < upper else x <= upper
  if (!lo_ok || !hi_ok) {
    stop(sprintf("%s: `%s` must be in %s%s, %s%s; got %g.",
                 context, name,
                 if (lower_open) "(" else "[", format(lower),
                 format(upper), if (upper_open) ")" else "]",
                 x), call. = FALSE)
  }
  invisible(as.numeric(x))
}

#' Check a numeric vector is all finite, naming the offending elements
#'
#' Shared by [validate_visits()] and `check_dropout_values()` in design.R,
#' which both stop with the same "must be finite; element(s) ... are not"
#' diagnosis over a vector rather than a scalar -- `check_scalar()` covers the
#' single-value case, this the vector one.
#' @noRd
check_finite_vector <- function(x, name, ctx) {
  if (any(!is.finite(x))) {
    stop(sprintf("%s: `%s` must be finite; element(s) %s are not.",
                 ctx, name, paste(which(!is.finite(x)), collapse = ", ")), call. = FALSE)
  }
  invisible(x)
}

#' Validate a probability in (0, 1)
#' @noRd
check_probability <- function(x, name, context) {
  check_scalar(x, name, context, lower = 0, upper = 1,
               lower_open = TRUE, upper_open = TRUE)
}

#' Validate a variance component (strictly positive)
#' @noRd
check_variance <- function(x, name, context) {
  check_scalar(x, name, context, lower = 0, upper = Inf, lower_open = TRUE)
}

#' Is a symmetric matrix positive definite?
#' @noRd
is_positive_definite <- function(m, tol = 1e-10) {
  if (!is.matrix(m) || nrow(m) != ncol(m)) return(FALSE)
  if (any(!is.finite(m))) return(FALSE)
  ev <- tryCatch(eigen(m, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  !is.null(ev) && all(ev > tol * max(1, abs(ev[1L])))
}

#' Check that a random-effects covariance triple is positive definite
#'
#' Shared by [new_slope_params()] -- the single validation point at
#' construction, for both routes into the class -- and [check_params()] in
#' power.R, which re-checks a hand-built object at the trust boundary every
#' stage-two calculation funnels through (a hand-built object can bypass the
#' constructor entirely, and the marginal covariance built in `sigma_at()` can
#' stay positive definite even when this one is not, because a large enough
#' `sigma2_residual` masks it). One helper means the matrix and the message can
#' never drift between the two call sites.
#' @noRd
check_re_covariance <- function(sigma2_intercept, sigma2_slope, sigma_cov, context) {
  G <- matrix(c(sigma2_intercept, sigma_cov, sigma_cov, sigma2_slope), 2L, 2L)
  if (!is_positive_definite(G)) {
    stop(sprintf(paste0("%s: the implied random-effects covariance matrix is not ",
                        "positive definite (var_int = %g, var_slope = %g, cov = %g)."),
                 context, sigma2_intercept, sigma2_slope, sigma_cov), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate a scalar is bounded and a whole number
#'
#' `check_scalar()` handles the bound; this adds the integrality check shared
#' by [run_bootstrap()]'s `R` and [solve_slope()]'s `n`, which differ only in
#' the bound and the noun used to name what is being counted.
#' @noRd
check_whole_number <- function(x, name, noun, context, lower, upper = Inf, lower_open = FALSE) {
  check_scalar(x, name, context, lower = lower, upper = upper, lower_open = lower_open)
  if (x != floor(x)) {
    stop(sprintf("%s: `%s` must be a whole number of %s; got %g.", context, name, noun, x),
         call. = FALSE)
  }
  invisible(as.numeric(x))
}

#' The two-sided critical value, computed exactly as the Stata original does
#'
#' `qnorm(1 - alpha / 2)` rather than the numerically preferable
#' `qnorm(alpha / 2, lower.tail = FALSE)`, because `slopepower.ado:634` writes
#' `invnormal(1 - 0.5 * alpha)` and the two differ in the last ulp for most
#' alphas -- enough, with `ceiling()` downstream, to step a published N by one.
#' Parity wins; this function exists so that the one place it is written is
#' also the one place its failure mode is handled.
#'
#' That failure mode: `1 - alpha / 2` rounds to exactly 1 once alpha falls
#' below about 2.2e-16, and `qnorm(1)` is `Inf`, so `size_per_arm()` returned
#' `n = Inf` -- "no sample size is enough" -- where the true requirement is
#' finite and unremarkable (7,584 at alpha = 1e-16 on the reference parameters,
#' against 712 at alpha = 0.05). Stata degenerates the same way and reports
#' N as missing; per CONTRACT.md section 6 the port says so instead.
#' [slope_sample_size_floor.slope_result()] guards the identical `qnorm(1)`
#' degeneracy on the power side.
#' @noRd
z_alpha <- function(alpha, context) {
  z <- stats::qnorm(1 - alpha / 2)
  if (!is.finite(z)) {
    stop(sprintf(paste0(
      "%s: `alpha` = %g is too small to compute a critical value with: 1 - alpha/2\n",
      "  is not distinguishable from 1 in double precision, so the two-sided z value\n",
      "  is infinite and no finite sample size can be reported. The Stata original\n",
      "  returns a missing N here. Use an alpha above about 1e-15."),
      context, alpha), call. = FALSE)
  }
  z
}

#' Extract the subject grouping from a `outcome ~ time | subject` formula
#'
#' Returns a list with `outcome`, `time` and `subject` language objects. The bar
#' is optional only when `subject` is supplied separately by the caller.
#' @noRd
parse_slope_formula <- function(formula, context) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop(sprintf("%s: `formula` must be two-sided, e.g. `outcome ~ time | subject`.",
                 context), call. = FALSE)
  }
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]
  subject <- NULL
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    subject <- rhs[[3L]]
    rhs <- rhs[[2L]]
  }

  # The subject side is held to the same rule as the time side below, and for
  # the same reason: it is *evaluated*, not expanded as a model formula, so
  # `id + site` becomes the arithmetic sum and `site/id` the quotient, and
  # `factor(as.character(...))` of either is a set of fabricated participant
  # identifiers that groups unrelated rows together. The fit then succeeds and
  # returns a plausible-looking slope -- on `slpower1`, `sdmt ~ visit | site/id`
  # gives -0.75 against the true -1.67. The hazard is not hypothetical: nested
  # and crossed groupings are exactly how `nlme` and `lme4` users write a
  # multi-level random-effects specification, and there is no such thing here
  # (see "further levels of clustering" in `?slope_params`), so the only
  # possible meanings of such an expression are a mistake.
  combining <- c("+", "-", "*", "/", ":", "^", "%in%", "offset")
  if (is.call(subject) && as.character(subject[[1L]])[1L] %in% combining) {
    stop(sprintf(paste0(
      "%s: the subject identifier must be a single grouping term, but `%s`\n",
      "  combines terms with `%s`. There is exactly one level of clustering here --\n",
      "  the participant -- so a nested or crossed grouping such as `site/id` or\n",
      "  `id + site` cannot be expressed; it would be evaluated arithmetically and\n",
      "  the result fitted as though it identified participants. Name the\n",
      "  participant identifier alone, e.g. `outcome ~ time | id`."),
      context, paste(deparse(subject), collapse = " "),
      as.character(subject[[1L]])[1L]), call. = FALSE)
  }

  # The time side must be exactly one model term. A transformation is fine --
  # `I(vdate / 365)` is one term -- but `time + age` is two, and without this
  # check it would be evaluated as the arithmetic sum and fitted as if it were
  # the time variable, silently returning a meaningless slope. The method of
  # Nash et al. models the outcome as a linear function of time alone; there is
  # no covariate adjustment, so extra terms can only be a mistake.
  # Two checks, because they catch different things. `terms()` applies formula
  # semantics while the expression is later *evaluated* with arithmetic
  # semantics, and the two disagree: `time - age` and `offset(age)` each yield a
  # single term label, so a count alone would let them through to be evaluated
  # as arithmetic and fitted as though the result were the time variable. So
  # reject any right-hand side whose head is an operator that combines or
  # removes model terms. A transformation stays legal because `I()`, `log()` and
  # friends are ordinary calls, not combining operators. `combining` is the list
  # shared with the subject-side check above.
  if (is.call(rhs) && as.character(rhs[[1L]])[1L] %in% combining) {
    single_term_error(context, rhs,
                      sprintf("combines terms with `%s`", as.character(rhs[[1L]])[1L]))
  }

  labels <- tryCatch(
    attr(stats::terms(stats::as.formula(paste("~", paste(deparse(rhs), collapse = " ")),
                                        env = baseenv())),
         "term.labels"),
    error = function(e) NULL)
  if (!is.null(labels) && length(labels) != 1L) {
    single_term_error(context, rhs,
                      sprintf("has %d: %s", length(labels),
                              paste(sQuote(labels), collapse = ", ")))
  }

  list(outcome = lhs, time = rhs, subject = subject)
}

#' Reject a right-hand side that is not a single time term
#'
#' Two checks in [parse_slope_formula()] reach the same conclusion by different
#' routes -- an operator that combines model terms, and a `terms()` expansion
#' with more than one label -- and owe the user the same explanation of why the
#' formula is narrow. `complaint` supplies only the clause that differs.
#' @noRd
single_term_error <- function(context, rhs, complaint) {
  stop(sprintf(paste0(
    "%s: the right-hand side must be a single time term, but `%s` %s.\n",
    "  This port models the outcome as a linear function of time only, as in\n",
    "  Nash et al. (2021); there is no covariate adjustment, and the Stata\n",
    "  original refuses such a formula at parse time. To transform time, wrap\n",
    "  the arithmetic in I(), e.g. `outcome ~ I(vdate / 365) | subject`.\n",
    "  Baseline in particular needs no adjustment: it is modelled as a\n",
    "  correlated outcome with a single intercept for both arms (paper\n",
    "  section 2.1), rather than entered as a covariate."),
    context, paste(deparse(rhs), collapse = " "), complaint), call. = FALSE)
}

#' Format a numeric vector compactly for messages and printing
#' @noRd
fmt_num <- function(x) {
  vapply(x, function(v) format(v, trim = TRUE, drop0trailing = TRUE), character(1L))
}

#' Render a numeric vector as a comma-separated list
#'
#' One rule for every place a vector of times or proportions is shown to the
#' user -- error messages, grid labels and the print methods alike -- so that
#' the same vector cannot render two ways. `fmt_call_vec()` is the deliberately
#' different rule that wraps the result in `c(...)` to make a copy-pasteable
#' call.
#' @noRd
label_numeric <- function(x) paste(fmt_num(x), collapse = ", ")

#' Find a fixed-effect name accounting for `lme`'s two possible interaction
#' spellings
#'
#' `lme` names an interaction after the order its components appear in the
#' model formula, so `sp_case * sp_time` yields `sp_case:sp_time` while
#' `sp_time * sp_case` yields `sp_time:sp_case`. The two fits are identical --
#' only the label differs -- so extraction must not depend on which spelling
#' arose.
#'
#' @return The matching name in `names(b)`, or `NA_character_` if neither
#'   spelling is present.
#' @noRd
resolve_fixef_name <- function(b, parts) {
  cand <- unique(c(paste(parts, collapse = ":"), paste(rev(parts), collapse = ":")))
  hit <- cand[cand %in% names(b)]
  if (length(hit)) hit[1L] else NA_character_
}

#' Add `effectiveness` to an argument list, unless target = "observed" implies it
#'
#' `target = "observed"` fixes `effectiveness` at 1 internally, and
#' [check_target_effectiveness()] rejects the two supplied together. Anything
#' that rebuilds a stage-two call must therefore omit `effectiveness` under
#' that target rather than supply it -- the rule belongs here once, rather than
#' being re-expressed at each call site: [slope_bootstrap()]'s `resolve_args()`,
#' [slopepower()] and both grid functions use this directly.
#' @noRd
maybe_add_effectiveness <- function(args, effectiveness, target) {
  if (!identical(target, "observed")) args$effectiveness <- effectiveness
  args
}

#' Reject `effectiveness` alongside target = "observed"
#'
#' The one place this rule is enforced, called by each of the four entry points
#' that can be handed both -- [slope_sample_size()], [slope_power()] and the two
#' grid wrappers -- immediately after `match.arg()`ing `target`. It has to sit at
#' that boundary rather than deeper in the calculation: "did the caller type an
#' `effectiveness`?" is a `missing()` question, and `missing()` can only be asked
#' of the function whose argument it is.
#' @noRd
check_target_effectiveness <- function(target, supplied, context) {
  if (identical(target, "observed") && isTRUE(supplied)) {
    stop(sprintf(paste0(
      "%s: supply only one of `effectiveness` and target = \"observed\".\n",
      "  target = \"observed\" reuses the treatment effect observed in the ",
      "previous trial, which fixes effectiveness at 1."), context), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate that a value is a single column name present in `data`
#'
#' Shared by the two places [slopepower()], the Stata-interface wrapper, takes
#' a column name as a string rather than a bare name: `depvar`/`subject`/`time`
#' and, separately, whichever of `casecon`/`treat` applies to the chosen model.
#' @noRd
check_column_name <- function(value, name, data, context) {
  if (!is.character(value) || length(value) != 1L) {
    stop(sprintf("%s: `%s` must be a single column name.", context, name), call. = FALSE)
  }
  if (!value %in% names(data)) {
    stop(sprintf("%s: `%s` = \"%s\" is not a column of `data`.", context, name, value),
         call. = FALSE)
  }
  invisible(value)
}

#' Format a labelled value the way the Stata command does, for print methods
#' @noRd
fmt_line <- function(label, value, width = 39L, digits = 3L) {
  val <- if (is.character(value)) {
    value
  } else if (is.na(value)) {
    "."
  } else if (digits == 0L && value == round(value) && abs(value) < 1e9) {
    # `scientific = FALSE` because this branch exists to print a whole number as
    # a whole number, and bare format() does not: it switches to scientific
    # notation whenever that is the shorter string, so a sample size of exactly
    # 100000 printed as "1e+05" and an 800-row model as "800" but a 100000-row
    # one as "1e+05". Stata's %5.0f never does that, and this block is a
    # transcription of Stata's output.
    format(value, scientific = FALSE)
  } else {
    formatC(value, format = "f", digits = digits)
  }
  sprintf("%s = %s", formatC(label, width = width), val)
}
