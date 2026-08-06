# Shared internal helpers. Read-only for the layer implementations.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Tolerance used when comparing dropout proportions against 1
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
  # friends are ordinary calls, not combining operators.
  combining <- c("+", "-", "*", "/", ":", "^", "%in%", "offset")
  if (is.call(rhs) && as.character(rhs[[1L]])[1L] %in% combining) {
    stop(sprintf(paste0(
      "%s: the right-hand side must be a single time term, but `%s` combines\n",
      "  terms with `%s`. This port models the outcome as a linear function of\n",
      "  time only, as in Nash et al. (2021); there is no covariate adjustment,\n",
      "  and the Stata original refuses such a formula at parse time.\n",
      "  To transform time, wrap the arithmetic in I(), e.g.\n",
      "  `outcome ~ I(vdate / 365) | subject`.\n",
      "  Baseline in particular needs no adjustment: it is modelled as a\n",
      "  correlated outcome with a single intercept for both arms (paper\n",
      "  section 2.1), rather than entered as a covariate."),
      context, paste(deparse(rhs), collapse = " "),
      as.character(rhs[[1L]])[1L]), call. = FALSE)
  }

  labels <- tryCatch(
    attr(stats::terms(stats::as.formula(paste("~", paste(deparse(rhs), collapse = " ")),
                                        env = baseenv())),
         "term.labels"),
    error = function(e) NULL)
  if (!is.null(labels) && length(labels) != 1L) {
    stop(sprintf(paste0(
      "%s: the right-hand side must be a single time term, but `%s` has %d: %s.\n",
      "  This port models the outcome as a linear function of time only, as in\n",
      "  Nash et al. (2021); there is no covariate adjustment, and the Stata\n",
      "  original refuses such a formula at parse time. A transformation of\n",
      "  time is fine, e.g. `outcome ~ I(vdate / 365) | subject`.\n",
      "  Baseline in particular needs no adjustment: it is modelled as a\n",
      "  correlated outcome with a single intercept for both arms (paper\n",
      "  section 2.1), rather than entered as a covariate."),
      context, paste(deparse(rhs), collapse = " "), length(labels),
      paste(sQuote(labels), collapse = ", ")), call. = FALSE)
  }

  list(outcome = lhs, time = rhs, subject = subject)
}

#' Format a numeric vector compactly for messages and printing
#' @noRd
fmt_num <- function(x) {
  vapply(x, function(v) format(v, trim = TRUE, drop0trailing = TRUE), character(1L))
}

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

#' Reject `effectiveness` alongside target = "observed"
#'
#' The one place this rule is enforced: called directly by `effect_components()`
#' and, separately, by the grid wrappers before their loop -- the grid loop
#' omits `effectiveness` from the call it builds when target = "observed",
#' which would otherwise bypass the check inside `effect_components()` itself.
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

#' Format a labelled value the way the Stata command does, for print methods
#' @noRd
fmt_line <- function(label, value, width = 39L, digits = 3L) {
  val <- if (is.character(value)) {
    value
  } else if (is.na(value)) {
    "."
  } else if (value == round(value) && abs(value) < 1e9 && digits == 0L) {
    format(value)
  } else {
    formatC(value, format = "f", digits = digits)
  }
  sprintf("%s = %s", formatC(label, width = width), val)
}
