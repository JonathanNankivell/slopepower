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
  labels <- tryCatch(
    attr(stats::terms(stats::as.formula(paste("~", paste(deparse(rhs), collapse = " ")),
                                        env = baseenv())),
         "term.labels"),
    error = function(e) NULL)
  if (!is.null(labels) && length(labels) != 1L) {
    stop(sprintf(paste0(
      "%s: the right-hand side must be a single time term, but `%s` has %d: %s.\n",
      "  This port models the outcome as a linear function of time only, as in\n",
      "  Nash et al. (2021); there is no covariate adjustment. A transformation\n",
      "  of time is fine, e.g. `outcome ~ I(vdate / 365) | subject`."),
      context, paste(deparse(rhs), collapse = " "), length(labels),
      paste(sQuote(labels), collapse = ", ")), call. = FALSE)
  }

  list(outcome = lhs, time = rhs, subject = subject)
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
