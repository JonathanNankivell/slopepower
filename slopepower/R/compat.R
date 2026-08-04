# Layer 4 -- a translation layer for the Stata command's interface.
#
# This exists so that an existing Stata script can be ported mechanically, and so
# that the worked examples in Nash et al. (2021) can be run verbatim as a parity
# check. New code should prefer slope_params() + trial_design() +
# slope_sample_size() or slope_power().

#' Sample size or power using the Stata command's interface
#'
#' A direct translation of the Stata `slopepower` command of Nash et al. (2021)
#' onto this package's pipeline. Argument names, defaults and conventions follow
#' the Stata original: `schedule` lists follow-up visits only with baseline assumed
#' at time 0, `scale` divides the time variable, `dropouts` aligns element by
#' element with `schedule`, and the data type is declared with `obs`/`rct` rather
#' than by which grouping variable is supplied.
#'
#' New code should use [slope_params()], [trial_design()] and then
#' [slope_sample_size()] or [slope_power()] directly. Those separate the two
#' stages of the calculation, so the mixed model is fitted once and any number of
#' designs evaluated against it; they express the visit schedule in real time
#' rather than through the `scale` workaround; and they ask one question each,
#' rather than switching between sample size and power depending on which
#' argument was left out.
#'
#' @section Differences from the Stata command:
#'
#' The interface is faithful, but three defects in the original are deliberately
#' not reproduced. The `dropouts` length check is enforced here; in the Stata code
#' the counter that drives it is built with a malformed macro reference, and an
#' over-long list can silently reduce the completer weight and inflate the
#' estimated sample size. The check that dropout proportions do not exceed 100% is
#' applied with a tolerance; the original accumulates by repeated subtraction, so
#' `dropouts(0.3 0.3 0.4)` fails despite summing to exactly 1. And a malformed
#' string in the original aborts the command instead of warning when `treat()` is
#' given with observational data; here it warns, as intended.
#'
#' The `schedule` restriction to ascending integers of at least 1 is also lifted,
#' since this port builds the covariance matrix at the requested times rather than
#' on a unit-integer grid. `scale` is therefore never necessary, and a message
#' notes as much when it is used.
#'
#' @param data A data frame in long format with repeated measurements.
#' @param depvar,subject,time Column names, as strings.
#' @param schedule Numeric vector of follow-up visit times for the planned trial,
#'   in units of `time / scale`. Baseline at time 0 is implicit.
#' @param obs,rct Declare the data as observational or as a previous trial.
#'   Exactly one must be `TRUE`.
#' @param nocontrols Use with `obs` when every subject has the condition.
#' @param casecon Column name identifying cases (1) and healthy controls (0). Used
#'   with `obs`.
#' @param treat Column name identifying the treated arm (1) and control arm (0).
#'   Used with `rct`.
#' @param dropouts Numeric vector of incremental dropout proportions aligned with
#'   `schedule`.
#' @param scale Number of `time` units in one `schedule` unit. Defaults to 1.
#' @param alpha Two-sided significance level. Defaults to 0.05.
#' @param power,n Supply one, as in Stata. `n` gives the power that sample size
#'   achieves; `power` gives the sample size that reaches it. `power` defaults to
#'   0.8 when neither is given.
#' @param effectiveness Proportion of the slope difference the treatment should
#'   remove. Defaults to 0.25. Mutually exclusive with `usetrt`.
#' @param usetrt Target the treatment effect observed in the supplied trial data.
#'   Requires `rct`.
#' @param ... Passed to [slope_params()], e.g. `common_variance`.
#'
#' @return A `slope_sample_size` object, or a `slope_power` object when `n` is
#'   supplied. Both inherit from `slope_result`.
#'
#' @examples
#' # The paper's first worked example, p.588. Observational data with no
#' # comparison group, follow-up at 1 and 2 years, powering against a treatment
#' # expected to remove a third of the decline. Gives N = 712.
#' slopepower(slpower1, "sdmt", "id", "visit", schedule = c(1, 2),
#'            obs = TRUE, nocontrols = TRUE, effectiveness = 0.33)
#'
#' # p.594. A completed trial, powered against the effect it actually observed
#' # rather than an assumed proportion. The uneven schedule draws the message
#' # about `scale` described above; the Stata original would reject it outright.
#' slopepower(slpower3, "sdmt", "id", "visit", schedule = c(0.5, 2),
#'            rct = TRUE, treat = "treat", usetrt = TRUE)
#'
#' @seealso [slope_sample_size()], [slope_power()]
#' @export
slopepower <- function(data, depvar, subject, time, schedule,
                       obs = FALSE, rct = FALSE, nocontrols = FALSE,
                       casecon = NULL, treat = NULL,
                       dropouts = NULL, scale = 1,
                       alpha = 0.05, power = NULL, n = NULL,
                       effectiveness = NULL, usetrt = FALSE, ...) {
  context <- "slopepower()"

  if (!is.data.frame(data)) data <- as.data.frame(data)
  for (nm in c("depvar", "subject", "time")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L) {
      stop(sprintf("%s: `%s` must be a single column name.", context, nm), call. = FALSE)
    }
    if (!v %in% names(data)) {
      stop(sprintf("%s: `%s` = \"%s\" is not a column of `data`.", context, nm, v),
           call. = FALSE)
    }
  }

  # ---- model selection, following slopepower.ado lines 48-64 ----------------
  if (obs && rct) {
    stop(sprintf("%s: you cannot specify both `obs` and `rct`.", context), call. = FALSE)
  }
  if (!obs && !rct) {
    stop(sprintf("%s: no model specified; set either `obs = TRUE` or `rct = TRUE`.", context),
         call. = FALSE)
  }
  if (rct && nocontrols) {
    stop(sprintf("%s: you cannot specify `nocontrols` with `rct`.", context), call. = FALSE)
  }
  model <- if (obs && !nocontrols) 1L else if (obs) 2L else 3L

  if (model == 1L && is.null(casecon)) {
    stop(sprintf(paste0("%s: observational data with healthy controls needs `casecon`; ",
                        "for data in which every subject has the condition set ",
                        "`nocontrols = TRUE`."), context), call. = FALSE)
  }
  if (model == 3L && is.null(treat)) {
    stop(sprintf("%s: `rct = TRUE` requires `treat`.", context), call. = FALSE)
  }
  if (!is.null(casecon) && model != 1L) {
    warning(sprintf(paste0("%s: `casecon` applies only to observational data with both cases ",
                           "and healthy controls. It will be ignored."), context), call. = FALSE)
    casecon <- NULL
  }
  if (!is.null(treat) && model != 3L) {
    warning(sprintf("%s: `treat` applies only to trial data (`rct = TRUE`). It will be ignored.",
                    context), call. = FALSE)
    treat <- NULL
  }
  if (usetrt && model != 3L) {
    warning(sprintf("%s: `usetrt` applies only to trial data (`rct = TRUE`). It will be ignored.",
                    context), call. = FALSE)
    usetrt <- FALSE
  }
  if (usetrt && !is.null(effectiveness)) {
    stop(sprintf("%s: supply only one of `effectiveness` and `usetrt`, not both.", context),
         call. = FALSE)
  }

  group_arg <- if (model == 1L) "healthy" else if (model == 3L) "treated" else NULL
  group_col <- if (model == 1L) casecon else if (model == 3L) treat else NULL
  if (!is.null(group_col)) {
    if (!is.character(group_col) || length(group_col) != 1L) {
      stop(sprintf("%s: `%s` must be a single column name.", context,
                   if (model == 1L) "casecon" else "treat"), call. = FALSE)
    }
    if (!group_col %in% names(data)) {
      stop(sprintf("%s: \"%s\" is not a column of `data`.", context, group_col), call. = FALSE)
    }
  }

  # ---- schedule and scale --------------------------------------------------
  check_scalar(scale, "scale", context, lower = 0, upper = Inf, lower_open = TRUE)
  if (!is.numeric(schedule) || length(schedule) < 1L || any(!is.finite(schedule))) {
    stop(sprintf("%s: `schedule` must be a numeric vector of at least one follow-up visit time.",
                 context), call. = FALSE)
  }
  if (any(schedule <= 0)) {
    stop(sprintf(paste0("%s: `schedule` lists follow-up visits only; baseline at time 0 is ",
                        "implicit, so every entry must be positive."), context), call. = FALSE)
  }
  if (is.unsorted(schedule, strictly = TRUE)) {
    stop(sprintf("%s: `schedule` must be strictly increasing.", context), call. = FALSE)
  }
  if (any(schedule != floor(schedule))) {
    message(sprintf(paste0("%s: `schedule` contains non-integer visit times. The Stata command ",
                           "requires ascending integers of at least 1 and provides `scale` to ",
                           "compensate; this port evaluates the covariance at the requested ",
                           "times directly, so fractional visits are supported and `scale` is ",
                           "unnecessary."), context))
  }
  if (scale != 1) {
    message(sprintf(paste0("%s: `scale` = %s divides `time` before fitting. In new code, express ",
                           "`visits` directly in the units you want and omit `scale`."),
                    context, format(scale)))
  }

  # ---- fit stage one -------------------------------------------------------
  work <- data
  work[[".slopepower_time"]] <- as.numeric(coerce_time(data[[time]], context)) / scale

  fml <- stats::as.formula(sprintf("`%s` ~ `.slopepower_time` | `%s`", depvar, subject))

  env <- new.env(parent = parent.frame())
  assign(".slopepower_data", work, envir = env)
  call_args <- list(quote(slope_params), formula = fml, data = quote(.slopepower_data))
  if (!is.null(group_arg)) call_args[[group_arg]] <- as.name(group_col)
  params <- eval(as.call(c(call_args, list(...))), env)

  # ---- stage two -----------------------------------------------------------
  design <- trial_design(c(0, schedule), dropout = dropouts)

  args <- list(params = params, design = design, alpha = alpha)
  if (usetrt) {
    args$target <- "observed"
  } else {
    args$target <- "effectiveness"
    args$effectiveness <- effectiveness %||% 0.25
  }

  # The Stata command takes n() or power(), never both, and picks the calculation
  # from which was given. That single bimodal interface is the thing this port
  # splits in two, so the branch lives here -- in the compatibility shim whose
  # job is to mirror Stata -- rather than in the functions it delegates to.
  if (!is.null(n) && !is.null(power)) {
    stop(sprintf(paste0(
      "%s: supply only one of `n` and `power`.\n",
      "  `n` gives the power that sample size achieves; `power` gives the sample\n",
      "  size that reaches it."), context), call. = FALSE)
  }
  if (!is.null(n)) {
    do.call(slope_power, c(args, list(n = n)))
  } else {
    do.call(slope_sample_size, c(args, list(power = power %||% 0.8)))
  }
}
