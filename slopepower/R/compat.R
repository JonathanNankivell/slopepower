# Layer 4 -- a translation layer for the Stata command's interface.
#
# This exists so that an existing Stata script can be ported mechanically, and so
# that the worked examples in Nash et al. (2021) can be run verbatim as a parity
# check. New code should prefer slope_params() + trial_design() +
# slope_sample_size() or slope_power().

#' Warn and reset an argument that does not apply to the chosen model
#'
#' Shared shape of `slopepower()`'s three "you supplied this, but the model
#' you selected doesn't use it" checks: if `condition` holds, warn in
#' `message` (one `%s` for `context`) and answer `off_value`; otherwise leave
#' `value` untouched. `condition` already encodes both "was this supplied"
#' and "does the chosen model use it", since the two differ by argument
#' (`!is.null(x)` for `casecon`/`treat`, the bare flag for `usetrt`).
#' @noRd
warn_unused_arg <- function(value, condition, off_value, message, context) {
  if (!condition) return(value)
  warning(sprintf(message, context), call. = FALSE)
  off_value
}

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
#' The interface is faithful, but two defects in the original are deliberately
#' not reproduced. The check that dropout proportions do not exceed 100% is
#' applied to their sum with a tolerance; the original accumulates by repeated
#' subtraction, and `1 - .3 - .3 - .4` is `-5.551e-17`, so `dropouts(.3 .3 .4)`
#' trips a bare `< 0` guard despite summing to exactly 1. And the original's
#' warning for `treat()` given with observational data contains an unbalanced
#' macro quote, so what it prints -- or whether it aborts where a warning was
#' meant -- is anyone's guess; here it warns, as intended.
#'
#' The `dropouts` length check is *not* one of them: the counter driving it in
#' the Stata code is built with a space inside a macro name, which looks like it
#' would make the guard dead, but Stata trims the name and the guard fires
#' correctly. This wrapper enforces the same rule, and enforces it before the
#' mixed model is fitted, as the original does.
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
#' @param dropouts Numeric vector of incremental dropout proportions, one per
#'   entry of `schedule`: element `j` is the proportion whose **first missed
#'   visit is `schedule[j]`**, so their last attended visit is the one before it
#'   --- `schedule[j - 1]`, or the implicit baseline at time 0 when `j` is 1.
#'   `dropouts[1]` therefore describes participants seen at baseline only, who
#'   contribute nothing to a slope; that alignment is the Stata original's, and
#'   [trial_design()] prints both readings side by side. They are handled by the
#'   Dawson and Lagakos (1991, 1993) pattern mixture, as in the Stata original;
#'   see [trial_design()].
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
    check_column_name(get(nm), nm, data, context)
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
  casecon <- warn_unused_arg(casecon, !is.null(casecon) && model != 1L, NULL,
    paste0("%s: `casecon` applies only to observational data with both cases ",
          "and healthy controls. It will be ignored."), context)
  treat <- warn_unused_arg(treat, !is.null(treat) && model != 3L, NULL,
    "%s: `treat` applies only to trial data (`rct = TRUE`). It will be ignored.", context)
  usetrt <- warn_unused_arg(usetrt, usetrt && model != 3L, FALSE,
    "%s: `usetrt` applies only to trial data (`rct = TRUE`). It will be ignored.", context)
  if (usetrt && !is.null(effectiveness)) {
    stop(sprintf("%s: supply only one of `effectiveness` and `usetrt`, not both.", context),
         call. = FALSE)
  }

  group_arg <- if (model == 1L) "healthy" else if (model == 3L) "treated" else NULL
  group_col <- if (model == 1L) casecon else if (model == 3L) treat else NULL
  if (!is.null(group_col)) {
    check_column_name(group_col, if (model == 1L) "casecon" else "treat", data, context)
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

  # Built now, before the (potentially slow) stage-one fit below, so that a
  # purely syntactic problem with `dropouts` -- the wrong length, or a total
  # over 1 -- is reported without paying for a REML fit first.
  design <- trial_design(c(0, schedule), dropout = dropouts)

  # ---- fit stage one -------------------------------------------------------
  work <- data
  work[[".slopepower_time"]] <- coerce_time(data[[time]], context) / scale

  # env = baseenv(): every symbol this formula names (the backtick-quoted
  # column names) is resolved against `data` inside slope_params(), never
  # against the formula's own environment, so that environment need not be
  # this frame. Left at the default, as.formula() would environment it in
  # this call's own frame instead, and slope_params()'s match.call() stores
  # this exact formula object in `params$call`, a field kept for the life of
  # every slope_sample_size/slope_power result. That would pin slopepower()'s
  # whole frame -- `data` and `work`, a full second copy, included -- in
  # memory for as long as the result exists: the same leak
  # fit_none_model()/fit_treated_model() in params.R were written to avoid
  # for `params$fit`, one layer up.
  fml <- stats::as.formula(sprintf("`%s` ~ `.slopepower_time` | `%s`", depvar, subject),
                           env = baseenv())

  # Parented on this frame, which encloses the package namespace -- not on
  # parent.frame(), which roots the lookup in the caller's environment instead.
  # `slope_params` is a free symbol in the call built below, so from the
  # caller's environment it is found only when the package happens to be
  # attached: `slopepower::slopepower(...)` on its own failed with "could not
  # find function \"slope_params\"", and a user object of that name in the
  # global environment would have been picked up in preference to ours.
  # Nothing here needs the caller's environment: `...` has already been forced
  # by list(), so the only free symbols left are `slope_params`,
  # `.slopepower_data` and the group column name, which slope_params() resolves
  # against `data`.
  env <- new.env(parent = environment())
  assign(".slopepower_data", work, envir = env)
  call_args <- list(quote(slope_params), formula = fml, data = quote(.slopepower_data))
  if (!is.null(group_arg)) call_args[[group_arg]] <- as.name(group_col)
  params <- eval(as.call(c(call_args, list(...))), env)

  # ---- stage two -----------------------------------------------------------
  args <- list(params = params, design = design, alpha = alpha,
               target = if (usetrt) "observed" else "effectiveness")
  # Whether `effectiveness` may be passed alongside the chosen target is
  # maybe_add_effectiveness()'s rule to state, not this shim's to restate.
  args <- maybe_add_effectiveness(args, effectiveness %||% 0.25, args$target)

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
