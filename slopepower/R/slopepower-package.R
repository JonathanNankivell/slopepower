#' @keywords internal
#'
#' @section Imports:
#' Every call in `R/` is written `nlme::lme()`, `stats::qnorm()` and so on, so
#' these `@importFrom` directives are not needed to make the code run. They are
#' kept because they are the package's declared surface onto its dependencies:
#' if `nlme` drops or renames one of these, loading `slopepower` fails
#' immediately with the offending name, rather than at the moment a user happens
#' to reach the branch that calls it.
#'
#' @importFrom nlme fixef getData getVarCov lme lmeControl
#' @importFrom nlme pdBlocked pdIdent pdSymm varIdent
#' @importFrom stats as.formula ave coef na.omit pnorm qnorm
#' @importFrom stats quantile setNames sigma terms vcov
"_PACKAGE"
