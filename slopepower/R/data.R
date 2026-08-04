# Documentation for the three packaged example datasets. The data themselves are
# built by data-raw/make-data.R from the Stata .dta files in the repository root.

#' Single-group longitudinal data
#'
#' The first worked example of Nash et al. (2021), used on pp. 588-589 and for
#' Table 1 on p. 595. Two hundred people with a progressive condition, each seen
#' four times a year apart, with no comparison group: the target trial is one
#' whose treatment aims to halt decline altogether, so the comparator slope is
#' zero rather than estimated.
#'
#' @format A data frame with 800 rows (200 subjects x 4 visits) and 3 columns:
#' \describe{
#'   \item{id}{Subject identifier, 1 to 200.}
#'   \item{visit}{Time of the visit in years: 0, 1, 2 or 3.}
#'   \item{sdmt}{Symbol Digit Modalities Test score, 0 to 68. Lower is worse, so
#'     the fitted slope is negative.}
#' }
#'
#' @source Simulated data distributed with the Stata `slopepower` command. The
#'   generating code is given in the appendix of Nash et al. (2021).
#'
#' @references
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#'
#' @examples
#' # The paper's p.588 fit: slope -1.672 points per year.
#' slope_params(sdmt ~ visit | id, data = slpower1)
#'
#' @seealso [slpower2] and [slpower3] for the two-group examples.
"slpower1"

#' Case-control longitudinal data
#'
#' The second worked example of Nash et al. (2021), pp. 590-593. Two hundred and
#' fifty cases and 250 healthy controls, each seen four times. Because the
#' controls are healthy rather than untreated patients, the comparator slope is
#' estimated from them and the trial is powered on the difference.
#'
#' Visit times are recorded as calendar dates rather than as a visit number, so
#' the time variable has to be built before fitting; see the examples. The paper
#' works in years, dividing the elapsed days by 365.
#'
#' @format A data frame with 2000 rows (500 subjects x 4 visits) and 4 columns:
#' \describe{
#'   \item{id}{Subject identifier, 1 to 500.}
#'   \item{case}{1 for a case, 0 for a healthy control. Constant within subject.}
#'   \item{vdate}{Date of the visit, between 2009-02-05 and 2012-12-22.}
#'   \item{sdmt}{Symbol Digit Modalities Test score, 0 to 83.}
#' }
#'
#' @source Simulated data distributed with the Stata `slopepower` command. The
#'   generating code is given in the appendix of Nash et al. (2021).
#'
#' @references
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#'
#' @examples
#' # Dates to years, as on p.590, then the paper's fit: slope -1.715.
#' d <- slpower2
#' d$time <- as.numeric(d$vdate) / 365
#' slope_params(sdmt ~ time | id, data = d, healthy = case)
#'
#' @seealso [slpower1], [slpower3].
"slpower2"

#' Randomised trial data
#'
#' The third worked example of Nash et al. (2021), p. 594. One hundred and fifty
#' people randomised between a treated and a control arm, each seen three times
#' over two years. Both slopes come from a completed trial, so the observed
#' treatment effect can be used directly as the effect to power against, via
#' `target = "observed"`, instead of assuming a proportion of the decline is
#' removed.
#'
#' The visit schedule is deliberately uneven -- 0, 6 months, 2 years -- which is
#' what makes this the example where the spacing of visits matters most.
#'
#' @format A data frame with 450 rows (150 subjects x 3 visits) and 4 columns:
#' \describe{
#'   \item{id}{Subject identifier, 1 to 150.}
#'   \item{treat}{1 for the treated arm, 0 for the control arm, 75 subjects each.
#'     Constant within subject.}
#'   \item{visit}{Time of the visit in years: 0, 0.5 or 2.}
#'   \item{sdmt}{Symbol Digit Modalities Test score, 2 to 58.}
#' }
#'
#' @source Simulated data distributed with the Stata `slopepower` command. The
#'   generating code is given in the appendix of Nash et al. (2021).
#'
#' @references
#' Nash, S., K. E. Morgan, C. Frost, and A. Mulick. 2021. Power and sample-size
#' calculations for trials that compare slopes over time: Introducing the
#' slopepower command. \emph{Stata Journal} 21(3): 575--601.
#'
#' @examples
#' # The paper's p.594 fit: treated slope -1.852, control -1.104.
#' slope_params(sdmt ~ visit | id, data = slpower3, treated = treat)
#'
#' @seealso [slpower1], [slpower2].
"slpower3"
