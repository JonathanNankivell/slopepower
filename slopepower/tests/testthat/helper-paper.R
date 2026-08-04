# Helpers for reproducing the published results of Nash et al. (2021),
# Stata Journal 21(3): 575-601.
#
# The paper's three simulated datasets ship with the package as slpower1,
# slpower2 and slpower3, so these tests run from a built tarball and under
# R CMD check rather than skipping. The .dta files they were built from still
# live in the parent of the package root; test-packaged-data.R checks the two
# against each other whenever those files are present.

# Locates the repository root, for the things that are *not* packaged: the raw
# .dta files and, as a fallback, ../stata-reference/. Returns NULL from a built
# tarball, where neither is available.
paper_data_dir <- function() {
  candidates <- c(
    file.path(testthat::test_path(), "..", "..", ".."),
    file.path(getwd(), "..", "..", ".."),
    file.path(getwd(), "..")
  )
  for (d in candidates) {
    if (file.exists(file.path(d, "slpower1.dta"))) return(normalizePath(d))
  }
  NULL
}

#' Load and prepare one of the paper's datasets
#'
#' Adds the `time` column in the units the paper works in, so that visit
#' schedules are expressed in the same units as the fitted time variable:
#'   slpower1  time = visit  (years)
#'   slpower2  time = vdate converted from days to years
#'   slpower3  time = visit  (years)
#'
#' The per-subject re-origining of time is left to `slope_params(origin =)`.
load_paper_data <- function(which = c("slpower1", "slpower2", "slpower3")) {
  which <- match.arg(which)
  d <- switch(which,
              slpower1 = slopepower::slpower1,
              slpower2 = slopepower::slpower2,
              slpower3 = slopepower::slpower3)
  d$time <- if (which == "slpower2") as.numeric(d$vdate) / 365 else d$visit
  d
}

# ---------------------------------------------------------------------------
# Published values. Page numbers refer to the Stata Journal article.
# ---------------------------------------------------------------------------

# Fitted slopes reported in the "Data characteristics" block of each example.
paper_slopes <- list(
  slpower1 = list(slope = -1.672, page = 588),
  slpower1_halfyear = list(slope = -0.836, page = 589),  # same fit, half-year axis
  slpower2 = list(slope = -1.715, comparator = 0.975,
                  difference = -2.690, page = 590),
  slpower3 = list(slope = -1.852, comparator = -1.104,
                  difference = -0.747, page = 594)
)

# Sample-size and power examples, sections 4.1.1 - 4.1.3.
paper_examples <- list(
  list(id = "p588_annual",   data = "slpower1", page = 588,
       visits = c(0, 1, 2), dropout = NULL, effectiveness = 0.33,
       expect = list(tte = 0.552, n = 712, n_per_arm = 356)),
  list(id = "p588_extended", data = "slpower1", page = 588,
       visits = c(0, 1, 2, 5), dropout = c(0, 0, 0.1), effectiveness = 0.33,
       expect = list(tte = 0.552, n = 328, n_per_arm = 164)),
  list(id = "p589_sixmonth", data = "slpower1", page = 589,
       visits = c(0, 0.5, 1, 1.5, 2), dropout = NULL, effectiveness = 0.33,
       expect = list(n = 620, n_per_arm = 310)),
  list(id = "p590_casecontrol", data = "slpower2", page = 590,
       visits = c(0, 1, 2), dropout = NULL, effectiveness = 0.33,
       comparator = "healthy",
       expect = list(tte = 0.888, n = 296, n_per_arm = 148)),
  list(id = "p593_power", data = "slpower2", page = 593,
       visits = c(0, 1, 2), dropout = c(0.05, 0.05), effectiveness = 0.33,
       comparator = "healthy", n = 200,
       expect = list(tte = 0.888, power = 0.597)),
  list(id = "p594_rct", data = "slpower3", page = 594,
       visits = c(0, 2, 3), dropout = c(0.2, 0.1), target = "observed",
       comparator = "treated",
       expect = list(tte = 0.747, n = 318, n_per_arm = 159))
)

# Table 1, p.595: slpower1, n = 450, effectiveness = 0.33.
paper_table1 <- data.frame(
  design  = rep(c("final_only", "annual", "six_month"), times = 3),
  dropout = rep(c("none", "5pc", "10pc"), each = 3),
  power   = c(0.798, 0.817, 0.865,
              0.732, 0.771, 0.828,
              0.648, 0.716, 0.783),
  stringsAsFactors = FALSE
)

table1_visits <- list(
  final_only = c(0, 3),
  annual     = 0:3,
  six_month  = seq(0, 3, by = 0.5)
)

# Per-year dropout converted to per-visit incremental proportions, exactly as
# the paper's appendix code does it.
table1_dropout <- list(
  none       = list(final_only = NULL, annual = NULL, six_month = NULL),
  `5pc`      = list(final_only = 0.15,
                    annual     = rep(0.05, 3),
                    six_month  = rep(0.025, 6)),
  `10pc`     = list(final_only = 0.30,
                    annual     = rep(0.10, 3),
                    six_month  = rep(0.05, 6))
)
