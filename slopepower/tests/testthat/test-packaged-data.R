# The packaged datasets are a second copy of the .dta files in the repository
# root, converted by data-raw/make-data.R. Everything else in the suite now
# reads the packaged copies, so this file is what stops the two drifting: it
# rebuilds the conversion from the .dta files and demands an exact match.
#
# This is the one place that still needs the .dta files and haven, so it is also
# the one place that still skips -- from a built tarball there is nothing to
# compare against, and the rest of the suite runs regardless.

skip_without_dta <- function() {
  testthat::skip_if_not(
    !is.null(paper_data_dir()) && requireNamespace("haven", quietly = TRUE),
    "slpower*.dta not found, or haven not installed"
  )
}

read_dta_source <- function(which) {
  d <- as.data.frame(haven::read_dta(
    file.path(paper_data_dir(), paste0(which, ".dta"))))
  d$id   <- as.integer(d$id)
  d$sdmt <- as.numeric(d$sdmt)
  if (which == "slpower2") {
    d$case  <- as.integer(haven::zap_labels(d$case))
    d$vdate <- as.Date(d$vdate)
  } else {
    d$visit <- as.numeric(d$visit)
    if (which == "slpower3") d$treat <- as.integer(haven::zap_labels(d$treat))
  }
  d
}

test_that("the packaged datasets are identical to the .dta files they came from", {
  skip_without_dta()
  for (which in c("slpower1", "slpower2", "slpower3")) {
    expect_equal(get(which, envir = asNamespace("slopepower")),
                 read_dta_source(which),
                 info = which)
  }
})

# Shape and coding, checked without reference to the .dta files so that they
# hold in a tarball too. These are the properties the rest of the suite assumes.
test_that("the packaged datasets have the shape the paper describes", {
  expect_equal(dim(slopepower::slpower1), c(800L, 3L))
  expect_equal(dim(slopepower::slpower2), c(2000L, 4L))
  expect_equal(dim(slopepower::slpower3), c(450L, 4L))

  expect_named(slopepower::slpower1, c("id", "visit", "sdmt"))
  expect_named(slopepower::slpower2, c("id", "case", "vdate", "sdmt"))
  expect_named(slopepower::slpower3, c("id", "treat", "visit", "sdmt"))

  expect_equal(length(unique(slopepower::slpower1$id)), 200L)
  expect_equal(length(unique(slopepower::slpower2$id)), 500L)
  expect_equal(length(unique(slopepower::slpower3$id)), 150L)

  expect_equal(sort(unique(slopepower::slpower1$visit)), c(0, 1, 2, 3))
  expect_equal(sort(unique(slopepower::slpower3$visit)), c(0, 0.5, 2))
  expect_s3_class(slopepower::slpower2$vdate, "Date")

  # Group columns are plain 0/1 integers and constant within subject -- the
  # coding slope_params() requires, and the reason the labels are dropped.
  for (nm in c("slpower2", "slpower3")) {
    d <- get(nm, envir = asNamespace("slopepower"))
    g <- d[[if (nm == "slpower2") "case" else "treat"]]
    expect_type(g, "integer")
    expect_equal(sort(unique(g)), c(0L, 1L), info = nm)
    expect_true(all(tapply(g, d$id, function(x) length(unique(x))) == 1L), info = nm)
  }
})

test_that("the packaged data reproduces the paper's fitted slopes", {
  # The reason the datasets are worth shipping at all: pp. 588, 590 and 594.
  expect_equal(paper_fit("slpower1")$slope, -1.672, tolerance = 5e-4)
  expect_equal(paper_fit("slpower2")$slope, -1.715, tolerance = 5e-4)
  expect_equal(paper_fit("slpower3")$slope, -1.852, tolerance = 5e-4)
})
