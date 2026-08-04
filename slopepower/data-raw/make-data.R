# Regenerate data/slpower1.rda, data/slpower2.rda, data/slpower3.rda from the
# Stata .dta files in the repository root.
#
# Run from the package root:  Rscript data-raw/make-data.R
#
# The .dta files are the simulated datasets of Nash et al. (2021), whose
# generating code is given in the paper's appendix. Only storage types are
# changed here -- every value is carried across untouched, so that fits on the
# packaged data reproduce the published numbers exactly. Stata's value labels
# are dropped (`case` and `treat` become plain 0/1 integers) and Stata's
# elapsed-day `vdate` becomes a `Date`.

src <- file.path("..", c("slpower1.dta", "slpower2.dta", "slpower3.dta"))
stopifnot(all(file.exists(src)))

read_one <- function(path) as.data.frame(haven::read_dta(path))

slpower1 <- read_one(src[1])
slpower1$id    <- as.integer(slpower1$id)
slpower1$visit <- as.numeric(slpower1$visit)
slpower1$sdmt  <- as.numeric(slpower1$sdmt)

slpower2 <- read_one(src[2])
slpower2$id    <- as.integer(slpower2$id)
slpower2$case  <- as.integer(haven::zap_labels(slpower2$case))
slpower2$sdmt  <- as.numeric(slpower2$sdmt)
# vdate arrives from haven as a Date already; keep it as one.
slpower2$vdate <- as.Date(slpower2$vdate)

slpower3 <- read_one(src[3])
slpower3$id    <- as.integer(slpower3$id)
slpower3$treat <- as.integer(haven::zap_labels(slpower3$treat))
slpower3$visit <- as.numeric(slpower3$visit)
slpower3$sdmt  <- as.numeric(slpower3$sdmt)

dir.create("data", showWarnings = FALSE)
save(slpower1, file = "data/slpower1.rda", version = 2, compress = "xz")
save(slpower2, file = "data/slpower2.rda", version = 2, compress = "xz")
save(slpower3, file = "data/slpower3.rda", version = 2, compress = "xz")

# Guard: the packaged data must still reproduce the paper's fitted slopes.
# slpower1 p.588, slpower2 p.590, slpower3 p.594.
if (requireNamespace("slopepower", quietly = TRUE)) {
  fit <- function(d, ...) slopepower::slope_params(sdmt ~ time | id, d, ...)
  slpower1$time <- slpower1$visit
  slpower2$time <- as.numeric(slpower2$vdate) / 365
  slpower3$time <- slpower3$visit
  got <- c(
    slpower1 = fit(slpower1)$slope,
    slpower2 = suppressMessages(fit(slpower2, healthy = case))$slope,
    slpower3 = fit(slpower3, treated = treat)$slope
  )
  want <- c(slpower1 = -1.672, slpower2 = -1.715, slpower3 = -1.852)
  stopifnot(all(abs(got - want) < 5e-4))
  cat("paper slopes reproduced:", paste(sprintf("%.3f", got), collapse = " "), "\n")
}
