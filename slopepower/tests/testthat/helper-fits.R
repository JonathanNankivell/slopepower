# Cached stage-one fits for the paper's three datasets.
#
# The mixed-model fits are the slow part of the suite (the cases-and-controls
# model in particular), and several test files need the same three objects. This
# is also the architectural point of the two-stage split the port introduces:
# fit stage one once, then price out as many trial designs as you like. In Stata
# every row of Table 1 refits the model; here the fit is reused.

sp_fit_cache <- function() {
  if (!exists(".sp_fit_cache", envir = globalenv(), inherits = FALSE)) {
    assign(".sp_fit_cache", new.env(parent = emptyenv()), envir = globalenv())
  }
  get(".sp_fit_cache", envir = globalenv(), inherits = FALSE)
}

paper_fit <- function(which) {
  cache <- sp_fit_cache()
  if (!is.null(cache[[which]])) return(cache[[which]])
  d <- load_paper_data(which)
  fit <- switch(
    which,
    slpower1 = slope_params(sdmt ~ time | id, d),
    slpower2 = suppressMessages(slope_params(sdmt ~ time | id, d, comparator = "healthy", group = case)),
    slpower3 = slope_params(sdmt ~ time | id, d, comparator = "treated", group = treat)
  )
  cache[[which]] <- fit
  fit
}

# Ten subjects, a small model that fits quickly and (at ci_method = "bca") needs
# only ten leave-one-subject-out refits rather than slpower1's two hundred.
# Built on first use and cached, so it is fitted once for the whole suite
# rather than once per test, and not at all by files that never ask.
small_fit <- function() {
  cache <- sp_fit_cache()
  if (is.null(cache$small)) {
    set.seed(4)
    subj <- data.frame(id = 1:10, a = rnorm(10, 50, 8), b = rnorm(10, -2, 0.6))
    d <- merge(subj, data.frame(visit = 0:3))
    d$y <- d$a + d$b * d$visit + rnorm(nrow(d), 0, 2)
    cache$small <- slope_params(y ~ visit | id, d)
  }
  cache$small
}

# The main bootstrapped quantity of a single-design bootstrap, under short
# generic names, for tests that check its numbers rather than its field names
# (test-output-names.R checks those). A sample size is read on its stored
# total basis.
boot_main <- function(b) {
  q <- slopepower:::stat_field(b$statistic)
  list(observed = b[[q]],
       replicates = b[[paste0(q, "_replicates")]],
       ci = c(b[[paste0(q, "_lower")]], b[[paste0(q, "_upper")]]),
       boot_mean = b[[paste0(q, "_mean")]],
       boot_sd = b[[paste0(q, "_sd")]],
       ci_method = b[[paste0(b$statistic, "_ci_method")]])
}
