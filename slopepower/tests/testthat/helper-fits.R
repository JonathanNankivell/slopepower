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
    slpower2 = suppressMessages(slope_params(sdmt ~ time | id, d, healthy = case)),
    slpower3 = slope_params(sdmt ~ time | id, d, treated = treat)
  )
  cache[[which]] <- fit
  fit
}
