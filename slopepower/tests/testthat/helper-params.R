# A reference parameter set, shared by every file that needs a valid
# `slope_params` object without caring what is in it.
#
# The values are the simulation truth used to generate slpower1 (paper
# appendix), with the slope the paper reports for it. Tests that assert what the
# constructor *stores* keep their literals -- there the numbers are the point --
# but everything that just needs a well-formed object calls this, so a change to
# the class's requirements is one edit rather than nine.
ref_params <- function(comparator = "none", slope = -1.672,
                       slope_comparator = NA_real_) {
  slope_params_manual(
    slope            = slope,
    sigma2_intercept = 100,
    sigma2_slope     = 2,
    sigma_cov        = 5,
    sigma2_residual  = 10,
    slope_comparator = slope_comparator,
    comparator       = comparator
  )
}
