*! Grid 5 --- the fitted variance components themselves.
*!
*! Why this file exists. The first run of grids 1-3 came back with the fixed
*! effects agreeing to 1e-15 and var_tte disagreeing by about 1e-5 relative.
*! That is the signature of two REML optimisers stopping at slightly different
*! points, not of a formula error --- slpower1 is balanced, so the GLS slope is
*! the OLS slope regardless of the covariance structure, which is exactly why
*! the slope can agree to the last bit while the variance components do not.
*! Tightening nlme's control moved the R side from 8.36388645 to 8.36394022
*! against Stata's 8.36394172, which confirms it.
*!
*! But "confirms it" was an argument, not a measurement, because r(table) does
*! not carry the variance components. This file posts them, which buys two
*! things the other grids cannot give:
*!
*!   1. A direct comparison of the four components, so the size of the REML
*!      disagreement is a number rather than an inference.
*!   2. A tolerance-free test of the algebra: feed Stata's OWN components into
*!      slope_params_manual() on the R side and the two implementations must
*!      then agree on var_tte to machine precision. Any residual difference at
*!      that point is a real defect in Sigma*, X* or F*.
*!
*! It also settles open question 2 for every model at once. For each fit the
*! file posts the whole of e(b) with its column names, so the R side can read
*! the parameterisation off the names, and separately posts the components as
*! the .ado's own positional formulas compute them. If the two disagree, the
*! .ado's indexing is wrong; if they agree, :371 and :397 are correct.
*!
*! Run with:  do 05_variance_components.do

clear all
set more off
version 13.0

local data ".."

capture program drop vcpost
program define vcpost
	* Reads e(b) from the fit just run and posts one row.
	* Globals: VTAG VDS VMODEL VSCALE VSUBSET VCONTVAR
	tempname b
	matrix `b' = e(b)
	local k = colsof(`b')
	local names : colfullnames `b'

	forvalues i = 1/20 {
		local p`i' = .
	}
	forvalues i = 1/`k' {
		local p`i' = `b'[1, `i']
	}

	local nobs = e(N)
	tempname ng
	matrix `ng' = e(N_g)
	local nsub = `ng'[1,1]
	local conv = e(converged)

	* The .ado's own positional formulas, reproduced exactly. Line references
	* are to slopepower.ado v2.1.
	local slope     = .
	local slopecomp = .
	local vslope    = .
	local vint      = .
	local cov       = .
	local vres      = .

	if "$VMODEL" == "obs_cases" & $VCONTVAR == 0 {          // :364-:371
		local slopecomp = `p3'
		local slope     = `p3' + `p5'
		local vslope    = (exp(`p7'))^2
		local vint      = (exp(`p8'))^2
		local cov       = tanh(`p9') * exp(`p7') * exp(`p8')
		local vres      = (exp(`p13' + `p14'))^2
	}
	if "$VMODEL" == "obs_cases" & $VCONTVAR == 1 {          // :390-:397
		local slopecomp = `p3'
		local slope     = `p3' + `p5'
		local vslope    = (exp(`p7'))^2
		local vint      = (exp(`p8'))^2
		local cov       = tanh(`p9') * exp(`p7') * exp(`p8')
		local vres      = (exp(`p11' + `p12'))^2
	}
	if "$VMODEL" == "obs_nocont" {                          // :419-:425
		local slope  = `p1'
		local vslope = (exp(`p3'))^2
		local vint   = (exp(`p4'))^2
		local cov    = tanh(`p5') * exp(`p3') * exp(`p4')
		local vres   = (exp(`p6'))^2
	}
	if "$VMODEL" == "rct" {                                 // :449-:461
		local slopecomp = `p1'
		local slope     = `p1' + `p3'
		local vslope    = (exp(`p5'))^2
		local vint      = (exp(`p6'))^2
		local cov       = tanh(`p7') * exp(`p5') * exp(`p6')
		local vres      = (exp(`p8'))^2
	}

	* The competing reading of the two residual parameters under
	* residuals(independent, by()): two free log-SDs rather than base+offset.
	local vres_alt = .
	if "$VMODEL" == "obs_cases" & $VCONTVAR == 0 local vres_alt = (exp(`p14'))^2
	if "$VMODEL" == "obs_cases" & $VCONTVAR == 1 local vres_alt = (exp(`p12'))^2
	* And the other group's residual variance under the base+offset reading.
	local vres_ref = .
	if "$VMODEL" == "obs_cases" & $VCONTVAR == 0 local vres_ref = (exp(`p13'))^2
	if "$VMODEL" == "obs_cases" & $VCONTVAR == 1 local vres_ref = (exp(`p11'))^2

	post V ("$VTAG") ("$VDS") ("$VMODEL") ($VSCALE) ("$VSUBSET") ($VCONTVAR) ///
	       (`conv') (`nobs') (`nsub') (`k') ("`names'") ///
	       (`slope') (`slopecomp') (`vint') (`vslope') (`cov') (`vres') ///
	       (`vres_alt') (`vres_ref') ///
	       (`p1') (`p2') (`p3') (`p4') (`p5') (`p6') (`p7') ///
	       (`p8') (`p9') (`p10') (`p11') (`p12') (`p13') (`p14')

	display as text "  [`conv'] $VTAG  (`k' parameters)"
end

* Prepare the data exactly as slopepower does (ado :305-:332), then fit.
capture program drop vcfit
program define vcfit
	syntax , OUTcome(varname) SUBJect(varname) TIMe(varname) ///
	         [IFexp(string) GRoup(varname) SCale(real 1)]

	preserve
		if "`ifexp'" != "" keep if `ifexp'
		quietly replace `time' = `time' / `scale'
		quietly drop if missing(`outcome')
		tempvar first
		quietly bysort `subject' (`time') : gen double `first' = `time'[1]
		quietly bysort `subject' : replace `time' = `time' - `first'

		if "$VMODEL" == "obs_nocont" {
			quietly mixed `outcome' c.`time' || `subject': `time', ///
				cov(uns) reml iterate(16000)
		}
		if "$VMODEL" == "obs_cases" {
			quietly drop if missing(`group')
			tempvar cse ctl tc tk
			quietly gen byte `cse' = (`group' != 0)
			quietly gen byte `ctl' = (`group' == 0)
			quietly gen double `tc' = `time' * `cse'
			quietly gen double `tk' = `time' * `ctl'
			if $VCONTVAR == 0 {
				quietly mixed `outcome' `cse'##c.`time' ///
					|| `subject': `tc' `cse', cov(uns) nocons ///
					|| `subject': `tk' `ctl', cov(uns) nocons ///
					res(ind, by(`cse')) reml iterate(16000)
			}
			else {
				quietly mixed `outcome' `cse'##c.`time' ///
					|| `subject': `tc' `cse', cov(uns) nocons ///
					|| `subject': `ctl', cov(id) nocons ///
					res(ind, by(`cse')) reml iterate(16000)
			}
		}
		if "$VMODEL" == "rct" {
			quietly drop if missing(`group')
			tempvar plac
			quietly gen byte `plac' = (`group' == 0)
			quietly mixed `outcome' `time' `plac'#c.`time' ///
				|| `subject': `time', cov(uns) reml iterate(16000)
		}

		vcpost
	restore
end

capture program drop vcreset
program define vcreset
	global VTAG ""
	global VDS ""
	global VMODEL ""
	global VSCALE 1
	global VSUBSET ""
	global VCONTVAR 0
end

postfile V str40 vtag str10 dataset str12 model double scale str24 subset ///
          byte contvar byte converged double obs double subjects int nparam ///
          str500 colnames ///
          double slope double slope_comp double v_intercept double v_slope ///
          double cov_islope double v_residual ///
          double v_residual_alt double v_residual_ref ///
          double b1 double b2 double b3 double b4 double b5 double b6 ///
          double b7 double b8 double b9 double b10 double b11 double b12 ///
          double b13 double b14 ///
          using "05_variance_components.dta", replace

* ---------------------------------------------------------------------------
* Every distinct fit the other four grids rely on
* ---------------------------------------------------------------------------

* --- slpower1, single group, across the scale sweep and the subsets ---------
use "`data'/slpower1.dta", clear

foreach sc in 0.2 0.25 0.5 1 2 3 1000 0.001 {
	vcreset
	global VDS "slpower1"
	global VMODEL "obs_nocont"
	global VSCALE `sc'
	global VTAG "F1-sc`sc'"
	vcfit, outcome(sdmt) subject(id) time(visit) scale(`sc')
}

local subs `" "id<=50" "id<=100" "id>100" "visit<3" "visit>0" "visit!=1" "visit==0 | visit==3" "'
foreach sub of local subs {
	vcreset
	global VDS "slpower1"
	global VMODEL "obs_nocont"
	global VSUBSET "`sub'"
	global VTAG "F1-`sub'"
	vcfit, outcome(sdmt) subject(id) time(visit) ifexp(`sub')
}

* --- slpower2, cases and controls, both variance structures -----------------
use "`data'/slpower2.dta", clear

foreach sc in 91.25 182.5 365 730 {
	forvalues cv = 0/1 {
		vcreset
		global VDS "slpower2"
		global VMODEL "obs_cases"
		global VSCALE `sc'
		global VCONTVAR `cv'
		global VTAG "F2-sc`sc'-cv`cv'"
		vcfit, outcome(sdmt) subject(id) time(vdate) group(case) scale(`sc')
	}
}

* The same data pooled through the single-group model, and each group alone.
foreach spec in "" "case==1" "case==0" {
	local subname "`spec'"
	if "`spec'" == "" local subname "all"
	vcreset
	global VDS "slpower2"
	global VMODEL "obs_nocont"
	global VSCALE 365
	global VSUBSET "`spec'"
	global VTAG "F2-pooled-`subname'"
	if "`spec'" == "" {
		vcfit, outcome(sdmt) subject(id) time(vdate) scale(365)
	}
	else {
		vcfit, outcome(sdmt) subject(id) time(vdate) scale(365) ifexp(`spec')
	}
}

* --- slpower3, randomised trial ---------------------------------------------
use "`data'/slpower3.dta", clear

foreach sc in 0.25 0.5 1 2 {
	vcreset
	global VDS "slpower3"
	global VMODEL "rct"
	global VSCALE `sc'
	global VTAG "F3-sc`sc'"
	vcfit, outcome(sdmt) subject(id) time(visit) group(treat) scale(`sc')
}

local subs3 `" "id<=40 | id>110" "id>40 & id<=110" "visit<2" "visit>0" "'
foreach sub of local subs3 {
	vcreset
	global VDS "slpower3"
	global VMODEL "rct"
	global VSUBSET "`sub'"
	global VTAG "F3-`sub'"
	vcfit, outcome(sdmt) subject(id) time(visit) group(treat) ifexp(`sub')
}

foreach spec in "" "treat==1" "treat==0" {
	local subname "`spec'"
	if "`spec'" == "" local subname "all"
	vcreset
	global VDS "slpower3"
	global VMODEL "obs_nocont"
	global VSUBSET "`spec'"
	global VTAG "F3-pooled-`subname'"
	if "`spec'" == "" {
		vcfit, outcome(sdmt) subject(id) time(visit)
	}
	else {
		vcfit, outcome(sdmt) subject(id) time(visit) ifexp(`spec')
	}
}

postclose V

use "05_variance_components.dta", clear
export delimited using "05_variance_components.csv", replace
display as text "wrote 05_variance_components.csv, " _N " rows"

* ---------------------------------------------------------------------------
* Read this before trusting any model-1 row anywhere
* ---------------------------------------------------------------------------
list vtag v_residual v_residual_alt v_residual_ref if model == "obs_cases", ///
	noobs sepby(vtag) abbrev(16)

display as text _n "The R port estimates 10.354 for the cases' residual variance"
display as text "and 10.699 for the controls'. On the scale(365) rows above:"
display as text "  v_residual     is exp(b13+b14)^2, what slopepower.ado:371 computes"
display as text "  v_residual_alt is exp(b14)^2, the two-free-log-SDs reading"
display as text "  v_residual_ref is exp(b13)^2, the base alone"
display as text "Whichever column lands on 10.354 is the parameterisation, and"
display as text "if it is not v_residual then :371 is wrong."
