*! Grid 1 --- the stage-two arithmetic, densely.
*!
*! Everything downstream of the mixed-model fit is pure arithmetic: the Sigma*
*! construction, the two-person F* inverse, the Dawson-Lagakos dropout
*! weighting, the effectiveness rescaling, the two var_tte branches.  None of
*! it depends on which dataset produced the variance components, so it only
*! needs to be swept once --- and it should be swept on slpower1, whose
*! single-group model (Stata model 2) is by far the fastest to refit.  Stata
*! refits on every call, so this file is the one place a big cross is affordable.
*!
*! Two sweeps rather than one full cross, because the full cross is mostly
*! redundant:
*!
*!   A. schedule x dropout pattern, at fixed effectiveness/alpha/power.
*!      This is where the matrix algebra lives.  Irregular and single-visit
*!      schedules are included on purpose: the R port builds Sigma at arbitrary
*!      real times while Stata selects rows from a unit-integer grid, and the
*!      two are only guaranteed to agree at integer times.
*!
*!   B. effectiveness x alpha x (power|n), at two fixed designs, one with
*!      dropout and one without.  This is where the scaling lives --- in
*!      particular the effectiveness factor that appears both inside tte and
*!      again in the sample size formula, and the two different var_tte
*!      branches under dropout.
*!
*! Run with:  do 01_grid_arithmetic.do

clear all
set more off
version 13.0

do "_sprun.do"

local data ".."

use "`data'/slpower1.dta", clear

spopen "01_grid_arithmetic"

* ---------------------------------------------------------------------------
* Sweep A: schedule x dropout pattern
* ---------------------------------------------------------------------------
*
* Schedules chosen to cover: a single follow-up visit (the smallest design that
* identifies a slope), even spacing, uneven spacing, a long gap after a short
* one and the reverse, a schedule that skips baseline-adjacent times entirely,
* and a long dense schedule.

local scheds `" "1" "2" "3" "1 2" "1 3" "2 3" "1 2 3" "1 2 5" "1 2 3 4" "1 4 9" "2 4 6" "1 2 3 4 5 6" "'

local runA = 0

foreach s of local scheds {

	local K : word count `s'

	forvalues p = 1/6 {

		* ---- build the dropout list for this pattern and this length ----
		local d ""
		local pname ""

		if `p' == 1 {
			local pname "none"
		}
		if `p' == 2 {
			local pname "flat05"
			forvalues i = 1/`K' {
				local d "`d' 0.05"
			}
		}
		if `p' == 3 {
			* All dropout at the first follow-up visit.  These participants
			* attend baseline only, contribute no slope information, and are
			* the stratum both implementations must skip rather than treat as
			* an infinite sample size.
			local pname "baseonly20"
			local d "0.2"
			forvalues i = 2/`K' {
				local d "`d' 0"
			}
		}
		if `p' == 4 {
			* All dropout at the last visit: everyone who drops out does so
			* after seeing every visit but the final one.
			local pname "late20"
			forvalues i = 1/`K' {
				if `i' < `K' local d "`d' 0"
				else         local d "`d' 0.2"
			}
		}
		if `p' == 5 {
			* Monotone increasing hazard.
			local pname "graded"
			forvalues i = 1/`K' {
				local vs = string(0.02 * `i', "%6.4f")
				local d "`d' `vs'"
			}
		}
		if `p' == 6 {
			* Almost everybody drops out: the completer stratum carries 1% of
			* the weight and the strata carry the rest.  Exercises the far end
			* of the pattern-mixture sum.
			local pname "near99"
			local vs = string(0.99 / `K', "%8.6f")
			forvalues i = 1/`K' {
				local d "`d' `vs'"
			}
		}

		local d = trim("`d'")

		* pattern 3 is meaningless when K == 1 (it is then just "late")
		if `p' == 3 & `K' == 1 continue

		* ---- solve for N, then solve for power on the same design ----
		forvalues m = 1/2 {

			spreset
			global DS     "slpower1"
			global MODEL  "obs_nocont"
			global SCHED  "`s'"
			global DROPS  "`d'"
			global EFF    "0.33"
			global ALPHA  0.05

			local dopt ""
			if "`d'" != "" local dopt "dropouts(`d')"

			if `m' == 1 {
				global MODE "power"
				global PIN  0.8
				global TAG  "A-`s'-`pname'-N"
				global CMD  "slopepower sdmt, subject(id) time(visit) obs nocont schedule(`s') `dopt' effectiveness(0.33) alpha(0.05) power(0.8)"
			}
			else {
				global MODE "n"
				global NIN  450
				global TAG  "A-`s'-`pname'-P"
				global CMD  "slopepower sdmt, subject(id) time(visit) obs nocont schedule(`s') `dopt' effectiveness(0.33) alpha(0.05) n(450)"
			}

			sprun
			local runA = `runA' + 1
		}
	}
}

display as text _n "sweep A: `runA' runs"

* ---------------------------------------------------------------------------
* Sweep B: effectiveness x alpha x (power | n)
* ---------------------------------------------------------------------------
*
* Two designs only.  The no-dropout design exercises the closed-form var_tte
* (which is just F*[3,3] and should be invariant to alpha, power, n and
* effectiveness --- a strong internal check the R side can assert directly).
* The dropout design exercises the two back-solved var_tte branches, which are
* deliberately NOT the same number: solving for N carries the ceiling()
* rounding, solving for power does not.

local runB = 0

forvalues design = 1/2 {

	if `design' == 1 {
		local s "1 2 3"
		local d ""
		local dname "nodrop"
	}
	else {
		local s "1 2 3"
		local d "0.05 0.10 0.15"
		local dname "drop"
	}

	local dopt ""
	if "`d'" != "" local dopt "dropouts(`d')"

	foreach e in 0.05 0.25 0.33 0.5 1 {
		foreach a in 0.01 0.05 0.10 {

			foreach pw in 0.5 0.8 0.9 {
				spreset
				global DS     "slpower1"
				global MODEL  "obs_nocont"
				global SCHED  "`s'"
				global DROPS  "`d'"
				global EFF    "`e'"
				global ALPHA  `a'
				global MODE   "power"
				global PIN    `pw'
				global TAG    "B-`dname'-e`e'-a`a'-pw`pw'"
				global CMD    "slopepower sdmt, subject(id) time(visit) obs nocont schedule(`s') `dopt' effectiveness(`e') alpha(`a') power(`pw')"
				sprun
				local runB = `runB' + 1
			}

			foreach nn in 30 200 450 5000 {
				spreset
				global DS     "slpower1"
				global MODEL  "obs_nocont"
				global SCHED  "`s'"
				global DROPS  "`d'"
				global EFF    "`e'"
				global ALPHA  `a'
				global MODE   "n"
				global NIN    `nn'
				global TAG    "B-`dname'-e`e'-a`a'-n`nn'"
				global CMD    "slopepower sdmt, subject(id) time(visit) obs nocont schedule(`s') `dopt' effectiveness(`e') alpha(`a') n(`nn')"
				sprun
				local runB = `runB' + 1
			}
		}
	}
}

display as text _n "sweep B: `runB' runs"

spclose
