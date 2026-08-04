*! Grid 2 --- the comparator branches.
*!
*! Grid 1 sweeps the arithmetic on a single-group fit.  What it cannot reach is
*! how the reference slope is chosen, and that is the part of the port with the
*! most branches per line of code:
*!
*!   model 1 (obs + casecon)  reference = the healthy controls' slope
*!   model 2 (obs + nocont)   reference = 0
*!   model 3 (rct)            reference = 0        when effectiveness is used
*!                            reference = treated  when usetrt is used
*!
*! The third row is the one worth hammering.  With rct data and no usetrt,
*! Stata IGNORES the observed treated-arm slope and measures the effect toward
*! zero; the treated slope is estimated, reported in the on-screen block, and
*! then dropped from r(table).  That is deliberate, it is the default, and it
*! is exactly the kind of thing a port silently "corrects".
*!
*! Also covered: nocontvar (model 1 only --- the controls' random-effects block
*! is collapsed to a single intercept variance), and the interaction of usetrt
*! with effectiveness, which Stata reports as missing rather than 1.
*!
*! These fits are slower than slpower1's --- model 1 in particular is a
*! four-variance-component model on 500 subjects, and Stata refits it on every
*! call --- so the crosses here are deliberately narrower.
*!
*! Run with:  do 02_grid_comparators.do

clear all
set more off
version 13.0

do "_sprun.do"

local data ".."

spopen "02_grid_comparators"

* ---------------------------------------------------------------------------
* slpower2 --- model 1, cases and healthy controls
* ---------------------------------------------------------------------------
*
* time is a date in days; scale(365) puts the fit on the year axis the paper
* uses, which is also what the R helper does with vdate / 365.  The pairing of
* scale() with a date variable is itself worth a reference row: it is the only
* place in the paper where the time variable is not already in the units the
* schedule is written in.

use "`data'/slpower2.dta", clear

local scheds2 `" "1" "1 2" "1 2 3" "1 2 5" "2 4 6" "'

foreach s of local scheds2 {

	local K : word count `s'

	forvalues p = 1/2 {

		local d ""
		local pname "none"
		if `p' == 2 {
			local pname "flat05"
			forvalues i = 1/`K' {
				local d "`d' 0.05"
			}
			local d = trim("`d'")
		}

		local dopt ""
		if "`d'" != "" local dopt "dropouts(`d')"

		forvalues cv = 0/1 {

			local cvopt ""
			local cvname "contvar"
			if `cv' == 1 {
				local cvopt "nocontvar"
				local cvname "nocontvar"
			}

			* solve for N
			spreset
			global DS      "slpower2"
			global MODEL   "obs_cases"
			global SCALE   365
			global SCHED   "`s'"
			global DROPS   "`d'"
			global EFF     "0.33"
			global CONTVAR `cv'
			global MODE    "power"
			global PIN     0.8
			global TAG     "C2-`s'-`pname'-`cvname'-N"
			global CMD     "slopepower sdmt, subject(id) time(vdate) obs casecon(case) scale(365) schedule(`s') `dopt' effectiveness(0.33) alpha(0.05) power(0.8) `cvopt'"
			sprun

			* solve for power
			spreset
			global DS      "slpower2"
			global MODEL   "obs_cases"
			global SCALE   365
			global SCHED   "`s'"
			global DROPS   "`d'"
			global EFF     "0.33"
			global CONTVAR `cv'
			global MODE    "n"
			global NIN     200
			global TAG     "C2-`s'-`pname'-`cvname'-P"
			global CMD     "slopepower sdmt, subject(id) time(vdate) obs casecon(case) scale(365) schedule(`s') `dopt' effectiveness(0.33) alpha(0.05) n(200) `cvopt'"
			sprun
		}
	}
}

* ---------------------------------------------------------------------------
* slpower3 --- model 3, randomised trial data
* ---------------------------------------------------------------------------
*
* Three targets per design:
*   eff033 / eff1  measure toward zero, ignoring the treated arm
*   usetrt         measure toward the observed treated slope, effectiveness 1
*
* eff1 and usetrt are the informative pair: they differ ONLY in the reference
* slope, so any disagreement between Stata and R isolates cleanly to the
* reference-slope rule rather than to the scaling.

use "`data'/slpower3.dta", clear

local scheds3 `" "1 2" "2 3" "1 2 3" "1 2 5" "'

foreach s of local scheds3 {

	local K : word count `s'

	forvalues p = 1/2 {

		local d ""
		local pname "none"
		if `p' == 2 {
			local pname "front"
			local d "0.2"
			forvalues i = 2/`K' {
				local d "`d' 0.1"
			}
			local d = trim("`d'")
		}

		local dopt ""
		if "`d'" != "" local dopt "dropouts(`d')"

		forvalues t = 1/3 {

			if `t' == 1 {
				local topt "effectiveness(0.33)"
				local tname "eff033"
				local effval "0.33"
				local ut 0
			}
			if `t' == 2 {
				local topt "effectiveness(1)"
				local tname "eff1"
				local effval "1"
				local ut 0
			}
			if `t' == 3 {
				local topt "usetrt"
				local tname "usetrt"
				local effval ""
				local ut 1
			}

			* solve for N
			spreset
			global DS     "slpower3"
			global MODEL  "rct"
			global SCHED  "`s'"
			global DROPS  "`d'"
			global EFF    "`effval'"
			global USETRT `ut'
			global MODE   "power"
			global PIN    0.8
			global TAG    "C3-`s'-`pname'-`tname'-N"
			global CMD    "slopepower sdmt, subject(id) time(visit) rct treat(treat) schedule(`s') `dopt' `topt' alpha(0.05) power(0.8)"
			sprun

			* solve for power
			spreset
			global DS     "slpower3"
			global MODEL  "rct"
			global SCHED  "`s'"
			global DROPS  "`d'"
			global EFF    "`effval'"
			global USETRT `ut'
			global MODE   "n"
			global NIN    318
			global TAG    "C3-`s'-`pname'-`tname'-P"
			global CMD    "slopepower sdmt, subject(id) time(visit) rct treat(treat) schedule(`s') `dopt' `topt' alpha(0.05) n(318)"
			sprun
		}
	}
}

spclose
