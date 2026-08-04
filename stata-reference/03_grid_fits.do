*! Grid 3 --- everything that changes the fitted parameters.
*!
*! Grids 1 and 2 hold the fit fixed and vary the trial design.  This one does
*! the opposite.  It matters most, because the fit is where the port genuinely
*! diverges: Stata's -mixed- against nlme::lme, Stata's log-SD/atanh
*! parameterisation against nlme's, and Stata's scale() against the R port's
*! decision to drop scale() entirely and build Sigma at arbitrary real times.
*!
*! Three levers, all of which change the estimated variance components:
*!
*!   scale()   Stata divides time by scale BEFORE fitting, then works on a unit
*!             integer grid.  So scale(s) with schedule(j1 ... jK) is the only
*!             way Stata can express a visit at a non-integer time, and it is
*!             the exact Stata counterpart of the R port's
*!                 visits = s * c(0, j1, ..., jK)
*!             on the unscaled time axis.  In exact arithmetic the two agree;
*!             the interesting question is whether they agree numerically once
*!             two different REML optimisers have been through them, and
*!             whether the agreement survives s < 1 (where the slope variance
*!             is inflated by 1/s^2) as well as s > 1.
*!
*!   if        A subset is a free extra dataset.  Each one gives a fresh set of
*!             variance components to push through the whole pipeline, at no
*!             cost in new data.  Subsets that drop the baseline visit also
*!             exercise the per-subject re-origining of time, which slopepower
*!             does silently apart from a warning, and which the R port exposes
*!             as an argument.
*!
*!   nocont / casecon / rct on the SAME data
*!             slpower2 and slpower3 both have a group variable that can be
*!             ignored.  Fitting slpower2 with nocont pools cases and controls
*!             into one group; fitting slpower3 with nocont pools the arms.
*!             Neither is a sensible analysis, and both are a clean check that
*!             the model-2 branch is reached identically from different data.
*!
*! Run with:  do 03_grid_fits.do

clear all
set more off
version 13.0

do "_sprun.do"

local data ".."

spopen "03_grid_fits"

* ---------------------------------------------------------------------------
* slpower1 --- scale sweep
* ---------------------------------------------------------------------------
*
* visit runs 0..3, so scale(s) with schedule(1 ... K) reaches visit times
* s, 2s, ... Ks.  The pairs below put visits at halves, quarters, fifths,
* doubles and triples of a year.  The R side must reproduce each with
* visits = s * c(0, schedule) and NO rescaling of its own time column, and
* again with the time column divided by s and integer visits; contract section
* 7 requires both routes to agree.

use "`data'/slpower1.dta", clear

local n1 = 0
foreach sc in 0.2 0.25 0.5 1 2 3 {
	local scheds `" "1" "1 2" "1 2 3" "1 3" "1 2 3 4" "'
	foreach s of local scheds {
		spreset
		global DS     "slpower1"
		global MODEL  "obs_nocont"
		global SCALE  `sc'
		global SCHED  "`s'"
		global EFF    "0.33"
		global MODE   "power"
		global PIN    0.8
		global TAG    "S1-sc`sc'-`s'"
		global CMD    "slopepower sdmt, subject(id) time(visit) obs nocont scale(`sc') schedule(`s') effectiveness(0.33) alpha(0.05) power(0.8)"
		sprun
		local n1 = `n1' + 1
	}
}
display as text _n "slpower1 scale sweep: `n1' runs"

* ---------------------------------------------------------------------------
* slpower1 --- subsets
* ---------------------------------------------------------------------------
*
* if visit>0 is the important one: it removes every subject's baseline, so
* slopepower shifts each subject's first visit to zero and prints a warning.
* The remaining design is three equally spaced visits, but the intercept
* variance is now the variance at what used to be year 1.

* Written so that every condition is valid R as well as valid Stata --- the R
* harness evaluates these strings verbatim against the same columns.
local subs `" "id<=50" "id<=100" "id>100" "visit<3" "visit>0" "visit!=1" "visit==0 | visit==3" "'

local n2 = 0
foreach sub of local subs {
	foreach s in "1" "1 2" {
		spreset
		global DS     "slpower1"
		global MODEL  "obs_nocont"
		global SUBSET "`sub'"
		global SCHED  "`s'"
		global EFF    "0.33"
		global MODE   "power"
		global PIN    0.8
		global TAG    "U1-`sub'-`s'"
		global CMD    "slopepower sdmt if `sub', subject(id) time(visit) obs nocont schedule(`s') effectiveness(0.33) alpha(0.05) power(0.8)"
		sprun
		local n2 = `n2' + 1
	}
}
display as text _n "slpower1 subsets: `n2' runs"

* ---------------------------------------------------------------------------
* slpower2 --- scale sweep and the pooled-group fit
* ---------------------------------------------------------------------------
*
* vdate is in days.  scale(365) is the paper's year axis; 730 makes the unit
* two years; 182.5 and 91.25 put the visits at half- and quarter-year spacing,
* which is where the R port's real-valued Sigma has to earn its keep.
*
* Note that scale here is doing double duty --- it converts the units AND sets
* the grid --- so an error in either shows up as the same disagreement.  The
* subset rows disambiguate: they change the fit without touching scale.

use "`data'/slpower2.dta", clear

local n3 = 0
foreach sc in 91.25 182.5 365 730 {
	foreach s in "1" "1 2" "1 2 3" {
		spreset
		global DS     "slpower2"
		global MODEL  "obs_cases"
		global SCALE  `sc'
		global SCHED  "`s'"
		global EFF    "0.33"
		global MODE   "power"
		global PIN    0.8
		global TAG    "S2-sc`sc'-`s'"
		global CMD    "slopepower sdmt, subject(id) time(vdate) obs casecon(case) scale(`sc') schedule(`s') effectiveness(0.33) alpha(0.05) power(0.8)"
		sprun
		local n3 = `n3' + 1
	}
}

* The same data through the model-2 branch: cases and controls pooled.  Also
* the cases alone, which is the analysis the single-group model is actually
* for, and which should land somewhere between the pooled fit and slpower1.
foreach spec in "" "case==1" "case==0" {
	local subname "`spec'"
	if "`spec'" == "" local subname "all"
	local ifpart ""
	if "`spec'" != "" local ifpart "if `spec'"
	foreach s in "1 2" "1 2 3" {
		spreset
		global DS     "slpower2"
		global MODEL  "obs_nocont"
		global SCALE  365
		global SUBSET "`spec'"
		global SCHED  "`s'"
		global EFF    "0.33"
		global MODE   "power"
		global PIN    0.8
		global TAG    "P2-`subname'-`s'"
		global CMD    "slopepower sdmt `ifpart', subject(id) time(vdate) obs nocont scale(365) schedule(`s') effectiveness(0.33) alpha(0.05) power(0.8)"
		sprun
		local n3 = `n3' + 1
	}
}
display as text _n "slpower2: `n3' runs"

* ---------------------------------------------------------------------------
* slpower3 --- scale sweep, subsets, and the pooled-arm fit
* ---------------------------------------------------------------------------
*
* visit takes only 0, 0.5 and 2, so this is the dataset whose observed visit
* times are already irregular and non-integer.  scale(0.5) puts the fitted time
* unit at half a year, which is the finest unit that lands every observed visit
* on an integer; scale(0.25) subdivides further.  "if visit<2" leaves only two
* timepoints per subject, which may or may not support a random-slope model ---
* a non-zero rc there is a result, not a failure of the harness.

use "`data'/slpower3.dta", clear

local n4 = 0
foreach sc in 0.25 0.5 1 2 {
	foreach s in "1 2" "1 2 4" {
		spreset
		global DS     "slpower3"
		global MODEL  "rct"
		global SCALE  `sc'
		global SCHED  "`s'"
		global EFF    "0.33"
		global MODE   "power"
		global PIN    0.8
		global TAG    "S3-sc`sc'-`s'"
		global CMD    "slopepower sdmt, subject(id) time(visit) rct treat(treat) scale(`sc') schedule(`s') effectiveness(0.33) alpha(0.05) power(0.8)"
		sprun
		local n4 = `n4' + 1
	}
}

* treat is assigned by id block in slpower3 (1-75 control, 76-150 treated), so
* any subset of the form id<=k leaves a single arm and slopepower refuses it
* with "Treatment variable must have exactly two levels". These straddle the
* boundary instead, and stay balanced.
local subs3 `" "id<=40 | id>110" "id>40 & id<=110" "visit<2" "visit>0" "'
foreach sub of local subs3 {
	spreset
	global DS     "slpower3"
	global MODEL  "rct"
	global SUBSET "`sub'"
	global SCHED  "1 2"
	global EFF    "0.33"
	global MODE   "power"
	global PIN    0.8
	global TAG    "U3-`sub'"
	global CMD    "slopepower sdmt if `sub', subject(id) time(visit) rct treat(treat) schedule(1 2) effectiveness(0.33) alpha(0.05) power(0.8)"
	sprun
	local n4 = `n4' + 1
}

* Arms pooled through the model-2 branch, and each arm alone.
foreach spec in "" "treat==1" "treat==0" {
	local subname "`spec'"
	if "`spec'" == "" local subname "all"
	local ifpart ""
	if "`spec'" != "" local ifpart "if `spec'"
	spreset
	global DS     "slpower3"
	global MODEL  "obs_nocont"
	global SUBSET "`spec'"
	global SCHED  "1 2"
	global EFF    "0.33"
	global MODE   "power"
	global PIN    0.8
	global TAG    "P3-`subname'"
	global CMD    "slopepower sdmt `ifpart', subject(id) time(visit) obs nocont schedule(1 2) effectiveness(0.33) alpha(0.05) power(0.8)"
	sprun
	local n4 = `n4' + 1
}
display as text _n "slpower3: `n4' runs"

spclose
