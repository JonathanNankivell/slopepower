*! Grid 4 --- boundaries, guards and the places the two implementations were
*! always going to disagree.
*!
*! Rows here are as much about _rc as about the numbers.  The R port turns
*! several of Stata's silent behaviours into errors and one of Stata's errors
*! into an accepted input, and each of those decisions should be recorded
*! against what Stata actually does rather than against what the .ado looks
*! like it does.
*!
*! The rows to read first:
*!
*!   SAT-*   Power saturation with dropout.  When solving for power, Stata
*!           reports var_tte by inverting the sample size formula using the
*!           power it just computed.  Once power hits exactly 1 in double
*!           precision, invnormal(1) is +inf and the reported var_tte collapses
*!           to 0.  The R port uses the algebraically equivalent closed form,
*!           which stays finite.  These rows pin down where Stata tips over.
*!
*!   SUM1-*  A dropout list summing to exactly 1 in decimal.  Stata accumulates
*!           1 - 0.3 - 0.3 - 0.4 and tests < 0; in binary that is -5.6e-17.  If
*!           it errors, the R port's explicit 1e-8 tolerance is a real fix.
*!
*!   ODD-*   n(451) must come back as 450: 2*floor(n/2).  n(3) must come back
*!           as 2.  Both implementations claim to do this; only one has been
*!           checked.
*!
*!   IGN-*   Options Stata warns about and then ignores rather than refusing.
*!           The R port has no equivalent of a silently ignored argument, so
*!           these rows document a deliberate divergence.
*!
*! Run with:  do 04_edge_cases.do

clear all
set more off
version 13.0

do "_sprun.do"

local data ".."

spopen "04_edge_cases"

* A small helper so each row is one readable line.
capture program drop edge
program define edge
	gettoken tag rest : 0
	global TAG "`tag'"
	global CMD "`rest'"
	sprun
end

* ---------------------------------------------------------------------------
* slpower1 --- the fast dataset carries most of the edge cases
* ---------------------------------------------------------------------------

use "`data'/slpower1.dta", clear

local base "slopepower sdmt, subject(id) time(visit) obs nocont"

* --- power saturation, with and without dropout -----------------------------
foreach nn in 450 2000 10000 40000 200000 1000000 {
	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCHED "1 2 3"
	global EFF "0.33"
	global MODE "n"
	global NIN `nn'
	edge "SAT-nodrop-`nn'" `base' schedule(1 2 3) effectiveness(0.33) n(`nn')

	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCHED "1 2 3"
	global DROPS "0.05 0.05 0.05"
	global EFF "0.33"
	global MODE "n"
	global NIN `nn'
	edge "SAT-drop-`nn'" `base' schedule(1 2 3) dropouts(0.05 0.05 0.05) effectiveness(0.33) n(`nn')
}

* --- dropout sums at and over the boundary ----------------------------------
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2 3"
global DROPS "0.3 0.3 0.4"
global EFF "0.33"
global PIN 0.8
edge "SUM1-exact" `base' schedule(1 2 3) dropouts(0.3 0.3 0.4) effectiveness(0.33) power(0.8)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2 3"
global DROPS "0.25 0.25 0.5"
global EFF "0.33"
global PIN 0.8
edge "SUM1-exact-b" `base' schedule(1 2 3) dropouts(0.25 0.25 0.5) effectiveness(0.33) power(0.8)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2 3"
global DROPS "0.4 0.4 0.4"
global EFF "0.33"
global PIN 0.8
edge "SUM1-over" `base' schedule(1 2 3) dropouts(0.4 0.4 0.4) effectiveness(0.33) power(0.8)

* Everyone drops out at the first follow-up: every stratum but the completers
* is the baseline-only stratum, which carries no information, and the completer
* stratum has weight 0.  The effect size should be 0 and N undefined.
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2 3"
global DROPS "1 0 0"
global EFF "0.33"
global PIN 0.8
edge "ZERO-allbaseline" `base' schedule(1 2 3) dropouts(1 0 0) effectiveness(0.33) power(0.8)

* --- n parity and validity --------------------------------------------------
foreach nn in 2 3 4 5 451 999 {
	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCHED "1 2"
	global EFF "0.33"
	global MODE "n"
	global NIN `nn'
	edge "ODD-`nn'" `base' schedule(1 2) effectiveness(0.33) n(`nn')
}

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global EFF "0.33"
global MODE "n"
edge "ODD-fractional" `base' schedule(1 2) effectiveness(0.33) n(450.5)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global EFF "0.33"
global MODE "n"
edge "ODD-one" `base' schedule(1 2) effectiveness(0.33) n(1)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global EFF "0.33"
global MODE "n"
edge "ODD-zero" `base' schedule(1 2) effectiveness(0.33) n(0)

* --- alpha, power and effectiveness at their boundaries ---------------------
foreach a in 0 0.0001 0.5 0.999 1 {
	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCHED "1 2"
	global EFF "0.33"
	global ALPHA `a'
	global PIN 0.8
	edge "ALPHA-`a'" `base' schedule(1 2) effectiveness(0.33) alpha(`a') power(0.8)
}

foreach pw in 0 0.001 0.5 0.999 0.99999 1 {
	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCHED "1 2"
	global EFF "0.33"
	global PIN `pw'
	edge "POW-`pw'" `base' schedule(1 2) effectiveness(0.33) power(`pw')
}

foreach e in 0 0.0001 1 1.0001 2 {
	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCHED "1 2"
	global EFF "`e'"
	global PIN 0.8
	edge "EFF-`e'" `base' schedule(1 2) effectiveness(`e') power(0.8)
}

* --- both / neither of power and n ------------------------------------------
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
global NIN 450
edge "BOTH-power-n" `base' schedule(1 2) effectiveness(0.33) power(0.8) n(450)

* Neither: Stata falls back to power(0.8).  The R port refuses to guess in
* slope_power() and only keeps this default in the slopepower() wrapper.
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global EFF "0.33"
edge "NEITHER-default" `base' schedule(1 2) effectiveness(0.33)

* Nor effectiveness: the documented default is 0.25.
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
edge "DEFAULT-eff" `base' schedule(1 2)

* --- schedule shapes --------------------------------------------------------
*
* A far-out final visit forces Stata to build a (max+1)-square covariance
* matrix and select three rows from it, while the R port builds a 3x3 directly.
* If the two agree at schedule(1 100) they agree on the construction, not just
* on the small cases.
foreach s in "1 100" "1 2 100" "1 50 100" {
	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCHED "`s'"
	global EFF "0.33"
	global PIN 0.8
	edge "SCHED-far-`s'" `base' schedule(`s') effectiveness(0.33) power(0.8)
}

* Rejected by numlist rather than by slopepower: non-integer, non-ascending,
* zero, negative, and a repeated visit.
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "0 1 2"
edge "SCHED-zero" `base' schedule(0 1 2) effectiveness(0.33) power(0.8)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 1.5 2"
edge "SCHED-fractional" `base' schedule(1 1.5 2) effectiveness(0.33) power(0.8)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "3 2 1"
edge "SCHED-descending" `base' schedule(3 2 1) effectiveness(0.33) power(0.8)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2 2"
edge "SCHED-repeat" `base' schedule(1 2 2) effectiveness(0.33) power(0.8)

* --- scale boundaries -------------------------------------------------------
foreach sc in 0 -1 0.001 1000 {
	spreset
	global DS "slpower1"
	global MODEL "obs_nocont"
	global SCALE `sc'
	global SCHED "1 2"
	global EFF "0.33"
	global PIN 0.8
	edge "SCALE-`sc'" `base' scale(`sc') schedule(1 2) effectiveness(0.33) power(0.8)
}

* --- options that are warned about and then ignored -------------------------
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "IGN-nocontvar-model2" `base' schedule(1 2) effectiveness(0.33) power(0.8) nocontvar

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "IGN-usetrt-model2" `base' schedule(1 2) effectiveness(0.33) power(0.8) usetrt

* --- the schedule/dropout length guard --------------------------------------
* Duplicated from 00_open_questions.do so the answer lands in the CSV too.
spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2 3"
global DROPS "0.05 0.05"
global EFF "0.33"
global PIN 0.8
edge "LEN-short" `base' schedule(1 2 3) dropouts(0.05 0.05) effectiveness(0.33) power(0.8)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2"
global DROPS "0.05 0.05 0.05"
global EFF "0.33"
global PIN 0.8
edge "LEN-long" `base' schedule(1 2) dropouts(0.05 0.05 0.05) effectiveness(0.33) power(0.8)

spreset
global DS "slpower1"
global MODEL "obs_nocont"
global SCHED "1 2 3"
global DROPS "-0.1 0.2 0.2"
global EFF "0.33"
global PIN 0.8
edge "LEN-negative" `base' schedule(1 2 3) dropouts(-0.1 0.2 0.2) effectiveness(0.33) power(0.8)

* ---------------------------------------------------------------------------
* slpower3 --- the option conflicts that only exist for rct data
* ---------------------------------------------------------------------------

use "`data'/slpower3.dta", clear

local base3 "slopepower sdmt, subject(id) time(visit) rct treat(treat)"

spreset
global DS "slpower3"
global MODEL "rct"
global SCHED "1 2"
global EFF "0.33"
global USETRT 1
global PIN 0.8
edge "CONF-usetrt-eff" `base3' schedule(1 2) effectiveness(0.33) usetrt power(0.8)

* usetrt alone: Stata sets the reported effectiveness to missing rather than 1.
spreset
global DS "slpower3"
global MODEL "rct"
global SCHED "1 2"
global USETRT 1
global PIN 0.8
edge "CONF-usetrt-alone" `base3' schedule(1 2) usetrt power(0.8)

spreset
global DS "slpower3"
global MODEL "rct"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "IGN-casecon-model3" `base3' schedule(1 2) effectiveness(0.33) power(0.8) casecon(treat)

spreset
global DS "slpower3"
global MODEL "rct"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "IGN-nocontvar-model3" `base3' schedule(1 2) effectiveness(0.33) power(0.8) nocontvar

* rct without treat(), and obs+casecon without casecon(): both must refuse.
spreset
global DS "slpower3"
global MODEL "rct"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "MISS-treat" slopepower sdmt, subject(id) time(visit) rct schedule(1 2) effectiveness(0.33) power(0.8)

spreset
global DS "slpower3"
global MODEL "obs_cases"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "MISS-casecon" slopepower sdmt, subject(id) time(visit) obs schedule(1 2) effectiveness(0.33) power(0.8)

* obs and rct together, and neither.
spreset
global DS "slpower3"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "MISS-both-models" slopepower sdmt, subject(id) time(visit) obs rct treat(treat) schedule(1 2) effectiveness(0.33) power(0.8)

spreset
global DS "slpower3"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "MISS-no-model" slopepower sdmt, subject(id) time(visit) schedule(1 2) effectiveness(0.33) power(0.8)

* nocont with rct: refused outright rather than warned about.
spreset
global DS "slpower3"
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "MISS-rct-nocont" slopepower sdmt, subject(id) time(visit) rct nocont treat(treat) schedule(1 2) effectiveness(0.33) power(0.8)

* ---------------------------------------------------------------------------
* slpower2 --- a group variable that is not 0/1
* ---------------------------------------------------------------------------

use "`data'/slpower2.dta", clear
gen byte case12 = case + 1
gen byte case3  = mod(id, 3)

spreset
global DS "slpower2"
global MODEL "obs_cases"
global SCALE 365
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "GRP-12" slopepower sdmt, subject(id) time(vdate) obs casecon(case12) scale(365) schedule(1 2) effectiveness(0.33) power(0.8)

spreset
global DS "slpower2"
global MODEL "obs_cases"
global SCALE 365
global SCHED "1 2"
global EFF "0.33"
global PIN 0.8
edge "GRP-3levels" slopepower sdmt, subject(id) time(vdate) obs casecon(case3) scale(365) schedule(1 2) effectiveness(0.33) power(0.8)

spclose
