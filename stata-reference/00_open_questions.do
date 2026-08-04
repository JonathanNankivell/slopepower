*! Two questions about slopepower.ado that only a Stata licence can settle.
*!
*! Neither produces a grid.  Run this first: question 2 decides whether the
*! model-1 (cases-and-controls) reference numbers in every other file are
*! trustworthy, and it is cheap.
*!
*! Run with:  do 00_open_questions.do
*! Then paste the whole log into the R port's notes.

clear all
set more off
version 13.0

local data ".."          // where slpower*.dta live, relative to this directory

capture log close _all
log using "00_open_questions.log", replace text

* ---------------------------------------------------------------------------
* QUESTION 1 --- is the schedule/dropout length guard vacuous?
*
* slopepower.ado:239 and :276 both read
*
*     local sched_length = `sched_length ' + 1
*
* with a space inside the macro name.  If Stata expands `sched_length ' to
* empty rather than trimming to `sched_length', the expression evaluates as
* " + 1" = 1 on every pass, both counters are stuck at 1, and the guard at :279
*
*     if `sched_length' != `drop_length'  ->  1 != 1  ->  never fires
*
* is dead code.  A dropouts() list of the wrong length would then be accepted
* silently, and the missing strata would be treated as completers --- which
* inflates N without warning.
*
* Test: give a 3-visit schedule a 2-element dropout list.  It must error 198
* with "Dropout list must correspond with visit schedule".  If it prints a
* result instead, the guard is vacuous and the R port's strict length check is
* a genuine fix rather than a restatement.
* ---------------------------------------------------------------------------

use "`data'/slpower1.dta", clear

display as text _n "{hline 70}"
display as text "Q1a: schedule(1 2 3) with dropouts(0.05 0.05) --- SHORT by one"
display as text "{hline 70}"
capture noisily slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2 3) dropouts(0.05 0.05) effectiveness(0.33) power(0.8)
display as text "Q1a _rc = " _rc

display as text _n "{hline 70}"
display as text "Q1b: schedule(1 2) with dropouts(0.05 0.05 0.05) --- LONG by one"
display as text "{hline 70}"
capture noisily slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2) dropouts(0.05 0.05 0.05) effectiveness(0.33) power(0.8)
display as text "Q1b _rc = " _rc

* The direct probe, independent of slopepower: does a macro name with a
* trailing space still expand?
*
* Written without backticks in the display strings: a backtick-doublequote pair
* inside a display argument opens a compound quote and takes the rest of the
* file with it, which is how the first run of this file died at this line and
* left Q2 unreached.
display as text _n "{hline 70}"
display as text "Q1c: does a macro name with a trailing space expand?"
display as text "{hline 70}"
local counter = 7
local probe = `counter ' + 1
display as text "counter = 7, then evaluating counter-with-a-trailing-space + 1 gives " ///
                as result `probe'
display as text "(8 => Stata trims the name, guard is live." _n ///
                " 1 => the name does not expand, guard is vacuous.)"

* For reference: what does the correct 3-element dropout list give?
display as text _n "{hline 70}"
display as text "Q1d: the well-formed call, for comparison"
display as text "{hline 70}"
slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2 3) dropouts(0.05 0.05 0.05) effectiveness(0.33) power(0.8)
return list

* ---------------------------------------------------------------------------
* QUESTION 2 --- how does residuals(independent, by(case)) store its two
* variances?
*
* slopepower.ado:371 computes the CASES' residual variance as
*
*     var_res = (exp(b[1,13] + b[1,14]))^2
*
* which is correct only if mixed stores a base log-SD in b[13] and a
* group-specific OFFSET in b[14].  If instead it stores two independent
* log-SDs, the right answer is exp(b[1,14])^2 and the shipped code is
* multiplying the two standard deviations together.
*
* The R port sidesteps this by extracting variance components by name, and it
* gets 10.354 for cases against 10.699 for controls --- close enough that the
* two readings differ by a factor of ~10.5 in variance, which would be
* impossible to miss in N.  Since the port reproduces the paper's N = 296
* exactly using 10.354, base+offset is the likely answer, but "likely" is not
* "checked".
*
* ANSWERED, 2026-08-04: base+offset, so :371 is CORRECT.
*
* The names are not the ones guessed above. mixed stores the two residual
* parameters as
*     lnsig_e:_cons      the base log-SD, which is the REFERENCE group (cse=0,
*                        the healthy controls)
*     r_lns2ose:_cons    a log-RATIO offset -- the "r_lns" prefix is what makes
*                        base+offset the right reading
* so exp(b13 + b14) is sd_controls * ratio = sd_cases. The arithmetic below
* lands on 10.354254 against mixed's own printed "1: var(e) = 10.35425" for the
* cases, and exp(b13)^2 = 10.698895 against "0: var(e) = 10.6989" for the
* controls. Same story in the nocontvar branch at b[11], b[12].
* ---------------------------------------------------------------------------

use "`data'/slpower2.dta", clear

* Reproduce slopepower's own data preparation exactly (ado lines 305-325).
replace vdate = vdate / 365
drop if missing(sdmt)
tempvar first_time
bysort id (vdate): gen `first_time' = vdate[1]
bysort id: replace vdate = vdate - `first_time'
drop if missing(case)
gen byte cse      = (case != 0)
gen byte ctl      = (case == 0)
gen double tcase  = vdate * cse
gen double tctl   = vdate * ctl

display as text _n "{hline 70}"
display as text "Q2: model-1 mixed fit --- the exact call slopepower makes"
display as text "{hline 70}"

mixed sdmt cse##c.vdate ///
	|| id: tcase cse, cov(uns) nocons ///
	|| id: tctl ctl, cov(uns) nocons ///
	res(ind, by(cse)) reml iterate(16000)

display as text _n "--- e(b), with names: look at the last two columns ---"
matrix b = e(b)
matrix list b, format(%12.8f)
display as text _n "--- column names ---"
local names : colfullnames b
display as text "`names'"

display as text _n "--- estat sd: the variances on the natural scale ---"
capture noisily estat sd

display as text _n "--- the two readings of the cases' residual variance ---"
display as text "  exp(b13 + b14)^2 = " as result %12.6f (exp(b[1,13] + b[1,14]))^2 ///
                as text "   <- what :371 computes"
display as text "  exp(b14)^2       = " as result %12.6f (exp(b[1,14]))^2
display as text "  exp(b13)^2       = " as result %12.6f (exp(b[1,13]))^2
display as text _n "The R port reports 10.354 for cases and 10.699 for controls."
display as text "Whichever line above lands on 10.354 is the parameterisation."

* Same question for the nocontvar branch, which reads b[11] and b[12].
display as text _n "{hline 70}"
display as text "Q2b: the nocontvar branch --- b[11], b[12]"
display as text "{hline 70}"

mixed sdmt cse##c.vdate ///
	|| id: tcase cse, cov(uns) nocons ///
	|| id: ctl, cov(id) nocons ///
	res(ind, by(cse)) reml iterate(16000)

matrix b2 = e(b)
matrix list b2, format(%12.8f)
local names2 : colfullnames b2
display as text "`names2'"
display as text _n "  exp(b11 + b12)^2 = " as result %12.6f (exp(b2[1,11] + b2[1,12]))^2
display as text "  exp(b12)^2       = " as result %12.6f (exp(b2[1,12]))^2

log close
