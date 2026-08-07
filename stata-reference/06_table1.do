*! Grid 6 --- the paper's table 1, at full precision.
*!
*! Table 1 (p.595) is the only published result the port checks against a
*! *rounded* figure and nothing else: the paper prints nine powers as
*! percentages to one decimal place, so test-paper-parity.R can only assert
*! them to +/- 5e-4.  Three of the nine already appear in grid 1 by accident of
*! its sweep (schedule 3 and 1 2 3, no dropout and flat 5%), but the other six
*! do not, and they are the interesting ones:
*!
*!   - the three six-month cells need scale(0.5), and grid 1 is entirely at
*!     scale 1.  This is the one design in the paper where Stata's integer
*!     schedule and the R port's real visit times are forced to meet through
*!     scale(), so it is exactly where the two could quietly disagree.
*!   - the two "baseline and final visit only" dropout cells use a
*!     single-element dropouts() list, dropouts(0.15) and dropouts(0.3).  A
*!     one-visit schedule with one dropout stratum is the smallest case the
*!     Dawson-Lagakos weighting has to handle, and nothing else generates it at
*!     these magnitudes.
*!
*! The nine calls in block A are transcribed from the appendix of the paper
*! verbatim (p.601, "commands to produce table 1"), including the conversion of
*! "5% per year" into 0.15 over three years for the final-visit-only design and
*! 0.05 per visit for annual visits.  That transcription is itself worth
*! pinning: if the R port's reading of the dropout column is wrong, the printed
*! percentages are the only thing that would have caught it, and only to 3 d.p.
*!
*! Block B re-runs the same nine designs in the other direction --- solve for N
*! at 80% power instead of power at N = 450 --- which the paper never does.  It
*! costs nine more fits of the cheapest dataset and gives every cell a second,
*! independent comparison, on the mode where Stata and the R port define
*! var_tte the same way (CONTRACT.md 5.6), so var_tte is comparable there and
*! not in block A.
*!
*! No new fits are needed on the R side: grid 5 already carries F1-sc1 and
*! F1-sc0.5, which are the two fits every row here uses, so the tolerance-free
*! test joins onto them and pins all nine powers to machine precision.
*!
*! Run with:  do 06_table1.do

clear all
set more off
version 13.0

do "_sprun.do"

local data ".."

use "`data'/slpower1.dta", clear

spopen "06_table1"

* ---------------------------------------------------------------------------
* The three planned trial designs, all three years long.
* ---------------------------------------------------------------------------
*
* Stata's schedule() takes ascending positive integers and multiplies them by
* scale() to get times, with baseline at 0 implied.  So the six-month design is
* schedule(1 2 3 4 5 6) scale(0.5), which is visits at 0.5, 1, ..., 3 years.
* scale(1) is the default and is passed explicitly on the other two so that
* every row of the CSV records the scale it ran at rather than inheriting it.

local nm1  "final"
local des1 "3"
local sc1  "1"

local nm2  "annual"
local des2 "1 2 3"
local sc2  "1"

local nm3  "sixmonth"
local des3 "1 2 3 4 5 6"
local sc3  "0.5"

* ---------------------------------------------------------------------------
* The three dropout scenarios, per design.
* ---------------------------------------------------------------------------
*
* dropouts() is a list of *incremental* proportions, one per follow-up visit:
* the fraction of the cohort whose last attended visit is the one before.  The
* paper quotes dropout per year and converts, so the three designs get lists of
* different lengths that describe the same annual rate:
*
*   5% per year   final only  0.15 over the whole three years
*                 annual      0.05 at each of three visits
*                 six-month   0.025 at each of six visits
*
* This is the paper's own arithmetic, not a reinterpretation of it: the lists
* below are copied from the appendix.

local lv1   "none"
local dr1_1 ""
local dr2_1 ""
local dr3_1 ""

local lv2   "5pc"
local dr1_2 "0.15"
local dr2_2 "0.05 0.05 0.05"
local dr3_2 "0.025 0.025 0.025 0.025 0.025 0.025"

local lv3   "10pc"
local dr1_3 "0.3"
local dr2_3 "0.1 0.1 0.1"
local dr3_3 "0.05 0.05 0.05 0.05 0.05 0.05"

* ---------------------------------------------------------------------------
* Block A: power at n = 450 --- the nine cells of table 1 as published.
* Block B: N at 80% power --- the same nine designs, the other direction.
* ---------------------------------------------------------------------------

local runs = 0

forvalues i = 1/3 {
	forvalues j = 1/3 {

		local s   "`des`i''"
		local sc  "`sc`i''"
		local nm  "`nm`i''"
		local d   "`dr`i'_`j''"
		local lv  "`lv`j''"

		local dopt ""
		if "`d'" != "" local dopt "dropouts(`d')"

		local opts "subject(id) time(visit) obs nocontrols schedule(`s') scale(`sc') `dopt' effectiveness(0.33) alpha(0.05)"

		* ---- A: n(450), solve for power.  This is the published cell. ----
		spreset
		global DS      "slpower1"
		global MODEL   "obs_nocont"
		global SCALE   `sc'
		global SCHED   "`s'"
		global DROPS   "`d'"
		global EFF     "0.33"
		global ALPHA   0.05
		global MODE    "n"
		global NIN     450
		global TAG     "T1-`nm'-`lv'"
		global CMD     "slopepower sdmt, `opts' n(450)"
		sprun
		local runs = `runs' + 1

		* ---- B: power(0.8), solve for N.  Not in the paper. ----
		spreset
		global DS      "slpower1"
		global MODEL   "obs_nocont"
		global SCALE   `sc'
		global SCHED   "`s'"
		global DROPS   "`d'"
		global EFF     "0.33"
		global ALPHA   0.05
		global MODE    "power"
		global PIN     0.8
		global TAG     "T1N-`nm'-`lv'"
		global CMD     "slopepower sdmt, `opts' power(0.8)"
		sprun
		local runs = `runs' + 1
	}
}

display as text _n "table 1: `runs' runs (expected 18)"

spclose
