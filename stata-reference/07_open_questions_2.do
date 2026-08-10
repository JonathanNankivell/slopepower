*! Two more questions about slopepower.ado that only a Stata licence can settle.
*!
*! BOTH ARE NOW ANSWERED (2026-08-10), and like questions 1 and 2 in 00, both
*! came back in the .ado's favour. See the per-question notes below. The file
*! is kept because it is the evidence, and it re-checks in a few seconds.
*!
*!   3. The dropout-total guard does NOT reject a legal list. Stata's locals
*!      round-trip through a decimal string, so the subtractive residue is
*!      exactly 0, not -5.551e-17. dropouts(0.3 0.3 0.4) is accepted and gives
*!      N = 1418, which is what the R port gives. No divergence.
*!   4. The unbalanced macro quote at :196 prints the backtick literally and
*!      carries on: _rc = 0, N = 712, same as the well-formed sibling at :132.
*!      Answer (a) below. The R port's warn-and-continue is faithful.
*!
*! Run with:  do 07_open_questions_2.do
*!
*! Question 4c deliberately runs a malformed line. It is written to a one-line
*! do-file and executed from there precisely so that, if the unbalanced macro
*! quote swallows what follows it, the damage is confined to that file rather
*! than to the rest of this one. See the note in README.md about the first run
*! of 00 dying on a stray backtick -- and the note at Q4c itself, which did
*! exactly that on ITS first run, for a reason worth reading.

clear all
set more off
version 13.0

local data ".."          // where slpower*.dta live, relative to this directory

capture log close _all
log using "07_open_questions_2.log", replace text

* ---------------------------------------------------------------------------
* QUESTION 3 --- does the dropout-total guard reject a legal dropout list?
*
* slopepower.ado:259-271 checks "dropouts cannot exceed 100%" by accumulating
* by repeated subtraction rather than by summing:
*
*     local dfrac_complete = 1
*     foreach fr of numlist `dropouts' {
*         ...
*         local dfrac_complete = `dfrac_complete' - `fr'
*     }
*     if `dfrac_complete' < 0 {
*         display as error _n "Dropouts cannot exceed 100%"
*         exit 198
*     }
*
* In IEEE double precision 1 - 0.3 - 0.3 - 0.4 is -5.551e-17, not 0. If Stata's
* locals carry the full double through each step, dropouts(0.3 0.3 0.4) -- a
* perfectly legal list, summing to exactly 1 in decimal, meaning "nobody
* completes" -- trips that guard and is rejected. If Stata instead round-trips
* each intermediate through a shorter string representation (1 - 0.3 stored as
* "0.7", then 0.7 - 0.3 stored as "0.4"), the residue is exactly 0 and there is
* no defect here at all.
*
* This matters because the R port sums once and compares against 1 with a 1e-8
* tolerance, and DIVERGENCES.md claimed that as a divergence on the strength of
* reading the source, not of running it.
*
* ANSWERED, 2026-08-10: there is no defect. Stata ACCEPTS dropouts(0.3 0.3 0.4)
* and returns N = 1418, the same as the R port.
*
* Q3d below is the mechanism, and it is the interesting part. A Stata local
* does not carry the full double from step to step -- it round-trips through a
* decimal string representation:
*
*     1   - 0.3  stored as ".7"
*     0.7 - 0.3  stored as ".4"    <- the double here is 0.39999999999999997
*     0.4 - 0.4  stored as "0"
*
* so the accumulated residue is exactly 0, while the same arithmetic in one
* expression is -5.551e-17. The two %21x lines printed by Q3d differ for that
* reason: +0.0000000000000X-3ff against -1.0000000000000X-036. The rounding
* that looks like sloppiness is what makes the guard safe.
*
* So the R port's tolerance restates Stata rather than repairing it, and
* DIVERGENCES.md now records this under "Claims checked and rejected".
*
* What the R port does with the same three calls, for comparison:
*
*     dropouts(0.3 0.3 0.4)   accepted, N = 1418   (var_tte 27.51684853)
*     dropouts(0.1 0.2 0.7)   accepted, N =  928   (var_tte 18.00820553)
*     dropouts(0.5 0.3 0.3)   rejected: sums to 1.1
*
* Note the second row. 1 - 0.1 - 0.2 - 0.7 is exactly 0 in double precision, so if
* Stata rejects 0.3 0.3 0.4 but accepts 0.1 0.2 0.7, the guard is not uniformly
* broken -- it fails on some legal lists and not others, which is worse to
* document and better to know.
* ---------------------------------------------------------------------------

use "`data'/slpower1.dta", clear

display as text _n "{hline 70}"
display as text "Q3a: dropouts(0.3 0.3 0.4) --- legal, sums to exactly 1 in decimal"
display as text "{hline 70}"
capture noisily slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2 3) dropouts(0.3 0.3 0.4) effectiveness(0.33) power(0.8)
display as text "Q3a _rc = " _rc
display as text "(0 with N = 1418 => no defect, the R port merely restates Stata."
display as text " 198 'Dropouts cannot exceed 100%' => the divergence is real.)"

display as text _n "{hline 70}"
display as text "Q3b: dropouts(0.1 0.2 0.7) --- also sums to 1, but residue is 0"
display as text "{hline 70}"
capture noisily slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2 3) dropouts(0.1 0.2 0.7) effectiveness(0.33) power(0.8)
display as text "Q3b _rc = " _rc
display as text "(Expect 0 with N = 928 on either answer to Q3a.)"

display as text _n "{hline 70}"
display as text "Q3c: dropouts(0.5 0.3 0.3) --- genuinely over 1, must still error"
display as text "{hline 70}"
capture noisily slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2 3) dropouts(0.5 0.3 0.3) effectiveness(0.33) power(0.8)
display as text "Q3c _rc = " _rc
display as text "(Anything other than 198 would mean the guard is dead, not merely"
display as text " fragile, and every dropout figure in the port needs re-checking.)"

* The direct probe, independent of slopepower: does a local round-trip a
* double exactly? This is the whole question behind Q3a, isolated.
display as text _n "{hline 70}"
display as text "Q3d: what does a local actually store after subtraction?"
display as text "{hline 70}"
local step1 = 1 - 0.3
local step2 = `step1' - 0.3
local step3 = `step2' - 0.4
display as text "  1 - 0.3           stored as: " as result "`step1'"
display as text "  that - 0.3        stored as: " as result "`step2'"
display as text "  that - 0.4        stored as: " as result "`step3'"
display as text _n "  as a double, that final residue is:"
display as text "    %18.0g " as result %18.0g `step3'
display as text "    %21x   " as result %21x   `step3'
display as text "  and the guard's own test, residue < 0, is: " ///
                as result cond(`step3' < 0, "TRUE (rejects)", "FALSE (accepts)")
display as text _n "  for reference, the same three steps in one expression:"
display as text "    %21x   " as result %21x   (1 - 0.3 - 0.3 - 0.4)
display as text "(If the two %21x lines differ, Stata is losing precision in the"
display as text " macro round-trip, which is what would make the guard survive.)"

* Same probe for the 0.1 0.2 0.7 ordering, which should land on exactly 0.
local alt = 1 - 0.1
local alt = `alt' - 0.2
local alt = `alt' - 0.7
display as text _n "  0.1 0.2 0.7 residue: " as result %21x `alt' ///
                as text "   < 0 is " as result cond(`alt' < 0, "TRUE", "FALSE")

* ---------------------------------------------------------------------------
* QUESTION 4 --- what does the unbalanced macro quote at :196 do?
*
* slopepower.ado:195-197 is the warn-and-continue guard for treat() supplied
* with observational data:
*
*     if (`model'!=3) & ("`treat'"!="") {
*         dis as error "Treatment can only be specified with RCT data<BT>. Treatment variable will be ignored."
*         }
*
* where <BT> is a literal backtick -- written here as a placeholder because
* this file must not itself contain an unmatched one, for the same reason 00
* keeps backticks out of its display strings. It opens a macro reference that
* is never closed. The three sibling guards (:132 casecon, :164 usetrt, :232
* nocontvar) are all well formed and all warn-and-continue.
*
* Three possible answers, and they are not equally comfortable:
*
*   (a) Stata prints something mangled and carries on. The R port's warn-and-
*       continue is then a faithful reading of the author's intent, and
*       DIVERGENCES.md section 19 is right as written.
*   (b) Stata aborts with a non-zero _rc. Then supplying treat() with obs data
*       is an error in the original and a warning in the port -- a real
*       behavioural divergence that section 19 currently only speculates about.
*   (c) The unmatched backtick consumes the following lines of the ado-file,
*       in which case the guard at :198 and everything after it is affected
*       too, and the blast radius needs mapping before anything is claimed.
*
* The R port warns and continues, returning the model-2 answer: for the call
* below, N = 712 (the paper's p.588 result, which does not involve treat() at
* all -- the point is that the ignored option leaves it unchanged).
*
* ANSWERED, 2026-08-10: answer (a). Stata prints the stray backtick as a
* literal character --
*
*     Treatment can only be specified with RCT data<BT>. Treatment variable will be ignored.
*
* -- and carries straight on: _rc = 0, N = 712, identical to the well-formed
* sibling in Q4b. So the malformed string is cosmetic, nothing downstream is
* consumed, and the R port's warn-and-continue is a faithful reading rather
* than a repair.
* ---------------------------------------------------------------------------

use "`data'/slpower1.dta", clear
gen byte fake_treat = mod(id, 2)      // any varname will do; :196 fires before
gen byte fake_cc    = mod(id, 2)      // the variable's contents are looked at

display as text _n "{hline 70}"
display as text "Q4a: obs nocont with treat() --- exercises :196 in situ"
display as text "{hline 70}"
capture noisily slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2) treat(fake_treat) effectiveness(0.33) power(0.8)
display as text "Q4a _rc = " _rc
display as text "(0 with N = 712 => answer (a), warns and continues."
display as text " non-zero    => answer (b), the malformed string aborts.)"

* The well-formed sibling at :132, for contrast: same shape of mistake by the
* user, same intent by the author, but a string with no stray backtick.
display as text _n "{hline 70}"
display as text "Q4b: obs nocont with casecon() --- the well-formed sibling, :132"
display as text "{hline 70}"
capture noisily slopepower sdmt, subject(id) time(visit) obs nocont ///
	schedule(1 2) casecon(fake_cc) effectiveness(0.33) power(0.8)
display as text "Q4b _rc = " _rc
display as text "(Expect 0 with N = 712 and a WARNING line. This is what :196"
display as text " is trying to do; Q4a says whether it manages it.)"

* The isolated probe: the offending line on its own, followed by a sentinel.
* Written to a one-line do-file at run time so that THIS file never itself
* contains an unmatched backtick.
*
* Q4a above has already settled the question -- the command completed and
* returned N = 712, so nothing downstream was consumed and answer (c) is out.
* This is kept as the direct demonstration, and as a warning.
*
* On its first run it was written as
*
*     local bt = char(96)
*     file write `fh' `"... RCT data`bt'. Treatment ..."' _n
*
* on the theory that macro expansion is a single left-to-right pass, so a
* backtick substituted INTO the line would land as inert text. It does not:
* Stata rejected the line with "invalid syntax", r(198), which aborted the
* do-file before `log close' and left a zero-byte probe behind. The expanded
* backtick is re-exposed to the parser.
*
* The fix is to keep the character out of the command text altogether.
* `file write' evaluates a parenthesised argument as an expression, so
* (char(96)) writes a backtick without one ever appearing in the line Stata
* has to parse. Do not "simplify" this back into the string.
display as text _n "{hline 70}"
display as text "Q4c: the offending line alone, in a throwaway do-file"
display as text "{hline 70}"

* Named, not a tempfile: `do' assumes a .do suffix when the name it is given
* has none, so a tempfile -- which never has one -- would send it looking for a
* file that does not exist. Erased at the end instead.
local probe "_q4c_probe.do"
tempname fh
capture erase "`probe'"
file open `fh' using "`probe'", write text replace
file write `fh' `"display as error "Treatment can only be specified with RCT data"' ///
                (char(96)) `". Treatment variable will be ignored.""' _n
file write `fh' `"display as text "SENTINEL: the line after the malformed one still ran.""' _n
file close `fh'

display as text "--- contents written to the probe file ---"
type "`probe'"
display as text "--- running it ---"
capture noisily do "`probe'"
display as text "Q4c _rc = " _rc
display as text "(SENTINEL printed => the backtick did not eat the next line,"
display as text " so answer (c) is ruled out and the blast radius is one line.)"
capture erase "`probe'"

log close
