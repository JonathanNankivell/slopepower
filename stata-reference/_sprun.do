*! Shared harness for generating golden-reference tables from slopepower.ado
*!
*! Every grid do-file in this directory opens a postfile called H with the
*! column layout below, sets a handful of globals describing one run, and calls
*! sprun.  The globals carry the *inputs* (so the R side can reconstruct the
*! design without parsing Stata syntax) and r(table) carries the *outputs*.
*!
*! Column layout, fixed for all files:
*!
*!   tag      short label for the run, unique within a file
*!   dataset  slpower1 | slpower2 | slpower3
*!   model    obs_nocont | obs_cases | rct   (Stata models 2, 1, 3)
*!   scale    the scale() argument
*!   subset   the if-condition, "" for none
*!   sched    the schedule() numlist, space separated
*!   drops    the dropouts() numlist, "" for none
*!   effin    the effectiveness() argument, "" when omitted or usetrt
*!   usetrt   1 if usetrt was given, else 0
*!   contvar  1 if nocontvar was given, else 0
*!   alpha    the alpha() argument
*!   mode     "power" (solve for N) or "n" (solve for power)
*!   pin      the power() argument, . when solving for power
*!   nin      the n() argument, . when solving for N
*!   rc       _rc from the call: 0 on success, the Stata error code otherwise
*!   o_*      the ten columns of r(table), plus the two count scalars
*!
*! r(table) columns are, in order:
*!   1 alpha  2 power  3 N  4 N1  5 N2  6 effectiveness
*!   7 tte    8 var_tte  9 slope_cases/untreated  10 slope_controls/treated
*! Column 10 is missing by construction for model 2, and for model 3 without
*! usetrt, even though a comparator slope was in fact estimated.

capture program drop sprun
program define sprun

	capture noisily $CMD

	local rc = _rc

	local o1 = .
	local o2 = .
	local o3 = .
	local o4 = .
	local o5 = .
	local o6 = .
	local o7 = .
	local o8 = .
	local o9 = .
	local o10 = .
	local onobs = .
	local onsub = .

	if `rc' == 0 {
		tempname R
		matrix `R' = r(table)
		local o1  = `R'[1,1]
		local o2  = `R'[1,2]
		local o3  = `R'[1,3]
		local o4  = `R'[1,4]
		local o5  = `R'[1,5]
		local o6  = `R'[1,6]
		local o7  = `R'[1,7]
		local o8  = `R'[1,8]
		local o9  = `R'[1,9]
		local o10 = `R'[1,10]
		local onobs = r(obs_in_model)
		local onsub = r(subjects_in_model)
	}

	post H ("$TAG") ("$DS") ("$MODEL") ($SCALE) ("$SUBSET") ("$SCHED") ///
	       ("$DROPS") ("$EFF") ($USETRT) ($CONTVAR) ($ALPHA) ("$MODE") ///
	       ($PIN) ($NIN) (`rc') ///
	       (`o1') (`o2') (`o3') (`o4') (`o5') (`o6') (`o7') (`o8') ///
	       (`o9') (`o10') (`onobs') (`onsub')

	display as text "  [`rc'] $TAG"

end

* Open the postfile.  Call as:  spopen "name-of-output"
capture program drop spopen
program define spopen
	args stem
	global SPSTEM "`stem'"
	postfile H str40 tag str10 dataset str12 model double scale str24 subset ///
	         str40 sched str64 drops str8 effin byte usetrt byte contvar ///
	         double alpha str8 mode double pin double nin int rc ///
	         double o_alpha double o_power double o_n double o_n1 double o_n2 ///
	         double o_effectiveness double o_tte double o_var_tte ///
	         double o_slope double o_slope_comp double o_obs double o_subjects ///
	         using "`stem'.dta", replace
end

* Close it and write the CSV the R suite reads.
capture program drop spclose
program define spclose
	postclose H
	preserve
		use "${SPSTEM}.dta", clear
		export delimited using "${SPSTEM}.csv", replace
		display as text "wrote ${SPSTEM}.csv, " _N " rows"
	restore
end

* Reset the run description to neutral defaults before each call.
capture program drop spreset
program define spreset
	global TAG     ""
	global DS      ""
	global MODEL   ""
	global SCALE   1
	global SUBSET  ""
	global SCHED   ""
	global DROPS   ""
	global EFF     ""
	global USETRT  0
	global CONTVAR 0
	global ALPHA   0.05
	global MODE    "power"
	global PIN     .
	global NIN     .
	global CMD     ""
end
