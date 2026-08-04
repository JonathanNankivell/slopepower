{smcl}
{* *! version 2.1)}{...}
{vieweralsosee "[R] mixed" "help mixed"}{...}
{viewerjumpto "Title" "slopepower##title"}{...}
{viewerjumpto "Syntax" "slopepower##syntax"}{...}
{viewerjumpto "Description" "slopepower##description"}{...}
{viewerjumpto "Options" "slopepower##options"}{...}
{viewerjumpto "Examples" "slopepower##examples"}{...}
{viewerjumpto "Stored results" "slopepower##results"}{...}
{viewerjumpto "Authors" "slopepower##authors"}{...}
{viewerjumpto "Also see" "slopepower##alsosee"}{...}
{cmd:help slopepower}{right: ({browse "https://doi.org/10.1177/1536867X211045512":SJ21-3: st0647})}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 19 21 2}{...}
{p2col :{cmd:slopepower} {hline 2}}Sample-size and power calculator for
outcomes analyzed using a slope (that is, repeated measures across multiple
timepoints); a linear mixed model is used on data in memory to obtain
estimates for slopes and variances among people with (and possibly without)
the condition of interest{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:slopepower}
{it:{help depvar}}
{ifin}{cmd:,}
{opth subj:ect(varname)}
{opth tim:e(varname)}
{opth sched:ule(numlist)}
{c -(}{cmd:obs}|{cmd:rct}{c )-}
[{it:{help slopepower##options_table:options}}]

{synoptset 20 tabbed}{...}
{marker options_table}{...}
{synopthdr}
{synoptline}
{syntab :Options for data in memory}
{p2coldent :* {opth subj:ect(varname)}}variable defining each subject{p_end}
{p2coldent :* {opth tim:e(varname)}}variable defining the time of each measurement{p_end}
{p2coldent :+ {opt obs}}data in memory are from an observational study{p_end}
{p2coldent :+ {opt rct}}data in memory are from a randomized controlled trial
(RCT){p_end}
{synopt :{opt nocont:rols}}observational data in memory contain no healthy
controls; use only with option {cmd:obs}{p_end}
{synopt :{opth case:con(varname)}}variable defining if a person is a case or
healthy control; required if observational data are used{p_end}
{synopt :{opth tr:eat(varname)}}variable defining if a person received the
intervention; required if RCT data are used{p_end}

{syntab :Options for planned trial}
{p2coldent :* {opth sched:ule(numlist)}}visit times for the proposed study{p_end}
{synopt :{opth drop:outs(numlist)}}proportion of dropouts at each visit; must correspond to the schedule list{p_end}
{synopt :{opt sca:le(#)}}ratio between the time and schedule timescales{p_end}
{synopt :{opt a:lpha(#)}}significance level; default is {cmd:alpha(0.05)}{p_end}
{synopt :{opt pow:er(#)}}power; default is {cmd:power(0.8)}; required to compute sample size{p_end}
{synopt :{opt n(#)}}total sample size; required to compute power{p_end}
{synopt :{opt eff:ectiveness(#)}}target effectiveness of the treatment to be
trialed; default is {cmd:effectiveness(0.25)}{p_end}
{synopt :{opt uset:rt}}use the observed effectiveness of the treatment; use
only with option {opt rct}{p_end}

{syntab :Model options}
{synopt :{opt iter:ate(#)}}maximum number of iterations allowed in the {cmd:mixed} model{p_end}
{synopt :{opt nocontv:ar}}omit the random slope variance and covariance for healthy controls in the {cmd:mixed} model{p_end}
{synoptline}
{pstd}* {cmd:subject()}, {cmd:time()}, and {cmd:schedule()} are
required.{p_end}
{pstd}+ Exactly one of {cmd:obs} or {opt rct} must be specified.
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:slopepower} performs a sample-size or power calculation for a proposed
two-arm parallel group randomized clinical trial where the outcome of interest
is a slope measured over time and the treatment is hoped to alter this slope
to be more similar to the slope of people without the condition.  The
calculations are based partly upon data the user has read into Stata's memory
and partly on user input.  A linear mixed model is run (using {helpb mixed})
on the data in memory to estimate a plausible treatment-effect variance and
the slope in those who are untreated; the remaining parameters are specified
by the user.  The data should come from either an observational study or a
similar clinical trial, and they should contain repeated measurements of the
outcome in long format (see {helpb reshape} for more details).  The data in
memory will not be altered by this command.


{marker options}{...}
{title:Options}

    {title:Options for data in memory}

{phang}
{opth subject(varname)} is the unique identifier for participants in the
user-supplied dataset.  {cmd:subject()} is required.

{phang}
{opt time(varname)} is the time variable of visits in the dataset.  This can
be in any units (for example, days, months, years).  It is assumed to be time
since start of observation for each individual.  If it is not (for instance,
if it is an actual calendar date), {cmd:slopepower} will issue a warning and
rescale it accordingly.  {cmd:time()} is required.

{phang}
{opt obs} and {opt rct} tell Stata the nature of the data in memory.
{cmd:obs} should generally be used for observational data and {cmd:rct} for
previously collected trial data (with an exception mentioned below).  Exactly
one of {cmd:obs} or {opt rct} must be specified.

{phang}
{cmd:nocontrols} should be used with {opt obs} if all the subjects in your
observational data have the condition of interest (that is, if there are no
healthy controls).

{phang}
{opt casecon(varname)} specifies the variable used to identify cases in
observational data; it can be used only with {opt obs}.  It must be a binary
0/1 variable with the cases coded as 1.

{phang}
{opth treat(varname)} specifies the treatment variable when you are using RCT
data; it can be specified only with {opt rct}.  It must be a binary 0/1
variable, with the experimental group coded as 1.

    {title:Options for planned trial}

{phang}
{opth schedule(numlist)} specifies the visit times for the proposed trial.  A
baseline visit at time 0 is assumed; this list should describe subsequent
visits in whole-number units of time.  The default is to use the same time unit
as the time variable in the dataset.  To use a different timescale, specify how
many {opt time()} units make one {opt schedule()} unit in the {cmd:scale()}
option. {cmd:schedule()} is required.

{phang}
{opth dropouts(numlist)} specifies the estimated proportion of dropouts you
expect at each study visit.  It must correspond exactly to the schedule list.
Each number in the list is a proportion between 0 and 1; this is the fraction
of subjects (of those who start the trial) you estimate will fail to attend
that visit.  We follow the pattern-mixture method of Dawson and Lagakos (1991,
1993).

{phang}
{opt scale(#)} specifies the ratio between the time and visit timescales.  For
instance, if the time variable in your dataset is in days and you wish to have
visits annually for three years, you would specify {cmd:scale(365)} and
{cmd:schedule(1 2 3)}.

{phang}
{opt alpha(#)} sets the significance level (also known as type I error rate)
to be used in the planned study.  The default is {cmd:alpha(0.05)}.

{phang}
{opt power(#)} sets the power for the planned study.  The default is
{cmd:power(0.8)}.  This option is required to compute the sample size.

{phang}
{opt n(#)} specifies the total number of participants who will be in the
trial.  If an odd number is given, (n-1) will be used to allow equal numbers
per arm.  This option is required to compute the power.  Only one of
{cmd:power()} or {opt n()} may be specified.

{phang}
{opt effectiveness(#)} and {opt usetrt} specify the effect size you would like
to be able to detect in the future trial.  {cmd:effectiveness()} specifies this effect
size as a proportion of the difference between cases and healthy controls in
the observational data in memory.  If RCT data, or observational data with no
healthy controls, are used, {cmd:effectiveness()} is a proportion of the
difference toward a slope of 0.  This must be a number between 0 and 1; the
default is {cmd:effectiveness(0.25)}.  {cmd:usetrt} specifies that, when RCT
data are used, the planned study is targeting the same effect size as observed
from the previous dataset.  You can specify only one of {cmd:effectiveness()}
or {cmd:usetrt}.

    {title:Model options}

{phang}
{opt iterate(#)} is used as an option in the {cmd:mixed} command, which
specifies the maximum number of iterations allowed in the mixed model.

{phang}
{opt nocontvar} specifies that the mixed model should not estimate a
random-slopes variance parameter or the covariance between random slopes and
intercepts for healthy controls.  This is applicable only when you are using
observational data with healthy controls.  Ignoring this variance and
covariance may help the model to converge.


{marker examples}{...}
{title:Examples}

{pstd}
Load example observational data with no healthy controls ({cmd:slpower1.dta},
available as an ancillary file):{p_end}
{phang2}{bf:{stata "use slpower1":. use slpower1}}{p_end}

{pstd} 
Use {cmd:slopepower} to give the sample size for an RCT with annual visits
over 2 years, assuming no dropouts, with 80% power to detect a treatment
effect that will eliminate one-third of the slope in the outcome {cmd:sdmt}.
The default value of 5% type I error rate is used:{p_end}
{phang2}{bf:{stata "slopepower sdmt, schedule(1 2) subject(id) time(visit) obs nocontrols effectiveness(0.33)":. slopepower sdmt, schedule(1 2) subject(id) time(visit) obs nocontrols effectiveness(0.33)}}{p_end}

{pstd} 
Extend the trial to five years, with no additional interim visits, and
assuming that 10% of participants will be lost to follow-up between the visit
at year two and the final visit:{p_end}
{phang2}{bf:. {stata slopepower sdmt, schedule(1 2 5) subject(id) time(visit) obs nocontrols effectiveness(0.33) dropouts(0 0 0.1)}}{p_end}

{pstd} 
Schedule visits every six months in a two-year trial by using the
{cmd:scale()} option:{p_end}
{phang2}{bf:{stata "slopepower sdmt, schedule(1 2 3 4) scale(0.5) subject(id) time(visit) obs nocontrols effectiveness(0.33)":. slopepower sdmt, schedule(1 2 3 4) scale(0.5) subject(id) time(visit) obs nocontrols effectiveness(0.33)}}{p_end}

{pstd}
Load example observational data with a healthy control group
({cmd:slpower2.dta}, available as an ancillary file):{p_end}
{phang2}{bf:{stata "use slpower2, clear":. use slpower2, clear}}{p_end}

{pstd} 
Use {cmd:slopepower} to give the sample size for an RCT with annual visits
over 2 years, assuming no dropouts, with 80% power to detect a treatment
effect that will eliminate one-third of the slope.  The {cmd:scale()} option
is used to tell Stata that the time variable is a date (recorded in
days):{p_end}
{phang2}{bf:{stata "slopepower sdmt, schedule(1 2) scale(365) subject(id) time(vdate) obs casecon(case) effectiveness(0.33)":. slopepower sdmt, schedule(1 2) scale(365) subject(id) time(vdate) obs casecon(case) effectiveness(0.33)}}{p_end}

{pstd} 
Calculate the power for a 200-person trial, assuming a dropout rate of 5% per
year of those who start the trial:{p_end}
{phang2}{bf:. {stata slopepower sdmt, schedule(1 2) scale(365) subject(id) time(vdate) obs casecon(case) effectiveness(0.33) n(200) dropouts(0.05 0.05)}}{p_end}

{pstd}
Load example RCT data with treated and untreated groups ({cmd:slpower3.dta},
available as an ancillary file):{p_end}
{phang2}{bf:{stata "use slpower3, clear":. use slpower3, clear}}{p_end}

{pstd} 
Use {cmd:slopepower} to give the sample size for a 3-year RCT with an interim
visit at 2 years, with 80% power to detect the same treatment effect as
observed in the previous RCT.  Assume a dropout rate of 10% per year:{p_end}
{phang2}{bf:{stata "slopepower sdmt, schedule(2 3) subject(id) time(visit) rct treat(treat) usetrt dropout(0.2 0.1)":. slopepower sdmt, schedule(2 3) subject(id) time(visit) rct treat(treat) usetrt dropout(0.2 0.1)}}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:slopepower} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 24 2: Scalars}{p_end}
{synopt:{cmd:r(subjects_in_model)}}number of subjects that were used in the mixed model{p_end}
{synopt:{cmd:r(obs_in_model)}}number of observations that were used in the mixed model{p_end}
{synopt:{cmd:r(alpha)}}type I error rate (significance level) of the planned trial{p_end}
{synopt:{cmd:r(power)}}power of the planned trial{p_end}
{synopt:{cmd:r(fupvisits)}}specified number of follow-up visits{p_end}
{synopt:{cmd:r(sampsize)}}total sample size{p_end}
{synopt:{cmd:r(effectiveness)}}specified target effectiveness of the proposed treatment{p_end}
{synopt:{cmd:r(tte)}}target treatment effect{p_end}
{synopt:{cmd:r(var_tte)}}variance of the target treatment effect for a hypothetical future two-person trial{p_end}
{synopt:{cmd:r(slope_cases)}}observed slope in the cases (if calculated){p_end}
{synopt:{cmd:r(slope_controls)}}observed slope in the healthy controls (if calculated){p_end}
{synopt:{cmd:r(slope_untreated)}}observed slope in the control arm of the previous RCT (if calculated){p_end}
{synopt:{cmd:r(slope_treated)}}observed slope in the experimental arm of the previous RCT (if calculated){p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 24 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}table of results{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Dawson, J. D., and S. W. Lagakos. 1991. Analyzing laboratory marker changes in
AIDS clinical trials. {it:Journal of Acquired Immune Deficiency Syndromes} 4:
667-676.

{phang}
------. 1993. Size and power of two-sample tests of repeated measures 
data. {it:Biometrics} 49: 1022-1032.
{browse "https://doi.org/10.2307/2532244"}.


{marker authors}{...}
{title:Authors}

{pstd}
Stephen Nash, Katy E. Morgan, Amy Mulick{break}
London School of Hygiene and Tropical Medicine{break}
London, UK{break}
katy.morgan@lshtm.ac.uk


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 21, number 3: {browse "https://doi.org/10.1177/1536867X211045512":st0647}{p_end}
