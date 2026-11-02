*============================
* Alise Ruža || CODING SAMPLE
*============================
* Estimation of first childbirth effect on mothers' earnings ("child penalty") with staggered treatment timing.
* Steps: 1. Data descriptives
* 		 2. Showing the "forbidden comparison" problem of TWFE in a two-period example, and decompose it (Goodman-Bacon)
*		 3. Building a stacked dataset with a clean comparison groups
*		 4. Estimate a stacked event study and compare it with conventional TWFE
*		 5. Express the effect in % of counterfactual earnings (bootstrap SEs)

version 18.5
clear all
set more off

* Required packages (install once from SSC before running):
* ssc install reghdfe   
* ssc install ftools
* ssc install estout
* ssc install palettes
* ssc install colrspace 
* ssc install coefplot
* ssc install bacondecomp
* See README.txt for versions that are used

*========================================================
* Paths: run this file from the project folder. All outputs are written to /output.
*========================================================
global root "`c(pwd)'"
capture confirm file "$root/data/dataset.dta"
if _rc {
    di as error "dataset.dta not found. Set the working directory to the project folder."
    exit 601
}

global root "`c(pwd)'"
global data "$root/data"
global out  "$root/output"
cap mkdir "$out"

cap log close
log using "$out/panel_assignment.log", replace text
cd "$out"
*========================================================
* 1. Uploading the dataset and conducting descriptive analysis
*========================================================
use "$data/dataset.dta", clear

label variable pid         "Unique individual identifier"
label variable agefirst    "Age at first childbirth"
label variable age         "Individual's age at observation"
label variable reltime     "Time relative to first childbirth (age - (age at first birth))"
label variable labearnings "Annual labour earnings from employment"
label variable treated     "Treatment status: 1 if individual has experienced first childbirth, 0 otherwise"
label define treated_lbl 0 "Not yet treated" 1 "Treated"
label values treated treated_lbl


* Check the panel structure (there is one row per person and age)
isid pid age

* Check treatment timing 
assert agefirst <= age if treated == 1 & !missing(agefirst, age)
* Check treatment status
assert treated == (age >= agefirst) if !missing(agefirst)

*** 1a) Plotting number of individuals by the timing of treatment, so that to assess the distribution of the groups that recieved treatment at different points in time. (How many individuals are treated at each age?)
histogram agefirst, discrete freq title("Distribution of the groups that recieved treatment") color("eltblue") xtitle("Timing of treatment (age at first childbirth)") ytitle("Number of individuals")
graph export "1a.png", replace

**** Defining a helper program for plotting
* Plots one line per group (e.g. per age at first birth). Used for the descriptive graphs in 1b, 1c and 2b.
capture program drop plot_by_group
program define plot_by_group
    syntax varlist(min=2 max=2), group(varname) [title(string) xtitle(string) ///
        ytitle(string) xline(string) prefix(string) file(string) cols(integer 6)]
    gettoken y x : varlist
    levelsof `group', local(levels)
    local n : word count `levels'
    colorpalette hcl, n(`n') nograph
    local i 0
    local cmd ""
    local leg ""
    foreach g of local levels {
        local ++i
        local col `r(p`i')'
        local cmd `"`cmd' (connected `y' `x' if `group' == `g', lwidth(medium) msymbol(D) msize(vsmall) lcolor("`col'") mcolor("`col'"))"'
        local leg `"`leg' `i' "`prefix'`g'""'
    }
    if "`xline'" != "" local xl xline(`xline', lcolor(black) lpattern(dash))
    twoway `cmd', title("`title'") xtitle("`xtitle'") ytitle("`ytitle'") `xl' ///
        legend(order(`leg') size(vsmall) cols(`cols') pos(6))
    if "`file'" != "" graph export "$out/`file'", replace
end


*** 1b) Plot average outcome over time dimension by timing of treatment
preserve
collapse (mean) labearnings, by(age agefirst)
plot_by_group labearnings age, group(agefirst) prefix("Age ") ///
    title("Average annual earnings by age and timing of first birth") ///
    xtitle("Age") ytitle("Annual earnings") file("1b.png")
restore

*** 1c) Plot average outcome over relative event time by timing of treatment
* Lets me check whether groups follow similar earnings paths before birth and how earnings change once each group is treated.
preserve
collapse (mean) labearnings, by(reltime agefirst)
plot_by_group labearnings reltime, group(agefirst) prefix("Age ") xline(0) ///
    title("Average annual earnings by event time and timing of first birth") ///
    xtitle("Time relative to first childbirth") ytitle("Annual earnings") file("1c.png")
restore

*========================================================
*** 2. Conducting analysis of potential problem of forbidden comparison.
*========================================================
* With staggered timing, TWFE also compares newly treated units with units that were treated earlier. Their earnings are still changing because of their own treatment, so they are not a valid control (this is called "forbidden comparison").
* I illustrate this with two ages (28 and 29) and 4 cohorts:
* One group that recieved treatment in 1st period == 28; in 2nd period == 29; much earlier == 22; much later == 33
*   - treated at 29: switches into treatment in this window   -> the "treated" group
*   - treated at 33: not yet treated  -> correct, clean control
*   - treated at 22: treated long before -> already-treated control
*   - treated at 28: treated just before  -> already-treated control

*** 2a) Analyze the problem in static case (with just two periods). 
*** Restrict the sample to this window. (two consequitve time periods 28 and 29)
preserve
keep if age >= 28 & age <= 29
*** Restrict the sample to few groups in terms of treatment timing. 
keep if agefirst == 22 | agefirst == 28 | agefirst == 29 | agefirst == 33

*** 2b) Adjust outcomes of each group such that the trends from 28 to 29 are isolated for each group by subtracting the levels at 28, and plot these outcomes over 28 to 29 in one graph. 
collapse (mean) labearnings, by(age agefirst)
bysort agefirst (age): gen d_earn = labearnings - labearnings[1]   // change since age 28
plot_by_group d_earn age, group(agefirst) prefix("Age at birth ") cols(2) ///
    title("Earnings change since age 28 by timing of first birth") ///
    xtitle("Age") ytitle("Earnings change since age 28") file("2b.png")
restore

*** 2c) Make a table of treatment status by group and age
preserve
keep if age >= 28 & age <= 29
keep if agefirst == 22 | agefirst == 28 | agefirst == 29 | agefirst == 33
table agefirst age, statistic(mean treated) nformat(%9.2f)

*** 2d) Run the regression model (1) 𝑦_𝑖𝑎 = 𝛼 + 𝛽𝑡𝑟𝑒𝑎𝑡𝑒𝑑_𝑖𝑎 + 𝛾_𝑖 + 𝜃_𝑎 + 𝜖_𝑖𝑎 (1) on the restricted dataset with the groups chosen in 2a). (TWFE regression with individual and age FE)
xtset pid age
xtreg labearnings treated i.age, fe vce(cluster pid)
restore

*** 2e) Which comparisons influence the TWFE estimate?
* Run the regression (1) further restricting the sample to include:
* 1. Treated at p1 (29) vs later-treated (33)
preserve
keep if age >= 28 & age <= 29
keep if agefirst == 29 | agefirst == 33
xtset pid age
xtreg labearnings treated i.age, fe vce(cluster pid)
restore

* 2. Treated at p1 (29) vs all earlier-treated (22 and 28)
preserve
keep if age >= 28 & age <= 29
keep if agefirst == 29 | agefirst == 22 | agefirst == 28
xtset pid age
xtreg labearnings treated i.age, fe vce(cluster pid)
restore

* 3a. Treated at p1 (29) vs earlier-treated at 22
preserve
keep if age >= 28 & age <= 29
keep if agefirst == 29 | agefirst == 22
xtset pid age
xtreg labearnings treated i.age, fe vce(cluster pid)
restore

* 3b. Treated at p1 (29) vs earlier-treated at 28
preserve
keep if age >= 28 & age <= 29
keep if agefirst == 29 | agefirst == 28
xtset pid age
xtreg labearnings treated i.age, fe vce(cluster pid)
restore

*** 2f) Bacon decomposition on the 2a restricted sample
* Splits the TWFE estimate from 2d into its underlying 2x2 comparisons and
* shows the weight each one receives, including the forbidden ones.
preserve
keep if age >= 28 & age <= 29
keep if agefirst == 22 | agefirst == 28 | agefirst == 29 | agefirst == 33
xtset pid age
bacondecomp labearnings treated, ddetail
restore

*========================================================
*** 3. Constructing the stacked dataset
*========================================================
* For each cohort s (age at first birth), build a "sub-event" containing:
*   - treated: individuals with first birth at age s
*   - controls: the next G cohorts (first birth at s+1,...,s+G), kept only while they are not yet treated
* Each sub-event spans P years before and G years after the birth.
* Stacking the sub-events gives a dataset in which no already-treated unit is ever used as a control.

*** 3a) Specify the local parameters
* Define the number of next groups included in each control group (e.g G=5)
* Number of pre-periods for pre-trends test (e.g P=3)
use "$data/dataset.dta", clear
summarize age
local min_age = r(min)
local max_age = r(max)
summarize agefirst
local min_af = r(min)
local max_af = r(max)

local G = 5  // number of later cohorts used as controls (= years after birth)
local P = 3  // number of pre-periods
local s_min = `min_af' + `P'  // first cohort with P pre-periods in the data
local s_max = `max_af' - `G'  // last cohort with G later cohorts available as controls
local l_min = -`P' // event timeframe minimum
local l_max = `G' - 1 // event timeframe maximum

* Store as globals so Section 5 can use them
global P     = `P'
global s_min = `s_min'
global s_max = `s_max'
global l_min = `l_min'
global l_max = `l_max'


*** 3b) Create empty dataset to store stacked data
local stacked_data "stacked_data.dta"
clear
save `stacked_data', emptyok replace     // empty file; variables come from the appends

*** 3c) Create a loop over the specified values of treatment timing s
forvalues s = `s_min'/`s_max' {
    di _newline "age at first birth = `s'"
    use "$data/dataset.dta", clear
    local age_start = `s' - `P'
    local age_end = `s' + `G' - 1
    keep if (agefirst == `s' & age >= `age_start' & age <= `age_end') | ///
            (agefirst >= `s' + 1 & agefirst <= `s' + `G' & age >= `age_start' & age <= `age_end')
    drop if agefirst > `s' & age >= agefirst
    gen subevent = `s'
    gen treatedsubevent = (agefirst == `s')
    gen evtime = age - `s'
    tab evtime treatedsubevent
    * Ploting avg outcomes for T and C in this sub-event
    preserve
    collapse (mean) avg_labearnings = labearnings, by(evtime treatedsubevent)
    twoway (connected avg_labearnings evtime if treatedsubevent == 1, ///
                msymbol(O) lcolor(purple) mcolor(purple)) ///
           (connected avg_labearnings evtime if treatedsubevent == 0, ///
                msymbol(T) lcolor(eltblue) mcolor(eltblue)), ///
        xtitle("Event Time") ytitle("Average yearly earnings") ///
        title("Subevent `s'") ///
        legend(order(1 "Treatment group (agefirst=`s')" 2 "Control group(agefirst=`=`s'+1' to `=`s'+`G'')") ///
            size(small) pos(6)) ///
        xline(0, lcolor(black) lpattern(dash)) ///
        xlabel(`l_min'(1)`l_max') ///
        graphregion(color(white)) bgcolor(white)
    graph export "3d_sub_`s'.png", replace width(1000)
    restore
    append using `stacked_data'
    save `stacked_data', replace
}


*========================================================
*** 4. Comparison between conventional TWFE and stacked event study design
*========================================================
* 4a) Saturated TWFE regression fully interacted with sub-event identifier, to the stacked dataset
* One dummy per sub-event x event time for the treated (t = -1 is the reference).
* Fixed effects are interacted with the sub-event, so each sub-event is its DiD. SEs are clustered by person, because the same person can appear as a control in several sub-events.
use `stacked_data', clear
assert !missing(pid, age, subevent, treatedsubevent, evtime)
local n_sub = `s_max' - `s_min' + 1   // number of sub-events

forvalues s = `s_min'/`s_max' {
    forvalues l = 0/`=`l_max'-`l_min'' {
        local actual_l = `l' + `l_min'
        if `actual_l' != -1 {
            local lname = `actual_l' + `P'
            gen d_s`s'_l`lname' = (subevent == `s') & (treatedsubevent == 1) & (evtime == `actual_l')
        }
    }
}
reghdfe labearnings d_s*_l*, absorb(pid#subevent age#subevent) vce(cluster pid)

*** 4b) Weight each sub-event by its share of treated individuals, so the average reflects the effect for the typical treated mother.
keep if treatedsubevent == 1
bysort subevent pid: keep if _n == 1 // one row per treated person per sub-event
collapse (count) n_treated = pid, by(subevent)
egen total_treated = total(n_treated)
gen weight = n_treated / total_treated
matrix W = J(1, `n_sub', .)
forvalues s = `s_min'/`s_max' {
    local col = `s' - `s_min' + 1
    quietly su weight if subevent == `s'
    matrix W[1, `col'] = r(mean)
}

*** Storing the sub-event estimates and 95% CIs in a dataset for plotting
local n_obs = `n_sub' * (`l_max' - `l_min')
clear
set obs `n_obs'
gen subevent = .
gen evtime   = .
gen beta     = .
gen ci_lo    = .
gen ci_hi    = .
gen wavg     = .
gen wavg_lo  = .
gen wavg_hi  = .

local row = 0
forvalues s = `s_min'/`s_max' {
    forvalues l = 0/`=`l_max'-`l_min'' {
        local actual_l = `l' + `l_min'
        if `actual_l' != -1 {
            local row = `row' + 1
            local lname = `actual_l' + `P'
            quietly replace subevent = `s' in `row'
            quietly replace evtime  = `actual_l' in `row'
            quietly replace beta = _b[d_s`s'_l`lname'] in `row'
            quietly replace ci_lo = _b[d_s`s'_l`lname'] - 1.96*_se[d_s`s'_l`lname'] in `row'
            quietly replace ci_hi  = _b[d_s`s'_l`lname'] + 1.96*_se[d_s`s'_l`lname'] in `row'
        }
    }
}

* Add the reference period (t = -1, effect normalised to 0) for each sub-event
local n_obs2 = `row' + `n_sub'
set obs `n_obs2'
forvalues s = `s_min'/`s_max' {
    local row = `row' + 1
    quietly replace subevent = `s' in `row'
    quietly replace evtime = -1 in `row'
    quietly replace beta  = 0 in `row'
    quietly replace ci_lo = 0 in `row'
    quietly replace ci_hi = 0 in `row'
}

* Weighted average across sub-events for each event time; lincom gives the SE of the weighted sum, accounting for covariance between estimates.
forvalues l = 0/`=`l_max'-`l_min'' {
    local actual_l = `l' + `l_min'
    if `actual_l' != -1 {
        local lname = `actual_l' + `P'
        local lincom_cmd ""
        forvalues s = `s_min'/`s_max' {
            local col = `s' - `s_min' + 1
            local w = W[1, `col']
            if "`lincom_cmd'" == "" local lincom_cmd "`w'*d_s`s'_l`lname'"
            else local lincom_cmd "`lincom_cmd' + `w'*d_s`s'_l`lname'"
        }
        quietly lincom `lincom_cmd'
        quietly replace wavg  = r(estimate) if evtime == `actual_l' & subevent == `s_min'
        quietly replace wavg_lo = r(estimate) - 1.96*r(se) if evtime == `actual_l' & subevent == `s_min'
        quietly replace wavg_hi = r(estimate) + 1.96*r(se) if evtime == `actual_l' & subevent == `s_min'
    }
    else {
        quietly replace wavg = 0 if evtime == -1 & subevent == `s_min'
        quietly replace wavg_lo = 0 if evtime == -1 & subevent == `s_min'
        quietly replace wavg_hi = 0 if evtime == -1 & subevent == `s_min'
    }
}

sort subevent evtime

* Plot 1 Group-specific (estimates for each sub-event separately (are effects similar across cohorts?))
colorpalette hcl, n(`n_sub') nograph
local i = 0
local plotcmd ""
local legorder ""
forvalues s = `s_min'/`s_max' {
    local i = `i' + 1
    local col `r(p`i')'
    local plotcmd `"`plotcmd' (rcap ci_lo ci_hi evtime if subevent == `s', lcolor("`col'") lwidth(vvthin))"'
    local plotcmd `"`plotcmd' (connected beta evtime if subevent == `s', lwidth(thin) msymbol(O) msize(vsmall) lcolor("`col'") mcolor("`col'"))"'
    local legorder `"`legorder' `=(`i'-1)*2+1' "s=`s'""'
}
twoway `plotcmd', ///
    xtitle("Event time (wrt first birth)") ///
    ytitle("Effect on yearly earnings") ///
    title("Stacked DiD: Group-Specific Estimates") ///
    legend(order(`legorder') size(vsmall) cols(4) pos(6)) ///
    xline(0, lcolor(black) lpattern(dash)) ///
    yline(0, lcolor(gray) lpattern(dot)) ///
    xlabel(`l_min'(1)`l_max')
graph export "4b_group_specific.png", replace

* Plot 2 Weighted average (across sub-events (main stacked estimate))
twoway (rarea wavg_lo wavg_hi evtime if subevent == `s_min', color("255 182 193%50")) ///
       (connected wavg evtime if subevent == `s_min', lwidth(thick) msymbol(D) msize(medium) lcolor("219 112 147") mcolor("219 112 147")), ///
    xtitle("Event time (wrt first birth)") ///
    ytitle("Effect on yearly earnings") ///
    title("Stacked DiD: Weighted-Average Estimates") ///
    legend(order(2 "Weighted avg") size(small) pos(6)) ///
    xline(0, lcolor(gs8) lpattern(dash)) ///
    yline(0, lcolor(gs12) lpattern(dot)) ///
    xlabel(`l_min'(1)`l_max')
graph export "4b_weighted_avg.png", replace


* Save stacked wavg for comparison plot (in 4c)
preserve
keep if subevent == `s_min'
keep evtime wavg wavg_lo wavg_hi
save "stack_wavg.dta", replace
restore

*** 4c) Conventional TWFE event study with restrictions rfrom s_min to s_max
*** I am using the same cohorts and same event-time window as the stack, so that differences come only from the forbidden comparisons TWFE makes, not from the sample.
use "$data/dataset.dta", clear
keep if agefirst >= `s_min' & agefirst <= `s_max'
keep if inrange(reltime, `l_min', `l_max')   // only -1 omitted, same window as the stack

forvalues l = 0/`=`l_max'-`l_min'' {
    local actual_l = `l' + `l_min'
    if `actual_l' != -1 {
        local lname = `actual_l' + `P'
        gen d_l`lname' = (reltime == `actual_l')
    }
}

reghdfe labearnings d_l*, absorb(pid age) vce(cluster pid)

* Collect estimates and 95% CIs for plotting
local n_obs = `l_max' - `l_min'
clear
set obs `n_obs'
gen evtime = .
gen beta   = .
gen ci_lo  = .
gen ci_hi  = .

local row = 0
forvalues l = 0/`=`l_max'-`l_min'' {
    local actual_l = `l' + `l_min'
    if `actual_l' != -1 {
        local row = `row' + 1
        local lname = `actual_l' + `P'
        quietly replace evtime = `actual_l' in `row'
        quietly replace beta = _b[d_l`lname'] in `row'
        quietly replace ci_lo = _b[d_l`lname'] - 1.96*_se[d_l`lname'] in `row'
        quietly replace ci_hi = _b[d_l`lname'] + 1.96*_se[d_l`lname'] in `row'
    }
}


* Adding the reference period (t = -1)
local row = `row' + 1
set obs `row'
quietly replace evtime = -1 in `row'
quietly replace beta  = 0 in `row'
quietly replace ci_lo  = 0 in `row'
quietly replace ci_hi  = 0 in `row'

sort evtime

twoway (rarea ci_lo ci_hi evtime, color(gs12%40)) ///
       (connected beta evtime, lwidth(thick) msymbol(D) msize(medium) lcolor(cranberry) mcolor(cranberry)), ///
    xtitle("Event time (wrt first birth)") ///
    ytitle("Effect on yearly earnings") ///
    title("Conventional Event Study") ///
    legend(order(2 "TWFE estimate") size(small) pos(6)) ///
    xline(0, lcolor(black) lpattern(dash)) ///
    yline(0, lcolor(gray) lpattern(dot)) ///
    xlabel(`l_min'(1)`l_max')
graph export "4c.png", replace

* Comparison plot between stacked DiD vs TWFE. The difference shows the bias from forbidden comparisons in the TWFE estimator.
save "twfe_ests.dta", replace
merge 1:1 evtime using "stack_wavg.dta", nogen
rename (wavg wavg_lo wavg_hi) (beta_stack lo_stack hi_stack)
rename (beta ci_lo ci_hi) (beta_twfe lo_twfe hi_twfe)

twoway (rarea lo_stack hi_stack evtime, color("255 182 193%40")) ///
       (connected beta_stack evtime, lwidth(thick) msymbol(D) msize(medium) lcolor("219 112 147") mcolor("219 112 147")) ///
       (rarea lo_twfe hi_twfe evtime, color(purple%20)) ///
       (connected beta_twfe evtime, lwidth(thick) msymbol(O) msize(medium) lcolor(purple) mcolor(purple)), ///
    xtitle("Event time (wrt first birth)") ///
    ytitle("Effect on yearly earnings") ///
    title("Stacked DiD vs TWFE") ///
    legend(order(2 "Stacked DiD (weighted avg)" 4 "TWFE") size(small) pos(6)) ///
    xline(0, lcolor(black) lpattern(dash)) ///
    yline(0, lcolor(gray) lpattern(dot)) ///
    xlabel(`l_min'(1)`l_max')
graph export "4c_comparison.png", replace


*==============================================================================
*** 5. Child penalty in % of counterfactual earnings (bootstrap SEs)
*==============================================================================
* For comparability purpouses,  we express earnings loss as a % of what treated mothers would have earned without the birth.
* Counterfactual earnings of the treated = observed mean - estimated effect, so penalty(t) = sum_s w_s*b_st / sum_s w_s*(ybar_st - b_st), w_s = share of treated people.
* The ratio has no analytic SE -> cluster bootstrap over individuals.
* Index k = event time + 3 (k = 2, i.e. t = -1, is the reference period).

*** 5a) Stacked data with the same event-time dummies as in 4a
use "$out/stacked_data.dta", clear
forvalues s = $s_min/$s_max {
    foreach k in 0 1 3 4 5 6 7 {
        gen d_s`s'_l`k' = subevent == `s' & treatedsubevent == 1 & evtime == `k' - $P
    }
}
*** 5b) Define a program that runs the full estimation once and returns the results
* (since bootstrap needs the estimation wrapped in a program so it can rerun it)
capture program drop child_penalty
program define child_penalty, rclass
    args id                                   // pid, or the bootstrap id
	* Weights count each treated person once per sub-event
    quietly reghdfe labearnings d_s*_l*, absorb(`id'#subevent age#subevent)
    tempvar first                             // each treated person once per sub-event
    quietly bysort subevent `id': gen `first' = _n == 1 & treatedsubevent == 1
    quietly count if `first'
    local N = r(N)
    foreach k in 0 1 3 4 5 6 7 {
        local num = 0
        local den = 0
        forvalues s = $s_min/$s_max {
            quietly count if `first' & subevent == `s'
            local w = r(N) / `N'
            quietly summarize labearnings if treatedsubevent == 1 & subevent == `s' & evtime == `k' - $P, meanonly
            local num = `num' + `w' * _b[d_s`s'_l`k']
            local den = `den' + `w' * (r(mean) - _b[d_s`s'_l`k'])
        }
        local b`k'   = `num'
        local pct`k' = 100 * `num' / `den'  // effect in % of counterfactual earnings
    }
    foreach k in 0 1 3 4 5 6 7 {             
        return scalar b`k'   = `b`k''
        return scalar pct`k' = `pct`k''
    }
end

*** 5c) Cluster bootstrap over individuals
* I resample people instead of rows since a person appears in several sub-events, so all their rows must be drawn together. idcluster() gives each drawn copy a new id.
gen bs_pid = pid
bootstrap b0=r(b0) b1=r(b1) b3=r(b3) b4=r(b4) b5=r(b5) b6=r(b6) b7=r(b7)             ///
          pct0=r(pct0) pct1=r(pct1) pct3=r(pct3) pct4=r(pct4) pct5=r(pct5)            ///
          pct6=r(pct6) pct7=r(pct7), reps(499) seed(2026) cluster(pid) idcluster(bs_pid): ///
          child_penalty bs_pid
estimates store boot
local reps = e(N_reps)

*** 5d) Graph and table of the child penalty in %
local lbl pct0="-3" pct1="-2" pct3="0" pct4="1" pct5="2" pct6="3" pct7="4"
coefplot boot, keep(pct*) vertical coeflabels(`lbl') recast(connected) msymbol(D)    ///
    yline(0, lpattern(dot)) xline(2.5, lpattern(dash))                                ///
    xtitle("Event time (wrt first birth)") ytitle("% of counterfactual earnings")     ///
    title("Child penalty (stacked DiD)")                                              ///
    note("Reference period is t = -1. 95% CIs from bootstrap (`reps' repetitions)," ///
         "these appear narrower than the markers; hence SEs reported in 5_child_penalty.tex.")
graph export "$out/5_child_penalty_pct.png", replace

esttab boot using "$out/5_child_penalty.tex", replace keep(pct*)                   ///
    coeflabels(pct0 "t = -3" pct1 "t = -2" pct3 "t = 0" pct4 "t = 1"               ///
               pct5 "t = 2" pct6 "t = 3" pct7 "t = 4")                             ///
    b(%9.1f) se(%9.2f) nonumbers booktabs mtitles("% of counterfactual earnings")
	
log close
