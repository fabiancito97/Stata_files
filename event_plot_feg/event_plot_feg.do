
*** Main program

capture program drop event_plot_feg
program event_plot_feg, eclass
	version 17.0
	
	syntax [namelist(name=mat)] ,   ///
		pre_cof(string asis) /// prefix of leads
		post_cof(string asis) /// prefix of lags
		levels(numlist min=1 max=2) /// ci levels to graph
		[ /// 
		compar(string asis) /// if we have a comparison period, defaul compar == -1
		zero /// set coefficient = 0 in comparison period
		add_note(string) /// add note to graph
		did_imputation(integer 0) /// if use did_imputation results to grapg
		dropline /// vertical line in comparison period
		perturbline(string asis) /// move dropline 
		n_size(string asis) /// specify N 
		uci /// add uniform confidence intervals 
		ciplot(string asis) /// ciplot type
		onlypre /// graph only pre- coefficients
		note_stats /// add note of stats
		detrend /// perform detrend 
		ltrend /// show linear pre-trend
		*] /// other options are twoway custom graphs


*** Initializate
tempname estim betas results uniform_ci se results tabl table
tempname save_estimates

estimates store `save_estimates'

preserve
clear

quietly {

local comand1 = e(cmd)
local comand2 = e(cmd2)

local df = e(df_r)

if "`comand1'" == "reghdfe"  | "`comand1'" == "reghdfejl" | "`comand1'" == "lpdid" | "`comand2'" == "xtevent" {

*if "`comand1'" == "reghdfe"  | "`comand1'" == "reghdfejl" | "`comand2'" == "xtevent" {

* Output
local N = trim(string(e(N),"%-9.0fc")) // observations
local level  = r(level) // confidence level

* Confidence intervals 90% and 95%
matrix define `betas' = e(b)'
mata st_matrix("`se'",sqrt(diagonal(st_matrix("e(V)")))) // matrix of se

*}

*if "`comand1'" == "lpdid" {
*
*matrix define `table' = e(results) 
*matrix define `betas' = `table'[1,....]
*
*----
*
*mata st_matrix("`se'",sqrt(diagonal(st_matrix("e(V)")))) // matrix of se
*
*
*}

matrix `estim' = `betas'

local levels: list sort levels

foreach level of local levels{ 

local alpha = 1 - (`level'/100)
matrix `estim' = `estim', `betas' + invt(`df', `alpha'/2)*`se', `betas' + invt(`df',1-`alpha'/2)*`se'
}

matrix check1= `estim'
matrix check2= `betas'
matrix check3= `se'


* Multiple test of coefficients
local pre_test 
local post_test 
local index_pre = 0
local index_post = 0

mat `results' = J(1,`=colsof(`estim')'+1, .)



if "`comand2'" != "xtevent" local vars = e(indepvars)
if "`comand2'" == "xtevent" local vars = e(inexog)
if "`comand1'" == "lpdid"{

	local vars 
	local post_win = e(post_window)
	local pre_win = e(pre_window)
	forvalues h = 2(1)`pre_win'{
		local vars `vars' `pre_cof'`h'
	}
	
	forvalues h = 0(1)`post_win'{
		local vars `vars' `post_cof'`h'
	}
}

foreach var of local vars {
	if strpos("`var'", "`pre_cof'") {
		local pre_test `pre_test' (`var'=0)
		matrix `results' = `results' \ -real(substr("`var'", strlen("`pre_cof'")+1, strlen("`var'"))), `estim'[rownumb(`estim',"`var'"), ....] 
	}
	
	if strpos("`var'", "`post_cof'") {
		local post_test `post_test' (`var'=0)
		matrix `results' = `results' \ real(substr("`var'", strlen("`post_cof'")+1, strlen("`var'"))), `estim'[rownumb(`estim',"`var'"), ....]
	}
}

*test `pre_test'
*local p_pre = trim(string(r(p),"%-9.3fc")) 

*test `post_test'
*local p_pos = trim(string(r(p),"%-9.3fc"))

if "`uci'" != ""{

uci 

matrix `uniform_ci' = r(rtable2)
if "`comand2'" != "xtevent" matrix `uniform_ci' = `uniform_ci'[....,4..5]
if "`comand2'" == "xtevent" matrix `uniform_ci' = `uniform_ci'[2...,4..5]
matrix `uniform_ci' = J(1,2, .) \ `uniform_ci'

matrix `results' = `results', `uniform_ci'
}

}

if e(cmd) == "csdid2" | e(cmd) == "csdid"{
* Output
*local N = trim(string(`n_size',"%-9.0fc"))

*estat event
if e(cmd) == "csdid2"{
 matrix define `betas' = e(b)
 matrix `betas' = `betas'[1,1..7]
 matrix `betas' = `betas''
 
 }
 
if e(cmd) == "csdid" matrix define `betas' = e(table_90)[1..7,1]

matrix `estim' = `betas'
mat list `estim'

local levels: list sort levels
foreach level of local levels{ 

*estat event, level(`level')
*matrix tabl =  r(table)'
*matrix tabl = tabl[...., 5..6]

matrix `tabl' = e(table_`level')
 matrix `tabl' = `tabl'[...., 5..6]
 
mat list `estim' 
 
matrix `estim' = `estim', `tabl'

}

local vars: rownames `estim'
mat `results' = J(1,`=colsof(`estim')'+1, .)



foreach var of local vars {
	if strpos("`var'", "`pre_cof'") {
		local pre_test `pre_test' (`var'=0)
		matrix `results' = `results' \ -real(substr("`var'", strlen("`pre_cof'")+1, strlen("`var'"))), `estim'[rownumb(`estim',"`var'"), ....] 
	}
	
	if strpos("`var'", "`post_cof'") {
		local post_test `post_test' (`var'=0)
		matrix `results' = `results' \ real(substr("`var'", strlen("`post_cof'")+1, strlen("`var'"))), `estim'[rownumb(`estim',"`var'"), ....]
	}
}

if e(cmd) != "csdid"{
test `pre_test'
local p_pre = trim(string(r(p),"%-9.3fc")) 

test `post_test'
local p_pos = trim(string(r(p),"%-9.3fc"))

}

}


*** Prepare to make plot ----------------------------------------------------------------

mat list `results'

*matrix check =  `results'

svmat `results'
drop if `results'1 == .




*** Only pre graph

if "`onlypre'"!="" keep if `results'1 < -1

*** Zero, if require

if "`compar'" == "" local compar = -1

if "`zero'" != "" { 
	
	count 
    local obs = r(N) + 1
    set obs `obs'
    replace `results'1 = `compar' in `obs'
	replace `results'2 = 0   in `obs'
    replace `results'3 = 0   in `obs'
    replace `results'4 = 0   in `obs'
    
}

if "`comand2'" == "xtevent" {
	count 
    local obs = r(N) + 1
    set obs `obs'
    replace `results'1 = -2 in `obs'
	replace `results'2 = 0   in `obs'
    replace `results'3 = 0   in `obs'
    replace `results'4 = 0   in `obs'

}


*** Calculate pre-trend
reg `results'2 `results'1 if `results'1<=`compar'

tempvar pre_trend
predict `pre_trend', xb

*** Make adjustment
if "`detrend'"!=""{

replace `results'2=`results'2-`pre_trend'
*** first level
replace `results'3=`results'3-`pre_trend'
replace `results'4=`results'4-`pre_trend'
*** second level
replace `results'5=`results'5-`pre_trend'
replace `results'6=`results'6-`pre_trend'

sum `results'2 if `results'1==`compar'
local shift = r(min)
replace `results'2=`results'2-`shift'
*** first level
replace `results'3=`results'3-`shift'
replace `results'4=`results'4-`shift'
*** second level
replace `results'5=`results'5-`shift'
replace `results'6=`results'6-`shift'

}


*** Create labels and axis

sort `results'1

* x-axis
summarize `results'1 
local max_x = r(max)
local min_x = r(min)

* y-axis
*ds results1, not
*local variables = r(varlist)
*local variables : subinstr local variables " " ",", all
*summarize results3 
*local l_minimun = r(min)
*quietly: summarize results4 
*local u_maximun = r(max)

local vars=`=colsof(`results')'

local var_min

forvalues r=2(1)`vars'{
	local var_min `var_min' `results'`r'
}

tempvar min max

egen `min' = rowmin(`var_min')

egen `max' = rowmax(`var_min')

sum `min'
local l_minimun = r(min)
sum `max' 
local u_maximun = r(max)

limits, lims(`u_maximun' `l_minimun')
local max_y  = r(max_val)
local min_y  = r(min_val)
local delta  = r(delta)
local factor = r(factor)
local format = r(format)

tempvar range 


local labels
local val = `min_y' - `delta'
while `val' < `max_y'{
	local val = `val' + `delta'
*forvalues val = `min_y'(`delta')`max_y'{
	if round(`val', 1/`factor') != 0 {
		local newlab = trim(string(`val',"`format'"))
		local labels `"`labels' `val' "`newlab'""'
	}
	
}

local max_val = `val'

generate `range' = `max_val'

*display as text "`labels'"

*local labels `"`labels' 0 "0" "'
*local labels ylabel(`labels')

*local labels `"`labels' 0 "0" "'

*local labels ylabel(`labels')
local labels ylabel(`min_y'(`delta')`max_y')
*local labels ""

*** Indicate significant coeffficient
tempvar significant 
generate `significant' = (`results'3 > 0 & `results'4 > 0) | (`results'3 < 0 & `results'4 < 0) // 

*** Pre-trend
if "`detrend'"=="" & "`ltrend'" != "" local pre_trend_line (line `pre_trend' `results'1, lcolor(blue) lpattern("-"))

*** Line in comparison, if required
if "`dropline'" != "" {
	tempvar line1
	if "`perturbline'" == "" generate `line1' = -1  
	if "`perturbline'" != "" generate `line1' = -1 + `perturbline' 
	
	local vert_line (dropline `range' `line1', lcolor(black) lwidth(vthin) lpattern(dash) mcolor(none) base(`min_y'))

}

*** CI plot
if "`ciplot'" == "" local ciplot_cmd (rcap `results'3 `results'4 `results'1, lcolor(black)) 
if "`ciplot'" == "rcap" local ciplot_cmd (`ciplot' `results'3 `results'4 `results'1, lcolor(black))
if "`ciplot'" == "rarea" local ciplot_cmd (`ciplot' `results'3 `results'4 `results'1,  color(gs9))
if "`ciplot'" == "line" local ciplot_cmd (`ciplot' `results'3 `results'1,  lcolor(black) lpattern(dash)) (`ciplot' `results'4 `results'1,  lcolor(black) lpattern(dash))

*** CI plot, second level
if wordcount("`levels'") == 2{
	
	if "`ciplot'" == "" local ciplot2_cmd (rspike `results'5 `results'6 `results'1, lcolor(black) mlwidth(thin) msize(vsmall))
	if "`ciplot'" == "rcap" local ciplot2_cmd (rspike `results'5 `results'6 `results'1, lcolor(black) mlwidth(thin) msize(vsmall))
	if "`ciplot'" == "rarea" local ciplot2_cmd (`ciplot' `results'5 `results'6 `results'1,  color(gs11))
	if "`ciplot'" == "line" local ciplot2_cmd (`ciplot' `results'5 `results'6 `results'1,  lcolor(black) lpattern(dash)) (`ciplot' `results'5 `results'6 `results'1,  lcolor(black) lpattern(dash))
 
 }

*** CI plot, uniform CI
if "`uci'" != ""{
	if "`ciplot'" == "" local uci_graph (rspike `results'5 `results'6 `results'1, lcolor(black) mlwidth(thin) msize(vsmall))
	if "`ciplot'" == "rcap" local uci_graph (rspike `results'5 `results'6 `results'1, lcolor(black) mlwidth(thin) msize(vsmall))
	if "`ciplot'" == "rarea" local uci_graph (`ciplot' `results'5 `results'6 `results'1,  color(gs11))
	if "`ciplot'" == "line" local uci_graph (`ciplot' `results'5 `results'6 `results'1,  lcolor(black) lpattern(dash)) (`ciplot' `results'5 `results'6 `results'1,  lcolor(black) lpattern(dash))
 
}

*** Point of coefficient estimations
local point_estim (scatter `results'2 `results'1, mcolor(white) mlcolor(black) mlwidth(medthin))

*** Line indicator of 0
local yzero yline(0, lcolor(red))

*** Note N and p-values
if "`note_stats'" != "" local note_stats note("N = `N'" "p-value pre = `p_pre'" "p-value post = `p_pos'", size(medium))

*** Command 

*display as text "LABELS ----------------> `labels'"

local graph_run twoway `vert_line' `pre_trend_line' `ciplot2_cmd' `ciplot_cmd' `uci_graph' `point_estim', `labels' `yzero' `note_stats' `options' 


*** Run graph
`graph_run'

restore

estimates restore `save_estimates'

}
	
end

*#### Aux functions -----------------------------------------------------------------------------------------

*** To calculate the limits in axis of graph
capture program drop limits
program define limits, rclass
version 17.0
	syntax, lims(numlist min=2 max=2)
	
	*local lims `namelist'
	
	tokenize "`lims'"
	
	if !(`1'*(-1)>0 & `2'*(-1)>0) local lims: list sort lims
	
	tokenize "`lims'"
	
if `2' >= abs(`1') {

    local maximun = max(`2', abs(`1'))
 	
	up_lim, value(`maximun')
	local max_y = r(value)/r(factor)
	
	local format = "%9.2f"
	if r(factor) == 1000 local format = "%9.3f"
	local factor = r(factor)
	
	local delta = `max_y'/2
    
	while `max_y' < `maximun' {
		local max_y = `max_y' + `delta'
	}
	
    local min_y = `max_y' - `delta'
    while `min_y' > `1'{
    	local min_y = `min_y' - `delta'
    }
    
    if `min_y' >= 0  local min_y = - `delta'
	
	*** CORRECTION TO VERY LONG SCALE
	
	if `max_y' >= `maximun' + `delta' local max_y = `max_y' - `delta'
	

}

if `2' < abs(`1') {

    local minimun = `1' 

	local up_lim = -`minimun'
	
	up_lim, value(`up_lim')
	local min_y = -r(value)/r(factor)
	
	local factor = r(factor)
	local format = "%9.2f"
	if r(factor) == 1000 local format = "%9.3f"
	
	local delta = -`min_y'/2
	while `min_y' > `minimun' {
		local min_y = `min_y' - `delta'
	}
    
    local max_y = `min_y' + `delta'
    while `max_y' < `2'{
    	local max_y = `max_y' + `delta'
    }
    
    if `max_y' <= 0  local max_y = `delta'
	
	*** CORRECTION TO VERY LONG SCALE
	
	if `min_y' <= `minimun' - `delta' local min_y = `min_y' + `delta'
	*if `max_y' >= `maximun' + `delta' local max_y = `max_y' - `delta'
	

}
return scalar max_val = `max_y'
return scalar min_val = `min_y'
return scalar delta   = `delta'
return scalar factor  = `factor'
return local format  "`format'"

end

*** To calculate the min-max in axis of graph

capture program drop up_lim
program up_lim, rclass 
version 17.0
	syntax, value(numlist min=1 max=1)
	
    local factor = 1
    while `value' <= 10{
        local factor = `factor'*10
    	display as text "`factor'"
        local value = `value'*10 
    	display as text "`value'"
    }
    
    forvalues j = 0(20)80{
        
    if `value' >= `j' & `value' < `j' + 20 {
        local value = min(abs(`j'-`value') , abs(`j'+20-`value'))*(-(abs(`j'-`value') < abs(`j'+20-`value')) + (abs(`j'-`value') > abs(`j'+20-`value'))) + `value'
        break
        }
    }
	return scalar value = `value'
	return scalar factor = `factor'

end	


	

