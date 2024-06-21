
*** Main program

capture program drop event_plot_feg
program event_plot_feg, eclass
	version 17.0
	
	*local cmdline: copy local 0
	*gettoken mat 0 : 0, parse(",[") match(paren)
	
	*if "`mat'" == "," local mat = ""
	
	*display as text "`mat'"
	
	syntax [namelist(name=mat)] ,   ///
		pre_cof(string asis) ///
		post_cof(string asis) ///
		levels(numlist min=1 max=2) ///
		[ /// 
		compar(string asis) /// if we have a comparison period
		zero /// set coefficient = 0 in comparison period
		add_note(string) /// add note to graph
		did_imputation(integer 0) /// if use did_imputation results to grapg
		dropline /// vertical line
		perturbline(string asis) /// perturbation of line 
		n_size(string asis) /// specify N 
		uci /// add uniform confidence intervals 
		ciplot(string asis) /// ciplot type
		*] /// other options are twoway custom graphs



preserve
clear

quietly {

local comand1 = e(cmd)
local comand2 = e(cmd2)

if "`comand1'" == "reghdfe"  | "`comand1'" == "reghdfejl" | "`comand2'" == "xtevent" {

* Output
local N = trim(string(e(N),"%-9.0fc")) // observations
local level  = r(level) // confidence level

* Confidence intervals 90% and 95%
matrix define betas = e(b)'
mata st_matrix("se",sqrt(diagonal(st_matrix("e(V)")))) // matrix of se

matrix estim = betas

local levels: list sort levels

foreach level of local levels{ 

local alpha = 1 - (`level'/100)
matrix estim = estim, betas + invt(e(df_r), `alpha'/2)*se, betas + invt(e(df_r),1-`alpha'/2)*se
}

* Multiple test of coefficients
local pre_test 
local post_test 
local index_pre = 0
local index_post = 0

mat results = J(1,`=colsof(estim)'+1, .)

if "`comand2'" != "xtevent" local vars = e(indepvars)
if "`comand2'" == "xtevent" local vars = e(inexog)
display as text "`vars'"


foreach var of local vars {
	if strpos("`var'", "`pre_cof'") {
		local pre_test `pre_test' (`var'=0)
		matrix results = results \ -real(substr("`var'", strlen("`pre_cof'")+1, strlen("`var'"))), estim[rownumb(estim,"`var'"), ....] 
	}
	
	if strpos("`var'", "`post_cof'") {
		local post_test `post_test' (`var'=0)
		matrix results = results \ real(substr("`var'", strlen("`post_cof'")+1, strlen("`var'"))), estim[rownumb(estim,"`var'"), ....]
	}
}

test `pre_test'
local p_pre = trim(string(r(p),"%-9.3fc")) 

test `post_test'
local p_pos = trim(string(r(p),"%-9.3fc"))

if "`uci'" != ""{

uci 

matrix uniform_ci = r(rtable2)
if "`comand2'" != "xtevent" matrix uniform_ci = uniform_ci[....,4..5]
if "`comand2'" == "xtevent" matrix uniform_ci = uniform_ci[2...,4..5]
matrix uniform_ci = J(1,2, .) \ uniform_ci

matrix results = results, uniform_ci
}

}

if e(cmd) == "csdid2" {
* Output
local N = trim(string(`n_size',"%-9.0fc"))

*estat event
matrix define betas = e(b)'

matrix estim = betas

local levels: list sort levels
foreach level of local levels{ 

*estat event, level(`level')
*matrix tabl =  r(table)'
*matrix tabl = tabl[...., 5..6]

matrix tabl = e(table_`level')
 matrix tabl = tabl[...., 5..6]
 
matrix estim = estim, tabl

}

local vars: rownames estim
mat results = J(1,`=colsof(estim)'+1, .)



foreach var of local vars {
	if strpos("`var'", "`pre_cof'") {
		local pre_test `pre_test' (`var'=0)
		matrix results = results \ -real(substr("`var'", strlen("`pre_cof'")+1, strlen("`var'"))), estim[rownumb(estim,"`var'"), ....] 
	}
	
	if strpos("`var'", "`post_cof'") {
		local post_test `post_test' (`var'=0)
		matrix results = results \ real(substr("`var'", strlen("`post_cof'")+1, strlen("`var'"))), estim[rownumb(estim,"`var'"), ....]
	}
}

test `pre_test'
local p_pre = trim(string(r(p),"%-9.3fc")) 

test `post_test'
local p_pos = trim(string(r(p),"%-9.3fc"))



}


clear
svmat results

*drop if results1 == .

if "`zero'" != "" & "`compar'" =="" { 
	count 
    local obs = r(N) + 1
    set obs `obs'
    replace results1 = -1 in `obs'
	replace results2 = 0   in `obs'
    replace results3 = 0   in `obs'
    replace results4 = 0   in `obs'
    
}

if "`zero'" != "" & "`compar'" !="" { 
	count 
    local obs = r(N) + 1
    set obs `obs'
    replace results1 = -`compar' in `obs'
	replace results2 = 0   in `obs'
    replace results3 = 0   in `obs'
    replace results4 = 0   in `obs'
    
}

if "`comand2'" == "xtevent" {
	count 
    local obs = r(N) + 1
    set obs `obs'
    replace results1 = -2 in `obs'
	replace results2 = 0   in `obs'
    replace results3 = 0   in `obs'
    replace results4 = 0   in `obs'

}

sort results1

tempvar range line1 significant // for vertical lines

* x-axis
summarize results1 
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

local vars=`=colsof(results)'

local var_min

forvalues r=2(1)`vars'{
	local var_min `var_min' results`r'
}

*ds results1, not
*local var_min = r(varlist)

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

generate `range' = `max_y'

* labels
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
local labels `"`labels' 0 "0" "'

local labels ylabel(`labels')
local labels ylabel(`min_y'(`delta')`max_y')


*if "`compar'" != "" generate `line1' = `compar' + `perturbline' //+ 0.5

generate `significant' = (results3 > 0 & results4 > 0) | (results3 < 0 & results4 < 0) // significant coef

*format `format'  results2 results3 results4 `range' 

if "`dropline'" != "" {
	if "`perturbline'" == "" generate `line1' = -1  
	if "`perturbline'" != "" generate `line1' = -1 + `perturbline' 
	
	local vert_line (dropline `range' `line1', lcolor(black) lwidth(vthin) lpattern(dash) mcolor(none) base(`min_y'))

}


if wordcount("`levels'") == 2 local sencond_level (rspike results5 results6 results1, lcolor(black) mlwidth(thin) msize(vsmall))

if "`uci'" != "" local uci_graph (rspike results5 results6 results1, lcolor(black) mlwidth(thin) msize(vsmall))

if "`ciplot'" == "" local ciplot rcap
if "`ciplot'" == "rcap" local ciplot_cmd (`ciplot' results3 results4 results1, lcolor(black))
if "`ciplot'" == "rarea" local ciplot_cmd (`ciplot' results3 results4 results1,  color(gs9))
if "`ciplot'" == "line" local ciplot_cmd (`ciplot' results3 results1,  lcolor(black) lpattern(dash)) (`ciplot' results4 results1,  lcolor(black) lpattern(dash))


twoway  `vert_line' /// vertical line 
	   `ciplot_cmd' ///
	   `sencond_level' ///
	   `uci_graph' ///
       /// (line results2 results1, lcolor(black))  ///
	   /// (scatter results2 results1 if `significant' == 0, mcolor(white) mlcolor(black) mlwidth(medthin)) ///
	   /// (scatter results2 results1 if `significant' == 1, mcolor(red) mlcolor(black) mlwidth(medthin)) ///
	   (scatter results2 results1, mcolor(white) mlcolor(black) mlwidth(medthin)) ///
	   , ///
	   /// ylabel(`min_y'(`delta')`max_y') ///
	   `labels' ///
	   yline(0, lcolor(red)) ///
	   note("N = `N'" "p-value pre = `p_pre'" "p-value post = `p_pos'", size(medium)) ///
	   `options'
	
restore

*drop results*

}
	
end

*#### -----------------------------------------------------------------------------------------

capture program drop heter_plot_feg
program heter_plot_feg, eclass
	version 17.0
	
	*local cmdline: copy local 0
	*gettoken mat 0 : 0, parse(",[") match(paren)
	
	*if "`mat'" == "," local mat = ""
	
	*display as text "`mat'"
	
	syntax namelist(name=mat), [  ///
		single ///
		*] /// other options are twoway custom graphs
	

quietly {

preserve
clear
matrix define results = `mat'

svmat results
drop if results1 == .

* y-axis
ds results1, not
local var_min = r(varlist)

tempvar min max

egen `min' = rowmin(`var_min')

egen `max' = rowmax(`var_min')

sum `min'
local l_minimun = r(min)
sum `max' 
local u_maximun = r(max)

limits, lims( `l_minimun' `u_maximun')
local max_y  = r(max_val)
local min_y  = r(min_val)
local delta  = r(delta)
local factor = r(factor)
local format = r(format)

* labels
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
local labels `"`labels' 0 "0" "'

if "`single'" == "" {

twoway (rcap    results3 results4 results1 if results1 == 2, lcolor(blue)) ///
       (rspike  results5 results6 results1 if results1 == 2, lcolor(blue)) /// 
	   (scatter results2 results1 if results1 == 2, mcolor(blue) mlcolor(blue) mlwidth(medthin) msymbol(T)) ///
	   (rcap   results3 results4 results1 if results1 == 1, lcolor(green)) ///
	   (rspike results5 results6 results1 if results1 == 1, lcolor(green)) ///
	   (scatter results2 results1 if results1 == 1, mcolor(green) mlcolor(green) mlwidth(medthin) msymbol(S)) ///
	   (rcap   results3 results4 results1 if results1 == 0, lcolor(orange)) ///
	   (rspike results5 results6 results1 if results1 == 0, lcolor(orange)) ///
	   (scatter results2 results1 if results1 == 0, mcolor(orange) mlcolor(orange) mlwidth(medthin)) ///
	   , ///
	   ylabel(`labels') ///
	   yline(0, lcolor(red)) ///
	   `options'
}

if "`single'" != "" {

twoway (rcap results3 results4 results1, lcolor(green)) ///
	  (rspike  results5 results6 results1, lcolor(green)) /// 
	   (scatter results2 results1, mcolor(green) mlcolor(green) mlwidth(medthin) msymbol(S)) ///
	   , ///
	   ylabel(`labels') ///
	   yline(0, lcolor(red)) ///
	   `options'
}	

restore
}
	
end

*#### -----------------------------------------------------------------------------------------

capture program drop heter_plot_mult
program heter_plot_mult, eclass
	version 17.0
	
	*local cmdline: copy local 0
	*gettoken mat 0 : 0, parse(",[") match(paren)
	
	*if "`mat'" == "," local mat = ""
	
	*display as text "`mat'"
	
	syntax namelist(name=mat), [  ///
		single ///
		*] /// other options are twoway custom graphs
	

quietly {

preserve
clear
matrix define results = `mat'

svmat results
drop if results1 == .

* y-axis
ds results1 results2, not
local var_min = r(varlist)


tempvar min max
egen `min' = rowmin(`var_min')

egen `max' = rowmax(`var_min')

sum `min'
local l_minimun = r(min)
sum `max' 
local u_maximun = r(max)

limits, lims( `l_minimun' `u_maximun')
local max_y  = r(max_val)
local min_y  = r(min_val)
local delta  = r(delta)
local factor = r(factor)
local format = r(format)

* labels
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
local labels `"`labels' 0 "0" "'
local labels ylabel(`labels')

local text = `max_y'
 
twoway (rcap    results4 results5 results1 if results1 == 2 & results2 == 1, lcolor(blue)) ///
       (rspike  results6 results7 results1 if results1 == 2 & results2 == 1, lcolor(blue)) /// 
	   (scatter results3 results1 if results1 == 2 & results2 == 1, mcolor(blue) mlcolor(blue) mlwidth(medthin) msymbol(T)) ///
	   (rcap    results4 results5 results1 if results1 == 1 & results2 == 1, lcolor(green)) ///
	   (rspike  results6 results7 results1 if results1 == 1 & results2 == 1, lcolor(green)) ///
	   (scatter results3 results1 if results1 == 1 & results2 == 1, mcolor(green) mlcolor(green) mlwidth(medthin) msymbol(S)) ///
	   (rcap    results4 results5 results1 if results1 == 0 & results2 == 1, lcolor(orange)) ///
	   (rspike  results6 results7 results1 if results1 == 0 & results2 == 1, lcolor(orange)) ///
	   (scatter results3 results1 if results1 == 0 & results2 == 1, mcolor(orange) mlcolor(orange) mlwidth(medthin)) ///
	   ///
	   (rcap    results4 results5 results1 if results1 == 6 & results2 == 2, lcolor(blue)) ///
       (rspike  results6 results7 results1 if results1 == 6 & results2 == 2, lcolor(blue)) /// 
	   (scatter results3 results1 if results1 == 6 & results2 == 2, mcolor(blue) mlcolor(blue) mlwidth(medthin) msymbol(T)) ///
	   (rcap    results4 results5 results1 if results1 == 5 & results2 == 2, lcolor(green)) ///
	   (rspike  results6 results7 results1 if results1 == 5 & results2 == 2, lcolor(green)) ///
	   (scatter results3 results1 if results1 == 5 & results2 == 2, mcolor(green) mlcolor(green) mlwidth(medthin) msymbol(S)) ///
	   (rcap    results4 results5 results1 if results1 == 4 & results2 == 2, lcolor(orange)) ///
	   (rspike  results6 results7 results1 if results1 == 4 & results2 == 2, lcolor(orange)) ///
	   (scatter results3 results1 if results1 == 4 & results2 == 2, mcolor(orange) mlcolor(orange) mlwidth(medthin)) ///
	   , `options' ///
	   `labels' ///
	   yline(0, lcolor(red)) ///
	   text(`text' 1 "OLS", box bcolor(white%100) fcolor(white%100) ) ///
	   text(`text' 5 "CSDID", box bcolor(white%100) fcolor(white%100))

	
restore
}
	
end

*#### -----------------------------------------------------------------------------------------

capture program drop event_make_matrix
program define event_make_matrix
version 16.0
/*
 From a regression (with reghdfe), make the matrix to plot an event study
*/
	syntax namelist(name=vars)
	
	tokenize "`vars'"
	
	matrix estim = `1'
	*matrix estim = r(table)
	
	matrix results = .,.,.,.
	
foreach var of local vars  {
	matrix results =  results \ ///
	estim[rownumb(estim,"b"),  colnumb(estim, "`var'")], /// coefficient
    estim[rownumb(estim,"ll"), colnumb(estim, "`var'")], /// low limit ci
	estim[rownumb(estim,"ul"), colnumb(estim, "`var'")], /// upper limit ci
	real(substr("`var'", strlen("`var'") - 3, strlen("`var'")))
}
	
end


********** To calculate the limits in axis of graph

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

********** To calculate the min-max in axis of graph

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


	

