	
*#### Main program -----------------------------------------------------------------------

capture program drop table_regs_feg

program table_regs_feg, eclass byable(recall)
    
	version 16.0

	local cmdline : copy local 0
	// local 0 save the cmdline
		local stop 0
		local tot_eq = 0
		while !`stop' { 
			gettoken eqn 0 : 0, parse(" ,[") match(paren)
			/*
			gettoken first 0 : 0
			  which obtains the first token from `0' and saves the rest back in `0'
			   parse(" ,[") parsing characters
			   match(paren) specifies that parentheses be matched in determining the token
			*/
			
			
			
			IsStop stop : `"`eqn'"'
			if !`stop' {
                local tot_eq = `tot_eq' + 1 // sum equations
				if "`paren'" != "" {		
					local alleqn "`alleqn' (`eqn')"
				}
				else    local alleqn `alleqn' `eqn' 
			}
		}
		local 0 `"`eqn' `0'"'
		
		
		*#### Syntax of command
		
		syntax, ///
		    contfe(string asis)    /// contros-FE and labels
			treatment(string asis) /// treatment variable/s
		    [ ///  /* optionals */
				fetitle(string)   /// title for FE section
			    fortable(string asis)  /// format table 
			    foroutput(string asis) /// format table final output
				varlabels /// use var treatment labels
				ctitles(string asis)   /// titles of table
			    save_file(string asis) /// path save table
				verbose 			   /// display intermediate tables
				by(varlist)            /// varlist of categorical variables to define groups
			* ] // other options for outreg	
		
		
		
		*** Treatment vars and labels
		if "`varlabels'"!=""{
		local t_treat=0 //tot treatments
		foreach word of local treatment{
		    local ++t_treat
			local label: variable label `word'
			if "`label'" != "" local lab_treat`t_treat' `label'
			if "`label'" == "" local lab_treat`t_treat' `word'
			 
			
        }
		}
		else{
		local t_treat=0
		foreach word of local treatment{
		    local ++t_treat
			local lab_treat`t_treat' `word'
			
        }
		}
		*/
		
		if "`verbose'"  == "" set output error
		{
		
		
		
		if "`by'" != ""{ // create groups in order to make different estimations
		    
			
			local j = 0
			local vars_group
			foreach var of local by{  
			    local j = `j' + 1
				
			    tempvar g_`j'
				clonevar `g_`j'' = `var'
				
				quietly summarize `var'
				local max = r(max)
				local min = r(min)
				forvalues v = `min'(1)`max'{
				    local labels `labels' `v' "`var' == `v'"
				}
				
				capture label drop g_`j'
				label define g_`j' `labels'
				label value `g_`j'' g_`j'
				
				local vars_group `vars_group' `g_`j''
				
			}
			
			tempvar groups
		    egen `groups' = group(`vars_group'), label
            quietly summarize `groups'
			local min_groups = r(min)
			local max_groups = r(max)
			
		}
		
		if "`fetitle'" == "" {
			local add_rows = `"" \noalign{\medskip} \textit{Controles y efectos fijos}" "'
		}
		if "`fetitle'" != "" {
			local add_rows = `" "`fetitle'" "'
		}
		
		local t = 0 // tot of variables and labels
		local l1 = 0 // number of label
		local l2 = 0 // number of variable
		
		foreach word of local contfe{
		    local t = `t' + 1
			if  mod(`t',2) == 0{
			    local l1 = `l1' + 1
				local lab_`l1' `word'
			}
			else {
			    local l2 = `l2' + 1
				local var_`l2' `word'
			}
        }		
		local tot_labels = `t'/2
		
		forvalues i = 1(1)`tot_labels'{
		    local row`i' = `" "`lab_`i''" "'
		}
		
		forvalues k = 1(1)`tot_eq'{
		    gettoken eqn alleqn : alleqn, parse(" ,[") match(paren)
			
			`eqn'
			outreg, keep(`treatment') store(`k') `options'
			
			if "`by'" != ""{
			forvalues g = `min_groups'(1)`max_groups'{
			preserve
				keep if `groups' == `g'
			    `eqn'
				outreg, keep(`treatment') store(`k'_`g') `options'
			restore
			}
			
			}
			
			*** Return all controls variables and fixed effects
			if e(cmd) == "reghdfe"{
			    local all_vars  `e(indepvars)' `e(absvars)' 
			}
			
			if e(cmd) == "regress"{
			    local all_vars `e(cmdline)'
			}
			
			local stop 0
			local fe_vars 
			while `stop' != 1{
			    gettoken var all_vars : all_vars, parse(" ,[")
			    IsStop stop : `"`var'"'
				if `stop' != 1 & ///
				   !inlist("`var'", "`e(cmd)'", "`e(depvar)'", "`treatment'", "_cons") & ///
				   strpos("`treatment'", "`var'") == 0{
				    local fe_vars `fe_vars' `var' 
				}
			}
			
			*** Detect and complete rows
			forvalues f = 1(1)`tot_labels'{
			    local yes_no = "No"
				local rep = 1
				while "`yes_no'" == "No" & `rep' <= wordcount("`fe_vars'"){
					local fe_var : word `rep' of `fe_vars'
					
	                if strmatch("`var_`f''", "`fe_var'") > 0 local yes_no = "Si"
					
					local rep = `rep' + 1
				}
				local row`f' = `" `row`f'' , "`yes_no'" "'
            }
			
			if `k' != 1{
			    outreg, replay(1) merge(`k') store(1)
				
				if "`by'" != ""{
				forvalues g = `min_groups'(1)`max_groups'{
				    outreg, replay(1_`g') merge(`k'_`g') store(1_`g')
				}
				}
			}
			
		}
		
		*** Final table
		if "`by'" != ""{
			forvalues g = `min_groups'(1)`max_groups'{
				outreg, replay(1) append(1_`g') store(1)
		    }
		}
		
		forvalues f = 1(1)`tot_labels'{
	        local add_rows = `" `add_rows' \ `row`f'' "'
        }
		
		*if "`ctitles'" != "" {
		*    local ctitles = "~" // next create titles of columns
        *    forvalues c = 1(1)`tot_eq'{
		*    	local ctitles = `" `ctitles', (`c')"'
		*    }	
		*}
		
		local ctitles `ctitles'
		
		
		if `"`ctitles'"' == "" {
		    local ctitles = "~" // next create titles of columns
            forvalues c = 1(1)`tot_eq'{
		    	local ctitles = `" `ctitles', (`c')"'
		    }	
		}
		
	
		

		if "`by'" == ""{
		local labtreat `""`lab_treat1' ""' 
		if `t_treat' > 1{
			forvalue k = 2(1)`t_treat'{
				local labtreat `"`labtreat' \ "" \ "`lab_treat`k''""'
			
			}
		}
			
		    local row_titles = `" `labtreat'   \ "" \ "R2"\"N" "'
		    *local row_titles = `"\noalign{\medskip} `treatment'"   \ "" \ "\noalign{\medskip}R2"\"N"\ "'
		}
		
		if "`by'" != ""{
		
		if `t_treat' > 1{
		local labtreat `" `lab_treat1'  "' 
			forvalue k = 2(1)`t_treat'{
				local labtreat `"`labtreat' " \ "" \ "`lab_treat`k''""'
			
			}
		}
		if `t_treat' == 1 local labtreat `"  `lab_treat1'  " "' 
		
		    local row_titles = `""\noalign{\medskip} \textbf{Panel: Todos} \\ `labtreat' \ "" \ "R2"\"N""'
			
			forvalues g = `min_groups'(1)`max_groups'{
			    local lab_value = "`:label (`groups') `g''"
				local row_titles = `"`row_titles' \ "\noalign{\medskip} \textbf{Panel: `lab_value'} \\ `labtreat' \ "" \ "R2"\"N""'
			}
		
		}
		
		if "`save_file'" != "" local save_file using `save_file'
		
		
		}
		
		if "`verbose'"  == "" set output proc
		
		outreg `save_file', replay(1) statfont(fs12) nolegend  /// 
		ctitle(`ctitles') ///
		rtitle(`row_titles') ///
		addrows(`add_rows') ///
		`options' 		 
end

*#### Auxiliary functions ---------------------------------------------------------

capture program drop IsStop
program define IsStop // this functiond detect a character and stop the loop
	args stop colon token
	if 	     `"`token'"' == "[" /*
		*/ | `"`token'"' == "," /*
		*/ | `"`token'"' == "if" /*
		*/ | `"`token'"' == "in" /*
		*/ | `"`token'"' == "" {
		c_local `stop' 1
	}
	else	c_local `stop' 0
end

