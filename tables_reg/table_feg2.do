
capture program drop table_feg2
program table_feg2, eclass byable(recall)
{
    version 17.0
	
	local cmdline : copy local 0
	// local 0 save the cmdline
		local stop 0
		local tot_betas = 0
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
                local tot_betas = `tot_betas' + 1 // sum equations
				if "`paren'" != "" {		
					local alleqn "`alleqn' (`eqn')"
				}
				else    local alleqn `alleqn' `eqn' 
			}
		}
		local 0 `"`eqn' `0'"'

local options `0'

syntax , ///
		    [ ///  /* optionals */
				pval /// present pval of coefficients
				pval_precoef /// for precoef, only present pval
				addstats(string asis) /// add stats to table
				store(string) /// name of store
				beta_matrix(string) /// matrix with betas, e(b) is default
				se_matrix(string) /// matrix with V, e(V) is default
				*]
			
tempname starts starts_post starts_att starts_pre
tempname stats p betas se stable main 

tempname stable_post stable_att stable_pre

capture frmttable, clear(`store')

qui{

	forvalues k = 1(1)`tot_betas'{
		    gettoken eqn alleqn : alleqn, parse(" ,[") match(paren)
			gettoken beta other : eqn, parse(" ,[") match(paren)

		tempname `beta'_tab
		*maketable `eqn' coeff_matrix(`betas') se_matrix(`se') store(``beta'_tab')
		maketable `eqn' store(``beta'_tab')
		
		if `k' == 1 frmttable, replay(``beta'_tab') store(`store')
		if `k' != 1 frmttable, replay(`store') append(``beta'_tab') store(`store')
	}	
		
*** Added stats

if `"`addstats'"' != ""{
	addstatscmd `addstats' nametab(`stats')
	frmttable, replay(`store') append(`stats') store(`store')
}	

}

frmttable, replay(`store')

}
end

*#### Auxiliary programs ------------------------------------------------

capture program drop addstatscmd
program addstatscmd, eclass byable(recall)
{
	version 17.0
	
	syntax  anything, [ ///
			nametab(string) /// name of the store tab
			*] //
				
	local l1 = 0 
	foreach word of local anything{
		 local l1 = `l1' + 1
		 if `l1' == 1 local added "`word'"
		 if `l1' != 1 local added `" `added' \ `word' "'
	}	
	tempname stats_mat
	matrix `stats_mat' = `added'	
	frmttable, statmat(`stats_mat') store(`nametab') `options'
}
end


capture program drop maketable
program maketable, eclass byable(recall)
{
    version 17.0
	
	syntax namelist(min=1 max=1), ///
		    [ ///  /* optionals */
				substats(string asis) /// sub stats of coefficient
				coeffdec(string) /// decimals of coefficient
				starts(string) /// position for stars, = 0 if not whant stars. default == 1
				starts_symbol(string asis) /// symbols for starts, default = *,**,***
				coefflabel(string asis) /// label of coefficient, default = cmd line
				brackets(string asis) /// brackest to present coeff + substats, default ="","" \ (,) \ [,]
				coeff_matrix(string asis) /// matrix that store the coefficients, default previous regression
				se_matrix(string asis) /// matrix that store the se
				v_matrix(string) /// covariance matrix, default previous regression (doesnt work if se_matrix is declared)
				store(string) /// table where store
				*]
	
	*** default values
	
	if "`starts_symbol'" == "" local starts_symbol `"*,**,***"'
	if "`brackets'" == "" local brackets `""",""\(,)\[,]"'
	if `""`coefflabel'""' == "" local coefflabel `namelist'
	if "`starts'" == "" local starts=1
	
	if "`coeff_matrix'"=="" {
		tempname coeff_matrix
		matrix `coeff_matrix'= e(b)' 
	}		
	
	if "`v_matrix'"=="" {
		tempname v_matrix
		matrix `v_matrix' = e(V)
	}
	
	if "`se_matrix'"=="" {
		tempname se_matrix
		mata st_matrix("`se_matrix'",sqrt(diagonal(st_matrix("`v_matrix'"))))
		local rownames : rownames `coeff_matrix'
		matrix rownames `se_matrix' = `rownames'
	}

	
	local coeff_name `namelist'
	
	local coeff = `coeff_matrix'[rownumb(`coeff_matrix',"`coeff_name'"), 1]
	local se = `se_matrix'[rownumb(`se_matrix',"`coeff_name'"), 1]
	local p = 2*(normal(-(abs(`coeff')/`se')))
		
	local matrix_cmd `coeff'
	local tot_stats = 1
	local rtitle "`coefflabel'"
	
	foreach stat of local substats{
		
		local ++tot_stats
		
		if "`stat'"=="se" {
			local matrix_cmd `matrix_cmd', `se'
			local rtitle `"`rtitle' \ "" "'
		}
		if "`stat'"=="p" {
			local matrix_cmd `matrix_cmd', `p'
			local rtitle `"`rtitle' \ "" "'
		}
		
		
			
			
	}
	tempname matrix_store starts_matrix
	matrix `matrix_store' = `matrix_cmd'  
			
	matrix `starts_matrix'=J(1,`tot_stats',0)
			
	
			
	if `starts' != 0 matrix `starts_matrix'[1,`starts'] = (`p'<0.10)+(`p'<0.05)+(`p'<0.01)
	
	local tot_substats = `tot_stats'-1
	
	mat list `starts_matrix'	
	
	*display `starts_symbol'
	
	frmttable, statmat(`matrix_store') substat(`tot_substats') sdec(`coeffdec') store(`store') annotate(`starts_matrix') asymbol(`starts_symbol')	brackets(`brackets') rtitle(`rtitle') 


}
end	


capture program drop IsStop
program define IsStop // this functiond detect a character and stop the loop
{
	args stop colon token
	if 	     `"`token'"' == "[" /*
		*/ | `"`token'"' == "," /*
		*/ | `"`token'"' == "if" /*
		*/ | `"`token'"' == "in" /*
		*/ | `"`token'"' == "" {
		c_local `stop' 1
	}
	else	c_local `stop' 0
}	
end


