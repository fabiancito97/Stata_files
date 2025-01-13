
global intent .19040466

capture program drop table_feg2

program table_feg2, eclass byable(recall)
    
	version 17.0
	
	syntax  namelist(min=1 max=1), ///
		    [ ///  /* optionals */
				pval /// present pval of coefficients
				pval_precoef /// for precoef, only present pval
				addstats(string asis) /// add stats to table
				addstats_labs(string asis) /// labels to added stats
				*]
			
			

qui{
			
tokenize `namelist'	

tempname starts starts_post starts_att starts_pre
tempname stats p betas se stable main 

tempname stable_post stable_att stable_pre

if "`pval'"!="" {
	matrix `starts_post'=J(1,3,0)
	matrix `starts_att'=J(1,1,0)
	
	
	if "`pval_precoef'"!="" matrix `starts_pre'=J(1,1,0)
	if "`pval_precoef'"=="" matrix `starts_pre'=J(1,3,0)
	
	local substat_sett substat(2) 
	local rtitle_cof_itt rtitle("\noalign{\medskip} ITT" \"" \ "") 
	if "`pval_precoef'"=="" local rtitle_cof_pre rtitle("\noalign{\medskip} Pre coef." \"" \"" )
	if "`pval_precoef'"!="" local rtitle_cof_pre rtitle("\noalign{\medskip} Pre p-val." )
	
}
if "`pval'"=="" {
	matrix `starts_post'=J(1,2,0)
	matrix `starts_att'=J(1,1,0)
	
	if "`pval_precoef'"!="" matrix `starts_pre'=J(1,1,0)
	if "`pval_precoef'"=="" matrix `starts_pre'=J(1,2,0)
	
	local substat_sett substat(1) 
	local cut_stable matrix `stable' = `stable'[....,1..2]
	local rtitle_cof_itt  rtitle("\noalign{\medskip} ITT" \"") 
	if "`pval_precoef'"=="" local rtitle_cof_pre rtitle("\noalign{\medskip} Pre coef." \"" )
	
}

*** Added stats

*matrix `stats' = e(mean) \ e(N) \ e(N_clust)


		

if "`addstats'" != "" {
	*local t = 0 // tot of variables and labels
	local l1 = 0 // number of stats
	local l2 = 0 // number of label


	foreach word of local addstats{
		 local l1 = `l1' + 1
		 if `l1' == 1 local added "`word'"
		 if `l1' != 1 local added `" `added' \ `word' "'
	}	
	
	foreach word of local addstats_labs{
		local l2 = `l2' + 1
		if `l2' == 1 local added_lab `""`word'""'
		if `l2' != 1 local added_lab `" `added_lab' \ "`word'" "'
	}
	
	
		*local tot_labels = `t'/2
		*
		*forvalues i = 1(1)`tot_labels'{
		*    local row`i' = `" "`lab_`i''" "'
		*}
	
	
	matrix `stats' = `added'
	frmttable, statmat(`stats') store(`stats') sdec(3 \ 0 \ 0) sfmt(fc) rtitle(`added_lab')
	
}



*matrix define `p' = e(table)' // p-values matrix
*matrix `p' =  `p'[...., colnumb(`p', "pvalue")] 
matrix define `betas' = e(b)' // betas matrix
mata st_matrix("`se'",sqrt(diagonal(st_matrix("e(V)"))))

*** p-value

test Dpre1
matrix define `p' = r(p) // betas matrix
test Dpost0
matrix define `p' = `p' \ r(p) // betas matrix

*** ATT
lincom _b[Dpost0]/${intent}
matrix `stable_att' = r(estimate)


matrix `stable_post' = (`betas'[2,1], `se'[2,1], `p'[2,1])  

if "`pval_precoef'"!="" matrix `stable_pre'  =`p'[1,1]
if "`pval_precoef'"=="" matrix `stable_pre'  =(`betas'[1,1], `se'[1,1], `p'[1,1])


*** Starts for table

matrix `starts_post'[1,1] = (`stable_post'[1,3]<0.10)+(`stable_post'[1,3]<0.05)+(`stable_post'[1,3]<0.01)
matrix `starts_att'[1,1] = (`stable_post'[1,3]<0.10)+(`stable_post'[1,3]<0.05)+(`stable_post'[1,3]<0.01)
if "`pval_precoef'"== "" matrix `starts_pre'[1,1] = (`stable_pre'[1,3]<0.10)+(`stable_pre'[1,3]<0.05)+(`stable_pre'[1,3]<0.01)
*if "`pval_precoef'"== ""

mat list `stable_post' 
mat list `stable_pre'


*matrix `stable' = `stable_post' \ `stable_pre'

*mat list `stable'


*matrix `starts' = `starts_post' \ `starts_pre'


*** Make table
`cut_stable' 

frmttable, statmat(`stable_att') store(`stable_att') sdec(3) annotate(`starts_att') asymbol(*,**,***) rtitle("\noalign{\medskip} ATT")
frmttable, statmat(`stable_post') `substat_sett' sdec(3) store(`stable_post') annotate(`starts_post') asymbol(*,**,***) `rtitle_cof_itt'
if "`pval_precoef'" == "" frmttable, statmat(`stable_pre') `substat_sett' sdec(3) store(`stable_pre') annotate(`starts_pre') asymbol(*,**,***) `rtitle_cof_pre'
if "`pval_precoef'" != "" frmttable, statmat(`stable_pre') sdec(3) store(`stable_pre') annotate(`starts_pre') asymbol(*,**,***) brackets([,]) `rtitle_cof_pre'



frmttable, replay(`stable_att') append(`stable_post') store(`main') 
frmttable, replay(`main') append(`stable_pre') store(`main') 


}

frmttable, replay(`main') append(`stats') store(`1') 




end

