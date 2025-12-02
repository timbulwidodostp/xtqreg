*! v.1.0.0. N.Orsini, M Bottai 31mar2017
* First created on 22oct2015
* Updated by Giola Santoni on 04apr2016
* Updated by Matteo Bottai on 31mar2017

capture program drop xtqreg
program xtqreg, eclass properties(mi) byable(onecall)
version 13

if _by() {
		local BY `"by `_byvars'`_byrc0':"'
}
	
`version' `BY' _vce_parserun xtqreg, mark(OFFset CLuster) : `0'

if "`s(exit)'" != "" {
		version 10: ereturn local cmdline `"xtqreg `0'"'
		exit
	}

if replay() {
		if ("`e(cmd)'"!="xtqreg")  error 301  
		Replay `0'
	}
else `version' `BY' Estimate `0'
end

capture program drop Estimate 
version 12, missing

program Estimate, eclass byable(recall)
	syntax varlist  [if] [in]  ///
	    [,  Level(integer $S_level) ///
	   Quantiles(numlist) ///
	   Random(string)  ///
	   Reps(string) ///
	   SEED(string) ///
	   Method(string) ///
	   COVariance(string) ///
	   NOConstant /// 
	   pathr(string) ///
	   Version(string) ///
	   DIRectory(string) ///
	 * ]
	   
    local cmdline : copy local 0
	
	preserve 
		
	marksample touse 
	
*************Start Added********************************************************
	if "`version'"=="" local version "3.1.3"
	if "`pathr'"=="" & "$Rterm_path"==""{
		local pathr "C:\Program Files\R\R-`version'\bin\x64\R.exe"
		noi di as text "Path for R.exe: " as res "`pathr'"
		}
	else if "`pathr'"=="" local pathr "$Rterm_path"

	// Check validity of directory for saving and reading data
	local cwd `c(pwd)'
	if substr("`c(pwd)'",1,2)=="\\"{
		if "`directory'"==""{
			di  as err "Directory on the server." _n "Please specify a directory on the computer"
			exit
			}
		else local cwd "`directory'"
		}
  	if "`c(os)'" == "Windows" local cwd : subinstr local cwd "\" "/" , all  
*************End Added**********************************************************

	// Get the panel variable from xtset or tsset
	capture quietly xtset 
	if _rc != 0 {
		di as err _n "panel variable not set; xtset or tsset the data"
		exit 198
	}
	local group `r(panelvar)'
	tempvar uniquev
	qui by `touse' `group', sort: gen `uniquev' = _n==1
	qui replace `uniquev' = sum(`uniquev') if `touse'
	qui sum `uniquev' if `touse', meanonly
	local ngroups = r(max)
	
	// Check specified quantiles 

		tempname quantlist
		SetQ `quantiles'
		local listquant "`r(quants)'"
		tokenize "`r(quants)'"
		
		local nq 1
		while "``nq''" != "" {
			if (`nq' == 1) {
			                mat `quantlist' = ``nq''
							local listqr "c(``nq''"
			}						
			else {
					mat `quantlist' = (`quantlist' , ``nq'')
					local listqr "`listqr', ``nq''"
			}
		
			local q`nq' ``nq''
			local nq = `nq' + 1
		}
		local nq = `nq' - 1
		local listqr "`listqr')"
 
// set the seed
			if "`seed'" == "" {
				local tmp = round(runiform()*100000000)
				local seed `tmp'
			}

// Check number of bootstrap replications
// Set asymptotic variance estimator as default
 
 		if "`reps'" == "" {
							local nreps = 20
		}
		else {
				if `reps'<2  {
						di  as err "specify more than 2 replications"
						exit  
				}
				local nreps = `reps'
		}
		
// Check number of observations 

		quietly count if `touse'
		if r(N)<4 { 
			di in red "insufficient observations"
			exit 2001
		}

// Remove collinear variables

		gettoken depv indepv : varlist		
		_rmcoll `indepv' [`weight'`exp'] 
		local indepv `r(varlist)'  
	    qui regress `depv' `indepv' if `touse'
 		local nobs `e(N)'
		
		tokenize "`indepv'"
		local i 1
		while "``i''" != "" {
			if _se[``i''] == 0 {
				di in blu /*
			*/ "note: ``i'' dropped because of collinearity"
				local `i' " "
			}
			local i = `i' + 1
		}
		
// Check random effects

		_rmcoll `random' [`weight'`exp'] 
		local randomeff `r(varlist)'  

	foreach v of local randomeff {
		if regexm("`indepv'", "`v'") != 1 {
			di as err "specify a random effect that is included in the fixed effect"
			exit 198
		}
	}	
// Check quantiles and get a list 

		local pow = 10
		local done = 0
		while !`done' {
			local pow = `pow'*10
			local done = 1
			forvalues k=1/`nq' {
				local q = (`q`k'')*`pow'
				local done = `done' &  (reldif(`q',trunc(`q'))<c(epsfloat))			
			}
		}
	    
		local pow = 10
		local pow = `pow'*10
		forvalues k=1/`nq' {
			local q = round(`q`k''*`pow')
			if length(string(`q')) == 1 local myeq "q0`q'"
			else local myeq "q`q'"
			
			local eqnames `eqnames' `myeq'
		}
 
		local eqnames : list uniq eqnames
		local k : word count `eqnames'
	
		if `k' < `nq' {
			di as err "only `k' out of `nq' quantiles are " /*
			 */ "unique within a relative precision of c(epsfloat)"
			exit 498
		}
				
 	    if "`noconstant'" == "" local conams "`indepv' _cons"
 	    else local conams "`indepv'"

		tokenize "`eqnames'"
		forv i = 1/`nq' {
			foreach v of local conams {
		    	local eqnams "`eqnams' ``i''"
			}
 		}
        
        local getit "`conams'"
 		forv i = 2/`nq' {
 			local conams "`conams' `getit'"
  		}
 
//  Dimension of the beta vector
		
	local indepv "`indepv'"
	local ncols  : word count `eqnams'
	if "`noconstant'" != ""  local ncov = `: word count `indepv''  
	else local ncov = `: word count `indepv'' + 1
	tempname coefs

// Check the optimization method 

	if "`method'" == "" local mopt "gs"
	else local mopt "`method'"
	
	if inlist("`mopt'", "gs", "df") != 1 {
			di as err `"`method'" not allowed"' 	
			exit 198
			}
			
// Create the R formula for fixed and random effects

	local xfixed "1"
    if "`noconstant'" != "" local xfixed "0"
	tokenize "`indepv'"
		local i 1
		while "``i''" != "" {
			local xfixed "`xfixed' + ``i''"
			local i = `i' + 1
		}
	
	local xrandom "1"
    if "`noconstant'" != "" local xrandom "0"

		tokenize "`randomeff'"
		local i 1
		while "``i''" != "" {
			local xrandom "`xrandom' + ``i''"
			local i = `i' + 1
		}

* Eventually erase created datasets*********************************************
quietly {
		capture rm `"`cwd'/d_xtqreg.dta"'
		capture rm `"`cwd'/bootf_xtqreg.dta"'
		capture rm `"`cwd'/b_xtqreg.dta"'
		capture rm `"`cwd'/bootr_xtqreg.dta"'
		capture rm `"`cwd'/vr_xtqreg.dta"'
} // end quietly****************************************************************
		
	
// save a dataset in txt format
 
quietly keep  if `touse' == 1
quietly keep `depv' `indepv' `group' `touse'

cap quietly saveold `"`cwd'/d_xtqreg"' , version(11)  replace nolabel
if _rc!=0 quietly saveold `"`cwd'/d_xtqreg"' , replace nolabel

// make sure -rsource- is installed 

quietly capture which rsource
if _rc != 0 qui ssc install rsource

// Covariance of the random effect 

if "`covariance'"  == "" {
	 local covr = "pdDiag"
	 local covartype "Independent"
}
if substr("`covariance'",1,2) == "in" {
     local covr = "pdDiag"
	 local covartype "Independent"
}
if substr("`covariance'",1,2) == "un" {
	local covr = "pdSymm"
	 local covartype "Unstructured"
}
if substr("`covariance'",1,2) == "ex" {
	 local covr = "pdCompSymm"
	 local covartype "Exchangeable"
}
if substr("`covariance'",1,2) == "id" {
	local covr = "pdIdent"
	 local covartype "Identity"
}

if "`covr'" == "" {
		di as err `""`covariance'" covariance not allowed"'
		exit 198
}
		
// Define labels for the random effects

local varre "_cons"
if "`noconstant'" != "" local varre ""

foreach v of local randomeff {
	local varre "`varre' `v'"
}

local nre : word count `varre'

// If unstructured, exchangeable (add the covariances)

local conamsre "`varre'"
local covtype "long"

if inlist("`covr'", "pdSymm", "pdCompSymm")==1 {
	local covre ""
	tokenize "`varre'"
	forv a = 1/`=`nre'-1' {
   		forv b = `=`a'+1'/`nre' {    
        	local covre "`covre' ``a''_``b''"
			local i = `i'+ 1
		}
		local conamsre "`varre' `covre'"
		local covtype "wide"
	}
}

local eqnamsre ""
local nvcre : word count `conamsre'

forv k = 1/`nq' {
	local myeq : word `k' of `eqnames'
	foreach v of local conamsre {
		local eqnamsre "`eqnamsre' `myeq'"
	}
}

local conamsrefull ""
forv k = 1/`nq' {
	local conamsrefull "`conamsrefull' `conamsre'"
}
 
// create a file with the R syntax and get results
 
  	//local cwd `c(pwd)'
  	//if "`c(os)'" == "Windows" local cwd : subinstr local cwd "\" "/" , all  
 quietly {
    tempname cmdr
    file open 	`cmdr'  using "`cwd'/xtqreg_to_r.R",  write text replace all   
	file write  `cmdr'  `"setwd("`cwd'")"' _n
	file write  `cmdr'  "library(lqmm)" _n	
	file write  `cmdr'  "library(foreign)"  _n
	file write  `cmdr'  `"mydata <- read.dta("d_xtqreg.dta")"' _n
    file write  `cmdr'   "datalqmm <- data.frame(mydata)" _n
	file write  `cmdr'  `"fit.lqmm <- lqmm(`depv' ~ `xfixed', "' _n
	file write  `cmdr'  `"random = ~ `xrandom', group = `group',"' _n
	file write  `cmdr'  `"control = list(method = "`mopt'"),"' _n
	file write  `cmdr'  `"tau = `listqr', covariance = "`covr'", data = datalqmm)"' _n
	file write  `cmdr'  `"b = data.frame(coef(fit.lqmm))"' _n
	file write  `cmdr'  `"write.dta(b, "b_xtqreg.dta")"' _n
	file write  `cmdr'  `"vr = data.frame(VarCorr(fit.lqmm))"' _n
	file write  `cmdr'  `"write.dta(vr, "vr_xtqreg.dta")"' _n
	file write  `cmdr'  `"bootstrap <- boot.lqmm(fit.lqmm, R=`nreps', seed = `seed')"' _n
	file write  `cmdr'  `"bootf = data.frame(lqmm:::extractBoot(bootstrap, which = "fixed"))"' _n
	file write  `cmdr'  `"write.dta(bootf, "bootf_xtqreg.dta")"' _n
	file write  `cmdr'  `"bootr = data.frame(lqmm:::extractBoot(bootstrap, which = "random"))"' _n
	file write  `cmdr'  `"write.dta(bootr, "bootr_xtqreg.dta")"' _n    
	file close  `cmdr' 
   rsource using "`cwd'/xtqreg_to_r.R" ,  lsource roptions(--slave) rpath("`pathr'") //noloutput 
  }


* viewsource xtqreg_to_r.R
  
//  Prepare b of the fixed effects
 
capture confirm file "`cwd'/b_xtqreg.dta"
if _rc != 0 {
    	      //di as err _n "the R script did not run, check the R path"
			  di as err _n "Please check if the R.exe directory is `pathr' and R version is `version'" _n ///
			   "If the directory is not correct use option pathr." _n ///
			   "If only the version is not correct use option version." _n ///
			   "If both the directory and the version are correct," _n ///
			   "something went wrong during the run of the R script."
			  exit 198   
   		 }

quietly  {
	use `"`cwd'/b_xtqreg"', clear
    
  	if "`noconstant'" == "" { 
   	 	set obs `=`c(N)'+1'
		foreach v of varlist *  {
    	            replace `v' = `v'[1] in l
   		 }
		drop in 1
  	}
	 
	tempname bq
	local i = 1 
	foreach v of varlist *  {
                mkmat `v' , matrix(`bq')
                mat list `bq'
                if `i' == 1 mat `coefs' = `bq'' 
                else mat `coefs' = `coefs' , `bq''
                local i = `i' + 1
     }
     mat colnames `coefs' = `conams'
	 mat coleq `coefs' = `eqnams'
} // end quietly
	
//  Prepare V  of the fixed effects
		
quietly  {
	use `"`cwd'/bootf_xtqreg"', clear
	drop scale*
	ds *
	list
	local listR  "`r(varlist)'"
    local a = 1 
    local b = `ncov'
    forv i = 1/`nq' {
    	local firstname : word `a' of `listR'
    	local lastname : word `b' of `listR'
    	order `firstname' , after(`lastname')
    	local a = `a' + `ncov'
    	local b = `b' + `ncov'
    }
    
    tempname VCE
    mat accum `VCE' = *, dev nocons
	mat rownames `VCE' = `conams'
	mat roweq `VCE' = `eqnams'
	mat colnames `VCE' = `conams'
    mat coleq `VCE' = `eqnams'
		
	mat `VCE'=`VCE'*(1/(`nreps'-1))
    
} // end quietly

*   Prepare b of the random effects
 
quietly  {
	use `"`cwd'/vr_xtqreg"', clear
	clist
	tempname coefsr bqr
 	if ("`covtype'" == "long") {
		local i = 1 
	    foreach v of varlist *  {
                mkmat `v' , matrix(`bqr')
                 if `i' == 1 mat `coefsr' = `bqr'' 
                else mat `coefsr' = `coefsr' , `bqr''
                local i = `i' + 1
       }
	}

 	if ("`covtype'" == "wide") {
 		tempname getvc getall one vd 
 		mkmat _all, matrix(`getall')
 		local c = 1
		forv i = 1/`nq' {
				 tempname getvc`i'
 	 			 mat `getvc`i'' = `getall'[1..`nre', `c'..`=`c'+`nre'-1']
	 			 mat `vd' = vecdiag(`getvc`i'')

	 		     if `c' == 1 mat `coefsr' = `vd' 
                 else mat `coefsr' = `coefsr' , `vd'
                 
	 			 forv a = 1/`=`nre'-1' {
   						forv b = `=`a'+1'/`nre' {    
        					mat `one' = `getvc`i''[`a', `b']
        					mat `coefsr' = `coefsr' , `one'
						}
				}
	 			 
	 			 local c = `c' + `nre'			
			}	
		}

     mat colnames `coefsr' = `conamsrefull'
	 mat coleq `coefsr' = `eqnamsre'

} // end quietly
   	
	//  Prepare V  of the random effects

quietly  {
	
	use `"`cwd'/bootr_xtqreg"', clear
    tempname VCER
	
	if ("`covtype'" == "long") {
   	 mat accum `VCER' = *, dev nocons
   	 mat list `VCER'
	 mat rownames `VCER' = `conamsre'
	 mat roweq `VCER' = `eqnamsre'
	 mat colnames `VCER' = `conamsre'
     mat coleq `VCER' = `eqnamsre'
	 mat `VCER'=`VCER'*(1/(`nreps'-1))
	}
 	
 	 if ("`covtype'" == "wide") {
   	 mat accum `VCER' = *, dev nocons
   	 mat list `VCER'
	 mat rownames `VCER' = `conamsre'
	 mat roweq `VCER' = `eqnamsre'
	 mat colnames `VCER' = `conamsre'
     mat coleq `VCER' = `eqnamsre'
	 mat `VCER'=`VCER'*(1/(`nreps'-1))
	}
 } // end quietly

// Display results 
 	
 		qui use `"`cwd'/d_xtqreg"' , clear
  		ereturn post `coefs' `VCE', obs(`nobs')  depn(`depv') 

		ereturn local qlist "`listquant'"

		ereturn repost, esample(`touse')
		
		ereturn local depvar "`depv'"

		ereturn scalar reps = `nreps'
		ereturn scalar seed = `seed'
		ereturn scalar N = `nobs'
		ereturn scalar groups = `ngroups'
		ereturn scalar n_q = `nq'
		ereturn local eqnames "`eqnames'"
		ereturn local method "`mopt'"
		ereturn local vcetype "Bootstrap"
		ereturn local predict "sqreg_p"
		ereturn local cmdline `"xtqreg `cmdline'"'
		ereturn local cmd "xtqreg"
		ereturn local ivars "`group'"
		ereturn local redim "`nre'"
		ereturn local vartypes "`covartype'"
		ereturn local revars "`randomeff'"
	    ereturn local depvar "`depv'"

		di _n in gr "Mixed-effects Quantile regression"       _c
		di in gr _col(53) "No. of obs    = " in ye %10.0g e(N)
		if "`e(method)'" == "df" di in gr "Optimization: " in y "Nelder-Mead" in gr _col(53) "No. of groups = " in y %10.0f  e(groups) 
		else di in gr "Optimization: " in y "Gradient-Search" in gr _col(53) "No. of groups = " in y %10.0f  e(groups)
	    di  in gr "Group variable: " as res e(ivars)   in gr _col(50) "No. of bootstrap = " in y %10.0f e(reps)   
		di in gr "Covariance random-effects: " in y "`e(vartypes)'" 
		ereturn display, level(`level')  

		di as txt  "Random-effects Parameters" 
		 _coef_table, neq(`e(n_q)') level(`level')  bmatrix(`coefsr') vmatrix(`VCER')  coeftitle("Estimate") 
		ereturn matrix br = `coefsr'
		ereturn matrix Vr = `VCER'
		/* Eventually erase created datasets
		quietly {
		rm `"d_xtqreg.dta"'
		rm `"boot_xtqreg.dta"'
		rm `"b_xtqreg.dta"'
		rm `"vr_xtqreg.dta"'
		} // end quietly
		*/
end


capture program drop Replay
program Replay
	syntax [, Level(cilevel) ]
	ereturn display, level(`level')  
end

// Sub-programs similar to sqreg

capture program drop SetQ
program define SetQ /* <nothing> | # [,] # ... */ , rclass
	version 6.0
		if "`*'"=="" {
		ret local quants ".5"
		exit
	}
	local orig "`*'"
	tokenize "`*'", parse(" ,")

	while "`1'" != "" {
		FixNumb "`orig'" `1'
		ret local quants "`return(quants)' `r(q)'"
		mac shift 
		if "`1'"=="," {
			mac shift
		}
	}
end

capture program drop FixNumb
program define FixNumb /* # */ , rclass
	version 6.0
	local orig "`1'"
	mac shift
	capture confirm number `1'
	if _rc {
		Invalid "`orig'" "`1' not a number"
	}
	if `1' >= 1 {
		ret local q = `1'/100
	}
	else 	ret local q `1'
	if `return(q)'<=0 | `return(q)'>=1 {
		Invalid "`orig'" "`return(q)' out of range"
	}
end
		

capture program drop Invalid 
program define Invalid /* "<orig>" "<extra>" */
	version 6.0
	di in red "quantiles(`1') invalid"
	if "`2'" != "" {
		di in red "`2'"
	}
	exit 198
end

capture program drop PrForm
program define PrForm /* # */ , rclass
	local aa : di %8.2f `1'
	ret local pr `aa'
	if substr("`return(pr)'",1,1)=="0" {
		ret local pr = substr("`return(pr)'",2,.)
	}
end

exit
 /*
use http://www.imm.ki.se/biostatistics/data/wtloss, clear
reshape long y, i(id) j(month)
gen inter = month*prog
xtset id month
* 1. Set the R path once for all
//global Rterm_path "C:\Program Files\R\R-4.5.0\bin\R.exe"
xtqreg y  prog month 
exit
set trace off
xtqreg y  prog month
exit
xtqreg y prog month inter ,  q(25 75)  pathr("/usr/bin/r")
 * 2. Set the R path as option
xtqreg y  prog month    , pathr("/usr/bin/r") //pathr("C:\Program Files\R\R-4.5.0\bin\R.exe")
  */
