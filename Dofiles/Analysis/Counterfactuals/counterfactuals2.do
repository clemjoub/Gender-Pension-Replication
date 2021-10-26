	**PRELIMINARIES******************************************************************************

	clear
	clear matrix
	clear mata
	set mem 30000000
	set seed 1
	

	* path
	 if c(username)=="WB452275" {
		global path "C:\Users\WB452275\Documents\GitHub\gender-pension\"
    global datapath "C:\Users\wb452275\Dropbox\DATA\"
	 }
	 else if c(username)=="clement" {
		global path "C:\Users\clement\Dropbox\Projects\Genderequity"
	 }
		else{
		global path "C:\Genderequity"
	 }
	 
	global dopath "$path\Dofiles\Analysis\Counterfactuals"
	global resultspath = "$path\Counterfactuals"
	*global graphpath = "C:\Users\clement\Dropbox\Projects\Genderequity\draft\draft032016\Graphs"
	*global tablepath = "C:\Users\clement\Dropbox\Projects\Genderequity\draft\draft032016\Tables"
	global graphpath= "$path\Figures\"
	global tablepath = "$path\Tables\"

	* variables and parameters
	global varlist1 "id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption wdeath hdeath halive walive"

	*global year = "year>2008 & year<2020"
	global Tub = 75
	global Tlb = 36
	global graphformat = "png"
	global exportoptions = "width(800) height(600)"
	global counterfactuals = "0 1 3 4 5 6 7 8 10"

	global load = 0
	
	* define samples
*	global samplenewretirees="hhagegroup5==65 & year>2009 & year<2029" //impacts on variables at retirement age
*	global sampleatret="hhage==64 & year>2009 & year<2029" //impacts on variables at retirement age
*	global samplenearretirement="hhagegroup5==60 & year>2009 & year<2029" //impacts on behavior near retirement
*	global samplecosts="hhage>64 & year>2009 & year<2029" 		 //impacts on all retirees


	
	global samplegendergap="hhagegroup5==65 & year>2009 & year<2029" 				
	global samplestocks="$samplegendergap" 									
	global samplelaborsupply="hhage<65 & hhage>55 & year>2009 & year<2029" 	
	global samplecosts="hhage>=60 & year>2009 & year<2029" 		 		
	global samplesupplement="year>2007 & year<2029" 		 		
	

	global outputfile = "\Tables_Draft"
	cd "$path"

	**DATA******************************************************************************

	
	
	***load data

	if $load ==1 {
		do "$dopath\LoadCounterfactuals.do"
	}

	***create variables 
	use "$resultspath\AppendedCounterfactuals", clear
		rename id folio
		capture drop _merge

		if c(username)=="WB452275" {
    merge m:1 folio using "${datapath}\Chile\EPSsecondary\EPSsamples\EPSweights.dta"
    }
    else
    {
    		merge m:1 folio using "C:\Genderequity\DATA\RawData\EPSweights.dta"
     }

		*ta _merge
		keep if inlist(_merge, 2,3)
		*codebook weights04
		drop if counterfactual==.
		rename folio id
		
		drop if inlist(counterfactual,2,9)

		do $path\dofiles\variabledefs.do

if c(username)=="WB452275" {
    do $path\dofiles\counterfactuals\variabledefs_counterfactuals.do
} 
else
{
		do $path\dofiles\Analysis\counterfactuals\variabledefs_counterfactuals.do
}

	***label counterfactuals
	label define counterfactual 0 "Old System" 1 "Baseline reform" 2 "14% contributions" 3 "Equal retirement age" 4 "Contribution split with non-working wife" 5 "Extra contribution for non-working wife" 6 "Gender neutral life tables" 7 "Reform without child bonus" 8 "Reform without divorce rule" 9 "Reform without APS" 10 "All features" 11 "Equal earning opportunities"
	label values counterfactual counterfactual
		

	***show that demographics of projections make sense
	*drop if $year==0

	*ta year hhagegroup5  	 [aw=weights06] if clone==1 & counterfactual==0 & inlist(year,2004,2014,2024) & hhage>64, row nofreq

	*table year agein2004 [aw=weights06] if clone==1 & counterfactual==0, c(n married mean hhage)
	*table year agein2004 [aw=weights06] if clone==1 & counterfactual==0 & female==1 & married==0, c(n female  mean wdeath  mean walive)
	*table year married HH_edugrp [aw=weights06] if counterfactual==0 & hhagegroup5==60 & year>2008 & year<2029 & female==1, c(mean LFP)

	***look at individual histories
	*sort  folio clone year counterfactual
	*br counterfactual folio clone year type married female wage hage shw sww sformh sformw pensionw pensionh hfexper wfexper earnh earnw sdiv wdeath hdeath supplement if inlist(counterfactual,1,2)

	***compare two counterfactuals
	/*local c1=1
	local c2=2

	foreach var of varlist pension penwithdraw hhfexper wealth {
		capture drop `var'1 
		capture drop `var'2
		capture drop `var'diff 
		capture drop `var'diffpct 
		capture drop temp1
		capture drop temp2
		g temp1=`var' if counterfactual==`c1'
		g temp2=`var' if counterfactual==`c2'
		bysort folio clone year : egen `var'1=max(temp1)
		bysort folio clone year : egen `var'2=max(temp2)

		g `var'diff=`var'2-`var'1
		g `var'diffpct=`var'diff/`var'1
		

	}
	su pensiondiff if year==2028 & hhagegroup5==65 & counterfactual==1, d
	su pensiondiffpct if year==2028 & hhagegroup5==65 & counterfactual==1, d

	su penwithdrawdiff if year==2028 & hhagegroup5==65 & counterfactual==1, d
	su penwithdrawdiffpct if year==2028 & hhagegroup5==65 & counterfactual==1, d

	su hhfexperdiff if year==2028 & hhagegroup5==65 & counterfactual==1, d
	su hhfexperdiffpct if year==2028 & hhagegroup5==65 & counterfactual==1, d

	su wealthdiff if year==2028 & hhagegroup5==65 & counterfactual==1, d
	su wealthdiffpct if year==2028 & hhagegroup5==65 & counterfactual==1, d	
	*/
	
	* create gender specific pension quartiles
	save "$resultspath\CounterfactualsFigures.dta", replace
	
**OUTPUTS******************************************************************************


	use "$resultspath\CounterfactualsFigures.dta", clear	
	
**Tabulate 99******************************************************************************	
use "$resultspath\99.dta", clear

tabulate W_edugrp, generate(W_edugrp)
tabulate H_edugrp, generate(H_edugrp)
tabstat married numk W_edugrp1 W_edugrp2 W_edugrp3 W_edugrp4 H_edugrp1 H_edugrp2 H_edugrp3 H_edugrp4 sww shw, by(type) s(mean) save
mat T1 = r(Stat1)
mat T2 = r(Stat2)
mat T3 = r(Stat3)
mat T4 = r(Stat4)
putexcel set "Summary.xlsx"
putexcel A1 = matrix(T1), names
putexcel A3 = matrix(T2), names
putexcel A5 = matrix(T3), names
putexcel A7 = matrix(T4), names

tab2xl W_edugrp  type using $tablepath$outputfile, col(1) row(1)
tab2xl H_edugrp  type using $tablepath$outputfile, col(1) row(8)
tab2xl HH_edugrp  type using $tablepath$outputfile, col(1) row(15)
tab2xl married  type using $tablepath$outputfile, col(1) row(22)
tab2xl numk  type using $tablepath$outputfile, col(1) row(27)
tab2xl LFP  type using $tablepath$outputfile, col(1) row(36)

	
	replace numk=1 if numk==2
	replace numk=3 if numk>3 & numk!=.
	replace hhfexpergrp=20 if hhfexpergrp>20 & hhfexpergrp!=.
	
	keep if inrange(hhfexpergrp,0,40) 
	keep if inrange(year,2008,2029)
	g ones=1
	
	
	
	* Impacts on the Gender Gap
	
	g pensionquartile = .
 
 foreach count in $counterfactuals { 
	*creates quartiles that are specific to each gender in each counterfactual

		foreach gender in 0 1 {
			di "`count' `gender'"
			capture drop temp
			xtile temp = penwithdraw  		if $samplegendergap & female==`gender' & counterfactual==`count' [aw=weights06], n(3) 
			replace pensionquartile= temp 	if $samplegendergap & female==`gender' & counterfactual==`count'
		}
	}	
		
	putexcel set "$tablepath$outputfile", modify sheet(CounterfGap)
	
	local var = "penwithdraw"	
	* outcome of interest
		
	local varlist = "ones pensionquartile hhfexpergrp HH_edugrp married numk" 
	* subgroups
	
	foreach var2 of varlist `varlist' {
		preserve	
		levelsof `var2', local(var2groups) 
		* values of subgroup
		display "subgroup `var2' takes values `var2groups'"
		
		collapse (mean) `var' = `var' (count) nobs = penwithdraw if $samplegendergap [aw=weights06], by (counterfactual female `var2') 
		* computes the mean pension benefit for each subgroup
		codebook female `var2'
		keep if `var2'!=.
		g group_female=`var2'*1000+female 
		* split subgroups by gender

		drop female `var2'
		reshape wide `var' nobs, i(counterfactual) j(group_female)
		local varlist = ""
		
		foreach group in `var2groups' { 
		* for each value of the subgroup variable
			local suffixmale= `group'*1000+0
			local suffixfemale= `group'*1000+1
			g genderratio_`group'=(`var'`suffixmale'-`var'`suffixfemale')/`var'`suffixmale'  
			*computes the % difference between the mean female and male pension for each subgroup

			local varlist = "`varlist'" + " genderratio_`group'"
		}
		sort counterfactual
		mkmat  `varlist', mat(`var'_`var2')
		restore

	}
	
	putexcel C5  = matrix(penwithdraw_ones')
	putexcel C7  = matrix(penwithdraw_pensionquartile')
	putexcel C11 = matrix(penwithdraw_HH_edugrp')
	putexcel C16 = matrix(penwithdraw_hhfexpergrp')
	putexcel C20 = matrix(penwithdraw_married')
	putexcel C23 = matrix(penwithdraw_numk')


	* labor supply before retirement
	
	foreach var of varlist LFP {
		preserve
		putexcel set "$tablepath$outputfile", modify sheet(CounterfBehav)
			
		collapse `var' if $samplelaborsupply [aw=weights06], by(counterfactual sex)
		g group=sex
		levelsof group, local(var2groups)
		drop sex
		reshape wide `var', i(counterfactual) j(group)	
		
		local varlist =""
		foreach group in `var2groups' {
			local varlist = "`varlist'" + " `var'`group'"

		}
		
		di "`varlist'"
		mkmat  `varlist', mat(`var'mat)			
		restore

	}
	putexcel M6 = matrix(LFPmat')
	*putexcel M12 = matrix(sformmat')

	* experience and pension wealth accumulation
	
	foreach var of varlist hhfexper pension {
	preserve
	putexcel set "$tablepath$outputfile", modify sheet(CounterfBehav)
		
	collapse `var' if $samplestocks [aw=weights06], by(counterfactual sex)
	g group=sex
	levelsof group, local(var2groups)
	drop sex
	reshape wide `var', i(counterfactual) j(group)	
	
	local varlist =""
	foreach group in `var2groups' {
		local varlist = "`varlist'" + " `var'`group'"
	}
	
	di "`varlist'"
	mkmat  `varlist', mat(`var'mat)			

	restore
	}
	putexcel M9 = matrix(hhfexpermat')
	putexcel M12 = matrix(pensionmat')
	
	* non-pension wealth accumulation
	
	foreach var of varlist wealth {
	preserve
	putexcel set "$tablepath$outputfile", modify sheet(CounterfBehav)
	collapse `var' if $samplestocks [aw=weights06], by(counterfactual ones)
	g group=ones
	levelsof group, local(var2groups)
	reshape wide `var', i(counterfactual) j(group)	
	
	local varlist =""
	foreach group in `var2groups' {
		local varlist = "`varlist'" + " `var'`group'"
	}
	
	di "`varlist'"
	mkmat  `varlist', mat(`var'mat)			
	restore
	}
	

	putexcel M15 = matrix(wealthmat')


	* program costs over all retirees
	local l=18
	foreach var of varlist  spsmcostopt   spsfcostopt   fmpgcostopt   mmpgcostopt  pasiscostopt  {
		preserve
		putexcel set "$tablepath$outputfile", modify sheet(CounterfBehav)

		collapse (sum) `var' if $samplecosts [aw=weights06], by(counterfactual ones)
		g group=ones
		levelsof group, local(var2groups)
		reshape wide `var', i(counterfactual) j(group)	
		
		local varlist =""
		foreach group in `var2groups' {
			local varlist = "`varlist'" + " `var'`group'"

		}
		
		di "`varlist'"
		mkmat  `varlist', mat(`var'mat)			
		putexcel M`l' = matrix(`var'mat')
		restore
		local l=`l'+1

	}

foreach var of varlist supplement  {
	preserve
	putexcel set "$tablepath$outputfile", modify sheet(CounterfBehav)

		collapse (sum) `var' if $samplesupplement [aw=weights06], by(counterfactual ones)
		g group=ones
		levelsof group, local(var2groups)
		reshape wide `var', i(counterfactual) j(group)	
		
		local varlist =""
		foreach group in `var2groups' {
		local varlist = "`varlist'" + " `var'`group'"
		}
		
		di "`varlist'"
		mkmat  `varlist', mat(`var'mat)			
	putexcel M`l' = matrix(`var'mat'/1000)
	restore
	local l=`l'+1
	}
	/*
	
	foreach var of varlist  systemcost {
	table counterfactual sex married    if $samplecosts [aw=weights06], c(sum `var')
	}
	
	foreach var of varlist systemcost {
	preserve
	putexcel set "$tablepath/$outputfile", modify sheet(`var')
		putexcel B1 = "Male"
		putexcel C1 = "Female"
	
		collapse (sum) `var' if $samplecosts [aw=weights06], by(counterfactual sex)
		
		g group=sex
		levelsof group, local(var2groups)
		drop sex
		reshape wide `var', i(counterfactual) j(group)	
		
		local varlist =""
		foreach group in `var2groups' {
		local varlist = "`varlist'" + " `var'`group'"
		}
		
		di "`varlist'"
		mkmat counterfactual `varlist', mat(`var'mat)
		putexcel A2 = matrix(`var'mat)
				
		levelsof counterfactual, local(levels)
		local lbe : value label counterfactual
		local ct=2
				
		foreach l of local levels {
			local temp : label `lbe' `l'
			putexcel A`ct'="`temp'"
			local ct=`ct'+1
		}
	restore
	}
	
	* compare counterfactuals
	/*table female  counterfactual  if year==2028 & hhagegroup5==65 & inlist(counterfactual,1,2), c(p10 penwithdraw p25 penwithdraw p50 penwithdraw p75 penwithdraw)
	table female  counterfactual  if year==2028 & hhagegroup5==65 & inlist(counterfactual,1,2), c(p10 hhfexper p25 hhfexper p50 hhfexper p75 hhfexper p90 hhfexper)
	table female  counterfactual  if year==2028 & hhagegroup5==65 & hhfexper>35 & inlist(counterfactual,1,2), c(p10 pension p25 pension p50 pension p75 pension p90 pension)


	table  female married counterfactual  if year==2028 & hhagegroup5==65, c(p10 penwithdraw p25 penwithdraw p50 penwithdraw p75 penwithdraw)
	table  female married counterfactual  if year==2010 & hhagegroup5==65, c(p10 penwithdraw p25 penwithdraw p50 penwithdraw p75 penwithdraw)
	table female married counterfactual  if year==2028 & hhagegroup5==65, c(p10 hhfexper p25 hhfexper p50 hhfexper p75 hhfexper p90 hhfexper)
	table female married counterfactual  if year==2028 & hhagegroup5==65, c(p10 wealth p25 wealth p50 wealth p75 wealth p90 wealth)
	table female married counterfactual  if year==2028 & hhagegroup5==65, c(mean spsmcostopt mean spsfcostopt mean fmpgcostopt mean mmpgcostopt mean pasiscostopt)
	table female married counterfactual  if year>2009 & hhage>59 & hhage<66, c(mean LFP)
	table female married counterfactual  if year>2009 & hhage>65 & hhage<71, c(mean LFP)
	table female married counterfactual  if year>2009 & hhage>59 & hhage<66, c(mean sform)
	table female married counterfactual  if year>2009 & hhage>59 & hhage<66, c(mean tax_paidopt)

	table female married year  if  hhagegroup5==65 & year<2029, by(counterfactual) c(mean penwithdraw)
	table female married year  if  hhagegroup5==60 & year<2029, by(counterfactual) c(mean LFP)
	table female married year  if  hhagegroup5==65 & year<2029, by(counterfactual) c(mean LFP)
	table female married year  if  hhagegroup5==65 & year<2029, by(counterfactual) c(mean wealth)
	 

	table female married hhagegroup5 if year>2009, by(counterfactual) c(mean penwithdraw)

	
	table female married year  if  hhage>59 & hhage<66,by(counterfactual) c(mean spsmcostopt mean spsfcostopt mean fmpgcostopt mean mmpgcostopt mean pasiscostopt)
	/*
	

