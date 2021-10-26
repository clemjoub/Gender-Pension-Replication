clear
clear matrix
clear mata

 * path
	 if c(username)=="WB452275" {
		global path "C:\Users\WB452275\Documents\GitHub\gender-pension\"
	 global tablepath = "C:\Users\WB452275\Documents\GitHub\gender-pension-draft\Tables"
}
	 else if c(username)=="clement" {
		global path "C:\Users\clement\Dropbox\Projects\Genderequity"
	 }
		else{
		global path "C:\Genderequity"
	 global tablepath = "C:\Genderequity\Tables"
}
 
global dopath "$path\Dofiles\Analysis\ReformImpacts"

global oldsyspath = "$path\Counterfactuals\0\outputs"
global resultspath = "$path\Counterfactuals\1\outputs"
*global oldsyspath = "$path\reformimp\outputs"
global graphpath= "$resultspath\..\reformimp\"

global weights ="[aw=weights06]"



* Generate file with data and simulations
global modelfit 0
global stderrors 0
global SRimpacts 1
global counterfactuals 0
global loaddata 1

if $loaddata ==1 {
	do $path\Dofiles\comparedatasets.do
}

use "$resultspath\comparedataSRI.dta", clear

sort folio
capture drop _merge
merge m:1 folio using "$path\WorkingDatasets\EPSweights.dta"

keep if _merge==3


cd "$resultspath"

global fittables="$tablepath\fittables"

do $dopath\variabledefs.do


global year = "year>2009 & year<2016"
global Tub = 75
global Tlb = 36
global graphformat = "png"
global exportoptions = "width(1600) height(1200)"

********************************************************************************

* Transfer receipts
/*
table hfexpergrp H_edugrp married if $year & hagegrp>64 & female==0 & hage>$Tlb+year-2004 & inrange(hfexpergrp,0,40) , c(mean mpasisqual)
table hfexpergrp H_edugrp married if $year & hagegrp>64 & female==0 & hage>$Tlb+year-2004 & inrange(hfexpergrp,0,40) , c(mean mmpgqual)
table hfexpergrp H_edugrp married if $year & hagegrp>64 & female==0 & hage>$Tlb+year-2004 & inrange(hfexpergrp,0,40) , c(mean mpbsqual)

table hfexpergrp H_edugrp if $year & hagegrp>64 & female==0 & hage>$Tlb+year-2004 & inrange(hfexpergrp,0,40) , c(mean mpenwithdraw mean mpenwithdraw_sim)
table hfexpergrp H_edugrp married if $year & hagegrp>64 & female==0 & hage>$Tlb+year-2004 & inrange(hfexpergrp,0,40) , c(mean transfers mean transfers_sim)

table wfexpergrp W_edugrp married if $year & wagegrp>64 & female==1 & wage>$Tlb+year-2004 & inrange(wfexpergrp,0,40) , c(mean fpasisqual)
table wfexpergrp W_edugrp married if $year & wagegrp>64 & female==1 & wage>$Tlb+year-2004 & inrange(wfexpergrp,0,40) , c(mean fmpgqual)
table wfexpergrp W_edugrp married if $year & wagegrp>64 & female==1 & wage>$Tlb+year-2004 & inrange(wfexpergrp,0,40) , c(mean fpbsqual)

table wfexpergrp W_edugrp  if $year & wagegrp>64 & female==1 & wage>$Tlb+year-2004 & inrange(wfexpergrp,0,40) , c(mean fpenwithdraw mean fpenwithdraw_sim)
table wfexpergrp W_edugrp married if $year & wagegrp>64 & female==1 & wage>$Tlb+year-2004 & inrange(wfexpergrp,0,40) , c(mean transfers mean transfers_sim)
*/
*drop if $year==0

/*
g penwithdraw = fpenwithdrawtc 
replace penwithdraw = mpenwithdrawtc if female==0
g penwithdraw_sim = fpenwithdrawtc_sim
replace penwithdraw_sim = mpenwithdrawtc_sim if female==0
g LFP= sww2
replace LFP = shw if female==0
g LFP_sim= sww2_sim
replace LFP_sim = shw_sim if female==0
g hhage=wage
replace hhage = hage if female==0

table female if year>2008 & year<2020 & hhage>64 & hhage<90, c(p10 penwithdraw p25 penwithdraw p50 penwithdraw p75 penwithdraw p90 penwithdraw)
table female if year>2008 & year<2020 & hhage>59 & hhage<66, c(mean LFP)

table female if year>2008 & year<2020 & hhage>64 & hhage<90, c(p10 penwithdraw_sim p25 penwithdraw_sim p50 penwithdraw_sim p75 penwithdraw_sim p90 penwithdraw_sim)
table female if year>2008 & year<2020 & hhage>59 & hhage<66, c(mean LFP_sim)

*/

do $dopath\graphspensions3.do
*do $dopath\graphscoverage.do
*do $dopath\graphscoveragechge.do
*do $dopath\graphsLFP.do
*do $dopath\graphsLPFchge.do
*do $dopath\graphssavings.do
do $dopath\graphsinequality.do



/*


*Savings
********

latabstat wealth wealth_sim if $year, by(hhagegrp2) s(mean) tf($tablepath/SRImpacts.tex) replace


local restriction2="(wage>45 & wage<70 & female==1)"
local restriction3 ="(W_edugrp==1 | W_edugrp==2)"
tabstat sww2  if `restriction2' & `restriction3' & numk==0, by(wage) s(n mean)

g mmpg_dum = mmpgcostopt>0 & mmpgcostopt!=.
ta mmpg_dum year if hage>64

g mmpg_dum_sim = mmpgcostopt_sim>0 & mmpgcostopt_sim!=.
ta mmpg_dum_sim year if hage_sim>64

g qual=hfexper>19
latabstat mmpg_dum qual pensionh if hage>64 & female==0, by(year) 

*/

*rm "$resultspath\comparedataSRI.dta"
