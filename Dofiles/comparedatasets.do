 clear
 clear matrix
 clear mata
 
 version 10
 capture log close
 set more off 
  set mem 3000000
/*
 global modelfit 0
 global stderrors 0
 global SRimpacts 1
 global counterfactuals 0

 
 
 if c(username)=="WB452275" {
 
 global path "C:\Users\wb452275\Dropbox\Projects\Genderequity"
 }
 else {
 global path "C:\Users\clement\Dropbox\Projects\Genderequity"
 }
 */

* FILES TO BE COMPARED

if $modelfit==1{
	*global varlist1 "id year numk kidnow married shw sww hexper wexper hfexper hinformexp wfexper winformexp hdeath wdeath hdead wdead sdiv mardur sformh sformw hed wed hage wage female wealth pensionh pensionw earnh earnw lformh lformw lhw lww singlesamp"
	*global varlist1 "id year married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed"
	global varlist1 ///
	                "id year married female wage hage shw sww sformh sformw singlesamp pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw  hed wed sdiv TaxableIncomeh TaxableIncomew "
	*global varlist1 "id year married female wage hage shw sww sformh sformw numk kidnow wed hed savrate pensionw pensionh wealth hexper hfexper wexper wfexper earnh earnw  mpenwithdraw fpenwithdraw"  

	global path1 ///
	"$modelfitpath/../inputs"
	global file1 ///
	"dataformoments.txt"

	global path2 ///
	"$simpath"
	global file2 ///
	"simdata.asc"
	*global varlist2 "id year married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv"
	global varlist2 ///
	"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption hdeath wdeath hdead wdead"
	*global varlist2 ///
	*"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption"
	*global varlist2 "id year clone married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv util consumption hdeath wdeath"
	global outputfile "$resultspath\comparedataMF.dta"
	global year "year>2003 & year<2009"
	global varsort "folio year"
}

if ($SRimpacts==1) {
	global path1 ///
	"$oldsyspath"
	global file1 ///
	"simdata.asc"
	*"simdatanorefo.asc"
	*global varlist1 "id year married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv"
	*global varlist1 ///
*"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption"
	*global varlist1 "id year clone married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv util consumption hdeath wdeath"
	global varlist1 ///
	"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption hdeath wdeath hdead wdead"

	global path2 ///
	"$resultspath"
	global file2 ///
	"simdata.asc"
	*"simdatarefo.asc"

	*global varlist1 "id year married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv"
	*global varlist2 ///
	*"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption"
	*global varlist1 "id year clone married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv util consumption hdeath wdeath"
	global varlist2 ///
	"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption hdeath wdeath hdead wdead"
	global outputfile "$resultspath\comparedataSRI.dta"
	global year "year>2003 & year<2020"
	global varsort "folio clone year"
 }
 
 if ($counterfactuals==1) {
	global path1 ///
	"$resultspath"
	global file1 ///
	"simdata.asc"
	*global varlist1 "id year married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv"
	global varlist1 ///
	"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption"
	*global varlist1 "id year clone married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv util consumption hdeath wdeath"

	global path2 ///
	"$resultspath"
	global file2 ///
	"simdatacounter.asc"
	*global varlist1 "id year married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv"
	global varlist2 ///
	"id year clone type married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv	tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt supplement util consumption"
	*global varlist1 "id year clone married female wage hage shw sww sformh sformw savrate pensionw pensionh wealth numk kidnow hexper hfexper wexper wfexper earnh earnw mpenwithdraw fpenwithdraw hed wed sdiv util consumption hdeath wdeath"
 	global outputfile "$resultspath\comparedataSRI.dta"
	global year "year>2003 & year<2075"
	global varsort "folio clone year"
}

 if ($stderrors==1) {
	global year "year==2005|year==2006"
	}
********************************************************************************

	* INFILE FILE 1
	
	infile $varlist1 using "$path1/$file1", clear

	if $stderrors==1 {
	replace female=female-1
	}
	
	do $path\dofiles\variabledefs.do
	
	keep if $year
	sort $varsort
	
	save "$path1/file1.dta", replace

	save "$path\temp.dta", replace

********************************************************************************

	* INFILE FILE 2
	
	infile $varlist2 using "$path2/$file2", clear
	do $path\dofiles\variabledefs.do
	*codebook folio clone year married female hed wed numk shw sww2 hcovered wcovered earnh earnw wealth if married==1
	
	renvars _all, postfix(_sim)
	
	rename folio_sim folio
	rename year_sim year
	rename clone_sim clone
	sort folio clone year
	save "$path2/file2.dta", replace

********************************************************************************

	* MERGE FILES
	keep if $year
	sort $varsort	
	
	if $modelfit==1 {
	merge m:1 $varsort using "$path\temp.dta"
	}
	else if $stderrors==1 {
	merge m:1 $varsort using "$path\temp.dta"
	}
	else {
	merge 1:1 $varsort using "$path\temp.dta"
	}
	
	keep if _merge ==3  /* note: lose 200 folios */

/*
	foreach var in "numk" "shw" "sww2" "hcovered" "wcovered" "earnh" "earnw" "wealth" {
	replace `var'_sim=. if `var'==.
	replace `var'=. if `var'_sim==.
	}
*/
********************************************************************************


	g hedgrp=int(hed/5)*5
	g wedgrp=int(wed/5)*5 /*XXXalready have hedugrp!*/
	g dfactor=0.94

	*drop if inlist(folio,719939,1040190,1430126,2411979,2418290,2419130)

save "$outputfile", replace
