
* Coverage
******


sort folio year clone


*putexcel set "$tablepath/Tables_Draft", modify sheet(Covchge)
putexcel set "$tablepath\Tables_Draft", modify sheet(CondCovchge)


*local restriction2="(hagegrp>49 & hagegrp<65 & female==0 & hage>$Tlb+year-2004)"
*local restriction3="(wagegrp>49 & wagegrp<65 & female==1 & wage>$Tlb+year-2004)"
local restriction2="(hagegrp>49 & hagegrp<65 & female==0 & hage>$Tlb+year-2004 & shw==1)"
local restriction3="(wagegrp>49 & wagegrp<65 & female==1 & wage>$Tlb+year-2004 & sww2==1)"
local tabstatoptions = "s(mean) f(%5.2f)  save nototal"

table hfexpergrp H_edugrp married if $year & `restriction2' & inrange(hfexpergrp,0,40) , c(mean sformhchge)
table wfexpergrp W_edugrp married if $year & `restriction3' & inrange(wfexpergrp,0,40) , c(mean sformwchge)



local i=5
local var agegrp
tabstat sformh sformh_sim sformhchge if $year & `restriction2' , by(h`var') `tabstatoptions'
tabstatmat A
putexcel D`i' = matrix(A)
tabstat sformw sformw_sim sformwchge if $year & `restriction3' , by(w`var') `tabstatoptions'
tabstatmat A
putexcel I`i' = matrix(A)
local i=`i'+5


local var edugrp
tabstat sformh sformh_sim sformhchge if $year & `restriction2' , by(H_`var') `tabstatoptions'
tabstatmat A
putexcel D`i' = matrix(A)
tabstat sformw sformw_sim sformwchge if $year & `restriction3' , by(W_`var') `tabstatoptions'
tabstatmat A
putexcel I`i' = matrix(A)
local i=`i'+6


local var married
tabstat sformh sformh_sim sformhchge if $year & `restriction2' , by(`var') `tabstatoptions'
tabstatmat A
putexcel D`i' = matrix(A)
tabstat sformw sformw_sim sformwchge if $year & `restriction3' , by(`var') `tabstatoptions'
tabstatmat A
putexcel I`i' = matrix(A)


/*


local titlemacro = `"title("Reform impact on formal work") subtitle("Simulated change in age profiles") xtitle("Age")"'
**************************************************************

local outfile = "\Coveragechge"
local notemacro = `"note("Fraction of sample that accepts a formal job offer not conditional on working")"'
local outcome = "sformhchge"
local xvariable = "hage"
local outcome2 = "sformwchge"
local xvariable2 = "wage"
local restriction1 = `"$year"'
local restriction2="(hagegrp>$Tlb & hagegrp<$Tub & female==0 & hage>$Tlb+year-2004)"
local restriction3="(wagegrp>$Tlb & wagegrp<$Tub & female==1 & wage>$Tlb+year-2004)"
local option1 = " "
local option2 = " "
local option3 = `"legend(order(1 "Men" 2 "Women") pos(6) ring(0) col(2))"'
*local yrange = `"ysc(r(0(0.02)0.3)) ylabel(0(0.02)0.3)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2', by(`xvariable')
egen mean2 = mean(`outcome2') if $year & `restriction3', by(`xvariable2')
capture drop rank
g rank = 0
replace rank= `xvariable' + 10000 if `restriction2'
replace rank = `xvariable2' + 20000 if `restriction3'
sort rank

twoway 	line 	mean1 `xvariable' if $year & `restriction2',`option1' `yrange' || ///
		line	mean2 `xvariable2' if $year & `restriction3',`option1' `yrange' ///
		`option3'`titlemacro' `notemacro' ///
		
graph export $graphpath/`outfile'.$graphformat, replace



local titlemacro = `"title("Reform impact on formal work by marital status") subtitle("Simulated change in age profiles")"'
*************************************************************************

	local outfile = "\Coveragechgemarital"
	local notemacro = `"note("Fraction of sample that accepts a formal job offer not conditional on working")"'
	local restriction1 = `"$year"'
	local option1 = " "
	local option2 = " "
	local yrange = `"ysc(r(-0.1(0.05)0.15)) ylabel(-0.1(0.05)0.15)"'

	local outcome = "sformhchge"
	local xvariable = "hage"
	local xtitlemacro= `"xtitle("Age")"'
	local restriction2="(hage>$Tlb & hage<$Tub & female==0 & hage>$Tlb+year-2004)"
	local restriction3 ="(married==1)"
	local restriction4 ="(married==0)"
	local legendmacro = `"legend(order(1 "Married men" 2 "Single men")  col(2))"'
	capture drop mean1
	capture drop mean2
	egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
	egen mean2 = mean(`outcome') if $year & `restriction2' & `restriction4', by(`xvariable')
	sort `xvariable'
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || ///
			line	mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange'  ///
			  `legendmacro' `xtitlemacro' saving(temp1, replace)
			
	local outcome = "sformwchge"
	local xvariable = "wage"
	local xtitlemacro= `"xtitle("Age")"'
	local restriction2="(wage>$Tlb & wage<$Tub & female==1 & wage>$Tlb+year-2004)"
	local restriction3 ="(married==1)"
	local restriction4 ="(married==0)"
	local legendmacro = `"legend(order(1 "Married women" 2 "Single women") col(2))"'
	capture drop mean1
	capture drop mean2
	egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
	egen mean2 = mean(`outcome') if $year & `restriction2' & `restriction4', by(`xvariable')
	sort `xvariable'
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || ///
			line	mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange'  ///
			`legendmacro' `xtitlemacro' saving(temp2, replace)
			
		
	graph combine temp1.gph temp2.gph, r(1) c(2) `titlemacro' `notemacro'  
	graph export $graphpath/`outfile'.$graphformat, replace
	


local titlemacro = `"title("Reform impact on formal work by schooling") subtitle("Simulated change in age profiles")"'
*************************************************************************

	local outfile = "\Coveragechgeschooling"
local notemacro = `"note("Fraction of sample that accepts a formal job offer not conditional on working")"'
	local restriction1 = `"$year"'
	local option1 = " "
	local option2 = " "
	local yrange = `"ysc(r(-0.1(0.05)0.15)) ylabel(-0.1(0.05)0.15)"'

	local outcome = "sformhchge"
	local xvariable = "hage"
	local xtitlemacro= `"xtitle("Age")"'
	local restriction2="(hage>$Tlb & hage<$Tub & female==0 & hage>$Tlb+year-2004)"
	local restriction3 ="(H_edugrp==1 | H_edugrp==2)"
	local restriction4 ="(H_edugrp==3 | H_edugrp==4)"
	local legendmacro = `"legend(order(1 "Men < HS grad" 2 "Men >= HS grad") rows(2) col(1))"'
	capture drop mean1
	capture drop mean2
	egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
	egen mean2 = mean(`outcome') if $year & `restriction2' & `restriction4', by(`xvariable')
	sort `xvariable'
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || ///
			line	mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange'  ///
			  `legendmacro' `xtitlemacro' saving(temp1, replace)
			
	local outcome = "sformwchge"
	local xvariable = "wage"
	local xtitlemacro= `"xtitle("Age")"'
	local restriction2="(wage>$Tlb & wage<$Tub & female==1 & wage>$Tlb+year-2004)"
local restriction3 ="(W_edugrp==1 | W_edugrp==2)"
	local restriction4 ="(W_edugrp==3 | W_edugrp==4)"
	local legendmacro = `"legend(order(1 "Women < HS grad" 2 "Women >= HS grad") rows(2) col(1))"'
	capture drop mean1
	capture drop mean2
	egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
	egen mean2 = mean(`outcome') if $year & `restriction2' & `restriction4', by(`xvariable')
	sort `xvariable'
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || ///
			line	mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange'  ///
			  `legendmacro' `xtitlemacro' saving(temp2, replace)
			
		
	graph combine temp1.gph temp2.gph, r(1) c(2) `titlemacro' `notemacro'
	graph export $graphpath/`outfile'.$graphformat, replace

/*
local titlemacro = `"title("Short- vs. Long-run impacts on formal work") subtitle("Men") xtitle("Year")"'
**************************************************************

local outfile = "\covchgebyyearmen"
local notemacro = `"note(" ")"'
local outcome = "sformhchge"
local xvariable = "year"
local outcome2 = "sformhchge"
local xvariable2 = "year"
local outcome3 = "sformhchge"
local xvariable3 = "year"
local restriction1 = "year>2006"
local restriction2="(hage>60 & hage<70 & female==0)"
local restriction3="(hage>50 & hage<=60 & female==0)"
local restriction4="(hage>40 & hage<=50 & female==0)"
local option1 = " "
local option2 = " "
local option3 = `"legend(order(1 "Ages 60-70" 2 "Ages 50-60" 3 "Ages 40-50") pos(6) ring(0) col(2))"'
*local yrange = `"ysc(r(0(0.02)0.3)) ylabel(0(0.02)0.3)"'
local xrange = `"xsc(r(2007(2)2020)) xline(2008,lp(dash)) xlabel(2008 "Reform" 2007(2)2020, angle(90))"'
capture drop mean1
capture drop mean2
capture drop mean3
egen mean1 = mean(`outcome') if $year & `restriction2', by(`xvariable')
egen mean2 = mean(`outcome2') if $year & `restriction3', by(`xvariable2')
egen mean3 = mean(`outcome3') if $year & `restriction4', by(`xvariable3')
capture drop rank
g rank = 0
replace rank= `xvariable' + 10000 if `restriction2'
replace rank = `xvariable2' + 20000 if `restriction3'
replace rank = `xvariable3' + 30000 if `restriction4'
sort rank

twoway 	line 	mean1 `xvariable' if `restriction1' & `restriction2',`option1' `yrange' || ///
		line	mean2 `xvariable2' if `restriction1' & `restriction3',`option1' `yrange' || ///
		line	mean3 `xvariable3' if `restriction1' & `restriction4',`option1' `yrange' ///
		`option3'`titlemacro' `notemacro' `xrange' ///
		
graph export $graphpath/`outfile'.$graphformat, replace



local titlemacro = `"title("Short- vs. Long-run impacts on formal work") subtitle("Women") xtitle("Year")"'
**************************************************************

local outfile = "\covchgebyyearwomen"
local notemacro = `"note(" ")"'
local outcome = "sformwchge"
local xvariable = "year"
local outcome2 = "sformwchge"
local xvariable2 = "year"
local outcome3 = "sformwchge"
local xvariable3 = "year"
local restriction1 = "year>2006"
local restriction2="(wage>60 & wage<70 & female==1)"
local restriction3="(wage>50 & wage<=60 & female==1)"
local restriction4="(wage>40 & wage<=50 & female==1)"
local option1 = " "
local option2 = " "
local option3 = `"legend(order(1 "Ages 60-70" 2 "Ages 50-60" 3 "Ages 40-50") pos(6) ring(0) col(2))"'
*local yrange = `"ysc(r(0(0.02)0.3)) ylabel(0(0.02)0.3)"'
local xrange = `"xsc(r(2007(2)2020)) xline(2008,lp(dash)) xlabel(2008 "Reform" 2007(2)2020, angle(90))"'
capture drop mean1
capture drop mean2
capture drop mean3
egen mean1 = mean(`outcome') if $year & `restriction2', by(`xvariable')
egen mean2 = mean(`outcome2') if $year & `restriction3', by(`xvariable2')
egen mean3 = mean(`outcome3') if $year & `restriction4', by(`xvariable3')
capture drop rank
g rank = 0
replace rank= `xvariable' + 10000 if `restriction2'
replace rank = `xvariable2' + 20000 if `restriction3'
replace rank = `xvariable3' + 30000 if `restriction4'
sort rank

twoway 	line 	mean1 `xvariable' if `restriction1' & `restriction2',`option1' `yrange' || ///
		line	mean2 `xvariable2' if `restriction1' & `restriction3',`option1' `yrange' || ///
		line	mean3 `xvariable3' if `restriction1' & `restriction4',`option1' `yrange' ///
		`option3'`titlemacro' `notemacro' `xrange' ///
		
graph export $graphpath/`outfile'.$graphformat, replace
