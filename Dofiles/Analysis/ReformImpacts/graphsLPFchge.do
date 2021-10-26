

* LFP
******

sort folio year clone


putexcel set "$tablepath/Tables_Draft", modify sheet(LFPchge)

local restriction2="(hagegrp>49 & hagegrp<75 & female==0 & hage>$Tlb+year-2004)"
local restriction3="(wagegrp>49 & wagegrp<$Tub & female==1 & wage>$Tlb+year-2004)"
local tabstatoptions = "s(mean) f(%5.2f)  save nototal"

local i=5
local var agegrp
tabstat shw shw_sim shwchge if $year & `restriction2' , by(h`var') `tabstatoptions'
tabstatmat A
putexcel D`i' = matrix(A)
tabstat sww2 sww2_sim sww2chge if $year & `restriction3' , by(w`var') `tabstatoptions'
tabstatmat A
putexcel I`i' = matrix(A)
local i=`i'+7


local var edugrp
tabstat shw shw_sim shwchge if $year & `restriction2' , by(H_`var') `tabstatoptions'
tabstatmat A
putexcel D`i' = matrix(A)
tabstat sww2 sww2_sim sww2chge if $year & `restriction3' , by(W_`var') `tabstatoptions'
tabstatmat A
putexcel I`i' = matrix(A)
local i=`i'+6


local var married
tabstat shw shw_sim shwchge if $year & `restriction2' , by(`var') `tabstatoptions'
tabstatmat A
putexcel D`i' = matrix(A)
tabstat sww2 sww2_sim sww2chge if $year & `restriction3' , by(`var') `tabstatoptions'
tabstatmat A
putexcel I`i' = matrix(A)






local titlemacro = `"title("Reform impact on labor force participation") subtitle("Simulated change in age profiles") xtitle("Age")"'
**************************************************************

local outfile = "\LFPchge"
local notemacro = `"note(" ")"'
local outcome = "shwchge"
local xvariable = "hagegrp3"
local outcome2 = "sww2chge"
local xvariable2 = "wagegrp3"
local restriction1 = `"$year"'
local restriction2="(hagegrp>$Tlb & hagegrp<$Tub & female==0 & hage>$Tlb+year-2004)"
local restriction3="(wagegrp>$Tlb & wagegrp<$Tub & female==1 & wage>$Tlb+year-2004)"
local option1 = " "
local option2 = " "
local option3 = `"legend(order(1 "Men" 2 "Women") pos(6) ring(0) col(2))"'
*local yrange = `"ysc(r(0(0.02)0.3)) ylabel(0(0.02)0.3)"'

tabstat shwchge if $year & `restriction2' & married==1, by(hagegrp)
tabstat sww2chge if $year & `restriction3' & married==1, by(wagegrp)
tabstat shwchge if $year & `restriction2' & married==0, by(hagegrp)
tabstat sww2chge if $year & `restriction3' & married==0, by(wagegrp)


capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2', by(`xvariable')
egen mean2 = mean(`outcome2') if $year & `restriction3', by(`xvariable2')
capture drop rank
g rank = 0
replace rank= `xvariable' + 10000 if `restriction2'
replace rank = `xvariable2' + 20000 if `restriction3'
sort rank



twoway 	line  mean1 `xvariable' if $year & `restriction2',`option1' `yrange' || line mean2 `xvariable2' if $year & `restriction3',`option1' `yrange'	`option3'`titlemacro' `notemacro'

		
graph export $graphpath/`outfile'.$graphformat, replace $exportoptions



local titlemacro = `"title("Reform impact on LFP by marital status") subtitle("Simulated change in age profiles")"'
*************************************************************************

	local outfile = "\LFPchgemarital"
	local notemacro = `"note(" ")"'
	local restriction1 = `"$year"'
	local option1 = " "
	local option2 = " "
	local yrange = `"ysc(r(-0.1(0.05)0.05)) ylabel(-0.1(0.05)0.05)"'

	local outcome = "shwchge"
	local xvariable = "hagegrp3"
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
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || line mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange'  `notemacro' `legendmacro' `xtitlemacro' saving(temp1, replace)

			
	local outcome = "sww2chge"
	local xvariable = "wagegrp3"
	local xtitlemacro= `"xtitle("Age")"'
	local restriction2="(wage>$Tlb & wage<$Tub & female==1 & wage>$Tlb+year-2004)"
	local restriction3 ="(married==1)"
	local restriction4 ="(married==0)"
	local legendmacro = `"legend(order(1 "Married wom." 2 "Single wom.") col(2))"'
	capture drop mean1
	capture drop mean2
	egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
	egen mean2 = mean(`outcome') if $year & `restriction2' & `restriction4', by(`xvariable')
	sort `xvariable'
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || line mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange'  `notemacro' `legendmacro' `xtitlemacro' saving(temp2, replace)

			
		
	graph combine temp1.gph temp2.gph, r(1) c(2) `titlemacro' 
	graph export $graphpath/`outfile'.$graphformat, replace $exportoptions
	


local titlemacro = `"title("Reform impact on LFP by schooling") subtitle("Simulated change in age profiles")"'
*************************************************************************

	local outfile = "\LFPchgeschooling"
	local notemacro = `"note(" ")"'
	local restriction1 = `"$year"'
	local option1 = " "
	local option2 = " "
	local yrange = `"ysc(r(-0.1(0.05)0.05)) ylabel(-0.1(0.05)0.05)"'

	local outcome = "shwchge"
	local xvariable = "hagegrp3"
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
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || line mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange' `notemacro' `legendmacro' `xtitlemacro' saving(temp1, replace)

			
	local outcome = "sww2chge"
	local xvariable = "wagegrp3"
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
	twoway 	line 	mean1 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' || line mean2 `xvariable' if $year & `restriction2' & `restriction4',`option1' `yrange' `notemacro' `legendmacro' `xtitlemacro' saving(temp2, replace)

			
		
	graph combine temp1.gph temp2.gph, r(1) c(2) `titlemacro' 
	graph export $graphpath/`outfile'.$graphformat, replace $exportoptions
	
	
	/*

	local titlemacro = `"title("Short- vs. Long-run impacts on labor force participation") subtitle("Men") xtitle("Year")"'
**************************************************************

local outfile = "\LFPchgebyyearmen"
local notemacro = `"note(" ")"'
local outcome = "shwchge"
local xvariable = "year"
local outcome2 = "shwchge"
local xvariable2 = "year"
local outcome3 = "shwchge"
local xvariable3 = "year"
local restriction1 = "year>2006"
local restriction2="(hage>60 & hage<70 & female==0)"
local restriction3="(hage>50 & hage<=60 & female==0)"
local restriction4="(hage>40 & hage<=50 & female==0)"
local option1 = " "
local option2 = " "
local option3 = `"legend(order(1 "Ages 60-70" 2 "Ages 50-60" 3 "Ages 40-50")  col(3))"'
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
		
graph export $graphpath/`outfile'.$graphformat, replace $exportoptions

* LFP
******

local titlemacro = `"title("Short- vs. Long-run impacts on labor force participation") subtitle("Women") xtitle("Year")"'
**************************************************************

local outfile = "\LFPchgebyyearwomen"
local notemacro = `"note(" ")"'
local outcome = "sww2chge"
local xvariable = "year"
local outcome2 = "sww2chge"
local xvariable2 = "year"
local outcome3 = "sww2chge"
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
		
graph export $graphpath/`outfile'.$graphformat, replace $exportoptions

