

* LFP
******

local titlemacro = `"title("Male labor force participation") subtitle("Simulated with and without reform") xtitle("Age")"'
**************************************************************

local outfile = "$graphpath\LFPmale"
local notemacro = `"note(" ")"'
local outcome = "shw"
local xvariable = "hage"
local restriction1 = `"year"'
local restriction2="(hagegrp>54 & hagegrp<75 & female==0 & hage>$Tlb+year-2004)"
local option1 = "graphregion(color(white)) bgcolor(white) "
local option2 = " "
local option3 = `"legend(order(1 "without reform" 2 "with reform") pos(2) ring(0) col(1))"'
*local yrange = `"ysc(r(0(0.2)1)) ylabel(0(0.2)1)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2', by(`xvariable')
sort `xvariable'
twoway line mean1 mean2 `xvariable' if $year & `restriction2' $weights,`option1' `yrange'  `option3' `titlemacro' `notemacro' 

graph export `outfile'.$graphformat, replace $exportoptions


local titlemacro = `"title("Female labor force participation") subtitle("Simulated with and without reform") xtitle("Age")"'
**************************************************************

local outfile = "$graphpath\LFPfemale"
local notemacro = `"note(" ")"'
local outcome = "sww2"
local xvariable = "wage"
local restriction1 = `"year"'
local restriction2="(wagegrp>$Tlb & wagegrp<$Tub & female==1 & wage>$Tlb+year-2004)"
local option1 = "graphregion(color(white)) bgcolor(white) "
local option2 = " "
local option3 = `"legend(order(1 "without reform" 2 "with reform") pos(2) ring(0) col(1))"'
*local yrange = `"ysc(r(0(0.2)1)) ylabel(0(0.2)1)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' $weights,`option1' `yrange' `option3'`titlemacro' `notemacro'

graph export `outfile'.$graphformat, replace $exportoptions

/*

local titlemacro = `"title("Labor force participation by marital status") subtitle("Simulated with and without reform")"'
*************************************************************************

local outfile = "$graphpath\LFPmarital"
local notemacro = `"note(" ")"'
local restriction1 = `"year"'
local option1 = " "
local option2 = " "
*local yrange = `"ysc(r(0(0.2)1)) ylabel(0(0.2)1)"'

local outcome = "shw"
local xvariable = "hage"
local xtitlemacro= `"xtitle("Married men")"'
local restriction2="(hage>$Tlb & hage<$Tub & female==0 & hage>$Tlb+year-2004)"
local restriction3 ="(married==1)"
local legendmacro =`"legend(off)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		 `notemacro' `legendmacro' `xtitlemacro' saving(temp1, replace)
		

local outcome = "shw"
local xvariable = "hage"
local xtitlemacro= `"xtitle("Single men")"'
local restriction2="(hage>$Tlb & hage<$Tub & female==0 & hage>$Tlb+year-2004)"
local restriction3 ="(married==0)"
local legendmacro = `"legend(off)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		 `notemacro' `legendmacro' `xtitlemacro' saving(temp2, replace)
		

local outcome = "sww2"
local xvariable = "wage"
local xtitlemacro= `"xtitle("Married women")"'
local restriction2="(wage>$Tlb & wage<$Tub & female==1 & wage>$Tlb+year-2004)"
local restriction3 ="(married==1)"
local legendmacro =`"legend(off)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		`notemacro' `legendmacro' `xtitlemacro' saving(temp3, replace)
		

local outcome = "sww2"
local xvariable = "wage"
local xtitlemacro= `"xtitle("Single women")"'
local restriction2="(wage>$Tlb & wage<$Tub & female==1  & wage>$Tlb+year-2004)"
local restriction3 ="(married==0)"
local legendmacro = `"legend(order(1 "without reform" 2 "with reform") pos(1) ring(0) col(1))"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		 `notemacro' `legendmacro' `xtitlemacro' saving(temp4, replace)

		
graph combine temp1.gph temp2.gph temp3.gph temp4.gph, r(2) c(2) `titlemacro' ycommon
graph export `outfile'.$graphformat, replace $exportoptions



local titlemacro = `"title("Labor force participation by schooling") subtitle("Simulated with and without reform")"'
*************************************************************************

local outfile = "$graphpath\LFPschooling"
local notemacro = `"note(" ")"'
local restriction1 = `"year"'
local option1 = " "
local option2 = " "
*local yrange = `"ysc(r(0(0.2)1)) ylabel(0(0.2)1)"'

local outcome = "shw"
local xvariable = "hage"
local xtitlemacro= `"xtitle("Men < HS grad")"'
local restriction2="(hage>$Tlb & hage<$Tub & female==0  & hage>$Tlb+year-2004)"
local restriction3 ="(H_edugrp==1 | H_edugrp==2)"
local legendmacro =`"legend(off)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		 `notemacro' `legendmacro' `xtitlemacro' saving(temp1, replace)
		

local outcome = "shw"
local xvariable = "hage"
local xtitlemacro= `"xtitle("Men >= HS grad")"'
local restriction2="(hage>$Tlb & hage<$Tub & female==0  & hage>$Tlb+year-2004)"
local restriction3 ="(H_edugrp==3 | H_edugrp==4)"
local legendmacro = `"legend(off)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		 `notemacro' `legendmacro' `xtitlemacro' saving(temp2, replace)
		

local outcome = "sww2"
local xvariable = "wage"
local xtitlemacro= `"xtitle("Women < HS grad")"'
local restriction2="(wage>$Tlb & wage<$Tub & female==1  & wage>$Tlb+year-2004)"
local restriction3 ="(W_edugrp==1 | W_edugrp==2)"
local legendmacro =`"legend(off)"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		`notemacro' `legendmacro' `xtitlemacro' saving(temp3, replace)
		

local outcome = "sww2"
local xvariable = "wage"
local xtitlemacro= `"xtitle("Women >= HS grad")"'
local restriction2="(wage>$Tlb & wage<$Tub & female==1  & wage>$Tlb+year-2004)"
local restriction3 ="(W_edugrp==3 | W_edugrp==4)"
local legendmacro = `"legend(order(1 "without reform" 2 "with reform") pos(1) ring(0) col(1))"'
capture drop mean1
capture drop mean2
egen mean1 = mean(`outcome') if $year & `restriction2' & `restriction3', by(`xvariable')
egen mean2 = mean(`outcome'_sim) if $year & `restriction2' & `restriction3', by(`xvariable')
sort `xvariable'
twoway 	line 	mean1 mean2 `xvariable' if $year & `restriction2' & `restriction3',`option1' `yrange' ///
		 `notemacro' `legendmacro' `xtitlemacro' saving(temp4, replace)

		
graph combine temp1.gph temp2.gph temp3.gph temp4.gph, r(2) c(2) `titlemacro' ycommon
graph export `outfile'.$graphformat, replace $exportoptions

