
global width=0.33333333333

* Pensions
***********

*Common graph options

local notemacro = `"note("Note 1: Pensions are topcoded at 3 million pesos for legibility")"'
local restriction1 = `"$year"'
local option1 =`"start(-0.05) width($width) fraction bfcolor(none) blcolor(navy) blwidth(medthick)"'
local option2 =`"start(-0.05) width($width) fraction bfcolor(none) blcolor(maroon) blwidth(medthick) "'
local option3 = `"legend(order(1 "without reform" 2 "with reform") pos(2) ring(0) col(1))xline(0.565 1.13 0.9,lc(grey) lp(dash)) xlabel(0.565 "PASIS" 1.13 "MPG" 0.9 "PBS" 0(1)3, angle(90) labs(small)) graphregion(color(white)) bgcolor(white)"'
local yrange = `"ysc(r(0(0.2)0.7)) ylabel(0(0.2)0.7)"'
local xtitlemacro= `"xtitle("Annual pensions (million pesos)")"'

local titlemacro = `"title("Simulated pensions for individuals over 64")"'
**************************************************************
local outfile = "$graphpathpensionsmaried"

local legendmacro =`"legend(off)"'


local outcome = "mpenwithdrawtc"
local restriction2="(female==0) & hage>64"
local restriction3 ="(married==1)"


local xtitlemacro= `"xtitle("") title("Married") ytitle("Without Reform")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp1,replace)
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("With Reform")"'
twoway 	histogram 	`outcome'_sim if $year & `restriction2' & `restriction3',`option2' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp1_sim,replace)

		
		
local outcome = "mpenwithdrawtc"
local restriction2="(female==0) & hage>64"
local restriction3 ="(married==0)"

local xtitlemacro= `"xtitle("") title("Single") ytitle("")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp2,replace)
		
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("")"'
twoway 	histogram `outcome'_sim if $year & `restriction2' & `restriction3',`option2'	`yrange' `option3' `legendmacro' `xtitlemacro' saving(temp2_sim,replace)


local outcome = "fpenwithdrawtc"
local restriction2="(female==1) & wage>64"
local restriction3 ="(married==1)"

local xtitlemacro= `"xtitle("") title("Married") ytitle("")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1'`yrange' `option3' `legendmacro' `xtitlemacro' saving(temp3,replace)
		
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("")"'
twoway 	histogram `outcome'_sim if $year & `restriction2' & `restriction3',`option2' `yrange'	`option3' `legendmacro' `xtitlemacro' saving(temp3_sim,replace)

		

local outcome = "fpenwithdrawtc"
local restriction2="(female==1) & wage>64"
local restriction3 ="(married==0)"
 
local xtitlemacro= `"xtitle("") title("Single") ytitle("")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp4,replace)
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("")"'
twoway 	histogram `outcome'_sim if $year & `restriction2' & `restriction3',`option2'	`yrange' `option3' `legendmacro' `xtitlemacro' saving(temp4_sim,replace)

		


graph combine temp1.gph temp2.gph temp1_sim.gph temp2_sim.gph, r(2) c(2) graphregion(color(white)) saving(tempmen, replace) title(Men)
graph combine temp3.gph temp4.gph  temp3_sim.gph temp4_sim.gph, r(2) c(2) graphregion(color(white)) saving(tempwomen, replace) title(Women)
graph combine tempmen.gph tempwomen.gph, r(1) c(2) graphregion(color(white)) 
graph export `outfile'Figure6.$graphformat, replace $exportoptions

/*
*graph combine temp1.gph temp2.gph temp3.gph temp4.gph temp1_sim.gph temp2_sim.gph temp3_sim.gph temp4_sim.gph, r(2) c(4) graphregion(color(white)) 
*graph export `outfile'.$graphformat, replace $exportoptions

*/



local titlemacro = `"title("Simulated pensions for individuals over 64")"'
**************************************************************

local outfile = "$graphpathpensionsschool"

local legendmacro =`"legend(off)"'
local notemacro = `"note("Note 1: Pensions are topcoded at 3 million pesos for legibility")"'

local outcome = "mpenwithdrawtc"
local restriction2="(female==0) & hage>64"
local restriction3 ="(H_edugrp==1 | H_edugrp==2)"

local xtitlemacro= `"xtitle("") title("Less than HS") ytitle("Without Reform")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp1,replace)
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("With Reform")"'
twoway 	histogram 	`outcome'_sim if $year & `restriction2' & `restriction3',`option2' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp1_sim,replace)


local outcome = "mpenwithdrawtc"
local restriction2="(female==0) & hage>64"
local restriction3 ="(H_edugrp==3 | H_edugrp==4)"

local xtitlemacro= `"xtitle("") title("HS grads") ytitle("")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp2,replace)
		
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("")"'
twoway 	histogram `outcome'_sim if $year & `restriction2' & `restriction3',`option2'	`yrange' `option3' `legendmacro' `xtitlemacro' saving(temp2_sim,replace)


local outcome = "fpenwithdrawtc"
local restriction2="(female==1) & wage>64"
local restriction3 ="(W_edugrp==1 | W_edugrp==2)"

local xtitlemacro= `"xtitle("") title("Less than HS") ytitle("")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1'`yrange' `option3' `legendmacro' `xtitlemacro' saving(temp3,replace)
		
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("")"'
twoway 	histogram `outcome'_sim if $year & `restriction2' & `restriction3',`option2' `yrange'	`option3' `legendmacro' `xtitlemacro' saving(temp3_sim,replace)

		

local outcome = "fpenwithdrawtc"
local restriction2="(female==1) & wage>64"
local restriction3 ="(W_edugrp==3 | W_edugrp==4)"

local xtitlemacro= `"xtitle("") title("HS grads") ytitle("")"'
twoway 	histogram 	`outcome' if $year & `restriction2' & `restriction3',`option1' `yrange' `option3' `legendmacro' `xtitlemacro' saving(temp4,replace)
local xtitlemacro= `"xtitle("Million CLP") title("") ytitle("")"'
twoway 	histogram `outcome'_sim if $year & `restriction2' & `restriction3',`option2'	`yrange' `option3' `legendmacro' `xtitlemacro' saving(temp4_sim,replace)

		


graph combine temp1.gph temp2.gph temp1_sim.gph temp2_sim.gph, r(2) c(2) graphregion(color(white)) saving(tempmen, replace) title(Men)
graph combine temp3.gph temp4.gph  temp3_sim.gph temp4_sim.gph, r(2) c(2) graphregion(color(white)) saving(tempwomen, replace) title(Women)
graph combine tempmen.gph tempwomen.gph, r(1) c(2) graphregion(color(white)) 
graph export `outfile'.$graphformat, replace $exportoptions



