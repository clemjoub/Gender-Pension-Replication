*Inequality
***********


recode earnh earnw wealth mpenwithdraw fpenwithdraw earnh_sim earnw_sim wealth_sim mpenwithdraw_sim fpenwithdraw_sim (.=0)


**household pensions
g hhpenwithdraw=mpenwithdraw+fpenwithdraw
g hhpenwithdraw_sim=mpenwithdraw_sim+fpenwithdraw_sim

*household income (no asset returns)
g hhinc = earnh +earnw + mpenwithdraw + fpenwithdraw
g hhinc_sim = earnh_sim +earnw_sim  + mpenwithdraw_sim + fpenwithdraw_sim

*household income (with asset returns)
g hhinc2 = earnh +earnw + wealth*0.05+ mpenwithdraw + fpenwithdraw
g hhinc2_sim = earnh_sim +earnw_sim + wealth*0.05 + mpenwithdraw_sim + fpenwithdraw_sim

*pension of the head of household
g penwithdraw=mpenwithdraw
replace penwithdraw=fpenwithdraw if female==1
g penwithdraw_sim=mpenwithdraw_sim
replace penwithdraw_sim=fpenwithdraw_sim if female==1

save temp, replace
use temp, clear

*sample restriction
global population = "$year & hhagegrp2>64"


local graphoptions= `" graphregion(color(white))"'
***
local outfile = "$graphpathlorenzpens"

local titlemacro = `"title("Lorenz curve of invidual pensions") subtitle("Individuals over 65")"'

preserve
stack penwithdraw  penwithdraw_sim if $population, into(penwithdraw) clear
rename _stack reform
label var penwithdraw "Pension level"

label define reform 1 "No reform" 2 "Reform"
label var reform reform
replace reform=reform-1
glcurve penwithdraw, by(reform) split lorenz legend(label(1 "No reform") label(2 "Reform")) ytitle("Lorenz Curve") title("Individual Pension Income") `graphoptions' saving(tempind,replace)
restore


***
local titlemacro = `"title("Lorenz curve of total household pensions") subtitle("household heads over 65")"'


preserve
stack hhpenwithdraw  hhpenwithdraw_sim if $population, into(hhpenwithdraw) clear
rename _stack reform
replace reform=reform-1
glcurve hhpenwithdraw, by(reform) split lorenz legend(label(1 "No reform") label(2 "Reform"))  ytitle("Lorenz Curve") title("Household Pension Income") `graphoptions' saving(temphh,replace)
restore


graph combine tempind.gph temphh.gph, r(1) c(2) graphregion(color(white)) 
graph export `outfile'Figure7.$graphformat, replace $exportoptions


/*
***
local outfile = "$graphpath\lorenzhhinc1"
local titlemacro = `"title("Lorenz curves of household income (excluding asset returns)") subtitle("household head over 65")"'


preserve
stack hhinc  hhinc_sim if $population, into(hhinc) clear
rename _stack reform
replace reform=reform-1
glcurve hhinc, by(reform) split  lorenz legend(label(1 "No reform") label(2 "Reform")) `graphoptions'
graph export `outfile'.$graphformat, replace $exportoptions

restore

***
local outfile = "$graphpath\lorenzhhinc2"
local titlemacro = `"title("Lorenz curves of household income (including asset returns)") subtitle("household head over 65")"'


preserve
stack hhinc2  hhinc2_sim if $population, into(hhinc) clear
rename _stack reform
replace reform=reform-1
glcurve hhinc, by(reform) split lorenz legend(label(1 "No reform") label(2 "Reform")) `graphoptions' 
graph export `outfile'.$graphformat, replace $exportoptions

restore

***
local outfile = "$graphpath\lorenzhhcons"
local titlemacro = `"title("Lorenz curves of household consumption") subtitle("household head over 65")"'

preserve
stack consumption  consumption_sim if $population, into(consumption) clear
rename _stack reform
replace reform=reform-1
glcurve consumption, by(reform) split  lorenz legend(label(1 "No reform") label(2 "Reform")) `graphoptions'
graph export `outfile'.$graphformat, replace $exportoptions

restore
*/

local outfile = "$graphpathlorenz"


grc1leg temp1.gph temp2.gph, r(2) c(2) leg(temp2.gph) graphregion(color(white)) 
graph export `outfile'.$graphformat, replace $exportoptions

/*
local outfile = "$graphpath\lorenzhhutil"
local titlemacro = `"title("Lorenz curves of household flow utility") subtitle("household head over 65")"'

preserve
stack util  util_sim if $population, into(util) clear
rename _stack reform
replace reform=reform-1
glcurve util, by(reform) split  lorenz legend(label(1 "No reform") label(2 "Reform")) `graphoptions'
graph export `outfile'.$graphformat, replace $exportoptions

restore
