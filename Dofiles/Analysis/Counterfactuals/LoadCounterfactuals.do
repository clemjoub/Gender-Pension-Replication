
infile $varlist1 using "$resultspath\0\outputs\simdata.asc", clear

g counterfactual=0
save "$resultspath\AppendedCounterfactuals", replace

use "$resultspath\AppendedCounterfactuals", clear

forvalues c=1(1)13 {  	

if (`c'==2 | `c'==9) {
continue
}
display "appending counterfactual `c'"
infile $varlist1 using "$resultspath\\`c'\outputs\simdata.asc", clear


g counterfactual=`c'
append using "$resultspath\AppendedCounterfactuals"
save "$resultspath\AppendedCounterfactuals", replace
}
save "$resultspath\AppendedCounterfactuals", replace

