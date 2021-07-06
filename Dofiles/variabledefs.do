*codebook hage shw sformh pensionh hexper hfexper earnh mpenwithdraw hed  if married==0 & female==1
*codebook wage sww sformw pensionw wexper wfexper earnw fpenwithdraw wed  if married==0 & female==0
recode _all (-99=.)
rename id folio

g hagegrp=5*int(hage/5) 
g wagegrp=5*int(wage/5)
g hagegrp2=10*int((hage+5)/10)-5 
g wagegrp2=10*int((wage+5)/10)-5
g hhagegrp2=hagegrp2 
replace hhagegrp2=wagegrp2 if female==1

g W_edugrp = 1*(wed<8)+2*(wed>=8 & wed<12)+3*(wed>=12 & wed<16)+4*(wed>=16)
g H_edugrp = 1*(hed<8)+2*(hed>=8 & hed<12)+3*(hed>=12 & hed<16)+4*(hed>=16)

g HH_edugrp = H_edugrp 
replace HH_edugrp = W_edugrp if female==1

replace earnh =. if shw ==0 
replace earnw =. if sww ==0
replace hed =. if married==0 & female==1
g sww2 =int((sww +1)/2)

replace wealth=wealth/1000000 
replace earnh=earnh/1000000 
replace earnw=earnw/1000000

g hcovered  = sformh 
replace hcovered  =. if shw ==0
g wcovered  = sformw 
replace wcovered  =. if sww ==0
g bothcovered=(shw==1 & sformh==1 & sww>0 & sformw==1)
g bothuncovered=(shw==1 & sformh==0 & sww>0 & sformw==0)
g twoearners=(shw==1 & sww>0) 
g noearner=(shw==0 & sww==0)

* sample composition

replace wage=. if married ==0 & female ==0
replace hage =. if married ==0 & female ==1
replace wed =. if married ==0 & female ==0
replace hed =. if married ==0 & female ==1
replace sww =. if married ==0 & female ==0
replace sww2 =. if married ==0 & female ==0
replace shw =. if married ==0 & female ==1
replace sformw =. if married ==0 & female ==0
replace sformh =. if married ==0 & female ==1
replace wexper =. if married ==0 & female ==0
replace hexper =. if married ==0 & female ==1
replace wfexper =. if married ==0 & female ==0
replace hfexper =. if married ==0 & female ==1

