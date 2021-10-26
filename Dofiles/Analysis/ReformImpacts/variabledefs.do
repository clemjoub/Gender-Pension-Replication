foreach suf in "" "_sim" {
g transfers`suf'=pasiscostopt`suf'+spsmcostopt`suf'+spsfcostopt`suf'+ fmpgcostopt`suf'+ mmpgcostopt`suf'
}

g shwchge=shw_sim-shw
g sww2chge=sww2_sim-sww2
g hcoveredchge=hcovered_sim-hcovered
g wcoveredchge=wcovered_sim-wcovered
g sformhchge=sformh_sim-sformh
g sformwchge=sformw_sim-sformw
g wealthchge=wealth_sim-wealth
g parttime=(sww==2) if sww>0 & sww!=.
g parttime_sim=(sww_sim==2) if sww_sim>0 & sww!=.
g parttimechge=parttime_sim-parttime
g mpenwithdrawchge=mpenwithdraw_sim-mpenwithdraw
g fpenwithdrawchge=fpenwithdraw_sim-fpenwithdraw
g savratechge=savrate_sim-savrate

g mlowpen=mpenwithdraw<1000000
g mmpgqual=(mpenwithdraw==1130472)
g fmpgqual=(fpenwithdraw==1130472)
g mpasisqual=(mpenwithdraw==565236)
g fpasisqual=(fpenwithdraw==565236)
g fpasisqual2=(pasiscostopt>0 & pasiscostopt!=. & (married==1 | female==1)) 

* slight difference with other definition

/* run: br folio clone year female married fpasisqual* pasiscostopt fpenwithdraw  mpenwithdraw if fpasisqual!=fpasisqual2
*/

g mpbsqual=(mpenwithdraw_sim==900000)
g fpbsqual=(fpenwithdraw_sim==900000)

/*
ta mmpgqual if  year==2014 & hagegrp==65
ta fmpgqual if  year==2014 & wagegrp==65
ta mpasisqual if  year==2014 & hagegrp==65
ta fpasisqual if  year==2014 & wagegrp==65
ta mpbsqual if  year==2014 & hagegrp==65
ta fpbsqual if  year==2014 & wagegrp==65
*/

g hcontdens = hfexper/(hexper)
g wcontdens = wfexper/(wexper)

g hcontdensgr=1*(hcontdens<0.25)+2*(hcontdens>=0.25 & hcontdens<0.85)+3*(hcontdens>=0.85)
g wcontdensgr=1*(wcontdens<0.25)+2*(wcontdens>=0.25 & wcontdens<0.85)+3*(wcontdens>=0.85)
g W_edugrp2 = 1*(wed<12)+2*(wed>=12)
g H_edugrp2 = 1*(hed<12)+2*(hed>=12)

g fpenwithdrawtc=fpenwithdraw/1000000
replace fpenwithdrawtc=3 if fpenwithdrawtc>3 & fpenwithdrawtc!=.
g fpenwithdrawtc_sim=fpenwithdraw_sim/1000000
replace fpenwithdrawtc_sim=3 if fpenwithdrawtc_sim>3 & fpenwithdrawtc_sim!=.

g mpenwithdrawtc=mpenwithdraw/1000000
replace mpenwithdrawtc=3 if mpenwithdrawtc>3 & mpenwithdrawtc!=.
g mpenwithdrawtc_sim=mpenwithdraw_sim/1000000
replace mpenwithdrawtc_sim=3 if mpenwithdrawtc_sim>3 & mpenwithdrawtc_sim!=.

g wealthtc = wealth
replace wealthtc=10 if wealthtc>10 & wealth!=.
g wealthtc_sim = wealth_sim
replace wealthtc_sim=10 if wealthtc_sim>10 & wealth_sim!=.

g lowwealth= (wealth<5 & wealth!=.)
g lowwealth_sim= (wealth_sim<5 & wealth_sim!=.)

/*
g hcontdens = hfexper/(hexper)
g wcontdens = wfexper/(wexper)

g hcontdensgr=1*(hcontdens<0.25)+2*(hcontdens>=0.25 & hcontdens<0.85)+3*(hcontdens>=0.85)
g wcontdensgr=1*(wcontdens<0.25)+2*(wcontdens>=0.25 & wcontdens<0.85)+3*(wcontdens>=0.85)
g W_edugrp2 = 1*(wed<12)+2(wed>=12)
g H_edugrp2 = 1*(hed<12)+2(hed>=12)
*/

g hagegrp3=3*int(hage/3) if female==0 | married==1
g wagegrp3=3*int(wage/3) if female==1 | married==1

g hfexpergrp=10*int(hfexper/10) if female==0 | married==1
g wfexpergrp=10*int(wfexper/10) if female==1 | married==1

