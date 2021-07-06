
! this subroutine considers all the choices at a given point in time
! and solves for the maximum utility (for each given vector of shocks)

subroutine maxu(o,ct,clone,thisyear,umaxmc,nmc,fixvar,ivart,rvart,yevec,alpha,thetap,Tbar,numz,ktp,totp,mptype,intmat,optch,optsafe,optrisky,earnh,earnw,kidmax,experh,experw,maxearnf,maxearnh,numa1cat,deathprob,cnu,optpenw,optpenm,probvec,maxsafe,mpenwithdraw,fpenwithdraw,reform,actret,nompg,mpg, nopasis, pasis,nobono, bono, nodivrule, pmax,pmin,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt, optc, tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual,numint)
			 
	implicit none

	! external u,TermDiv,TermHDec,TermWDec,ran1

	integer, parameter :: ncoeff=150
	integer, parameter :: numper=75
	integer lastage,lastworkage,counterfactual,numint 
	integer mptype,intmat(75,75), numk, clone, thisyear,hageret,wageret,wageret2
	integer nmc,ktp,totp,tt,wlim,hlim,wllim,hllim,reform,nompg, nopasis, nobono, nodivrule
	integer t1,t2,ct,i,o,kk,optch(nmc,7),numch
	real(8) yevec(nmc*7),yyevec(7),c,deathprob(91,3), mpg, pasis
	real(8) umaxmc(nmc),thetap(3,numz),cnu(76,4),utilm,utilf
	integer fert,hw,ww,Tbar,numz,formh,formw,kidmax
	real(8) experw,experh,safe,risky,ret1,maxsafe(numper)
	real(8) util,term4
	real(8) alpha(ncoeff,ktp),earnh(nmc),earnw(nmc),ttemp1,ttemp2
	real(8) probcofferf,probcofferm,r1(1),r2(1),earnhf,earnhi,earnwf,earnwi
	real(8) earnhfp,earnhip,earnwfp,earnwip,optpenw(nmc),optpenm(nmc), npensionh,npensionw
	real(8) probvec(4),pmax,pmin,actret, nwealth,tao, fexperh, fexperw

	real(8) rvart(6,1,numper)      ! row 1 = wealth - safe
								   ! row 2 = wealth - risky
								   ! row 3 = rate of return safe
								   ! row 4 = rate of return risky
								   ! row 5 = womens pension wealth
								   ! row 6 = mens pension wealth
	integer ivart(9,1,numper)      ! row 1 = children (n)
								   ! row 2 = married
								   ! row 3 = work h
								   ! row 4 = work w
								   ! row 5 = survivh
								   ! row 6 = survivw
								   ! row 7 = marriage dur
								   ! row 8 = formal h
								   ! row 9 = formal w
	integer  fixvar(4,1)           ! row 1 = hed
								   ! row 2 = wed
								   ! row 3 = hage
								   ! row 4 = female indicator 

	integer married,hage,wage,a1,a2
	real(8) nsafe,nrisky,optsafe(nmc),optrisky(nmc),optc(nmc)
	integer female
	real(8) maxearnf,maxearnh,fpenwithdraw,mpenwithdraw,tmpenwithdraw,tfpenwithdraw
	integer hed,wed,first,idum,formlimf,formlimm
	integer numa1cat
	real(8) tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost,term4lag
	real(8) tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt, bono

	idum = 11111 
	nsafe = 0.0d-0
	nrisky = 0.0d-0
	yyevec = 0.0d-0
	maxearnf = 0.0d-0
	maxearnh = 0.0d-0
	tax_paid = 0.0d-0
	spsmcost = 0.0d-0
	spsfcost = 0.0d-0
	fmpgcost = 0.0d-0
	mmpgcost = 0.0d-0
	pasiscost = 0.0d-0           
	util = 0.0d-0
	term4 = 0.0d-0
	term4lag = 0.0d-0

	!---------------------------------------------------------------------------------------------
	!---------------------------------------------------------------------------------------------
	! state variables

	hed = fixvar(1,1)
	wed = fixvar(2,1)
	female = fixvar(4,1)
	married = ivart(2,1,ct)
	wage = ct+15
	hage = min0(fixvar(3,1)+ct-1,90)
	safe = rvart(1,1,ct-1)

	!number of kids
	numk = 0
	numk = sum(ivart(1,1,1:ct))

	!experience
	experh = 0.0d-0
	experw = 0.0d-0
	ttl1:   do tt = 1,(ct-1) 
		if (ivart(3,1,tt).eq.1) then
			 experh = experh + 1.d-0
		end if
		if (ivart(4,1,tt).eq.1) then
			 experw = experw + 1.0d-0
		else if (ivart(4,1,tt).eq.2) then
			 experw = experw + 0.5d-0
		end if
	end do ttl1

	!permanent component of earnings that does not depend on shocks
	earnhfp = alpha(16,mptype)*dexp(alpha(17,1)*(hed-12)+ alpha(18,1)*experh+alpha(19,1)*experh*experh+alpha(20,1)*hed*hed)*1000000.0d-0 
	earnhip = alpha(21,mptype)*dexp(alpha(22,1)*(hed-12)+ alpha(23,1)*experh+alpha(24,1)*experh*experh+alpha(25,1)*hed*hed)*1000000.0d-0 
	earnwfp = alpha(26,mptype)*dexp(alpha(27,1)*(wed-12)+ alpha(28,1)*experw+alpha(29,1)*experw*experw+alpha(30,1)*wed*wed)*1000000.0d-0
	earnwip = alpha(31,mptype)*dexp(alpha(32,1)*(wed-12)+ alpha(33,1)*experw+alpha(34,1)*experw*experw+alpha(35,1)*wed*wed)*1000000.0d-0

	if (reform.eq.1.and.counterfactual.eq.11) then ! equate male and female earning opportunities in this counterfactual
		earnwfp=earnhfp
		earnwip=earnhip
	endif

	maxearnf = max(earnwfp,earnwip)
	maxearnh = max(earnhfp,earnhip)
	  
	wlim = 2            
	hlim = 1
	hllim = 0
	wllim = 0
	formlimf = 0
	formlimm = 0

	! Determine whehter they get an offer from the formal sector
	probcofferm = dexp(alpha(36,1)+alpha(37,1)*dble(hed)+alpha(38,1)*dble(ivart(8,1,(ct-1)))+alpha(39,1)*dble(hage))
	probcofferm = probcofferm/(1.0d-0+probcofferm)
	probcofferf = dexp(alpha(40,1)+alpha(41,1)*dble(wed)+alpha(42,1)*dble(ivart(9,1,(ct-1)))+alpha(43,1)*dble(wage))
	probcofferf = probcofferf/(1.0d-0+probcofferf)

	r1(1) = 0.0d-0
	r2(1) = 0.0d-0
	call ran1(idum,r1,1) 
	call ran1((idum+2),r2,1) 
	if (r1(1).le.probcofferm) then 
		formlimm = 1
	end if
	if (r2(1).le.probcofferf) then
		formlimf = 1
	end if

	!male and female formal sector work experience (needed to qualify for minimum pension benefit guarantee)

	fexperh = 0
	fexperw = 0
	if (ct.gt.1) then
		fexperh = sum(ivart(8,1,1:(ct-1)))   !male formal sector work experience
		fexperw = sum(ivart(9,1,1:(ct-1)))   !woman formal sector work experience
	end if

	!---------------------------------------------------------------------------------------------
	!---------------------------------------------------------------------------------------------
	!--------------------------------------------------------------------------------------------- 
	! Loop through montecarlo shock draws
	kkloop: do kk = 1,nmc

		first = 1
		t1 = (kk-1)*7+1
		t2 = t1+6
		yyevec(1:7) = yevec(t1:t2)

		earnhf = earnhfp*dexp(yyevec(1))
		earnhi = earnhip*dexp(yyevec(2)) 
		earnwf = earnwfp*dexp(yyevec(3)) 
		earnwi = earnwip*dexp(yyevec(4))

		util = 0.0d-0
		ttemp1 = 0.0d-0
		ttemp2 = 0.0d-0
		numch = 0
			
			
		!---------------------------------------------------------------------------------------------
		!---------------------------------------------------------------------------------------------
		!---------------------------------------------------------------------------------------------
		!Loop across all possible choices for a given monte carlo draw

		hhloop: do hw = hlim,hllim,-1             ! husband work
			wwloop: do ww = wlim,wllim,-1              ! wife work (1 is full time, 2 is part time)
				if (hage.ge.lastworkage.and.hw.gt.0) cycle hhloop
				if (wage.ge.lastworkage.and.ww.gt.0) cycle wwloop
				if (married.eq.0.and.female.eq.1.and.hw.gt.0) cycle hhloop
				if (married.eq.0.and.female.eq.0.and.ww.gt.0) cycle wwloop
				formhloop: do formh = 0,formlimm   ! husband formal or informal
					formwloop: do formw = 0,formlimf   ! wife formal or informal
						if (ww.eq.0.and.formw.eq.1) cycle formwloop
						if (hw.eq.0.and.formh.eq.1) cycle formhloop
						if (married.eq.0.and.female.eq.1.and.formh.gt.0) cycle hhloop
						if (married.eq.0.and.female.eq.0.and.formw.gt.0) cycle wwloop
						a1loop: do a1 = 1,numa1cat       !asset type 1
							a2loop: do a2 = 0,0    !asset type 2 (unused) 
								   
								numch = numch+1

								!if (thisyear.ge.2008.and.thisyear.lt.2010) then
								!continue
								!endif
								!if (ct.gt.49.and.a1.gt.1) cycle a1loop
								!if (hw.gt.0) cycle hhloop
								!if (ww.gt.0) cycle wwloop
								!term4lag=term4
		   
								call u(util,ct,o,yyevec,fixvar,ivart,rvart,fert,hw,ww,Tbar,numz,thetap, alpha,ktp,totp,numper,mptype,hage,wage,experh,experw,fexperh,fexperw,a1,a2,intmat, nsafe,nwealth, nrisky,ttemp1,ttemp2,c,formh,formw,married,female,earnhf,earnhi,earnwf,earnwi, maxearnf,maxearnh,term4,numa1cat,nmc,deathprob,cnu,npensionw,npensionh,probvec,maxsafe, tmpenwithdraw,tfpenwithdraw,reform,actret, nompg, mpg, nopasis, pasis, nobono, bono,nodivrule, pmax,pmin,tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost,utilm,utilf, numk,tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual,numint) 

								!if (a1.gt.1.and.term4.gt.term4lag) then
								!print*, "nonmonotonic: period",ct,"state draw",o,"hw",hw,"ww",ww,"formh",formh,"formw",formw,"a1",a1
								!continue
								!endif

								!XXX TESTING
								!
								!if (thisyear.ge.2008.and.thisyear.lt.2010) then
								!write(9999,*) o,ct,clone, mptype, thisyear, wage, hage, ww,hw,formh,formw,a1, married,female,numk, reform, safe, tmpenwithdraw, tfpenwithdraw, tax_paid, nwealth, nsafe, c,ttemp1,ttemp2,npensionw,npensionh, util, utilf, utilm,term4,optch(kk,4), optch(kk,1), optch(kk,2),optch(kk,6),optch(kk,7),yyevec(1:5),tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost
								!endif
			
								! Compare utilities from each subsequent choice and store the highest              
								if (first.eq.1.or.util.gt.umaxmc(kk)) then 
									optch(kk,1) = ww
									optch(kk,2) = hw
									optch(kk,3) = fert
									optch(kk,4) = a1
									optch(kk,5) = a2
									optch(kk,6) = formh 
									optch(kk,7) = formw 
									optsafe(kk) = nsafe
									optrisky(kk) = nrisky
									earnh(kk) = ttemp1
									earnw(kk) = ttemp2
									optc(kk) = c
									umaxmc(kk) = util
									optpenw(kk) = npensionw  !end-of-period pension balance at the optimal choice
									optpenm(kk) = npensionh
									mpenwithdraw = tmpenwithdraw
									fpenwithdraw = tfpenwithdraw
									tax_paidopt=tax_paid/1000000
									spsmcostopt=spsmcost/1000000
									spsfcostopt=spsfcost/1000000
									fmpgcostopt=fmpgcost/1000000
									mmpgcostopt=mmpgcost/1000000
									pasiscostopt=pasiscost/1000000
								end if
								
								if (first.eq.1) first=0

							end do a2loop
						end do a1loop
					end do formwloop
				end do formhloop
			end do wwloop
		end do hhloop
	end do kkloop
end subroutine maxu

!--------------------------------------------------------------------------------------------------
!       This subroutine calculates the utility associated with a given choice, which
!       is designated by the combination (fert,hw,ww,formh,formw) 
!
!       In doing so, it updates the state space for each set of choices
!       being considered
!---------------------------------------------------------------------------------------------------

subroutine u(util,ct,o,yyevec,fixvar,ivart,rvart,fert,hw,ww,Tbar,numz,thetap, &
        alpha,ktp,totp,numper,mptype,hage,wage,experh,experw,fexperh,fexperw,a1,a2,intmat, &
        nsafe, nwealth, nrisky,earnh,earnw,c,formh,formw,married,female,earnhf,earnhi,earnwf,earnwi, &
        maxearnf,maxearnh,term4,numa1cat,nmc,deathprob,cnu,npensionw,npensionh,probvec,maxsafe, &
        mpenwithdraw,fpenwithdraw,reform,actret, nompg,mpg, nopasis, pasis, nobono, bono, nodivrule, pmax,pmin,&
        tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost,utilm,utilf, numk, tao,hageret,wageret,wageret2,&
	lastage,lastworkage, counterfactual,numint)

	implicit none

	! external ProbDiv,TermDiv,TermHDec,TermWDec,tax
	integer, parameter :: ncoeff=150
	integer, parameter :: nbracket=8

	integer hageret,wageret,wageret2,intmat(75,75),female,nmc,reform, nompg, nopasis, nobono, nodivrule, numint
	integer numa1cat
	integer ktp,totp,numz,numper,hage,wage
	real(8) experh,experw,earnhf,earnhi,earnwf,earnwi,deathprob(91,3), bono
	real(8) thetap(3,numz),probvec(4),pmax,pmin,actret, supplement
	real(8) ProbDiv,cnu(76,4),maxsafe(numper)
	real(8) rvart(6,1,numper)      ! row 1 = wealth - safe
								   ! row 2 = wealth - risky
								   ! row 3 = rate of return safe
								   ! row 4 = blank 
								   ! row 5 = womens pension wealth 
								   ! row 6 = mens pension wealth 
	integer ivart(9,1,numper)      ! row 1 = children (n)
								   ! row 2 = married
								   ! row 3 = work h
								   ! row 4 = work w
								   ! row 5 = survivh
								   ! row 6 = survivw
								   ! row 7 = marriage dur
								   ! row 8 = formal h
								   ! row 9 = formal w
	integer fixvar(4,1)            ! row 1 = hed
								   ! row 2 = wed
								   ! row 3 = hage
								   ! row 4 = female 
	real(8) nrvart(6,1,numper)     ! row 1 = wealth - safe
								   ! row 2 = wealth - risky
								   ! row 3 = rate of return safe
								   ! row 4 = rate of return risky
								   ! row 5 = womens pension wealth 
								   ! row 6 = mens pension wealth 
	integer nivart(9,1,numper)     ! row 1 = children (n)
								   ! row 2 = married
								   ! row 3 = work h
								   ! row 4 = work w
								   ! row 5 = survivh
								   ! row 6 = survivw
								   ! row 7 = marriage dur
								   ! row 8 = formal h
								   ! row 9 = formal w
	integer nfixvar(4,1)           ! row 1 = hed
								   ! row 2 = wed
								   ! row 3 = hage
								   ! row 4 = female
	real(8) maxearnf,maxearnh
	real(8) tao
	integer hed,wed,fert,hw,ww,mptype,lastage,lastworkage
	integer married,exph,expw
	integer ct,o,ctt,numkid
	real(8) yyevec(7),safe,risky
	integer a1,a2,counterfactual
	real(8) term4
	real(8) util,alpha(ncoeff,ktp)
	integer Tbar,numk
	real(8) earnh,earnw
	real(8) c, utilf,utilm
	real(8) pdiv,phdec,pwdec,cmin,mpg
	real(8) beta,ret1,nwealth,nsafe,nrisky
	integer ind65h,ind60w,formh,formw,ind70h
	real(8) probfert,temp22,pensionw,pensionh,npensionw,npensionh
	integer mswitch,fswitch,mreentry,freentry,fqualearly,mqualearly
	real(8) earlycut(20)
	real(8) ubequest
	
	!! variables related to taxes
	integer :: tax_toggle,tempw,tempm
	real(8) bracket_min(nbracket)    !bracket lower bounds 
	real(8) bracket_max(nbracket)   !bracket upper bounds
	real(8) bracket_rate(nbracket)   !marginal tax rates
	real(8) bracket_adj(nbracket)    !taxes due are obtained by multiplying
									 !the taxable income with the marginal tax rate
									 !and subtracting the backet-specific adjustment
	real(8) taxable_wages,tax_paid,taxw,taxm,tax,pasis,fpenwithdraw,mpenwithdraw,faminc,incthreshold
	real(8) pasisqualm,pasisqualf,mpgqualf,mpgqualm,spsm,spsf,fexperh,fexperw
	real(8) temp111m,temp111f,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost, consgrid(numa1cat)

	! Initialize

	data bracket_min(:) /0 ,4.909d-0,10.910d-0,18.184d-0,25.458d-0,32.732d-0,43.643d-0,54.554d-0 /&
	bracket_max(:) /4.909d-0,10.910d-0,18.184d-0,25.458d-0,32.732d-0,43.643d-0,54.554d-0,9999999/ &
	bracket_rate(:) /0,0.05d-0,0.1d-0,0.15d-0,0.25d-0,0.32d-0,0.37d-0,0.4d-0 / &
	bracket_adj(:) /0,0.245d-0,0.791d-0,1.7d-0,4.246d-0,6.537d-0,8.719d-0,10.356d-0 /  

	beta = alpha(15,1)
	cmin = alpha(85,1)*1000000.0d-0    ! set a minimum consumption level
	ret1 = 0.05d-0
	tax_toggle = 2
	util = -9999999999.0d-0
	ctt=ct+1
	pensionw = 0.0d-0
	pensionh = 0.0d-0
	pasisqualf = 0
	pasisqualm = 0
	mpgqualf = 0
	mpgqualm = 0
	temp111f = 0.0d-0
	temp111m = 0.0d-0
	c = 0
	fpenwithdraw = 0.0d-0
	mpenwithdraw = 0.0d-0
	npensionh = 0.0d-0
	npensionw = 0.0d-0
	safe = 0.0d-0
	risky = 0.0d-0
	mswitch = 0
	fswitch = 0
	mreentry = 0
	freentry = 0
	ind65h = 0
	ind70h = 0
	ind60w = 0
	earnh = 0.0d-0
	earnw = 0.0d-0
	taxable_wages = 0.0d-0
	taxm = 0.0d-0
	taxw = 0.0d-0
	tax_paid = 0.0d-0
	fqualearly = 0
	mqualearly = 0

	earlycut(20) = 13682941.0d-0
	earlycut(19) = 14021373.0d-0
	earlycut(18) = 14650655.0d-0
	earlycut(17) = 14943005.0d-0
	earlycut(16) = 15221433.0d-0
	earlycut(15) = 15486603.0d-0
	earlycut(14) = 15739146.0d-0
	earlycut(13) = 15979663.0d-0
	earlycut(12) = 15979663.0d-0
	earlycut(11) = 16208727.0d-0
	earlycut(10) = 16426882.0d-0
	earlycut(9) = 16634650.0d-0
	earlycut(8) = 16832524.0d-0
	earlycut(7) = 17020975.0d-0
	earlycut(6) = 17200452.0d-0
	earlycut(5) = 17371383.0d-0
	earlycut(4) = 17534174.0d-0
	earlycut(3) = 17689214.0d-0
	earlycut(2) = 17836870.0d-0
	earlycut(1) = 17977495.0d-0

	! Create some interim variables
	hed = fixvar(1,1)
	wed = fixvar(2,1)

	wage = ct+15
	married = ivart(2,1,ct)

	if (ct.gt.1) then
		safe = rvart(1,1,(ct-1))    
		risky = rvart(2,1,(ct-1))
		pensionw = rvart(5,1,(ct-1))
		pensionh = rvart(6,1,(ct-1))
	end if

	if (hage.ge.65) then
	   ind65h = 1
	end if
	if (hage.ge.lastworkage) then
	   ind70h = 1
	end if
	if (wage.ge.60) then
	   ind60w = 1
	end if

	if (ct.gt.1) then
	  if ((ivart(8,1,(ct-1)).eq.1.and.(formh.eq.0.and.hw.eq.1)).or.&
		   ((ivart(8,1,(ct-1)).eq.0.and.ivart(3,1,(ct-1)).eq.1).and.formh.eq.1)) then
			mswitch = 1
	  end if
	  if ((ivart(9,1,(ct-1)).eq.1.and.(formw.eq.0.and.ww.ge.1)).or.&
		   ((ivart(9,1,(ct-1)).eq.0.and.ivart(4,1,(ct-1)).ge.1).and.formw.eq.1)) then
			fswitch = 1
	  end if
	  if (ivart(3,1,(ct-1)).eq.0.and.hw.eq.1) then  ! determine if reentering the workforce
			mreentry = 1
	  end if
	  if (ivart(4,1,(ct-1)).eq.0.and.ww.ge.1) then  ! determine if reentering the workfoce
			freentry = 1
	  end if
	end if

	!--------------------------------------------------------------------------------------------
	!--------------------------------------------------------------------------------------------
	!--------------------------------------------------------------------------------------------
	!      calculate male and female earnings and taxes given their work choices 

	if (hw.eq.1.and.formh.eq.1) then
	  earnh = (1-tao)*earnhf  
	  pensionh = pensionh + tao*earnhf    ! 10% of formal earnings goes to pension
	  taxm = (1-tao)*earnh   !do not get taxed on contributions
	end if
	if (hw.eq.1.and.formh.eq.0) then
	  earnh = earnhi
	end if
	if (ww.eq.1.and.formw.eq.1) then
	  earnw = (1-tao)*earnwf 
	  pensionw = pensionw + tao*earnwf      ! 10% of formal earnings goes to pension
	  taxw = (1-tao)*earnw  !do not get taxed on contributions
	end if
	if (ww.eq.2.and.formw.eq.1) then
	  earnw = (1-tao)*0.5d-0*earnwf   !earn half of a full-time worker
	  pensionw = pensionw + tao*0.5d-0*earnwf      ! 10% of formal earnings goes to pension
	  taxw = (1-tao)*0.5d-0*earnw   !do not get taxes on contributions
	end if
	if (ww.eq.1.and.formw.eq.0) then
	  earnw = earnwi 
	end if
	if (ww.eq.2.and.formw.eq.0) then
	  earnw = 0.5d-0*earnwi 
	end if

	if (reform.eq.1.and.(counterfactual.eq.4).and.hw.eq.1.and.formh.eq.1.and.formw.eq.0) then !counterfactual: hb contributes to non-contributing wife's account
	  pensionh = pensionh-(tao/2)*earnhf ! pensionsplit counterfactual: half of husband's contribution goes to wife
	  pensionw = pensionw+(tao/2)*earnhf
	endif

	if (reform.eq.1.and.(counterfactual.eq.5.or.counterfactual.eq.10).and.hw.eq.1.and.formh.eq.1.and.formw.eq.0) then !counterfactual: hb contribute to non-contributing wife's account
	  earnh = earnh - tao*earnhf
	  pensionw = pensionw+tao*earnhf
	endif

	if (married.eq.1) then
	  taxable_wages = taxm+taxw
	else if (married.eq.0.and.female.eq.1) then
		 taxable_wages = taxw
	else if (married.eq.0.and.female.eq.0) then
		 taxable_wages=taxm
	end if

	if (taxable_wages.gt.0.0d-0) then
	   tax_paid = tax((taxable_wages+0.05d-0*safe),tax_toggle,bracket_min,bracket_max,bracket_rate,bracket_adj)
	end if

	!----------------------------------------------------------------------------------------------------------
	!       Compute programmed withdrawal pension benefits and  next period's balance
	!   assets earn interest

	if (nmc.gt.1) then 
		pensionh = pensionh*(1+actret+yyevec(5))
		pensionw = pensionw*(1+actret+yyevec(5))
	else
		pensionh = pensionh*(1+actret)   !no shock for simulation when using actual return
		pensionw = pensionw*(1+actret) 
	end if

	safe=safe*(1+ret1)

	!   determine if male or female qualify for early retirement

	if (hage.ge.45.and.hage.le.64) then
	   if (pensionh.gt.earlycut(hage-44)) then
		  mqualearly = 1
	   end if
	end if
	if (wage.ge.45.and.wage.le.59) then
		  if (pensionw.gt.earlycut(wage-44)) then
			fqualearly = 1
		  end if
	end if


	tempw = wage-34
	tempm = hage-34
	if (tempw.lt.1) then
		tempw = 1
	else if (tempw.gt.76) then
		tempw = 76
	end if 
	if (tempm.lt.1) then
		tempm = 1
	else if (tempm.gt.76) then
		tempm = 76
	end if 

	if (married.eq.1) then
		if (hage.ge.hageret.or.mqualearly.eq.1) then
		  mpenwithdraw=pensionh/cnu(tempm,4)  !determine programmed withdrawal amount
		end if
		if (wage.ge.wageret.or.fqualearly.eq.1) then
		  fpenwithdraw=pensionw/cnu(tempw,3)  !determine programmed withdrawal amount
		end if
		npensionh = dmax1((pensionh-mpenwithdraw),0.0d-0) !npensionh is remaining balance after pension withdrawal
		npensionw = dmax1((pensionw-fpenwithdraw),0.0d-0)
	end if
	if (married.eq.0.and.female.eq.1) then
		if (wage.ge.wageret.or.fqualearly.eq.1) then
		  fpenwithdraw=pensionw/cnu(tempw,3)  !determine programmed withdrawal amount
		end if
		npensionw = dmax1((pensionw-fpenwithdraw),0.0d-0)
	end if
	if (married.eq.0.and.female.eq.0) then
		if (hage.ge.hageret.or.mqualearly.eq.1) then
		  mpenwithdraw=pensionh/cnu(tempm,2)  !determine programmed withdrawal amount
		end if
		npensionh = dmax1((pensionh-mpenwithdraw),0.0d-0)
	end if

	!----------------------------------------------------------------------------------------------------------------       
	!   Add minimum pension benefits to the baseline pension withdrawals

	faminc=(earnh+mpenwithdraw)*(1-female*(1-married))+(earnw+fpenwithdraw)*(1-(1-female)*(1-married))+safe*0.05 ! family income for purpose of eligibility to SPS
	incthreshold=12.0d-0*162084.0d-0 !annual income of the 60th percentile Chilean household in 2009 (VII Encuesta de presupuestos familiares, ine.cl)
	spsm = 0.0d-0   ! solidarity pension benefit
	spsf = 0.0d-0   ! solidarity pension benefit

	if (reform.eq.0) then
	
		if (married.eq.1) then
		
			! A. MPG - married couple
			if (wage.ge.wageret.and.fexperw.ge.20) then
				mpgqualf = 1   ! wife qualifies
				fmpgcost=dmax1(mpg-fpenwithdraw,0.0d-0)
				fpenwithdraw = dmax1(mpg,fpenwithdraw)  !supplement pension
			end if
			if (hage.ge.hageret.and.fexperh.ge.20) then
				mpgqualm = 1   ! husband qualifies
				mmpgcost=dmax1(mpg-mpenwithdraw,0.0d-0)
				mpenwithdraw = dmax1(mpg,mpenwithdraw)  !supplement pension
			end if

			! B. PASIS - married couple
			! only can get pasis if neither spouse qualifies for mpg and average income for household is less than pasis (XXX age restriction?)
			if (mpgqualf.eq.0.and.mpgqualm.eq.0.and.(faminc/2.0d-0).lt.pasis.and.&
					(wage.ge.wageret2.or.hage.ge.hageret)) then
				pasisqualf = 1  !only one spouse can qualify for pasis
				pasiscost=dmax1(pasis-fpenwithdraw,0.0d-0)
				fpenwithdraw = dmax1(pasis,fpenwithdraw)   !only one spouse can get pasis: give it to woman
			end if
			
		end if  ! married.eq.1

		if (married.eq.0.and.female.eq.1) then   ! get the MPG pension if beyond retirement age, PASIS if older than 65
			
			! A. MPG - single woman
			if (fexperw.ge.20.and.wage.ge.wageret) then
				mpgqualf = 1
				fmpgcost=dmax1(mpg-fpenwithdraw,0.0d-0)
				if (fmpgcost.gt.0.0d-0) then
					continue
				endif
				fpenwithdraw = dmax1(mpg,fpenwithdraw)
			end if 
			
			! B. PASIS - single woman
			if (mpgqualf.eq.0.and.(faminc).lt.pasis.and.wage.ge.wageret2) then
				pasisqualf = 1
				pasiscost=dmax1(pasis-fpenwithdraw,0.0d-0)
				fpenwithdraw = dmax1(pasis,fpenwithdraw)
			end if
			
		end if
		
		if (married.eq.0.and.female.eq.0) then   ! get the PASIS or MPG pension if beyond retirement age

		! A. MPG - single man
			if (fexperh.ge.20.and.hage.ge.hageret) then
				mpgqualm = 1
				mmpgcost=dmax1(mpg-mpenwithdraw,0.0d-0)
				mpenwithdraw = dmax1(mpg,mpenwithdraw)
			end if

		! B. PASIS - single man
			if (mpgqualm.eq.0.and.(faminc).lt.pasis.and.hage.ge.hageret) then
				pasisqualm = 1
				pasiscost=dmax1(pasis-mpenwithdraw,0.0d-0)
				mpenwithdraw = dmax1(pasis,mpenwithdraw)
			end if
			
		end if
		
	endif

	if (reform.eq.1) then

		!C. Solidarity pension/APS

		if (married.eq.1) then

			! -  married couple
			if (wage.ge.wageret2.and.fpenwithdraw.lt.pmax) then
				temp111f = fpenwithdraw/pmax
				spsf = pmin*(1.0d-0-temp111f)
				if (faminc.gt.incthreshold) spsf=0.0
				if (faminc.le.incthreshold) then
					continue
				endif
				fpenwithdraw = fpenwithdraw+spsf
				spsfcost=spsf
			end if
			  
			if (hage.ge.hageret.and.mpenwithdraw.lt.pmax) then
				temp111m = mpenwithdraw/pmax
				spsm = pmin*(1.0d-0-temp111m)
				if (faminc.gt.incthreshold) spsm=0.0
				mpenwithdraw = mpenwithdraw+spsm
				spsmcost=spsm			
			end if		
		
		end if

		! - single individuals
		if (married.eq.0) then

			!- female
			if (female.eq.1.and.wage.ge.wageret2.and.fpenwithdraw.lt.pmax) then
				temp111f = fpenwithdraw/pmax
				spsf = pmin*(1.0d-0-temp111f)
				if (faminc.gt.incthreshold) spsf=0.0
				fpenwithdraw = fpenwithdraw+spsf
				spsfcost=spsf
			end if

			! - male
			if (female.eq.0.and.hage.ge.hageret.and.mpenwithdraw.lt.pmax) then
				temp111m = mpenwithdraw/pmax
				spsm = pmin*(1.0d-0-temp111m)
				if (faminc.gt.incthreshold) spsm=0.0
				mpenwithdraw = mpenwithdraw+spsm
				spsmcost=spsm
			end if
				
		end if
	
	endif  

	!----------------------------------------------------------------------------------------------------------------
	!      calculate and store today's consumption and tomorrow's non-pension assets

	nwealth = safe + earnh + earnw + mpenwithdraw + fpenwithdraw - tax_paid !      nwealth is disposable wealth this period
	nwealth=dmax1(nwealth,cmin*(1+married)) !implement consumption floor as a wealth floor (hh can't save out of the wealth floor given my consumption decision set up)
	if (numa1cat.eq.30) then
		consgrid=(/0.50,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,12.0,15.0,20.0,25.0,30.0,50.0,75.0,100.0/)
	elseif (numa1cat.eq.15) then
		consgrid=(/0.50,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.0,15.0,30.0/)
	elseif (numa1cat.eq.10) then
		consgrid=(/0.50,1.0,1.5,2.0,3.0,5.0,7.0,10.0,15.0,30.0/)
	endif 


	c=consgrid(a1)*1000000*(1+married)

	if (c.gt.nwealth) then
		if(a1.eq.1) then            
			c=nwealth
		elseif (a1.gt.1) then
			util=-99999999999999.0d-0
			goto 100
		endif
	endif

	nsafe=nwealth-c        
	if ((ct+1).le.numper) then
		nsafe = dmin1(nsafe,maxsafe(ct+1)) 
	end if


	!----------------------------------------------------------------------------------------------------------------
	!       Compute probabilities of divorce, fertility and death for next period

	pdiv = 0.0d-0
	if (married.eq.1.and.ct.lt.45) then   !no divorce after woman age 60
		pdiv = ProbDiv(ncoeff,alpha,hage,wage,ktp,hed)
	end if
	probvec(1) = pdiv  !store divorce probability

	temp22 = 0.0d-0
	probfert = 0.0d-0
	!can have a child if a married couple or a single female
	if ((female.eq.1.or.married.eq.1).and.ct.le.25) then
		temp22 = dexp(alpha(44,1)+alpha(45,1)*(ct+15)+alpha(46,1)*dble(married) + alpha(47,1)*dble(numk)) + &
				  alpha(67,1)*wed + alpha(68,1)*dble(married)*dble(numk)
		probfert = 1.0d-0 - (temp22 / (1.0d-0+temp22) )
	end if
	probvec(2) = probfert

	if (hage.ge.90) then   !people only live to age 90
		phdec = 1.0d-0
	elseif (hage.ge.0) then
		phdec = deathprob((hage-19),2)
	else
		phdec = 0.0d-0
	end if
	if (wage.eq.90) then
		pwdec = 1.0d-0
	elseif (wage.ge.0) then
		pwdec = deathprob((wage-19),3)
	else
		pwdec = 0.0d-0
	end if
	probvec(3)=phdec
	probvec(4)=pwdec

	call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm, ubequest)

	!------------------------------------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------------------------------------
	!------------------------------------------------------------------------------------------------------------------
	!        Compute EMAX terms  if NOT last period

	term4=0.0d-0

	if ((ct+15).lt.Tbar.and.(hage.lt.90.or.married==1.or.female==1)) then 

	!  create temporary state variable arrays that can be modified for each possible choice (nxxxvart)
		nrvart(:,1,ct) = rvart(:,1,ct)
		nivart(:,1,ct:ct+1) = ivart(:,1,ct:ct+1)
		nfixvar = fixvar
		
		nivart(1,1,ct) = 0 
		nivart(3,1,ct) = hw 
		nivart(4,1,ct) = ww 
		nivart(8,1,ct) = formh
		nivart(9,1,ct) = formw
	  
	   call contvalue(married,female,numper,term4,thetap,reform,nodivrule,nivart,nrvart,nfixvar,ct,o,phdec,pwdec,pdiv,probfert,nsafe,npensionh,npensionw,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,numz,totp,intmat,pasis,bono,cnu,lastage,lastworkage,numint)
	   
	   util = util+beta*term4
	   
	end if !ct+15 

	! --------------------------------------------------------------------------------------------------------------------
	! --------------------------------------------------------------------------------------------------------------------
	! --------------------------------------------------------------------------------------------------------------------
	!    output choice utilities

	!if (reform.eq.0.and.o==2045.and.(ct==74)) then
	!open(unit=590,position='append',file='emaxhat.txt')
	!write(590,*) reform,o, ct, hw,formh,ww,formw,a1,util,term4,wsum1,married,female,hed,wed,numk,experw,experh,fexperw,fexperh,safe,pensionw,pensionh
	!close(590)
	!endif

	if ((o.eq.781.or.o.eq.960).and.ct==25) then
	continue
	endif

100 end subroutine u 

function ProbDiv(ncoeff,alpha,hage,wage,ktp,hed)

	implicit none

	integer ncoeff,hage,wage,ktp
	real(8) alpha(ncoeff,ktp),xb
	real(8) ProbDiv
	integer hed,agediff2,wage2,agediff

	agediff = hage-wage
	agediff2 = (hage-wage)**2
	wage2= wage**2
	! specify prob of divorce
	xb = alpha(48,1)+ alpha(49,1)*dble(wage)+alpha(50,1)*dble(wage2)+ &
		  alpha(51,1)*dble(hed)+alpha(52,1)*dble(agediff)+alpha(53,1)*dble(agediff2)
	ProbDiv = 1-dexp(xb)/(1+dexp(xb))

end function ProbDiv

!-------------------------------------------------------------------------------------------------
! The probabilities of dying are taken from life tables and gets the programmed
! withdrawal amounts

Subroutine ProbDec(deathprob,cnu)
	
	!---------------------------------------------------------------
	! This subroutine takes in age of a person, sex and race and spits
	! out a probability of death;
	! age = 1,2, etc. in years
	! sex = 1 Male; 2 Female
	!--------------------------------------------------
	IMPLICIT NONE
	REAL(8) deathprob(91,3)  !first col age, second men, third women
	integer age ,kk,vv1
	real(8) maleprob, femaleprob,cnu(76,4)
	!--------------------- map actual age in years into intervals
	!

	maleprob  = 0.0d-0
	femaleprob = 0.0d-0
	open(unit=444,file='./inputs/Lifetable_men.txt',position='rewind')
	open(unit=445,file='./inputs/Lifetable_women.txt',position='rewind')
	do kk = 20,110
		read(fmt=*,unit=444) age,maleprob
		read(fmt=*,unit=445) age,femaleprob
		deathprob((kk-19),1)=age
		deathprob((kk-19),2)=maleprob
		deathprob((kk-19),3)=femaleprob
	end do

	close(unit=444)
	close(unit=445)

	open(unit=446,file='inputs/CNU.txt',position='rewind')
	read(fmt=*,unit=446) (cnu(vv1,1),vv1=1,76) 
	read(fmt=*,unit=446) (cnu(vv1,2),vv1=1,76) 
	read(fmt=*,unit=446) (cnu(vv1,3),vv1=1,76) 
	read(fmt=*,unit=446) (cnu(vv1,4),vv1=1,76) 
	close(unit=446)
	
end subroutine ProbDec

function tax(taxable,tax_toggle,bracket_min,bracket_max,bracket_rate,bracket_adj)

	implicit none

	double precision, parameter:: tax_rate = 0.2d-0       ! rate used in progressive taxation

	integer, parameter :: nbracket=8
	integer :: tax_toggle
	integer :: br               !index for tax brackets
	double precision :: tax     !taxes due
	double precision :: taxable !taxable income
	double precision :: bracket_min(nbracket) 
	double precision :: bracket_max(nbracket) 
	double precision :: bracket_adj(nbracket) 
	double precision :: bracket_rate(nbracket) 

	tax = 0.0d-0
	if (tax_toggle==1) then
		tax = taxable*tax_rate
		elseif (tax_toggle==2) then
			do br = 1,nbracket
				if ((bracket_min(br).lt.taxable).and.(bracket_max(br).ge.taxable)) then
					tax = bracket_rate(br)*taxable-bracket_adj(br)
					exit
				endif
			end do
	endif

end function tax
