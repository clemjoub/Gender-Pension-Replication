
subroutine paramsinterp(wageparams,typeparams,prefparams,segmentparams,alpha,ncoeff,ktp,Tbar,theta,totp,alabel)

	integer,parameter:: nbins=10
	integer,parameter:: binsize=500000
	integer,parameter:: nwdraws=100
	integer,parameter:: ntypesim=4
	integer,parameter:: lastage = 90
	integer,parameter:: lastworkage = 70
	integer::hedvec(ntypesim),wedvec(ntypesim),cohortvec(ntypesim),smarriedvec(ntypesim),sfemalevec(ntypesim),hagevec(ntypesim),wagevec(ntypesim)
	real(8)::temp111,typeprob(ntypesim),ret
	integer:: ncoeff,ktp,totp,ind65h,ind60w,wage,hage,Tbar,ct,mptype,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,parttime,fulltime
	real(8):: nsafe,ubequest,ubequest1,ubequest2,ubequest1k,alpha(ncoeff,ktp),npensionh,npensionw,beta,phdec,pwdec,yyevec(7),util,c,utilf,utilm
	real(8):: MUC,MUC1,MUC2,MUC3,MUC1k,MUbq,MULh,MULw,MULw3, MRSh,MRSw,MRSwk,MRSw2,MRSw3,cbase,switchcosth,switchcostw,util1,util1k,util2,util3,util4,util5,nonpech,nonpecw
	real(8):: entrycosth,entrycostw,nsafebase
	real(8):: theta(75,3,106),wdraws(2,2,nwdraws)
	integer:: experh,experw,wageparams, typeparams,prefparams,segmentparams,bequestparams
	integer::i,wage_sim(2,2,2,nbins),bin,sec,gen,typ,draw,hed,wed,binubound,binlbound,idum
	real(8):: earn(2,2,2),earnp(2,2,2),wage_dist(2,2,2,nbins),var(2,2),v(1,1)
	character*50::alabel(300)
	real(8):: mswitchcost, fswitchcost, mreentrycost, freentrycost
	real(8):: deathprob(91,3),cnu(76,4)
	real :: leisuredraws(2,1000)
	real(8) :: leisdraws(2,1000)
	idum=dint(-1.0d-0*abs(548997))

	if (typeparams.eq.1) then
		hedvec=(/4,4,12,16/)
		wedvec=(/4,12,12,16/)
		smarriedvec=(/1,0,1,0/)
		hagevec=(/35,45,65,35/)
		wagevec=(/35,45,65,35/)
		sfemalevec=(/0,1,0,1/)
		cohortvec=(/5,4,2,5/)
		do i=1,ntypesim
			temp111 =   alpha(75,1)+&
				alpha(78,1)*dble(smarriedvec(i))+&
				alpha(79,1)*dble(hagevec(i))*dble(1-sfemalevec(i))+&
				alpha(79,1)*dble(wagevec(i))*dble(sfemalevec(i))+&
				alpha(76,1)*dble(wedvec(i))*dble(smarriedvec(i))+ &
				alpha(76,1)*dble(wedvec(i))*dble(sfemalevec(i))*(1.0d-0-dble(smarriedvec(i))) + &
				alpha(77,1)*dble(hedvec(i))*dble(smarriedvec(i)) + &
				alpha(77,1)*dble(hedvec(i))*(1.0d-0-dble(sfemalevec(i)))*(1.0d-0-dble(smarriedvec(i)))
			typeprob(i) = dexp(temp111)/(1.0d-0+dexp(temp111))
		 enddo
		 write(*,*) "type probabilities"
		 write(*,*)
		 write(*,*) "constant",alpha(75,1)
		 write(*,*) "married",alpha(78,1)
		 write(*,*) "sampled age", alpha(79,1)
		 write(*,*) "female schooling", alpha(76,1)
		 write(*,*) "male schooling", alpha(77,1)
		 write(*,*) 
		 write(*,fmt='(A10,4I6)') "hed",hedvec
		 write(*,fmt='(A10,4I6)') "wed",wedvec
		 write(*,fmt='(A10,4I6)') "married",smarriedvec
		 write(*,fmt='(A10,4I6)') "hagevec",hagevec
		 write(*,fmt='(A10,4I6)') "wagevec",wagevec
		 write(*,fmt='(A10,4I6)') "female",sfemalevec
		 write(*,fmt='(A10,4I6)') "cohortvec",cohortvec
		 write(*,fmt='(A10,4f6.3)') "prob typ1",typeprob
		 write(*,*) 
	endif
	if (wageparams.eq.1) then
		
		write(*,*) "Returns to schooling"
		do hed=0,20,5
			write(*,*)
			ret=alpha(17,1)+2*hed*alpha(20,1)
			write(*,*) "ret to school, formal (man w",hed,"yrs of school):",ret        
			ret=alpha(22,1)+2*hed*alpha(25,1)
			write(*,*) "returns to schooling, informal sector (man with",hed,"yrs of schooling):",ret        
			ret=alpha(27,1)+2*hed*alpha(30,1)
			write(*,*) "returns to schooling, formal sector (woman with",hed,"yrs of schooling):",ret        
			ret=alpha(32,1)+2*hed*alpha(35,1)
			write(*,*) "returns to schooling, informal sector (woman with",hed,"yrs of schooling):",ret        
	   enddo 
		
		write(*,*)  "Returns to experience:"
		do hed=0,45,10
			write(*,*)
			ret=alpha(18,1)+2*hed*alpha(19,1)
			write(*,*) "returns to exp, formal sector (man with",hed,"yrs of exp):",ret        
			ret=alpha(23,1)+2*hed*alpha(24,1)
			write(*,*) "returns to exp, informal sector (man with",hed,"yrs of exp):",ret        
			ret=alpha(28,1)+2*hed*alpha(29,1)
			write(*,*) "returns to exp, formal sector (woman with",hed,"yrs of exp):",ret        
			ret=alpha(33,1)+2*hed*alpha(34,1)
			write(*,*) "returns to exp, informal sector (woman with",hed,"yrs of exp):",ret        
		enddo 
		
		! Generate distribution of wage offers
		hed= 10
		wed= 14
		experh=35
		experw=18
		wage_sim=0.0d-0
		earnp(1,1,1) = alpha(16,1)*dexp(alpha(17,1)*(hed-12)+ alpha(18,1)*experh+alpha(19,1)*experh**2.0d-0+alpha(20,1)*hed*hed)*1000000.0d-0 
		earnp(2,1,1) = alpha(21,1)*dexp(alpha(22,1)*(hed-12)+ alpha(23,1)*experh+alpha(24,1)*experh**2.0d-0+alpha(25,1)*hed*hed)*1000000.0d-0 
		earnp(1,2,1) = alpha(26,1)*dexp(alpha(27,1)*(wed-12)+ alpha(28,1)*experw+alpha(29,1)*experw**2.0d-0+alpha(30,1)*wed*wed)*1000000.0d-0
		earnp(2,2,1) = alpha(31,1)*dexp(alpha(32,1)*(wed-12)+ alpha(33,1)*experw+alpha(34,1)*experw**2.0d-0+alpha(35,1)*wed*wed)*1000000.0d-0
		earnp(1,1,2) = alpha(16,2)*dexp(alpha(17,1)*(hed-12)+ alpha(18,1)*experh+alpha(19,1)*experh**2.0d-0+alpha(20,1)*hed*hed)*1000000.0d-0 
		earnp(2,1,2) = alpha(21,2)*dexp(alpha(22,1)*(hed-12)+ alpha(23,1)*experh+alpha(24,1)*experh**2.0d-0+alpha(25,1)*hed*hed)*1000000.0d-0 
		earnp(1,2,2) = alpha(26,2)*dexp(alpha(27,1)*(wed-12)+ alpha(28,1)*experw+alpha(29,1)*experw**2.0d-0+alpha(30,1)*wed*wed)*1000000.0d-0
		earnp(2,2,2) = alpha(31,2)*dexp(alpha(32,1)*(wed-12)+ alpha(33,1)*experw+alpha(34,1)*experw**2.0d-0+alpha(35,1)*wed*wed)*1000000.0d-0
		var=0.0d-0
		var(1,1)=alpha(105,1)
		var(2,1)=alpha(106,1)
		var(1,2)=alpha(107,1)
		var(2,2)=alpha(108,1)
			
		do draw=1,nwdraws
			do sec=1,2
				do gen=1,2
					v(1,1)=var(sec,gen)
					call get_eps(wdraws(sec,gen,:),var(sec,gen),1,nwdraws,idum,1)
				
					do typ=1,2
							earn(sec,gen,typ)=earnp(sec,gen,typ)*dexp(wdraws(sec,gen,draw))
							bin=min0(int(earn(sec,gen,typ)/binsize)+1,nbins)
							wage_sim(sec,gen,typ,bin)=wage_sim(sec,gen,typ,bin)+1
					enddo
				enddo
			enddo
		enddo
		
		write(*,*) "wage offers"
		write(*,*)
		write(*,fmt='(A5,A1,A5,A1,8A9)') "","","","","tp 1","","","","tp 2","","",""
		write(*,fmt='(A5,A1,A5,A1,8A9)') "","","","","men","","wom","","men","","wom",""
		write(*,fmt='(A5,A1,A5,A1,8A9)') "","","","","for","infor","for","infor","for","infor","for","infor"
		write(*,*)
		write(*,*) "- permanent component"
		do i=1,5
		   write(*,fmt='(A5,A1,A5,A1,8f9.4)')"","","","",alpha(15+i,1),alpha(20+i,1),alpha(25+i,1),alpha(30+i,1),alpha(15+i,2),alpha(20+i,2),alpha(25+i,2),alpha(30+i,2)
		enddo
		write(*,*) 
		write(*,fmt='(A5,A1,A5,A1,8f9.0)')"","","","",earnp(1,1,1),earnp(2,1,1),earnp(1,2,1),earnp(2,2,1),earnp(1,1,2),earnp(2,1,2),earnp(1,2,2),earnp(2,2,2)    
		write(*,*)
		write(*,*) "- variance"
		write(*,fmt='(A5,A1,A5,A1,8f9.4)')"","","","",alpha(105,1),alpha(106,1),alpha(107,1),alpha(108,1),alpha(105,1),alpha(106,1),alpha(107,1),alpha(108,1)
		write(*,*)
		wage_dist=100*wage_sim/nwdraws
		do bin=1,nbins
			binubound=bin*binsize/1000
			binlbound=(bin-1)*binsize/1000
			write(*,fmt='(I5,A1,I5,A1,8f9.1)') binlbound,"-",binubound,":", wage_dist(1,1,1,bin),wage_dist(2,1,1,bin),wage_dist(1,2,1,bin),wage_dist(2,2,1,bin), &
					& wage_dist(1,1,2,bin),wage_dist(2,1,2,bin),wage_dist(1,2,2,bin),wage_dist(2,2,2,bin)
		enddo
	endif

	if (prefparams.eq.1) then

		female=1
		married = 1
		ct = 40
		hage = 55
		wage = ct+15
		mptype = 1
		numk = 0

		ind65h = 0
		ind60w = 0

		hw = 1
		ww = 1
		formh = 1
		formw = 1
		fswitch = 0
		freentry = 0 
		mswitch = 0
		mreentry = 0

		nsafebase = 5000000.0d-0
		nsafe=nsafebase
		npensionh = 10000000.0d-0
		npensionw = 10000000.0d-0
		cbase =2500000.0d-0
		c=cbase

		yyevec = 0.0d-0


		call ProbDec(deathprob,cnu)
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

		beta=alpha(15,1)

		write(*,*)
		write(*,*) "Preference parameters"
		write(*,*)
		write(*,5) alpha(1,1),"",alabel(1)
		write(*,5) alpha(2,1),"",alabel(2)
		write(*,5) alpha(3,1),"",alabel(3)
		write(*,5) alpha(4,1),"",alabel(4)
		write(*,5) alpha(5,1),"",alabel(5)
		write(*,5) alpha(6,1),"",alabel(6)
		write(*,5) alpha(6,2),"",alabel(6)
		write(*,5) alpha(7,1),"",alabel(7)
		write(*,5) alpha(7,2),"",alabel(7)
		write(*,5) alpha(8,1),"",alabel(8)
		write(*,5) alpha(9,1),"",alabel(9)
		write(*,5) alpha(10,1),"",alabel(10)
		write(*,5) alpha(11,1),"",alabel(11)
		write(*,5) alpha(13,1),"",alabel(13)
		write(*,5) alpha(14,1),"",alabel(14)
		write(*,5) alpha(15,1),"",alabel(15)
		write(*,5) alpha(62,1),"",alabel(62)
		write(*,5) alpha(63,1),"",alabel(63)
		write(*,5) alpha(64,1),"",alabel(64)
		write(*,5) alpha(65,1),"",alabel(65)
		write(*,5) alpha(70,1),"",alabel(70)
		write(*,5) alpha(70,2),"",alabel(70)
		write(*,5) alpha(71,1),"",alabel(71)
		write(*,5) alpha(71,2),"",alabel(71)
		write(*,5) alpha(81,1),"",alabel(81)
		write(*,5) alpha(83,1),"",alabel(83)
		write(*,5) alpha(110,1),"",alabel(110)
		write(*,5) alpha(111,1),"",alabel(111)
		write(*,*)
		5 format(f28.15,A5,A50)

		do mptype=1,2

			write(*,*)
			write(*,*) "type:",mptype
			write(*,*)

			! Baseline utility: married couple, both working formally, no kids
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm, ubequest)

			util1 = util
			ubequest1 = ubequest

			! Baseline utility: married couple, both working formally, two kids
			numk=2 
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			util1k = util
			ubequest1k=ubequest
			numk=0

			! Baseline utility: single men working formally, no kids

			married=0
			female=0
			c=cbase/2.0d-0    
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			util2 = util
			c=cbase
			married=1
			female=1

			! Baseline utility: single women working formally, no kids

			married=0
			c=cbase/2.0d-0    
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			util3 = util
			c=cbase    
			married=1
			
			! Get marginal utility of consumption
			
			write(*,*) "Marginal utility of consumption"
			write(*,*)
			c=cbase+1.0d-0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MUC1=(util-util1)/(c-cbase)
			write(*,6) "MUC married couple:","",MUC1
			c=cbase

			numk=2
			c=cbase+1.0d-0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MUC1k=(util-util1k)/(c-cbase)
			write(*,6) "MUC married couple with child:","",MUC1k
			numk=0
			c=cbase

			married=0
			c=cbase/2.0d-0+1.0d-0    
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MUC3=(util-util3)/(c-cbase/2.0d-0)
			write(*,6) "MUC single women:","",MUC3
			c=cbase    
			married=1

			! marginal utility of bequest
			
			write(*,*) "Marginal rates of substitution"
			write(*,*)
			nsafe=nsafebase+1.0d-0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MUbq=phdec*pwdec*(ubequest-ubequest1)/(nsafe-nsafebase) 
			MUbq=MUbq/MUC1
			write(*,6) "MU of bequest:","",MUbq
			MUbq=(ubequest-ubequest1)/(nsafe-nsafebase) 
			MUbq=MUbq/MUC1
			write(*,6) "MU of bequest last period:", "", MUbq
			nsafe=nsafebase

			nsafe=nsafebase+1.0d-0
			numk=2
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MUbq=phdec*pwdec*(ubequest-ubequest1k)/(nsafe-nsafebase) 
			MUbq=MUbq/MUC1k
			write(*,6) "MU of bequest, 2 kids:","",MUbq
			MUbq=(ubequest-ubequest1k)/(nsafe-nsafebase) 
			MUbq=MUbq/MUC1k
			write(*,6) "MU of bequest last period, 2 kids:", "", MUbq
			nsafe=nsafebase
			numk=0

			! marginal rate of substitution for leisure

			!   draw leisure shocks
			v(1,1)=alpha(110,1)*alpha(110,1)
			call get_eps(leisdraws(1,:),v,1,1000,idum,1)
			leisuredraws(1,:)=leisdraws(1,:)
			call piksrt(1000,leisuredraws(1,:))
			v(1,1)=alpha(111,1)*alpha(111,1)
			call get_eps(leisdraws(2,:),v,1,1000,idum,1)
			leisuredraws(2,:)=leisdraws(2,:)
			call piksrt(1000,leisuredraws(2,:))
			
			hw=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MULh=util-util1
			MRSh=MULh/MUC1
			write(*,6) "MRS married men, no kids:","",MRSh
			hw=1

			yyevec(6)=leisuredraws(1,250)
			hw=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MULh=util-util1
			MRSh=MULh/MUC1
			write(*,6) "MRS married men, no kids p25:","",MRSh
			hw=1
			yyevec(6)=0.0d-0
		   
			yyevec(6)=leisuredraws(1,750)
			hw=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MULh=util-util1
			MRSh=MULh/MUC1
			write(*,6) "MRS married men, no kids p75:","",MRSh
			hw=1
			

			yyevec(7)=0.0d-0
			ww=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MULw=util-util1
			MRSw=MULw/MUC1
			!write(*,6) "MU of married women's leisure:","",MULw
			write(*,6) "MRS married women, no kids:","",MRSw
			ww=1

			yyevec(7)=leisuredraws(2,250)
			ww=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MULw=util-util1
			MRSw=MULw/MUC1
			!write(*,6) "MU of married women's leisure:","",MULw
			write(*,6) "MRS married women, no kids p25:","",MRSw
			ww=1
			yyevec(7)=0.0d-0
		   
			yyevec(7)=leisuredraws(2,750)
			ww=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MULw=util-util1
			MRSw=MULw/MUC1
			!write(*,6) "MU of married women's leisure:","",MULw
			write(*,6) "MRS married women, no kids p75:","",MRSw
			ww=1
			yyevec(7)=0.0d-0

			ww=0
			married=0
			c=c/2.0d-0    
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MULw3=util-util3
			MRSw3=MULw3/MUC3 
			!write(*,6) "MU of single women's leisure:","",MULw3
			write(*,6) "MRS Single women:","",MRSw3
			married=1
			c=c*2.0d-0    
			ww=1

			numk=2
			ww=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			MRSwk=(util-util1k)/MUC1k
			write(*,6) "MRS of married women with kid's leisure:","" ,MRSwk
			numk=0
			ww=1

			! non-pecuniary benefits

			formh=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			nonpech=(util-util1)/MUC1
			write(*,6)"MU of men's formal work rel to consumption:","", nonpech
			formh=1

			formw=0
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			nonpecw=(util-util1)/MUC1
			write(*,6)"MU of men's formal work rel to consumption:","", nonpecw
			formw=1

			! switching costs

			mswitch=1
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			mswitchcost=(util-util1)/MUC1
			write(*,6)"Men's switch. costs rel to consumption:","" ,mswitchcost
			mswitch=0

			fswitch=1
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			fswitchcost=(util-util1)/MUC1
			write(*,6)"Women's switch. costs rel to consumption:","", fswitchcost
			fswitch=0

			! reentry costs

			mreentry=1
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			mreentrycost=(util-util1)/MUC1
			write(*,6)"Men's rentry costs rel to consumption:","", mreentrycost
			mreentry=0

			freentry=1
			call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
			freentrycost=(util-util1)/MUC1
			write(*,6)"Women's rentry costs rel to consumption:","", freentrycost
			freentry=0
			
		enddo
		
	endif

	6 format(A50,A5,f28.15)
	
	!----------------------------------------------------------------------------------------------------
	!numper=75
	!ctt=50
	!o=1
	!phdec
	!pwdec
	!pdiv
	!probfert
	!numz=106
	!intmat
	!cnu
	!maxearnh
	!maxearnf
	!
	!tao=0.10d-0
	!pasis = 565236.0d-0 ! this correspond to the value of PASIS for 70-75 years old in august 2007
	!mpg = 2*pasis
	!pmin = 12.0d-0*500000.0d-0
	!pmax = pmin
	!bono=tao*1.5d-0*165000.0d-0*12.0d-0
	!cmin = alpha(85,1)*1000000.0d-0    ! set a minimum consumption level
	!
	!
	!reform=0
	!nodivrule=0
	!
	!married
	!female
	!nsafe
	!npensionh
	!npensionw
	!experh
	!experw
	!numk
	!fexperh
	!fexperw
	!
	!nivart=
	!nfixvar=
	!nrvart=
	!
	!call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
	!call contvalue(married,female,numper,term4,ubequest,theta,reform,nodivrule,nivart,nrvart,nfixvar,ct,o,pdiv,probfert,nsafe,npensionh,npensionw,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,numz,totp,intmat,pmax,pmin,pasis,bono,cmin,cnu,lastage,lastworkage)

endsubroutine paramsinterp
