PROGRAM savingsmod 

	implicit none

	!external DrawHist,ReadIn,ran1,EVMPI,EVFUN,createinputfile

	integer retVal

	CALL MainProg(retVal,1)

27	END PROGRAM savingsmod

SUBROUTINE MainProg(Out_Arg, In_Arg)
	implicit none
	
	integer Out_Arg, In_Arg  
 
	integer, parameter :: ncoeff=150
	integer, parameter :: maxnumtypes=4
	integer ktp 
	integer:: silent,reform
	character*200::pathoutput
	real(8) thetamtype(2,75,3,106),thetamreftype(2,75,3,106),thetamreftype1(75,3,106),thetamtype1(75,3,106)
	integer totp,kidmax,numper
	real(8) vseed(3)
	real(8),allocatable :: parms(:)
	integer numiter,numz,maxit
	real(8) lval
	real(8) alphaa(ncoeff,maxnumtypes),abumpp(ncoeff,maxnumtypes),rtol,ftol,tausm
	integer aiterr(ncoeff,maxnumtypes)
	real(8) alphaloww(ncoeff,maxnumtypes),alphahighh(ncoeff,maxnumtypes)
	real(8), allocatable :: alpha(:,:),abump(:,:)
	integer, allocatable ::  aiter(:,:)
	real(8), allocatable ::  alphalow(:,:),alphahigh(:,:)
	integer Sbar,Tbar,T0,nmc,nsimul,simul,dnumdraw,snumdraw
	real(8),allocatable :: rvart(:,:,:)
	integer,allocatable :: ivart(:,:,:)   
	integer, allocatable :: fixvar(:,:)        
	real(8), allocatable :: maxsafe(:)    
	integer regind(75),intmat(75,75),numint,i,nedraw,ncounter
	CHARACTER inarg*100,outarg*100,tag1*100,temp3
	REAL(8) argc,cbar
	integer temp2,jj,kk,counter,mtype,estimation,wageparams,segmentparams,typeparams,prefparams
	character*50::alabel(300)	
	integer :: ibump,  bumpup
	real :: bumpsize,range,bump

	!! Clock variables
	integer::						time_array_0(8), time_array_1(8)
	real::							start_time, finish_time

	!!	Start clock
	call date_and_time(values=time_array_0)
	start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 + time_array_0 (7) + 0.001 * time_array_0 (8)
	
	estimation = 1
	
	wageparams = 0
	typeparams = 0
	prefparams = 0
	segmentparams = 0
	
	ncounter=1
	if (estimation.eq.1) ncounter=1
	
	COUNTERLOOP: do counter=1,ncounter
		
		!Get file containing parameter estimates (infile.asc)
		argc=IARGC()
		CALL getarg(1,inarg)
		CALL getarg(2,outarg)
		CALL getarg(3,tag1)
		if (estimation.eq.0) then
			inarg='infile.asc'
			outarg='outfile.asc'
		endif
		
		!Initialize variables
		totp = 0
		kidmax = 0
		numper = 0
		vseed = 0
		numiter = 0
		alphaa = 0.0d-0
		alphaloww = 0.0d-0
		alphahighh = 0.0d-0
		abumpp = 0.0d-0
		aiterr = 0
		Tbar = 0
		ktp = 0
		T0 = 36 
		Sbar = 0
			
		!Read parameters in (estimated and not estimated) from input.txt
		call Readin(nedraw,totp,kidmax,ktp,numper,vseed,alphaa,abumpp,aiterr,numiter,Sbar,Tbar,nmc,&
				regind,intmat,numint,nsimul,simul,dnumdraw,snumdraw,tausm,rtol,ftol,maxit,alphaloww,alphahighh)

		!Store in allocatable arrays	
		allocate (alpha(ncoeff,ktp))
		allocate (alphalow(ncoeff,ktp))
		allocate (alphahigh(ncoeff,ktp))
		allocate (aiter(ncoeff,ktp))
		allocate (abump(ncoeff,ktp))
		alpha = 0.0d-0
		alphalow = 0.0d-0
		alphahigh = 0.0d-0
		abump = 0.0d-0
		aiter = 0
		do i=1,ktp
			alpha(:,i)=alphaa(:,i)
			alphalow(:,i)=alphaloww(:,i)
			alphahigh(:,i)=alphahighh(:,i)
			abump(:,i)=abumpp(:,i)
			aiter(:,i)=aiterr(:,i)
		enddo

		numz = 75+numint
		numiter = sum(aiter(1:ncoeff,1:ktp))

		allocate(parms(numiter))
		allocate(rvart(6,totp,numper))      
		allocate(ivart(9,totp,numper))   
		allocate(fixvar(4,totp))
		allocate(maxsafe(numper))

		!Read latest parameters
		open(unit=37,file=inarg,position='rewind')
		read(37,*) temp3
		read(37,*) temp2

		!read in the initial parameter values 
		do 100 i = 1,numiter
			read(37,*) parms(i)
		parms(i) = parms(i)*1.0d-0
		100 continue 
		close(37)

		!Draw histories for states at which to evaluate Emax
		call DrawHist(totp,kidmax,ktp,numper,vseed,fixvar,ivart,rvart,maxsafe)
		lval = 0.0d-0
		
		!Replace estimable parameters by the last estimates and bump up/down if computing std errors
		ibump = 0 
		bump = 0.0d-0
		bumpsize = 0.0d-0
		bumpup = 1	
		bumpup=bumpup - 1 ! recoding 0 into -1, 1 to 0 and 2 to 1	
		i = 0
		if (sum(aiter).gt.0) then
			jj1: do jj = 1,ncoeff
				kk1: do kk = 1,ktp
					if (aiter(jj,kk).eq.1) then
						i = i + 1
						range = alphahigh(jj,kk)-alphalow(jj,kk)
						bump = range*bumpup*bumpsize
					if (i.ne.ibump) then
								bump = 0 
					else if (i.eq.ibump) then
						write(*,*) "bumping alpha ",jj,",",kk," by ",bump
					endif
					alpha(jj,kk)=parms(i)+bump
					end if
				end do	kk1 !jj
			end do	jj1 !kk
		end if !iterating on values of alpha

		!normalize some preference parameter values by the marginal utility of consumption evaluated at cbar
		cbar=3.0d-0
		do mtype=1,ktp
			alpha(6,mtype)=alpha(6,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alpha(7,mtype)=alpha(7,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alpha(70,mtype)=alpha(70,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alpha(71,mtype)=alpha(71,mtype)*(cbar**(alpha(1,1)-1.0d-0))
		enddo
		alpha(62,1)=alpha(62,1)*(cbar**(alpha(1,1)-1.0d-0))
		alpha(63,1)=alpha(63,1)*(cbar**(alpha(1,1)-1.0d-0))
		alpha(64,1)=alpha(64,1)*(cbar**(alpha(1,1)-1.0d-0))
		alpha(65,1)=alpha(65,1)*(cbar**(alpha(1,1)-1.0d-0))
		
		do mtype=1,ktp
			alphahigh(6,mtype)=alphahigh(6,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alphahigh(7,mtype)=alphahigh(7,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alphahigh(70,mtype)=alphahigh(70,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alphahigh(71,mtype)=alphahigh(71,mtype)*(cbar**(alpha(1,1)-1.0d-0))
		enddo
		alphahigh(62,1)=alphahigh(62,1)*(cbar**(alpha(1,1)-1.0d-0))
		alphahigh(63,1)=alphahigh(63,1)*(cbar**(alpha(1,1)-1.0d-0))
		alphahigh(64,1)=alphahigh(64,1)*(cbar**(alpha(1,1)-1.0d-0))
		alphahigh(65,1)=alphahigh(65,1)*(cbar**(alpha(1,1)-1.0d-0))
		do mtype=1,ktp
			alphalow(6,mtype)=alpha(6,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alphalow(7,mtype)=alpha(7,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alphalow(70,mtype)=alpha(70,mtype)*(cbar**(alpha(1,1)-1.0d-0))
			alphalow(71,mtype)=alpha(71,mtype)*(cbar**(alpha(1,1)-1.0d-0))
		enddo
		alphalow(62,1)=alphalow(62,1)*(cbar**(alpha(1,1)-1.0d-0))
		alphalow(63,1)=alphalow(63,1)*(cbar**(alpha(1,1)-1.0d-0))
		alphalow(64,1)=alphalow(64,1)*(cbar**(alpha(1,1)-1.0d-0))
		alphalow(65,1)=alphalow(65,1)*(cbar**(alpha(1,1)-1.0d-0))		
		
		!parameter interpretation
		if (wageparams.eq.1.or.typeparams.eq.1.or.prefparams.eq.1.or.segmentparams.eq.1) then
			reform=0
			pathoutput="./outputs/"
			call readtheta2(1,pathoutput,reform,thetamtype,thetamreftype,75,106,ktp)
			call createinputfile(parms,numiter,pathoutput,alabel,ktp)
			call paramsinterp(wageparams,typeparams,prefparams,segmentparams,alpha,ncoeff,ktp,Tbar,thetamtype1,totp,alabel)
			go to 34
		endif	  
		
		!Solve for Emax and simulate moments
		call EVMPI(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed, &
			fixvar,ivart,rvart,lval,alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc, &
			regind,intmat,numint,nsimul,maxsafe,pathoutput,counter,estimation)
		
		!Output MSM criterion
		open(unit=11,file=outarg)
		write(11,*) 1 
		write(11,*) lval 
		if (estimation.eq.0) write(*,*) 'MSM criterion==',lval
		close(unit=11)

		!Create new input file recording parameter values used in this run	
		pathoutput="./outputs/"
		if (estimation.ne.1) call createinputfile(parms,numiter,pathoutput,alabel,ktp)
			
		deallocate(parms)
		deallocate(rvart)      
		deallocate(ivart)   
		deallocate(fixvar)
		deallocate(maxsafe)
		deallocate(alpha)
		deallocate(alphalow)
		deallocate(alphahigh)
		deallocate(abump)
		deallocate(aiter)

		if (estimation.eq.0) print*, "counterloop", counter
	
		call date_and_time(values=time_array_1)
		finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60+ time_array_1 (7) + 0.001 * time_array_1 (8)
		
		if (estimation.eq.0) then
			write(*,*) "starting time", time_array_0(5),":", time_array_0(6),":",time_array_0(7) 
			write(*,*) "finish time"  , time_array_1(5),":", time_array_1(6),":",time_array_1(7) 
			write (*,*) 'elapsed wall clock time:',(finish_time - start_time)/60   
		endif
		
	enddo COUNTERLOOP
 	
	34  Out_Arg = 222
	return
 
end subroutine MainProg

!   27     END PROGRAM savingsmod 

real(8) function logic2dbl(a)

	!   This function added for gfortran compatibility changes
	!   elemental pure double precision function logic2dbl(a)
	logical, intent(in) :: a

	if (a) then
		logic2dbl = -1.d0
	else
		logic2dbl = 0.d0
	end if
	
end function logic2dbl
	
