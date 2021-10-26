  
  PROGRAM savingsmod 

	implicit none
	
 !       external DrawHist,ReadIn,ran1,EVMPI,EVFUN,dp,createinputfile
 	
	integer retVal
	
	CALL MainProg(retVal,1)

27	END PROGRAM savingsmod

SUBROUTINE MainProg(Out_Arg, In_Arg)
	implicit none
	
	integer Out_Arg, In_Arg     
 
        integer, parameter :: ncoeff=150
        integer, parameter :: maxnumtypes=4
        integer ktp 
        integer:: silent,reform,group,xx,zz
        character*200::pathoutput
        real(8) thetamtype(2,75,3,106),thetamreftype(2,75,3,106),thetamtype1(75,3,106),thetamtype2(75,3,106),thetamreftype1(75,3,106),thetamreftype2(75,3,106)

        integer totp,kidmax,numper
        real(8) vseed(3)
        real(8),allocatable :: workh(:,:)
        real(8),allocatable :: workw(:,:)
        integer,allocatable :: hed(:)
        integer,allocatable :: wed(:)
        integer,allocatable :: n(:,:)
        integer,allocatable :: hage(:)
        integer,allocatable :: typ(:)
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
        real(8) dpmat(11,3,2)
        integer regind(75),intmat(75,75),numint,j,i,nedraw,ncounter
        CHARACTER inarg*100,outarg*100,tag1*100,tt,temp3
        REAL(8) temp1,argc,cbar
        integer temp2,jj,kk,counter,mtype,estimation,interpparams,wageparams,segmentparams,typeparams,prefparams
        character*50::alabel(300)
	
	integer :: ibump,  bumpup
	real :: bumpsize,range,bump

        !! Clock variables
        integer::						time_array_0(8), time_array_1(8)
        real::							start_time, finish_time

        !!	Start clock
        call date_and_time(values=time_array_0)
        start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 + time_array_0 (7) + 0.001 * time_array_0 (8)
        
        estimation = 0
        
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
        	call DP(dpmat) !!XXX

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
		bumpsize = 0.001d-0
		bumpup = 0
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
        
        
!        	!parameter interpretation
        	if (wageparams.eq.1.or.typeparams.eq.1.or.prefparams.eq.1.or.segmentparams.eq.1) then
            		reform=0
            		pathoutput="./outputs/"
            		call readtheta2(1,pathoutput,reform,thetamtype,thetamtype1,thetamtype2,thetamreftype,thetamreftype1,thetamreftype2,75,106,ktp)
            		call createinputfile(parms,numiter,pathoutput,alabel,ktp)
            		call paramsinterp(wageparams,typeparams,prefparams,segmentparams,alpha,ncoeff,ktp,Tbar,thetamtype1,totp,alabel)
            		go to 34
        	endif
      
        
        	!Solve for Emax and simulate moments
        	call EVMPI(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed, &
              	fixvar,ivart,rvart,lval,alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc,dpmat, &
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
 	
34    Out_Arg = 222
return
 
end subroutine MainProg

!   27     END PROGRAM savingsmod 

real(8) function logic2dbl(a)

 !   This function added for gfortran compatibility changes
 !    elemental pure double precision function logic2dbl(a)
         logical, intent(in) :: a
 
         if (a) then
         	logic2dbl = -1.d0
         else
         	logic2dbl = 0.d0
         end if
     end function logic2dbl
	

         SUBROUTINE Readin(nedraw,totp,kidmax,ktp,numper,vseed,alphaa,abumpp,aiterr,numiter, &
              Sbar,Tbar,nmc,regind,intmat,numint,nsimul,simul,dnumdraw,snumdraw,tausm,rtol,ftol,maxit,alphaloww,alphahighh)


!	This subroutine reads in the parameter values
!     from the file input.txt.

        implicit none

	integer, parameter :: ncoeff=150
	integer, parameter :: maxnumtypes=4


	integer ktp
        integer totp,Tbar,Sbar,kidmax
	integer nedraw,nmc,FTbar,vl
	integer simul,nsimul
	integer dnumdraw,snumdraw
        real(8) vseed(3)
	real(8) disc
	real(8) ftol,rtol,temp1
	integer maxit,temp2
        real(8) alphaa(ncoeff,maxnumtypes),abumpp(ncoeff,maxnumtypes)
        integer aiterr(ncoeff,maxnumtypes)
        real(8) alphaloww(ncoeff,maxnumtypes),alphahighh(ncoeff,maxnumtypes)
        integer i,j,numper,numiter
        real(8) tausm,tempvar
        integer checkp,numint,f1ind,sind,regind(75),intmat(75,75)

        alphaa = 0.0
        abumpp = 0.0
        alphahighh = 0.0d-0
        alphaloww = 0.0d-0
        aiterr = 0
        Sbar = 0
        Tbar = 0
        kidmax = 0
        nedraw = 0
        nmc = 0
        i = 0
        j = 0 



        open(unit=100, file='./inputs/input.txt', position='rewind')

	
	read(100,*) totp
	read(100,*) Tbar
	read(100,*) Sbar


        numper = Tbar-Sbar+1
	FTbar = Tbar-Sbar-2
	vl = FTbar*totp
	read(100,*) kidmax
	read(100,*) nedraw
	read(100,*) nmc
	read(100,*) vseed(1)
	read(100,*) vseed(2)
	read(100,*) vseed(3)
	read(100,*) ktp
	read(100,*) disc
        read(100,*) temp2
	read(100,*) simul
	read(100,*) nsimul
	read(100,*) dnumdraw
	read(100,*) snumdraw
	read(100,*) tausm
	read(100,*) temp1  
	read(100,*) ftol
	read(100,*) rtol
	read(100,*) maxit
      

        
	
	Loop4: do i=1,132
	     if (i.eq.1.or.i.eq.6.or.i.eq.7.or.i.eq.16.or.i.eq.21.or.i.eq.26.or.i.eq.31.or.i.eq.70.or.i.eq.71.or.i.eq.75.or.i.eq.76.or.i.eq.77.or.i.eq.78.or.i.eq.79) then  
 					
               loopj: do j = 1,maxnumtypes		  
                 read(unit=100,fmt=*) alphaa(i,j),alphaloww(i,j),alphahighh(i,j),abumpp(i,j),aiterr(i,j)
	       end do loopj
             else
               read(unit=100,fmt=*) alphaa(i,1),alphaloww(i,1),alphahighh(i,1),abumpp(i,1),aiterr(i,1)
             end if
	end do Loop4

	read(100,*) tempvar
	read(100,*) tempvar
	read(100,*) tempvar
	read(100,*) tempvar
	read(100,*) tempvar
		  
	!read in the indicators for whether to include each regressor

	regrd: do i = 1,75
	   read(100,*) regind(i)
	end do regrd

	
	read(100,*) tempvar
	read(100,*) tempvar
	read(100,*) tempvar
	read(100,*) tempvar
	read(100,*) tempvar


	! read in the location of the interactions

	intmat = 0
	checkp = 0
	numint = 0
	do while(checkp.eq.0)
	read(100,*) f1ind, sind
	if (f1ind.ne.-9) then
	intmat(f1ind,sind) = 1
	intmat(sind,f1ind) = 1
 	numint = numint+1
	else
	checkp = 1
	end if
	end do



	close(unit=100)


	return
	end subroutine ReadIn


SUBROUTINE dp(X)
!---------------------------------------------------------------
! This subroutine reads in the info in DeathProbabilities.txt
! and with with fills an array of dimension (11,3,2); where
! the first dimension is for ages, the second for race (W,B,H);
! and the third for gender (male, female);
!--------------------------------------------------
IMPLICIT NONE
REAL(8), DIMENSION(11,3,2), INTENT(OUT) :: X
INTEGER :: i,j
OPEN(UNIT=1,FILE='inputs/dp.txt',POSITION='REWIND')
READ(1,FMT=*) 
READ(1,FMT=*) 
DO i = 1,11
	READ(1,FMT=*) (X(i,j,1), j=1,3), (X(i,j,2), j=1,3)
END DO
CLOSE(UNIT=1) 
RETURN
END SUBROUTINE dp 




	SUBROUTINE DrawHist(totp,kidmax,ktp,numper,vseed,fixvar,ivart,rvart,maxsafe)
	
!       This subroutine draws simulated state points at every marriage duration (=period of choice)
!
!       state variables:
!         total wealth  (wealth)  [-100,00 to +500,000]
!         workexp of male, female - fexp,mexp  
!         education levels of male and female (years of education)
!         number of children (up to kidmax) - n
!         current ages of spuses - hage, wage
!         rates of return: assume fixed so for now we don't draw these

        IMPLICIT NONE
  !      external ran1

	    integer totp,kidmax,ktp,numper
        integer i,o,idum,j
        real(8) cut1,cut2,cut3,cut4,r1(1),r2(1),r3(1),r4(1),r5(1),r6(1),r7(1),r20(1),r21(1)
        real(8) r100(1),r101(1),r102(1),r103(1),r104(1),r8(1),r9(1),r10(1),r11(1),r12(1),r13(1),r14(1),r15(1),r16(1), r17(1),r18(1),r19(1)
        real(8) maxsafe(numper)
        integer n(totp,numper), fert(totp,numper)!n is number of kids, fert is current fetility decision
        integer hages(totp),married(totp,numper)
        real(8) safe(totp,numper),vseed(3),risky(totp,numper),pensionw(totp,numper),pensionh(totp,numper)
        integer initmarried
        integer sum1,firstchild, hage_cur(totp,numper)
        integer hed(totp),wed(totp),mardur(totp,numper),mindur,maxdur
        integer workh(totp,numper),workw(totp,numper)
        real(8) rvart(6,totp,numper)   ! row 1 = safe
                                       ! row 2 = risky 
                                       ! row 3 = rate of return safe
                                       ! row 4 = rate of return risky
                                       ! row 5 - womens pension wealth
                                       ! row 6 - mens pension wealth
        integer ivart(9,totp,numper)   ! row 1 = children (n)
                                       ! row 2 = married
                                       ! row 3 = work h
                                       ! row 4 = work w
                                       ! row 5 = whether male surviving
                                       ! row 6 = whether female surviving
                                       ! row 7 = marriage duration 
                                       ! row 8 = formal h
                                       ! row 9 = formal w
        integer fixvar(4,totp)         ! row 1 = hed
                                       ! row 2 = wed
                                       ! row 3 = hage
                                       ! row 4 = female (needed for singles for divorced - whether male or female) 
        integer survivh(totp,numper),survivw(totp,numper)
        integer formh(totp,numper), formw(totp,numper),female(totp)

        mardur = 0 
        idum = dint(-1.0d-0*abs(vseed(2)))
	    cut1 = 0.2d-0
 	    cut2 = 2.00d-0
	    cut3 = 0.40d-0
	    cut4 = 0.85d-0
        n = 0
        survivh = 1
        survivw = 1
 	    r1 = 0.0d-0
	    r2 = 0.0d-0
	    r3 = 0.0d-0
	    r4 = 0.0d-0
	    r5 = 0.0d-0
	    r6 = 0.0d-0
	    r7 = 0.0d-0
	    r9 = 0.0d-0       
	    fert = 0
        
	! loop through simulated histories loop (will draw totp histories)
    
    OLoop: do o = 1,totp
 

        firstchild = -99
        sum1=0
             

    ! Draw gender, which is needed for singles and to know the gender for the person being followed in the case of divorce

	    female(o) = 1
	    if (o.le.int(totp/4).or.o.gt.int(3*totp/4)) then 
		    female(o) = 0
	    endif

	
    ! Draw marital status    

        married(o,:) = 1
        !call ran1((idum+6),r3,1)
        if (o.gt.int(totp/2)) then 
            married(o,:) = 0
        end if
        !draw husband age
        if (female(o).eq.1.or.married(o,1).eq.1) then
            call ran1(idum,r5,1)
            hages(o) = dint(r5(1)*10.0d-0+16.0d-0)
        else
            hages(o) = 16   !for men who never marry, make ct index go along with their agem
        end if
        
    ! Draw education levels for male and female, ranging from 4 to 20 years
    ! men:

        r7 = 0.0d-0
        call ran1((idum+22),r7,1)
	    hed(o)= dint(r7(1)*16+4)

    ! women: 
        r8 = 0.0d-0
        call ran1((idum+24),r8,1)
	    wed(o)= dint(r8(1)*16+4)
       
    ! Draw fraction of time individual will spend working and workink formally/Informally

	    call ran1((idum+26),r100,1)
	    call ran1((idum+27),r101,1)
	    call ran1((idum+28),r102,1)
	    call ran1((idum+29),r103,1)

    ! Draw expected number of children between 0 and 4 and compute the per period probability that will yield that number
    
        call ran1((idum+30),r104,1)
        cut2=int(5.0d-0*r104(1))/25.0d-0
       
       

    ! Loop over ages 15 to 90
    !   latest date of preg is 40  (even though model will begin at age 35)
    !   Draw number of kids in every time period, taking into account that in period 1 no kids, in period 2 at most 1... etc.

        Loop1: do i = 1,numper   ! ages 15 through 90 
                     
            hage_cur(o,i)=min0(hages(o)+i-1,90)    ! get mans current age, cap at 90
            
     ! randomly draw fertility histories 
     
            if (i.le.25) then   !if woman less than 40
                call ran1((idum+2),r1,1)       
                if (i.eq.1) then
                    sum1 = 0
                end if
                if ((r1(1).lt.cut2).and.sum1.lt.kidmax) then     
	                fert(o,i) = 1
                    ivart(1,o,i) = fert(o,i)
	            end if
                sum1=sum(fert(o,1:(i)))
                n(o,i)=sum1
            else !woman older than 40
                fert(o,i)=0
                ivart(1,o,i)=fert(o,i)
                sum1=sum(fert(o,1:(i-1)))
                n(o,i)=sum1
            end if 

            firstchild = 0 
            if (sum1.eq.0.and.fert(o,i).eq.1) then
                firstchild = i+15   !I set firstchild as age at which woman was pregnant for the first time
            end if

                
    ! for some singles draw whether separated/divorced from previous marriage

            call ran1((idum+6),r3,1)

                
            if (married(o,i).eq.0.and.r3(1).lt.0.003d-0) then   !with some prob get divorced/widowed
                mardur(o,i)=int(i/2) !give mariage duration equal to half the number of elapsed periods
            end if
            
            if (married(o,i).eq.1.and.r3(1).lt.0.8) then
                mardur(o,i)=mardur(o,i)+1   !update marriage duration (do not always update to create variation in mariage duration)
            end if


    ! draw wealth levels between 1 and 500,000 - using safe assets

            call ran1((idum+10),r9,1)
            call ran1((idum),r10,1)
            safe(o,i) = (r9(1)**4.0d-0)**2.0d-0*dble(i+15)*1000000.0d-0
            risky(o,i)=0.0d-0
            maxsafe(i) = maxval(safe(:,i))
            if (r10(1).lt.0.2d-0) then
                safe(o,i) = 0.0d-0
            end if 

        ! draw work pattern - men can only work fulltime or not at all
        ! women can work also part time 
        ! men:
            
            call ran1((idum+14),r10,1)   
            call ran1((idum+16),r11,1)
            workh(o,i)=0
            workw(o,i)=0
            
            if (r10(1).lt.r100(1)) workh(o,i) = 1                    
            if (r11(1).lt.r101(1)) workw(o,i) = 1      
            
        ! draw whether working part-time
!      	        
!            call ran1((idum+28),r11,1)    
!            if (workw(o,i).eq.1.and.r11(1).lt.0.2d-0) then
!                workw(o,i) = 2
!  	        end if

        ! draw whether working in the formal or iformal sector

            call ran1((idum+18),r16,1)
            call ran1((idum+20),r17,1)
            formw(o,i) = 0
            formh(o,i) = 0
            
            if (r16(1).lt.r102(1).and.workh(o,i).eq.1) formh(o,i)=1                    
            if (r17(1).lt.r103(1).and.workw(o,i).eq.1) formw(o,i)=1 
                               
        ! draw women and mens pension levels between 1 and 500,000
         
            call ran1((idum+12),r18,1)
            call ran1((idum+22),r19,1)
            call ran1((idum+24),r20,1)
            call ran1((idum+26),r21,1)

            if (o.ge.(totp-2)) then
                r18(1) = 1.0d-0
                r19(1) = 1.0d-0
            end if

            pensionw(o,i) = ((r18(1))**4)*(3000000.0d-0)*sum(formw(o,1:i))
            pensionh(o,i) = ((r18(1))**4)*(3000000.0d-0)*sum(formh(o,1:i))

        ! Handle always single case
        ! if single, female at ct 35, then set husband pension equal to 0
            
            if (female(o).eq.1.and.married(o,i).eq.0) then
                pensionh(o,i)=0.0d-0
            end if   
                  
        ! if single, male at ct 35, then set wife pension equal to 0
        
            if (female(o).eq.0.and.married(o,i).eq.0) then
                pensionw(o,i)=0.0d-0
            end if
            
	    end do Loop1

         
	
	end do OLoop 


    ivart(1,:,:) = fert(:,:)       !this is not number of kids, it is fertility each period. Will use this to infer 
                                   !   numkids in subroutine u
    ivart(2,:,:) = married(:,:)    !whether married or not
    ivart(3,:,:) = workh(:,:)      !current work choice husband, will use this to compute experience
    ivart(4,:,:) = workw(:,:)      !current work chocie wife, will use this to compute ecperience
    ivart(5,:,:) = survivh(:,:)
    ivart(6,:,:) = survivw(:,:)
    ivart(7,:,:) = mardur(:,:)
    ivart(8,:,:) = formh(:,:)
    ivart(9,:,:) = formw(:,:)
    fixvar(1,:) = hed(:)
    fixvar(2,:) = wed(:)
    fixvar(3,:) = hages(:)
    fixvar(4,:) = female(:)                     
    rvart(1,:,:) = safe(:,:)
    rvart(2,:,:) = risky(:,:) 
    rvart(3,:,:) = 1.05d-0 
    rvart(4,:,:) = 0.0d-0              !not using risky, so set return equal to 0
    rvart(5,:,:) = pensionw(:,:)
    rvart(6,:,:) = pensionh(:,:)

      
        end subroutine DrawHist

subroutine paramsinterp(wageparams,typeparams,prefparams,segmentparams,alpha,ncoeff,ktp,Tbar,theta,totp,alabel)

integer,parameter:: nbins=10
integer,parameter:: binsize=500000
integer,parameter:: nwdraws=100
integer,parameter:: ntypesim=4
integer,parameter:: lastage = 90
integer,parameter:: lastworkage = 70
integer::bidule, hedvec(ntypesim),wedvec(ntypesim),cohortvec(ntypesim),smarriedvec(ntypesim),sfemalevec(ntypesim),hagevec(ntypesim),wagevec(ntypesim)
real(8)::temp111,typeprob(ntypesim),ret
integer:: ncoeff,ktp,totp,ind65h,ind60w,wage,hage,Tbar,ct,mptype,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,parttime,fulltime
real(8):: nsafe,ubequest,ubequest1,ubequest2,ubequest1k,alpha(ncoeff,ktp),npensionh,npensionw,beta,phdec,pwdec,yyevec(7),util,c,utilf,utilm,cf,cm
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
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    util1 = util
    ubequest1 = ubequest

    ! Baseline utility: married couple, both working formally, two kids
    numk=2 
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    util1k = util
    ubequest1k=ubequest
    numk=0

    ! Baseline utility: single men working formally, no kids

    married=0
    female=0
    c=cbase/2.0d-0    
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    util2 = util
    c=cbase
    married=1
    female=1

    ! Baseline utility: single women working formally, no kids

    married=0
    c=cbase/2.0d-0    
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    util3 = util
    c=cbase    
    married=1

    
    ! Get marginal utility of consumption
    
    write(*,*) "Marginal utility of consumption"
    write(*,*)
    c=cbase+1.0d-0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MUC1=(util-util1)/(c-cbase)
    write(*,6) "MUC married couple:","",MUC1
    c=cbase

    numk=2
    c=cbase+1.0d-0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MUC1k=(util-util1k)/(c-cbase)
    write(*,6) "MUC married couple with child:","",MUC1k
    numk=0
    c=cbase

    married=0
    c=cbase/2.0d-0+1.0d-0    
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MUC3=(util-util3)/(c-cbase/2.0d-0)
    write(*,6) "MUC single women:","",MUC3
    c=cbase    
    married=1
    


    ! marginal utility of bequest
    
    write(*,*) "Marginal rates of substitution"
    write(*,*)
    nsafe=nsafebase+1.0d-0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MUbq=phdec*pwdec*(ubequest-ubequest1)/(nsafe-nsafebase) 
    MUbq=MUbq/MUC1
    write(*,6) "MU of bequest:","",MUbq
    MUbq=(ubequest-ubequest1)/(nsafe-nsafebase) 
    MUbq=MUbq/MUC1
    write(*,6) "MU of bequest last period:", "", MUbq
    nsafe=nsafebase

    nsafe=nsafebase+1.0d-0
    numk=2
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
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
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MULh=util-util1
    MRSh=MULh/MUC1
    write(*,6) "MRS married men, no kids:","",MRSh
    hw=1

    yyevec(6)=leisuredraws(1,250)
    hw=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MULh=util-util1
    MRSh=MULh/MUC1
    write(*,6) "MRS married men, no kids p25:","",MRSh
    hw=1
    yyevec(6)=0.0d-0
   
    yyevec(6)=leisuredraws(1,750)
    hw=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MULh=util-util1
    MRSh=MULh/MUC1
    write(*,6) "MRS married men, no kids p75:","",MRSh
    hw=1
    

    yyevec(7)=0.0d-0
    ww=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MULw=util-util1
    MRSw=MULw/MUC1
    !write(*,6) "MU of married women's leisure:","",MULw
    write(*,6) "MRS married women, no kids:","",MRSw
    ww=1

    yyevec(7)=leisuredraws(2,250)
    ww=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MULw=util-util1
    MRSw=MULw/MUC1
    !write(*,6) "MU of married women's leisure:","",MULw
    write(*,6) "MRS married women, no kids p25:","",MRSw
    ww=1
    yyevec(7)=0.0d-0
   
    yyevec(7)=leisuredraws(2,750)
    ww=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MULw=util-util1
    MRSw=MULw/MUC1
    !write(*,6) "MU of married women's leisure:","",MULw
    write(*,6) "MRS married women, no kids p75:","",MRSw
    ww=1
    yyevec(7)=0.0d-0

    ww=0
    married=0
    c=c/2.0d-0    
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MULw3=util-util3
    MRSw3=MULw3/MUC3 
    !write(*,6) "MU of single women's leisure:","",MULw3
    write(*,6) "MRS Single women:","",MRSw3
    married=1
    c=c*2.0d-0    
    ww=1

    numk=2
    ww=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    MRSwk=(util-util1k)/MUC1k
    write(*,6) "MRS of married women with kid's leisure:","" ,MRSwk
    numk=0
    ww=1

    ! non-pecuniary benefits

    formh=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    nonpech=(util-util1)/MUC1
    write(*,6)"MU of men's formal work rel to consumption:","", nonpech
    formh=1

    formw=0
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    nonpecw=(util-util1)/MUC1
    write(*,6)"MU of men's formal work rel to consumption:","", nonpecw
    formw=1

    ! switching costs

    mswitch=1
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    mswitchcost=(util-util1)/MUC1
    write(*,6)"Men's switch. costs rel to consumption:","" ,mswitchcost
    mswitch=0

    fswitch=1
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    fswitchcost=(util-util1)/MUC1
    write(*,6)"Women's switch. costs rel to consumption:","", fswitchcost
    fswitch=0

    ! reentry costs

    mreentry=1
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
    mreentrycost=(util-util1)/MUC1
    write(*,6)"Men's rentry costs rel to consumption:","", mreentrycost
    mreentry=0

    freentry=1
    call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
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
!
!nivart=
!nfixvar=
!nrvart=
!
!call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)
!call contvalue(married,female,numper,term4,ubequest,theta,reform,nodivrule,nivart,nrvart,nfixvar,ct,o,phdec,pwdec,pdiv,probfert,nsafe,npensionh,npensionw,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,numz,totp,intmat,pmax,pmin,pasis,bono,cmin,cnu,lastage,lastworkage)

endsubroutine paramsinterp

!This function will be called by the simplex program
!Given a set of parameters, it will calculate
!the emax coefficients and return the likelihood value


        SUBROUTINE EVMPI(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed,fixvar,ivart,rvart, &
          lval,alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc,dpmat,regind,intmat,numint,nsimul,maxsafe,pathoutput,counter,estimation)

	! specify the number of columns of the theta matrix
	! and the number of time periods in which decisions are
	! being made (from loop runs from age 59 to age 13)

        implicit none

!        external EVFUN


        integer, parameter :: ncoeff=150
        integer, parameter :: nmoments=258        !number of moments used in the estimation

        integer mptype,numiter,i,nsimul,nedraw,estimation
        real(8) loss(1,1:nmoments)
        integer totp,kidmax,ktp,numper,numint,regind(75),intmat(75,75)
        real(8) vseed(3)
        real(8) lval 
        real(8) parms(numiter)
        real(8) alpha(ncoeff,ktp),alphahigh(ncoeff,ktp),alphalow(ncoeff,ktp),abump(ncoeff,ktp)
        integer aiter(ncoeff)
        real(8) flik(totp),tprob(totp)
        integer Sbar,Tbar,T0,numz,nmc,counter
        real(8) rvart(6,totp,numper)   ! row 1 = safe assets
                                       ! row 2 = risky 
                                       ! row 3 = rate of return safe
                                       ! row 4 = rate of return risky
                                       ! row 5 = womens pension wealth 
                                       ! row 6 = mens pension wealth 
        integer ivart(9,totp,numper)   ! row 1 = children (n)
                                       ! row 2 = married
                                       ! row 3 = work h
                                       ! row 4 = work w
                                       ! row 5 = husband survival
                                       ! row 6 = wife survival
                                       ! row 7 = marriage duration
                                       ! row 8 = formal h
                                       ! row 9 = formal w
        integer fixvar(4,totp)         ! row 1 = hed
                                       ! row 2 = wed
                                       ! row 3 = hage
                                       ! row 4 = type 
        real(8) dpmat(11,3,2),maxsafe(numper)
        character*200::pathoutput



        mptype = 1
        lval = 0.0d-0
print*,"pre evfun"        
        
        
        call EVFUN(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed, &
           fixvar,ivart,rvart,mptype,flik,tprob, &
           alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc,dpmat,regind,numint,intmat,nsimul,&
           maxsafe,lval,pathoutput, counter,estimation,nmoments,loss)
print*,"past evfun"      
  end SUBROUTINE EVMPI

SUBROUTINE EVFUN(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed, &
  fixvar,ivart,rvart,mptype,flik,tprob, &
  alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc,dpmat,regind,numint,intmat,nsimul, &
  maxsafe,criterion,pathoutput,counter,estimation,nmoments,loss)

!XXX should more rigorously separate incidental variables (smarried) from history vectors (snivart(2))

implicit none


! external extract_var,getemax,getstate2,maxu,u,get_eps,GetTheta,GetData,ProbDec,ran1,MSM

integer, parameter :: ncoeff=150
integer, parameter :: numfam=6109
!integer, parameter :: numfam=5314
integer, parameter :: nvar1 = 24 
integer, parameter :: nvar2 = 12 
integer, parameter :: lastage = 90
integer, parameter :: lastworkage = 75
integer, parameter :: numa1cat = 15
integer:: nmoments,wageret2,wageret,wageret_ref,wageret_baseline,hageret,ntyp
real(8) loss(1,1:nmoments), fractionsimtype(ktp)
integer flag,counter,nedraw
integer mptype,ct,o,s1,s2,mpt,ctt,skiploop,numint,idum,ss1,ss2,ii,bidule
integer totp,kidmax,ktp,numper,numiter,nmc,regind(75),intmat(75,75),intmatbis(75,75),tt
integer j,wage,numvalid,update,reform,ttt,ftt,tempref, nompg, nopasis
real(8) vseed(3),maxsafe(numper),typeprob(ktp),rr1(1), supp
real(8) deathprob(91,3)
real(8) getemax,actret,cnu(76,4),cnu_baseline(76,4),cnu_ref(76,4)
real(8) r1(1),r2(1),r3(1),r4(1), divorce_increase
integer kk,jj,xx,zz,nosimul,optch(nmc,7),optch2(1,7),jjj,kkk,nsimul
real(8) flik(totp),tprob(totp)
real(8) parms(numiter),probvec(4)
integer snumk,skidnow
integer vv,typ(totp),ftyp,yearmat(numfam,2),thisyear
real(8) alpha(ncoeff,ktp),alphalow(ncoeff,ktp),alphahigh(ncoeff,ktp),abump(ncoeff,ktp),nalpha(ncoeff,ktp)
real(8) supplement,pmin,pmax,temp111(ktp),temp112,temp113
integer aiter(ncoeff,ktp)
real(8) v(7,7),vbump(7,7),rhomat(7,7),newv(7,7),wexper,fexperw,hexper,fexperh
integer viter(7,7),gg
real(8) temprho,dpmat(11,3,2),experh,experw,mpenwithdraw,fpenwithdraw
integer Sbar,Tbar,T0,i,d(4),numz
integer togv(numz),numv,agev(8),zz1,zz2,maxfam
integer sfert,smarried,sfemale,sdiv,shdec,swdec,totp1,totp2,totp3,nedraw1,nedraw2,nedraw3
real(8) thetap(3,numz),thetam(numper,3,numz),maxearnf,maxearnh
real(8) thetamref(numper,3,numz),thetamreftype(ktp,numper,3,numz),thetamreftype1(numper,3,numz),thetamreftype2(numper,3,numz)
real(8) thetamtype(ktp,numper,3,numz),thetamtype1(numper,3,numz),mintypeprobtemp,mintypeprob,maxtypeprob 
real(8) thetamtype2(numper,3,numz), avgtypeprob
real(8) evec(numper*nmc*7),rsqv(numper,1),rsqcvv(numper,1),rsq,rsqcv
real(8) yevec(nmc*7),thetav(1,numz),sevv(1,numz),sevm(numper,numz)
real(8) umaxmc(nmc),Ztemp(1,numz), optc(nmc)
real(8),allocatable :: emaxv(:,:),Z(:,:),Zbis(:,:),emaxvbis(:,:)
real(8) rvart(6,totp,numper)   ! row 1 = wealth -safe (end of period)
                               ! row 2 = wealth - risky 
                               ! row 3 = rate of return safe
                               ! row 4 = rate of return risky 
                               ! row 5 = womens pension wealth
                               ! row 6 = mens pension wealth
integer ivart(9,totp,numper)   ! row 1 = children (n)
                               ! row 2 = married
                               ! row 3 = work h
                               ! row 4 = work w
                               ! row 5 = survivh
                               ! row 6 = survivw
                               ! row 7 = marriage duration
                               ! row 8 = formal h
                               ! row 9 = formal w
integer fixvar(4,totp)         ! row 1 = hed
                               ! row 2 = wed
                               ! row 3 = hage (age difference with wife)
                               ! row 4 = race 
real(8) nrvart(6,1,numper)     ! row 1 = wealth - safe assets (end of period)
                               ! row 2 = wealth - risky assets
                               ! row 3 = rate of return safe
                               ! row 4 = rate of return risky
                               ! row 5 = womens pension 
                               ! row 6 = mens pension wealth 
integer nivart(9,1,numper)      ! row 1 = children (n)
                               ! row 2 = married
                               ! row 3 = work h
                               ! row 4 = work w
                               ! row 5 = survivh
                               ! row 6 = survivw
                               ! row 7 = marriage duration 
                               ! row 8 = formal h
                               ! row 9 = formal w
integer nfixvar(4,1)           ! row 1 = hed
                               ! row 2 = wed
                               ! row 3 = hage
                               ! row 4 = female sampled 
integer sivart(9,numfam,numper),sfixvar(4,numfam),hed,wed
integer snivart(9,1,numper),snfixvar(4,1),tempnk,templsh,templsw,year(numfam),syear
integer allstate(numfam,numper),found,idvec(numfam)
real(8) srvart(6,numfam,numper),sumax(1),snrvart(6,1,numper),hhinc, bono,tao,tao_baseline, tao_ref
real(8) sevec(numper*7),ssevec(7),optsafe(nmc),optrisky(nmc),soptsafe(1),soptrisky(1)
real(8) searnh(1),searnw(1),earnh(nmc),earnw(nmc),simpensionw,simpensionh,simwealth
integer earnhi,earnwi
integer hage,loopstart
real(8) optpenw(nmc),soptpenw(1)
real(8) optpenm(nmc),soptpenm(1)
integer, allocatable :: simdata1(:,:)
real(8), allocatable :: simdata2(:,:)
integer, allocatable :: tsimdata1(:,:)
real(8), allocatable :: tsimdata2(:,:)
integer, allocatable :: simdiv(:),tsimdiv(:),simhdeath(:),simwdeath(:),simtyp(:)
logical, allocatable :: masksimtyp(:)
integer :: clone, ncounter 
integer :: noAPS, nobono, nodivrule, halfPBS,noPBS
character*6:: spec
real(8) :: pasis, mpg,rsqmat(numper,ktp,3)
real(8) criterion                         ! the MSM criterion function returned by subroutine MSM given current parameter vector
integer nsim                              !number of year-household-clone simulated
real(8) tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt
character*200::pathoutput,command,path
character*1::tag
character*2::counterchar
integer:: debug, generous,reformimp,newbest,formalh,formalw,informalh,informalw,married,evermarried,female,numk,silent,hw,ww,n(numper),group,ind
real(8) pensionh,pensionw,hcaph,hcapw,wealth,cmin        
integer::counterfactual,time_array_0(8),time_array_1(8),parabump,estimation,nyrs,highcontrib,highPBS,univPBS,hretage
real::	start_time, finish_time

!---------------------------------------------------------------------------------------------
! Specifications
!---------------------------------------------------------------------------------------------
print*,"evfun"
reform = 1
counterfactual = 7 

silent = 0
skiploop = 0
nosimul = 0
nyrs = 26				!number of calendar years to simulate ahead 
maxfam = 9999999
!nmc = 2
!totp = 300
!nedraw = 300

ncounter = 1
spec = ""
path = ""
pathoutput=adjustl(trim(path))//"outputs/"


! Policy parameters

	tao_baseline = 0.10d-0
	tao_ref = 0.10d-0
	pasis = 565236.0d-0 ! this correspond to the value of PASIS for 70-75 years old in august 2007
	mpg = 2*pasis
	pmin = 12.0d-0*75000.0d-0
	pmax = 12.0d-0*255000.0d-0
	hageret = 65
	wageret_baseline = 60
	wageret_ref = 60
	wageret2 = 65
	bono = tao_ref*1.5d-0*165000.0d-0*12.0d-0
	
	cmin = alpha(85,1)*1000000.0d-0    ! set a minimum consumption level
        actret = alpha(104,1)		   ! interest rate used in the model solution
	nodivrule = 0
	nobono = 1
	noAPS = 0

if (counter.gt.ncounter) goto 88



! Initialize and allocate variables

	write(counterchar,fmt='(I2)') counter
print*,"numfam=",numfam,"nyrs=",nyrs,"nsimul=",nsimul,"nvar1=",nvar1	
	allocate(simdata1((numfam*nyrs*nsimul),nvar1))
print*,"simdata1"	
	allocate(simdata2((numfam*nyrs*nsimul),nvar2))
print*,"simdatar2"	
	allocate(simdiv(numfam*nyrs*nsimul))
print*,"simdiv"	
	allocate(simhdeath(numfam*nyrs*nsimul))
	allocate(simwdeath(numfam*nyrs*nsimul))
print*,"simdeath"	
	allocate(simtyp(numfam*nyrs*nsimul))
print*,"simtyp"	
	allocate(masksimtyp(numfam*nyrs*nsimul))
print*,"allocated sim arrays"
	simtyp = 0 
	simdiv = 0
	flag = 0 
	clone=1
	thisyear=-99
	i = 0
	bidule = 0
	nalpha = alpha
	newv = v
	thetam = 0.0d-0
	thetamtype = 0.0d-0
	thetamreftype = 0.0d-0
	Rsqmat=1.0d-0

	call ProbDec(deathprob,cnu_baseline)
	cnu_ref=cnu_baseline
	if (counterfactual.eq.6.or.counterfactual.eq.10) cnu_ref(:,3)=cnu_baseline(:,2)  
	vv=23

	call extract_var(v,vbump,viter,alpha,abump,aiter,vv,ktp)
	! store interactions in a more compact form, keep original interaction array as "intmatbis", to be used in GetTogV
	intmatbis=intmat 
	ind=0
	intmat=-9
	do i=1,75
	do j=i,75
		if (intmatbis(i,j)==1) then
    			ind=ind+1
    			intmat(ind,1)=i
    			intmat(ind,2)=j
		endif
	enddo 
	enddo 
   
! Adjust parameters for counterfactuals
 
	if (reform.eq.0) then
		tao=tao_baseline
		wageret=wageret_baseline
		cnu=cnu_baseline
	elseif (reform.eq.1) then
		tao=tao_ref
		wageret=wageret_ref
		cnu=cnu_ref
	endif





! Draw shocks

	idum = dint(-1.0d-0*abs(vseed(3)))
	evec = 0.0d-0
	call get_eps(evec,v,7,nmc,idum,numper)
	write(tag,fmt='(I1)') counter
	allocate(emaxv(totp,1))
	allocate(Z(totp,numz))
print*,"allocated emax"


!----------------------------------------------------------------------
! Model solution (stored in thetam)
!----------------------------------------------------------------------


thetam=-99.0d-0 		!This period's emax coefficients
thetap=0.0d-0			!Last period's emax coefficients


if (skiploop.eq.0) then

	if (silent.ne.1) open(unit=981,action='write',file=adjustl(trim(pathoutput))//'R2'//adjustl(trim(counterchar))//".txt",position='rewind')



    	TTTLOOP: do ttt = 1,ktp  ! Loop over types

    	mptype = ttt 			! XXX difference between mptype and ttt??
    
     
        YRLOOP: do i = Tbar,T0,-1 

            if (silent.ne.1) write(*,*) 'working on age',i
            ct = i-15

            ! Initialize variables for Emax regression
            ! if not the last period, then get the thetas that were calculated last period 

            if (i.lt.Tbar) then
               thetap(:,1:numz) = thetam((ct+1),:,1:numz) 
            else
               thetap(:,1:numz)=0.0d-0
            end if

            Z = 0.0d-0
            emaxv = 0.0d-0
            maxearnf = 0.0d-0
            maxearnh = 0.0d-0
 
            FAMLOOP: do o = 1,totp 
            ! if (o.gt.15) cycle famloop


		! Initialize variables

                probvec = 0.0d-0
                mpenwithdraw = 0.0d-0
                fpenwithdraw = 0.0d-0
                umaxmc=0.0d-0
                hage = fixvar(3,o)+ct-1
                nfixvar(1:4,1) = fixvar(1:4,o)				!reduce arrays of state variables to only the current family 
                nivart(1:9,1,1:numper) = ivart(1:9,o,1:numper)
                nrvart(1:6,1,1:numper) = rvart(1:6,o,1:numper)
                optch = 0
                optsafe = 0.0d-0
                optrisky = 0.0d-0
                Ztemp = 0.0d-0


                ! specify subvectors and matrices containing the epsilons
                ! for current history under consideration
                ! will use same draws
                ! for all histories, but different draws every time period
                ! in order to impose iid assumption

                s1 = (nmc*7)*(ct-1)+1	  ! indices for monte carlo draws
                s2 = (s1-1)+nmc*7
                yevec = evec(s1:s2)

                ! call the subroutine that gets the max utility for
                ! all the possible choices and all the epsilon MC draws
                      		
                call maxu(o,ct,clone, thisyear, umaxmc,nmc,nfixvar,nivart,nrvart,yevec,nalpha,thetap,Tbar,numz,ktp,totp,mptype,dpmat,intmat,optch,optsafe,optrisky,earnh,earnw,kidmax,experh,experw,maxearnf,maxearnh,numa1cat,deathprob,cnu,optpenw,optpenm,probvec,maxsafe,mpenwithdraw,fpenwithdraw,reform,actret, nompg, mpg, nopasis, pasis,nobono,bono,nodivrule, pmax,pmin,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt,optc,tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual)

                !  Integrate (get emaxv) by taking the mean over the estimated umax values
                !  for the different monte carlo draws for the epsilons	for the different
                !  types
                
                emaxv(o,1) = sum(umaxmc)/(1.0d-0*nmc)
    
                ! get vector of emax regressors: functions of state variables at beginning of "ct" period
                
                fexperh = 0.0d-0
                fexperw = 0.0d-0
                if (ct.gt.1) then
                   fexperh = sum(nivart(8,1,1:ct))*1.0d-0   !male formal sector work experience
                   fexperw = sum(nivart(9,1,1:ct))*1.0d-0   !woman formal sector work experience
                end if
                numk = sum(nivart(1,1,1:ct))*1.0d-0 
                call getstate2(o,ct,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 

                Z(o,1:numz)=Ztemp(1,1:numz)


            end do FAMLOOP

            !After processed all the histories, run regression to approximate 
            !the emax functions.

            thetav = 0.0d-0      ! vector that will contain estimated thetas
            sevv = 0.0d-0        ! vector that will contain estimated theta standard 
                                 ! errors
            togv = 0   	         ! vector that indicators whether the column of theta
                                 ! is nonzero
            rsq = 0.0d-0	     ! variable containing R-squared from the regression
            rsqcv = 0.0d-0	     ! variable containing R-squared from the cross-val.
                                 ! regression
            
            do group=1,3
           ! figure out which Z variables to turn on or off in the emax
           ! regression

           totp1=int(totp/2)
           totp2=int(totp/4)                 
           totp3=totp-totp2-totp1
           nedraw1=int(totp1*nedraw/totp)
           nedraw2=int(totp2*nedraw/totp)
           nedraw3=nedraw-nedraw1-nedraw2
          

           ! Run the regression, store results in thetav,sevv,rsq,rsqcv
           ! Separate regressions are run for each type (but state variables
           ! will be the same for all types as long as do not include 
           ! type-specific intercepts)
       
           if (group==1) then !married
               allocate(Zbis(totp1,numz))
               allocate(emaxvbis(totp1,1))
               Zbis=Z(1:totp1,:)
               emaxvbis(:,1)=emaxv(1:totp1,1)
               call GetTogV(ct,Zbis,togv,numz,intmatbis,regind,totp1,numint,numv,silent)
  	           call GetTheta(o,ct,Zbis,emaxvbis,thetav(1,:),togv,sevv,rsq,nedraw1,totp1,numz,numv,flag,pathoutput,group,ttt,reform,silent)
          elseif (group==2) then !single female
  	           allocate(Zbis(totp2,numz))
               allocate(emaxvbis(totp1,1))
               Zbis=Z(totp1+1:totp1+totp2,:)
               emaxvbis(:,1)=emaxv(totp1+1:totp1+totp2,1)
               call GetTogV(ct,Zbis,togv,numz,intmatbis,regind,totp2,numint,numv,silent)
  	           call GetTheta(o,ct,Zbis,emaxvbis,thetav(1,:),togv,sevv,rsq,nedraw2,totp2,numz,numv,flag,pathoutput,group,ttt,reform,silent)
           elseif  (group==3) then !single male
  	           allocate(Zbis(totp3,numz))
               allocate(emaxvbis(totp3,1))
               Zbis=Z(totp1+totp2+1:totp,:)
               emaxvbis(:,1)=emaxv(totp1+totp2+1:totp,1)
               call GetTogV(ct,Zbis,togv,numz,intmatbis,regind,totp3,numint,numv,silent)
 	           call GetTheta(o,ct,Zbis,emaxvbis,thetav(1,:),togv,sevv,rsq,nedraw3,totp3,numz,numv,flag,pathoutput,group,ttt,reform,silent)
          endif
       deallocate(Zbis)
       deallocate(emaxvbis)
       if (silent.ne.1) write(*,*)'ct',counter, 'age',(ct+15), 'tp',ttt,'grp',group,'ref',reform
       if (silent.ne.1) write(981,fmt=36) 'age',(ct+15), 'R-square type',ttt,'group',group,'reform',reform,':',rsq
   36 format(A3,I4,A15,I4,A7,I4,A7,I4,A3,f10.7)   
       thetam(ct,group,:) = thetav(1,:) 
       Rsqmat(ct,ttt,group)=rsq
        enddo

        !   sevm(ct,:) = sevv(1,:)
        !   rsqv(ct,1)=rsq
        !   rsqcvv(ct,1)=rsqcv
        
        call date_and_time(values=time_array_1)
        finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60+ time_array_1 (7) + 0.001 * time_array_1 (8)
      !  write (*,*) 'elapsed wall clock time:',finish_time    

 
        end do YRLOOP
       
        ! store the thetas for each type
    
        if (reform.eq.0) then
            thetamtype(ttt,:,:,:) = thetam
	    if (ttt.eq.1) then
               thetamtype1 = thetam
            else
               thetamtype2 = thetam
            end if 
        end if
        if (reform.eq.1) then
            thetamreftype(ttt,:,:,:) = thetam
            if (ttt.eq.1) then
               thetamreftype1 = thetam
            else
               thetamreftype2 = thetam
            end if 
        end if
        if (estimation.ne.1) call ReportRes(thetam,sevm,rsqv,rsqcvv,mptype,Sbar,Tbar,numz,numint,totp,totp,intmatbis,reform,pathoutput,ktp)

    end do TTTLOOP

deallocate(emaxv)
deallocate(Z)

close(981)
end if !skiploop eq 0

!--------------------------------------------------------
! read in the
! thetam matrix from a previous run of the program and do
! not recalculate them
if (estimation.ne.1) call readtheta2(silent,pathoutput,reform,thetamtype,thetamtype1,thetamtype2,thetamreftype,thetamreftype1,thetamreftype2,numper,numz,ktp)

print*,"readtheta"
!----------------------------------------------------------------------
! Simulate sample decisions

!      open (file='checksav.txt', unit=9999)

nmc = 1
update = 1 
mintypeprob=0.5
maxtypeprob=0.5

if (nosimul.eq.0) then

    loopstart = 0
    
    ! write to separate files depending on whether simulating
    ! with or without the reform

    if (silent.eq.0) then
        open(unit=700,file=adjustl(trim(pathoutput))//'simdata'//adjustl(trim(spec))//'.asc',RECL=300,status='REPLACE',form=&
        &'FORMATTED', action='WRITE')
        rewind(unit=700)
    end if

    
    ! initialize idum for drawing the random shocks
    idum = dint(-1.0d-0*abs(vseed(2)))
    
    ! call the subroutine that reads in the data and keeps track of what years have complete state variables
    srvart = 0.0d-0
    sivart = 0 
    sfixvar = 0
    allstate = 0
    numvalid = 0     !number of valid cases in dataset that will be processed
    call GetData(sivart,srvart,sfixvar,allstate,idvec,yearmat,numvalid,kidmax,T0,Tbar) 
    

    s1loop: do s1 = 1,numfam

    if (s1.gt.maxfam) cycle s1loop
    
    if (silent.ne.1) write(*,*) 'family #',s1
         
    ! create family specific state vectors

    snfixvar(1:4,1) = sfixvar(1:4,s1)
    snivart(1:9,1,1:numper) = sivart(1:9,s1,1:numper)  
    snrvart(1:6,1,1:numper) = srvart(1:6,s1,1:numper)     

    cloneloop: do clone=1,nsimul

        ! Find initial year: move forward to a time period (wife's age -15 or single man's age -15) when the state variables are complete
        jj = 1
        do while(allstate(s1,jj).ne.1.and.jj.lt.(Tbar-15)) 
            jj = jj + 1 
        end do  
        ct = jj ! ct is the initial household age -15 
              
        ! Determine the family's permanent or initial characteristics 

        hed = snfixvar(1,1)
        wed = snfixvar(2,1)
        sfemale = snfixvar(4,1)        
        smarried = snivart(2,1,ct)
        if (smarried.eq.0) then
            if (sfemale.eq.1)  snivart(5,1,ct:Tbar-15)=0 !If single, set spouse as dead
            if (sfemale.eq.0)  snivart(6,1,ct:Tbar-15)=0
        endif
        if (smarried.eq.-99) smarried=0
        if (smarried==1) then 
            wage = ct+15
            hage = min0(snfixvar(3,1)+ct-1,90)
        elseif (smarried==0.and.sfemale==0) then
            wage=-99
            hage= ct+15
        elseif (smarried==0.and.sfemale==1) then
            wage=ct+15
            hage=-99
        endif

        ! Compute type probability for this family => type should be included in the fixvar vector XXX
        
        call ran1(idum,rr1,1)
        
	do ftt=1,ktp
        	temp111(ftt) =  alpha(75,ftt)+&
                    		alpha(78,ftt)*dble(smarried)+&
                    		alpha(79,ftt)*dble(hage)*dble(1-sfemale)+&
                    		alpha(79,ftt)*dble(wage)*dble(sfemale)+&
                    		alpha(76,ftt)*dble(wed)*dble(smarried)+ &
                        	alpha(76,ftt)*dble(wed)*dble(sfemale)*(1.0d-0-dble(smarried)) + &
                        	alpha(77,ftt)*dble(hed)*dble(smarried) + &
                    		alpha(77,ftt)*dble(hed)*(1.0d-0-dble(sfemale))*(1.0d-0-dble(smarried))
        	temp111(ftt)= dexp(temp111(ftt))
	enddo
	
	temp112=0.0d-0
	do ftt=1,ktp
		temp112 = temp112 + temp111(ftt)
	enddo
	do ftt=1,ktp
		typeprob(ftt) = temp111(ftt)/temp112
	enddo


        mintypeprobtemp=minval(typeprob(:))
        mintypeprob=dmin1(mintypeprobtemp,mintypeprob)
        
	! Simulate this family's type
        temp113=0.0d-0
	do ftt=1,ktp 
        	temp113=temp113+typeprob(ftt)
		if (rr1(1).le.temp113) then
            		mptype=ftt
        	exit
		endif
	enddo      
        

        ! Draw a lifetime of random shocks for this family 

        sevec = 0.0d-0
        call get_eps(sevec,v,7,1,idum,numper)
        ssevec = 0.0d-0
        
        ! Initialize optimal choice vector
        
        optch2 = -99  
                    !                    optch(kk,1) = ww 
                    !                    optch(kk,2) = hw
                    !                    optch(kk,3) = fert
                    !                    optch(kk,4) = a1
                    !                    optch(kk,5) = a2
                    !                    optch(kk,6) = formh 
                    !                    optch(kk,7) = formw 

        ! Loop over the time periods
        
        ! for each time period when the state space is complete, need to
        ! make predictions until the next time period the state space is complete.

        if (allstate(s1,jj).eq.1) then 

            kkloop: do kk = ct,ct+nyrs-1    ! do 2004-2014. kk is the current hh age -15
                
                thisyear = (kk-yearmat(s1,2)+yearmat(s1,1))      
                               
               ! cycle cohorts that are less than T0 in 2004 
		
		if (ct+15.lt.T0) cycle kkloop

		 ! cycle dead households (i.e. the drawn individual in the household is dead)
                
                if (kk.gt.75) cycle kkloop
                if ((sfemale.eq.1.and.snivart(6,1,kk).eq.0).or.(sfemale.eq.0.and.snivart(5,1,kk).eq.0)) cycle kkloop

                ! only simulate the last cohort (birthyear >1963) after 2019
                
                !if (thisyear-kk-15.lt.1963.and.thisyear.gt.2019) cycle kkloop
                
                ! Set year specific variables:
                
                smarried=snivart(2,1,kk)
                
                !tempref denotes whether the reform is currently in place
                if (thisyear.lt.2009.or.reform.eq.0) then 
                      tempref = 0
                else 
                      tempref = 1 
                end if

             
		if (tempref==1) then
                	select case (thisyear) 
                	 case (2009)
                	    pmin=12.0d-0*67500.0d-0
                	    pmax=12.0d-0*102500.0d-0
                	 case (2010)
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*227500.0d-0
               		  case (2011)
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*175000.0d-0
                	 case (2012,2013,2014,2015,2016)
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*255000.0d-0
                	end select

                	if (thisyear.ge.2016) then
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*255000.0d-0
                	endif
			
			wageret=wageret_ref
			tao=tao_ref
			cnu=cnu_ref
		elseif (tempref==0) then
			wageret=wageret_baseline
			tao=tao_baseline
			cnu=cnu_baseline
		endif
                
                if (thisyear.eq.2004) then
                  actret = 0.1055d-0
                else if (thisyear.eq.2005) then
                  actret = 0.0886d-0
                else if (thisyear.eq.2006) then
                  actret = 0.0458d-0
                else if (thisyear.eq.2007) then
                  actret = 0.1577d-0
                else if (thisyear.eq.2008) then
                  actret = -0.1998d-0
                else if (thisyear.eq.2009) then
                  actret = 0.1770d-0
                else if (thisyear.eq.2010) then
                  actret = 0.063d-0-0 
                end if 

                ! get the appropriate thetap vector and optimal choices
                ! if pre 2009 and simulating with reform, using thetam, otherwise use thetamref


                    if (tempref.eq.0) then
                        if ((kk+15).lt.Tbar) then
                            thetap(:,1:numz) = thetamtype(mptype,(kk+1),:,1:numz)
                        else
                            thetap(:,1:numz)=0.0d-0
                        end if  !(kk+15) lt Tbar
                    else      ! use post-reform thetas for emax 
                        if ((kk+15).lt.Tbar) then
                            thetap(:,1:numz) = thetamreftype(mptype,(kk+1),:,1:numz)
                        else
                            thetap(:,1:numz)=0.0d-0
                        end if  !(kk+15) lt Tbar
                    end if


                ! Extract montecarlo draws

                ss1 = 7*(kk-1)+1	  ! indices for monte carlo draws
                ss2 = ss1+6   
                ssevec(1:7)=sevec(ss1:ss2)

                soptsafe = 0.0d-0     ! initialize variables that will store optimal choices of safe and
                soptrisky = 0.0d-0    ! risky
                call maxu(s1,kk,clone,thisyear, sumax,nmc,snfixvar,snivart,snrvart,ssevec,nalpha,thetap,Tbar,numz,ktp,numfam,mptype,dpmat,intmat,optch2,soptsafe,soptrisky,searnh,searnw,kidmax,experh,experw,maxearnf,maxearnh,numa1cat,deathprob,cnu,soptpenw,soptpenm,probvec,maxsafe,mpenwithdraw,fpenwithdraw,tempref,actret, nompg, mpg, nopasis,pasis, nobono,bono, nodivrule, pmax,pmin,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt,optc,tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual)

                ! update the state variables to reflect the optimal choices
                
                if (thisyear.eq.2004) then
                    hexper = 0.0d-0
                    fexperh = 0.0d-0
                    wexper = 0.0d-0
                    fexperw = 0.0d-0
                    ggloop: do gg = 1,(kk-1)
                         hexper = hexper+snivart(3,1,gg)
                         fexperh = fexperh+snivart(8,1,gg)
                         if (snivart(4,1,gg).eq.1) then
                           wexper = wexper+1.0d-0
                           fexperw = fexperw+snivart(9,1,gg)
                         elseif (snivart(4,1,gg).eq.2) then
                           wexper = wexper+0.5d-0
                           fexperw = fexperw+0.5d-0*snivart(9,1,gg)
                         else
                           wexper = wexper+0.0d-0
                           fexperw = fexperw+0.0d-0
                         end if
                    end do ggloop
                 end if
                
                if (sfemale.eq.0.or.smarried.eq.1) then
                    hexper = hexper+optch2(1,2)
                    fexperh = fexperh+optch2(1,2)*optch2(1,6)
                else
                    hexper = -99
                    fexperh= -99
                endif
                
                if (sfemale.eq.1.or.smarried.eq.1) then
                if (optch2(1,1).eq.1) then
                    wexper = wexper+optch2(1,1)
                    fexperw = fexperw+optch2(1,7)*optch2(1,1)
                else if (optch2(1,1).eq.2) then
                    wexper = wexper+0.5d-0
                    fexperw = fexperw+optch2(1,7)*0.5d-0
                else
                    wexper = wexper+0.0d-0
                    fexperw = fexperw+0.0d-0
                end if
                else
                    wexper = -99
                    fexperw = -99
                endif 

                snivart(3,1,kk)  = optch2(1,2)   !work status of husband
                snivart(4,1,kk)  = optch2(1,1)   !work status of wife
                snivart(8,1,kk)  = optch2(1,6)   !formal work status of husband
                snivart(9,1,kk)  = optch2(1,7)   !formal work status of wife
                
                ! set next period as having a complete state space (since simulating many years forward)
                ! XXX figure this out better
 
                if ((kk+1).lt.65) then
                    allstate(s1,(kk+1)) = 1
                end if 

                simpensionw = soptpenw(1) !end of pension balance at the optimal choice
                simpensionh = soptpenm(1)
                simwealth = soptsafe(1) 

                ! Simulate stochastic transitions
                
		if ((kk+1).le.Tbar-15) then 
		if (sfemale.eq.1.and.snivart(5,1,kk).eq.1) snivart(6,1,kk+1)=1 ! by default, sampled individual stays alive until next period
		if (sfemale.eq.0.and.snivart(6,1,kk).eq.1) snivart(5,1,kk+1)=1 
		endif

		! Absent divorce or death the household structure is preserved
		if (smarried.eq.1.and.(kk+1).le.(Tbar-15)) then
                    snivart(2,1,(kk+1))=1 ! still married
		    snivart(5,1,(kk+1))=1 ! still live husband		
		    snivart(6,1,(kk+1))=1 ! still live wife 		
                end if
                
                ! simulate divorce
                
		sdiv = 0
                if (smarried.eq.1.and.wage.lt.60) then !only allow divorce until wife turns 60
                    r1=0.0d-0
		    call ran1(idum,r1,1)
					divorce_increase=0.00d-0
					if (thisyear.gt.2008.and.reform.eq.1) divorce_increase=0.02d-0
                    if (r1(1).lt.probvec(1)+divorce_increase) then !adjust in simulation for the fact that divorce rates increased (unanticipated)
                         sdiv = 1
                         snivart(2,1,(kk+1):(Tbar-15))=0 ! if divorce, household becomes unmarried
                         if (sfemale.eq.1) snivart(5,1,(kk+1):(Tbar-15))=0 !if divorce and female sampled, hh no longer contains live husband
                         if (sfemale.eq.0) snivart(6,1,(kk+1):(Tbar-15))=0 !if divorce and male sampled, hh no longer contains live wife
                         simwealth = 0.5d-0*simwealth   !get half of private assets
                         if (tempref.eq.1.and.nodivrule.eq.0) then !if divrule applies
                            simpensionh=(simpensionw+simpensionh)/2.0d-0  !after reform pension is split between spouses
                            simpensionw=(simpensionw+simpensionh)/2.0d-0  !after reform pension is split between spouses
                         end if
                    end if !r1
                end if !smarried

                ! Simulate whether a member of the household dies
                
                shdec = 0
                swdec = 0
                r1=0.0d-0 
                call ran1(idum,r1,1)
                r2=0.0d-0 
                call ran1(idum,r2,1)

                ! Husband dies
                if (snivart(5,1,kk).eq.1) then ! alive male in the household
                if (r1(1).lt.probvec(3).or.kk.eq.(Tbar-15)) then ! male dies
                    shdec = 1
                    snivart(5,1,(kk+1):Tbar-15)=0   !male will be dead
                    snivart(2,1,(kk+1):Tbar-15)=0   !household will be unmarried
                    if (smarried.eq.1) then
                        simpensionw = simpensionw+0.6d-0*simpensionh  !if husband dies, wife inherits 0.6 of his pension
                        simpensionh = 0.0d-0
                    endif
                else !r1
                    shdec = 0 
                    snivart(5,1,kk+1) = 1
                endif
                endif
                
		! Wife dies
                if (snivart(6,1,kk).eq.1) then !alive female in the household
                if (r2(1).lt.probvec(4).or. kk.eq.(Tbar-15)) then
                    swdec = 1
                    snivart(6,1,(kk+1):Tbar-15)=0   !female will be dead
                    snivart(2,1,(kk+1):Tbar-15)=0   !household will be unmarried
                    if (smarried.eq.1) then                             
                        simpensionh = simpensionh  !if wife dies, husband does not inherit her pension (gets own)
                        simpensionw = 0.0d-0
                    endif
                else
                    swdec = 0 
                    snivart(6,1,kk+1) = 1
                end if  !r2
                endif

                !draw fertility 
                snumk=sum(snivart(1,1,1:kk))
                skidnow = 0   ! simulated fertility
                if ((sfemale.eq.1.or.smarried.eq.1).and.kk.le.25.and.snumk.lt.kidmax) then
                   r1=0.0d-0
                   call ran1(idum,r1,1)
                   if (r1(1).lt.probvec(2)) then
                        snivart(1,1,(kk+1)) = 1
                        skidnow = 1 
                   else 
                        snivart(1,1,(kk+1)) = 0
                        skidnow = 0 
                   end if  !r1
                end if  !sfemale 


                ! Give out bono por hijo if eligible
                supplement = 0.0d-0
                if (reform.eq.1.and.thisyear.eq.2008.and.nobono.eq.0) then
                    supplement = bono*dble(snumk)   !pension reform fertility supplement
                end if            
                if (reform.eq.1.and.thisyear.gt.2008.and.nobono.eq.0) then
                    supplement = bono*dble(skidnow)
                end if 

                simpensionw = simpensionw + supplement
                
                ! fill in values of simdata1 and simdata2 and output simulation to text file
                if (loopstart.eq.0) then
                
                	simdata1=-99
                	simdata2=-99
                	simdiv=-99
                	simtyp=-99
                	simhdeath=-99
                	simwdeath=-99
     
                	nsim = 1
                	!     write(700,*) 'idvec thisyear married female wage hage shw sww sformh sformw savrate simpensionw simpensionh simwealth numk kidnow hexper fexperh wexper fexperw earnh earnw mpenwithdraw fpenwithdraw wed hed sdiv tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt'
                end if

                simwealth=dble(int(simwealth)*1000.0d-0)/1000.0d-0
                optc(1)=optc(1)/1000000
                supp=supplement/1000.0d-0
                                          
                simdata1(nsim,1) = idvec(s1)
                simdata1(nsim,2) = thisyear
                simdata1(nsim,3) = smarried 
                simdata1(nsim,4) = sfemale 
                simdata1(nsim,5) = wage
                simdata1(nsim,6) = hage
                simdata1(nsim,7) = optch2(1,2)
                simdata1(nsim,8) = optch2(1,1)
                simdata1(nsim,9) = optch2(1,6)
                simdata1(nsim,10) = optch2(1,7)
                simdata2(nsim,1) = optch2(1,4)*1.0d-0 
                simdata2(nsim,2) = simpensionw
                simdata2(nsim,3) = simpensionh
                simdata2(nsim,4) = simwealth
                simdata1(nsim,11) = snumk
                simdata1(nsim,12) = skidnow
                simdata2(nsim,5) = hexper 
                simdata2(nsim,6) =  fexperh
                simdata2(nsim,7) =  wexper
                simdata2(nsim,8) =  fexperw
                simdata2(nsim,9) =  searnh(1)
                simdata2(nsim,10) = searnw(1)
                simdata2(nsim,11) =  mpenwithdraw
                simdata2(nsim,12) =  fpenwithdraw
                simdata1(nsim,13) = hed
                simdata1(nsim,14) = wed
                simdiv(nsim) = sdiv
                simtyp(nsim) = mptype
		simhdeath(nsim) = shdec
		simwdeath(nsim) = swdec
                

                if (silent.eq.0.and.nsim>2) write(700,fmt=26) idvec(s1),thisyear,clone, mptype, smarried,sfemale, wage,hage,optch2(1,2),optch2(1,1),optch2(1,6), &
                  optch2(1,7),optch2(1,4),simpensionw,simpensionh,simwealth,snumk,skidnow,hexper,fexperh,wexper,fexperw, &
                searnh(1),searnw(1),mpenwithdraw,fpenwithdraw,hed,wed,sdiv,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt, &
                  pasiscostopt, supp, sumax(1), optc(1),swdec,shdec ,snivart(5,1,kk),snivart(6,1,kk)

26      format(i7,i6, i3, i3, i5,i5,i5,i5,i2,i2,i2,i2,i4,3f14.2,i5,i5,4f8.2,4f14.2,i5,i5,i5,7f9.2,f16.2, f8.2)

                ! Update ages, husband stops aging at 90
                
                if (sfemale.eq.1.or.smarried.eq.1) wage = kk+1+15
                if (sfemale.eq.0.or.smarried.eq.1) hage = min0(snfixvar(3,1)+kk-1+1,90)
                
		! Update state variables XXX eventually do all the state variables here and only use hte temp variables elsewhere in the kkloop
		! Will create halive and walive for snivart(5 and snivart(6
	
                snrvart(5,1,kk) = simpensionw 
                snrvart(6,1,kk) = simpensionh
                snrvart(1,1,kk) = simwealth
		

		loopstart = 1
                nsim = nsim + 1

            end do kkloop

        end if  !allstate 
        end do cloneloop
        
   ! Update full state variable vector using family specific state vectors

   sfixvar(1:4,s1) = snfixvar(1:4,1)
   sivart(1:9,s1,1:numper) = snivart(1:9,1,1:numper)  
   srvart(1:6,s1,1:numper) = snrvart(1:6,1,1:numper)
   
   end do s1loop
        

   ! Rule out near-degenerate type mixtures

   if (mintypeprob.gt.0.95.or.maxtypeprob.lt.0.05) flag=1
   do ttt=1,ktp
	do i=1,nsim
 	masksimtyp(i)=simtyp(i).eq.ttt
	enddo
	ntyp=count(masksimtyp)
	fractionsimtype(ttt)=(ntyp*1.0d-0)/((nsim-1)*1.0d-0)
   	if (silent.eq.0) write(*,*)"Fraction of type",ttt,":",fractionsimtype(ttt)
	if (fractionsimtype(ttt).lt.0.1.or.fractionsimtype(ttt).gt.1.9) flag=1
   enddo
   !avgtypeprob=sum(simtyp(1:(nsim-1)))*1.0d-0
   !avgtypeprob=avgtypeprob/(nsim-1)*1.0d-0
   !if (avgtypeprob.lt.1.1.or.avgtypeprob.gt.1.9) flag=1     
     

close(unit=700)

nsim=nsim-1

allocate(tsimdata1(nsim,nvar1))
allocate(tsimdata2(nsim,nvar2))
allocate(tsimdiv(nsim))
print*,"allocated tsimdata"
tsimdata1(1:nsim,1:nvar1) = simdata1(1:nsim,1:nvar1)
tsimdata2(1:nsim,1:nvar2) = simdata2(1:nsim,1:nvar2)
tsimdiv(1:nsim)=simdiv(1:nsim)

tsimdata2(:,4)=tsimdata2(:,4)/1000000
tsimdata2(:,9)=tsimdata2(:,9)/1000000
tsimdata2(:,10)=tsimdata2(:,10)/1000000

call MSM(nsim,nmoments,nvar1,nvar2,criterion,tsimdata1,tsimdata2,tsimdiv,loss,pathoutput,estimation,counterchar)
!if (silent.ne.1) print*,"criterion", criterion
!if (silent.ne.1.or.estimation.eq.1) print*,minval(Rsqmat(20:75,:,:))
if (flag.eq.1) criterion=9999999999999

deallocate(tsimdata1)
deallocate(tsimdata2)
deallocate(tsimdiv)
end if  !simul eq 1


deallocate(simdata2)
deallocate(simdiv)
deallocate(simhdeath)
deallocate(simwdeath)
deallocate(simdata1)
deallocate(simtyp)
deallocate(masksimtyp)
print*,"end of evfun"


88      end subroutine EVFUN


!------------------------------------------------------------------------------

subroutine extract_var(v,vbump,viter,alpha,abump,aiter,vv,ktp)

integer, parameter :: ncoeff=150

integer ktp
real(8) v(7,7),vbump(7,7),alpha(ncoeff,ktp)
real(8) rhomat(7,7),abump(ncoeff,ktp)
integer vv, i,aiter(ncoeff,ktp),viter(7,7),j



vv = 105 
v(1,1) = alpha(vv,1)
v(2,2) = alpha(vv+1,1)
v(3,3) = alpha(vv+2,1)
v(4,4) = alpha(vv+3,1)
v(5,5) = alpha(vv+4,1)
v(6,6) = alpha(vv+5,1)
v(7,7) = alpha(vv+6,1)
vbump(1,1) = abump(vv,1)
vbump(2,2) = abump(vv+1,1)
vbump(3,3) = abump(vv+2,1)
vbump(4,4) = abump(vv+3,1)
vbump(5,5) = abump(vv+4,1)
vbump(6,6) = abump(vv+5,1)
vbump(7,7) = abump(vv+6,1)
viter(1,1) = aiter(vv,1)
viter(2,1) = aiter(vv+1,1)
viter(3,1) = aiter(vv+2,1)
viter(4,1) = aiter(vv+3,1)
viter(5,1) = aiter(vv+4,1)
viter(6,1) = aiter(vv+5,1)
viter(7,1) = aiter(vv+6,1)
rhomat(1,2) = alpha(vv+7,1)
rhomat(1,3) = alpha(vv+8,1)
rhomat(1,4) = alpha(vv+9,1)
rhomat(1,5) = alpha(vv+10,1)
rhomat(1,6) = alpha(vv+11,1)
rhomat(1,7) = alpha(vv+12,1)
rhomat(2,3) = alpha(vv+13,1)
rhomat(2,4) = alpha(vv+14,1)
rhomat(2,5) = alpha(vv+15,1)
rhomat(2,6) = alpha(vv+16,1)
rhomat(2,7) = alpha(vv+17,1)
rhomat(3,4) = alpha(vv+18,1)
rhomat(3,5) = alpha(vv+19,1)
rhomat(3,6) = alpha(vv+20,1)
rhomat(3,7) = alpha(vv+21,1)
rhomat(4,5) = alpha(vv+22,1)
rhomat(4,6) = alpha(vv+23,1)
rhomat(4,7) = alpha(vv+24,1)
rhomat(5,6) = alpha(vv+25,1)
rhomat(5,7) = alpha(vv+26,1)
rhomat(6,7) = alpha(vv+27,1)
viter(1,2) = aiter(vv+7,1)
viter(1,3) = aiter(vv+8,1)
viter(1,4) = aiter(vv+9,1)
viter(1,5) = aiter(vv+10,1)
viter(1,6) = aiter(vv+11,1)
viter(1,7) = aiter(vv+12,1)
viter(2,3) = aiter(vv+13,1)
viter(2,4) = aiter(vv+14,1)
viter(2,5) = aiter(vv+15,1)
viter(2,6) = aiter(vv+16,1)
viter(2,7) = aiter(vv+17,1)
viter(3,4) = aiter(vv+18,1)
viter(3,5) = aiter(vv+19,1)
viter(3,6) = aiter(vv+20,1)
viter(3,7) = aiter(vv+21,1)
viter(4,5) = aiter(vv+22,1)
viter(4,6) = aiter(vv+23,1)
viter(4,7) = aiter(vv+24,1)
viter(5,6) = aiter(vv+25,1)
viter(5,7) = aiter(vv+26,1)
viter(6,7) = aiter(vv+27,1)

v(1,2) = rhomat(1,2)*dsqrt(v(1,1))*dsqrt(v(2,2))
v(1,3) = rhomat(1,3)*dsqrt(v(1,1))*dsqrt(v(3,3))
v(1,4) = rhomat(1,4)*dsqrt(v(1,1))*dsqrt(v(4,4))
v(1,5) = rhomat(1,5)*dsqrt(v(1,1))*dsqrt(v(5,5))
v(1,6) = rhomat(1,6)*dsqrt(v(1,1))*dsqrt(v(6,6))
v(1,7) = rhomat(1,7)*dsqrt(v(1,1))*dsqrt(v(7,7))
v(2,1)=v(1,2)
v(3,1)=v(1,3)
v(4,1)=v(1,4)
v(5,1)=v(1,5)
v(6,1)=v(1,6)
v(7,1)=v(1,7)
v(2,3) = rhomat(2,3)*dsqrt(v(2,2))*dsqrt(v(3,3))
v(2,4) = rhomat(2,4)*dsqrt(v(2,2))*dsqrt(v(4,4))
v(2,5) = rhomat(2,5)*dsqrt(v(2,2))*dsqrt(v(5,5))
v(2,6) = rhomat(2,6)*dsqrt(v(2,2))*dsqrt(v(6,6))
v(2,7) = rhomat(2,7)*dsqrt(v(2,2))*dsqrt(v(7,7))
v(3,2)=v(2,3)
v(4,2)=v(2,4)
v(5,2)=v(2,5)
v(6,2)=v(2,6)
v(7,2)=v(2,7)
v(3,4) = rhomat(3,4)*dsqrt(v(3,3))*dsqrt(v(4,4))
v(3,5) = rhomat(3,5)*dsqrt(v(3,3))*dsqrt(v(5,5))
v(3,6) = rhomat(3,6)*dsqrt(v(3,3))*dsqrt(v(6,6))
v(3,7) = rhomat(3,7)*dsqrt(v(3,3))*dsqrt(v(7,7))
v(4,5) = rhomat(4,5)*dsqrt(v(4,4))*dsqrt(v(5,5))
v(4,6) = rhomat(4,6)*dsqrt(v(4,4))*dsqrt(v(6,6))
v(4,7) = rhomat(4,7)*dsqrt(v(4,4))*dsqrt(v(7,7))
v(5,4)=v(4,5)
v(5,3)=v(3,5)
v(4,3)=v(3,4)
v(5,3)=v(3,5)
v(6,3)=v(3,6)
v(7,3)=v(3,7)
v(5,4)=v(4,5)
v(6,4)=v(4,6)
v(7,4)=v(4,7)


i5: do i = 1,7
rhomat(i,i) = 1.0d-0
end do i5
rhomat(2,1)=rhomat(1,2)
rhomat(3,1)=rhomat(1,3)
rhomat(4,1)=rhomat(1,4)
rhomat(5,1)=rhomat(1,5)
rhomat(3,2)=rhomat(2,3)
rhomat(4,2)=rhomat(2,4)
rhomat(5,2)=rhomat(2,5)
rhomat(4,3)=rhomat(3,4)
rhomat(5,3)=rhomat(3,5)
rhomat(5,4)=rhomat(4,5)
rhomat(6,5)=rhomat(5,6)
rhomat(7,5)=rhomat(5,7)
rhomat(7,6)=rhomat(6,7)

end subroutine extract_var

        subroutine GetData(nivart,nrvart,nfixvar,allstate,idvec,yearmat,numvalid,kidmax,T0,Tbar)

        implicit none

! XXX clean up so that counter corresponds to age of sampled individual rather than female age ?

!       Data elements

! household id
! year of observation (e.g. 2004 or 2006)
! number of kids not including current year 
! whether had a kid this period
! whether married (0=no, 1 = yes)
! whether husband working (0 = no, 1 = yes) - code retired as 0
! whether wife working (0=no, 1 = yes, fulltime,2=year, parttime)
! husband work experience in years (coming into period, not including current year)
! wife work experience in years, where part-time work counts as 0.5  (coming into period, not including 2004)
! indicator for whether husband died
! indicator for whether wife died
! indicator for whether divorced
! marriage duration 
! whether husband works in formal sector (1 if yes, 0 if no)
! whether wife works in formal sector  (1 if yes, 0 if no)
! husband's education in years
! wife's education in years
! husband's age in years
! woman's age in years
! whether female (for single households or households that become divorced/widowed)
! nonpension wealth level  (in pesos) - coming into this period
! pension wealth for woman  (in pesos, can be zero if this is a single man HH) - coming into this period
! pension wealth for man  (in pesos, can be zero if this is a single woman HH) - coming into this period
! male earnings
! female earnings
! male formal work last period
! female formal work last period
! male work last period (0 if no, 1 if yes) 
! female work last period (0 if no, 1 if fulltime, 2 if parttime)

        integer, parameter :: numfam = 6109 ! number of households
        integer, parameter :: maxl = 12218 ! number of lines in the input data set 
        !integer, parameter :: numfam = 5314 ! number of households
        !integer, parameter :: maxl = 10628 ! number of lines in the input data set 
        integer, parameter :: numper = 75 


        integer id,nobs,year,hed,wed,race,hage,wage,ss
        integer kidmax
        integer preg,married,singlesamp
        real(8) totwealth,temp,totexph,totexpw
        integer agekid(8),age1,age2,age3,age4,age5,age6,age7,age8,zz
        integer idvec(numfam),numvalid
        integer yearmat(numfam,2)  !observe families up to 3 years - need to keep track
        integer nlines ,i,lastid,j,ctt,k,ct,kk,selfam
        real(8) nrvart(6,numfam,numper)     ! row 1 = wealth - safe assets
                                            ! row 2 = wealth - risky assets
                                            ! row 3 = rate of return safe
                                            ! row 4 = rate of return risky
                                            ! row 5 = women's pension wealth
                                            ! row 6 = men's pension wealth
        integer nivart(9,numfam,numper)   ! row 1 = children (n)
                                       ! row 2 = married
                                       ! row 3 = work h
                                       ! row 4 = work w
                                       ! row 5 = whether male surviving
                                       ! row 6 = whether female surviving
                                       ! row 7 = marriage duration 
                                       ! row 8 = formal h
                                       ! row 9 = formal w
        integer nfixvar(4,numfam)      ! row 1 = hed
                                       ! row 2 = wed
                                       ! row 3 = hage
                                       ! row 4 = female (needed for singles for divorced - whether male or female) 
        integer allstate(numfam,numper)  ! keep track of what years the state variables are
                                         ! complete (no missing values)
        integer indskip
	    integer hdead,wdead
        integer numk, hw,ww,formh,formw,survivh,survivf,div,mardur,kidnow,totk
        integer female,lhw,lww,lformh,lformw
        real(8) earnm,earnw
        real(8) hexper,wexper,pensionw,pensionh,wealth,thexper,twexper
        real(8) winformexp,hinformexp,tfhexper,tfwexper
        integer hdeath,wdeath,T0,Tbar
	real(8) hformexp,wformexp
        

	!set the variable that will count the number of lines
	!read in from the input file. Will read until end of file
	!marker
		i = 1

	!initialize variables
		lastid = -1
        	j = 0 ! indexes the households
		yearmat = 0
        	allstate = 0 
        	nivart(1,1:numfam,1:numper) = 0
        	nivart(2,1:numfam,1:numper) = 0
        	nivart(5,1:numfam,1:numper) = 0
        	nivart(6,1:numfam,1:numper) = 0
        	nivart(7,1:numfam,1:numper) = 0
        	nivart(8,1:numfam,1:numper) = 0
        	nivart(9,1:numfam,1:numper) = 0
  
        open(unit=100,file='./inputs/Initialconditions_new.txt',position='rewind')

	do while(i.le.maxl) 

	    	read(fmt=*,unit=100,end=500) id,year,numk,kidnow,married,hw,ww,hexper,wexper, &
          		hformexp,hinformexp,wformexp,winformexp,hdeath,wdeath,hdead,wdead,div,mardur, &
          		formh,formw,hed,wed,hage,wage,female, &
          		wealth,pensionh,pensionw,earnm,earnw,lformh,lformw,lhw,lww,singlesamp

		!deal with data anomalies
		if (id.eq.77025) then
			female = 0
		endif        	
		!recode missing or topcoded variables
		if (wage.eq.-99) then 
           		wage = hage  !deal with case of a single male or married couple with missing wife age
        	end if
            	if (wage.eq.-99.or.(female.eq.0.and.singlesamp.eq.1.and.hage.eq.-99)) then
              		write(*,*) 'problem, missing age***'
            		pause
		end if
        	if (married.eq.1.and.wage.gt.hage) hage = wage 	! hard to deal with women who are older than men (but rare)
        	if (numk.gt.kidmax) then
         		numk=kidmax
        	endif

        	if (wage.lt.Tbar.and.wage.ge.T0) then 
            		if (hed.eq.-99.and.wed.eq.-99) then
                		write(*,*) 'id, ed missing for both',id
                		pause
            		end if 
            		if (hed.eq.-99) then
                		hed = 0
            		end if
            		if (wed.eq.-99) then
                		wed = 0
            		end if
            		if (pensionw.eq.-99.0d-0) then
                		pensionw = 0.0d-0
            		end if
            		if (pensionh.eq.-99.0d-0) then
                		pensionh = 0.0d-0
            		end if
            		if (female.eq.2) then 
                		female = 1
            		else
                		female = 0
            		end if
        		if (married.eq.0) then
            			if (hw.eq.-99) then
               				hw = 0
            			end if
            			if (ww.eq.-99) then
               				ww = 0
            			end if
        		end if
            		if (hw.eq.0) then
                		earnm = 0.0d-0
            		end if
            		if (ww.eq.0) then
                		earnw = 0.0d-0
            		end if
            		if (kidnow.eq.-99) then
                		write(*,*) 'kidnow:',kidnow
            		end if
	
            
	    		!ct = model period : corresponds to woman's age, or if single male, to man's age
			if (female.eq.1.or.singlesamp.eq.0) then
	            		ct = wage - 15 
	    		elseif (female.eq.0.and.singlesamp.eq.1) then
		    		ct = hage - 15
	    		endif

			!update fixed household variables if a new id is read
        		if(lastid.ne.id) then
              			j=j+1
              			lastid = id
              			yearmat(j,1)=year ! first year in which the household is observed
              			yearmat(j,2)=ct	! first model period in which the household is observed
				nfixvar(1,j) = hed
            			nfixvar(2,j) = wed
            			nfixvar(3,j) = hage - ct + 1 ! husband's age when counter is 1 (ie. the "household was 16" )
            			if (hage.eq.-99) then
                			nfixvar(3,j) = -99
            			end if
            			nfixvar(4,j) = female
            			idvec(j) = id
            		end if


			!determine if all state variables are non missing for that period
           		allstate(j,ct) = 1 
            		if (married.eq.1.and.(hexper.eq.-99.or.wexper.eq.-99.or.wage.eq.-99.or.hage.eq.-99.or.lww.eq.-99.or.lhw.eq.-99)) then
                		allstate(j,ct) = 0 ! married couples
            		end if
            		if (married.eq.0.and.female.eq.1.and.(wexper.eq.-99.or.wage.eq.-99.or.lww.eq.-99)) then
                		allstate(j,ct) = 0 ! single female
            		end if
            		if (married.eq.0.and.female.eq.0.and.(hexper.eq.-99.or.hage.eq.-99.or.lhw.eq.-99)) then
                		allstate(j,ct) = 0 ! single male
            		end if
            		if (wealth.eq.-99.0d-0) then
                		allstate(j,ct) = 0
            		end if       
			
			!other reasons for not simulating this household:
            		if (wage.ge.Tbar) then
                		allstate(j,ct) = 1
            		end if
            		if (married.eq.0.and.female.eq.0.and.hage.ge.Tbar) then
                		allstate(j,ct) = 0 
            		end if
		
			!update variable household characteristics	
        		nivart(3,j,ct) = hw 
        		nivart(4,j,ct) = ww 
        		if (ct.gt.1) then      ! lagged work sector
          			nivart(3,j,(ct-1)) = kidnow   ! if had kid this period, then preg last period
          			nivart(3,j,(ct-1)) = lhw 
          			nivart(4,j,(ct-1)) = lww 
        		end if
        		nivart(2,j,ct) = married 
        		nivart(5,j,ct) = 1-hdeath                    
        		nivart(6,j,ct) = 1-wdeath 

        		if (hdeath.eq.1.and.wdeath.eq.1) then 
            			write(*,*) 'both husband and wife died'
            			write(*,*) 'folio:',id
            			allstate(j,ct) = 0 
            			pause 
        		end if
        		if (year.eq.2004.and.ct.lt.numper) then   !set survival and marital status to be same in 2005 as in 2004
          			nivart(2,j,(ct+1)) =  nivart(2,j,ct)                    
          			nivart(7,j,(ct+1)) =  nivart(7,j,ct)+nivart(2,j,ct)  !marriage duration                    
         	 		nivart(5,j,(ct+1)) = 1-hdeath                    
          			nivart(6,j,(ct+1)) = 1-wdeath 
        		end if


 
        		!XXXif did not die this period, then need to assign some time
        		!whether woman or man died previous period
        		!for now I am treating them as being not married

        		if (hdeath.eq.0.and.hdead.eq.1) then
            			married = 0
        		end if
        		if (wdeath.eq.0.and.wdead.eq.1) then
            			married = 0
        		end if
        		if (ct.gt.1.and.hdeath.eq.0.and.hdead.eq.1) then
           			nivart(5,j,(ct-1)) = 1-hdead             
        		end if      
        		if (ct.gt.1.and.wdeath.eq.0.and.wdead.eq.1) then
           			nivart(6,j,(ct-1)) = 1-wdead             
        		end if      
                 

        		! fill in some values in work so that experience 
        		! will add to the right number
        		thexper = hexper
        		twexper = wexper
        		tfhexper = hformexp
        		tfwexper = wformexp
        		if (ct.gt.1) then
        		exploop: do k = (ct-1),1,-1
          			if (thexper.gt.0) then
             				nivart(3,j,k) =  1
             				thexper = thexper-1.0d-0
          			else
             				nivart(3,j,k) = 0 
          			end if
          			if (tfhexper.gt.0) then
             				nivart(8,j,k) =  1
             				tfhexper = tfhexper-1.0d-0
          			else
             				nivart(8,j,k) = 0 
          			end if
          			if (twexper.gt.0) then
             				nivart(4,j,k) =  1
             				twexper = twexper-1.0d-0
          			else
             				nivart(4,j,k) = 0 
          			end if
          			if (tfwexper.gt.0) then
             				nivart(9,j,k) =  1
             				tfwexper = tfwexper-1.0d-0
          			else
             				nivart(9,j,k) = 0 
          			end if
        		end do exploop
        		end if
			
			! create fertility history
        		totk = 0
        		totk = numk-sum(nivart(1,j,1:ct))   ! number of pregnancies woman needs to still have recorded
                                            ! in nivart

        		kidloop: do k = min(ct-1,25),1,-1   ! starting at ct and going backwards
               				if (totk.gt.0.and.nivart(1,j,k).ne.1) then  
                     				nivart(1,j,k) = 1  ! record a pregnancy
                     				totk = totk-1      ! move indicator backwards
               				end if
        		end do kidloop

       
        		if (ct.gt.1) then
          			nrvart(1,j,(ct-1)) = wealth            ! put all wealth in safe asset,
          			nrvart(5,j,(ct-1)) = pensionw            
          			nrvart(6,j,(ct-1)) = pensionh            
        		end if

        		nrvart(3,j,1:numper) = 1.05d-0       ! XXX returns (ensure consistency with rest of code)
        		nrvart(4,j,1:numper) = 0.0d-0        ! returns


        		i = i+1

        	end if !wage between T0 and Tbar 

        end do  !end do while loop

500        close(unit=100)

        numvalid = j

        return
	end subroutine GetData


! this subroutine considers all the choices at a given point in time
! and solves for the maximum utility (for each given vector of shocks)

subroutine maxu(o,ct,clone,thisyear,umaxmc,nmc,fixvar,ivart,rvart,yevec,alpha,thetap,Tbar,numz,ktp,totp,mptype,dpmat,intmat,optch,optsafe,optrisky,earnh,earnw,kidmax,experh,experw,maxearnf,maxearnh,numa1cat,deathprob,cnu,optpenw,optpenm,probvec,maxsafe,mpenwithdraw,fpenwithdraw,reform,actret,nompg,mpg, nopasis, pasis,nobono, bono, nodivrule, pmax,pmin,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt, optc, tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual)
         
implicit none

! external u,TermDiv,TermHDec,TermWDec,ran1

integer, parameter :: ncoeff=150
integer, parameter :: numper=75
integer lastage,lastworkage,counterfactual 
integer mptype,intmat(75,75), numk, clone, thisyear,hageret,wageret,wageret2
integer nmc,ktp,totp,jj,tt,wlim,hlim,wllim,hllim,reform,nompg, nopasis, nobono, nodivrule
integer vv,t1,t2,ct,i,o,kk,temp1,temp2,optch(nmc,7),numch
real(8) yevec(nmc*7),yyevec(7),c,cvec(nmc),deathprob(91,3), mpg, pasis
real(8) umaxmc(nmc),thetap(3,numz),wsum,wsumopt,cnu(76,4),ubequest,utilm,utilf
integer fert,hw,ww,perc,Tbar,numz,ii,formh,formw,kidmax
real(8) experw,experh,safe,risky,ret1,ret2,maxsafe(numper)
real(8) util,maxutil,TermDiv,TermHDec,TermWDec,term4
real(8) alpha(ncoeff,ktp),dpmat(11,3,2),earnh(nmc),earnw(nmc),ttemp1,ttemp2
real(8) probcofferf,probcofferm,r1(1),r2(1),earnhf,earnhi,earnwf,earnwi
real(8) earnhfp,earnhip,earnwfp,earnwip,optpenw(nmc),optpenm(nmc), npensionh,npensionw
real(8) penw,penm,probvec(4),pmax,pmin,actret, nwealth,tao, fexperh, fexperw

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

 integer fertlim,deceased,married,terminal,dechus,decwife,hage,wage,a1,a2
 real(8) termval,nsafe,nrisky,optsafe(nmc),optrisky(nmc),optc(nmc)
 integer female
 real(8) maxearnf,maxearnh,fpenwithdraw,mpenwithdraw,tmpenwithdraw,tfpenwithdraw
 integer hed,wed,temp,first,idum,formlimf,formlimm,evermarried
 integer numa1cat
 real(8) tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost,term4lag
 real(8) tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt  , bono


idum = 11111 
wsum = 0.0d-0
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
evermarried = ivart(2,1,20)
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
  

fertlim = 0  !fertility not a choice in this model
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
    t2 = t1+4
    yyevec(1:7) = yevec(t1:t2)

    earnhf = earnhfp*dexp(yyevec(1))
    earnhi = earnhip*dexp(yyevec(2)) 
    earnwf = earnwfp*dexp(yyevec(3)) 
    earnwi = earnwip*dexp(yyevec(4))

    perc = 0.5d-0
    util = 0.0d-0
    ttemp1 = 0.0d-0
    ttemp2 = 0.0d-0
    maxutil = 0.0d-0
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

   
     call u(util,ct,o,yyevec,fixvar,ivart,rvart,fert,hw,ww,perc,Tbar,numz,thetap, alpha,ktp,totp,numper,mptype,hage,wage,experh,experw,fexperh,fexperw,dpmat,a1,a2,intmat, nsafe,nwealth, nrisky,ttemp1,ttemp2,c,wsum,formh,formw,married,female,earnhf,earnhi,earnwf,earnwi, maxearnf,maxearnh,evermarried,term4,numa1cat,nmc,deathprob,cnu,npensionw,npensionh,probvec,maxsafe, tmpenwithdraw,tfpenwithdraw,reform,actret, nompg, mpg, nopasis, pasis, nobono, bono,nodivrule, pmax,pmin,tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost,ubequest,utilm,utilf, numk,tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual) 

!if (a1.gt.1.and.term4.gt.term4lag) then
!print*, "nonmonotonic: period",ct,"state draw",o,"hw",hw,"ww",ww,"formh",formh,"formw",formw,"a1",a1
!continue
!endif

!XXX TESTING
!
!if (thisyear.ge.2008.and.thisyear.lt.2010) then
!write(9999,*) o,ct,clone, mptype, thisyear, wage, hage, ww,hw,formh,formw,a1, married,female,numk, reform, safe, tmpenwithdraw, tfpenwithdraw, tax_paid, nwealth, nsafe, c,ttemp1,ttemp2,npensionw,npensionh, util, utilf, utilm,term4,ubequest,optch(kk,4), optch(kk,1), optch(kk,2),optch(kk,6),optch(kk,7),yyevec(1:5),tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost
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
        wsumopt = wsum
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
!       is designated by the combination (fert,hw,ww,perc,formh,formw) 
!
!       In doing so, it updates the state space for each set of choices
!       being considered
!---------------------------------------------------------------------------------------------------


subroutine u(util,ct,o,yyevec,fixvar,ivart,rvart,fert,hw,ww,perc,Tbar,numz,thetap, &
        alpha,ktp,totp,numper,mptype,hage,wage,experh,experw,fexperh,fexperw,dpmat,a1,a2,intmat, &
        nsafe, nwealth, nrisky,earnh,earnw,c,wsum,formh,formw,married,female,earnhf,earnhi,earnwf,earnwi, &
        maxearnf,maxearnh,evermarried,term4,numa1cat,nmc,deathprob,cnu,npensionw,npensionh,probvec,maxsafe, &
        mpenwithdraw,fpenwithdraw,reform,actret, nompg,mpg, nopasis, pasis, nobono, bono, nodivrule, pmax,pmin,&
        tax_paid,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost,ubequest,utilm,utilf, numk, tao,hageret,wageret,wageret2,&
	lastage,lastworkage, counterfactual)


implicit none

! external ProbDiv,TermDiv,TermHDec,TermWDec,tax
integer, parameter :: ncoeff=150
integer, parameter :: nbracket=8

integer hageret,wageret,wageret2,intmat(75,75),female,evermarried,nmc,reform, nompg, nopasis, nobono, nodivrule
integer numa1cat
integer ktp,totp,numz,numper,hhs,hsomec,hcol,whs,wsomec,wcol,hage,wage
real(8) experh,experw,earnhf,earnhi,earnwf,earnwi,deathprob(91,3), bono,B(numz,1),A(numz)
real(8) thetap(3,numz),dpmat(11,3,2),probvec(4),pmax,pmin,actret, supplement
real(8) ProbDiv,TermDiv,TermHDec,TermWDec,cnu(76,4),maxsafe(numper),term4bis
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
real(8) curwealth,prevwealth
real(8) retsafe,maxearnf,maxearnh
real(8) retrisky,tao
integer hed,wed,fert,hw,ww,perc,mptype,lastage,lastworkage
integer n(numper),kk,jj
integer married,exph,expw
integer ct,o,ctt,numkid,schoolh,schoolw
real(8) yyevec(7),safe,risky,wsum1,wsum2,wsum,wsum1d,wsum1fd,wsum2d
real(8) wsum1f,wsum2f,wsumf
integer numsvar,nextct,a1,a2,counterfactual
real(8) Ztemp(1,numz),term4,term41,term42
real(8) util,alpha(ncoeff,ktp),femutil,malutil
integer Tbar,ii,numk
real(8) sum1,earnh,earnw,sum2,sum3,sum4
real(8) c, cf,cm,utilf,utilm,hcaph,hcapw
real(8) pdiv,phdec,pwdec,term1,term2,term3,cmin,mpg
real(8) beta,ret1,ret2,nwealth,nsafe,nrisky,temp,nsafe_inv,nrisky_inv
integer mardur,fulltime,ind2024,ind2529,ind3035,ind3639,ind4045,parttime
integer ind65h,ind60w,formh,formw,ind70h
real(8) probfert,temp22,pensionw,pensionh,npensionw,npensionh
integer eswitch,mswitch,fswitch,mreentry,freentry,fqualearly,mqualearly
real(8) earlycut(20)
real(8) ubequest, ubequest1, ubequest2, cbis, cdif,util2,ubequestbis

!! variables related to taxes
integer :: tax_toggle,tempw,tempm
real(8) bracket_min(nbracket)    !bracket lower bounds 
real(8) bracket_max(nbracket)   !bracket upper bounds
real(8) bracket_rate(nbracket)   !marginal tax rates
real(8) bracket_adj(nbracket)    !taxes due are obtained by multiplying
                                 !the taxable income with the marginal tax rate
                                 !and subtracting the backet-specific adjustment
real(8) taxable_wages,tax_paid,taxw,taxm,tax,pasis,fpenwithdraw,mpenwithdraw,faminc,incthreshold
real(8) pasisqualm,pasisqualf,mpgqualf,mpgqualm,spsm,spsf,fexperh,fexperw,utilfbis,utilbis,utilmbis
real(8) temp111m,temp111f,spsmcost,spsfcost,fmpgcost,mmpgcost,pasiscost, earnpot,consgrid(numa1cat)
integer group

! Initialize


data bracket_min(:) /0 ,4.909d-0,10.910d-0,18.184d-0,25.458d-0,32.732d-0,43.643d-0,54.554d-0 /&
bracket_max(:) /4.909d-0,10.910d-0,18.184d-0,25.458d-0,32.732d-0,43.643d-0,54.554d-0,9999999/ &
bracket_rate(:) /0,0.05d-0,0.1d-0,0.15d-0,0.25d-0,0.32d-0,0.37d-0,0.4d-0 / &
bracket_adj(:) /0,0.245d-0,0.791d-0,1.7d-0,4.246d-0,6.537d-0,8.719d-0,10.356d-0 /  

beta = alpha(15,1)
cmin = alpha(85,1)*1000000.0d-0    ! set a minimum consumption level
ret1 = 0.05d-0
ret2 = 0.0d-0
tax_toggle = 2
util = -9999999999.0d-0
ctt=ct+1
ind2024 = 0
ind2529 = 0
ind3035 = 0
ind3639 = 0
ind4045 = 0
fulltime = 0
parttime = 0
pensionw = 0.0d-0
pensionh = 0.0d-0
pasisqualf = 0
pasisqualm = 0
mpgqualf = 0
mpgqualm = 0
temp111f = 0.0d-0
temp111m = 0.0d-0
c = 0
cf = 0
cm = 0
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

if (ww.eq.1) then 
        fulltime = 1
else if (ww.eq.2) then
        parttime = 1
end if

wage = ct+15
mardur = ivart(7,1,ct)
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
  taxm = earnh   !do not get taxed on contributions
end if
if (hw.eq.1.and.formh.eq.0) then
  earnh = earnhi
end if
if (ww.eq.1.and.formw.eq.1) then
  earnw = (1-tao)*earnwf 
  pensionw = pensionw + tao*earnwf      ! 10% of formal earnings goes to pension
  taxw = earnw  !do not get taxed on contributions
end if
if (ww.eq.2.and.formw.eq.1) then
  earnw = (1-tao)*0.5d-0*earnwf   !earn half of a full-time worker
  pensionw = pensionw + tao*0.5d-0*earnwf      ! 10% of formal earnings goes to pension
  taxw = earnw   !do not get taxes on contributions
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
temp111m = 0.0d-0
temp111f = 0.0d-0

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
    pdiv = ProbDiv(ct,ncoeff,alpha,hage,wage,numk,ktp,mardur,hed,wed)
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

call utilfunction2(alpha,ncoeff,ktp,ct,mptype,c,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafe,npensionh,npensionw,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)

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

  
   call contvalue(married,female,numper,term4,ubequest,thetap,reform,nodivrule,nivart,nrvart,nfixvar,ct,o,phdec,pwdec,pdiv,probfert,nsafe,npensionh,npensionw,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,numz,totp,intmat,pmax,pmin,pasis,bono,cmin,cnu,lastage,lastworkage)
   
   util = util+beta*term4
   
end if !ct+15 

! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
!    output choice utilities


!if (reform.eq.0.and.o==2045.and.(ct==74)) then
!open(unit=590,position='append',file='emaxhat.txt')
!write(590,*) reform,o, ct, hw,formh,ww,formw,a1,util,term4,wsum1,ubequest,married,female,hed,wed,numk,experw,experh,fexperw,fexperh,safe,pensionw,pensionh
!close(590)
!endif

if ((o.eq.781.or.o.eq.960).and.ct==25) then
continue
endif

100        end subroutine u 


























function ProbDiv(ct,ncoeff,alpha,hage,wage,numk,ktp,mardur,hed,wed)
implicit none

integer ncoeff,hage,wage,ktp,numk
real(8) alpha(ncoeff,ktp),temp,xb
real(8) ProbDiv,pnodiv
integer mardur,ct,wed,hed,agediff2,wage2,agediff

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


subroutine utilfunction2(alpha,ncoeff,ktp,ct,mptype,cbis,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafebis,npensionhbis,npensionwbis,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)

integer:: ncoeff,ktp,ind65h,ind60w,hage,Tbar,ct,mptype,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,parttime,fulltime
real(8):: nsafe,ubequest,ubequest1,ubequest2,alpha(ncoeff,ktp),npensionh,npensionw,beta,phdec,pwdec,yyevec(7),util,c,utilf,utilm,cf,cm
real(8):: MUCshiftf,MUCshiftm,MULf,MULm,switchf,switchm,nonpecf,nonpecm,cbis,nsafebis,npensionhbis,npensionwbis

c=0.0d-0
nsafe=0.0d-0
npensionh=0.0d-0
npensionw=0.0d-0
c=cbis/1000000.0d-0 !express consumption in million pesos
nsafe=nsafebis/1000000.0d-0
npensionh=npensionhbis/1000000.0d-0
npensionw=npensionwbis/1000000.0d-0

parttime=0
fulltime=0
if (ww.eq.1) then 
fulltime = 1
else if (ww.eq.2) then
parttime = 1
end if

util = 0.0d-0
utilm = 0.0d-0
utilf = 0.0d-0

cf = c/(1+married)
cm = c/(1+married)       

MUCshiftf=dexp(alpha(2,1)*dble(numk)+ alpha(4,1)*dble((1-parttime)*(1-fulltime)))
MULf=alpha(6,mptype)*dexp(alpha(10,1)*dble(numk) +alpha(11,1)*married + alpha(14,1)*dble(ind60w) + yyevec(7))
switchf=dble(fswitch)*alpha(63,1)+dble(freentry)*alpha(65,1)
nonpecf=dble(1-formw)*dble(parttime+fulltime)*alpha(70,mptype)

MUCshiftm=dexp(alpha(3,1)*dble(numk)+alpha(5,1)*dble(1-hw))
MULm=alpha(7,mptype)*dexp(alpha(13,1)*dble(ind65h) + yyevec(6))
switchm=dble(mswitch)*alpha(62,1)+dble(mreentry)*alpha(64,1)
nonpecm=dble(1-formh)*dble(hw)*alpha(71,mptype)


utilf = ((cf**alpha(1,1))/alpha(1,1)) * MUCshiftf &
+ (dble((1-parttime)*(1-fulltime))+alpha(8,1)*dble(parttime))* MULf &
- switchf + nonpecf   

utilm = ((cm**alpha(1,1))/alpha(1,1)) * MUCshiftm &
+ dble(1-hw)* MULm &
- switchm + nonpecm 

util = utilf*dble(married+(1-married)*female)+utilm*dble(married+(1-married)*(1-female))

!bequest

ubequest1=alpha(81,1)*dble(maxval((/numk,1/))) !number of periods and kids over which the bequest will be split
ubequest2=(nsafe+npensionh+npensionw)/(ubequest1+0.00001)+alpha(83,1) !split bequest + baseline income of heirs
ubequest=ubequest1*((ubequest2+0.00001)**alpha(1,1))/alpha(1,1) !utility of heirs'consumption

!add bequest if last period
if ((ct+15).eq.Tbar.or.(hage.eq.90.and.married.eq.0.and.female.eq.0)) util= util + beta*ubequest


endsubroutine

subroutine contvalue(married,female,numper,term4,ubequest,thetap,reform,nodivrule,nivart,nrvart,nfixvar,ct,o,phdec,pwdec,pdiv,probfert,nsafe,npensionh,npensionw,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,numz,totp,intmat,pmax,pmin,pasis,bono,cmin,cnu,lastage,lastworkage)

!!!INCORPORATE LASTAGE


integer:: married,female,group,ctt,numz,reform,nodivrule,ct,o,numper,numk,totp,intmat(75,75),nivart(9,1,numper),nfixvar(4,1),lastage,lastworkage
real(8):: ubequest,wsum,wsum1,wsum1d,wsum2,wsum2d,wsumf,wsum1f,wsum1fd,wsum2f,term4,term42,term41,phdec,pwdec,pdiv,bono,probfert
real(8):: Ztemp(1,numz), thetap(3,numz),nrvart(6,1,numper)
real(8):: experh,experw,fexperh,fexperw,maxearnh,maxearnf,pmax,pmin,pasis,cmin,cnu(76,4),npensionh,npensionw,nsafe

!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!        Compute EMAX terms  if NOT last period

 
    wsum = 0.0d-0
    wsum1 = 0.0d-0
    wsum1d = 0.0d-0
    wsum2 = 0.0d-0
    wsum2d = 0.0d-0
    wsumf = 0.0d-0
    wsum1f = 0.0d-0
    wsum1fd = 0.0d-0
    wsum2f = 0.0d-0
    
       ! CASE 1: couple's continuation value: no fertility => wsum
   group=1 !married couple
   ctt=ct+1
   nfixvar(4,1) = female                        ! household gender stays same
   nivart(1,1,ctt) = 0                          ! fert is 0
   nivart(2,1,ctt) = 1                          ! stay married, without fertility
   nivart(7,1,ctt) = nivart(7,1,ct)+1           ! mardur will be one year longer if stay married 
   nrvart(1,1,ct) = nsafe          
   nrvart(5,1,ct) = npensionw          
   nrvart(6,1,ct) = npensionh
   call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
   wsum = dot_product(Ztemp(1,:),thetap(group,:))
!   if(ct.eq.74.and.o.eq.1) then
!   print*,wsum
!   print*,Ztemp(1,:)
!   print*,thetap(group,:)
!   continue
!   endif


       ! CASE 2: single or divorced female continuation value: no fertility => wsum1
   group=2 !single female
   nfixvar(4,1) = 1                             ! household is female
   nivart(1,1,ctt) = 0                          ! fert is 0
   nivart(2,1,ctt) = 0                          ! become an unmarried female (divorce), without fertility
   nivart(7,1,ctt) = nivart(7,1,ct)             ! mardur will stay the same if get divorced/widowed        
   nrvart(1,1,ct)=0.5d-0*nsafe                  ! get half of private savings
   if (reform.eq.1.and.nodivrule.eq.0) then
     nrvart(5,1,ct)= (npensionh+npensionw)/2.0d-0    
     nrvart(6,1,ct)= (npensionh+npensionw)/2.0d-0    
   else
     nrvart(5,1,ct) = npensionw
     nrvart(6,1,ct) = npensionh 
   end if
   call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
   wsum1 = dot_product(Ztemp(1,:),thetap(group,:))
   

       ! CASE 3: wife's continuation value after husband dies => wsum1d
   group=2 !single female
   nfixvar(4,1) = 1                             ! female
   nivart(1,1,ctt) = 0                          ! fert is 0
   nivart(2,1,ctt) = 0                          ! become  unmarried (death of husband)
   nivart(7,1,ctt) = nivart(7,1,ct)             ! mardur will stay the same if get divorced/widowed   
   nrvart(1,1,ct) = nsafe                       ! keeps all the private savings          
   nrvart(5,1,ct)= npensionw+0.6d-0*npensionh   ! survivorship pension put in wifes account 
   nrvart(6,1,ct)= 0.0d-0     
   call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
   wsum1d = dot_product(Ztemp(1,:),thetap(group,:))


       ! CASE 4: single or divorced male's continuation value after divorce => wsum2
    group=3                                     ! single  male
    nfixvar(4,1) = 0                            ! male
    nivart(1,1,ctt) = 0                         ! fert is 0 
    nivart(2,1,ctt) = 0                         ! become an unmarried male (divorce), without fertility 
    nivart(7,1,ctt) = nivart(7,1,ct)            ! mardur will stay the same if get divorced/widowed        
    nrvart(1,1,ct) = 0.5d-0*nsafe               ! get half of private savings
    if (reform.eq.1.and.nodivrule.eq.0) then
     nrvart(5,1,ct)= (npensionh+npensionw)/2.0d-0    
     nrvart(6,1,ct)= (npensionh+npensionw)/2.0d-0    
   else
     nrvart(5,1,ct) = npensionw
     nrvart(6,1,ct) = npensionh
    end if
    call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
    wsum2 = dot_product(Ztemp(1,:),thetap(group,:))


       ! CASE 5: husband's continuation value after wife dies => wsum2d
    group=3         
    nivart(1,1,ctt) = 0                         ! fert is 0 
    nivart(2,1,ctt) = 0                         ! become an unmarried male (widowed), without fertility 
    nivart(7,1,ctt) = nivart(7,1,ct) 
    nfixvar(4,1) = 0
    nrvart(1,1,ct) = nsafe                      ! get all of private savings
    nrvart(5,1,ct)=0.0d-0                   
    nrvart(6,1,ct)= nrvart(6,1,ct)              !loses wifes pension 
    call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
    wsum2d = dot_product(Ztemp(1,:),thetap(group,:))

   ! fertility is an option if household includes  female< age 40
   if ((married.eq.1.or.female.eq.1).and.ct.le.25) then 
   
           ! CASE 6: couple's continuation value, with fertility => wsumf
       group=1 !married couple
       nfixvar(4,1)=female
       nivart(1,1,ctt) = 1                      ! fert - 1 
       nivart(2,1,ctt) = 1                      ! stay married, with fertility
       nivart(7,1,ctt) = nivart(7,1,ct)+1 
       nrvart(1,1,ct) = nsafe                      ! get all of private savings
       nrvart(5,1,ct)= npensionw   + bono  
       nrvart(6,1,ct)= npensionh
       call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk+1,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
       wsumf = dot_product(Ztemp(1,:),thetap(group,:))


            ! CASE 7: single or divorced woman continuation value, with fertility
        group=2                             !single female
        nfixvar(4,1) = 1
        nivart(1,1,ctt) = 1 
        nivart(2,1,ctt) = 0                 ! become an unmarried female(divorce), with fertility
        nivart(7,1,ctt) = nivart(7,1,ct) 
        nrvart(1,1,ct)=0.5d-0*nsafe
        if (reform.eq.1.and.nodivrule.eq.0) then
        nrvart(5,1,ct)= dmax1((npensionh+npensionw)/2.0d-0,npensionw)+bono 
        nrvart(6,1,ct) = 0.0d-0
        else
        nrvart(5,1,ct)= npensionw + bono
        nrvart(6,1,ct) = 0.0d-0
        end if
   	call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk+1,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
        wsum1f = dot_product(Ztemp(1,:),thetap(group,:))


            ! CASE 8: widow's continuation value, with fertility
        group=2 !single female
        nfixvar(4,1) = 1
        nivart(1,1,ctt) = 1 
        nivart(2,1,ctt) = 0          ! become an unmarried female(widowed), with fertility
        nivart(7,1,ctt) = nivart(7,1,ct) 
        nrvart(1,1,ct)= nsafe
        nrvart(5,1,ct)= npensionw+0.6d-0*npensionh  +bono  !get husbands pension put in wifes account 
        nrvart(6,1,ct)= 0.0d-0     
   	call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk+1,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
        wsum1fd = dot_product(Ztemp(1,:),thetap(group,:))

    end if  !married.eq.1.or... 


! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
!    add future values, taking account that there is some probability of
!    becoming divorced or deceased next period 


   ! if married, cont value does not depend on whether female sampled or not
   if (married.eq.1.and.ct.gt.25) then
   term4 = (1.0d-0 - phdec)*(1.0d-0 - pwdec)*(1.0d-0-pdiv)*wsum + &   !neither dies and no divorce
            phdec*wsum1d +   &                            !husband dies, with or without divorce
            (1.0d-0-phdec)*(1.0d-0 - pwdec)*pdiv*(0.5d-0*wsum1+0.5d-0*wsum2) + &     !neither dies, but get divorced
            pwdec*wsum2d          + &                         !wife dies, with or without divorce 
            pwdec*phdec*ubequest                            !both die, get utility from bequest
   end if         
   
   if (married.eq.1.and.ct.le.25) then
   term41 = (1.0d-0 - phdec)*(1.0d-0 - pwdec)*(1.0d-0-pdiv)*wsum + &   !neither dies and no divorce
            phdec*wsum1d +   &                             !husband dies, with or without divorce
            (1.0d-0-phdec)*(1.0d-0-pwdec)*pdiv*(0.5*wsum1+0.5*wsum2) + &      !neither dies, but get divorced 
            pwdec*wsum2d  +&
            pwdec*phdec*ubequest                            !both die, get utility from bequest
   term42 = (1.0d-0 - phdec)*(1.0d-0 - pwdec)*(1.0d-0-pdiv)*wsumf + &   !neither dies and no divorce
            phdec*wsum1fd +   &                              !husband dies, with or without divorce
            (1.0d-0-phdec)*(1.0d-0-pwdec)*pdiv*(0.5d-0*wsum1f+0.5d-0*wsum2) + &      !neither dies, but get divorced 
            pwdec*wsum2d    +&                               !assume no child if wife dies
               pwdec*phdec*ubequest                            !both die, get utility from bequest
   term4 = (1.0d-0-probfert)*term41+probfert*term42
   end if         

   if (married.eq.0.and.female.eq.1.and.ct.gt.25) then      !single or divorced female age gt 40
   term4 = (1.0d-0-pwdec)*wsum1  + pwdec * ubequest
   end if         
   if (married.eq.0.and.female.eq.1.and.ct.le.25) then      !single or divorced female age le 40
   term4 = (1.0d-0-probfert)*(1.0d-0-pwdec)*wsum1+probfert*(1.0d-0-pwdec)*wsum1f + pwdec* ubequest
   end if   
         
   if (married.eq.0.and.female.eq.0) then               !single or divorced male
   term4 = (1.0d-0-phdec)*wsum2 + phdec *ubequest
   end if         

endsubroutine contvalue


    ! -----------------------------------------------------------------------------------------------------
    !       This subroutines constructs the state variables (at beginning of "ctt" period) from the information contained in
    !       fixvar,ivart,rvart - it returns the Z vector 
    !       
    SUBROUTINE getstate2(o,ctt,fixvar,ivart,rvart,experh,experw,numk,fexperh,fexperw,hcaph,hcapw,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage)


    implicit none

    integer, parameter:: nseg=7,nseg2=6,nseg3=5
    integer totp,numper,o,ct,wage,hage,k,j,kk,singlemale,singlefemale,i1,i2, lastage,lastworkage
    real(8) rvart(6,1,numper)      ! row 1 = wealth - safe
                                   ! row 2 = wealth - risky
                                   ! row 3 = rate of return safe
                                   ! row 4 = rate of return risky
                                   ! row 5 = pension wealth women
                                   ! row 6 = pension wealth men
    integer ivart(9,1,numper)      ! row 1 = children (n)
                                   ! row 2 = married
                                   ! row 3 = work h
                                   ! row 4 = work w
                                   ! row 5 = survivh
                                   ! row 6 = survivw
                                   ! row 7 = marriage duration
                                   ! row 8 = formal h 
                                   ! row 9 = formal w
    integer fixvar(4,1)            ! row 1 = hed
                                   ! row 2 = wed
                                   ! row 3 = hage
                                   ! row 4 = female 


    integer numz,n(numper),i,full,part,reform,postref
    integer hed,wed,intmat(75,75),ind,married,ctt
    real(8) Ztemp(1,numz),experw,experh,wealth,indw1,indw2,c1,c2,hcaph,hcapw,cmin,fexperh,fexperw
    integer tt,female,changestat,ind65h,ind70h,evermarried
    integer dech,decw,pos1,temp1,temp2,nowealth,mardur,numk,ind80h,ind90h
    real(8) ret1, ret2,pensionw,pensionh,pmin,pmax,cnu(76,4)
    real(8) nsafe, nrisky,hw,ww,experw2,experh2,fpenwithdraw,mpenwithdraw,withdrawh(nseg2),withdraww(nseg2),withdsegh(nseg2),withdsegw(nseg2),withdh,withdw
    real(8) maxearnf,maxearnh,pasis,mpg,temp111f,spsf,temp111m,spsm,basecons,seg(nseg),apwseg(nseg),logapwreg(nseg),apwreg(nseg),crraapwreg(nseg),apw,logapw,crraapw
    integer formalh,formalw,informalh,informalw,male
    integer pasisqual,consperh,workperh,consperw,workperw
    real pseg(nseg2),mpgqualw,mpgqualh,xpseg(nseg3),hxpseg(nseg3),wxpseg(nseg3)


 !    real(16),external:: logic2dbl
 
     INTERFACE
         elemental pure double precision function logic2dbl(a)
             logical, intent(in) :: a
         end function logic2dbl
     END INTERFACE


    informalh=0
    informalw=0
    ct = ctt-1
    wealth = rvart(1,1,ct) !use counter ct because it is end of period wealth
    pensionw = rvart(5,1,ct)
    pensionh = rvart(6,1,ct)
    married = ivart(2,1,ctt)
    evermarried = ivart(2,1,20)   !if married at age 35, then evermarried equals 1
    hed = fixvar(1,1)
    wed = fixvar(2,1)
    female = fixvar(4,1)
    hage = max0(35,min0((fixvar(3,1)-1+ctt),90))  
    wage = min0(ctt+15,90)
    hw = ivart(3,1,ct)                 !employment status
    ww = ivart(4,1,ct)                 !employment status
    formalw = ivart(9,1,ct)   ! indicator for female worked in formal sector
    formalh = ivart(8,1,ct)   ! indicator for male worked in formal sector
    if (hw.eq.1.and.formalh.eq.0) then
         informalh = 1
    end if 
    if (ww.ge.1.and.formalw.eq.0) then
         informalw = 1
    end if 
    mpg = 2.0d-0*pasis
    male = 1-female
    singlemale=male*(1-married)
    singlefemale=female*(1-married)
    mpgqualh = dmax1(20.0d-0-fexperh,0.0d-0)
    mpgqualw = dmax1(20.0d-0-fexperw,0.0d-0)

! income effect regressors (based on "annualized permanent wealth")

    !consumption periods left
    consperh=max0(1,lastage-hage) 
    consperw=max0(1,lastage-wage) 

    !work periods left
    workperh=max0(0,lastworkage-hage)
    workperw=max0(0,lastworkage-wage)

    !wealth segments
    basecons=2000000.0d-0
    seg(1)=0.0
    seg(2)=0.5*basecons
    seg(3)=basecons
    seg(4)=1.5*basecons
    seg(5)=3*basecons
    seg(6)=5*basecons
    seg(7)=-9999999999999999999.0d-0

    apw=(wealth+pensionh+pensionw+workperh*hcaph*(1-singlefemale)+workperw*hcapw*(1-singlemale))/(consperh+consperw) !current wealth plus earnings from remaining work periods divided by remaining years

    do i=1,nseg-1
        apwseg(i)=logic2dbl(apw.ge.seg(i))*(-1.0d-0)
        apwreg(i)=apwseg(i)*(apw-seg(i))
    enddo


!implicit marginal tax regressors (based on pension wealth accumulations interacted with apw)

    pseg(1)=0.0
    pseg(2)=basecons/4
    pseg(3)=basecons/2
    pseg(4)=basecons
    pseg(5)=basecons*2
    pseg(6)=999999999999.0

    call piksrt(6,pseg)

    withdh=0.0d-0
    withdw=0.0d-0

    if (married.eq.1) then
    withdh=pensionh/dmax1(1.0d-0,(CNU(hage-34,4)-workperh*1.0d-0))
    withdw=pensionw/dmax1(1.0d-0,(CNU(wage-34,3)-workperw*1.0d-0))
    elseif (married.eq.0) then
    withdh=pensionh/dmax1(1.0d-0,(CNU(hage-34,2)-workperh*1.0d-0))
    withdw=pensionw/dmax1(1.0d-0,(CNU(wage-34,3)-workperw*1.0d-0))
    endif
    do i=1,nseg2-1 
    withdsegh(i)=logic2dbl(withdh>=pseg(i))*(-1.0d-0)
    withdrawh(i)=withdsegh(i)*(withdh-pseg(i))*apw
    withdsegw(i)=logic2dbl(withdw>=pseg(i))*(-1.0d-0)
    withdraww(i)=withdsegw(i)*(withdw-pseg(i))*apw
    enddo

!impact of years of contribution restrictions (based on number of years in formal sector)

    xpseg(1)=0.0
    xpseg(2)=10.0
    xpseg(3)=20.0
    xpseg(4)=30.0
    xpseg(5)=9999999999999.0

    do i=1,nseg3-1
    hxpseg(i)=logic2dbl(fexperh.ge.xpseg(i+1))*(-1.0d-0)
    wxpseg(i)=logic2dbl(fexperw.ge.xpseg(i+1))*(-1.0d-0)
    enddo
    
! fill the Z matrix (some vectors may be zero)

!linear terms
    Ztemp(1,1) = 1.0d-0
    Ztemp(1,2) = female*1.0d-0  !collinear
    Ztemp(1,3) = married*1.0d-0  !collinear
    Ztemp(1,4) = logic2dbl(numk==0)*1.0d-0  !collinear
    Ztemp(1,5) = max0(numk,1)*1.0d-0  !collinear
    Ztemp(1,6) = hed*1.0d-0 *(1-singlefemale)
    Ztemp(1,7) = wed*1.0d-0 *(1-singlemale)
    Ztemp(1,8) = experh*1.0d-0 *(1-singlefemale)
    Ztemp(1,9) = experw*1.0d-0 *(1-singlemale)
    Ztemp(1,10) = hcaph*1.0d-0 *(1-singlefemale)*logic2dbl(hage.le.lastworkage)
    Ztemp(1,11) = hcapw*1.0d-0 *(1-singlemale)*logic2dbl(wage.le.lastworkage)
    Ztemp(1,12) = formalh*1.0d-0 *(1-singlefemale)*logic2dbl(hage.lt.lastworkage)
    Ztemp(1,13) = formalw*1.0d-0 *(1-singlemale)*logic2dbl(wage.lt.lastworkage)
    Ztemp(1,14) = informalh*1.0d-0 *(1-singlefemale)*logic2dbl(hage.lt.lastworkage)
    Ztemp(1,15) = informalw*1.0d-0*(1-singlemale)*logic2dbl(wage.lt.lastworkage)
    Ztemp(1,16) = wealth*1.0d-0 !collinear
    Ztemp(1,17) = pensionh*1.0d-0 *(1-singlefemale)
    Ztemp(1,18) = pensionw*1.0d-0*(1-singlemale)
    Ztemp(1,19) = fexperh*1.0d-0*(1-singlefemale)
    Ztemp(1,20) = fexperw*1.0d-0*(1-singlemale)
    Ztemp(1,21) = logic2dbl(numk==1)*1.0d-0
    Ztemp(1,22) = logic2dbl(numk==2)*1.0d-0
    Ztemp(1,23) = logic2dbl(numk==3)*1.0d-0
    Ztemp(1,24) = logic2dbl(numk==4)*1.0d-0 !collinear
    Ztemp(1,25) = 0.0d-0
    Ztemp(1,26) = 0.0d-0
    Ztemp(1,27) = 0.0d-0
    Ztemp(1,28) = 0.0d-0
    Ztemp(1,29) = 0.0d-0
            
! income effect regressors 

    Ztemp(1,30) = apwreg(1)*1.0d-0 
    Ztemp(1,31) = apwreg(2)*1.0d-0 
    Ztemp(1,32) = apwreg(3)*1.0d-0 
    Ztemp(1,33) = apwreg(4)*1.0d-0 
    Ztemp(1,34) = apwreg(5)*1.0d-0 
    Ztemp(1,35) = apwseg(1)*1.0d-0  
    Ztemp(1,36) = apwseg(2)*1.0d-0 
    Ztemp(1,37) = apwseg(3)*1.0d-0 
    Ztemp(1,38) = apwseg(4)*1.0d-0 
            
! Age differences

    Ztemp(1,39) = (hage-wage) * married * logic2dbl(ctt.lt.72)* (-1.0d-0)
    Ztemp(1,40) = logic2dbl(hage==wage+1) * married * (-1.0d-0)        
    Ztemp(1,41) = logic2dbl(hage==wage+2) * married * (-1.0d-0)      
    Ztemp(1,42) = logic2dbl(hage==wage+3) * married * (-1.0d-0)        
    Ztemp(1,43) = logic2dbl(hage==wage+4) * married * (-1.0d-0)        
    Ztemp(1,44) = logic2dbl(hage==wage+5) * married * (-1.0d-0)         
    Ztemp(1,45) = logic2dbl(hage==wage+6) * married * (-1.0d-0)         
    Ztemp(1,46) = logic2dbl(hage==wage+7) * married * (-1.0d-0)        
    Ztemp(1,47) = logic2dbl(hage==wage+8) * married * (-1.0d-0)       
    Ztemp(1,48) = logic2dbl(hage==wage+9) * married * (-1.0d-0)       
    Ztemp(1,49) = logic2dbl(hage==wage+10) * married * (-1.0d-0)        

! implicit marginal tax regressors

    Ztemp(1,50) = withdrawh(1)*1.0d-0 *(1-singlefemale)
    Ztemp(1,51) = withdrawh(2)*1.0d-0 *(1-singlefemale)
    Ztemp(1,52) = withdrawh(3)*1.0d-0 *(1-singlefemale)
    Ztemp(1,53) = withdrawh(4)*1.0d-0 *(1-singlefemale)
    Ztemp(1,54) = withdrawh(5)*1.0d-0 *(1-singlefemale)
    Ztemp(1,55) = withdraww(1)*1.0d-0 *(1-singlemale)
    Ztemp(1,56) = withdraww(2)*1.0d-0 *(1-singlemale)
    Ztemp(1,57) = withdraww(3)*1.0d-0 *(1-singlemale)
    Ztemp(1,58) = withdraww(4)*1.0d-0 *(1-singlemale)
    Ztemp(1,59) = withdraww(5)*1.0d-0 *(1-singlemale)    
    Ztemp(1,60) = 0.0d-0
    Ztemp(1,61) = 0.0d-0
    Ztemp(1,62) = 0.0d-0
    Ztemp(1,63) = 0.0d-0
    Ztemp(1,64) = 0.0d-0
    
! Formal experience regressors

    Ztemp(1,65) = hxpseg(1)*1.0d-0*(1-singlefemale)*logic2dbl(ctt.gt.20)
    Ztemp(1,66) = hxpseg(2)*1.0d-0*(1-singlefemale)*logic2dbl(ctt.gt.30)
    Ztemp(1,67) = hxpseg(3)*1.0d-0*(1-singlefemale)*logic2dbl(ctt.gt.40)
    Ztemp(1,68) = 0.0d-0
    Ztemp(1,69) = 0.0d-0
    Ztemp(1,70) = wxpseg(1)*1.0d-0*(1-singlemale)*logic2dbl(ctt.gt.20)
    Ztemp(1,71) = wxpseg(2)*1.0d-0*(1-singlemale)*logic2dbl(ctt.gt.30)
    Ztemp(1,72) = wxpseg(3)*1.0d-0*(1-singlemale)*logic2dbl(ctt.gt.40)
    Ztemp(1,73) = 0.0d-0
    Ztemp(1,74) = 0.0d-0
    Ztemp(1,75) = 0.0d-0


    ind = 75

! make the interacted variables 
    loop2: do i = 1,50
           if (intmat(i,1)==-9) exit loop2
           i1=intmat(i,1)
           i2=intmat(i,2)
           Ztemp(1,ind+i)=Ztemp(1,i1)*Ztemp(1,i2)
    end do loop2

    return
    end




      !------------------------- REPORTRES ---------------------------------

	SUBROUTINE ReportRes(thetam,sevm,rsqv,rsqcvv, &
        mptype,Sbar,Tbar,numz,numint,totp,nedraw,intmat,reform,pathoutput,ktp)

	
	integer, parameter :: numper=75
        integer, parameter :: ncoeff=75

	integer Sbar,Tbar,kidmax,numz,totp,group,ktp
	real(8) vseed(3)
        integer vl,nedraw,nmc,FTbar
	integer simul,nsimul
	real(8) disc
        
        character*200::pathoutput

        integer myid,numprocs,mptype

	integer i,j,kk,cntr,intmat(75,75),numint,ind
	integer intv(numint,2),jj

	integer xx,zz,ww,yy,x,reform
	real(8) rsqv(numper,1),rsqcvv(numper,1)
	real(8) thetam(numper,3,numz),sevm(numper,numz)
	character*16 vn(ncoeff)

	! This subroutine generates two sets of output files
	! One set of files contains the estimated coefficients
	! for each type. The other contains the thetas, which can
	! be read into the program in later runs. 

	! first put the location of the interactions into intv
	! so they can be put to the output file
	
        write(*,*) 'writing coefficients to output files'

	intv=0
	ind = 0
	loopi: do i = 1,ncoeff
	  loopj: do j = i,ncoeff
	     if (intmat(i,j).eq.1) then
	        ind = ind +1
	        intv(ind,1)=i
	        intv(ind,2)=j
	     end if
	     
	  end do loopj
	end do loopi


	! specify labels of variables

	vn(1)=''
	vn(2)=''
	vn(3)=''
	vn(4)=''
	vn(5)=''
	vn(6)=''
	vn(7)=''
	vn(8)=''
	vn(9)=''
	vn(10)=''
	vn(11)=''
	vn(12)=''
	vn(13)=''
	vn(14)=''
	vn(15)=''
	vn(16)=''
	vn(17)=''
	vn(18)=''
	vn(19)=''
	vn(20)=''
	vn(21)=''
	vn(22)=''
	vn(23)=''
	vn(24)=''
	vn(25)=''
	vn(26)=''
	vn(27)=''
	vn(28)=''
	vn(29)=''
	vn(30)=''
	vn(31)=''
	vn(32)=''
	vn(33)=''
	vn(34)=''
	vn(35)=''
	vn(36)=''
	vn(37)=''
	vn(38)=''
	vn(39)=''
	vn(40)=''
	vn(41)=''
	vn(42)=''
	vn(43)=''
	vn(44)=''
	vn(45)=''
	vn(46)=''
	vn(47)=''
	vn(48)=''
	vn(49)=''
	vn(50)=''
	vn(51)=''
	vn(52)=''
	vn(53)=''
	vn(54)=''
	vn(55)=''
	vn(56)=''
	vn(57)=''
	vn(58)=''
	vn(59)=''
	vn(60)=''
	vn(61)=''
	vn(62)=''
	vn(63)=''
	vn(64)=''
	vn(65)=''
	vn(66)=''
	vn(67)=''
	vn(68)=''
	vn(69)=''
	vn(70)=''
	vn(71)=''
	vn(72)=''
	vn(73)=''
	vn(74)=''
	vn(75)=''

        ! write to different files depending on pre or post reform
        if (mptype.eq.1) then 
	if (reform.eq.0) then
	    open(unit=300,file=adjustl(trim(pathoutput))//'../inputs/pretheta',status='UNKNOWN',action='WRITE')
            rewind(unit=300)
          else
	    open(unit=300,file=adjustl(trim(pathoutput))//'../inputs/posttheta',status='UNKNOWN',action='WRITE')
            rewind(unit=300)
          end if
	endif

!	if (mptype.eq.1) then 
 !	  !open(unit=200,file=adjustl(trim(pathoutput))//'output.type1',status='UNKNOWN',action='WRITE')
  !        !rewind(unit=200)
   !       if (reform.eq.0) then
!	    open(unit=300,file=adjustl(trim(pathoutput))//'pretheta.type1',status='UNKNOWN',action='WRITE')
 !           rewind(unit=300)
  !        else
!	    open(unit=300,file=adjustl(trim(pathoutput))//'posttheta.type1',status='UNKNOWN',action='WRITE')
 !           rewind(unit=300)
  !        end if
!	end if

!	if (mptype.eq.2) then 
!	!open(unit=200,file=adjustl(trim(pathoutput))//'output.type2',status='UNKNOWN',action='WRITE')
 !       !rewind(unit=200)
  !        if (reform.eq.0) then
!	    open(unit=300,file=adjustl(trim(pathoutput))//'pretheta.type2',status='UNKNOWN',action='WRITE')
 !           rewind(unit=300)
  !        else
!	    open(unit=300,file=adjustl(trim(pathoutput))//'posttheta.type2',status='UNKNOWN',action='WRITE')
 !           rewind(unit=300)
  !        end if
!	end if

 !	if (mptype.eq.3) then 
!	!open(unit=200,file='output.type3',status='UNKNOWN',action='WRITE')
 !       !rewind(unit=200)
  !        if (reform.eq.0) then
!	     open(unit=300,file=adjustl(trim(pathoutput))//'pretheta.type3',status='UNKNOWN',action='WRITE')
 !            rewind(unit=300)
  !        else
!	     open(unit=300,file=adjustl(trim(pathoutput))//'posttheta.type3',status='UNKNOWN',action='WRITE')
 !            rewind(unit=300)
  !        end if
!	end if

!
!
!	write(200,*) 'total number of histories:',totp
!	write(200,*) 'number of state variables (max):',numz
!	write(200,*) 'number histories used in calculating emax', &
!        nedraw
!	write(200,*) '--------------------------'
!        write(200,*) 'Estimated thetas'
!
!        Loopth1: do i = (Tbar-15),(35-15),-1  
!
!	
!         	 write(200,*) '--------------------------'
!		 write(200,*) 'mothers age',i
! 
!         
!      loopc: do jj = 1,7
!	     kk = (jj-1)*4+1
!
!	     if (jj.le.6) then	    
!	     write(200,*) '      ',vn(kk),vn(kk+1),vn(kk+2),vn(kk+3)
!	     write(200,16) (thetam(i,j),j=kk,kk+3)
! 	     write(200,16) (sevm(i,j),j=kk,kk+3)
!	     else
!	     write(200,*) '      ',vn(kk)
!	     write(200,16) (thetam(i,j),j=kk,kk)
! 	     write(200,16) (sevm(i,j),j=kk,kk)
!		 end if
!	     write(200,*) ' '
!
!	end do loopc
!
!	if (ind.gt.0) then
!	write(200,*) 'INTERACTION TERMS'
!	kk = 1
!	jj = ncoeff+1
!	x = numint
!	do while (x.ge.4)
!	      write(200,21) (intv(j,1),intv(j,2),'     ',j=kk,kk+3)
!	      write(200,16) (thetam(i,j),j=jj,jj+3)
!	      write(200,16) (sevm(i,j),j=jj,jj+3)
!	      write(200,*) ' '
!	      x=x-4
!	      jj = jj+4
!	      kk=kk+4
!	end do
!
!	if (x.gt.0) then
! 	write(200,21) (intv(j,1),intv(j,2),'     ',j=kk,kk+x-1)
!	write(200,16) (thetam(i,j),j=kk,kk+x-1)
!	write(200,16) (sevm(i,j),j=kk,kk+x-1)
!	write(200,*) ' '
!	end if
!	    
!	end if
!
!	write(200,*) '------------------------ --'
!      write(200,*) 'R-square, CV R-square'
!	write(200,18) rsqv(i,1),rsqcvv(i,1)		
!      write(200,*) '--------------------------'
!	     
!	    
!	end do Loopth1
!
!	
!	close(unit=200)


!	write(*,*) 'finished writing output to output files'


20    format(5a16)
21    format(i10,i3,a10,i3,i3,a10,i3,i3,a10,i3,i3,a10,i3,i3,a10)
16    format(4f16.4)
17    format(3f16.4)
18    format(f16.4)
19    format(40f16.4)



	xxloop: do xx = 1,numper
	grouploop: do group =1,3
	zzloop: do zz = 1,numz

	    write(300,*) thetam(xx,group,zz)

    	end do zzloop 
    	end do grouploop
	end do xxloop
	
        if (mptype.eq.ktp) close(unit=300)


22      format(f24.8)


	return 
        end subroutine ReportRes
     
     


!----------------------------- GETTOGV --------------------------
! This subroutine creates a vector of indicator variables
! for whether the column of the Z matrix contains nonzero
! values (i.e. whether it should be used in the regression)


    subroutine GetTogV(ct,Z,togv,nnz,intmat,regind,totp,numint,numv,silent)

    implicit none

    integer silent,ct,i,j,nnz,mage,totp,numv
    integer intmat(75,75),ind,togv(nnz)
    integer regind(75),numint,A,B
    real(8) rk
    real(8) Z(totp,nnz)

    togv = 0
    rk = 0.0d-0

!Toggle on all non-zero columns

    Loop1: do j = 1,nnz
        i = 1
        do while (i.le.totp.and.togv(j).eq.0)  !this loop gives a value togv(j)=1 to any j column of Z where at least one element is nonzero 
            if (Z(i,j).ne.rk) then
                togv(j)=1 
            end if
            i = i+1
        end do
    end do Loop1

!Toggle off columns with only one non-zero element (in case they are later interacted this would create collinearity)
    
    do j=1,nnz
    A=count(Z(:,j).ne.0.0d-0)
    if(A.gt.0.and.A.lt.10) then
        togv(j)=0
        if (silent.ne.1) print*,"too few non-zeros for regressor",j
    endif
    enddo
    
! "manually" toggle off some colinear regressors:

    do j=1,4
    B=count(Z(:,35+j).ne.1.0d-0.and.Z(:,34+j).eq.1.0d-0)
    if(B.lt.10) then
        togv(30+j)=0
        togv(77+j)=0
        if (silent.ne.1) print*,"noone in apvsegment",30+j
    endif
    enddo
    
! Toggle off regressors as read from input file and stored in "regind"

    bloop: do i = 1,75
        if (regind(i).eq.0) then
            togv(i)=0
        end if
    end do bloop


! Compute number of regressors actually used in regression
    numv = sum(togv) 

16      format(7f8.4)   

    return
    end subroutine GetTogV


    !------------------------ GETTHETA ------------------------------------

     subroutine GetTheta(o,ct,Z,emaxv,theta,togv,sevec, Rsq,nedraw,totp,numz,numv,flag1,pathoutput,group,typ,ref,silent)

    ! This subroutine takes the state values and the calculated
    ! emaxv values and returns the thetam parameters. This theta 
    ! values are then used to construct the emax values for the
    ! previous time period

    implicit none

 !   external matpd,dmatinv 
    
    integer ct,totp,simul,nsimul,nedraw,group,typ,ref,silent
    real(8) disc
    character*200::pathoutput
    integer i,j,ind,togv(numz),ii,numz,numv,o
    integer tt,info,flag1
    real(8) emaxv(totp,1),rk
    real(8) theta(1,numz),semaxv(nedraw,1)
    real(8) predval(nedraw,1),res(nedraw,1)
    real(8) sigmasq,Rsq
    real(8) sevec(1,numz),SST
    real(8) SSE1,SST1,ybar1(nedraw,1)
    real(8) Z(totp,numz),ZS(nedraw,numz)
    real(8) xpxi(numv,numv)
    real(8) xpx(numv,numv),xpxt(numv,numv)
    real(8) xpy(numv,1),checki(numv,numv)
    real(8) temp(numv,1)
    real(8) X(nedraw,numv)
    real(8) XP(numv,nedraw)
    real(8) CX(nedraw,numv)
    real(8) sum1,dist1
    integer w(numv),k,jj
    character(len=1024) blah,blah2,blah3,blah4

    56  format(2f10.2)

! initialize

    X = 0.0d-0
    XP=0.0d-0
    CX = 0.0d-0
    temp = 0.0d-0
    xpx = 0.0d-0
    xpy = 0.0d-0
    xpxi = 0.0d-0
    checki = 0.0d-0
    theta = 0.0d-0
    sevec=0.0d-0

! subset of data to use for regression 

    ZS = Z(1:nedraw,1:numz) 
    

! Populate vector X with only non-excluded regressors:

    tt = 0	
    Loop3: do j = 1,numz
       if (togv(j).eq.1) then
          tt = tt+1
          if (tt.gt.numz) then 
             write(*,*) 'error - beyond length'
          end if
          X(1:nedraw,tt) =ZS(1:nedraw,j)
          XP(tt,1:nedraw)=ZS(1:nedraw,j)	
       end if
    end do Loop3 

     
! compute the X'X and X'Y matrices, passing only the nonzero
! columns (or rows in the case of X-prime)
! store results in xpx and xpy

    semaxv(1:nedraw,1)=emaxv(1:nedraw,1)

    call matpd(numv,nedraw,numv,XP,X,xpx)
    call matpd(numv,nedraw,1,XP,semaxv,xpy)

! solve for beta in (X'X) * beta = X'Y 

    temp = xpy
    xpxi=xpx
    xpxt = xpx
    call dmatinv(xpxt,xpxi,numv)	


! check that we got the inverse

    call matpd(numv,numv,numv,xpxi,xpx,checki)
    flag1 = 0
    iloop: do i = 1,numv
    dist1 = dabs(checki(i,i)-1.0d-0)
    if (dist1.gt.0.001d-0) then
        flag1 = 1 
        !write(*,*) '**** problem with the inverse', dist1
        cycle iloop
        ! Output emax and regressors (before exclusion of 0s and collinear vectors) for stata analysis 
    end if
    end do iloop
    
    if (flag1.eq.1) then
        write(blah,98) ct
        write(blah2,99) group
        write(blah3,99) ref
        write(blah4,99) typ
        open(unit=589,position='rewind',file=adjustl(trim(pathoutput))//'emax'//trim(blah)//trim(blah2)//trim(blah3)//trim(blah4)//'.txt')
        write(589,*)0,0,0, togv(:)
        do i=1,totp
            write(589,*)i,ct, emaxv(i,1),Z(i,:)
        enddo
        close(589)
        flag1=0
        98 format(I2)
        99 format(I1)	
    endif

! if we got the correct inverse, calculcate xpxi*xpy
! and store in temp

    call matpd(numv,numv,1,xpxi,xpy,temp)

! read estimated regression coefficients into the theta matrix,
! where now the elements correspond to the original vectors of Z
! and we assign a coefficient equal to 0 for the columns of Z
! that were all zero

    ind = 1
    sevec=0.0d-0
    Loop4: do i = 1,numz
    if (togv(i).eq.1) then
       theta(1,i)=temp(ind,1)
       sevec(1,i)=xpxi(ind,ind)
       ind = ind+1
    end if
    end do Loop4

! calculate the predicted values from the regression

    call matpd(nedraw,numv,1,X,temp,predval)

! calculate the residuals and the error and
! total sum of squares

    SSE1 = 0.0d-0
    SST1 = 0.0d-0
    ybar1(:,1)=sum(semaxv(:,1))/nedraw
    res = 0.0d-0
    Loop11: do i = 1,nedraw
     res(i,1)=semaxv(i,1)-predval(i,1)
     SSE1 = SSE1 + res(i,1)*res(i,1)
     SST1 = SST1 + (semaxv(i,1)-ybar1(i,1))*(semaxv(i,1)-ybar1(i,1))
    end do Loop11

!compute standard errors on the regression coefficients

    sigmasq = SSE1 / (nedraw-numv)
    Loop10: do i = 1,numz
       sevec(1,i) = dsqrt(sigmasq*sevec(1,i))
    end do Loop10

    Rsq=0.0d-0
    rk = 0.0d-0

    Rsq = 1.0d-0 - (SSE1/SST1)


    return
    end subroutine GetTheta



    !****************************************************
    ! Subroutine that does matrix multiplication
    !****************************************************


    subroutine matpd(n1,m1,k1,a1,b1,c1)

    implicit none

    integer i,j,ii,n1,m1,k1
    real(8) a1(n1,m1),b1(m1,k1),c1(n1,k1)
    real(8) s1

    do 91 i =1,n1
     do 90 j = 1,k1
      s1 = 0.0d-0
      do 80 ii=1,m1
         s1 = s1+a1(i,ii)*b1(ii,j)
    80        continue
      c1(i,j)=s1
    90       continue
    91    continue

    return
    end subroutine matpd

    
    subroutine MSM(nsim,nmoments,n1,n2,criterion,simdata1,simdata2,simdiv,loss,pathoutput,estimation,counterchar)
 
    implicit none
    integer, parameter:: noutcomes=8
    integer, parameter:: nconditions = 15
    integer, parameter:: numfam = 5314
    integer, parameter:: nyrs = 5
    
    character*6:: spec
    integer i,j,k,l                                                  !counter
    integer nsim                                                    !total number of year-household-clone simulations
    integer nclones                                                 !number of household clones used in the simulation
    integer nhh                                                     !number of households simulated
    integer nmoments                                                !number of moments used in estimation
    integer cohort(nsim)
    real(8) temp(nmoments,2*noutcomes+2*nconditions)
    logical temp6(1,noutcomes), temp5(1,nconditions), temp7(nsim), temp8(nsim)
    real(8) Y(nsim,noutcomes),outcomes(nsim,noutcomes),X(nsim), upperbound(nsim,nmoments,noutcomes), lowerbound(nsim,nmoments,noutcomes) 
    real(8) conditions(nsim,nconditions), conditions_ub(nsim,nmoments,nconditions), conditions_lb(nsim,nmoments,nconditions) 
    character*16:: condchar_lb(nmoments,nconditions)
    character*16:: condchar_ub(nmoments,nconditions)
    integer female(nsim),brou(nsim),bra(nmoments),bri(noutcomes,nsim)
    integer nyears                                                  !number of years during which households are simulated
    integer n1,n2                                             !number of integer and real variables in the simulated data
    real(8) simdata2(nsim,n2)                                    !array containing the simulated data (real variables)
    integer simdata1(nsim,n1)                                    !array containing the simulated data (integer variables)
    real(8) firstdiff(nsim,2)                                         !array containing the first differences in wages for men and women
    integer simhhage(nsim),simedugrp_w(nsim),simedugrp_h(nsim) ,simXPgr_w(nsim),simXPgr2_w(nsim),simXPgr_h(nsim)  !array containing the schooling and XP of the simulated agents grouped in  bins
    character(200) momentnames(nmoments)                            !array containing the description of the  moments used in estimation
    real(8) mdata(nmoments)                                         !array containing the data moments
    real(8) weights(nmoments)                                       !array containing the weights used in computing the MSM criterion function
    real(8) var(nmoments)                                           !array containing the variances of the data moments
    integer nobs(nmoments)                                          !array containing the number of observations contributing to a given moment
    real(8) msim(nmoments)                                          !array containing the simulated moments
    real(8) criterion                                               !contains the value of the MSM criterion function evaluated at current parameter values
    real(8) loss(1,nmoments),weightedloss(1,nmoments)               !array containing the contribution of each moment to the MSM criterion function
    logical mask(nsim),mask1(nsim),mask2(nsim),mask1b(nsim),mask2b(nsim)                                 !these arrays are used to count simulated households that contribute to a given moment
    integer ntrue(nmoments)                                                   !number of observations that contribute positively to a given moment
    real(8) sumcont(1)                                              !sum of contribution to a given moments (to be divided by number of contributing obs) ex: sum of earnings of males with no HS in 2006
    integer blah, temp4, temp3,flag                                                   !buffer
    character(4) bleh
    integer onlymfwage,onlywage, simdiv(nsim), simkidgrp(nsim), simwealthgrp(nsim), simdensgrp_w(nsim), simdensgrp_h(nsim)
    real(8) simdens_h(nsim), simdens_w(nsim),diffbin(5),total, earnbin(6)
    character*200::pathoutput
    integer offmom(nmoments), onmom(nmoments),estimation,equalweights,bin,g
    integer collmom(nmoments),meanmom(nmoments),nbinmom,nmeanmom, D(nsim)
    character*1::tag,counterchar
    character*200:: momgroup(nmoments)
    real(8) cloneaverage(nmoments)
real(8) contribution(nsim,nmoments)

    !note: simdata1 contains
!	1 folio
!	2 year
!	3 married
!	4 sexsampled
!	5 woman's age (x if household is a single male)
!	6 man's age (x if household is a single female)
!	7 husband lfp dummy (x if household is a single female)
!	8 wife's lfp (2 means part-time work) (x if household is a single male)
!	9 husband's formal work dummy (x if household is a single female)
!	10 woman's formal work dummy (x if household is a single male)
!	11 number of kids
!	12 arrival of a kid this period (dummy)
!	13 woman's years of schooling
!	14 man's years of schooling

    !note: simdata2 contains
!	1 saving rate
!	2 woman's pension savings (x if household is a single male)
!	3 man's pension savings (x if household is a single female)
!	4 household private savings
!	5 man's total experience (x if household is a single female)
!	6 man's formal sector experience (x if household is a single male)
!	7 woman's total experience (x if household is a single male)
!	8 woman's formal sector experience (x if household is a single male)
!	9 man's annual earnings (x if household is a single female)
!	10 woman's annual earnings (x if household is a single male)
!	11 man's withdrawal from pension fund (including governement benefits) (x if household is a single female)
!   12 woman's withdrawal from pension fund (including governement benefits) (x if household is a single male)


    counterchar="o"
    !0a. Equal weights?

    equalweights=0
    
    !0b. Moments to turn-off

    offmom=0
    onlymfwage=0
    onlywage=0
    if (onlymfwage.eq.1) then
         offmom(1:183)=(/(i,i=1,103),(i,i=107,115),(i,i=120,131),(i,i=137,165),167,169,171,(i,i=173,199)/)    
    endif
    if (onlywage.eq.1) then
         offmom(1:103)=(/(i,i=1,103)/)    
    endif
     
    !0c. Moments to emphasize
    onmom=0
    !onmom(1:32)=(/(i,i=92,103),(i,i=160,179)/)
     
	!1. Read text file with moment definitions

    temp=0.0d-0  
    open(file="./inputs/momdef.txt", unit=1)
		read(1,*)
		do i=1,nmoments
			!print*,i
			read(1,*) blah, momgroup(i),collmom(i), meanmom(i),temp(i,:)
			!print*, momgroup(i),collmom(i), meanmom(i), temp(i,:)
		enddo
    close(1)

    do i=1,nsim
	do j=1,noutcomes
		lowerbound(i,:,j)=temp(:,2*j-1)*1.0d-0
		upperbound(i,:,j)=temp(:,2*j)*1.0d-0
    enddo
	do j=1,nconditions
		conditions_lb(i,:,j)=temp(:,2*noutcomes+2*j-1)*1.0d-0
		conditions_ub(i,:,j)=temp(:,2*noutcomes+2*j)*1.0d-0
	enddo 
    enddo
	
	
	! Fill array of outcomes and array of conditions with either data or simulations stored in simdata arrays
    
	female=0
    female(:)=simdata1(:,4)
    
    outcomes=0.0d-0
    outcomes(:,1)=simdata2(:,4)						 								!household assets
    outcomes(:,2)=1.0d-0*(female(:)*simdata1(:,8)+(1-female(:))*simdata1(:,7)) 		!LFP (sampled hh mb) 
    outcomes(:,3)=1.0d-0*(female(:)*simdata1(:,10)+(1-female(:))*simdata1(:,9))		!Formal work (spled hh mb) uncond
    outcomes(:,4)=(female(:)*simdata2(:,10)+(1-female(:))*simdata2(:,9))			!Earnings (spled hh mb)
 
	do i=1,2 
		outcomes(i,5)=-9999															!lagged Earnings (spled hh mb)
		outcomes(i,6)=-9999															!Earnings second-differences	    
    enddo
    do i=3,nsim
	if (simdata1(i-2,1).eq.simdata1(i,1)) then
		outcomes(i,5)=(female(i)*simdata2(i-2,10)+(1-female(i))*simdata2(i-2,9))	!lagged Earnings (spled hh mb)
		outcomes(i,6)=outcomes(i,4)-outcomes(i,5)									!Earnings second-differences	    
	endif
    enddo
	
    outcomes(:,7)=1.0d-0*(simdata1(:,7)+simdata1(:,8))								!nb of earners in the household
    outcomes(:,8)=1.0d-0*(simdata1(:,9)+simdata1(:,10))   							!nb of formal earners in the household
    
    conditions=0.0d-0
    conditions(:,1)=1.0d-0*(female(:)*simdata1(:,5)+(1-female(:))*simdata1(:,6))	!Age (spled hh mb)
    conditions(:,2)=1.0d-0*(simdata1(:,3))											!Married
    conditions(:,3)=1.0d-0*(female(:))												!Female
    conditions(:,4)=1.0d-0*(simdata1(:,11))											!Nb kids
    conditions(:,5)=1.0d-0*(female(:)*simdata1(:,13)+(1-female(:))*simdata1(:,14))	!Education in years (spled hh mb)
    conditions(:,6)=female(:)*simdata2(:,7)+(1-female(:))*simdata2(:,5)				!Total experience (spled hh mb)
    conditions(:,7)=female(:)*simdata2(:,8)+(1-female(:))*simdata2(:,6)				!Formal experience (spled hh mb)
    conditions(:,8)=1.0d-0*(simdata1(:,2)-conditions(:,1)-1900)						!Cohort (spled hh mb)
    conditions(:,9)=outcomes(:,2)													!LFP
    conditions(:,10)=outcomes(:,3)													!Formal work uncond
    
	do i=1,2 
		conditions(i,11)=-9999	  													!lagged LFP (sampled hh mb) 
		conditions(i,12)=-9999  													!lagged Covered work (spled hh mb) uncond
    enddo
    do i=3,nsim
	if (simdata1(i-2,1).eq.simdata1(i,1)) then
		conditions(i,11)=1.0d-0*(female(i)*simdata1(i-2,8) +(1-female(i))*simdata1(i-2,7))  !lagged LFP (sampled hh mb) 
		conditions(i,12)=1.0d-0*(female(i)*simdata1(i-2,10)+(1-female(i))*simdata1(i-2,9)) 	!lagged Covered work (spled hh mb) uncond
	endif
    enddo
	
    conditions(:,13)=simdata2(:,4)													!HH assets
    conditions(:,14)=1.0d-0*(simdata1(:,2))											!Year
    conditions(:,15)=1.0d-0*(simdata1(:,2)-2*int(simdata1(:,2)/2))    				!Odd/Even year	



	! Create moment names 
	
    do k=1,nmoments
    do j=1,nconditions
    	write(condchar_lb(k,j),'(f16.0)') conditions_lb(1,k,j)
    	write(condchar_ub(k,j),'(f16.0)') conditions_ub(1,k,j)
    enddo
    enddo

    do k=1,nmoments
    	momentnames(k)=trim(adjustl(momgroup(k)))//  &
&		'     '// &
&		'_Age_'//trim(adjustl(condchar_lb(k,1)))//"_"//trim(adjustl(condchar_ub(k,1)))// &
&		'_Mar_'//trim(adjustl(condchar_lb(k,2)))//"_"//trim(adjustl(condchar_ub(k,2)))// &
&		'_Fem_'//trim(adjustl(condchar_lb(k,3)))//"_"//trim(adjustl(condchar_ub(k,3)))// &
&		'_Kid_'//trim(adjustl(condchar_lb(k,4)))//"_"//trim(adjustl(condchar_ub(k,4)))// &
&		'_Edu_'//trim(adjustl(condchar_lb(k,5)))//"_"//trim(adjustl(condchar_ub(k,5)))// &
&		'_XP_'//trim(adjustl(condchar_lb(k,6)))//"_"//trim(adjustl(condchar_ub(k,6)))// &
&		'_FXP_'//trim(adjustl(condchar_lb(k,7)))//"_"//trim(adjustl(condchar_ub(k,7)))// &
&		'_Coh_'//trim(adjustl(condchar_lb(k,8)))//"_"//trim(adjustl(condchar_ub(k,8)))// &
&		'_LFP_'//trim(adjustl(condchar_lb(k,9)))//"_"//trim(adjustl(condchar_ub(k,9)))// &
&		'_FLFP_'//trim(adjustl(condchar_lb(k,10)))//"_"//trim(adjustl(condchar_ub(k,10)))// &
&		'_Ass_'//trim(adjustl(condchar_lb(k,11)))//"_"//trim(adjustl(condchar_ub(k,11)))// & 
&		'_Year_'//trim(adjustl(condchar_lb(k,12)))//"_"//trim(adjustl(condchar_ub(k,12))) 

    enddo

        !3. Compute the simulated moments
    flag=0

    do k=1,nmoments
	if (meanmom(k).eq.0) then
		ntrue(k)=count(all(outcomes(:,:).ge.lowerbound(:,k,:).and.outcomes(:,:).le.upperbound(:,k,:),2).and.all(conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2)) 
		nobs(k)=count(all(    conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2)) 
		msim(k)=ntrue(k)*1.0d-0/nobs(k)*1.0d-0
		var(k)=msim(k)*(1.0d-0-msim(k))/nobs(k)
		D(:)=-1*	all(    conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2) 
		contribution(:,k)=-1*all(outcomes(:,:).ge.lowerbound(:,k,:).and.outcomes(:,:).le.upperbound(:,k,:),2).and.all(conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2) 
		contribution(:,k)=contribution(:,k)*D(:)

	elseif (meanmom(k).eq.1) then
		Y(:,:)=outcomes(:,:)*lowerbound(:,k,:) ! for moments that are averages (not proportions), lowerbound contains an indicator for the particular outcome relevant to that moment
		X(:)=sum(Y(:,:),2)
		mask(:)=all(conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2).and.X(:).ne.-99.0d-0
		nobs(k)=count(mask) 
		msim(k)=sum(X,mask)/nobs(k)
		var(k)=sum((X-msim(k))**2,mask)/nobs(k)
		D(:)=-1*all(    conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2) 
		contribution(:,k)=X(:)*D(:)
   	endif
	enddo     


	if (counterchar=="Z") then					!create the data moments
		open(unit=333, file='./inputs/mdata.txt')
		write(333,*) "	c1	  c2	  c3"
		do k=1,nmoments
			write(333,'(I5,2f15.5,I7,A360)') k, msim(k), var(k), nobs(k), momgroup(k)
		enddo
 		close(333)
		endif
		
	if (counterchar=="X") then	
		open(unit=334, file='./outputs/contribution.txt')
		open(unit=335, file='./outputs/nobs.txt')
		do i=1,nsim
			write(334,*) simdata1(i,1), simdata1(i,2)
			write(335,*) simdata1(i,1), simdata1(i,2)
			do k=1,nmoments
				write(334,'(f19.6)',advance="no") contribution(i,k)
			enddo
			do k=1,nmoments
				write(335,'(I7)',advance="no") nobs(k)
			enddo
		enddo
		close(334)
		close(335)
	endif
	
	!4. Read in the data moments
    

    open(unit=4002,file='./inputs/mdata.txt',position='rewind')
    	read(4002,*)
    
	do i=1,nmoments
    		read(4002,*) bleh, mdata(i), var(i)
    		if (var(i).eq.0.0d-0) var(i)=10000000
    	enddo
    close(4002)
    var(:)=var(:)
    weights(:)=1.0d-0/var(:)

	
	!4. Calculate the MSM objective function

    !! Turn-off/amplify some moments if needed

    do k=1,nmoments
        if (nobs(k).eq.0) then
          msim(k)=0.0d-0
        endif
    enddo

    if (equalweights==1) weights=1.0d-0
    
    do k=1,nmoments
    	if (offmom(k).gt.0.and.offmom(k).le.nmoments) weights(offmom(k))=0.0d-0
!       if (onmom(k).gt.0.and.onmom(k).le.nmoments)   weightedloss(1,onmom(k))=10*weightedloss(1,onmom(k))
    enddo
    
    !! Compute difference between data and simulated moments, store in loss()
    
    loss(1,:)= mdata(1:nmoments) - msim(1:nmoments)
    do i=1,nmoments
    	weightedloss(1,i)=loss(1,i)*loss(1,i)*weights(i)
    enddo

    criterion=sum(weightedloss) 
     
    if (estimation.ne.1) then
     
        !5. Output moment names, simulated moments, data moments and weights

    	open(unit=4004, file=trim(adjustl(pathoutput))//'criterion.txt',position='rewind')
        	write(4004,*) 'i' , ',' , ' m(i) ' , ',' , ' mdata(i) ' , ',' , ' weights(i) ' , ',' , ' loss(1,i) ' , ',' , 'nobs(i)' 
        	do i=1,nmoments
        		write(4004,'(I4,A1,f20.8,A1,f20.8,A1,f20.8,A1,f20.8,A1,I7,A1,A200)') i,",", msim(i),",", mdata(i),",", weights(i),",", loss(1,i),",", nobs(i),",",momentnames(i)
		enddo
    	close(4004)

	
    	open(unit=4004, file=trim(adjustl(pathoutput))//'moments.txt',position='rewind')
        	do i=1,nmoments
        		write(4004,*) msim(i)
		enddo
    	close(4004)
    endif
  


26      format(i7,i6,i5,i5,i5,i5,i2,i2,i2,i2,f6.2,3f16.2,i5,i5,4f8.2,4f16.2,i5,i5)
   
   endsubroutine
   

        subroutine createinputfile(parms,numiter,pathoutput,alabel,ktp)

        implicit none

                character*200::pathoutput

	integer, parameter :: maxnumtypes=4
        integer totp,Tbar,Sbar,kidmax
	integer nedraw,nmc,FTbar,vl,ncoeff2
	integer simul,nsimul,ncoeff
	integer dnumdraw,snumdraw,i,ktp
        real(8) vseed(3)
	real(8) disc
	real(8) ftol,rtol,temp1
	integer maxit,temp2
        real(8), allocatable :: alpha(:,:)
        real(8), allocatable :: alphalow(:,:)
        real(8), allocatable :: alphahigh(:,:)
        real(8), allocatable :: talpha(:,:)
        real(8), allocatable :: abump(:,:)
        integer, allocatable :: aiter(:,:)
        integer j,numper,numiter
        real(8) tausm,tempvar
        integer checkp,numint,f1ind,sind,regind(75),intmat(75,75)
        integer numz,jj,kk
        real(8) tparms(numiter)
        real(8) parms(numiter)
        real(8) scale(numiter)
        real(8),allocatable :: rvart(:,:,:)
        integer,allocatable :: ivart(:,:,:)   
        integer, allocatable :: fixvar(:,:)        
        real(8), allocatable :: maxsafe(:)        
        character*50 :: alabel(300)
	character*10 ::  macro(300,maxnumtypes)
 
        ncoeff = 150	


        allocate(alpha(ncoeff,maxnumtypes))      
        allocate(alphalow(ncoeff,maxnumtypes))      
        allocate(alphahigh(ncoeff,maxnumtypes))      
        allocate(talpha(ncoeff,maxnumtypes))      
        allocate(abump(ncoeff,maxnumtypes))      
        allocate(aiter(ncoeff,maxnumtypes))      

        alpha = 0.0d-0
        abump = 0.0d-0

        totp = 0
        kidmax = 0
        numper = 0
        vseed = 0
        numiter = 0
        alpha = 0.0d-0
        abump = 0.0d-0
        aiter = 0
        Tbar = 0
        Sbar = 0

        call Readin(nedraw,totp,kidmax,ktp,numper,vseed,alpha,abump,aiter,numiter,Sbar,Tbar,nmc,&
             regind,intmat,numint,nsimul,simul,dnumdraw,snumdraw,tausm,rtol,ftol,maxit,alphalow,alphahigh)


        numz = 75+numint

        allocate(rvart(6,totp,numper))      
        allocate(ivart(9,totp,numper))   
        allocate(fixvar(4,totp))
        allocate(maxsafe(numper))



        talpha = alpha

        ncoeff2 = ncoeff-numint
	


        ! update the values of the parameters to reflect latest values

	i = 0
	if (sum(aiter).gt.0) then
	jj1: do jj = 1,ncoeff
	 kk1: do kk = 1,ktp
	    if (aiter(jj,kk).eq.1) then
	       i = i + 1
	       talpha(jj,kk)=parms(i)
	    end if
	 end do	kk1 !jj
	end do	jj1 !kk
	end if !iterating on values of alpha


!       write out all the ingredients to create a new input file 

	macro(:,:)=""

        alabel(1) = 'alpha(1) - crra coefficient'; macro(1,:)="\crra"
        alabel(2) = 'alpha(2) - stock of children in female marg ut c (nu_0)'
        alabel(3) = 'alpha(3) - stock of children in male marg ut c (nu_0)'
        alabel(4) = 'alpha(4) - leisure female in marg ut c'
        alabel(5) = 'alpha(5) - leisure male in marg ut c'
        alabel(6) = 'alpha(6) - util of leisure female (delta_l)'
        alabel(7) = 'alpha(7) - util of leisure male (delta_l)'
        alabel(8) = 'alpha(8) - util of leisure if work parttime (delta_l)'
        alabel(9) = 'alpha(9) - util of leisure female*nkids*wooman age < 50'
        alabel(10) = 'alpha(10) - util of leisure female*nkids'
        alabel(11) = 'alpha(11) - util of leisure female*married'
        alabel(12) = 'alpha(12) - '
        alabel(13) = 'alpha(13) - man not working age 65+'
        alabel(14) = 'alpha(14) - women not working age 60+'
        alabel(15) = 'alpha(15) - discount rate'
        alabel(16) = 'alpha(16) - ** husbands log income formal - constant term'
        alabel(17) = 'alpha(17) - years of education'
        alabel(18) = 'alpha(18) - experience'
        alabel(19) = 'alpha(19) - experience squared'
        alabel(20) = 'alpha(20) - education squared'
        alabel(21) = 'alpha(21) - ** husbands log income informal - constant'
        alabel(22) = 'alpha(22) - years of education'
        alabel(23) = 'alpha(23) - experience'
        alabel(24) = 'alpha(24) - experience squared'
        alabel(25) = 'alpha(25) - education squared'
        alabel(26) = 'alpha(26) - ** female log income formal - constant'
        alabel(27) = 'alpha(27) - years of education'
        alabel(28) = 'alpha(28) - experience'
        alabel(29) = 'alpha(29) - experience squared'
        alabel(30) = 'alpha(30) - education squared'
        alabel(31) = 'alpha(31) - ** female log income informal - constant'
        alabel(32) = 'alpha(32) - years of education'
        alabel(33) = 'alpha(33) - experience'
        alabel(34) = 'alpha(34) - experience squared'
        alabel(35) = 'alpha(35) - education squared'
        alabel(36) = 'alpha(36) - ** prob covered sector male offer constant'
        alabel(37) = 'alpha(37) - years of schooling'
        alabel(38) = 'alpha(38) - prev period covered'
        alabel(39) = 'alpha(39) - male age'
        alabel(40) = 'alpha(40) - ** prob covered sector female offer constant'
        alabel(41) = 'alpha(41) - years of schooling'
        alabel(42) = 'alpha(42) - prev period covered'
        alabel(43) = 'alpha(43) - female age'
        alabel(44) = 'alpha(44) - ** prob no fert constant'
        alabel(45) = 'alpha(45) - age of woman'
        alabel(46) = 'alpha(46) - married'
        alabel(47) = 'alpha(47) - numk'
        alabel(48) = 'alpha(48) - ** prob no divorce'
        alabel(49) = 'alpha(49) - male age'
        alabel(50) = 'alpha(50) - hed'
        alabel(51) = 'alpha(51) - age diff'
        alabel(52) = 'alpha(52) - duration of marriage'
        alabel(53) = 'alpha(53) - duration squared'
        alabel(54) = 'alpha(54) - fem age'
        alabel(55) = 'alpha(55) - male age'
        alabel(56) = 'alpha(56) - blank'
        alabel(57) = 'alpha(57) - blank'
        alabel(58) = 'alpha(58) - blank'
        alabel(59) = 'alpha(59) - blank'
        alabel(60) = 'alpha(60) - blank'
        alabel(61) = 'alpha(61) - ** more utility parameters'
        alabel(62) = 'alpha(62) - male cost of switching sectors' 
        alabel(63) = 'alpha(63) - female cost of switching sectors'
        alabel(64) = 'alpha(64) - male cost of returning to work'
        alabel(65) = 'alpha(65) - female cost of returning to work'
        alabel(66) = 'alpha(66) - '
        alabel(67) = 'alpha(67) - fertility: woman years ed'
        alabel(68) = 'alpha(68) - fertility: married*kids'
        alabel(69) = 'alpha(69) - '
        alabel(70) = 'alpha(70) - non-pecuniary benefits informal sector (women)'
        alabel(71) = 'alpha(71) - non-pecuniary benefits informal sector (men)'
        alabel(72) = 'alpha(72) - '
        alabel(73) = 'alpha(73) - '
        alabel(74) = 'alpha(74) - '
        alabel(75) = 'alpha(75) - Type prob'
        alabel(76) = 'alpha(76) - wed' 
        alabel(77) = 'alpha(77) - hed '
        alabel(78) = 'alpha(78) - married'
        alabel(79) = 'alpha(79) - cohort'
        alabel(80) = 'alpha(80) - '
        alabel(81) = 'alpha(81) - bequest parameter 1'
        alabel(82) = 'alpha(82) - '
        alabel(83) = 'alpha(83) - bequest parameter 2'
        alabel(84) = 'alpha(84) - '
        alabel(85) = 'alpha(85) - consumption floor'
        alabel(86) = 'alpha(86) - '
        alabel(87) = 'alpha(87) - '
        alabel(88) = 'alpha(88) - '
        alabel(89) = 'alpha(89) - '
        alabel(90) = 'alpha(90) - '
        alabel(91) = 'alpha(91) - '
        alabel(92) = 'alpha(92) - '
        alabel(93) = 'alpha(93) - '
        alabel(94) = 'alpha(94) - '
        alabel(95) = 'alpha(95) - '
        alabel(96) = 'alpha(96) - '
        alabel(97) = 'alpha(97) - '
        alabel(98) = 'alpha(98) -  ***** borrowing constraint'
        alabel(99) = 'alpha(99) - '
        alabel(100) = 'alpha(100) - '
        alabel(101) = 'alpha(101) - '
        alabel(102) = 'alpha(102) - '
        alabel(103) = 'alpha(103) - '
        alabel(104) = 'alpha(104) - mean retun of asset'
        alabel(105) = 'alpha(105) - variance of shock e1: shock to husbands formal income'
        alabel(106) = 'alpha(106) - variance of shock e2: shock to husbands informal income'
        alabel(107) = 'alpha(107) - variance of shock e3: shock to wifes formal income'
        alabel(108) = 'alpha(108) - variance of shock e4: shock to wifes informal income'
        alabel(109) = 'alpha(109) - variance of shock e5: shock to asset'
        alabel(110) = 'alpha(110) - variance of shock e6: shock to hus leisure'
        alabel(111) = 'alpha(111) - variance of shock e5: shock to wife leisure'
        alabel(112) = 'alpha(112) - correlation of shocks e1 and e2'
        alabel(113) = 'alpha(113) - correlation of shocks e1 and e3'
        alabel(114) = 'alpha(114) - correlation of shocks e1 and r4'
        alabel(115) = 'alpha(115) - correlation of shocks e1 and e5'
        alabel(116) = 'alpha(116) - correlation of shocks e1 and e6'
        alabel(117) = 'alpha(117) - correlation of shocks e1 and e7'
        alabel(118) = 'alpha(118) - correlation of shocks e2 and e3'
        alabel(119) = 'alpha(119) - correlation of shocks e2 and e4'
        alabel(120) = 'alpha(120) - correlation of shocks e2 and e5'
        alabel(121) = 'alpha(121) - correlation of shocks e2 and e6'
        alabel(122) = 'alpha(122) - correlation of shocks e2 and e7'
        alabel(123) = 'alpha(123) - correlation of shocks e3 and e4'
        alabel(124) = 'alpha(124) - correlation of shocks e3 and e5'
        alabel(125) = 'alpha(125) - correlation of shocks e3 and e6'
        alabel(126) = 'alpha(126) - correlation of shocks e3 and e7'
        alabel(127) = 'alpha(127) - correlation of shocks e4 and e5'
        alabel(128) = 'alpha(128) - correlation of shocks e4 and e6'
        alabel(129) = 'alpha(129) - correlation of shocks e4 and e7'
        alabel(130) = 'alpha(130) - correlation of shocks e5 and e6'
        alabel(131) = 'alpha(131) - correlation of shocks e5 and e7'
        alabel(132) = 'alpha(132) - correlation of shocks e6 and e7'

        nedraw = 3000

        open(unit=50, file=adjustl(trim(pathoutput))//'parammacro.txt'	,status='UNKNOWN',action='WRITE',position='REWIND')
        open(unit=100,file=adjustl(trim(pathoutput))//'newinputfile.txt',status='UNKNOWN',action='WRITE',position='REWIND')

        write(100,*) totp,'Potential number of states' 
        write(100,*) Tbar,'T-bar' 
        write(100,*) Sbar,'S-bar' 
        write(100,*) kidmax,'Total number of kids possible' 
        write(100,*) nedraw,'Actual number of states used in calculating emax' 
        write(100,*) nmc,'Number of Monte Carlo draws used in emax' 
        write(100,*) vseed(1),'Seed for random epsilon draws' 
        write(100,*) vseed(2),'Seed for generating simulated people' 
        write(100,*) vseed(3),'Seed for other use' 
        write(100,*) ktp,'Number of types' 
        write(100,*) alpha(15,1),'discount rate' 
        write(100,*) 75,'max number of regressors in emax approx' 
        write(100,*) simul,'indicator for whether to do the simulation part' 
        write(100,*) nsimul,'number of people to simulate' 
        write(100,*) dnumdraw,'number of draws in the likelihood estimation part' 
        write(100,*) snumdraw,'number of draws in simulation based on data' 
        write(100,*) tausm,'tau smoothing parameter' 
        write(100,*) 1.0,'contraction factor (simplex parameter)' 
        write(100,*) ftol,'FTOL - simplex parameter' 
        write(100,*) rtol,'RTOL - simplex parameter' 
        write(100,*) maxit,'maximum number of simplex iterations' 


	Loop4: do i=1,132
	 							  
	   if (i.eq.1.or.i.eq.6.or.i.eq.7.or.i.eq.16.or.i.eq.21.or.i.eq.26.or.i.eq.31.or.i.eq.70.or.i.eq.71.or.i.eq.75.or.i.eq.76.or.i.eq.77.or.i.eq.78.or.i.eq.79) then  
	     Loop7: do j = 1,maxnumtypes	         
               write(100,18) talpha(i,j),alphalow(i,j),alphahigh(i,j),abump(i,j),aiter(i,j),alabel(i)
	       write(50,19)"\newcommand{",macro(i,j),"}{",talpha(i,j),"}  %",alabel(i)
	     end do Loop7
	     else
	       write(100,18) talpha(i,1),alphalow(i,1),alphahigh(i,1),abump(i,1),aiter(i,1),alabel(i)
	       write(50,19)"\newcommand{",macro(i,1),"}{",talpha(i,1),"}  %",alabel(i)
	   end if
	   
	end do Loop4
	

        write(100,*) 0.0,'*********************************'	
	write(100,*) 0.0,'*  specify regressors           *'
	write(100,*) 0.0,'*  to include in the emax       *'
	write(100,*) 0.0,'*  calculation (1=yes,0=no)     *'
        write(100,*) 0.0,'*********************************'	


	!write out the indicators for whether to include each regressor

        write(100,12) regind(1),'intercept                1'
        write(100,12) regind(2),'                         2'
        write(100,12) regind(3),'                         3'
        write(100,12) regind(4),'                         4'
        write(100,12) regind(5),'                         5'
        write(100,12) regind(6),'                         6'
        write(100,12) regind(7),'                         7'
        write(100,12) regind(8),'                         8'
        write(100,12) regind(9),'                         9'
        write(100,12) regind(10),'                       10'
        write(100,12) regind(11),'                       11'
        write(100,12) regind(12),'                       12'
        write(100,12) regind(13),'                       13'
        write(100,12) regind(14),'                       14'
        write(100,12) regind(15),'                       15'
        write(100,12) regind(16),'                       16'
        write(100,12) regind(17),'                       17'
        write(100,12) regind(18),'                       18'
        write(100,12) regind(19),'                       19'
        write(100,12) regind(20),'                       20'
        write(100,12) regind(21),'                       21'
        write(100,12) regind(22),'                       22'
        write(100,12) regind(23),'                       23'
        write(100,12) regind(24),'                       24'
        write(100,12) regind(25),'                       25'
        write(100,12) regind(26),'                       26'
        write(100,12) regind(27),'                       27'
        write(100,12) regind(28),'                       28'
        write(100,12) regind(29),'                       29'
        write(100,12) regind(30),'                       30'
        write(100,12) regind(31),'                       31'
        write(100,12) regind(32),'                       32'
        write(100,12) regind(33),'                       33'
        write(100,12) regind(34),'                       34'
        write(100,12) regind(35),'                       35'
        write(100,12) regind(36),'                       36'
        write(100,12) regind(37),'                       37'
        write(100,12) regind(38),'                       38'
        write(100,12) regind(39),'                       39'
        write(100,12) regind(40),'                       40'
        write(100,12) regind(41),'                       41'
        write(100,12) regind(42),'                       42'
        write(100,12) regind(43),'                       43'
        write(100,12) regind(44),'                       44'
        write(100,12) regind(45),'                       45' 
        write(100,12) regind(46),'                       46' 
        write(100,12) regind(47),'                       47' 
        write(100,12) regind(48),'                       48' 
        write(100,12) regind(49),'                       49' 
        write(100,12) regind(50),'                       50' 
        write(100,12) regind(51),'                       51' 
        write(100,12) regind(52),'                       52' 
        write(100,12) regind(53),'                       53' 
        write(100,12) regind(54),'                       54' 
        write(100,12) regind(55),'                       55' 
        write(100,12) regind(56),'                       56' 
        write(100,12) regind(57),'                       57' 
        write(100,12) regind(58),'                       58' 
        write(100,12) regind(59),'                       59' 
        write(100,12) regind(60),'                       60' 
        write(100,12) regind(61),'                       61' 
        write(100,12) regind(62),'                       62' 
        write(100,12) regind(63),'                       63' 
        write(100,12) regind(64),'                       64' 
        write(100,12) regind(65),'                       65' 
        write(100,12) regind(66),'                       66' 
        write(100,12) regind(67),'                       67' 
        write(100,12) regind(68),'                       68' 
        write(100,12) regind(69),'                       69' 
        write(100,12) regind(70),'                       70' 
        write(100,12) regind(71),'                       71' 
        write(100,12) regind(72),'                       72' 
        write(100,12) regind(73),'                       73' 
        write(100,12) regind(74),'                       74' 
        write(100,12) regind(75),'                       75' 

        write(100,*) 0.0,'**********************'	
	write(100,*) 0.0,'* interaction        *'
        write(100,*) 0.0,'* terms in emax      *'
        write(100,*) 0.0,'* approx             *'
        write(100,*) 0.0,'**********************'	


	
	! read in the location of the interactions

        intloop: do i = 1,75
          if (intmat(i,i).eq.1) then
             write(100,11) i,i
          end if
          int2loop: do j = (i+1),75
             if (intmat(i,j).eq.1) then
               write(100,11) i,j
             end if
          end do int2loop
        end do intloop
        write(100,11) -9,-9
	

	close(unit=100)
	close(unit=50)

11      format(i4,5x,i4)
12      format(i1,10x,a50)
13      format(i1,1x,i1,1x,i1,5x,a50)
14      format(f20.8,5x,a50)
15      format(i20,10x,a50)
16      format(f20.8,1x,f20.8,i4,5x,a50,a6,1x,i1)
17      format(f20.8,1x,f20.8,i4,5x,a50)
18      format(f20.8,1x,f20.8,1x,f20.8,1x,f20.8,i4,5x,a50)
19	format(A12,A10,A2,f10.3,A4,A50)
	return
	
        deallocate(alpha)      
        deallocate(alphalow)      
        deallocate(alphahigh)      
        deallocate(talpha)      
        deallocate(abump)      
        deallocate(aiter)      

        deallocate(rvart)      
        deallocate(ivart)   
        deallocate(fixvar)
        deallocate(maxsafe)
 
        end subroutine createinputfile





      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      implicit none

      REAL(8) A(NP,NP)
      REAL(8) VV(400),D,SUM,AAMAX,DUM
      REAL(8) TINY
      INTEGER INDX(N),N,NP,IMAX,NMAX
      INTEGER I,J,K


      NMAX = 400
      TINY = 1.0d-20 * 1.0d-0
      D=1.0d-0
      VV = 0.0d-0
      SUM = 0.0d-0

      DO 12 I=1,N
        AAMAX=0.0d-0
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.0d-0) then
        continue
!        PAUSE 'Singular matrix.'
        endif
        VV(I)=1.0d-0/AAMAX
12    CONTINUE


      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-1.0d-0 *A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.0d-0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-1.d-0 * D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.0d-0)A(J,J)=TINY
          DUM=1.0d-0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.0d-0)A(N,N)=TINY


      RETURN
      END

!---------------------------------------------------------

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      implicit none

      REAL(8) A(NP,NP)
      REAL(8) B(N),SUM
      INTEGER INDX(N),N,NP,I,II,LL,J,K

      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.0d-0) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
     
 
      RETURN
      END

!----------------------------------
       subroutine dmatinv(a,y,np)

       implicit none
 
       integer np
       real(8)  a(np,np),y(np,np)
       real(8), allocatable :: ac(:,:)
       real(8) d,b(np)
       integer indx(np),i,j,n

 !      external LUDCMP,LUBKSB

       ALLOCATE(ac(np,np))

!      make a copy of a so it will not be affected by 
!      the subroutines and initialize y

       y = 0.0d-0
       loop1: do i = 1,np
         ac(i,i) = a(i,i)
         y(i,i) = 1.0d-0
       loop2: do j = 1,np
         ac(i,j) = a(i,j)
       end do loop2
       end do loop1
       indx(1:np)=0
       n = np

       call LUDCMP(ac,np,np,indx,d)

       loop13: do j = 1,n
          call LUBKSB(ac,n,np,indx,y(1,j))
       end do loop13

       DEALLOCATE(ac)

       return
       end subroutine dmatinv


SUBROUTINE ran1sub(idum,x) 
IMPLICIT NONE
INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL(8) :: x 
INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647,IQ=127773,IR=2836
REAL(8), SAVE :: am
INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
if (idum <= 0 .or. iy < 0) then
     am=nearest(1.0d-0,-1.0d-0)/IM
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1
end if
ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/IQ
iy=IA*(iy-k*IQ)-IR*k
if (iy<0) iy=iy+IM
x=am*ior(iand(IM,ieor(ix,iy)),1)
END SUBROUTINE ran1sub


SUBROUTINE ran1(idum,vv,n)
 
implicit none

! external ran1sub

REAL(8) vv(n),x
INTEGER kk,n,idum

vvloop: do kk = 1,n
     call ran1sub(idum,x)
     vv(kk) = x
end do vvloop

idum = idum+2

END SUBROUTINE ran1
 


!------------------------


       SUBROUTINE GASDEV(idum,vv,n)

! This function returns a normally distributed deviate with
! zero mean and unit variance, using ran1(IDUM) as the
! source of uniform deviates

       implicit none

       integer iset,n,k,gotoc,idum
       real(8) v1,v2,r,fac,gset,gasdev1
       real(8) temp1(2),vv(n),ix1,ix2

!       external ran1

!       write(*,*) 'want',n,'random numbers'

       vvloop: do k = 1,n

       iset = 0
       v1 = 0.0d-0
       v2 = 0.0d-0
       r = 0.0d-0
       fac = 0.0d-0
       gset = 0.0d-0
       gasdev1 = 0.0d-0
       ix1=0.0d-0
       ix2=0.0d-0
       gotoc = 0

       if (iset.eq.0) then


1         call ran1(idum,temp1,2)

          v1 = (2.0d-0)*temp1(1)-1.0d-0 
          v2 = (2.0d-0)*temp1(2)-1.0d-0

!          write(*,*) 'v1',v1
!          write(*,*) 'v2',v2
!          write(*,*) 'idum',idum

          r = v1**2.0d-0+v2**2.0d-0
          if (r.ge.1.0d-0) then 
             gotoc = gotoc + 1
             if (gotoc.gt.n) then
                write(*,*) 'error in gasdev'
             end if
             go to 1
          end if

          fac = dsqrt(-2.0d-0*dlog(r)/r)
          gset = v1*fac
          gasdev1 = v2*fac
          iset = 1
       else
          gasdev1 = gset
          iset = 0
       end if
       
       vv(k)=gasdev1
!       write(*,*) 'random normal number',k,vv(k)

       end do vvloop

       return
       end SUBROUTINE GASDEV





         








 
         








 

!
!      ________________________________________________________
!     |                                                        |
!     |                INVERT A GENERAL MATRIX                 |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         V     --ARRAY CONTAINING MATRIX                |
!     |                                                        |
!     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
!     |                                                        |
!     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
!     |                                                        |
!     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
!     |                      ELEMENTS                          |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         V     --INVERSE                                |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: DABS                             |
!     |________________________________________________________|
!
      SUBROUTINE VERT(V,LV,N,W)

      implicit none

      REAL(8) V(LV,N),S,T
      INTEGER W(N),I,J,K,L,M,N,P,LV


      IF ( N .EQ. 1 ) GOTO 110
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
!     ---------------------------------------
!     |*** FIND PIVOT AND START ROW SWAP ***|
!     ---------------------------------------
      P = L
      IF ( M .GT. N ) GOTO 30
      S = DABS(V(L,L))
      DO 20 I = M,N
           T = DABS(V(I,L))
           IF ( T .LE. S ) GOTO 20
           P = I
           S = T
20    CONTINUE
      W(L) = P
30    S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0.0d-0 ) GOTO 120
!     -----------------------------
!     |*** COMPUTE MULTIPLIERS ***|
!     -----------------------------
      V(L,L) = -1.0d-0
      S = 1.0d-0/S
      DO 40 I = 1,N
40         V(I,L) = -S*V(I,L)
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0.0d-0 ) GOTO 50
!     ------------------------------
!     |*** ELIMINATE BY COLUMNS ***|
!     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO 60 I = 1,K
60         V(I,J) = V(I,J) + T*V(I,L)
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO 80 I = M,N
80         V(I,J) = V(I,J) + T*V(I,L)
      GOTO 50
!     -----------------------
!     |*** PIVOT COLUMNS ***|
!     -----------------------
90    L = W(K)
      DO 100 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,K)
100        V(I,K) = T
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
      RETURN
110   IF ( V(1,1) .EQ. 0.0d-0 ) GOTO 120
      V(1,1) = 1.0d-0/V(1,1)
      RETURN
120   WRITE(6,*) 'ERROR: MATRIX HAS NO INVERSE'
      STOP

      write(*,*) 'exiting vert'
      END
      ! cholesky decomposition routine

	SUBROUTINE choldc (a,n,p)

      implicit none 

      integer,intent(in)::n
      real(8),dimension(n,n),intent(in)::a
      real(8),dimension(n,n),intent(out)::p
      real(8),dimension(n,n)::aa
      real(8) sum
      integer i,j,k


      sum = 0.0d-0
      aa(1:n,1:n)=a(1:n,1:n)
      do 13 i = 1,n
         do 12 j = i,n
         sum=aa(i,j)
	   do 11 k = i-1,1,-1
	      sum = sum - aa(i,k)*aa(j,k)
11	   continue 
         if (i.eq.j) then
 !          if(sum.le.0.0d-0)pause 'choldc failed'
           p(i,i) = dsqrt(sum)
           else
             aa(j,i) = sum/p(i,i)
	       p(j,i) = aa(j,i)
         end if
12       continue
13    continue
      return

      END SUBROUTINE CHOLDC



	! ----------------------- GET_EPS --------------------------------

	SUBROUTINE get_eps(evec,evcv,n,nmc,idum,numper)

!	This subroutine draws the random variables epsilon from a multivariate 
!     random normal distribution with variance-covariance matrix evcv
!     This subroutine calls the cholesky decomposition subroutine

        implicit none

	integer nmc,vl,n,numper
	integer idum,i,j,k,l,s1
	real(8) evec(nmc*numper*n),evcv(n,n),p(n,n),eiid(n)
	real(8) r(nmc*numper*n),etemp(n),rtemp(n)

!	external choldc, ran1

	p = 0.0d-0
	eiid = 0.0d-0
	evec = 0.0d-0
	etemp=0.0d-0


	call choldc(evcv,n,p)
	call gasdev(idum,r,nmc*numper*n)

!     Reconstruct the epsilons

	Loop3: do k = 1,numper
	Loop4: do l = 1,nmc
	   s1 = (k-1)*n*nmc+(l-1)*n+1
	   rtemp=r(s1:s1+n-1)
	   etemp=0.0d-0

	Loop1: do i = 1,n
	Loop2: do j = 1,n
	   etemp(i) = etemp(i)+p(i,j)*rtemp(j)
	end do Loop2
	end do Loop1

	evec(s1:s1+n-1)=etemp

	end do Loop4
	end do Loop3

	return

	end subroutine get_eps



!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!* by straight insertion.                            *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table ARR                  *
!*          ARR	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: Straight insertion is a N² routine and      *
!*       should only be used for relatively small    *
!*       arrays (N<100).                             *
!*****************************************************         
SUBROUTINE PIKSRT(N,ARR)
   integer N,j,i
 real ARR(N),a
  
  do j=2, N
    a=ARR(j)
    do i=j-1,1,-1
      if (ARR(i)<=a) goto 10
      ARR(i+1)=ARR(i)
    end do
	i=0
10  ARR(i+1)=a
  end do
  return
END

!write table of size N to standard output
SUBROUTINE TWRIT(N,ARR)
integer I,N
real ARR(N)
  print *,' '
  WRITE(*,10) (ARR(I),I=1,N)
  return
10 FORMAT(5(' ',F6.1))
END

!end of file sort1.f90
subroutine readtheta2(silent,pathoutput,reform,thetamtype,thetamtype1,thetamtype2,thetamreftype,thetamreftype1,thetamreftype2,numper,numz,ktp)
integer:: silent,reform,group,xx,zz,numper,numz,ktp,ttt
character*200::pathoutput
real(8) thetamref(numper,3,numz),thetamreftype(ktp,numper,3,numz),thetamreftype1(numper,3,numz),thetamreftype2(numper,3,numz)
real(8) thetamtype(ktp,numper,3,numz),thetamtype1(numper,3,numz),thetamtype2(numper,3,numz)

open(unit=400,file=adjustl(trim(pathoutput))//'../inputs/pretheta', status='UNKNOWN',form='FORMATTED', &
        action='READ')
rewind(unit=400)
if (reform.eq.1) then
        open(unit=500,file=adjustl(trim(pathoutput))//'../inputs/posttheta', status='UNKNOWN',form='FORMATTED', &
            action='READ')
        rewind(unit=500)
end if

do ttt=1,ktp
    if (silent.ne.1) write(*,*) 'reading type', ttt
    
    xxloop: do xx = 1,numper
    grouploop: do group = 1,3
    zzloop: do zz = 1,numz
        read(400,*) thetamtype(ttt,xx,group,zz)
        if (reform.eq.1) then
            read(500,*) thetamreftype(ttt,xx,group,zz)
        else
                thetamreftype(ttt,xx,group,zz) = thetamtype(ttt,xx,group,zz)
        end if 
    end do zzloop 
    end do grouploop 
    end do xxloop

enddo

close(unit=400)
close(unit=500)
    33        format(f16.4)
    32        format(f24.8)

endsubroutine
