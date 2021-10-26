
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



