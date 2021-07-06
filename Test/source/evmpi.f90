
!This function will be called by the simplex program
!Given a set of parameters, it will calculate
!the emax coefficients and return the likelihood value

SUBROUTINE EVMPI(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed,fixvar,ivart,rvart, &
	lval,alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc,regind,intmat,numint,nsimul,maxsafe,pathoutput,counter,estimation)

	! specify the number of columns of the theta matrix
	! and the number of time periods in which decisions are
	! being made (from loop runs from age 59 to age 13)

    implicit none

	! external EVFUN

	integer, parameter :: ncoeff=150
	integer, parameter :: nmoments=258        !number of moments used in the estimation

	integer numiter,nsimul,nedraw,estimation
	real(8) loss(1,1:nmoments)
	integer totp,kidmax,ktp,numper,numint,regind(75),intmat(75,75)
	real(8) vseed(3)
	real(8) lval 
	real(8) parms(numiter)
	real(8) alpha(ncoeff,ktp),alphahigh(ncoeff,ktp),alphalow(ncoeff,ktp),abump(ncoeff,ktp)
	integer aiter(ncoeff)
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
	real(8) maxsafe(numper)
	character*200::pathoutput

	
	lval = 0.0d-0
	print*,"pre evfun"
        
    call EVFUN(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed, &
	   fixvar,ivart,rvart, &
	   alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc,regind,numint,intmat,nsimul,&
	   maxsafe,lval,pathoutput, counter,estimation,nmoments,loss)
	print*,"past evfun"
	
end SUBROUTINE EVMPI
