Program StandardErrors
implicit none

	integer, parameter :: ndropparams=4
	integer, parameter :: nyrs = 5
	integer, parameter :: numfam = 5314
	integer, parameter :: naltbump=23
        integer, parameter :: ncoeff=150
        integer, parameter :: nexcluded=0 !exclude collinear moments
	integer, parameter :: nmoments=258
	integer, parameter :: nparams=98-ndropparams
	integer, parameter :: nparams2=196- 2*ndropparams
        integer, parameter :: maxnumtypes=4
        
	integer mptype,mpt,i,simul,nsimul,jj,j
        integer kidmax,ktp,numper,numint,regind(75),intmat(75,75)
        real(8) lval,varcovar(nparams,nmoments) 
        real(8) params(nparams)
        integer Sbar,Tbar,numz,nmc,xx
        real(8) dpmat(11,3,2)
        integer gradient
        real(8) firstdiff_mat(nmoments,nparams2),gradient_mat(nmoments,nparams)
        real(8) h_mat(nparams),loss(1,nmoments),loss1(nmoments)
	character*10 temp2,temp3
	character*1 bit
	character*3 bit3
	integer dropparams(ndropparams+1)
        real(8) alphaa(ncoeff,maxnumtypes),abumpp(ncoeff,maxnumtypes),rtol,ftol,tausm
        integer aiterr(ncoeff,maxnumtypes)
        integer aiter(ncoeff*maxnumtypes)
        real(8) alphaloww(ncoeff,maxnumtypes),alphahighh(ncoeff,maxnumtypes)
        real(8) alphalow(ncoeff*maxnumtypes),alphahigh(ncoeff*maxnumtypes)
       
        integer nedraw,ncounter
        integer altbump(naltbump+1)
	integer totp
        real(8) vseed(3)
        integer numiter,maxit
        integer dnumdraw,snumdraw
	integer, allocatable :: idx(:)

	character*200 path	
	integer up,k
	real(8) bumpsize(nparams)

        altbump=(/26,27,34,35,42,43,50,51,73,74,75,77,78,79,81,82,83,85,86,87,89,90,91,99999/) !parameters to be bumpped based on predefined range 
	dropparams=(/85,86,87,98,99999/)	!parameters to not consider
	
	!Read parameter ranges
	call Readin(nedraw,totp,kidmax,ktp,numper,vseed,alphaa,abumpp,aiterr,numiter,Sbar,Tbar,nmc,&
          	regind,intmat,numint,nsimul,simul,dnumdraw,snumdraw,tausm,rtol,ftol,maxit,alphaloww,alphahighh)
	
	aiter=reshape(aiterr,(/ncoeff*maxnumtypes/))	
	alphalow=reshape(alphaloww,(/ncoeff*maxnumtypes/))	
	alphahigh=reshape(alphahighh,(/ncoeff*maxnumtypes/))	
		
	idx = PACK([(i, i=1,ncoeff*maxnumtypes)], aiter.eq.1)
	bumpsize(:)=alphahigh(idx)-alphalow(idx)
	
	!Read latest parameters
	jj=1
	open(unit=37,file="./infile.asc",position='rewind')
        	read(37,*) temp3
        	read(37,*) temp2
        
		!read in the initial parameter values 
        	do 100 i = 1,nparams+ndropparams
		if (i.eq.dropparams(jj).and.jj.le.ndropparams) then
			read(37,*)
			jj=jj+1
        	else
			read(37,*) params(i-jj+1)
        		params(i-jj+1) = params(i-jj+1)*1.0d-0
        	endif
       		100 enddo
	 close(37)

		
	j=1 	!counter for the column of firstdiff_mat to read: nparams cols for upward 
		!bumps and nparams for downward bumps. Updated everytime a new moments.txt is read
	do up=0,2,2
        	jj=1 	!counter for excluded parameters. Updated every time the counter kk hits an excluded parameter
        	do k=1,nparams+ndropparams
			if (k.ne.dropparams(jj)) then
				write(bit,'(I1)') up
				write(bit3,'(I3)') k
				path="../Gradient/"//trim(adjustl(bit3))//"/"//adjustl(trim(bit))//"/outputs/"
				open(unit=1,file=adjustl(trim(path))//"moments.txt")
				do i=1,nmoments
					read(1,*) firstdiff_mat(i,j)
				enddo
				j=j+1
	        		close(1)
			elseif (k.eq.dropparams(jj)) then
				jj=jj+1
			endif
		enddo
        enddo
        
	jj=1
        do j=1,nparams
		if (j.eq.altbump(jj).and.jj.le.naltbump) then
			jj=jj+1
        	else
        		bumpsize(j)=abs(params(j)*0.10d-0)
		endif
        	
		do i=1,nmoments
	 		gradient_mat(i,j)=(firstdiff_mat(i,j)-firstdiff_mat(i,(j+nparams)))/(2.0d-0*bumpsize(j))
        	enddo
        enddo
        close(33) 
       
        open(unit=12901, file='gradient.txt',position='rewind')
	write(12901,*) 'gradient_mat'
        do i=1,nmoments
        do j=1,nparams
            write(12901,*) gradient_mat(i,j)
        enddo
        enddo
        close(12901)
       
        write(*,*) 'calling varcov'
        call varcov(params,nmoments,nparams,nexcluded,bumpsize,varcovar)
	call outputstderrors(params, nparams,varcovar,nyrs*numfam,ndropparams,dropparams)
end program
		
	
subroutine varcov(params,nmoments,nparams,nexcluded,bumpsize,varcovar)

implicit none

!****************************************************************************************************************************************************************************
!	Computes varcovar=(1+1/S)*D'WiDi*DtWiVWiD*D'WiDi, where S is number of simulations
!****************************************************************************************************************************************************************************

!INCLUDE 'link_fnl_shared.h'
!use IMSL

!****************************************************************************************************************************************************************************
integer, parameter :: nskip=1																							!Number of lines and colums to skip if want to invert just a submatrix																	
	integer, parameter :: nclones = 10
double precision,allocatable::Vi(:,:),V(:,:),Wi(:,:),D(:,:),Dt(:,:),mat(:,:), mati(:,:), checki(:,:),Dtemp(:,:),DtWiD(:,:),DtWiDi(:,:),DtWiDit(:,:),exclude(:)
double precision,allocatable:: momvar(:), weights(:),DtWiVWiD(:,:),DtWiV(:,:),DtWiVWi(:,:),DtViD(:,:),temp5(:,:),DtWi(:,:),DtVi(:,:),temp6i(:,:),varcovar2(:,:)
double precision :: varcovar(nparams,nmoments)
character(len=40):: blah(nskip)
integer::nexcluded,nmoments,nparams,u,e,i,j
integer flag1,k
real(8) dist1
real(8) params(nparams),tstat(nparams),std(nparams)

real(8) bumpsize(nparams)
!****************************************************************************************************************************************************************************

allocate(DtWiVWiD(nparams,nparams))
allocate(DtWi(nparams,nmoments))
allocate(DtWiV(nparams,nmoments))
allocate(DtWiVWi(nparams,nmoments))
allocate(temp5(nparams,nparams))
allocate(DtVi(nparams,nmoments))
allocate(DtViD(nparams,nmoments))
!allocate(varcovar(nparams,nparams))
allocate(varcovar2(nparams,nparams))

print*, "working with", nmoments, "moments"

! 1. GET V 
!***********
! V is the varcovar of the difference between moments and their expectation. The expectation is approximated by simulating the moments using S draws of the random components of the model

u = nparams   !work on a subset of parameters
write(*,*), "u=" , u
!e=29	!xXX

allocate(mat(nmoments,nmoments))
allocate(Vi(nmoments,nmoments))
allocate(V(nmoments,nmoments))
allocate(checki(nmoments,nmoments))

!! Read in varcovar matrix

open(unit=444, file='./inputs/varcov.txt',position='rewind')

!do i=1,nskip
!	read(fmt=*,unit=444)	! This is if I want to try with a small submatrix first 
!enddo
read(fmt=*,unit=444)
!do i=1,nmoments
do j=1,nmoments
read(fmt=*,unit=444) blah(:),V(:,j)
enddo
!endo

close(444)

!For now use identity matrix
!V(:,:)=0.0d-0
!do j=1,nmoments
!	V(j,j)=1.0d-0
!enddo


! 2. GET D AND Wi
!******************

!! Read in gradient matrix with moments's numerical derivatives wrt paramters (Gradient_mat)
!! Exclude colinear moments and store in Gradient_mat_ex

allocate(D(nmoments,nparams))
allocate(Dt(nparams,nmoments))
allocate(Wi(nmoments,nmoments))
allocate(Dtemp(nmoments+nexcluded,nparams))
allocate(weights(nmoments))
allocate(momvar(nmoments))

!! Read weights from weights.txt

open (unit=777,file='./inputs/mdata.txt')
read(unit=777,fmt=*)
do i=1,nmoments+nexcluded
	read(unit=777,fmt=*) blah, blah,momvar(i) !XXX check that mdata contains variances not inverse variances
	if (momvar(i)>0) weights(i)=1/momvar(i)
enddo
close(777)

!! Read gradient into Gradient_mat array from gradient.txt

open (unit=222,file='gradient.txt')
read(fmt=*,unit=222)

do i=1,nmoments
	Dtemp(i,:)=0
	do j=1,nparams
		read(fmt=*,unit=222) Dtemp(i,j)
	enddo
enddo
close(222)

!! Exclude collinear moments

allocate(exclude(nexcluded+1))
!exclude=(/49,50,51,52,53,54,87,88,89,90,123,124,125,126,167,168,169,170,171,172,185,186,187,188,189,190,199,200,201,202,211,212,213,214,227,228,229,230,231,232,245,246,247,248,249,250,262,263,264,265,266,267,268,273,274,281,282,283,284,285,286,354,372,378,384,396,559,560,561,568,569,570,577,578,579,586,587,588,595,596,597,604,605,606,613,614,615,622,623,624,631,632,633,640,641,642,649,650,651,689,694,699,704,709,714,719,724,740,741,742,743,744,745,782,783,784,785,786,787,800,801,802,803,804,805,818,819,820,821,822,823,836,837,838,839,840,841,852,853,854,855,856,867,868,869,870,871,880,881,882,883,892,893,894,895,902,903,904,911,912,913,100000/)

k=1
D=0.0
Wi=0.0
do i=1,nmoments+nexcluded
	if (i==exclude(k)) then
		k=k+1
	else
		D(i-k+1,:)=Dtemp(i,:)
		Wi(i-k+1,i-k+1)=weights(i)
	endif
enddo

!do i=1,nmoments
!do j=1,nparams
!if (i.ne.j) D(i,j)=0.0
!enddo
!enddo

Dt=transpose(D)



! 3. GET DtWiDi
!***************


!
!V=0.0
!do i=1,nmoments
!V(i,i)=1.0
!enddo
!
!
!Wi=0.0
!do i=1,nmoments
!Wi(i,i)=1.0
!enddo

allocate(DtWiD(nparams,nparams))
allocate(DtWiDi(nparams,nparams))

call matpd(nparams,nmoments,nmoments,Dt,Wi,DtWi)
call matpd(nparams,nmoments,nparams,DtWi,D,DtWiD)


!Output diagonal, to detect large differences across parameters
do i=1,nparams
	write(*,*) i,DtWiD(i,i)
enddo

!Invert DtWiD term
call dmatinv (DtWiD(1:u,1:u),DtWiDi(1:u,1:u),u) 

	deallocate(checki)
	allocate(checki(nparams,nparams))
	call matpd(nparams,nparams,nparams,DtWiD,DtWiDi,checki)
	flag1=0
	CHECKINV2: do i= 1,nparams
			dist1=abs(checki(i,i)-1.0)
			if (dist1.gt.0.0001d-0) then 
				flag1=1
!				pause
			endif
	enddo CHECKINV2
	if (flag1.eq.1) then
		write(*,*)'****problem with the inverse'
		do i=1,nparams
		write(*,*) i,checki(i,i)
		enddo
	endif

	! Denormalize the gradient (XXX double check this)
!	do i=1,nparams
!		Dt(i,:)=Dt(i,:)/bumpsize(i)
!		D(:,i)=D(:,i)/bumpsize(i)
!		do j=1,nparams
!			DtWiDi(i,j)=DtWiDi(i,j)/(bumpsize(i)*bumpsize(j))
!		enddo
!	enddo

!! Multiply matrices to get parameters varcovars

! Build middle term DtWiVWiD

call matpd(nparams,nmoments,nmoments,Dt,Wi,DtWi)
call matpd(nparams,nmoments,nmoments,DtWi,V,DtWiV)
call matpd(nparams,nmoments,nmoments,DtWiV,Wi,DtWiVWi)
call matpd(nparams,nmoments,nparams,DtWiVWi,D,DtWiVWiD)

! Build and output varcovar matrix

call matpd(nparams,nparams,nparams,DtWiDi,DtWiVWiD,temp5)
call matpd(nparams,nparams,nparams,temp5,DtWiDi,varcovar)


! Check: build DtViDi (formula if we had used OWM)

!call dmatinv(V,Vi,nmoments)
!!call DLINRG (nmoments,V, nmoments, Vi, nmoments)
!call matpd(nparams,nparams,nparams,DtViD,varcovar2,checki)
!flag1=0
!CHECKINV3: do i= 1,nparams
!			dist1=abs(checki(i,i)-1.0)
!		if (dist1.gt.0.0001d-0) then 
!				flag1=1
!				pause
!			endif
!enddo CHECKINV3

!if (flag1.eq.1) then
!write(*,*)'****problem with the inverse'
!		do i=1,nparams
!		write(*,*) i,checki(i,i)
!		enddo
!endif
!call matpd(nparams,nmoments,nmoments,Dt,Vi,DtVi)
!call matpd(nparams,nmoments,nparams,DtVi,D,DtViD)

!call dmatinv(DtViD,varcovar2,nparams)
!call DLINRG (nparams,DtViD, nparams, varcovar2, nparams)
!call matpd(nparams,nparams,nparams,DtViD,varcovar2,checki)
!flag1=0
!CHECKINV4: do i= 1,nparams
!			dist1=abs(checki(i,i)-1.0)
!			if (dist1.gt.0.0001d-0) then 
!				flag1=1
!				pause
!		endif
!enddo CHECKINV4

!if (flag1.eq.1) then
!write(*,*)'****problem with the inverse'
!		do i=1,nparams
!		write(*,*) i,checki(i,i)
!		enddo
!endif

!print*, "DtWi", DtWi(1,:)
!print*, "DtWiVWiDt", DtWiVWiD(1,:)

deallocate(mat)
deallocate(Vi)
deallocate(checki)
deallocate(DtWiD)
deallocate(DtWiDi)
deallocate(exclude)
deallocate(DtWiVWiD)
deallocate(DtWi)
deallocate(DtWiV)
deallocate(DtWiVWi)
deallocate(temp5)
deallocate(DtViD)
!deallocate(varcovar)
deallocate(D)
deallocate(Dt)
deallocate(Wi)
deallocate(weights)
deallocate(momvar)




endsubroutine
