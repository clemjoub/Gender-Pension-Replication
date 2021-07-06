
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
	!		open(unit=200,file=adjustl(trim(pathoutput))//'output.type1',status='UNKNOWN',action='WRITE')
	!		rewind(unit=200)
	!       if (reform.eq.0) then
	!	    	open(unit=300,file=adjustl(trim(pathoutput))//'pretheta.type1',status='UNKNOWN',action='WRITE')
	!           rewind(unit=300)
	!       else
	!	    	open(unit=300,file=adjustl(trim(pathoutput))//'posttheta.type1',status='UNKNOWN',action='WRITE')
	!           	rewind(unit=300)
	!       end if
	!	end if

	!	if (mptype.eq.2) then 
	!	open(unit=200,file=adjustl(trim(pathoutput))//'output.type2',status='UNKNOWN',action='WRITE')
	!       rewind(unit=200)
	!       if (reform.eq.0) then
	!	    	open(unit=300,file=adjustl(trim(pathoutput))//'pretheta.type2',status='UNKNOWN',action='WRITE')
	!       	rewind(unit=300)
	!       else
	!	    	open(unit=300,file=adjustl(trim(pathoutput))//'posttheta.type2',status='UNKNOWN',action='WRITE')
	!           rewind(unit=300)
	!       end if
	!	end if

	!	if (mptype.eq.3) then 
	!		open(unit=200,file='output.type3',status='UNKNOWN',action='WRITE')
	!       	rewind(unit=200)
	!   	if (reform.eq.0) then
	!		     open(unit=300,file=adjustl(trim(pathoutput))//'pretheta.type3',status='UNKNOWN',action='WRITE')
	!   	        rewind(unit=300)
	!   	else
	!		     open(unit=300,file=adjustl(trim(pathoutput))//'posttheta.type3',status='UNKNOWN',action='WRITE')
	!       	    rewind(unit=300)
	!   	end if
	!	end if

	!
	!
	!	write(200,*) 'total number of histories:',totp
	!	write(200,*) 'number of state variables (max):',numz
	!	write(200,*) 'number histories used in calculating emax', &
	!        nedraw
	!	write(200,*) '--------------------------'
	!   write(200,*) 'Estimated thetas'
	!
	!   Loopth1: do i = (Tbar-15),(35-15),-1  
	!
	!	
	!   	write(200,*) '--------------------------'
	!		write(200,*) 'mothers age',i
	! 
	!         
	!       loopc: do jj = 1,7
	!
	!	    	kk = (jj-1)*4+1
	!
	!	    	if (jj.le.6) then	    
	!	    	write(200,*) '      ',vn(kk),vn(kk+1),vn(kk+2),vn(kk+3)
	!	    	write(200,16) (thetam(i,j),j=kk,kk+3)
	! 	    	write(200,16) (sevm(i,j),j=kk,kk+3)
	!	    	else
	!	    	write(200,*) '      ',vn(kk)
	!	    	write(200,16) (thetam(i,j),j=kk,kk)
	! 	    	write(200,16) (sevm(i,j),j=kk,kk)
	!			end if
	!
	!	    	write(200,*) ' '
	!
	!		end do loopc
	!
	!		if (ind.gt.0) then
	!
	!			write(200,*) 'INTERACTION TERMS'
	!			kk = 1
	!			jj = ncoeff+1
	!			x = numint
	!			do while (x.ge.4)
	!	    	    write(200,21) (intv(j,1),intv(j,2),'     ',j=kk,kk+3)
	!	    	    write(200,16) (thetam(i,j),j=jj,jj+3)
	!	      		write(200,16) (sevm(i,j),j=jj,jj+3)
	!	      		write(200,*) ' '
	!	      		x=x-4
	!	      		jj = jj+4
	!	      		kk=kk+4
	!			end do
	!
	!			if (x.gt.0) then
	! 				write(200,21) (intv(j,1),intv(j,2),'     ',j=kk,kk+x-1)
	!				write(200,16) (thetam(i,j),j=kk,kk+x-1)
	!				write(200,16) (sevm(i,j),j=kk,kk+x-1)
	!				write(200,*) ' '
	!			end if
	!	    
	!		end if
	!
	!		write(200,*) '------------------------ --'
	!       write(200,*) 'R-square, CV R-square'
	!		write(200,18) rsqv(i,1),rsqcvv(i,1)		
	!      	write(200,*) '--------------------------'
	!	     
	!	end do Loopth1
	!
	!	close(unit=200)
	!
	!	write(*,*) 'finished writing output to output files'

	20  format(5a16)
	21  format(i10,i3,a10,i3,i3,a10,i3,i3,a10,i3,i3,a10,i3,i3,a10)
	16  format(4f16.4)
	17  format(3f16.4)
	18  format(f16.4)
	19  format(40f16.4)

	xxloop: do xx = 1,numper
		grouploop: do group =1,3
			zzloop: do zz = 1,numz

				write(300,*) thetam(xx,group,zz)

			end do zzloop 
    	end do grouploop
	end do xxloop
	
    if (mptype.eq.ktp) close(unit=300)

	22  format(f24.8)

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

	16  format(7f8.4)   

    return
	
end subroutine GetTogV


!------------------------ GETTHETA ------------------------------------

subroutine GetTheta(o,ct,Z,emaxv,theta,togv,sevec, Rsq,nedraw,totp,numz,numv,flag1,pathoutput,group,typ,ref,silent)

    ! This subroutine takes the state values and the calculated
    ! emaxv values and returns the thetam parameters. This theta 
    ! values are then used to construct the emax values for the
    ! previous time period

    implicit none

	! external matpd,dmatinv 
    
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
		
    80  continue
	
      c1(i,j)=s1
	  
    90  continue
    91  continue

    return
	
end subroutine matpd

