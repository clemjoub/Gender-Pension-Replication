
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





         








 
         








 
