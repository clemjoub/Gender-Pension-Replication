
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