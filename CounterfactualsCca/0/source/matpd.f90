

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

    
 
