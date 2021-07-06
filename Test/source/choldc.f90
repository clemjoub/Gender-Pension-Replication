
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
	11 continue 
	if (i.eq.j) then
		!if(sum.le.0.0d-0)pause 'choldc failed'
		p(i,i) = dsqrt(sum)
	else
		aa(j,i) = sum/p(i,i)
		p(j,i) = aa(j,i)
	end if
	12 continue
	13 continue
	return

END SUBROUTINE CHOLDC

! ----------------------- GET_EPS --------------------------------

SUBROUTINE get_eps(evec,evcv,n,nmc,idum,numper)

	!This subroutine draws the random variables epsilon from a multivariate 
	!random normal distribution with variance-covariance matrix evcv
	!This subroutine calls the cholesky decomposition subroutine

	implicit none

	integer nmc,vl,n,numper
	integer idum,i,j,k,l,s1
	real(8) evec(nmc*numper*n),evcv(n,n),p(n,n),eiid(n)
	real(8) r(nmc*numper*n),etemp(n),rtemp(n)

	!external choldc, ran1

	p = 0.0d-0
	eiid = 0.0d-0
	evec = 0.0d-0
	etemp=0.0d-0

	call choldc(evcv,n,p)
	call gasdev(idum,r,nmc*numper*n)

	!Reconstruct the epsilons

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


