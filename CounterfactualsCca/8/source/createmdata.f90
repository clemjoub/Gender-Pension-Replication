program mdata
 
    implicit none
    integer nsim 		                                             !total number of year-household observations
    integer,parameter:: nmoments = 258					     !number of moments			
    integer,parameter:: n1 = 14   	                                     !number of integer variables in the data
    integer,parameter:: n2 = 12		                                     !number of real variables in the simulated data
    integer,parameter:: estimation = 0
    real(8) temp1(100000,n1)
    real(8) temp2(100000,n2)
    integer, allocatable:: data1(:,:)                                        !array containing the simulated data (integer variables)
    real(8), allocatable:: data2(:,:)           		             !array containing the simulated data (real variables)
    integer,allocatable:: data3(:)					     !unused data array
    real(8) criterion
    integer stat
    integer i,j,k,l                		                             !counters
    character*1::tag,counterchar
    character*200::pathoutput
    real(8):: loss(1,nmoments)
 
    pathoutput="./outputs/"
    counterchar="Z"
    
    ! Read the data
    
    open(file="./inputs/dataformoments13.txt",unit=1)
	i=0
	stat=1
	do while(stat.ge.0)
		i=i+1
		read(1,iostat=stat,fmt=*) temp1(i,:), temp2(i,:)
	enddo
	nsim=i-1
    close(1)
    

    ! Pass data to routine MSM that calculates moments
    
    allocate(data1(nsim,n1))
    allocate(data2(nsim,n2))
    allocate(data3(nsim))
    data1(:,:)=temp1(1:nsim,:)
    data2(:,:)=temp2(1:nsim,:)
    data3(:)=0.0d-0

    call MSM(nsim,nmoments,n1,n2,criterion,data1,data2,data3,loss,pathoutput,estimation,counterchar)

    deallocate(data1)
    deallocate(data2)
    deallocate(data3)

end program
