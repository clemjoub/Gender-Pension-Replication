    
    subroutine MSM(nsim,nmoments,n1,n2,criterion,simdata1,simdata2,simdiv,loss,pathoutput,estimation,counterchar)
 
    implicit none
    integer, parameter:: noutcomes=8
    integer, parameter:: nconditions = 15
    integer, parameter:: numfam = 5314
    integer, parameter:: nyrs = 5
    
    character*6:: spec
    integer i,j,k,l                                                 !counter
    integer nsim                                                    !total number of year-household-clone simulations
    integer nclones                                                 !number of household clones used in the simulation
    integer nhh                                                     !number of households simulated
    integer nmoments                                                !number of moments used in estimation
    integer cohort(nsim)
    real temp(nmoments,2*noutcomes+2*nconditions)
    logical temp6(1,noutcomes), temp5(1,nconditions), temp7(nsim), temp8(nsim)
    real(8) Y(nsim,noutcomes),outcomes(nsim,noutcomes),X(nsim) 
	real upperbound(nsim,nmoments,noutcomes), lowerbound(nsim,nmoments,noutcomes) 
    real(8) conditions(nsim,nconditions)
	real	conditions_ub(nsim,nmoments,nconditions), conditions_lb(nsim,nmoments,nconditions) 
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
	
	
	open(unit=101,action='write',file="1.txt",position='rewind')
	close (unit=101)
	
    do i=1,nsim
	do j=1,noutcomes
		lowerbound(i,:,j)=temp(:,2*j-1)
		upperbound(i,:,j)=temp(:,2*j)
    enddo
	do j=1,nconditions
		conditions_lb(i,:,j)=temp(:,2*noutcomes+2*j-1)
		conditions_ub(i,:,j)=temp(:,2*noutcomes+2*j)
	enddo 
    enddo

	open(unit=102,action='write',file="2.txt",position='rewind')
	close (unit=102)	
	
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
		D(:)=all(conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2)
		D(:)=-1*D(:) 
		contribution(:,k)=merge(-1.d0, 0.d0, all(outcomes(:,:).ge.lowerbound(:,k,:).and.outcomes(:,:).le.upperbound(:,k,:),2).and.all(conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2))
		contribution(:,k)=contribution(:,k)*D(:)

	elseif (meanmom(k).eq.1) then
		Y(:,:)=outcomes(:,:)*lowerbound(:,k,:) ! for moments that are averages (not proportions), lowerbound contains an indicator for the particular outcome relevant to that moment
		X(:)=sum(Y(:,:),2)
		mask(:)=all(conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2).and.X(:).ne.-99.0d-0
		nobs(k)=count(mask) 
		msim(k)=sum(X,mask)/nobs(k)
		var(k)=sum((X-msim(k))**2,mask)/nobs(k)
		D(:)=all(conditions(:,:).ge.conditions_lb(:,k,:).and.conditions(:,:).le.conditions_ub(:,k,:),2)
		D(:)=-1*D(:) 
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
   
