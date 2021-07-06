SUBROUTINE DrawHist(totp,kidmax,ktp,numper,vseed,fixvar,ivart,rvart,maxsafe)

	!This subroutine draws simulated state points at every marriage duration (=period of choice)
	!
	!state variables:
	!total wealth  (wealth)  [-100,00 to +500,000]
	!workexp of male, female - fexp,mexp  
	!education levels of male and female (years of education)
	!number of children (up to kidmax) - n
	!current ages of spuses - hage, wage
	!rates of return: assume fixed so for now we don't draw these
	
    IMPLICIT NONE
	
    !external ran1

	integer totp,kidmax,ktp,numper
	integer i,o,idum,j
	real(8) cut1,cut2,cut3,cut4,r1(1),r2(1),r3(1),r4(1),r5(1),r6(1),r7(1),r20(1),r21(1)
	real(8) r100(1),r101(1),r102(1),r103(1),r104(1),r8(1),r9(1),r10(1),r11(1),r12(1),r13(1),r14(1),r15(1),r16(1), r17(1),r18(1),r19(1)
	real(8) maxsafe(numper)
	integer n(totp,numper), fert(totp,numper)!n is number of kids, fert is current fetility decision
	integer hages(totp),married(totp,numper)
	real(8) safe(totp,numper),vseed(3),risky(totp,numper),pensionw(totp,numper),pensionh(totp,numper)
	integer initmarried
	integer sum1,firstchild, hage_cur(totp,numper)
	integer hed(totp),wed(totp),mardur(totp,numper),mindur,maxdur
	integer workh(totp,numper),workw(totp,numper)
	real(8) rvart(6,totp,numper)   ! row 1 = safe
								   ! row 2 = risky 
								   ! row 3 = rate of return safe
								   ! row 4 = rate of return risky
								   ! row 5 - womens pension wealth
								   ! row 6 - mens pension wealth
	integer ivart(9,totp,numper)   ! row 1 = children (n)
								   ! row 2 = married
								   ! row 3 = work h
								   ! row 4 = work w
								   ! row 5 = whether male surviving
								   ! row 6 = whether female surviving
								   ! row 7 = marriage duration 
								   ! row 8 = formal h
								   ! row 9 = formal w
	integer fixvar(4,totp)         ! row 1 = hed
								   ! row 2 = wed
								   ! row 3 = hage
								   ! row 4 = female (needed for singles for divorced - whether male or female) 
	integer survivh(totp,numper),survivw(totp,numper)
	integer formh(totp,numper), formw(totp,numper),female(totp)

	mardur = 0 
	idum = dint(-1.0d-0*abs(vseed(2)))
	cut1 = 0.2d-0
	cut2 = 2.00d-0
	cut3 = 0.40d-0
	cut4 = 0.85d-0
	n = 0
	survivh = 1
	survivw = 1
	r1 = 0.0d-0
	r2 = 0.0d-0
	r3 = 0.0d-0
	r4 = 0.0d-0
	r5 = 0.0d-0
	r6 = 0.0d-0
	r7 = 0.0d-0
	r9 = 0.0d-0       
	fert = 0
        
	! loop through simulated histories loop (will draw totp histories)
    
    OLoop: do o = 1,totp
        firstchild = -99
        sum1=0
		! Draw gender, which is needed for singles and to know the gender for the person being followed in the case of divorce

	    female(o) = 1
	    if (o.le.int(totp/4).or.o.gt.int(3*totp/4)) then 
		    female(o) = 0
	    endif
	
		! Draw marital status    

        married(o,:) = 1
        !call ran1((idum+6),r3,1)
        if (o.gt.int(totp/2)) then 
            married(o,:) = 0
        end if
        !draw husband age
        if (female(o).eq.1.or.married(o,1).eq.1) then
            call ran1(idum,r5,1)
            hages(o) = dint(r5(1)*10.0d-0+16.0d-0)
        else
            hages(o) = 16   !for men who never marry, make ct index go along with their agem
        end if
        
		! Draw education levels for male and female, ranging from 4 to 20 years
		! men:

        r7 = 0.0d-0
        call ran1((idum+22),r7,1)
	    hed(o)= dint(r7(1)*16+4)

		! women: 
        r8 = 0.0d-0
        call ran1((idum+24),r8,1)
	    wed(o)= dint(r8(1)*16+4)
       
		! Draw fraction of time individual will spend working and workink formally/Informally

	    call ran1((idum+26),r100,1)
	    call ran1((idum+27),r101,1)
	    call ran1((idum+28),r102,1)
	    call ran1((idum+29),r103,1)

		! Draw expected number of children between 0 and 4 and compute the per period probability that will yield that number
    
        call ran1((idum+30),r104,1)
        cut2=int(5.0d-0*r104(1))/25.0d-0      

		! Loop over ages 15 to 90
		!   latest date of preg is 40  (even though model will begin at age 35)
		!   Draw number of kids in every time period, taking into account that in period 1 no kids, in period 2 at most 1... etc.

        Loop1: do i = 1,numper   ! ages 15 through 90 
                     
            hage_cur(o,i)=min0(hages(o)+i-1,90)    ! get mans current age, cap at 90
            
			! randomly draw fertility histories 
     
            if (i.le.25) then   !if woman less than 40
                call ran1((idum+2),r1,1)       
                if (i.eq.1) then
                    sum1 = 0
                end if
                if ((r1(1).lt.cut2).and.sum1.lt.kidmax) then     
	                fert(o,i) = 1
                    ivart(1,o,i) = fert(o,i)
	            end if
                sum1=sum(fert(o,1:(i)))
                n(o,i)=sum1
            else !woman older than 40
                fert(o,i)=0
                ivart(1,o,i)=fert(o,i)
                sum1=sum(fert(o,1:(i-1)))
                n(o,i)=sum1
            end if 

            firstchild = 0 
            if (sum1.eq.0.and.fert(o,i).eq.1) then
                firstchild = i+15   !I set firstchild as age at which woman was pregnant for the first time
            end if

                
			! for some singles draw whether separated/divorced from previous marriage

            call ran1((idum+6),r3,1)

                
            if (married(o,i).eq.0.and.r3(1).lt.0.003d-0) then   !with some prob get divorced/widowed
                mardur(o,i)=int(i/2) !give mariage duration equal to half the number of elapsed periods
            end if
            
            if (married(o,i).eq.1.and.r3(1).lt.0.8) then
                mardur(o,i)=mardur(o,i)+1   !update marriage duration (do not always update to create variation in mariage duration)
            end if


			! draw wealth levels between 1 and 500,000 - using safe assets

            call ran1((idum+10),r9,1)
            call ran1((idum),r10,1)
            safe(o,i) = (r9(1)**4.0d-0)**2.0d-0*dble(i+15)*1000000.0d-0
            risky(o,i)=0.0d-0
            maxsafe(i) = maxval(safe(:,i))
            if (r10(1).lt.0.2d-0) then
                safe(o,i) = 0.0d-0
            end if 

			! draw work pattern - men can only work fulltime or not at all
			! women can work also part time 
			! men:
            
            call ran1((idum+14),r10,1)   
            call ran1((idum+16),r11,1)
            workh(o,i)=0
            workw(o,i)=0
            
            if (r10(1).lt.r100(1)) workh(o,i) = 1                    
            if (r11(1).lt.r101(1)) workw(o,i) = 1      
            
			! draw whether working part-time
			!      	        
			! call ran1((idum+28),r11,1)    
			! if (workw(o,i).eq.1.and.r11(1).lt.0.2d-0) then
			! 	workw(o,i) = 2
			! end if

			! draw whether working in the formal or iformal sector

            call ran1((idum+18),r16,1)
            call ran1((idum+20),r17,1)
            formw(o,i) = 0
            formh(o,i) = 0
            
            if (r16(1).lt.r102(1).and.workh(o,i).eq.1) formh(o,i)=1                    
            if (r17(1).lt.r103(1).and.workw(o,i).eq.1) formw(o,i)=1 
                               
			! draw women and mens pension levels between 1 and 500,000
         
            call ran1((idum+12),r18,1)
            call ran1((idum+22),r19,1)
            call ran1((idum+24),r20,1)
            call ran1((idum+26),r21,1)

            if (o.ge.(totp-2)) then
                r18(1) = 1.0d-0
                r19(1) = 1.0d-0
            end if

            pensionw(o,i) = ((r18(1))**4)*(3000000.0d-0)*sum(formw(o,1:i))
            pensionh(o,i) = ((r18(1))**4)*(3000000.0d-0)*sum(formh(o,1:i))

			! Handle always single case
			! if single, female at ct 35, then set husband pension equal to 0
            
            if (female(o).eq.1.and.married(o,i).eq.0) then
                pensionh(o,i)=0.0d-0
            end if   
                  
			! if single, male at ct 35, then set wife pension equal to 0
        
            if (female(o).eq.0.and.married(o,i).eq.0) then
                pensionw(o,i)=0.0d-0
            end if
            
	    end do Loop1

         
	
	end do OLoop 


    ivart(1,:,:) = fert(:,:)       !this is not number of kids, it is fertility each period. Will use this to infer 
                                   !numkids in subroutine u
    ivart(2,:,:) = married(:,:)    !whether married or not
    ivart(3,:,:) = workh(:,:)      !current work choice husband, will use this to compute experience
    ivart(4,:,:) = workw(:,:)      !current work chocie wife, will use this to compute ecperience
    ivart(5,:,:) = survivh(:,:)
    ivart(6,:,:) = survivw(:,:)
    ivart(7,:,:) = mardur(:,:)
    ivart(8,:,:) = formh(:,:)
    ivart(9,:,:) = formw(:,:)
    fixvar(1,:) = hed(:)
    fixvar(2,:) = wed(:)
    fixvar(3,:) = hages(:)
    fixvar(4,:) = female(:)                     
    rvart(1,:,:) = safe(:,:)
    rvart(2,:,:) = risky(:,:) 
    rvart(3,:,:) = 1.05d-0 
    rvart(4,:,:) = 0.0d-0              !not using risky, so set return equal to 0
    rvart(5,:,:) = pensionw(:,:)
    rvart(6,:,:) = pensionh(:,:)
      
end subroutine DrawHist
