
subroutine GetData(nivart,nrvart,nfixvar,allstate,idvec,yearmat,numvalid,kidmax,T0,Tbar)

	implicit none

	! XXX clean up so that counter corresponds to age of sampled individual rather than female age ?

	!       Data elements

	! household id
	! year of observation (e.g. 2004 or 2006)
	! number of kids not including current year 
	! whether had a kid this period
	! whether married (0=no, 1 = yes)
	! whether husband working (0 = no, 1 = yes) - code retired as 0
	! whether wife working (0=no, 1 = yes, fulltime,2=year, parttime)
	! husband work experience in years (coming into period, not including current year)
	! wife work experience in years, where part-time work counts as 0.5  (coming into period, not including 2004)
	! indicator for whether husband died
	! indicator for whether wife died
	! indicator for whether divorced
	! marriage duration 
	! whether husband works in formal sector (1 if yes, 0 if no)
	! whether wife works in formal sector  (1 if yes, 0 if no)
	! husband's education in years
	! wife's education in years
	! husband's age in years
	! woman's age in years
	! whether female (for single households or households that become divorced/widowed)
	! nonpension wealth level  (in pesos) - coming into this period
	! pension wealth for woman  (in pesos, can be zero if this is a single man HH) - coming into this period
	! pension wealth for man  (in pesos, can be zero if this is a single woman HH) - coming into this period
	! male earnings
	! female earnings
	! male formal work last period
	! female formal work last period
	! male work last period (0 if no, 1 if yes) 
	! female work last period (0 if no, 1 if fulltime, 2 if parttime)

	integer, parameter :: numfam = 6109 ! number of households
	integer, parameter :: maxl = 12218 ! number of lines in the input data set 
	!integer, parameter :: numfam = 5314 ! number of households
	!integer, parameter :: maxl = 10628 ! number of lines in the input data set 
	integer, parameter :: numper = 75 

	integer id,nobs,year,hed,wed,hage,wage
	integer kidmax
	integer preg,married,singlesamp
	integer idvec(numfam),numvalid
	integer yearmat(numfam,2)  !observe families up to 3 years - need to keep track
	integer i,lastid,j,k,ct
	real(8) nrvart(6,numfam,numper)     ! row 1 = wealth - safe assets
										! row 2 = wealth - risky assets
										! row 3 = rate of return safe
										! row 4 = rate of return risky
										! row 5 = women's pension wealth
										! row 6 = men's pension wealth
	integer nivart(9,numfam,numper)   ! row 1 = children (n)
								   ! row 2 = married
								   ! row 3 = work h
								   ! row 4 = work w
								   ! row 5 = whether male surviving
								   ! row 6 = whether female surviving
								   ! row 7 = marriage duration 
								   ! row 8 = formal h
								   ! row 9 = formal w
	integer nfixvar(4,numfam)      ! row 1 = hed
								   ! row 2 = wed
								   ! row 3 = hage
								   ! row 4 = female (needed for singles for divorced - whether male or female) 
	integer allstate(numfam,numper)  ! keep track of what years the state variables are
									 ! complete (no missing values)
	integer hdead,wdead
	integer numk, hw,ww,formh,formw,survivh,survivf,div,mardur,kidnow,totk
	integer female,lhw,lww,lformh,lformw
	real(8) earnm,earnw
	real(8) hexper,wexper,pensionw,pensionh,wealth,thexper,twexper
	real(8) winformexp,hinformexp,tfhexper,tfwexper
	integer hdeath,wdeath,T0,Tbar
	real(8) hformexp,wformexp
        
	!set the variable that will count the number of lines
	!read in from the input file. Will read until end of file
	!marker
	i = 1

	!initialize variables
	lastid = -1
	j = 0 ! indexes the households
	yearmat = 0
	allstate = 0 
	nivart(1,1:numfam,1:numper) = 0
	nivart(2,1:numfam,1:numper) = 0
	nivart(5,1:numfam,1:numper) = 0
	nivart(6,1:numfam,1:numper) = 0
	nivart(7,1:numfam,1:numper) = 0
	nivart(8,1:numfam,1:numper) = 0
	nivart(9,1:numfam,1:numper) = 0
  
	open(unit=100,file='./inputs/Initialconditions_new.txt',position='rewind')

	do while(i.le.maxl) 

		read(fmt=*,unit=100,end=500) id,year,numk,kidnow,married,hw,ww,hexper,wexper, &
			hformexp,hinformexp,wformexp,winformexp,hdeath,wdeath,hdead,wdead,div,mardur, &
			formh,formw,hed,wed,hage,wage,female, &
			wealth,pensionh,pensionw,earnm,earnw,lformh,lformw,lhw,lww,singlesamp

		!deal with data anomalies
		if (id.eq.77025) then
			female = 0
		endif   
     	
		!recode missing or topcoded variables
		if (wage.eq.-99) then 
				wage = hage  !deal with case of a single male or married couple with missing wife age
			end if
			if (wage.eq.-99.or.(female.eq.0.and.singlesamp.eq.1.and.hage.eq.-99)) then
				write(*,*) 'problem, missing age***'
				pause
		end if
		
		if (married.eq.1.and.wage.gt.hage) hage = wage 	! hard to deal with women who are older than men (but rare)
		if (numk.gt.kidmax) then
			numk=kidmax
		endif

		if (wage.lt.Tbar.and.wage.ge.T0) then 
				if (hed.eq.-99.and.wed.eq.-99) then
					write(*,*) 'id, ed missing for both',id
					pause
				end if 
				if (hed.eq.-99) then
					hed = 0
				end if
				if (wed.eq.-99) then
					wed = 0
				end if
				if (pensionw.eq.-99.0d-0) then
					pensionw = 0.0d-0
				end if
				if (pensionh.eq.-99.0d-0) then
					pensionh = 0.0d-0
				end if
				if (female.eq.2) then 
					female = 1
				else
					female = 0
				end if
			if (married.eq.0) then
					if (hw.eq.-99) then
						hw = 0
					end if
					if (ww.eq.-99) then
						ww = 0
					end if
			end if
			if (hw.eq.0) then
				earnm = 0.0d-0
			end if
			if (ww.eq.0) then
				earnw = 0.0d-0
			end if
			if (kidnow.eq.-99) then
				write(*,*) 'kidnow:',kidnow
			end if
		
			!ct = model period : corresponds to woman's age, or if single male, to man's age
			if (female.eq.1.or.singlesamp.eq.0) then
					ct = wage - 15 
			elseif (female.eq.0.and.singlesamp.eq.1) then
				ct = hage - 15
			endif

			!update fixed household variables if a new id is read
			if(lastid.ne.id) then
				j=j+1
				lastid = id
				yearmat(j,1)=year ! first year in which the household is observed
				yearmat(j,2)=ct	! first model period in which the household is observed
				nfixvar(1,j) = hed
				nfixvar(2,j) = wed
				nfixvar(3,j) = hage - ct + 1 ! husband's age when counter is 1 (ie. the "household was 16" )
				if (hage.eq.-99) then
					nfixvar(3,j) = -99
				end if
				nfixvar(4,j) = female
				idvec(j) = id
			end if

			!determine if all state variables are non missing for that period
			allstate(j,ct) = 1 
			if (married.eq.1.and.(hexper.eq.-99.or.wexper.eq.-99.or.wage.eq.-99.or.hage.eq.-99.or.lww.eq.-99.or.lhw.eq.-99)) then
				allstate(j,ct) = 0 ! married couples
			end if
			if (married.eq.0.and.female.eq.1.and.(wexper.eq.-99.or.wage.eq.-99.or.lww.eq.-99)) then
				allstate(j,ct) = 0 ! single female
			end if
			if (married.eq.0.and.female.eq.0.and.(hexper.eq.-99.or.hage.eq.-99.or.lhw.eq.-99)) then
				allstate(j,ct) = 0 ! single male
			end if
			if (wealth.eq.-99.0d-0) then
				allstate(j,ct) = 0
			end if       
		
			!other reasons for not simulating this household:
			if (wage.ge.Tbar) then
				allstate(j,ct) = 1
			end if
			if (married.eq.0.and.female.eq.0.and.hage.ge.Tbar) then
				allstate(j,ct) = 0 
			end if
	
			!update variable household characteristics	
			nivart(3,j,ct) = hw 
			nivart(4,j,ct) = ww 
			if (ct.gt.1) then      ! lagged work sector
				nivart(3,j,(ct-1)) = kidnow   ! if had kid this period, then preg last period
				nivart(3,j,(ct-1)) = lhw 
				nivart(4,j,(ct-1)) = lww 
			end if
			nivart(2,j,ct) = married 
			nivart(5,j,ct) = 1-hdeath                    
			nivart(6,j,ct) = 1-wdeath
			if (hdeath.eq.1.and.wdeath.eq.1) then 
					write(*,*) 'both husband and wife died'
					write(*,*) 'folio:',id
					allstate(j,ct) = 0 
					pause 
			end if
			if (year.eq.2004.and.ct.lt.numper) then   !set survival and marital status to be same in 2005 as in 2004
				nivart(2,j,(ct+1)) =  nivart(2,j,ct)                    
				nivart(7,j,(ct+1)) =  nivart(7,j,ct)+nivart(2,j,ct)  !marriage duration                    
				nivart(5,j,(ct+1)) = 1-hdeath                    
				nivart(6,j,(ct+1)) = 1-wdeath 
			end if

			!XXXif did not die this period, then need to assign some time
			!whether woman or man died previous period
			!for now I am treating them as being not married

			if (hdeath.eq.0.and.hdead.eq.1) then
					married = 0
			end if
			if (wdeath.eq.0.and.wdead.eq.1) then
					married = 0
			end if
			if (ct.gt.1.and.hdeath.eq.0.and.hdead.eq.1) then
				nivart(5,j,(ct-1)) = 1-hdead             
			end if      
			if (ct.gt.1.and.wdeath.eq.0.and.wdead.eq.1) then
				nivart(6,j,(ct-1)) = 1-wdead             
			end if      
			 
			! fill in some values in work so that experience 
			! will add to the right number
			thexper = hexper
			twexper = wexper
			tfhexper = hformexp
			tfwexper = wformexp
			if (ct.gt.1) then
				exploop: do k = (ct-1),1,-1
					if (thexper.gt.0) then
							nivart(3,j,k) =  1
							thexper = thexper-1.0d-0
					else
							nivart(3,j,k) = 0 
					end if
					if (tfhexper.gt.0) then
							nivart(8,j,k) =  1
							tfhexper = tfhexper-1.0d-0
					else
							nivart(8,j,k) = 0 
					end if
					if (twexper.gt.0) then
							nivart(4,j,k) =  1
							twexper = twexper-1.0d-0
					else
							nivart(4,j,k) = 0 
					end if
					if (tfwexper.gt.0) then
							nivart(9,j,k) =  1
							tfwexper = tfwexper-1.0d-0
					else
							nivart(9,j,k) = 0 
					end if
				end do exploop
			end if
		
			! create fertility history
			totk = 0
			totk = numk-sum(nivart(1,j,1:ct))   ! number of pregnancies woman needs to still have recorded
			
			! in nivart

			kidloop: do k = min(ct-1,25),1,-1   ! starting at ct and going backwards
						if (totk.gt.0.and.nivart(1,j,k).ne.1) then  
								nivart(1,j,k) = 1  ! record a pregnancy
								totk = totk-1      ! move indicator backwards
						end if
			end do kidloop
   
			if (ct.gt.1) then
				nrvart(1,j,(ct-1)) = wealth            ! put all wealth in safe asset,
				nrvart(5,j,(ct-1)) = pensionw            
				nrvart(6,j,(ct-1)) = pensionh            
			end if

			nrvart(3,j,1:numper) = 1.05d-0       ! XXX returns (ensure consistency with rest of code)
			nrvart(4,j,1:numper) = 0.0d-0        ! returns
			
			i = i+1

		end if !wage between T0 and Tbar 

	end do  !end do while loop

	500 close(unit=100)

    numvalid = j

    return

end subroutine GetData
