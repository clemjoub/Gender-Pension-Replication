subroutine createinputfile(parms,numiter,pathoutput,alabel,ktp)

    implicit none

    character*200::pathoutput

	integer, parameter :: maxnumtypes=4
    integer totp,Tbar,Sbar,kidmax
	integer nedraw,nmc,FTbar,vl,ncoeff2
	integer simul,nsimul,ncoeff
	integer dnumdraw,snumdraw,i,ktp
    real(8) vseed(3)
	real(8) disc
	real(8) ftol,rtol,temp1
	integer maxit,temp2
	real(8), allocatable :: alpha(:,:)
	real(8), allocatable :: alphalow(:,:)
	real(8), allocatable :: alphahigh(:,:)
	real(8), allocatable :: talpha(:,:)
	real(8), allocatable :: abump(:,:)
	integer, allocatable :: aiter(:,:)
	integer j,numper,numiter
	real(8) tausm,tempvar
	integer checkp,numint,f1ind,sind,regind(75),intmat(75,75)
	integer numz,jj,kk
	real(8) tparms(numiter)
	real(8) parms(numiter)
	real(8) scale(numiter)
	real(8),allocatable :: rvart(:,:,:)
	integer,allocatable :: ivart(:,:,:)   
	integer, allocatable :: fixvar(:,:)        
	real(8), allocatable :: maxsafe(:)        
	character*50 :: alabel(300)
	character*10 ::  macro(300,maxnumtypes)
 
	ncoeff = 150	

	allocate(alpha(ncoeff,maxnumtypes))      
	allocate(alphalow(ncoeff,maxnumtypes))      
	allocate(alphahigh(ncoeff,maxnumtypes))      
	allocate(talpha(ncoeff,maxnumtypes))      
	allocate(abump(ncoeff,maxnumtypes))      
	allocate(aiter(ncoeff,maxnumtypes))      

	alpha = 0.0d-0
	abump = 0.0d-0

	totp = 0
	kidmax = 0
	numper = 0
	vseed = 0
	numiter = 0
	alpha = 0.0d-0
	abump = 0.0d-0
	aiter = 0
	Tbar = 0
	Sbar = 0

	call Readin(nedraw,totp,kidmax,ktp,numper,vseed,alpha,abump,aiter,numiter,Sbar,Tbar,nmc,&
		 regind,intmat,numint,nsimul,simul,dnumdraw,snumdraw,tausm,rtol,ftol,maxit,alphalow,alphahigh)

	numz = 75+numint

	allocate(rvart(6,totp,numper))      
	allocate(ivart(9,totp,numper))   
	allocate(fixvar(4,totp))
	allocate(maxsafe(numper))

	talpha = alpha

	ncoeff2 = ncoeff-numint

	! update the values of the parameters to reflect latest values

	i = 0
	if (sum(aiter).gt.0) then
		jj1: do jj = 1,ncoeff
			 kk1: do kk = 1,ktp
				if (aiter(jj,kk).eq.1) then
				   i = i + 1
				   talpha(jj,kk)=parms(i)
				end if
			 end do	kk1 !jj
		end do	jj1 !kk
	end if !iterating on values of alpha

	!write out all the ingredients to create a new input file 

	macro(:,:)=""

	alabel(1) = 'alpha(1) - crra coefficient'; macro(1,:)="\crra"
	alabel(2) = 'alpha(2) - stock of children in female marg ut c (nu_0)'
	alabel(3) = 'alpha(3) - stock of children in male marg ut c (nu_0)'
	alabel(4) = 'alpha(4) - leisure female in marg ut c'
	alabel(5) = 'alpha(5) - leisure male in marg ut c'
	alabel(6) = 'alpha(6) - util of leisure female (delta_l)'
	alabel(7) = 'alpha(7) - util of leisure male (delta_l)'
	alabel(8) = 'alpha(8) - util of leisure if work parttime (delta_l)'
	alabel(9) = 'alpha(9) - util of leisure female*nkids*wooman age < 50'
	alabel(10) = 'alpha(10) - util of leisure female*nkids'
	alabel(11) = 'alpha(11) - util of leisure female*married'
	alabel(12) = 'alpha(12) - '
	alabel(13) = 'alpha(13) - man not working age 65+'
	alabel(14) = 'alpha(14) - women not working age 60+'
	alabel(15) = 'alpha(15) - discount rate'
	alabel(16) = 'alpha(16) - ** husbands log income formal - constant term'
	alabel(17) = 'alpha(17) - years of education'
	alabel(18) = 'alpha(18) - experience'
	alabel(19) = 'alpha(19) - experience squared'
	alabel(20) = 'alpha(20) - education squared'
	alabel(21) = 'alpha(21) - ** husbands log income informal - constant'
	alabel(22) = 'alpha(22) - years of education'
	alabel(23) = 'alpha(23) - experience'
	alabel(24) = 'alpha(24) - experience squared'
	alabel(25) = 'alpha(25) - education squared'
	alabel(26) = 'alpha(26) - ** female log income formal - constant'
	alabel(27) = 'alpha(27) - years of education'
	alabel(28) = 'alpha(28) - experience'
	alabel(29) = 'alpha(29) - experience squared'
	alabel(30) = 'alpha(30) - education squared'
	alabel(31) = 'alpha(31) - ** female log income informal - constant'
	alabel(32) = 'alpha(32) - years of education'
	alabel(33) = 'alpha(33) - experience'
	alabel(34) = 'alpha(34) - experience squared'
	alabel(35) = 'alpha(35) - education squared'
	alabel(36) = 'alpha(36) - ** prob covered sector male offer constant'
	alabel(37) = 'alpha(37) - years of schooling'
	alabel(38) = 'alpha(38) - prev period covered'
	alabel(39) = 'alpha(39) - male age'
	alabel(40) = 'alpha(40) - ** prob covered sector female offer constant'
	alabel(41) = 'alpha(41) - years of schooling'
	alabel(42) = 'alpha(42) - prev period covered'
	alabel(43) = 'alpha(43) - female age'
	alabel(44) = 'alpha(44) - ** prob no fert constant'
	alabel(45) = 'alpha(45) - age of woman'
	alabel(46) = 'alpha(46) - married'
	alabel(47) = 'alpha(47) - numk'
	alabel(48) = 'alpha(48) - ** prob no divorce'
	alabel(49) = 'alpha(49) - male age'
	alabel(50) = 'alpha(50) - hed'
	alabel(51) = 'alpha(51) - age diff'
	alabel(52) = 'alpha(52) - duration of marriage'
	alabel(53) = 'alpha(53) - duration squared'
	alabel(54) = 'alpha(54) - fem age'
	alabel(55) = 'alpha(55) - male age'
	alabel(56) = 'alpha(56) - blank'
	alabel(57) = 'alpha(57) - blank'
	alabel(58) = 'alpha(58) - blank'
	alabel(59) = 'alpha(59) - blank'
	alabel(60) = 'alpha(60) - blank'
	alabel(61) = 'alpha(61) - ** more utility parameters'
	alabel(62) = 'alpha(62) - male cost of switching sectors' 
	alabel(63) = 'alpha(63) - female cost of switching sectors'
	alabel(64) = 'alpha(64) - male cost of returning to work'
	alabel(65) = 'alpha(65) - female cost of returning to work'
	alabel(66) = 'alpha(66) - '
	alabel(67) = 'alpha(67) - fertility: woman years ed'
	alabel(68) = 'alpha(68) - fertility: married*kids'
	alabel(69) = 'alpha(69) - '
	alabel(70) = 'alpha(70) - non-pecuniary benefits informal sector (women)'
	alabel(71) = 'alpha(71) - non-pecuniary benefits informal sector (men)'
	alabel(72) = 'alpha(72) - '
	alabel(73) = 'alpha(73) - '
	alabel(74) = 'alpha(74) - '
	alabel(75) = 'alpha(75) - Type prob'
	alabel(76) = 'alpha(76) - wed' 
	alabel(77) = 'alpha(77) - hed '
	alabel(78) = 'alpha(78) - married'
	alabel(79) = 'alpha(79) - cohort'
	alabel(80) = 'alpha(80) - '
	alabel(81) = 'alpha(81) - bequest parameter 1'
	alabel(82) = 'alpha(82) - '
	alabel(83) = 'alpha(83) - bequest parameter 2'
	alabel(84) = 'alpha(84) - '
	alabel(85) = 'alpha(85) - consumption floor'
	alabel(86) = 'alpha(86) - '
	alabel(87) = 'alpha(87) - '
	alabel(88) = 'alpha(88) - '
	alabel(89) = 'alpha(89) - '
	alabel(90) = 'alpha(90) - '
	alabel(91) = 'alpha(91) - '
	alabel(92) = 'alpha(92) - '
	alabel(93) = 'alpha(93) - '
	alabel(94) = 'alpha(94) - '
	alabel(95) = 'alpha(95) - '
	alabel(96) = 'alpha(96) - '
	alabel(97) = 'alpha(97) - '
	alabel(98) = 'alpha(98) -  ***** borrowing constraint'
	alabel(99) = 'alpha(99) - '
	alabel(100) = 'alpha(100) - '
	alabel(101) = 'alpha(101) - '
	alabel(102) = 'alpha(102) - '
	alabel(103) = 'alpha(103) - '
	alabel(104) = 'alpha(104) - mean retun of asset'
	alabel(105) = 'alpha(105) - variance of shock e1: shock to husbands formal income'
	alabel(106) = 'alpha(106) - variance of shock e2: shock to husbands informal income'
	alabel(107) = 'alpha(107) - variance of shock e3: shock to wifes formal income'
	alabel(108) = 'alpha(108) - variance of shock e4: shock to wifes informal income'
	alabel(109) = 'alpha(109) - variance of shock e5: shock to asset'
	alabel(110) = 'alpha(110) - variance of shock e6: shock to hus leisure'
	alabel(111) = 'alpha(111) - variance of shock e5: shock to wife leisure'
	alabel(112) = 'alpha(112) - correlation of shocks e1 and e2'
	alabel(113) = 'alpha(113) - correlation of shocks e1 and e3'
	alabel(114) = 'alpha(114) - correlation of shocks e1 and r4'
	alabel(115) = 'alpha(115) - correlation of shocks e1 and e5'
	alabel(116) = 'alpha(116) - correlation of shocks e1 and e6'
	alabel(117) = 'alpha(117) - correlation of shocks e1 and e7'
	alabel(118) = 'alpha(118) - correlation of shocks e2 and e3'
	alabel(119) = 'alpha(119) - correlation of shocks e2 and e4'
	alabel(120) = 'alpha(120) - correlation of shocks e2 and e5'
	alabel(121) = 'alpha(121) - correlation of shocks e2 and e6'
	alabel(122) = 'alpha(122) - correlation of shocks e2 and e7'
	alabel(123) = 'alpha(123) - correlation of shocks e3 and e4'
	alabel(124) = 'alpha(124) - correlation of shocks e3 and e5'
	alabel(125) = 'alpha(125) - correlation of shocks e3 and e6'
	alabel(126) = 'alpha(126) - correlation of shocks e3 and e7'
	alabel(127) = 'alpha(127) - correlation of shocks e4 and e5'
	alabel(128) = 'alpha(128) - correlation of shocks e4 and e6'
	alabel(129) = 'alpha(129) - correlation of shocks e4 and e7'
	alabel(130) = 'alpha(130) - correlation of shocks e5 and e6'
	alabel(131) = 'alpha(131) - correlation of shocks e5 and e7'
	alabel(132) = 'alpha(132) - correlation of shocks e6 and e7'

	nedraw = 3000

	open(unit=50, file=adjustl(trim(pathoutput))//'parammacro.txt'	,status='UNKNOWN',action='WRITE',position='REWIND')
	open(unit=100,file=adjustl(trim(pathoutput))//'newinputfile.txt',status='UNKNOWN',action='WRITE',position='REWIND')

	write(100,*) totp,'Potential number of states' 
	write(100,*) Tbar,'T-bar' 
	write(100,*) Sbar,'S-bar' 
	write(100,*) kidmax,'Total number of kids possible' 
	write(100,*) nedraw,'Actual number of states used in calculating emax' 
	write(100,*) nmc,'Number of Monte Carlo draws used in emax' 
	write(100,*) vseed(1),'Seed for random epsilon draws' 
	write(100,*) vseed(2),'Seed for generating simulated people' 
	write(100,*) vseed(3),'Seed for other use' 
	write(100,*) ktp,'Number of types' 
	write(100,*) alpha(15,1),'discount rate' 
	write(100,*) 75,'max number of regressors in emax approx' 
	write(100,*) simul,'indicator for whether to do the simulation part' 
	write(100,*) nsimul,'number of people to simulate' 
	write(100,*) dnumdraw,'number of draws in the likelihood estimation part' 
	write(100,*) snumdraw,'number of draws in simulation based on data' 
	write(100,*) tausm,'tau smoothing parameter' 
	write(100,*) 1.0,'contraction factor (simplex parameter)' 
	write(100,*) ftol,'FTOL - simplex parameter' 
	write(100,*) rtol,'RTOL - simplex parameter' 
	write(100,*) maxit,'maximum number of simplex iterations' 

	Loop4: do i=1,132	 							  
	    if (i.eq.1.or.i.eq.6.or.i.eq.7.or.i.eq.16.or.i.eq.21.or.i.eq.26.or.i.eq.31.or.i.eq.70.or.i.eq.71.or.i.eq.75.or.i.eq.76.or.i.eq.77.or.i.eq.78.or.i.eq.79) then  
			Loop7: do j = 1,maxnumtypes	         
				write(100,18) talpha(i,j),alphalow(i,j),alphahigh(i,j),abump(i,j),aiter(i,j),alabel(i)
				write(50,19)"\newcommand{",macro(i,j),"}{",talpha(i,j),"}  %",alabel(i)
			end do Loop7
	    else
	       write(100,18) talpha(i,1),alphalow(i,1),alphahigh(i,1),abump(i,1),aiter(i,1),alabel(i)
	       write(50,19)"\newcommand{",macro(i,1),"}{",talpha(i,1),"}  %",alabel(i)
	   end if	   
	end do Loop4	

    write(100,*) 0.0,'*********************************'	
	write(100,*) 0.0,'*  specify regressors           *'
	write(100,*) 0.0,'*  to include in the emax       *'
	write(100,*) 0.0,'*  calculation (1=yes,0=no)     *'
    write(100,*) 0.0,'*********************************'

	!write out the indicators for whether to include each regressor
	write(100,12) regind(1),'intercept                1'
	write(100,12) regind(2),'                         2'
	write(100,12) regind(3),'                         3'
	write(100,12) regind(4),'                         4'
	write(100,12) regind(5),'                         5'
	write(100,12) regind(6),'                         6'
	write(100,12) regind(7),'                         7'
	write(100,12) regind(8),'                         8'
	write(100,12) regind(9),'                         9'
	write(100,12) regind(10),'                       10'
	write(100,12) regind(11),'                       11'
	write(100,12) regind(12),'                       12'
	write(100,12) regind(13),'                       13'
	write(100,12) regind(14),'                       14'
	write(100,12) regind(15),'                       15'
	write(100,12) regind(16),'                       16'
	write(100,12) regind(17),'                       17'
	write(100,12) regind(18),'                       18'
	write(100,12) regind(19),'                       19'
	write(100,12) regind(20),'                       20'
	write(100,12) regind(21),'                       21'
	write(100,12) regind(22),'                       22'
	write(100,12) regind(23),'                       23'
	write(100,12) regind(24),'                       24'
	write(100,12) regind(25),'                       25'
	write(100,12) regind(26),'                       26'
	write(100,12) regind(27),'                       27'
	write(100,12) regind(28),'                       28'
	write(100,12) regind(29),'                       29'
	write(100,12) regind(30),'                       30'
	write(100,12) regind(31),'                       31'
	write(100,12) regind(32),'                       32'
	write(100,12) regind(33),'                       33'
	write(100,12) regind(34),'                       34'
	write(100,12) regind(35),'                       35'
	write(100,12) regind(36),'                       36'
	write(100,12) regind(37),'                       37'
	write(100,12) regind(38),'                       38'
	write(100,12) regind(39),'                       39'
	write(100,12) regind(40),'                       40'
	write(100,12) regind(41),'                       41'
	write(100,12) regind(42),'                       42'
	write(100,12) regind(43),'                       43'
	write(100,12) regind(44),'                       44'
	write(100,12) regind(45),'                       45' 
	write(100,12) regind(46),'                       46' 
	write(100,12) regind(47),'                       47' 
	write(100,12) regind(48),'                       48' 
	write(100,12) regind(49),'                       49' 
	write(100,12) regind(50),'                       50' 
	write(100,12) regind(51),'                       51' 
	write(100,12) regind(52),'                       52' 
	write(100,12) regind(53),'                       53' 
	write(100,12) regind(54),'                       54' 
	write(100,12) regind(55),'                       55' 
	write(100,12) regind(56),'                       56' 
	write(100,12) regind(57),'                       57' 
	write(100,12) regind(58),'                       58' 
	write(100,12) regind(59),'                       59' 
	write(100,12) regind(60),'                       60' 
	write(100,12) regind(61),'                       61' 
	write(100,12) regind(62),'                       62' 
	write(100,12) regind(63),'                       63' 
	write(100,12) regind(64),'                       64' 
	write(100,12) regind(65),'                       65' 
	write(100,12) regind(66),'                       66' 
	write(100,12) regind(67),'                       67' 
	write(100,12) regind(68),'                       68' 
	write(100,12) regind(69),'                       69' 
	write(100,12) regind(70),'                       70' 
	write(100,12) regind(71),'                       71' 
	write(100,12) regind(72),'                       72' 
	write(100,12) regind(73),'                       73' 
	write(100,12) regind(74),'                       74' 
	write(100,12) regind(75),'                       75' 

	write(100,*) 0.0,'**********************'	
	write(100,*) 0.0,'* interaction        *'
	write(100,*) 0.0,'* terms in emax      *'
	write(100,*) 0.0,'* approx             *'
	write(100,*) 0.0,'**********************'
	
	! read in the location of the interactions
	intloop: do i = 1,75
	  if (intmat(i,i).eq.1) then
		 write(100,11) i,i
	  end if
	  int2loop: do j = (i+1),75
		 if (intmat(i,j).eq.1) then
		   write(100,11) i,j
		 end if
	  end do int2loop
	end do intloop
	
    write(100,11) -9,-9	

	close(unit=100)
	close(unit=50)

	11	format(i4,5x,i4)
	12  format(i1,10x,a50)
	13  format(i1,1x,i1,1x,i1,5x,a50)
	14  format(f20.8,5x,a50)
	15  format(i20,10x,a50)
	16  format(f20.8,1x,f20.8,i4,5x,a50,a6,1x,i1)
	17  format(f20.8,1x,f20.8,i4,5x,a50)
	18  format(f20.8,1x,f20.8,1x,f20.8,1x,f20.8,i4,5x,a50)
	19	format(A12,A10,A2,f10.3,A4,A50)
	
	return
	
	deallocate(alpha)      
	deallocate(alphalow)      
	deallocate(alphahigh)      
	deallocate(talpha)      
	deallocate(abump)      
	deallocate(aiter)      

	deallocate(rvart)      
	deallocate(ivart)   
	deallocate(fixvar)
	deallocate(maxsafe)
 
end subroutine createinputfile



