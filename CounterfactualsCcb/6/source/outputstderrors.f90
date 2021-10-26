
        subroutine outputstderrors(parms,nparms,varcovar,splsize,ndropparams,dropparams)

        implicit none

	integer, parameter :: nclones=10
	character*200::pathoutput
	
	real(8) tstat(nparms),std(nparms), varcovar(nparms,nparms)
	integer, parameter :: maxnumtypes=4
        integer splsize,ndropparams,dropparams(ndropparams),nparms,totp,Tbar,Sbar,kidmax
	integer nedraw,nmc,FTbar,vl,ncoeff2
	integer simul,nsimul,ncoeff,exc
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
        integer jj,kk
        real(8) tparms(nparms)
        real(8) parms(nparms)
        real(8) scale(nparms)
        character*100 :: alabel(300,maxnumtypes)
	character*3 ::stars(nparms)
	integer group(300)
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

	group(1:15)=1
	group(16:43)=2
	group(44:55)=3
	group(62:65)=1
	group(66:69)=3
	group(70:74)=1
	group(80:105)=1
	group(75:79)=4
	group(105:108)=2
	group(110:111)=1

        alabel(1,:) = 'CRRA coefficient                 	& $\sigma$ '
        alabel(2,:) = 'MUc - Stock of children (female) 	& $\nu^f_0$ '
        alabel(3,:) = 'MUc - Stock of children (male)   	& $\nu^m_0$ '
        alabel(4,:) = 'MUc - Leisure (female)           	& $\nu^f_1$ '
        alabel(5,:) = 'MUc - Leisure (male)             	& $\nu^m_1$ '
        alabel(6,1) = 'Ul  - female type 1 				& $\delta^f_I$ '
        alabel(6,2) = 'Ul  - female type 2 				& $\delta^f_I$ '
        alabel(6,3) = 'Ul  - female type 3 				& $\delta^f_I$ '
        alabel(6,4) = 'Ul  - female type 4 				& $\delta^f_I$ '
        alabel(7,1) = 'Ul  - male type 1 	 			& $\delta^m_I$ '
        alabel(7,2) = 'Ul  - male type 2	 			& $\delta^m_I$ '
        alabel(7,3) = 'Ul  - male type 3	 			& $\delta^m_I$ '
        alabel(7,4) = 'Ul  - male type 4	 			& $\delta^m_I$ '
        alabel(8,:) = 'Ul  - part-time				& $\delta^f_p$ '
        alabel(9,:) = 'Ul  - ??? '
        alabel(10,:) = 'Ul - number of children			& $\delta^f_n$ '
        alabel(11,:) = 'Ul - married				& $\delta^f_m$ '
        alabel(12,:) = '??? '
        alabel(13,:) = 'Ul - LRA (male)				& $\delta^m_C$ '
        alabel(14,:) = 'Ul - LRA (female)			& $\delta^f_C$ '
        alabel(15,:) = 'Discount rate				& $\beta$ '
        alabel(16,1) = 'Log formal earnings (male type 1) - constant		& $\theta^m_{0F}$ '
        alabel(16,2) = 'Log formal earnings (male type 2) - constant		& $\theta^m_{0F}$ '
        alabel(16,3) = 'Log formal earnings (male type 3) - constant		& $\theta^m_{0F}$ '
        alabel(16,4) = 'Log formal earnings (male type 4) - constant		& $\theta^m_{0F}$ '
        alabel(17,:) = 'Log formal earnings (male) - schooling		& $\theta^m_{1F}$ '
        alabel(18,:) = 'Log formal earnings (male) - experience		& $\theta^m_{3F}$ '
        alabel(19,:) = 'Log formal earnings (male) - quadratic exp.	& $\theta^m_{4F}$ '
        alabel(20,:) = 'Log formal earnings (male) - quadratic school.	& $\theta^m_{2F}$ '
        alabel(21,1) = 'Log informal earnings (male type 1) - constant		& $\theta^m_{0I}$ '
        alabel(21,2) = 'Log informal earnings (male type 2) - constant		& $\theta^m_{0I}$ '
        alabel(21,3) = 'Log informal earnings (male type 3) - constant		& $\theta^m_{0I}$ '
        alabel(21,4) = 'Log informal earnings (male type 4) - constant		& $\theta^m_{0I}$ '
        alabel(22,:) = 'Log informal earnings (male) - schooling	& $\theta^m_{1I}$ '
        alabel(23,:) = 'Log informal earnings (male) - experience	& $\theta^m_{3I}$ '
        alabel(24,:) = 'Log informal earnings (male) - quadratic exp.	& $\theta^m_{4I}$ '
        alabel(25,:) = 'Log informal earnings (male) - quadratic school.& $\theta^m_{2I}$ '
        alabel(26,1) = 'Log formal earnings (female type 1) - constant		& $\theta^m_{0F}$ '
        alabel(26,2) = 'Log formal earnings (female type 2) - constant		& $\theta^m_{0F}$ '
        alabel(26,3) = 'Log formal earnings (female type 3) - constant		& $\theta^m_{0F}$ '
        alabel(26,4) = 'Log formal earnings (female type 4) - constant		& $\theta^m_{0F}$ '
        alabel(27,:) = 'Log formal earnings (female) - schooling		& $\theta^f_{1F}$ '
        alabel(28,:) = 'Log formal earnings (female) - experience		& $\theta^f_{3F}$ '
        alabel(29,:) = 'Log formal earnings (female) - quadratic exp.		& $\theta^f_{4F}$ '
        alabel(30,:) = 'Log formal earnings (female) - quadratic school.	& $\theta^f_{2F}$ '
        alabel(31,1) = 'Log informal earnings (female type 1) - constant		& $\theta^f_{0I}$ '
        alabel(31,2) = 'Log informal earnings (female type 2) - constant		& $\theta^f_{0I}$ '
        alabel(31,3) = 'Log informal earnings (female type 3) - constant		& $\theta^f_{0I}$ '
        alabel(31,4) = 'Log informal earnings (female type 4) - constant		& $\theta^f_{0I}$ '
        alabel(32,:) = 'Log informal earnings (female) - schooling		& $\theta^f_{1I}$ '
        alabel(33,:) = 'Log informal earnings (female) - experience		& $\theta^f_{3I}$ '
        alabel(34,:) = 'Log informal earnings (female) - quadratic exp.		& $\theta^f_{4I}$ '
        alabel(35,:) = 'Log informal earnings (female) - quadratic school.	& $\theta^f_{2I}$ '
        alabel(36,:) = 'Formal offer probability (male) -  constant		& $\gamma^m_{0}$ '
        alabel(37,:) = 'Formal offer probability (male) - schooling		& $\gamma^m_{2}$ '
        alabel(38,:) = 'Formal offer probability (male) - formal		& $\gamma^m_{1}$ '
        alabel(39,:) = 'Formal offer probability (male) - age			& $\gamma^m_{3}$ '
        alabel(40,:) = 'Formal offer probability (female) -  constant		& $\gamma^f_{0}$ '
        alabel(41,:) = 'Formal offer probability (female) - schooling		& $\gamma^f_{2}$ '
        alabel(42,:) = 'Formal offer probability (female) - formal		& $\gamma^f_{1}$ '
        alabel(43,:) = 'Formal offer probability (female) - age			& $\gamma^f_{3}$ '
        alabel(44,:) = 'alpha(44) - ** prob no fert constant '
        alabel(45,:) = 'alpha(45) - age of woman'
        alabel(46,:) = 'alpha(46) - married'
        alabel(47,:) = 'alpha(47) - numk'
        alabel(48,:) = 'Separation probability - constant			& $\Pi_0$ '
        alabel(49,:) = 'Separation probability - age				& $\Pi_1$ '
        alabel(50,:) = 'Separation probability - schooling			& $\Pi_2$ '
        alabel(51,:) = 'Separation probability - age difference			& $\Pi_3$ '
        alabel(52,:) = 'Separation probability - ???				& $\Pi_4$ '
        alabel(53,:) = 'alpha(53) - duration squared '
        alabel(54,:) = 'alpha(54) - fem age'
        alabel(55,:) = 'alpha(55) - male age'
        alabel(56,:) = 'alpha(56) - blank'
        alabel(57,:) = 'alpha(57) - blank'
        alabel(58,:) = 'alpha(58) - blank'
        alabel(59,:) = 'alpha(59) - blank'
        alabel(60,:) = 'alpha(60) - blank'
        alabel(61,:) = 'alpha(61) - ** more utility parameters'
        alabel(62,:) = 'Switching costs (male)			& $\Phi^m_s$ ' 
        alabel(63,:) = 'Switching costs (female)		& $\Phi^f_s$ ' 
        alabel(64,:) = 'Entry costs (male)			& $\Phi^m_r$ '  
        alabel(65,:) = 'Entry costs (female)			& $\Phi^f_r$ ' 
        alabel(66,:) = 'alpha(66) - '
        alabel(67,:) = 'alpha(67) - fertility: woman years ed'
        alabel(68,:) = 'alpha(68) - fertility: married*kids'
        alabel(69,:) = 'alpha(69) - '
        alabel(70,1) = 'Non-pecuniary benefits informal sector (female type 1)		& $\Phi^f_2$ '
        alabel(70,2) = 'Non-pecuniary benefits informal sector (female type 2)		& $\Phi^f_2$ '
        alabel(70,3) = 'Non-pecuniary benefits informal sector (female type 3)		& $\Phi^f_2$ '
        alabel(70,4) = 'Non-pecuniary benefits informal sector (female type 4)		& $\Phi^f_2$ '
        alabel(71,1) = 'Non-pecuniary benefits informal sector (male type 1)		& $\Phi^m_2$ '
        alabel(71,2) = 'Non-pecuniary benefits informal sector (male type 2)		& $\Phi^m_2$ '
        alabel(71,3) = 'Non-pecuniary benefits informal sector (male type 3)		& $\Phi^m_2$ '
        alabel(71,4) = 'Non-pecuniary benefits informal sector (male type 4)		& $\Phi^m_2$ '
        alabel(72,:) = 'alpha(72) - '
        alabel(73,:) = 'alpha(73) - '
        alabel(74,:) = 'alpha(74) - '
        alabel(75,:) = 'Type logit - constant					& $\rho_0$ '
        alabel(76,:) = 'Type logit - schooling (female)				& $\rho_s^f$ ' 
        alabel(77,:) = 'Type logit - schooling (male)				& $\rho_s^m$ ' 
        alabel(78,:) = 'Type logit - married					& $\rho_m$ ' 
        alabel(79,:) = 'Type logit - cohort					& $\rho_c$ ' 
        alabel(80,:) = 'alpha(80) - '	
        alabel(81,:) = 'alpha(81) - bequest parameter 1'
        alabel(82,:) = 'alpha(82) - '
        alabel(83,:) = 'alpha(83) - bequest parameter 2'
        alabel(84,:) = 'alpha(84) - '
        alabel(85,:) = 'Consumption floor					& $C_{min}$ '
        alabel(86,:) = 'alpha(86) - '
        alabel(87,:) = 'alpha(87) - '
        alabel(88,:) = 'alpha(88) - '
        alabel(89,:) = 'alpha(89) - '
        alabel(90,:) = 'alpha(90) - '
        alabel(91,:) = 'alpha(91) - '
        alabel(92,:) = 'alpha(92) - '
        alabel(93,:) = 'alpha(93) - '
        alabel(94,:) = 'alpha(94) - '
        alabel(95,:) = 'alpha(95) - '
        alabel(96,:) = 'alpha(96) - '
        alabel(97,:) = 'alpha(97) - '
        alabel(98,:) = 'alpha(98) -  ***** borrowing constraint'
        alabel(99,:) = 'alpha(99) - '
        alabel(100,:) = 'alpha(100) - '
        alabel(101,:) = 'alpha(101) - '
        alabel(102,:) = 'alpha(102) - '
        alabel(103,:) = 'alpha(103) - '
        alabel(104,:) = 'alpha(104) - mean retun of asset'
        alabel(105,:) = 'Formal Earnings variance (male)			& $\sigma^m_F$ '
        alabel(106,:) = 'Informal Earnings variance (male)			& $\sigma^m_I$ '
        alabel(107,:) = 'Formal Earnings variance (female)			& $\sigma^f_F$ '
        alabel(108,:) = 'Informal Earnings variance (female)			& $\sigma^f_I$ '
        alabel(109,:) = 'alpha(109) - variance of shock e5: shock to asset'
        alabel(110,:) = 'Ul shock variance (male)				& $\sigma^m_H$ '
        alabel(111,:) = 'Ul shock variance (female)				& $\sigma^f_H$ '
        alabel(112,:) = 'alpha(112) - correlation of shocks e1 and e2'
        alabel(113,:) = 'alpha(113) - correlation of shocks e1 and e3'
        alabel(114,:) = 'alpha(114) - correlation of shocks e1 and r4'
        alabel(115,:) = 'alpha(115) - correlation of shocks e1 and e5'
        alabel(116,:) = 'alpha(116) - correlation of shocks e1 and e6'
        alabel(117,:) = 'alpha(117) - correlation of shocks e1 and e7'
        alabel(118,:) = 'alpha(118) - correlation of shocks e2 and e3'
        alabel(119,:) = 'alpha(119) - correlation of shocks e2 and e4'
        alabel(120,:) = 'alpha(120) - correlation of shocks e2 and e5'
        alabel(121,:) = 'alpha(121) - correlation of shocks e2 and e6'
        alabel(122,:) = 'alpha(122) - correlation of shocks e2 and e7'
        alabel(123,:) = 'alpha(123) - correlation of shocks e3 and e4'
        alabel(124,:) = 'alpha(124) - correlation of shocks e3 and e5'
        alabel(125,:) = 'alpha(125) - correlation of shocks e3 and e6'
        alabel(126,:) = 'alpha(126) - correlation of shocks e3 and e7'
        alabel(127,:) = 'alpha(127) - correlation of shocks e4 and e5'
        alabel(128,:) = 'alpha(128) - correlation of shocks e4 and e6'
        alabel(129,:) = 'alpha(129) - correlation of shocks e4 and e7'
        alabel(130,:) = 'alpha(130) - correlation of shocks e5 and e6'
        alabel(131,:) = 'alpha(131) - correlation of shocks e5 and e7'
        alabel(132,:) = 'alpha(132) - correlation of shocks e6 and e7'

        talpha = alpha

	varcovar=(1+1/nclones)*varcovar  !! Variance upper bound to account for simulations (Gourieroux Monfort)
	
	
	open(unit=2321, file='./outputs/stderrors-preferences.txt',position='rewind')
	open(unit=2322, file='./outputs/stderrors-earnings.txt',position='rewind')
	open(unit=2323, file='./outputs/stderrors-others.txt',position='rewind')
	open(unit=2324, file='./outputs/stderrors-types.txt',position='rewind')

		stars(:)="   "
		i = 0
		if (sum(aiter).gt.0) then
		exc=1
		jj1: do jj = 1,ncoeff
		kk1: do kk = 1,ktp
			if (aiter(jj,kk).eq.1.and.i.lt.nparms) then
		       		if (i.ne.dropparams(exc)) then
					i = i + 1
		       			talpha(jj,kk)=parms(i)
					std(i)=sqrt(varcovar(i,i)/splsize) !!
					tstat(i)=parms(i)/std(i)
					if (abs(tstat(i)).gt.1.645) stars(i)="*  "
					if (abs(tstat(i)).gt.1.960) stars(i)="** "
					if (abs(tstat(i)).gt.2.576) stars(i)="***"
				if (group(jj).eq.1) write(2321,'(A100,A3,E10.2,A3,A10,E10.2,A6)') adjustl(trim(alabel(jj,kk)))," & ", parms(i),stars(i), "& \textit{", std(i),"}  \\"
				if (group(jj).eq.2) write(2322,'(A100,A3,E10.2,A3,A10,E10.2,A6)') adjustl(trim(alabel(jj,kk)))," & ", parms(i),stars(i), "& \textit{", std(i),"}  \\"
				if (group(jj).eq.3) write(2323,'(A100,A3,E10.2,A3,A10,E10.2,A6)') adjustl(trim(alabel(jj,kk)))," & ", parms(i),stars(i), "& \textit{", std(i),"}  \\"
				if (group(jj).eq.4) write(2324,'(A100,A3,E10.2,A3,A10,E10.2,A6)') adjustl(trim(alabel(jj,kk)))," & ", parms(i),stars(i), "& \textit{", std(i),"}  \\"
		    		else if (i.eq.dropparams(exc)) then
					exc=exc+1
				end if
			end if
		end do	kk1 !jj
		end do	jj1 !kk
		end if !iterating on values of alpha

	close(2321)
	close(2322)
	close(2323)
	return
	
        deallocate(alpha)      
        deallocate(alphalow)      
        deallocate(alphahigh)      
        deallocate(talpha)      
        deallocate(abump)      
        deallocate(aiter)      

        end subroutine outputstderrors 



