
SUBROUTINE EVFUN(parms,numiter,nedraw,totp,kidmax,ktp,numper,vseed, &
  fixvar,ivart,rvart,mptype,flik,tprob, &
  alpha,alphalow,alphahigh,abump,aiter,Sbar,Tbar,T0,numz,nmc,dpmat,regind,numint,intmat,nsimul, &
  maxsafe,criterion,pathoutput,counter,estimation,nmoments,loss)

!XXX should more rigorously separate incidental variables (smarried) from history vectors (snivart(2))

implicit none


! external extract_var,getemax,getstate2,maxu,u,get_eps,GetTheta,GetData,ProbDec,ran1,MSM

integer, parameter :: ncoeff=150
integer, parameter :: numfam=6109
!integer, parameter :: numfam=5314
integer, parameter :: nvar1 = 24 
integer, parameter :: nvar2 = 12 
integer, parameter :: lastage = 90
integer, parameter :: lastworkage = 75
integer, parameter :: numa1cat = 10
integer:: nmoments,wageret2,wageret,wageret_ref,wageret_baseline,hageret,ntyp
real(8) loss(1,1:nmoments), fractionsimtype(ktp)
integer flag,counter,nedraw
integer mptype,ct,o,s1,s2,mpt,ctt,skiploop,numint,idum,ss1,ss2,ii,bidule
integer totp,kidmax,ktp,numper,numiter,nmc,regind(75),intmat(75,75),intmatbis(75,75),tt
integer j,wage,numvalid,update,reform,ttt,ftt,tempref, nompg, nopasis
real(8) vseed(3),maxsafe(numper),typeprob(ktp),rr1(1), supp
real(8) deathprob(91,3)
real(8) getemax,actret,cnu(76,4),cnu_baseline(76,4),cnu_ref(76,4)
real(8) r1(1),r2(1),r3(1),r4(1)
integer kk,jj,xx,zz,nosimul,optch(nmc,7),optch2(1,7),jjj,kkk,nsimul
real(8) flik(totp),tprob(totp)
real(8) parms(numiter),probvec(4)
integer snumk,skidnow
integer vv,typ(totp),ftyp,yearmat(numfam,2),thisyear
real(8) alpha(ncoeff,ktp),alphalow(ncoeff,ktp),alphahigh(ncoeff,ktp),abump(ncoeff,ktp),nalpha(ncoeff,ktp)
real(8) supplement,pmin,pmax,temp111(ktp),temp112,temp113
integer aiter(ncoeff,ktp)
real(8) v(7,7),vbump(7,7),rhomat(7,7),newv(7,7),wexper,fexperw,hexper,fexperh
integer viter(7,7),gg
real(8) temprho,dpmat(11,3,2),experh,experw,mpenwithdraw,fpenwithdraw
integer Sbar,Tbar,T0,i,d(4),numz
integer togv(numz),numv,agev(8),zz1,zz2,maxfam
integer sfert,smarried,sfemale,sdiv,shdec,swdec,totp1,totp2,totp3,nedraw1,nedraw2,nedraw3
real(8) thetap(3,numz),thetam(numper,3,numz),maxearnf,maxearnh
real(8) thetamref(numper,3,numz),thetamreftype(ktp,numper,3,numz),thetamreftype1(numper,3,numz),thetamreftype2(numper,3,numz)
real(8) thetamtype(ktp,numper,3,numz),thetamtype1(numper,3,numz),mintypeprobtemp,mintypeprob,maxtypeprob 
real(8) thetamtype2(numper,3,numz), avgtypeprob
real(8) evec(numper*nmc*7),rsqv(numper,1),rsqcvv(numper,1),rsq,rsqcv
real(8) yevec(nmc*7),thetav(1,numz),sevv(1,numz),sevm(numper,numz)
real(8) umaxmc(nmc),Ztemp(1,numz), optc(nmc)
real(8),allocatable :: emaxv(:,:),Z(:,:),Zbis(:,:),emaxvbis(:,:)
real(8) rvart(6,totp,numper)   ! row 1 = wealth -safe (end of period)
                               ! row 2 = wealth - risky 
                               ! row 3 = rate of return safe
                               ! row 4 = rate of return risky 
                               ! row 5 = womens pension wealth
                               ! row 6 = mens pension wealth
integer ivart(9,totp,numper)   ! row 1 = children (n)
                               ! row 2 = married
                               ! row 3 = work h
                               ! row 4 = work w
                               ! row 5 = survivh
                               ! row 6 = survivw
                               ! row 7 = marriage duration
                               ! row 8 = formal h
                               ! row 9 = formal w
integer fixvar(4,totp)         ! row 1 = hed
                               ! row 2 = wed
                               ! row 3 = hage (age difference with wife)
                               ! row 4 = race 
real(8) nrvart(6,1,numper)     ! row 1 = wealth - safe assets (end of period)
                               ! row 2 = wealth - risky assets
                               ! row 3 = rate of return safe
                               ! row 4 = rate of return risky
                               ! row 5 = womens pension 
                               ! row 6 = mens pension wealth 
integer nivart(9,1,numper)      ! row 1 = children (n)
                               ! row 2 = married
                               ! row 3 = work h
                               ! row 4 = work w
                               ! row 5 = survivh
                               ! row 6 = survivw
                               ! row 7 = marriage duration 
                               ! row 8 = formal h
                               ! row 9 = formal w
integer nfixvar(4,1)           ! row 1 = hed
                               ! row 2 = wed
                               ! row 3 = hage
                               ! row 4 = female sampled 
integer sivart(9,numfam,numper),sfixvar(4,numfam),hed,wed
integer snivart(9,1,numper),snfixvar(4,1),tempnk,templsh,templsw,year(numfam),syear
integer allstate(numfam,numper),found,idvec(numfam)
real(8) srvart(6,numfam,numper),sumax(1),snrvart(6,1,numper),hhinc, bono,tao,tao_baseline, tao_ref
real(8) sevec(numper*7),ssevec(7),optsafe(nmc),optrisky(nmc),soptsafe(1),soptrisky(1)
real(8) searnh(1),searnw(1),earnh(nmc),earnw(nmc),simpensionw,simpensionh,simwealth
integer earnhi,earnwi
integer hage,loopstart
real(8) optpenw(nmc),soptpenw(1)
real(8) optpenm(nmc),soptpenm(1)
integer, allocatable :: simdata1(:,:)
real(8), allocatable :: simdata2(:,:)
integer, allocatable :: tsimdata1(:,:)
real(8), allocatable :: tsimdata2(:,:)
integer, allocatable :: simdiv(:),tsimdiv(:),simhdeath(:),simwdeath(:),simtyp(:)
logical, allocatable :: masksimtyp(:)
integer :: clone, ncounter 
integer :: noAPS, nobono, nodivrule, halfPBS,noPBS
character*6:: spec
real(8) :: pasis, mpg,rsqmat(numper,ktp,3)
real(8) criterion                         ! the MSM criterion function returned by subroutine MSM given current parameter vector
integer nsim                              !number of year-household-clone simulated
real(8) tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt
character*200::pathoutput,command,path
character*1::tag
character*2::counterchar
integer:: debug, generous,reformimp,newbest,formalh,formalw,informalh,informalw,married,evermarried,female,numk,silent,hw,ww,n(numper),group,ind
real(8) pensionh,pensionw,hcaph,hcapw,wealth,cmin        
integer::counterfactual,time_array_0(8),time_array_1(8),parabump,estimation,nyrs,highcontrib,highPBS,univPBS,hretage
real::	start_time, finish_time

!---------------------------------------------------------------------------------------------
! Specifications
!---------------------------------------------------------------------------------------------
print*,"evfun"
reform = 0
counterfactual = 0 

silent = 1
skiploop = 1
nosimul = 0
nyrs = 25					!number of calendar years to simulate ahead 
maxfam = 9999999
!nmc = 2
!totp = 300
!nedraw = 300

ncounter = 1
spec = ""
path = ""
pathoutput=adjustl(trim(path))//"outputs/"


! Policy parameters

	tao_baseline = 0.10d-0
	tao_ref = 0.10d-0
	pasis = 565236.0d-0 ! this correspond to the value of PASIS for 70-75 years old in august 2007
	mpg = 2*pasis
	pmin = 12.0d-0*75000.0d-0
	pmax = 12.0d-0*255000.0d-0
	hageret = 65
	wageret_baseline = 60
	wageret_ref = 60
	wageret2 = 65
	bono = tao_ref*1.5d-0*165000.0d-0*12.0d-0
	
	cmin = alpha(85,1)*1000000.0d-0    ! set a minimum consumption level
        actret = alpha(104,1)		   ! interest rate used in the model solution
	nodivrule = 0
	nobono = 0
	noAPS = 0

if (counter.gt.ncounter) goto 88



! Initialize and allocate variables

	write(counterchar,fmt='(I2)') counter
print*,"numfam=",numfam,"nyrs=",nyrs,"nsimul=",nsimul,"nvar1=",nvar1	
	allocate(simdata1((numfam*nyrs*nsimul),nvar1))
print*,"simdata1"	
	allocate(simdata2((numfam*nyrs*nsimul),nvar2))
print*,"simdatar2"	
	allocate(simdiv(numfam*nyrs*nsimul))
print*,"simdiv"	
	allocate(simhdeath(numfam*nyrs*nsimul))
	allocate(simwdeath(numfam*nyrs*nsimul))
print*,"simdeath"	
	allocate(simtyp(numfam*nyrs*nsimul))
print*,"simtyp"	
	allocate(masksimtyp(numfam*nyrs*nsimul))
print*,"allocated sim arrays"
	simtyp = 0 
	simdiv = 0
	flag = 0 
	clone=1
	thisyear=-99
	i = 0
	bidule = 0
	nalpha = alpha
	newv = v
	thetam = 0.0d-0
	thetamtype = 0.0d-0
	thetamreftype = 0.0d-0
	Rsqmat=1.0d-0

	call ProbDec(deathprob,cnu_baseline)
	cnu_ref=cnu_baseline
	if (counterfactual.eq.6.or.counterfactual.eq.10) cnu_ref(:,3)=cnu_baseline(:,2)  
	vv=23

	call extract_var(v,vbump,viter,alpha,abump,aiter,vv,ktp)
	! store interactions in a more compact form, keep original interaction array as "intmatbis", to be used in GetTogV
	intmatbis=intmat 
	ind=0
	intmat=-9
	do i=1,75
	do j=i,75
		if (intmatbis(i,j)==1) then
    			ind=ind+1
    			intmat(ind,1)=i
    			intmat(ind,2)=j
		endif
	enddo 
	enddo 
   
! Adjust parameters for counterfactuals
 
	if (reform.eq.0) then
		tao=tao_baseline
		wageret=wageret_baseline
		cnu=cnu_baseline
	elseif (reform.eq.1) then
		tao=tao_ref
		wageret=wageret_ref
		cnu=cnu_ref
	endif
	If (noAPS.eq.1) pmax = pmin




! Draw shocks

	idum = dint(-1.0d-0*abs(vseed(3)))
	evec = 0.0d-0
	call get_eps(evec,v,7,nmc,idum,numper)
	write(tag,fmt='(I1)') counter
	allocate(emaxv(totp,1))
	allocate(Z(totp,numz))
print*,"allocated emax"


!----------------------------------------------------------------------
! Model solution (stored in thetam)
!----------------------------------------------------------------------


thetam=-99.0d-0 		!This period's emax coefficients
thetap=0.0d-0			!Last period's emax coefficients


if (skiploop.eq.0) then

	if (silent.ne.1) open(unit=981,action='write',file=adjustl(trim(pathoutput))//'R2'//adjustl(trim(counterchar))//".txt",position='rewind')



    	TTTLOOP: do ttt = 1,ktp  ! Loop over types

    	mptype = ttt 			! XXX difference between mptype and ttt??
    
     
        YRLOOP: do i = Tbar,T0,-1 

            if (silent.ne.1) write(*,*) 'working on age',i
            ct = i-15

            ! Initialize variables for Emax regression
            ! if not the last period, then get the thetas that were calculated last period 

            if (i.lt.Tbar) then
               thetap(:,1:numz) = thetam((ct+1),:,1:numz) 
            else
               thetap(:,1:numz)=0.0d-0
            end if

            Z = 0.0d-0
            emaxv = 0.0d-0
            maxearnf = 0.0d-0
            maxearnh = 0.0d-0
 
            FAMLOOP: do o = 1,totp 
            ! if (o.gt.15) cycle famloop


		! Initialize variables

                probvec = 0.0d-0
                mpenwithdraw = 0.0d-0
                fpenwithdraw = 0.0d-0
                umaxmc=0.0d-0
                hage = fixvar(3,o)+ct-1
                nfixvar(1:4,1) = fixvar(1:4,o)				!reduce arrays of state variables to only the current family 
                nivart(1:9,1,1:numper) = ivart(1:9,o,1:numper)
                nrvart(1:6,1,1:numper) = rvart(1:6,o,1:numper)
                optch = 0
                optsafe = 0.0d-0
                optrisky = 0.0d-0
                Ztemp = 0.0d-0


                ! specify subvectors and matrices containing the epsilons
                ! for current history under consideration
                ! will use same draws
                ! for all histories, but different draws every time period
                ! in order to impose iid assumption

                s1 = (nmc*7)*(ct-1)+1	  ! indices for monte carlo draws
                s2 = (s1-1)+nmc*7
                yevec = evec(s1:s2)

                ! call the subroutine that gets the max utility for
                ! all the possible choices and all the epsilon MC draws
                      		
                call maxu(o,ct,clone, thisyear, umaxmc,nmc,nfixvar,nivart,nrvart,yevec,nalpha,thetap,Tbar,numz,ktp,totp,mptype,dpmat,intmat,optch,optsafe,optrisky,earnh,earnw,kidmax,experh,experw,maxearnf,maxearnh,numa1cat,deathprob,cnu,optpenw,optpenm,probvec,maxsafe,mpenwithdraw,fpenwithdraw,reform,actret, nompg, mpg, nopasis, pasis,nobono,bono,nodivrule, pmax,pmin,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt,optc,tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual)

                !  Integrate (get emaxv) by taking the mean over the estimated umax values
                !  for the different monte carlo draws for the epsilons	for the different
                !  types
                
                emaxv(o,1) = sum(umaxmc)/(1.0d-0*nmc)
    
                ! get vector of emax regressors: functions of state variables at beginning of "ct" period
                
                fexperh = 0.0d-0
                fexperw = 0.0d-0
                if (ct.gt.1) then
                   fexperh = sum(nivart(8,1,1:ct))*1.0d-0   !male formal sector work experience
                   fexperw = sum(nivart(9,1,1:ct))*1.0d-0   !woman formal sector work experience
                end if
                numk = sum(nivart(1,1,1:ct))*1.0d-0 
                call getstate2(o,ct,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 

                Z(o,1:numz)=Ztemp(1,1:numz)


            end do FAMLOOP

            !After processed all the histories, run regression to approximate 
            !the emax functions.

            thetav = 0.0d-0      ! vector that will contain estimated thetas
            sevv = 0.0d-0        ! vector that will contain estimated theta standard 
                                 ! errors
            togv = 0   	         ! vector that indicators whether the column of theta
                                 ! is nonzero
            rsq = 0.0d-0	     ! variable containing R-squared from the regression
            rsqcv = 0.0d-0	     ! variable containing R-squared from the cross-val.
                                 ! regression
            
            do group=1,3
           ! figure out which Z variables to turn on or off in the emax
           ! regression

           totp1=int(totp/2)
           totp2=int(totp/4)                 
           totp3=totp-totp2-totp1
           nedraw1=int(totp1*nedraw/totp)
           nedraw2=int(totp2*nedraw/totp)
           nedraw3=nedraw-nedraw1-nedraw2
          

           ! Run the regression, store results in thetav,sevv,rsq,rsqcv
           ! Separate regressions are run for each type (but state variables
           ! will be the same for all types as long as do not include 
           ! type-specific intercepts)
       
           if (group==1) then !married
               allocate(Zbis(totp1,numz))
               allocate(emaxvbis(totp1,1))
               Zbis=Z(1:totp1,:)
               emaxvbis(:,1)=emaxv(1:totp1,1)
               call GetTogV(ct,Zbis,togv,numz,intmatbis,regind,totp1,numint,numv,silent)
  	           call GetTheta(o,ct,Zbis,emaxvbis,thetav(1,:),togv,sevv,rsq,nedraw1,totp1,numz,numv,flag,pathoutput,group,ttt,reform,silent)
          elseif (group==2) then !single female
  	           allocate(Zbis(totp2,numz))
               allocate(emaxvbis(totp1,1))
               Zbis=Z(totp1+1:totp1+totp2,:)
               emaxvbis(:,1)=emaxv(totp1+1:totp1+totp2,1)
               call GetTogV(ct,Zbis,togv,numz,intmatbis,regind,totp2,numint,numv,silent)
  	           call GetTheta(o,ct,Zbis,emaxvbis,thetav(1,:),togv,sevv,rsq,nedraw2,totp2,numz,numv,flag,pathoutput,group,ttt,reform,silent)
           elseif  (group==3) then !single male
  	           allocate(Zbis(totp3,numz))
               allocate(emaxvbis(totp3,1))
               Zbis=Z(totp1+totp2+1:totp,:)
               emaxvbis(:,1)=emaxv(totp1+totp2+1:totp,1)
               call GetTogV(ct,Zbis,togv,numz,intmatbis,regind,totp3,numint,numv,silent)
 	           call GetTheta(o,ct,Zbis,emaxvbis,thetav(1,:),togv,sevv,rsq,nedraw3,totp3,numz,numv,flag,pathoutput,group,ttt,reform,silent)
          endif
       deallocate(Zbis)
       deallocate(emaxvbis)
       if (silent.ne.1) write(*,*)'ct',counter, 'age',(ct+15), 'tp',ttt,'grp',group,'ref',reform
       if (silent.ne.1) write(981,fmt=36) 'age',(ct+15), 'R-square type',ttt,'group',group,'reform',reform,':',rsq
   36 format(A3,I4,A15,I4,A7,I4,A7,I4,A3,f10.7)   
       thetam(ct,group,:) = thetav(1,:) 
       Rsqmat(ct,ttt,group)=rsq
        enddo

        !   sevm(ct,:) = sevv(1,:)
        !   rsqv(ct,1)=rsq
        !   rsqcvv(ct,1)=rsqcv
        
        call date_and_time(values=time_array_1)
        finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60+ time_array_1 (7) + 0.001 * time_array_1 (8)
      !  write (*,*) 'elapsed wall clock time:',finish_time    

 
        end do YRLOOP
       
        ! store the thetas for each type
    
        if (reform.eq.0) then
            thetamtype(ttt,:,:,:) = thetam
	    if (ttt.eq.1) then
               thetamtype1 = thetam
            else
               thetamtype2 = thetam
            end if 
        end if
        if (reform.eq.1) then
            thetamreftype(ttt,:,:,:) = thetam
            if (ttt.eq.1) then
               thetamreftype1 = thetam
            else
               thetamreftype2 = thetam
            end if 
        end if
        if (estimation.ne.1) call ReportRes(thetam,sevm,rsqv,rsqcvv,mptype,Sbar,Tbar,numz,numint,totp,totp,intmatbis,reform,pathoutput,ktp)

    end do TTTLOOP

deallocate(emaxv)
deallocate(Z)

close(981)
end if !skiploop eq 0

!--------------------------------------------------------
! read in the
! thetam matrix from a previous run of the program and do
! not recalculate them
if (estimation.ne.1) call readtheta2(silent,pathoutput,reform,thetamtype,thetamtype1,thetamtype2,thetamreftype,thetamreftype1,thetamreftype2,numper,numz,ktp)

print*,"readtheta"
!----------------------------------------------------------------------
! Simulate sample decisions

!      open (file='checksav.txt', unit=9999)

nmc = 1
update = 1 
mintypeprob=0.5
maxtypeprob=0.5

if (nosimul.eq.0) then

    loopstart = 0
    
    ! write to separate files depending on whether simulating
    ! with or without the reform

    if (silent.eq.0) then
        open(unit=700,file=adjustl(trim(pathoutput))//'simdata'//adjustl(trim(spec))//'.asc',RECL=300,status='REPLACE',form=&
        &'FORMATTED', action='WRITE')
        rewind(unit=700)
    end if

    
    ! initialize idum for drawing the random shocks
    idum = dint(-1.0d-0*abs(vseed(2)))
    
    ! call the subroutine that reads in the data and keeps track of what years have complete state variables
    srvart = 0.0d-0
    sivart = 0 
    sfixvar = 0
    allstate = 0
    numvalid = 0     !number of valid cases in dataset that will be processed
    call GetData(sivart,srvart,sfixvar,allstate,idvec,yearmat,numvalid,kidmax,T0,Tbar) 
    

    s1loop: do s1 = 1,numfam

    if (s1.gt.maxfam) cycle s1loop
    
    if (silent.ne.1) write(*,*) 'family #',s1
         
    ! create family specific state vectors

    snfixvar(1:4,1) = sfixvar(1:4,s1)
    snivart(1:9,1,1:numper) = sivart(1:9,s1,1:numper)  
    snrvart(1:6,1,1:numper) = srvart(1:6,s1,1:numper)     

    cloneloop: do clone=1,nsimul

        ! Find initial year: move forward to a time period (wife's age -15 or single man's age -15) when the state variables are complete
        jj = 1
        do while(allstate(s1,jj).ne.1.and.jj.lt.(Tbar-15)) 
            jj = jj + 1 
        end do  
        ct = jj ! ct is the initial household age -15 
              
        ! Determine the family's permanent or initial characteristics 

        hed = snfixvar(1,1)
        wed = snfixvar(2,1)
        sfemale = snfixvar(4,1)        
        smarried = snivart(2,1,ct)
        if (smarried.eq.0) then
            if (sfemale.eq.1)  snivart(5,1,ct:Tbar-15)=0 !If single, set spouse as dead
            if (sfemale.eq.0)  snivart(6,1,ct:Tbar-15)=0
        endif
        if (smarried.eq.-99) smarried=0
        if (smarried==1) then 
            wage = ct+15
            hage = min0(snfixvar(3,1)+ct-1,90)
        elseif (smarried==0.and.sfemale==0) then
            wage=-99
            hage= ct+15
        elseif (smarried==0.and.sfemale==1) then
            wage=ct+15
            hage=-99
        endif

        ! Compute type probability for this family => type should be included in the fixvar vector XXX
        
        call ran1(idum,rr1,1)
        
	do ftt=1,ktp
        	temp111(ftt) =  alpha(75,ftt)+&
                    		alpha(78,ftt)*dble(smarried)+&
                    		alpha(79,ftt)*dble(hage)*dble(1-sfemale)+&
                    		alpha(79,ftt)*dble(wage)*dble(sfemale)+&
                    		alpha(76,ftt)*dble(wed)*dble(smarried)+ &
                        	alpha(76,ftt)*dble(wed)*dble(sfemale)*(1.0d-0-dble(smarried)) + &
                        	alpha(77,ftt)*dble(hed)*dble(smarried) + &
                    		alpha(77,ftt)*dble(hed)*(1.0d-0-dble(sfemale))*(1.0d-0-dble(smarried))
        	temp111(ftt)= dexp(temp111(ftt))
	enddo
	
	temp112=0.0d-0
	do ftt=1,ktp
		temp112 = temp112 + temp111(ftt)
	enddo
	do ftt=1,ktp
		typeprob(ftt) = temp111(ftt)/temp112
	enddo


        mintypeprobtemp=minval(typeprob(:))
        mintypeprob=dmin1(mintypeprobtemp,mintypeprob)
        
	! Simulate this family's type
        temp113=0.0d-0
	do ftt=1,ktp 
        	temp113=temp113+typeprob(ftt)
		if (rr1(1).le.temp113) then
            		mptype=ftt
        	exit
		endif
	enddo      
        

        ! Draw a lifetime of random shocks for this family 

        sevec = 0.0d-0
        call get_eps(sevec,v,7,1,idum,numper)
        ssevec = 0.0d-0
        
        ! Initialize optimal choice vector
        
        optch2 = -99  
                    !                    optch(kk,1) = ww 
                    !                    optch(kk,2) = hw
                    !                    optch(kk,3) = fert
                    !                    optch(kk,4) = a1
                    !                    optch(kk,5) = a2
                    !                    optch(kk,6) = formh 
                    !                    optch(kk,7) = formw 

        ! Loop over the time periods
        
        ! for each time period when the state space is complete, need to
        ! make predictions until the next time period the state space is complete.

        if (allstate(s1,jj).eq.1) then 

            kkloop: do kk = ct,ct+nyrs-1    ! do 2004-2014. kk is the current hh age -15
                
                thisyear = (kk-yearmat(s1,2)+yearmat(s1,1))      
                               
               ! cycle cohorts that are less than T0 in 2004 
		
		if (ct+15.lt.T0) cycle kkloop

		 ! cycle dead households (i.e. the drawn individual in the household is dead)
                
                if (kk.gt.75) cycle kkloop
                if ((sfemale.eq.1.and.snivart(6,1,kk).eq.0).or.(sfemale.eq.0.and.snivart(5,1,kk).eq.0)) cycle kkloop

                ! only simulate the last cohort (birthyear >1963) after 2019
                
                !if (thisyear-kk-15.lt.1963.and.thisyear.gt.2019) cycle kkloop
                
                ! Set year specific variables:
                
                smarried=snivart(2,1,kk)
                
                !tempref denotes whether the reform is currently in place
                if (thisyear.lt.2009.or.reform.eq.0) then 
                      tempref = 0
                else 
                      tempref = 1 
                end if

             
		if (tempref==1) then
                	select case (thisyear) 
                	 case (2009)
                	    pmin=12.0d-0*67500.0d-0
                	    pmax=12.0d-0*102500.0d-0
                	 case (2010)
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*227500.0d-0
               		  case (2011)
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*175000.0d-0
                	 case (2012,2013,2014,2015,2016)
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*255000.0d-0
                	end select

                	if (thisyear.ge.2016) then
                	    pmin=12.0d-0*75000.0d-0
                	    pmax=12.0d-0*255000.0d-0
                	endif
					If (noAPS.eq.1) pmax = pmin
			wageret=wageret_ref
			tao=tao_ref
			cnu=cnu_ref
		elseif (tempref==0) then
			wageret=wageret_baseline
			tao=tao_baseline
			cnu=cnu_baseline
		endif
                
                if (thisyear.eq.2004) then
                  actret = 0.1055d-0
                else if (thisyear.eq.2005) then
                  actret = 0.0886d-0
                else if (thisyear.eq.2006) then
                  actret = 0.0458d-0
                else if (thisyear.eq.2007) then
                  actret = 0.1577d-0
                else if (thisyear.eq.2008) then
                  actret = -0.1998d-0
                else if (thisyear.eq.2009) then
                  actret = 0.1770d-0
                else if (thisyear.eq.2010) then
                  actret = 0.063d-0-0 
                end if 

                ! get the appropriate thetap vector and optimal choices
                ! if pre 2009 and simulating with reform, using thetam, otherwise use thetamref


                    if (tempref.eq.0) then
                        if ((kk+15).lt.Tbar) then
                            thetap(:,1:numz) = thetamtype(mptype,(kk+1),:,1:numz)
                        else
                            thetap(:,1:numz)=0.0d-0
                        end if  !(kk+15) lt Tbar
                    else      ! use post-reform thetas for emax 
                        if ((kk+15).lt.Tbar) then
                            thetap(:,1:numz) = thetamreftype(mptype,(kk+1),:,1:numz)
                        else
                            thetap(:,1:numz)=0.0d-0
                        end if  !(kk+15) lt Tbar
                    end if


                ! Extract montecarlo draws

                ss1 = 7*(kk-1)+1	  ! indices for monte carlo draws
                ss2 = ss1+6   
                ssevec(1:7)=sevec(ss1:ss2)

                soptsafe = 0.0d-0     ! initialize variables that will store optimal choices of safe and
                soptrisky = 0.0d-0    ! risky
                call maxu(s1,kk,clone,thisyear, sumax,nmc,snfixvar,snivart,snrvart,ssevec,nalpha,thetap,Tbar,numz,ktp,numfam,mptype,dpmat,intmat,optch2,soptsafe,soptrisky,searnh,searnw,kidmax,experh,experw,maxearnf,maxearnh,numa1cat,deathprob,cnu,soptpenw,soptpenm,probvec,maxsafe,mpenwithdraw,fpenwithdraw,tempref,actret, nompg, mpg, nopasis,pasis, nobono,bono, nodivrule, pmax,pmin,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt,pasiscostopt,optc,tao,hageret,wageret,wageret2,lastage,lastworkage,counterfactual)

                ! update the state variables to reflect the optimal choices
                
                if (thisyear.eq.2004) then
                    hexper = 0.0d-0
                    fexperh = 0.0d-0
                    wexper = 0.0d-0
                    fexperw = 0.0d-0
                    ggloop: do gg = 1,(kk-1)
                         hexper = hexper+snivart(3,1,gg)
                         fexperh = fexperh+snivart(8,1,gg)
                         if (snivart(4,1,gg).eq.1) then
                           wexper = wexper+1.0d-0
                           fexperw = fexperw+snivart(9,1,gg)
                         elseif (snivart(4,1,gg).eq.2) then
                           wexper = wexper+0.5d-0
                           fexperw = fexperw+0.5d-0*snivart(9,1,gg)
                         else
                           wexper = wexper+0.0d-0
                           fexperw = fexperw+0.0d-0
                         end if
                    end do ggloop
                 end if
                
                if (sfemale.eq.0.or.smarried.eq.1) then
                    hexper = hexper+optch2(1,2)
                    fexperh = fexperh+optch2(1,2)*optch2(1,6)
                else
                    hexper = -99
                    fexperh= -99
                endif
                
                if (sfemale.eq.1.or.smarried.eq.1) then
                if (optch2(1,1).eq.1) then
                    wexper = wexper+optch2(1,1)
                    fexperw = fexperw+optch2(1,7)*optch2(1,1)
                else if (optch2(1,1).eq.2) then
                    wexper = wexper+0.5d-0
                    fexperw = fexperw+optch2(1,7)*0.5d-0
                else
                    wexper = wexper+0.0d-0
                    fexperw = fexperw+0.0d-0
                end if
                else
                    wexper = -99
                    fexperw = -99
                endif 

                snivart(3,1,kk)  = optch2(1,2)   !work status of husband
                snivart(4,1,kk)  = optch2(1,1)   !work status of wife
                snivart(8,1,kk)  = optch2(1,6)   !formal work status of husband
                snivart(9,1,kk)  = optch2(1,7)   !formal work status of wife
                
                ! set next period as having a complete state space (since simulating many years forward)
                ! XXX figure this out better
 
                if ((kk+1).lt.65) then
                    allstate(s1,(kk+1)) = 1
                end if 

                simpensionw = soptpenw(1) !end of pension balance at the optimal choice
                simpensionh = soptpenm(1)
                simwealth = soptsafe(1) 

                ! Simulate stochastic transitions
                
		if ((kk+1).le.Tbar-15) then 
		if (sfemale.eq.1.and.snivart(5,1,kk).eq.1) snivart(6,1,kk+1)=1 ! by default, sampled individual stays alive until next period
		if (sfemale.eq.0.and.snivart(6,1,kk).eq.1) snivart(5,1,kk+1)=1 
		endif

		! Absent divorce or death the household structure is preserved
		if (smarried.eq.1.and.(kk+1).le.(Tbar-15)) then
                    snivart(2,1,(kk+1))=1 ! still married
		    snivart(5,1,(kk+1))=1 ! still live husband		
		    snivart(6,1,(kk+1))=1 ! still live wife 		
                end if
                
                ! simulate divorce
                
		sdiv = 0
                if (smarried.eq.1.and.wage.lt.60) then !only allow divorce until wife turns 60
                    r1=0.0d-0
		    call ran1(idum,r1,1)
                    if (r1(1).lt.probvec(1)+0.02) then !adjust in simulation for the fact that divorce rates increased (unanticipated)
                         sdiv = 1
                         snivart(2,1,(kk+1):(Tbar-15))=0 ! if divorce, household becomes unmarried
                         if (sfemale.eq.1) snivart(5,1,(kk+1):(Tbar-15))=0 !if divorce and female sampled, hh no longer contains live husband
                         if (sfemale.eq.0) snivart(6,1,(kk+1):(Tbar-15))=0 !if divorce and male sampled, hh no longer contains live wife
                         simwealth = 0.5d-0*simwealth   !get half of private assets
                         if (tempref.eq.1.and.nodivrule.eq.0) then !if divrule applies
                            simpensionh=(simpensionw+simpensionh)/2.0d-0  !after reform pension is split between spouses
                            simpensionw=(simpensionw+simpensionh)/2.0d-0  !after reform pension is split between spouses
                         end if
                    end if !r1
                end if !smarried

                ! Simulate whether a member of the household dies
                
                shdec = 0
                swdec = 0
                r1=0.0d-0 
                call ran1(idum,r1,1)
                r2=0.0d-0 
                call ran1(idum,r2,1)

                ! Husband dies
                if (snivart(5,1,kk).eq.1) then ! alive male in the household
                if (r1(1).lt.probvec(3).or.kk.eq.(Tbar-15)) then ! male dies
                    shdec = 1
                    snivart(5,1,(kk+1):Tbar-15)=0   !male will be dead
                    snivart(2,1,(kk+1):Tbar-15)=0   !household will be unmarried
                    if (smarried.eq.1) then
                        simpensionw = simpensionw+0.6d-0*simpensionh  !if husband dies, wife inherits 0.6 of his pension
                        simpensionh = 0.0d-0
                    endif
                else !r1
                    shdec = 0 
                    snivart(5,1,kk+1) = 1
                endif
                endif
                
		! Wife dies
                if (snivart(6,1,kk).eq.1) then !alive female in the household
                if (r2(1).lt.probvec(4).or. kk.eq.(Tbar-15)) then
                    swdec = 1
                    snivart(6,1,(kk+1):Tbar-15)=0   !female will be dead
                    snivart(2,1,(kk+1):Tbar-15)=0   !household will be unmarried
                    if (smarried.eq.1) then                             
                        simpensionh = simpensionh  !if wife dies, husband does not inherit her pension (gets own)
                        simpensionw = 0.0d-0
                    endif
                else
                    swdec = 0 
                    snivart(6,1,kk+1) = 1
                end if  !r2
                endif

                !draw fertility 
                snumk=sum(snivart(1,1,1:kk))
                skidnow = 0   ! simulated fertility
                if ((sfemale.eq.1.or.smarried.eq.1).and.kk.le.25.and.snumk.lt.kidmax) then
                   r1=0.0d-0
                   call ran1(idum,r1,1)
                   if (r1(1).lt.probvec(2)) then
                        snivart(1,1,(kk+1)) = 1
                        skidnow = 1 
                   else 
                        snivart(1,1,(kk+1)) = 0
                        skidnow = 0 
                   end if  !r1
                end if  !sfemale 


                ! Give out bono por hijo if eligible
                supplement = 0.0d-0
                if (reform.eq.1.and.thisyear.eq.2008.and.nobono.eq.0) then
                    supplement = bono*dble(snumk)   !pension reform fertility supplement
                end if            
                if (reform.eq.1.and.thisyear.gt.2008.and.nobono.eq.0) then
                    supplement = bono*dble(skidnow)
                end if 

                simpensionw = simpensionw + supplement
                
                ! fill in values of simdata1 and simdata2 and output simulation to text file
                if (loopstart.eq.0) then
                
                	simdata1=-99
                	simdata2=-99
                	simdiv=-99
                	simtyp=-99
                	simhdeath=-99
                	simwdeath=-99
     
                	nsim = 1
                	!     write(700,*) 'idvec thisyear married female wage hage shw sww sformh sformw savrate simpensionw simpensionh simwealth numk kidnow hexper fexperh wexper fexperw earnh earnw mpenwithdraw fpenwithdraw wed hed sdiv tax_paidopt spsmcostopt spsfcostopt fmpgcostopt mmpgcostopt pasiscostopt'
                end if

                simwealth=dble(int(simwealth)*1000.0d-0)/1000.0d-0
                optc(1)=optc(1)/1000000
                supp=supplement/1000.0d-0
                                          
                simdata1(nsim,1) = idvec(s1)
                simdata1(nsim,2) = thisyear
                simdata1(nsim,3) = smarried 
                simdata1(nsim,4) = sfemale 
                simdata1(nsim,5) = wage
                simdata1(nsim,6) = hage
                simdata1(nsim,7) = optch2(1,2)
                simdata1(nsim,8) = optch2(1,1)
                simdata1(nsim,9) = optch2(1,6)
                simdata1(nsim,10) = optch2(1,7)
                simdata2(nsim,1) = optch2(1,4)*1.0d-0 
                simdata2(nsim,2) = simpensionw
                simdata2(nsim,3) = simpensionh
                simdata2(nsim,4) = simwealth
                simdata1(nsim,11) = snumk
                simdata1(nsim,12) = skidnow
                simdata2(nsim,5) = hexper 
                simdata2(nsim,6) =  fexperh
                simdata2(nsim,7) =  wexper
                simdata2(nsim,8) =  fexperw
                simdata2(nsim,9) =  searnh(1)
                simdata2(nsim,10) = searnw(1)
                simdata2(nsim,11) =  mpenwithdraw
                simdata2(nsim,12) =  fpenwithdraw
                simdata1(nsim,13) = hed
                simdata1(nsim,14) = wed
                simdiv(nsim) = sdiv
                simtyp(nsim) = mptype
		simhdeath(nsim) = shdec
		simwdeath(nsim) = swdec
                

                if (silent.eq.0.and.nsim>2) write(700,fmt=26) idvec(s1),thisyear,clone, mptype, smarried,sfemale, wage,hage,optch2(1,2),optch2(1,1),optch2(1,6), &
                  optch2(1,7),optch2(1,4),simpensionw,simpensionh,simwealth,snumk,skidnow,hexper,fexperh,wexper,fexperw, &
                searnh(1),searnw(1),mpenwithdraw,fpenwithdraw,hed,wed,sdiv,tax_paidopt,spsmcostopt,spsfcostopt,fmpgcostopt,mmpgcostopt, &
                  pasiscostopt, supp, sumax(1), optc(1),swdec,shdec ,snivart(5,1,kk),snivart(6,1,kk)

26      format(i7,i6, i3, i3, i5,i5,i5,i5,i2,i2,i2,i2,i4,3f14.2,i5,i5,4f8.2,4f14.2,i5,i5,i5,7f9.2,f16.2, f8.2)

                ! Update ages, husband stops aging at 90
                
                if (sfemale.eq.1.or.smarried.eq.1) wage = kk+1+15
                if (sfemale.eq.0.or.smarried.eq.1) hage = min0(snfixvar(3,1)+kk-1+1,90)
                
		! Update state variables XXX eventually do all the state variables here and only use hte temp variables elsewhere in the kkloop
		! Will create halive and walive for snivart(5 and snivart(6
	
                snrvart(5,1,kk) = simpensionw 
                snrvart(6,1,kk) = simpensionh
                snrvart(1,1,kk) = simwealth
		

		loopstart = 1
                nsim = nsim + 1

            end do kkloop

        end if  !allstate 
        end do cloneloop
        
   ! Update full state variable vector using family specific state vectors

   sfixvar(1:4,s1) = snfixvar(1:4,1)
   sivart(1:9,s1,1:numper) = snivart(1:9,1,1:numper)  
   srvart(1:6,s1,1:numper) = snrvart(1:6,1,1:numper)
   
   end do s1loop
        

   ! Rule out near-degenerate type mixtures

   if (mintypeprob.gt.0.95.or.maxtypeprob.lt.0.05) flag=1
   do ttt=1,ktp
	do i=1,nsim
 	masksimtyp(i)=simtyp(i).eq.ttt
	enddo
	ntyp=count(masksimtyp)
	fractionsimtype(ttt)=(ntyp*1.0d-0)/((nsim-1)*1.0d-0)
   	if (silent.eq.0) write(*,*)"Fraction of type",ttt,":",fractionsimtype(ttt)
	if (fractionsimtype(ttt).lt.0.1.or.fractionsimtype(ttt).gt.1.9) flag=1
   enddo
   !avgtypeprob=sum(simtyp(1:(nsim-1)))*1.0d-0
   !avgtypeprob=avgtypeprob/(nsim-1)*1.0d-0
   !if (avgtypeprob.lt.1.1.or.avgtypeprob.gt.1.9) flag=1     
     

close(unit=700)

nsim=nsim-1

allocate(tsimdata1(nsim,nvar1))
allocate(tsimdata2(nsim,nvar2))
allocate(tsimdiv(nsim))
print*,"allocated tsimdata"
tsimdata1(1:nsim,1:nvar1) = simdata1(1:nsim,1:nvar1)
tsimdata2(1:nsim,1:nvar2) = simdata2(1:nsim,1:nvar2)
tsimdiv(1:nsim)=simdiv(1:nsim)

tsimdata2(:,4)=tsimdata2(:,4)/1000000
tsimdata2(:,9)=tsimdata2(:,9)/1000000
tsimdata2(:,10)=tsimdata2(:,10)/1000000

call MSM(nsim,nmoments,nvar1,nvar2,criterion,tsimdata1,tsimdata2,tsimdiv,loss,pathoutput,estimation,counterchar)
!if (silent.ne.1) print*,"criterion", criterion
!if (silent.ne.1.or.estimation.eq.1) print*,minval(Rsqmat(20:75,:,:))
if (flag.eq.1) criterion=9999999999999

deallocate(tsimdata1)
deallocate(tsimdata2)
deallocate(tsimdiv)
end if  !simul eq 1


deallocate(simdata2)
deallocate(simdiv)
deallocate(simhdeath)
deallocate(simwdeath)
deallocate(simdata1)
deallocate(simtyp)
deallocate(masksimtyp)
print*,"end of evfun"


88      end subroutine EVFUN


!------------------------------------------------------------------------------

subroutine extract_var(v,vbump,viter,alpha,abump,aiter,vv,ktp)

integer, parameter :: ncoeff=150

integer ktp
real(8) v(7,7),vbump(7,7),alpha(ncoeff,ktp)
real(8) rhomat(7,7),abump(ncoeff,ktp)
integer vv, i,aiter(ncoeff,ktp),viter(7,7),j



vv = 105 
v(1,1) = alpha(vv,1)
v(2,2) = alpha(vv+1,1)
v(3,3) = alpha(vv+2,1)
v(4,4) = alpha(vv+3,1)
v(5,5) = alpha(vv+4,1)
v(6,6) = alpha(vv+5,1)
v(7,7) = alpha(vv+6,1)
vbump(1,1) = abump(vv,1)
vbump(2,2) = abump(vv+1,1)
vbump(3,3) = abump(vv+2,1)
vbump(4,4) = abump(vv+3,1)
vbump(5,5) = abump(vv+4,1)
vbump(6,6) = abump(vv+5,1)
vbump(7,7) = abump(vv+6,1)
viter(1,1) = aiter(vv,1)
viter(2,1) = aiter(vv+1,1)
viter(3,1) = aiter(vv+2,1)
viter(4,1) = aiter(vv+3,1)
viter(5,1) = aiter(vv+4,1)
viter(6,1) = aiter(vv+5,1)
viter(7,1) = aiter(vv+6,1)
rhomat(1,2) = alpha(vv+7,1)
rhomat(1,3) = alpha(vv+8,1)
rhomat(1,4) = alpha(vv+9,1)
rhomat(1,5) = alpha(vv+10,1)
rhomat(1,6) = alpha(vv+11,1)
rhomat(1,7) = alpha(vv+12,1)
rhomat(2,3) = alpha(vv+13,1)
rhomat(2,4) = alpha(vv+14,1)
rhomat(2,5) = alpha(vv+15,1)
rhomat(2,6) = alpha(vv+16,1)
rhomat(2,7) = alpha(vv+17,1)
rhomat(3,4) = alpha(vv+18,1)
rhomat(3,5) = alpha(vv+19,1)
rhomat(3,6) = alpha(vv+20,1)
rhomat(3,7) = alpha(vv+21,1)
rhomat(4,5) = alpha(vv+22,1)
rhomat(4,6) = alpha(vv+23,1)
rhomat(4,7) = alpha(vv+24,1)
rhomat(5,6) = alpha(vv+25,1)
rhomat(5,7) = alpha(vv+26,1)
rhomat(6,7) = alpha(vv+27,1)
viter(1,2) = aiter(vv+7,1)
viter(1,3) = aiter(vv+8,1)
viter(1,4) = aiter(vv+9,1)
viter(1,5) = aiter(vv+10,1)
viter(1,6) = aiter(vv+11,1)
viter(1,7) = aiter(vv+12,1)
viter(2,3) = aiter(vv+13,1)
viter(2,4) = aiter(vv+14,1)
viter(2,5) = aiter(vv+15,1)
viter(2,6) = aiter(vv+16,1)
viter(2,7) = aiter(vv+17,1)
viter(3,4) = aiter(vv+18,1)
viter(3,5) = aiter(vv+19,1)
viter(3,6) = aiter(vv+20,1)
viter(3,7) = aiter(vv+21,1)
viter(4,5) = aiter(vv+22,1)
viter(4,6) = aiter(vv+23,1)
viter(4,7) = aiter(vv+24,1)
viter(5,6) = aiter(vv+25,1)
viter(5,7) = aiter(vv+26,1)
viter(6,7) = aiter(vv+27,1)

v(1,2) = rhomat(1,2)*dsqrt(v(1,1))*dsqrt(v(2,2))
v(1,3) = rhomat(1,3)*dsqrt(v(1,1))*dsqrt(v(3,3))
v(1,4) = rhomat(1,4)*dsqrt(v(1,1))*dsqrt(v(4,4))
v(1,5) = rhomat(1,5)*dsqrt(v(1,1))*dsqrt(v(5,5))
v(1,6) = rhomat(1,6)*dsqrt(v(1,1))*dsqrt(v(6,6))
v(1,7) = rhomat(1,7)*dsqrt(v(1,1))*dsqrt(v(7,7))
v(2,1)=v(1,2)
v(3,1)=v(1,3)
v(4,1)=v(1,4)
v(5,1)=v(1,5)
v(6,1)=v(1,6)
v(7,1)=v(1,7)
v(2,3) = rhomat(2,3)*dsqrt(v(2,2))*dsqrt(v(3,3))
v(2,4) = rhomat(2,4)*dsqrt(v(2,2))*dsqrt(v(4,4))
v(2,5) = rhomat(2,5)*dsqrt(v(2,2))*dsqrt(v(5,5))
v(2,6) = rhomat(2,6)*dsqrt(v(2,2))*dsqrt(v(6,6))
v(2,7) = rhomat(2,7)*dsqrt(v(2,2))*dsqrt(v(7,7))
v(3,2)=v(2,3)
v(4,2)=v(2,4)
v(5,2)=v(2,5)
v(6,2)=v(2,6)
v(7,2)=v(2,7)
v(3,4) = rhomat(3,4)*dsqrt(v(3,3))*dsqrt(v(4,4))
v(3,5) = rhomat(3,5)*dsqrt(v(3,3))*dsqrt(v(5,5))
v(3,6) = rhomat(3,6)*dsqrt(v(3,3))*dsqrt(v(6,6))
v(3,7) = rhomat(3,7)*dsqrt(v(3,3))*dsqrt(v(7,7))
v(4,5) = rhomat(4,5)*dsqrt(v(4,4))*dsqrt(v(5,5))
v(4,6) = rhomat(4,6)*dsqrt(v(4,4))*dsqrt(v(6,6))
v(4,7) = rhomat(4,7)*dsqrt(v(4,4))*dsqrt(v(7,7))
v(5,4)=v(4,5)
v(5,3)=v(3,5)
v(4,3)=v(3,4)
v(5,3)=v(3,5)
v(6,3)=v(3,6)
v(7,3)=v(3,7)
v(5,4)=v(4,5)
v(6,4)=v(4,6)
v(7,4)=v(4,7)


i5: do i = 1,7
rhomat(i,i) = 1.0d-0
end do i5
rhomat(2,1)=rhomat(1,2)
rhomat(3,1)=rhomat(1,3)
rhomat(4,1)=rhomat(1,4)
rhomat(5,1)=rhomat(1,5)
rhomat(3,2)=rhomat(2,3)
rhomat(4,2)=rhomat(2,4)
rhomat(5,2)=rhomat(2,5)
rhomat(4,3)=rhomat(3,4)
rhomat(5,3)=rhomat(3,5)
rhomat(5,4)=rhomat(4,5)
rhomat(6,5)=rhomat(5,6)
rhomat(7,5)=rhomat(5,7)
rhomat(7,6)=rhomat(6,7)

end subroutine extract_var
