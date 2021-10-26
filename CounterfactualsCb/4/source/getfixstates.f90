

    ! -----------------------------------------------------------------------------------------------------
    !       This subroutines constructs the state variables (at beginning of "ctt" period) from the information contained in
    !       fixvar,ivart,rvart - it returns the Z vector 
    !       
    SUBROUTINE getfixstate(o,ctt,fixvar,ivart,rvart,experh,experw,numk,fexperh,fexperw,hcaph,hcapw,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage)


    implicit none

    integer, parameter:: nseg=7,nseg2=6,nseg3=5
    integer totp,numper,o,ct,wage,hage,k,j,kk,singlemale,singlefemale,i1,i2, lastage,lastworkage
    real(8) rvart(6,1,numper)      ! row 1 = wealth - safe
                                   ! row 2 = wealth - risky
                                   ! row 3 = rate of return safe
                                   ! row 4 = rate of return risky
                                   ! row 5 = pension wealth women
                                   ! row 6 = pension wealth men
    integer ivart(9,1,numper)      ! row 1 = children (n)
                                   ! row 2 = married
                                   ! row 3 = work h
                                   ! row 4 = work w
                                   ! row 5 = survivh
                                   ! row 6 = survivw
                                   ! row 7 = marriage duration
                                   ! row 8 = formal h 
                                   ! row 9 = formal w
    integer fixvar(4,1)            ! row 1 = hed
                                   ! row 2 = wed
                                   ! row 3 = hage
                                   ! row 4 = female 


    integer numz,n(numper),i,full,part,reform,postref
    integer hed,wed,intmat(75,75),ind,married,ctt
    real(8) Ztemp(1,numz),experw,experh,wealth,indw1,indw2,c1,c2,hcaph,hcapw,cmin,fexperh,fexperw
    integer tt,female,changestat,ind65h,ind70h,evermarried
    integer dech,decw,pos1,temp1,temp2,nowealth,mardur,numk,ind80h,ind90h
    real(8) ret1, ret2,pensionw,pensionh,pmin,pmax,cnu(76,4)
    real(8) nsafe, nrisky,hw,ww,experw2,experh2,fpenwithdraw,mpenwithdraw,withdrawh(nseg2),withdraww(nseg2),withdsegh(nseg2),withdsegw(nseg2),withdh,withdw
    real(8) maxearnf,maxearnh,pasis,mpg,temp111f,spsf,temp111m,spsm,basecons,seg(nseg),apwseg(nseg),logapwreg(nseg),apwreg(nseg),crraapwreg(nseg),apw,logapw,crraapw
    integer formalh,formalw,informalh,informalw,male
    integer pasisqual,consperh,workperh,consperw,workperw
    real pseg(nseg2),mpgqualw,mpgqualh,xpseg(nseg3),hxpseg(nseg3),wxpseg(nseg3)


    informalh=0
    informalw=0
    ct = ctt-1
    wealth = rvart(1,1,ct) !use counter ct because it is end of period wealth
    pensionw = rvart(5,1,ct)
    pensionh = rvart(6,1,ct)
    married = ivart(2,1,ctt)
    evermarried = ivart(2,1,20)   !if married at age 35, then evermarried equals 1
    hed = fixvar(1,1)
    wed = fixvar(2,1)
    female = fixvar(4,1)
    hage = max0(35,min0((fixvar(3,1)-1+ctt),90))  
    wage = min0(ctt+15,90)
    hw = ivart(3,1,ct)                 !employment status
    ww = ivart(4,1,ct)                 !employment status
    formalw = ivart(9,1,ct)   ! indicator for female worked in formal sector
    formalh = ivart(8,1,ct)   ! indicator for male worked in formal sector
    if (hw.eq.1.and.formalh.eq.0) then
         informalh = 1
    end if 
    if (ww.ge.1.and.formalw.eq.0) then
         informalw = 1
    end if 
    mpg = 2.0d-0*pasis
    male = 1-female
    singlemale=male*(1-married)
    singlefemale=female*(1-married)
    mpgqualh = dmax1(20.0d-0-fexperh,0.0d-0)
    mpgqualw = dmax1(20.0d-0-fexperw,0.0d-0)

! income effect regressors (based on "annualized permanent wealth")

    !consumption periods left
    consperh=max0(1,lastage-hage) 
    consperw=max0(1,lastage-wage) 

    !work periods left
    workperh=max0(0,lastworkage-hage)
    workperw=max0(0,lastworkage-wage)

    !wealth segments
    basecons=2000000.0d-0
    seg(1)=0.0
    seg(2)=0.5*basecons
    seg(3)=basecons
    seg(4)=1.5*basecons
    seg(5)=3*basecons
    seg(6)=5*basecons
    seg(7)=-9999999999999999999.0d-0

    apw=(wealth+pensionh+pensionw+workperh*hcaph*(1-singlefemale)+workperw*hcapw*(1-singlemale))/(consperh+consperw) !current wealth plus earnings from remaining work periods divided by remaining years

    do i=1,nseg-1
        apwseg(i)=(apw.ge.seg(i))*(-1.0d-0)
        apwreg(i)=apwseg(i)*(apw-seg(i))
    enddo


!implicit marginal tax regressors (based on pension wealth accumulations interacted with apw)

    pseg(1)=0.0
    pseg(2)=basecons/4
    pseg(3)=basecons/2
    pseg(4)=basecons
    pseg(5)=basecons*2
    pseg(6)=999999999999.0

    call piksrt(6,pseg)

    withdh=0.0d-0
    withdw=0.0d-0

    if (married.eq.1) then
    withdh=pensionh/dmax1(1.0d-0,(CNU(hage-34,4)-workperh*1.0d-0))
    withdw=pensionw/dmax1(1.0d-0,(CNU(wage-34,3)-workperw*1.0d-0))
    elseif (married.eq.0) then
    withdh=pensionh/dmax1(1.0d-0,(CNU(hage-34,2)-workperh*1.0d-0))
    withdw=pensionw/dmax1(1.0d-0,(CNU(wage-34,3)-workperw*1.0d-0))
    endif
    do i=1,nseg2-1 
    withdsegh(i)=(withdh>=pseg(i))*(-1.0d-0)
    withdrawh(i)=withdsegh(i)*(withdh-pseg(i))*apw
    withdsegw(i)=(withdw>=pseg(i))*(-1.0d-0)
    withdraww(i)=withdsegw(i)*(withdw-pseg(i))*apw
    enddo

!impact of years of contribution restrictions (based on number of years in formal sector)

    xpseg(1)=0.0
    xpseg(2)=10.0
    xpseg(3)=20.0
    xpseg(4)=30.0
    xpseg(5)=9999999999999.0

    do i=1,nseg3-1
    hxpseg(i)=(fexperh.ge.xpseg(i+1))*(-1.0d-0)
    wxpseg(i)=(fexperw.ge.xpseg(i+1))*(-1.0d-0)
    enddo
    
! fill the Z matrix (some vectors may be zero)

!linear terms
    Ztemp(1,1) = 1.0d-0
    Ztemp(1,2) = female*1.0d-0  !collinear
    Ztemp(1,3) = married*1.0d-0  !collinear
    Ztemp(1,4) = (numk==0)*1.0d-0  !collinear
    Ztemp(1,5) = max0(numk,1)*1.0d-0  !collinear
    Ztemp(1,6) = hed*1.0d-0 *(1-singlefemale)
    Ztemp(1,7) = wed*1.0d-0 *(1-singlemale)
    Ztemp(1,10) = hcaph*1.0d-0 *(1-singlefemale)*(hage.le.lastworkage)
    Ztemp(1,11) = hcapw*1.0d-0 *(1-singlemale)*(wage.le.lastworkage)
    !Ztemp(1,12) = formalh*1.0d-0 *(1-singlefemale)*(hage.lt.lastworkage)
    !Ztemp(1,13) = formalw*1.0d-0 *(1-singlemale)*(wage.lt.lastworkage)
    !Ztemp(1,14) = informalh*1.0d-0 *(1-singlefemale)*(hage.lt.lastworkage)
    !Ztemp(1,15) = informalw*1.0d-0*(1-singlemale)*(wage.lt.lastworkage)
    !Ztemp(1,16) = wealth*1.0d-0 !collinear
    !Ztemp(1,17) = pensionh*1.0d-0 *(1-singlefemale)
    !Ztemp(1,18) = pensionw*1.0d-0*(1-singlemale)
    !Ztemp(1,19) = fexperh*1.0d-0*(1-singlefemale)
    !Ztemp(1,20) = fexperw*1.0d-0*(1-singlemale)
    Ztemp(1,21) = (numk==1)*1.0d-0
    Ztemp(1,22) = (numk==2)*1.0d-0
    Ztemp(1,23) = (numk==3)*1.0d-0
    Ztemp(1,24) = (numk==4)*1.0d-0 !collinear
    Ztemp(1,25) = 0.0d-0
    Ztemp(1,26) = 0.0d-0
    Ztemp(1,27) = 0.0d-0
    Ztemp(1,28) = 0.0d-0
    Ztemp(1,29) = 0.0d-0
            
! income effect regressors 

    !Ztemp(1,30) = apwreg(1)*1.0d-0 
    !Ztemp(1,31) = apwreg(2)*1.0d-0 
    !Ztemp(1,32) = apwreg(3)*1.0d-0 
    !Ztemp(1,33) = apwreg(4)*1.0d-0 
    !Ztemp(1,34) = apwreg(5)*1.0d-0 
    !Ztemp(1,35) = apwseg(1)*1.0d-0  
    !Ztemp(1,36) = apwseg(2)*1.0d-0 
    !Ztemp(1,37) = apwseg(3)*1.0d-0 
    !Ztemp(1,38) = apwseg(4)*1.0d-0 
            
! Age differences

    Ztemp(1,39) = (hage-wage) * married * (ctt.lt.72)* (-1.0d-0)
    Ztemp(1,40) = (hage==wage+1) * married * (-1.0d-0)        
    Ztemp(1,41) = (hage==wage+2) * married * (-1.0d-0)      
    Ztemp(1,42) = (hage==wage+3) * married * (-1.0d-0)        
    Ztemp(1,43) = (hage==wage+4) * married * (-1.0d-0)        
    Ztemp(1,44) = (hage==wage+5) * married * (-1.0d-0)         
    Ztemp(1,45) = (hage==wage+6) * married * (-1.0d-0)         
    Ztemp(1,46) = (hage==wage+7) * married * (-1.0d-0)        
    Ztemp(1,47) = (hage==wage+8) * married * (-1.0d-0)       
    Ztemp(1,48) = (hage==wage+9) * married * (-1.0d-0)       
    Ztemp(1,49) = (hage==wage+10) * married * (-1.0d-0)        

! implicit marginal tax regressors

    !Ztemp(1,50) = withdrawh(1)*1.0d-0 *(1-singlefemale)
    !Ztemp(1,51) = withdrawh(2)*1.0d-0 *(1-singlefemale)
    !Ztemp(1,52) = withdrawh(3)*1.0d-0 *(1-singlefemale)
    !Ztemp(1,53) = withdrawh(4)*1.0d-0 *(1-singlefemale)
    !Ztemp(1,54) = withdrawh(5)*1.0d-0 *(1-singlefemale)
    !Ztemp(1,55) = withdraww(1)*1.0d-0 *(1-singlemale)
    !Ztemp(1,56) = withdraww(2)*1.0d-0 *(1-singlemale)
    !Ztemp(1,57) = withdraww(3)*1.0d-0 *(1-singlemale)
    !Ztemp(1,58) = withdraww(4)*1.0d-0 *(1-singlemale)
    !Ztemp(1,59) = withdraww(5)*1.0d-0 *(1-singlemale)    
    Ztemp(1,60) = 0.0d-0
    Ztemp(1,61) = 0.0d-0
    Ztemp(1,62) = 0.0d-0
    Ztemp(1,63) = 0.0d-0
    Ztemp(1,64) = 0.0d-0
    
! Formal experience regressors

    !Ztemp(1,65) = hxpseg(1)*1.0d-0*(1-singlefemale)*(ctt.gt.20)
    !Ztemp(1,66) = hxpseg(2)*1.0d-0*(1-singlefemale)*(ctt.gt.30)
    !Ztemp(1,67) = hxpseg(3)*1.0d-0*(1-singlefemale)*(ctt.gt.40)
    Ztemp(1,68) = 0.0d-0
    Ztemp(1,69) = 0.0d-0
    !Ztemp(1,70) = wxpseg(1)*1.0d-0*(1-singlemale)*(ctt.gt.20)
    !Ztemp(1,71) = wxpseg(2)*1.0d-0*(1-singlemale)*(ctt.gt.30)
    !Ztemp(1,72) = wxpseg(3)*1.0d-0*(1-singlemale)*(ctt.gt.40)
    Ztemp(1,73) = 0.0d-0
    Ztemp(1,74) = 0.0d-0
    Ztemp(1,75) = 0.0d-0


    ind = 75

! make the interacted variables 
    !loop2: do i = 1,50
     !      if (intmat(i,1)==-9) exit loop2
      !     i1=intmat(i,1)
       !    i2=intmat(i,2)
        !   Ztemp(1,ind+i)=Ztemp(1,i1)*Ztemp(1,i2)
    !end do loop2

    return
    end



