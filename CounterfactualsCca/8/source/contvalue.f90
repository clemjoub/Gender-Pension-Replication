
subroutine contvalue(married,female,numper,term4,ubequest,thetap,reform,nodivrule,nivart,nrvart,nfixvar,ct,o,phdec,pwdec,pdiv,probfert,nsafe,npensionh,npensionw,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,numz,totp,intmat,pmax,pmin,pasis,bono,cmin,cnu,lastage,lastworkage)

!!!INCORPORATE LASTAGE


integer:: married,female,group,ctt,numz,reform,nodivrule,ct,o,numper,numk,totp,intmat(75,75),nivart(9,1,numper),nfixvar(4,1),lastage,lastworkage
real(8):: ubequest,wsum,wsum1,wsum1d,wsum2,wsum2d,wsumf,wsum1f,wsum1fd,wsum2f,term4,term42,term41,phdec,pwdec,pdiv,bono,probfert
real(8):: Ztemp(1,numz), thetap(3,numz),nrvart(6,1,numper)
real(8):: experh,experw,fexperh,fexperw,maxearnh,maxearnf,pmax,pmin,pasis,cmin,cnu(76,4),npensionh,npensionw,nsafe

!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!        Compute EMAX terms  if NOT last period

 
    wsum = 0.0d-0
    wsum1 = 0.0d-0
    wsum1d = 0.0d-0
    wsum2 = 0.0d-0
    wsum2d = 0.0d-0
    wsumf = 0.0d-0
    wsum1f = 0.0d-0
    wsum1fd = 0.0d-0
    wsum2f = 0.0d-0
    
       ! CASE 1: couple's continuation value: no fertility => wsum
   group=1 !married couple
   ctt=ct+1
   nfixvar(4,1) = female                        ! household gender stays same
   nivart(1,1,ctt) = 0                          ! fert is 0
   nivart(2,1,ctt) = 1                          ! stay married, without fertility
   nivart(7,1,ctt) = nivart(7,1,ct)+1           ! mardur will be one year longer if stay married 
   nrvart(1,1,ct) = nsafe          
   nrvart(5,1,ct) = npensionw          
   nrvart(6,1,ct) = npensionh
   call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
   wsum = dot_product(Ztemp(1,:),thetap(group,:))
!   if(ct.eq.74.and.o.eq.1) then
!   print*,wsum
!   print*,Ztemp(1,:)
!   print*,thetap(group,:)
!   continue
!   endif


       ! CASE 2: single or divorced female continuation value: no fertility => wsum1
   group=2 !single female
   nfixvar(4,1) = 1                             ! household is female
   nivart(1,1,ctt) = 0                          ! fert is 0
   nivart(2,1,ctt) = 0                          ! become an unmarried female (divorce), without fertility
   nivart(7,1,ctt) = nivart(7,1,ct)             ! mardur will stay the same if get divorced/widowed        
   nrvart(1,1,ct)=0.5d-0*nsafe                  ! get half of private savings
   if (reform.eq.1.and.nodivrule.eq.0) then
     nrvart(5,1,ct)= (npensionh+npensionw)/2.0d-0    
     nrvart(6,1,ct)= (npensionh+npensionw)/2.0d-0    
   else
     nrvart(5,1,ct) = npensionw
     nrvart(6,1,ct) = npensionh 
   end if
   call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
   wsum1 = dot_product(Ztemp(1,:),thetap(group,:))
   

       ! CASE 3: wife's continuation value after husband dies => wsum1d
   group=2 !single female
   nfixvar(4,1) = 1                             ! female
   nivart(1,1,ctt) = 0                          ! fert is 0
   nivart(2,1,ctt) = 0                          ! become  unmarried (death of husband)
   nivart(7,1,ctt) = nivart(7,1,ct)             ! mardur will stay the same if get divorced/widowed   
   nrvart(1,1,ct) = nsafe                       ! keeps all the private savings          
   nrvart(5,1,ct)= npensionw+0.6d-0*npensionh   ! survivorship pension put in wifes account 
   nrvart(6,1,ct)= 0.0d-0     
   call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
   wsum1d = dot_product(Ztemp(1,:),thetap(group,:))


       ! CASE 4: single or divorced male's continuation value after divorce => wsum2
    group=3                                     ! single  male
    nfixvar(4,1) = 0                            ! male
    nivart(1,1,ctt) = 0                         ! fert is 0 
    nivart(2,1,ctt) = 0                         ! become an unmarried male (divorce), without fertility 
    nivart(7,1,ctt) = nivart(7,1,ct)            ! mardur will stay the same if get divorced/widowed        
    nrvart(1,1,ct) = 0.5d-0*nsafe               ! get half of private savings
    if (reform.eq.1.and.nodivrule.eq.0) then
     nrvart(5,1,ct)= (npensionh+npensionw)/2.0d-0    
     nrvart(6,1,ct)= (npensionh+npensionw)/2.0d-0    
   else
     nrvart(5,1,ct) = npensionw
     nrvart(6,1,ct) = npensionh
    end if
    call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
    wsum2 = dot_product(Ztemp(1,:),thetap(group,:))


       ! CASE 5: husband's continuation value after wife dies => wsum2d
    group=3         
    nivart(1,1,ctt) = 0                         ! fert is 0 
    nivart(2,1,ctt) = 0                         ! become an unmarried male (widowed), without fertility 
    nivart(7,1,ctt) = nivart(7,1,ct) 
    nfixvar(4,1) = 0
    nrvart(1,1,ct) = nsafe                      ! get all of private savings
    nrvart(5,1,ct)=0.0d-0                   
    nrvart(6,1,ct)= nrvart(6,1,ct)              !loses wifes pension 
    call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
    wsum2d = dot_product(Ztemp(1,:),thetap(group,:))

   ! fertility is an option if household includes  female< age 40
   if ((married.eq.1.or.female.eq.1).and.ct.le.25) then 
   
           ! CASE 6: couple's continuation value, with fertility => wsumf
       group=1 !married couple
       nfixvar(4,1)=female
       nivart(1,1,ctt) = 1                      ! fert - 1 
       nivart(2,1,ctt) = 1                      ! stay married, with fertility
       nivart(7,1,ctt) = nivart(7,1,ct)+1 
       nrvart(1,1,ct) = nsafe                      ! get all of private savings
       nrvart(5,1,ct)= npensionw   + bono  
       nrvart(6,1,ct)= npensionh
       call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk+1,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
       wsumf = dot_product(Ztemp(1,:),thetap(group,:))


            ! CASE 7: single or divorced woman continuation value, with fertility
        group=2                             !single female
        nfixvar(4,1) = 1
        nivart(1,1,ctt) = 1 
        nivart(2,1,ctt) = 0                 ! become an unmarried female(divorce), with fertility
        nivart(7,1,ctt) = nivart(7,1,ct) 
        nrvart(1,1,ct)=0.5d-0*nsafe
        if (reform.eq.1.and.nodivrule.eq.0) then
        nrvart(5,1,ct)= dmax1((npensionh+npensionw)/2.0d-0,npensionw)+bono 
        nrvart(6,1,ct) = 0.0d-0
        else
        nrvart(5,1,ct)= npensionw + bono
        nrvart(6,1,ct) = 0.0d-0
        end if
   	call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk+1,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
        wsum1f = dot_product(Ztemp(1,:),thetap(group,:))


            ! CASE 8: widow's continuation value, with fertility
        group=2 !single female
        nfixvar(4,1) = 1
        nivart(1,1,ctt) = 1 
        nivart(2,1,ctt) = 0          ! become an unmarried female(widowed), with fertility
        nivart(7,1,ctt) = nivart(7,1,ct) 
        nrvart(1,1,ct)= nsafe
        nrvart(5,1,ct)= npensionw+0.6d-0*npensionh  +bono  !get husbands pension put in wifes account 
        nrvart(6,1,ct)= 0.0d-0     
   	call getstate2(o,ctt,nfixvar,nivart,nrvart,experh,experw,numk+1,fexperh,fexperw,maxearnh,maxearnf,Ztemp,numz,totp,numper,intmat,pmax,pmin,pasis,cmin,cnu,lastage,lastworkage) 
        wsum1fd = dot_product(Ztemp(1,:),thetap(group,:))

    end if  !married.eq.1.or... 


! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------
!    add future values, taking account that there is some probability of
!    becoming divorced or deceased next period 


   ! if married, cont value does not depend on whether female sampled or not
   if (married.eq.1.and.ct.gt.25) then
   term4 = (1.0d-0 - phdec)*(1.0d-0 - pwdec)*(1.0d-0-pdiv)*wsum + &   !neither dies and no divorce
            phdec*wsum1d +   &                            !husband dies, with or without divorce
            (1.0d-0-phdec)*(1.0d-0 - pwdec)*pdiv*(0.5d-0*wsum1+0.5d-0*wsum2) + &     !neither dies, but get divorced
            pwdec*wsum2d          + &                         !wife dies, with or without divorce 
            pwdec*phdec*ubequest                            !both die, get utility from bequest
   end if         
   
   if (married.eq.1.and.ct.le.25) then
   term41 = (1.0d-0 - phdec)*(1.0d-0 - pwdec)*(1.0d-0-pdiv)*wsum + &   !neither dies and no divorce
            phdec*wsum1d +   &                             !husband dies, with or without divorce
            (1.0d-0-phdec)*(1.0d-0-pwdec)*pdiv*(0.5*wsum1+0.5*wsum2) + &      !neither dies, but get divorced 
            pwdec*wsum2d  +&
            pwdec*phdec*ubequest                            !both die, get utility from bequest
   term42 = (1.0d-0 - phdec)*(1.0d-0 - pwdec)*(1.0d-0-pdiv)*wsumf + &   !neither dies and no divorce
            phdec*wsum1fd +   &                              !husband dies, with or without divorce
            (1.0d-0-phdec)*(1.0d-0-pwdec)*pdiv*(0.5d-0*wsum1f+0.5d-0*wsum2) + &      !neither dies, but get divorced 
            pwdec*wsum2d    +&                               !assume no child if wife dies
               pwdec*phdec*ubequest                            !both die, get utility from bequest
   term4 = (1.0d-0-probfert)*term41+probfert*term42
   end if         

   if (married.eq.0.and.female.eq.1.and.ct.gt.25) then      !single or divorced female age gt 40
   term4 = (1.0d-0-pwdec)*wsum1  + pwdec * ubequest
   end if         
   if (married.eq.0.and.female.eq.1.and.ct.le.25) then      !single or divorced female age le 40
   term4 = (1.0d-0-probfert)*(1.0d-0-pwdec)*wsum1+probfert*(1.0d-0-pwdec)*wsum1f + pwdec* ubequest
   end if   
         
   if (married.eq.0.and.female.eq.0) then               !single or divorced male
   term4 = (1.0d-0-phdec)*wsum2 + phdec *ubequest
   end if         

endsubroutine contvalue
