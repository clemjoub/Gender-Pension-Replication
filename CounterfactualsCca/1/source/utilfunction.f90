

subroutine utilfunction2(alpha,ncoeff,ktp,ct,mptype,cbis,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,nsafebis,npensionhbis,npensionwbis,beta,phdec,pwdec,ind65h,ind60w,yyevec,hage,Tbar,util,utilf,utilm,ubequest)

integer:: ncoeff,ktp,ind65h,ind60w,hage,Tbar,ct,mptype,formh,formw,hw,ww,fswitch,freentry,mswitch,mreentry,married,female,numk,parttime,fulltime
real(8):: nsafe,ubequest,ubequest1,ubequest2,alpha(ncoeff,ktp),npensionh,npensionw,beta,phdec,pwdec,yyevec(7),util,c,utilf,utilm,cf,cm
real(8):: MUCshiftf,MUCshiftm,MULf,MULm,switchf,switchm,nonpecf,nonpecm,cbis,nsafebis,npensionhbis,npensionwbis

c=0.0d-0
nsafe=0.0d-0
npensionh=0.0d-0
npensionw=0.0d-0
c=cbis/1000000.0d-0 !express consumption in million pesos
nsafe=nsafebis/1000000.0d-0
npensionh=npensionhbis/1000000.0d-0
npensionw=npensionwbis/1000000.0d-0

parttime=0
fulltime=0
if (ww.eq.1) then 
fulltime = 1
else if (ww.eq.2) then
parttime = 1
end if

util = 0.0d-0
utilm = 0.0d-0
utilf = 0.0d-0

cf = c/(1+married)
cm = c/(1+married)       

MUCshiftf=dexp(alpha(2,1)*dble(numk)+ alpha(4,1)*dble((1-parttime)*(1-fulltime)))
MULf=alpha(6,mptype)*dexp(alpha(10,1)*dble(numk) +alpha(11,1)*married + alpha(14,1)*dble(ind60w) + yyevec(7))
switchf=dble(fswitch)*alpha(63,1)+dble(freentry)*alpha(65,1)
nonpecf=dble(1-formw)*dble(parttime+fulltime)*alpha(70,mptype)

MUCshiftm=dexp(alpha(3,1)*dble(numk)+alpha(5,1)*dble(1-hw))
MULm=alpha(7,mptype)*dexp(alpha(13,1)*dble(ind65h) + yyevec(6))
switchm=dble(mswitch)*alpha(62,1)+dble(mreentry)*alpha(64,1)
nonpecm=dble(1-formh)*dble(hw)*alpha(71,mptype)


utilf = ((cf**alpha(1,1))/alpha(1,1)) * MUCshiftf &
+ (dble((1-parttime)*(1-fulltime))+alpha(8,1)*dble(parttime))* MULf &
- switchf + nonpecf   

utilm = ((cm**alpha(1,1))/alpha(1,1)) * MUCshiftm &
+ dble(1-hw)* MULm &
- switchm + nonpecm 

util = utilf*dble(married+(1-married)*female)+utilm*dble(married+(1-married)*(1-female))

!bequest

ubequest1=alpha(81,1)*dble(maxval((/numk,1/))) !number of periods and kids over which the bequest will be split
ubequest2=(nsafe+npensionh+npensionw)/(ubequest1+0.00001)+alpha(83,1) !split bequest + baseline income of heirs
ubequest=ubequest1*((ubequest2+0.00001)**alpha(1,1))/alpha(1,1) !utility of heirs'consumption

!add bequest if last period
if ((ct+15).eq.Tbar.or.(hage.eq.90.and.married.eq.0.and.female.eq.0)) util= util + beta*ubequest


endsubroutine
