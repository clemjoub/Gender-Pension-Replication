
subroutine readtheta2(silent,pathoutput,reform,thetamtype,thetamtype1,thetamtype2,thetamreftype,thetamreftype1,thetamreftype2,numper,numz,ktp)
integer:: silent,reform,group,xx,zz,numper,numz,ktp,ttt
character*200::pathoutput
real(8) thetamref(numper,3,numz),thetamreftype(ktp,numper,3,numz),thetamreftype1(numper,3,numz),thetamreftype2(numper,3,numz)
real(8) thetamtype(ktp,numper,3,numz),thetamtype1(numper,3,numz),thetamtype2(numper,3,numz)

open(unit=400,file=adjustl(trim(pathoutput))//'../inputs/pretheta', status='UNKNOWN',form='FORMATTED', &
        action='READ')
rewind(unit=400)
if (reform.eq.1) then
        open(unit=500,file=adjustl(trim(pathoutput))//'../inputs/posttheta', status='UNKNOWN',form='FORMATTED', &
            action='READ')
        rewind(unit=500)
end if

do ttt=1,ktp
    if (silent.ne.1) write(*,*) 'reading type', ttt
    
    xxloop: do xx = 1,numper
    grouploop: do group = 1,3
    zzloop: do zz = 1,numz
        read(400,*) thetamtype(ttt,xx,group,zz)
        if (reform.eq.1) then
            read(500,*) thetamreftype(ttt,xx,group,zz)
        else
                thetamreftype(ttt,xx,group,zz) = thetamtype(ttt,xx,group,zz)
        end if 
    end do zzloop 
    end do grouploop 
    end do xxloop

enddo

close(unit=400)
close(unit=500)
    33        format(f16.4)
    32        format(f24.8)

endsubroutine
