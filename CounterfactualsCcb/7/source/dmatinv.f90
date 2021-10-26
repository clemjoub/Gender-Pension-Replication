

      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      implicit none

      REAL(8) A(NP,NP)
      REAL(8) VV(400),D,SUM,AAMAX,DUM
      REAL(8) TINY
      INTEGER INDX(N),N,NP,IMAX,NMAX
      INTEGER I,J,K


      NMAX = 400
      TINY = 1.0d-20 * 1.0d-0
      D=1.0d-0
      VV = 0.0d-0
      SUM = 0.0d-0

      DO 12 I=1,N
        AAMAX=0.0d-0
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.0d-0) then
        continue
!        PAUSE 'Singular matrix.'
        endif
        VV(I)=1.0d-0/AAMAX
12    CONTINUE


      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-1.0d-0 *A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.0d-0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-1.d-0 * D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.0d-0)A(J,J)=TINY
          DUM=1.0d-0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.0d-0)A(N,N)=TINY


      RETURN
      END

!---------------------------------------------------------

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      implicit none

      REAL(8) A(NP,NP)
      REAL(8) B(N),SUM
      INTEGER INDX(N),N,NP,I,II,LL,J,K

      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.0d-0) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
     
 
      RETURN
      END

!----------------------------------
       subroutine dmatinv(a,y,np)

       implicit none
 
       integer np
       real(8)  a(np,np),y(np,np)
       real(8), allocatable :: ac(:,:)
       real(8) d,b(np)
       integer indx(np),i,j,n

 !      external LUDCMP,LUBKSB

       ALLOCATE(ac(np,np))

!      make a copy of a so it will not be affected by 
!      the subroutines and initialize y

       y = 0.0d-0
       loop1: do i = 1,np
         ac(i,i) = a(i,i)
         y(i,i) = 1.0d-0
       loop2: do j = 1,np
         ac(i,j) = a(i,j)
       end do loop2
       end do loop1
       indx(1:np)=0
       n = np

       call LUDCMP(ac,np,np,indx,d)

       loop13: do j = 1,n
          call LUBKSB(ac,n,np,indx,y(1,j))
       end do loop13

       DEALLOCATE(ac)

       return
       end subroutine dmatinv

