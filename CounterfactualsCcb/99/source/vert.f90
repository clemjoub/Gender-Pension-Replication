
!
!      ________________________________________________________
!     |                                                        |
!     |                INVERT A GENERAL MATRIX                 |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         V     --ARRAY CONTAINING MATRIX                |
!     |                                                        |
!     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
!     |                                                        |
!     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
!     |                                                        |
!     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
!     |                      ELEMENTS                          |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         V     --INVERSE                                |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: DABS                             |
!     |________________________________________________________|
!
      SUBROUTINE VERT(V,LV,N,W)

      implicit none

      REAL(8) V(LV,N),S,T
      INTEGER W(N),I,J,K,L,M,N,P,LV


      IF ( N .EQ. 1 ) GOTO 110
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
!     ---------------------------------------
!     |*** FIND PIVOT AND START ROW SWAP ***|
!     ---------------------------------------
      P = L
      IF ( M .GT. N ) GOTO 30
      S = DABS(V(L,L))
      DO 20 I = M,N
           T = DABS(V(I,L))
           IF ( T .LE. S ) GOTO 20
           P = I
           S = T
20    CONTINUE
      W(L) = P
30    S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0.0d-0 ) GOTO 120
!     -----------------------------
!     |*** COMPUTE MULTIPLIERS ***|
!     -----------------------------
      V(L,L) = -1.0d-0
      S = 1.0d-0/S
      DO 40 I = 1,N
40         V(I,L) = -S*V(I,L)
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0.0d-0 ) GOTO 50
!     ------------------------------
!     |*** ELIMINATE BY COLUMNS ***|
!     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO 60 I = 1,K
60         V(I,J) = V(I,J) + T*V(I,L)
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO 80 I = M,N
80         V(I,J) = V(I,J) + T*V(I,L)
      GOTO 50
!     -----------------------
!     |*** PIVOT COLUMNS ***|
!     -----------------------
90    L = W(K)
      DO 100 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,K)
100        V(I,K) = T
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
      RETURN
110   IF ( V(1,1) .EQ. 0.0d-0 ) GOTO 120
      V(1,1) = 1.0d-0/V(1,1)
      RETURN
120   WRITE(6,*) 'ERROR: MATRIX HAS NO INVERSE'
      STOP

      write(*,*) 'exiting vert'
      END