

*----------------------------------------------------------

      SUBROUTINE LUDCMP (A,N,NP,INDX,D)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *

*----------------------------------------------------------

      implicit none

      INTEGER N,NP,INDX
      INTEGER NMAX,IMAX
      INTEGER I,J,K

      DOUBLE PRECISION A,D,VV,AAMAX,SUM,DUM,TINY

      PARAMETER (NMAX = 100,TINY = 1.D-20)

      DIMENSION A(NP,NP),INDX(N),VV(NMAX)

      D = 1.D0
      DO I = 1,N
         AAMAX = 0.D0
         DO J = 1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX = ABS(A(I,J))
c            if (i.eq.15) print *, 'j',j,'A', A(i,j)
         ENDDO
         if (AAMAX.eq.0.d0) then
            print *, 'i',i
            print *,'abs A = AAMAX', AAMAX
         endif
         IF (AAMAX.EQ.0.D0) STOP 'Singular matrix.'
         VV(I) = 1.D0/AAMAX
      ENDDO
      DO J = 1,N
         IF (J.GT.1) THEN
            DO I = 1,J-1
               SUM = A(I,J)
               IF (I.GT.1) THEN
                  DO K = 1,I-1
                     SUM = SUM-A(I,K)*A(K,J)
                  ENDDO
                  A(I,J) = SUM
               ENDIF
            ENDDO
         ENDIF
         AAMAX = 0.D0
         DO I = J,N
            SUM = A(I,J)
            IF (J.GT.1) THEN
               DO K = 1,J-1
                  SUM = SUM-A(I,K)*A(K,J)
               ENDDO
               A(I,J) = SUM
            ENDIF
            DUM = VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX = I
               AAMAX = DUM
            ENDIF
         ENDDO
         IF (J.NE.IMAX) THEN
            DO K = 1,N
               DUM = A(IMAX,K)
               A(IMAX,K) = A(J,K)
               A(J,K) = DUM
            ENDDO
            D = -D
            VV(IMAX) = VV(J)
         ENDIF
         INDX(J) = IMAX
         IF (J.NE.N) THEN
            IF (A(J,J).EQ.0.D0) A(J,J) = TINY
            DUM = 1.D0/A(J,J)
            DO I = J+1,N
               A(I,J) = A(I,J)*DUM
            ENDDO
         ENDIF
      ENDDO
      IF (A(N,N).EQ.0.D0) A(N,N) = TINY

      RETURN
      END
