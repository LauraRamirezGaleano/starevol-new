

*----------------------------------------------------------

      SUBROUTINE LUBKSB (A,N,NP,INDX,B)
*                                                                      *
* $LastChangedDate:: 2014-02-04 11:25:49 +0100 (Tue, 04 Feb 2014)    $ *
* $Author:: decress1                                                 $ *
* $Rev:: 5                                                           $ *
*                                                                      *

*----------------------------------------------------------

      implicit none

      INTEGER N,NP,INDX
      INTEGER I,II,J,LL

      DOUBLE PRECISION A,B,SUM

      DIMENSION A(NP,NP),INDX(N),B(N)

      II = 0
      DO I = 1,N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF (II.NE.0) THEN
            DO J = II,I-1
               SUM = SUM-A(I,J)*B(J)
            ENDDO
         ELSE IF (SUM.NE.0.D0) THEN
            II = I
         ENDIF
         B(I) = SUM
      ENDDO
      DO I = N,1,-1
         SUM = B(I)
         IF (I.LT.N) THEN
            DO J = I+1,N
               SUM = SUM-A(I,J)*B(J)
            ENDDO
         ENDIF
         B(I) = SUM/A(I,I)
      ENDDO

      RETURN
      END
