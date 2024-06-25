c      $Id$
      SUBROUTINE PCSHFT(A,B,D,N)
	implicit real*8 (a-h,o-z)
      DIMENSION D(N)
      CONST=2./(B-A)
      FAC=CONST
      DO 11 J=2,N
        D(J)=D(J)*FAC
        FAC=FAC*CONST
11    CONTINUE
      CONST=0.5*(A+B)
      DO 13 J=1,N-1
        DO 12 K=N-1,J,-1
          D(K)=D(K)-CONST*D(K+1)
12      CONTINUE
13    CONTINUE
      RETURN
      END
