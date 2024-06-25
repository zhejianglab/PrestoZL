c      $Id$
      SUBROUTINE CHEBPC(C,D,N)
	implicit real*8 (a-h,o-z)
      PARAMETER (NMAX=50)
      DIMENSION C(N),D(N),DD(NMAX)
      DO 11 J=1,N
        D(J)=0.
        DD(J)=0.
11    CONTINUE
      D(1)=C(N)
      DO 13 J=N-1,2,-1
        DO 12 K=N-J+1,2,-1
          SV=D(K)
          D(K)=2.*D(K-1)-DD(K)
          DD(K)=SV
12      CONTINUE
        SV=D(1)
        D(1)=-DD(1)+C(J)
        DD(1)=SV
13    CONTINUE
      DO 14 J=N,2,-1
        D(J)=D(J-1)-DD(J)
14    CONTINUE
      D(1)=-DD(1)+0.5*C(1)
      RETURN
      END
