c      $Id$
      FUNCTION DOT(A,B)
C **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C           COMPUTES DOT PRODUCT OF TWO VECTORS
      DIMENSION A(3),B(3)
      DOT=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      RETURN
      END
