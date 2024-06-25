c      $Id: cross.f,v 1.1 2002/11/25 15:03:41 wvanstra Exp $
      SUBROUTINE CROSS(A,B,R)
C **********************************************************************
      IMPLICIT NONE
C
C           COMPUTES CROSS PRODUCT OF TWO VECTORS (R=A X B)
      real*8 A(3),B(3),R(3)
      R(1)=A(2)*B(3)-A(3)*B(2)
      R(2)=A(3)*B(1)-A(1)*B(3)
      R(3)=A(1)*B(2)-A(2)*B(1)
      END
