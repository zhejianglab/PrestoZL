c      $Id$
c
      subroutine ell12bt()

C  Converts ELL1 parameters to BT parameters
C  N Wex, RNM May, 2000

      implicit real*8(a-h,o-z)
      parameter(twopi=6.28318530717958648d0)
      include 'dim.h'
      include 'orbit.h'

      e(1)   = dsqrt(eps1**2 + eps2**2)
      om     = datan2(eps1,eps2)
      omz(1) = om*360.d0/twopi
      t0(1)  = t0asc + pb(1)*om/twopi

      return
      end
