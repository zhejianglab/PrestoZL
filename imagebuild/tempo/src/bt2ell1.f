c      $Id$
      subroutine bt2ell1()

C  Converts orbital parameters for BT, DD, ... (a1,e,Pb,omz,T0,...)
C  into parameters for ELL1 model: eps1, eps2, T0 (new, defined with
C  respect to ascending node).

      implicit real*8(a-h,o-z)
      parameter(TWOPI=6.28318530717958648d0)
      include 'dim.h'
      include 'orbit.h'

      t0asc = t0(1)-omz(1)/360.d0*pb(1)

      om    = omz(1)*TWOPI/360.d0
      eps1  = e(1)*DSIN(om)
      eps2  = e(1)*DCOS(om)

      an    = TWOPI/pb(1) + omdot*TWOPI/360.d0/365.25d0
      pb(1) = TWOPI/an

      t0(1)=0.
      e(1)=0.
      omz(1)=0.

      return
      end

