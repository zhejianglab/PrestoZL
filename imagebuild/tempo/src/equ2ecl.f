c     $Id$
      subroutine equ2ecl(r)

c     rotate vector r(1..3) from equatorial to ecliptic coordinates
c     ==> rotate about "x" axis by angle epsilon

      implicit real*8 (a-h,o-z)

      include 'dim.h'
      include 'acom.h'

      real*8 r(3)
      real*8 tmpy, tmpz

c     historic content of this routine (until generalized for any value of 
c       epsilon, DJN, 2014-Jun-22:)
c     epsilon = 84381.412 arc sec at epoch of J2000.0 (IERS tech note 21, p. 19)
c     define:
c     ce = cos(epsilon), se = sin(epsilon)

c      parameter (CE=0.91748213149438d0)
c      parameter (SE=0.39777699580108d0)


c     n.b. r(1) remains unchanged

      tmpy = r(2)
      tmpz = r(3)

      r(2) = ceecl*tmpy + seecl*tmpz
      r(3) = -seecl*tmpy + ceecl*tmpz 

      return
      end

      

