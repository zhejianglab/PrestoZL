c      $Id$
      SUBROUTINE PRCNUT(fmjd, psid, epsd, ob, prn)

C     This sets up a matrix PRN(I,J) to precess coordinates to/from those
C     of 2000.0.  This is a modified version (DJN, July 1988) of the original
C     subroutine, which used 1950.0 coordinates.  (The old routine is now
C     called PRCNUT50.)

C     To precess coordinates of date to 2000.0, set RP(I)=PRN(I,J)*R(J);
C     Opposite transformation is RP(I)=PRN(J,I)*R(J)
c	(AWI 1995) changes:
c	independent variable is jd in TDB.
c	use call to precession and nutation to derive matrices.
c	n.b. prc and nut are transposes of the original prcnut routine
c	matrices.  This is handled properly when calculating prcnut
c	consistently with old convention.
c	psid = nutation in longitude in radians at jd
c	epsd = nutation in obliquity in radians at jd
c	ob = mean obliquity in radians

      implicit real*8 (A-H,O-Z)
      real*8 PRC(3,3),NUT(3,3),prn(3,3), dummy(3,3)
c  Compute precession matrix (corresponding to transpose of original 
c  prcnut version).
      call precession(fmjd, prc, dummy)
c  Compute nutation matrix (corresponding to transpose of original
c  prcnut version).
      call nutation(psid, epsd, ob, nut, dt)
      do i = 1,3
         do j = 1,3
c  n.b. RHS transforms J2000 positions to equator and equinox of date.
c  however, we take transpose on LHS so transforms from equator and
c  equinox of date to J2000.
	    prn(j,i) = nut(i,1)*prc(1,j) + nut(i,2)*prc(2,j) 
     +           + nut(i,3)*prc(3,j)
	enddo
      enddo
      end
