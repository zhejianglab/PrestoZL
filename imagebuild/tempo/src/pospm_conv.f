c      $Id$
	SUBROUTINE POSPM_CONV(alpha, delta, pmra, pmdec, r, u)

c       INPUT:	alpha: RA  in radians
c	        delta: DEC in radians
c	        pmra: proper motion in RA  in radians of arc per century
c	        pmdec: proper motion in DEC in radians of arc per century
c       OUTPUT: r(3): unit vector in the direction of the pulsar
c               u(3): vector in the direction of the transverse velocity of 
c                     the pulsar in radians of arc per century
c
c       N.W. 8/1997

	implicit none
	real*8 alpha,delta,pmra, pmdec
	real*8 r(3), u(3)
	real*8 ca,sa,cd,sd

	ca = DCOS(alpha)
	sa = DSIN(alpha)
	cd = DCOS(delta)
	sd = DSIN(delta)
	if(cd.lt.0.d0) stop 'space_motion2: |delta| > pi/2'

	r(1) = ca*cd
	r(2) = sa*cd
	r(3) = sd

	u(1) = -pmra*sa*cd - pmdec*ca*sd
	u(2) =  pmra*ca*cd - pmdec*sa*sd
	u(3) =  pmdec*cd

	RETURN
	END
