c      $Id$
c	$Id$
c	$Date$
c	%R%
c	%Y%
	subroutine nutation(delp, dele, eps, dnutat, dt)
c	use jpl ephemeris prediction of nutation.  this is special
c	tempo version that relies on outside call of JPL ephemeris to
c	provide delp, dele, and eps.
c	This result
c	has been checked extensively against previous nutate (written by
c	(Stephenson Yang) (involving huge IAU series).
c	result:  Got agreement to 10^-5 arcsec.
c	nutation matrix in dnutat.  x'(i) = dnutat(i,j) * x(j) transforms 
c	x (coordinates for mean equinox and equator of date) to x' (coordinates
c	of actual equinox and equator of date).
c	dt corresponds to equation of equinoxes in radians = delp*cos(eps+dele),
c	delp corresponds to the nutation in longitude in radians
c	dele corresponds to the nutation in obliquity in radians
c	eps corresponds to the mean obliquity in radians
	implicit none
	real*8 dnutat(3,3), dt, eps, ceps, seps, delp, dele
	ceps=cos(eps)
	seps=sin(eps)
c	this is only a first order correction according to AA explanatory  
c	Not clear if should use mean obliquity (eps) or true obliquity here, 
c	but we use mean obliquity.
	dnutat(1,1)=1.d0
	dnutat(1,2)=-delp*ceps
	dnutat(1,3)=-delp*seps
	dnutat(2,1)=-dnutat(1,2)
	dnutat(2,2)=1.d0
	dnutat(2,3)=-dele
	dnutat(3,1)=-dnutat(1,3)
	dnutat(3,2)=-dnutat(2,3)
	dnutat(3,3)=1.d0
c	this agrees with tables, but not clear that AA tables are correct.
	dt = delp*cos(eps+dele)
	end
