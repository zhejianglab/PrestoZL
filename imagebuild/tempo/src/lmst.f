c      $Id$
	subroutine lmst(nmjdu1,fmjdu1, olong, xlst, sdd)
c
c	calculate local mean sidereal time xlst (turns), and the derivative
c	of local mean sidereal time wrt mean solar time = ratio of 
c	sidereal to solar time sdd (sidereal time units/solar time units).
c
c	input: 
c	nmjdu1 + fmjdu1 = UT1 MJD
c	olong = WEST longitude in degrees.
	implicit none
c	tu is the interval of time measured in Julian centuries of
c	36525 days of UT (mean solar days) elapsed since epoch 2000,
c	January 1d 12h UT.  From Astronomical Almanac, 1991, B6.
c	tu0 is tu evaluated for previous midnight
c	dtu = tu-tu0

	integer nmjdu1
	real*8 fmjdu1
	real*8 seconds_per_jc, a, b, c, d, bprime, cprime, dprime,
     +	  tu, tu0, dtu, gmst0, gst
	real*8 olong, xlst, sdd
	parameter (seconds_per_jc = 86400.d0*36525d0)
	parameter (a = 24 110.548 41d0, b = 8640 184.812 866d0,
     +	  c = 0.093 104d0, d = -6.2d-6)
	parameter (bprime = 1.d0 + b/seconds_per_jc)
	parameter (cprime = 2.d0*c/seconds_per_jc)
	parameter (dprime = 3.d0*d/seconds_per_jc)

	tu0 = (dfloat(nmjdu1-51545)+0.5d0)/3.6525d4
	dtu = fmjdu1/3.6525d4
	tu = tu0 + dtu

c	gmst0 (gst of previous midnight) in days.
	gmst0 = (a + tu0*(b + tu0*(c + tu0*d)))/86400.d0
c	could add 86400*36525 coefficient to get correct total change of 
c	sidereal time in one solar day with zero point adjustment required 
c	to get exact agreement with current formula
c	evaluated at tjdut1 = 2451544.5, dtim = 0., tu = -0.5/36525.
c	gst=(aprime + tu*(
c	1 seconds_per_jc + b + tu*(c + d*tu)))/86400.d0,
c	where aprime = a + 43200 - c dtu' - d dtu'^2,
c	and dtu' = -0.5/36525
c	however, using this formula is not a good idea, because result is 
c	very large number with subsequent subtraction of aint leading to 
c	significance loss.  Instead solve for difference of gst from gmst0.
c	first take derivative of gst formula.
!	sdd=1.002737909350795d0 + tu*(5.9006q-11 - 5.9q-15*tu)
	sdd = bprime + tu*(cprime + tu*dprime)
c	expand differences of various powers of tu and tu0.
c	note gst ~ gmst0 + dtim*sdd (old formula).
	gst=gmst0+dtu*(seconds_per_jc + b + c*(tu+tu0) +
     +    d*(tu*tu + tu*tu0 + tu0*tu0))/86400.d0
	xlst=gst-olong/360.d0
	xlst=mod(xlst,1.d0)
	if(xlst.lt.0.d0)xlst=xlst+1.d0
	return
	end
