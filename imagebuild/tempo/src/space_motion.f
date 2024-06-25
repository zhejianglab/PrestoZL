c      $Id$
c	$Id$
c	$Date$
c	%R%
c	%Y%
	subroutine space_motion(alpha, delta, mualpha, mudelta, parallax, 
     +    v, r, rdot, ifangle_input)
c	if ifangle_input ge 1, then calculate space motion vectors r, and rdot
c	from alpha, delta, mualpha, mudelta, parallax, and v.
c	if ifangle_input lt 1 then calculate alpha, delta, etc. from r and rdot.
c	space motion vectors r(1)...r(3) = x, y, z, rdot(1)...rdot(3) = 
c	xdot, ydot, zdot.
c	n.b. all r and rdot coordinates are divided by a standard radial 
c	distance |r|, if ifangle_input = 1, thus,
c		r(1) to r(3) are dimensionless and rdot(1)....rdot(3) have units
c		of radians per century and measure xdot, etc. divided by the
c		standard distance.
c	n.b.	if ifangle_input gt 1, then r and rdot are not scaled by
c		standard radial distance, and units are au's for distance
c		and au/century for velocity.
c	alpha is the ra in radians
c	delta is the dec in radians
c	mualpha is the proper motion (change in ra) in radians of arc per
c		century.
c	mudelta is the proper motion (change in dec) in radians of arc per
c		century.
c	parallax is the parallax in radians.
c	v is the velocity in au per century.
c	A.E. p. B39.
	implicit none
	real*8 mualpha, mudelta
	real*8 r(3), rdot(3), pi, alpha, delta, parallax, v,
     +     ca, sa, cd, sd, dot, vparallax
	integer*4 ifangle_input, i
	parameter(pi=3.1415926535897932d0)
	if(ifangle_input.ge.1) then
		ca = cos(alpha)
		sa = sin(alpha)
		cd = cos(delta)
		if(cd.lt.0.d0) stop 'space_motion: |delta| > pi/2'
		sd = sin(delta)
		r(1) = ca*cd
		r(2) = sa*cd
		r(3) = sd
		vparallax = v*parallax
		rdot(1) = -mualpha*sa*cd - mudelta*ca*sd + vparallax*r(1)
		rdot(2) = mualpha*ca*cd - mudelta*sa*sd + vparallax*r(2)
		rdot(3) = mudelta*cd + vparallax*r(3)
		if(ifangle_input.gt.1) then
c			multiply through by distance in au's or divide by the
c			parallax in radians.
			do 10 i = 1,3
				r(i) = r(i)/parallax
10				rdot(i) = rdot(i)/parallax
		endif
	else
		parallax = 1.d0/sqrt(dot(r, r))
		delta = asin(parallax*r(3))
		alpha = atan2(r(2),r(1))
		if(alpha.lt.0.d0) alpha = alpha + 2.d0*pi
		ca = cos(alpha)
		sa = sin(alpha)
		cd = cos(delta)
		if(cd.le.0.d0) stop 'space_motion: cos dec .le. 0.d0'
		sd = sin(delta)
		mualpha = parallax*(-sa*rdot(1) + ca*rdot(2))/cd
		mudelta = parallax*(-ca*sd*rdot(1) - sa*sd*rdot(2) + cd*rdot(3))
		v = rdot(1)*ca*cd + rdot(2)*sa*cd + rdot(3)*sd
	endif
	end
