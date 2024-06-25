c      $Id$
	subroutine precession(fmjd, precess, dprecess)

c  calculate the precession matrix precess(3,3), and its derivative
c  wrt time in days dprecess(3,3) as a function of the julian
c  day number.  It is assumed the original coordinates being precessed
c  to rjd are epoch J2000.
c
c  from 1984 AA S19 with T = 0.d0.  note there are more significant
c  figures in these expressions than in the more standard AA B18.
c  These constants were originally given by 
c  lieske 1979 a&a 73, 282 and have been checked to give his
c  his precession matrix (also quoted in
c  Standish (1982 a&a 115, 20)) to ~16 figures.

	implicit none
	real*8 fmjd
	real*8 precess(3,3), dprecess(3,3),
     +     par_zeta(3), par_z(3), par_theta(3),
     +     pi, seconds_per_rad, t, czeta, zeta, dzeta, z, dz, theta,
     +     dtheta, szeta, cz, sz, dctheta, dczeta, dszeta, dcz, dsz, 
     +     ctheta, stheta, dstheta
	parameter (pi = 3.14159 26535 89793d0)
	parameter (seconds_per_rad = 3600.d0*180.d0/pi)
	data par_zeta /2306.2181d0, 0.30188d0, 0.017998d0/
	data par_z /2306.2181d0, 1.09468d0, 0.018203d0/
	data par_theta /2004.3109d0, -0.42665d0, -0.041833d0/

c  AA, B18 alternatives (not used).
c       data par_zeta /0.640 6161d0, 8.39d-5, 5.d-6/
c	data par_z /0.640 6161d0, 3.041d-4, 5.1d-6/
c	data par_theta /0.556 7530d0, -1.185d-4, -1.16d-5/
c  t is the change in time from J2000 in Julian centuries 
c  (J2000 = JD 2451545.0 = MJD 51544.5)

        t=(fmjd-51544.5d0)/36525.d0

        zeta=t*(par_zeta(1)+t*(par_zeta(2)+t*par_zeta(3)))/
     +     seconds_per_rad
	dzeta=(par_zeta(1)+t*(2.d0*par_zeta(2)+t*3.d0*par_zeta(3)))/
     +     seconds_per_rad/36525.d0
	z=t*(par_z(1)+t*(par_z(2)+t*par_z(3)))/
     +     seconds_per_rad
	dz=(par_z(1)+t*(2.d0*par_z(2)+t*3.d0*par_z(3)))/
     +     seconds_per_rad/36525.d0
	theta=t*(par_theta(1)+t*(par_theta(2)+t*par_theta(3)))/
     +     seconds_per_rad
	dtheta=(par_theta(1)+t*(2.d0*par_theta(2)
     +     +t*3.d0*par_theta(3)))/seconds_per_rad/36525.d0
     
	czeta = cos(zeta)
	szeta = sin(zeta)
	dczeta = -szeta*dzeta
	dszeta = czeta*dzeta
	cz = cos(z)
	sz = sin(z)
	dcz = -sz*dz
	dsz = cz*dz
	ctheta = cos(theta)
	stheta = sin(theta)
	dctheta = -stheta*dtheta
	dstheta = ctheta*dtheta

        precess(1,1) = czeta*ctheta*cz - szeta*sz
	precess(2,1) = czeta*ctheta*sz + szeta*cz
	precess(3,1) = czeta*stheta
	precess(1,2) = -szeta*ctheta*cz - czeta*sz
	precess(2,2) = -szeta*ctheta*sz + czeta*cz
	precess(3,2) = -szeta*stheta
	precess(1,3) = - stheta*cz
	precess(2,3) = -stheta*sz
	precess(3,3) = ctheta

        dprecess(1,1) = dczeta*ctheta*cz + czeta*dctheta*cz +
     +    czeta*ctheta*dcz - dszeta*sz - szeta*dsz
	dprecess(2,1) = dczeta*ctheta*sz + czeta*dctheta*sz +
     +    czeta*ctheta*dsz + dszeta*cz + szeta*dcz
	dprecess(3,1) = dczeta*stheta + czeta*dstheta
	dprecess(1,2) = - dszeta*ctheta*cz - szeta*dctheta*cz 
     +    - szeta*ctheta*dcz - dczeta*sz - czeta*dsz
	dprecess(2,2) = - dszeta*ctheta*sz - szeta*dctheta*sz 
     +    - szeta*ctheta*dsz + dczeta*cz + czeta*dcz
	dprecess(3,2) = - dszeta*stheta - szeta*dstheta
	dprecess(1,3) = - dstheta*cz - stheta*dcz
	dprecess(2,3) = - dstheta*sz - stheta*dsz
	dprecess(3,3) = dctheta

        end
