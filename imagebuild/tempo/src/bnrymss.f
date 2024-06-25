c      $Id$
	subroutine bnrymss(torb,fctn)

c  Model for main-sequence star binary pulsars (Wex 1998, astro-ph/9706086).
c  -> changes in x proportonal to Ae(u)
c  -> use of second time derivative in omega and x: om2dot, x2dot

c  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch.
c  Pulsar proper time is then TP=T+TORB.
c  Units are such that c=g=1. Thus masses have units of seconds, with
c  one solar mass = 4.925490947 usec.

c  Also computes the binary orbit-related values of fctn: partial
c  derivatives of each arrival time residual with respect to the model
c  parameters.

c  Initial guesses for all parameters must be placed in common/orbit/ by the
c  calling program.  

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*8 fctn(NPAP1),k,m2
	parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
	parameter (RAD=360.d0/twopi)
	include 'dp.h'
	include 'orbit.h'

	an=twopi/pb(1)
	ecc0=e(1)
	x0=a1(1)
	omega0=omz(1)/RAD
	k=omdot/an/(RAD*365.25d0*86400.d0)
	xi=xdot/an
	m2=am2*SUNMASS

	tt0=(ct-t0(1))*86400.d0
	orbits=tt0/pb(1) - 0.5d0*pbdot*(tt0/pb(1))**2
	norbits=orbits
	if(orbits.lt.0.d0) norbits=norbits-1
	phase=twopi*(orbits-norbits)

	ecc=ecc0 + edot*tt0
	er =ecc*(1.d0+dr)
	eth=ecc*(1.d0+dth)

C  Compute eccentric anomaly u by iterating Kepler's equation.
	u=phase+ecc*sin(phase)*(1+ecc*cos(phase))
10	du=(phase-(u-ecc*sin(u)))/(1.d0-ecc*cos(u))
	u=u+du
	if(dabs(du).gt.1.d-12) go to 10

c  DD equations 17b, 17c, 29, and 46 through 52
	su=sin(u)
	cu=cos(u)
	onemecu=1.d0-ecc*cu
	cae=(cu-ecc)/onemecu
	sae=sqrt(1.d0-ecc**2)*su/onemecu
	ae=atan2(sae,cae)
	if(ae.lt.0.0) ae=ae+twopi
	ae=twopi*orbits + ae - phase
	omega = omega0 +  k*ae + 0.5d0*om2dot*tt0**2 ! Wex 1998
	x     = x0     + xi*ae + 0.5d0* x2dot*tt0**2 ! Wex 1998
	sw=sin(omega)
	cw=cos(omega)
	alpha=x*sw
	beta=x*sqrt(1-eth**2)*cw
	bg=beta+gamma
	dre=alpha*(cu-er) + bg*su
	drep=-alpha*su + bg*cu
	drepp=-alpha*cu - bg*su
	anhat=an/onemecu

c  DD equations 26, 27, 57
	sqr1me2=sqrt(1-ecc**2)
	cume=cu-ecc
	brace=onemecu-si*(sw*cume+sqr1me2*cw*su)
	dlogbr=dlog(brace)
	ds=-2*m2*dlogbr
	da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw)

c  Now compute d2bar, the orbital time correction in DD equation 42
	d2bar=dre*(1-anhat*drep+(anhat**2)*(drep**2 + 0.5d0*dre*drepp -
     +    0.5d0*ecc*su*dre*drep/onemecu)) + ds + da
	torb=-d2bar

c  Partial derivatives, DD equations 62a - 62k
	csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu
	ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2
	cx=sw*cume+sqr1me2*cw*su
	comega=x*(cw*cume-sqr1me2*sw*su)
	cgamma=su
	cdth=-ecc*ecc*x*cw*su/sqr1me2
	cm2=-2*dlogbr
	csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace

	fctn( 9)=cx*f0
	fctn(10)=ce*f0
	fctn(11)=-csigma*an*f0*86400.d0
	fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
	fctn(13)=comega*f0
	fctn(14)=ae*fctn(13)/(an*RAD*365.25d0*86400.d0)
	fctn(15)=cgamma*f0
	fctn(18)=0.5d-6*tt0*fctn(12)
	fctn(20)=csi*f0
	fctn(22)=cm2*f0*SUNMASS
	fctn(23)=cdth*f0
	fctn(24)=ae*fctn(9)/an
	fctn(25)=tt0*fctn(10)
	fctn(39)=0.5d0*tt0**2*fctn(13) ! Wex 1998
	fctn(40)=0.5d0*tt0**2*fctn(9)  ! Wex 1998

	return
	end
