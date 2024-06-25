c      $Id$
	subroutine bnryddgr(ct,f0,n,torb,fctn)

c  Damour et Deruelle modele pour le chronometrage des temps d'arrive 
c  au barycentre du systeme solaire, a la premiere approximation
c  post-newtonienne de la Relativite Generale.

c  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch.
c  Pulsar proper time is then TP=T+TORB
c  Units are such that c=g=1. Thus masses have units of seconds, with
c  one solar mass = 4.925 490 947d-6 sec.

c  Also computes the binary orbit-related values of fctn: partial
c  derivatives of each arrival time residual with respect to the model
c  parameters.

c  Initial guesses for all parameters must be placed in common/orbit/ by the
c  calling program.  

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*8 fctn(NPAP1),k,m,m1,m2
	parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
	parameter (RAD=360.d0/twopi)
	parameter (ARRTOL=1.d-10)
	include 'orbit.h'

	tt0=(ct-t0(1))*86400.d0
	an=twopi/pb(1)
	x=a1(1)+xdot*tt0
	ecc=e(1)+edot*tt0
	m=am*SUNMASS
	m2=am2*SUNMASS
	m1=m-m2
	call mass2dd(am,am2,x,ecc,an,arr,ar,xk,si,gamma,pbdot)
	k=xk

c  DD equations 36, 37:
	dr=(3*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)
	er=ecc*(1+dr)
	dth=(3.5d0*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)
	eth=ecc*(1+dth)

c  Compute eccentric anomaly u by iterating Kepler's equation.

	orbits=tt0/pb(1) - 0.5d0*(pbdot+xpbdot)*(tt0/pb(1))**2
	norbits=orbits
	if(orbits.lt.0.d0) norbits=norbits-1
	phase=twopi*(orbits-norbits)
	u=phase+ecc*dsin(phase)*(1+ecc*dcos(phase))
10	fac=1/(1-ecc*dcos(u))
	du=fac*(phase-(u-ecc*dsin(u)))
	u=u+du
	if(dabs(du).gt.1.d-14) go to 10

C  DD equations 17a, 29:
	ae=2*datan(dsqrt((1+ecc)/(1-ecc))*dtan(0.5d0*u))
	if(ae.lt.0.0) ae=ae+twopi
	if(dabs(ae-phase).gt.3.14d0)  then
		print *,ae,phase,orbits,ct
		stop 'ae error'
		endif
	ae=twopi*orbits + ae - phase
	omega=omz(1)/RAD + (k + xomdot/(an*RAD*365.25d0*86400.d0))*ae

C  DD equations 46 through 52.
	su=dsin(u)
	cu=dcos(u)
	sw=dsin(omega)
	cw=dcos(omega)
	alpha=x*sw
	beta=x*dsqrt(1-eth**2)*cw
	bg=beta+gamma
	dre=alpha*(cu-er) + bg*su
	drep=-alpha*su + bg*cu
	drepp=-alpha*cu - bg*su
	onemecu=1-ecc*cu
	anhat=an/onemecu
c  DD equations 26, 27, 57:
	cume=cu-ecc
	sqr1me2=dsqrt(1-ecc**2)
	brace=onemecu-si*(sw*cume+sqr1me2*cw*su)
	dlogbr=dlog(brace)
	ds=-2*m2*dlogbr
C  These will be different if spin axis not aligned:
	a0aligned=an*ar/(twopi*f0*si*sqr1me2)
	a0=afac*a0aligned
	b0=0.d0
	da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw)
c  Now compute d2bar, the orbital time correction in DD equation 42.
	d2bar=dre*(1-anhat*drep+(anhat**2)*(drep**2 + 0.5d0*dre*drepp -
     +    0.5d0*ecc*su*dre*drep/onemecu)) + ds + da
	torb=-d2bar

	tmp=dre*(1-anhat*drep+(anhat**2)*(drep**2 + 0.5d0*dre*drepp -
     +    0.5d0*ecc*su*dre*drep/onemecu)) 

c  Now we need the partial derivatives. Use DD equations 62a - 62k.
	csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu
	ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2
	cx=sw*cume+sqr1me2*cw*su
	comega=x*(cw*cume-sqr1me2*sw*su)
      dtdm2=-2*dlogbr
      an0=dsqrt(m/arr**3)
      dnum=an0*(m-2*m2)/(2*arr*m)
      fact1=(m/(2*arr)) * ((m-m2)*m2/m**2 - 9)
      fact2=(3*m/(2*arr**4)) * (1.d0 + fact1)
      denom=an0*fact1/arr + fact2/an0
      darrdm2=dnum/denom

      cgamma=su
      dgmdm2=( (m+2*m2)/arr - (m2*(m+m2)*darrdm2/arr**2) )*ecc/
     1                          (an*m)
      csini=2*m2*(sw*cume+sqr1me2*cw*su)/brace
      dsidm2=-(m*x/(arr*m2)) * (1.d0/m2 + darrdm2/arr)
      ck=ae*comega
      dkdm2=-k*darrdm2/arr
      cdr=-ecc*x*sw
      ddrdm2=-dr*darrdm2/arr - 2*m2/(arr*m)
      cdth=-ecc**2*x*cw*su/sqr1me2
      dthdm2=-dth*darrdm2/arr - (m+m2)/(arr*m)
      cpbdot=-csigma*an*tt0**2/(2*pb(1))
      dpbdm2=pbdot/m2 - pbdot/(m-m2)
      cm2=dtdm2+cgamma*dgmdm2+csini*dsidm2+ck*dkdm2+
     1 cdr*ddrdm2+cdth*dthdm2+cpbdot*dpbdm2

      fact3=(m/2*arr) * (m2/m**2-2*(m-m2)*m2/m**3)
      denumm=(1+fact1)/(2*arr**3*an0) + an0*(fact1/m+fact3)
      fact4=(1+fact1)*3*m/(2*arr**4*an0)
      fact5=an0*fact1/arr
      denomm=fact4+fact5
      darrdm=denumm/denomm
      dkdm=k/m - k*darrdm/arr
      fact6=1.d0/(arr*m)
      fact7=-(m+m2)/(arr*m**2)
      fact8=-(m+m2)*darrdm/(arr**2*m)
      dgamdm=(ecc*m2/an)*(fact6+fact7+fact8)
      ddrdm=-dr/m -dr*darrdm/arr + 6/arr
      dthdm=-dth/m - dth*darrdm/arr + (7*m-m2)/(arr*m)
      dpbdm=pbdot/(m-m2) - pbdot/(3*m)
      dsidm=-(m*x/(arr*m2))*(-1.d0/m+darrdm/arr)
      cm=ck*dkdm+cgamma*dgamdm+cdr*ddrdm+cdth*dthdm+cpbdot*dpbdm+
     1     csini*dsidm

	fctn(9)=cx*f0
	fctn(10)=ce*f0
	fctn(11)=-csigma*an*f0*86400.d0
	fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
	fctn(13)=comega*f0
c following used to be fctn(14), omdot, but it is really xomdot.
	fctn(37)=ae*fctn(13)/(an*RAD*365.25d0*86400.d0)  
c following used to be fctn(18), pbdot, but it is really xpbdot
	fctn(38)=tt0*fctn(12)*0.5d-6
	fctn(21)=cm*f0*SUNMASS
	fctn(22)=cm2*f0*SUNMASS
	fctn(24)=tt0*fctn(9)
	fctn(25)=tt0*fctn(10)

	return
	end
