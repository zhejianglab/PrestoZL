c      $Id$
	subroutine bnryddt(ct,f0,torb,fctn)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*8 fctn(NPAP1),k,m,m1,m2,trb(3)
C	real*4 xbp(10),xbpp(10)
	parameter (TWOPI=6.28318530717958648d0,TSUN=4.925490947d-6)
	parameter (RAD=360.d0/TWOPI)
	include 'orbit.h'
C	data xbp /-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5/
C	data xbpp/-15.,-10.,-5.,0.,5.,10.,15.,10.,15.,20./

	tt0=(ct-t0(1))*86400.d0
	an=TWOPI/pb(1)
	x=a1(1)+xdot*tt0
	ecc=e(1)+edot*tt0
C  Compute eccentric anomaly u by iterating Kepler's equation.
	orbits=tt0/pb(1) - 0.5d0*pbdot*(tt0/pb(1))**2
	norbits=orbits
	if(orbits.lt.0.d0) norbits=norbits-1
	phase=TWOPI*(orbits-norbits)
	u=phase+ecc*sin(phase)*(1+ecc*cos(phase))
	fac=1/(1-ecc*cos(u))
10	du=fac*(phase-(u-ecc*sin(u)))
	u=u+du
	if(abs(du).gt.1.d-14) go to 10
	su=sin(u)
	cu=cos(u)

C  DD equations 17a, 29:
	ae=2*atan(sqrt((1+ecc)/(1-ecc))*tan(0.5d0*u))
	if(ae.lt.0.0) ae=ae+TWOPI
	if(abs(ae-phase).gt.3.14d0)  then
	  print *,ae,phase,orbits,ct
	  stop 'ae error'
	endif
	ae=TWOPI*orbits + ae - phase

C###	amz=am
	am1z=am
	am2z=am2
	deltam1=1.d-8
	deltam2=1.d-8

	do 20 i=1,4
	am1=am1z
	am2=am2z
	if(i.eq.1) am1=am1z+deltam1
	if(i.eq.2) am2=am2z+deltam2

C  Compute quantities depending on M, m2, bp, and bpp.
C###	am1=am-am2
	am=am1+am2
	m=am*TSUN
	m1=am1*TSUN
	m2=am2*TSUN				! r=m2
	aa=2.1569176d0
	bb=1.0261529d0
	cp=0.21d0
	c1=cp*am1
	c2=cp*am2
	a1a2=0.5d0*bp*bb*(c1**2 + c2**2)
	bk1=-c2-bb*c1**2 + (aa-3*bb)*c2**2 - 2*(aa-bb)*c1**2 * c2 +
     +    (2*aa**2 -7*aa*bb + 5*bb**2)*c1**2 * c2**2
	bk2=-3*c2**3 + 2*c1**2 * c2**2 + c2**4 + 0.5d0*c1**4 * c2 + 
     +    aa*c1**4 * c2**2
	a1b2a1=bp*bk1 + (bp*bb)**2 * bk2 + 0.5d0*bpp*bb*c2**2
	bk3=-c1-bb*c2**2 + (aa-3*bb)*c1**2 - 2*(aa-bb)*c2**2 * c1 +
     +    (2*aa**2 -7*aa*bb + 5*bb**2)*c2**2 * c1**2
	bk4=-3*c1**3 + 2*c2**2 * c1**2 + c1**4 + 0.5d0*c2**4 * c1 + 
     +    aa*c2**4 * c1**2
	a2b1a2=bp*bk3 + (bp*bb)**2 * bk4 + 0.5d0*bpp*bb*c1**2
	a0a2=0.5d0*bp*bb*c2**2

	brk1=(1.d0-a1a2/3.d0)*(1.d0+a1a2)**(-1.d0/3.d0) - 
     +    (1.d0/(6.d0*am)) *
     +    (am1*a1b2a1+am2*a2b1a2)*(1.d0+a1a2)**(-4.d0/3.d0)
	k=(3.d0/(1.d0-ecc**2)) * (an*TSUN*am)**(2.d0/3.d0) * brk1
	brk2=1.d0 + (am2/am)*(1.d0+a1a2)+a0a2
	gamma=(ecc/an)*(an*TSUN*am)**(2.d0/3.d0) * (am2/am) *
     +    (1.d0+a1a2)**(-1.d0/3.d0) * brk2
	si=x*an**(2.d0/3.d0) * (TSUN*am)**(-1.d0/3.d0) * (am/am2) *
     +    (1.d0+a1a2)**(-1.d0/3.d0)

	arr=(m/an**2)**(1.d0/3.d0)
	ar=arr*m2/m

C  DD equations 36, 37:
	dr=(3*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)
	er=ecc*(1+dr)
	dth=(3.5d0*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)
	eth=ecc*(1+dth)

C  DD equations 46 through 52.
	omega=omz(1)/RAD + k*ae
	sw=sin(omega)
	cw=cos(omega)
	alpha=x*sw
	beta=x*sqrt(1-eth**2)*cw
	bg=beta+gamma
	dre=alpha*(cu-er) + bg*su
	drep=-alpha*su + bg*cu
	drepp=-alpha*cu - bg*su
	onemecu=1-ecc*cu
	anhat=an/onemecu

C  DD equations 26, 27, 57:
	cume=cu-ecc
	sqr1me2=sqrt(1-ecc**2)
	brace=onemecu-si*(sw*cume+sqr1me2*cw*su)
	dlogbr=dlog(brace)
	ds=-2*m2*dlogbr

C  These will be different if spin axis not aligned:
	a0=an*ar/(TWOPI*f0*si*sqr1me2)
	b0=0.d0
	da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw)

C  Now compute d2bar, the orbital time correction in DD equation 42.
	d2bar=dre*(1-anhat*drep+(anhat**2)*(drep**2 + 0.5d0*dre*drepp -
     +    0.5d0*ecc*su*dre*drep/onemecu)) + ds + da
20	trb(i)=-d2bar
	torb=trb(3)

C  Now we need the partial derivatives. Use DD equations 62a - 62k.
	cx=sw*cume+sqr1me2*cw*su
	ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2
	csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu
	comega=x*(cw*cume-sqr1me2*sw*su)
	cm=-(trb(1)-torb)/deltam1
	cm2=-(trb(2)-torb)/deltam2

	fctn(9)=cx*f0
	fctn(10)=ce*f0
	fctn(11)=-csigma*an*f0*86400.d0
	fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
	fctn(13)=comega*f0
	fctn(18)=tt0*fctn(12)*0.5d-6
	fctn(21)=cm*f0
	fctn(22)=cm2*f0
	fctn(24)=tt0*fctn(9)
	fctn(25)=tt0*fctn(10)

	am=am1
	return
	end
