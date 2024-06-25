c      $Id$
	subroutine bnrydd(torb,fctn)

c  Damour et Deruelle modele pour le chronometrage des temps d'arrive 
c  au barycentre du systeme solaire, a la premiere approximation
c  post-newtonienne de la Relativite Generale.

c  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch.
c  Pulsar proper time is then TP=T+TORB.
c  Units are such that c=g=1. Thus masses have units of seconds, with
c  one solar mass = 4.925490947 usec.

c  Also computes the binary orbit-related values of fctn: partial
c  derivatives of each arrival time residual with respect to the model
c  parameters.

c  Initial guesses for all parameters must be placed in common/orbit/ by the
c  calling program.  

c  Variable name 'x' changed to 'xx' to avoid conflict with variable
c     'x' defined inacom.h

	implicit real*8 (a-h,o-z)
	include 'dim.h'
        include 'acom.h'
	real*8 fctn(NPAP1),k,m2
	parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
	parameter (RAD=360.d0/twopi)
	include 'dp.h'
	include 'orbit.h'
      common /CRDS/ RCB(6),RBE(6),RCE(6),RSE(3),RCA(3),RSA(3),RBA(3),
     +             REA(3),PSID,EPSD,PC,PS,TSID,PRN(3,3),ATUT,UT1UT,ETAT,
     +             DTGR,TDIS,BCLT      
	include 'trnsfr.h'

	an=twopi/pb(1)
	k=omdot/(RAD*365.25d0*86400.d0*an)
	m2=am2*SUNMASS

	tt0=(ct-t0(1))*86400.d0

	xx=a1(1)+xdot*tt0
	ecc=e(1)+edot*tt0
	er =ecc*(1.d0+dr)
	eth=ecc*(1.d0+dth)
		
	orbits=tt0/pb(1) - 0.5d0*(pbdot+xpbdot)*(tt0/pb(1))**2
	norbits=orbits
	if(orbits.lt.0.d0) norbits=norbits-1
	phase=twopi*(orbits-norbits)

c  Compute eccentric anomaly u by iterating Keplers equation.
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

        if (useannorb) then
          deltai0 = -rca(1)*sin(pra) + rca(2)*cos(pra)
          deltaj0 = (-rca(1))*sin(pdec)*cos(pra) +
     +              (-rca(2))*sin(pdec)*sin(pra) +
     +              ( rca(3))*cos(pdec)
          if (usefixeddist) then  ! use fixed distance in kpc
            dist = 499.00478364D0/(twopi/fixeddist/1000/3600/360)
          else                    ! use parallax distance
            dist = 499.00478364D0/(twopi*px/1000/3600/360)
          endif
          Omkopeikin = twopi/4-(PAAscNode/360*twopi)
          omegax = -(1.d0/si)/dist*
     +           (deltai0*cos(Omkopeikin)+deltaj0*sin(Omkopeikin))
          coti = -sqrt(1/si**2-1)  ! hack -- hard-wired for negative cot i
          xxterm = coti/dist * 
     +           (deltai0*sin(Omkopeikin)-deltaj0*cos(Omkopeikin))
          xx = xx*(1.d0+xxterm)
          omega=omz(1)/RAD + k*ae + omegax
        else
          omega=omz(1)/RAD + k*ae 
        endif
	sw=sin(omega)
	cw=cos(omega)
	alpha=xx*sw
	beta=xx*sqrt(1-eth**2)*cw
	bg=beta+gamma
	dre=alpha*(cu-er) + bg*su
	drep=-alpha*su + bg*cu
	drepp=-alpha*cu - bg*su
	anhat=an/onemecu
c  DD equations 26, 27, 57:
	sqr1me2=sqrt(1-ecc**2)
	cume=cu-ecc
	brace=onemecu-si*(sw*cume+sqr1me2*cw*su)
	dlogbr=dlog(brace)
	ds=-2*m2*dlogbr
	da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw)
c  Now compute d2bar, the orbital time correction in DD equation 42.
	d2bar=dre*(1-anhat*drep+(anhat**2)*(drep**2 + 0.5d0*dre*drepp -
     +    0.5d0*ecc*su*dre*drep/onemecu)) + ds + da
	torb=-d2bar

c  Now we need the partial derivatives. Use DD equations 62a - 62k.
	csigma=xx*(-sw*su+sqr1me2*cw*cu)/onemecu
	ce=su*csigma-xx*sw-ecc*xx*cw*su/sqr1me2
	cx=sw*cume+sqr1me2*cw*su
	comega=xx*(cw*cume-sqr1me2*sw*su)
	cgamma=su
	cdth=-ecc*ecc*xx*cw*su/sqr1me2
	cm2=-2*dlogbr
	csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace

	fctn(9)=cx*f0
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
	fctn(24)=tt0*fctn(9)
	fctn(25)=tt0*fctn(10)

	return
	end


