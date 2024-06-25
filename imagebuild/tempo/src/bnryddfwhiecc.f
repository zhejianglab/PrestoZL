c      $Id$
        subroutine bnryddfwhiecc(torb,fctn)
c   modified version of bnrydd (Damour-Deruelle), with Shapiro
c   terms re-parameterized as per Freire & Wex 2010, for
c   high eccentricity orbits.  Joel Weisberg, Yuping Huang, & Andrew Chael.
c   See Weisberg & Huang APJ 2016 for first use.
c   requires inputs of h3 and varsigma instead 
c   of m2 (aka r) and sini (aka s)

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

c	an is mean motion "n" in radians/s
	an=twopi/pb(1)
c	k is the MEAN input omdot (deg/yr), conv to  rad/orbit
c	as described below eq 16c in DD
	k=omdot/(RAD*365.25d0*86400.d0*an)
c	si, aka sin i, was taken in bnrydd (orig) from input. Now it is derived from inverting
c	the input FW parameter called varsigma (invert FW eq 12):
        varsigmasqp1 = varsigma**2.0 + 1.0

c	tt0 is elapsed time (s) since t0(1)
	tt0=(ct-t0(1))*86400.d0

c	note that x is a1(1) updated to epoch via input xdot
	xx=a1(1)+xdot*tt0
c	note that ecc is e(1) updated to epoch via input edot
	ecc=e(1)+edot*tt0
c	er is defined in terms of input "dr" aka delta_r (Eq 36)
	er =ecc*(1.d0+dr)
c	eth is defined in terms of input "dth" aka delta_theta (Eq 37)
	eth=ecc*(1.d0+dth)
		
c	orbits is orbits elapsed since tt0, including (pbdot+xpbdot)
	orbits=tt0/pb(1) - 0.5d0*(pbdot+xpbdot)*(tt0/pb(1))**2
c	norbits is next lowest integer part of orbits
	norbits=orbits
	if(orbits.lt.0.d0) norbits=norbits-1
c	phase is fractional part of orbit in radians, (0->2 pi). ie
c	it is mean anomaly in rad
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
c	Ae is true anomaly; cae and sae its cos and sine, see DD 17 a,b,c
	cae=(cu-ecc)/onemecu
	sae=sqrt(1.d0-ecc**2)*su/onemecu
	ae=atan2(sae,cae)
	if(ae.lt.0.0) ae=ae+twopi
c	add (2 pi *norbits)to ae to get total elapsed true anomaly since tt0:
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
c         need to calculate si from varsigma and h3 now
          si = 2.d0*varsigma/varsigmasqp1
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
c	DD 46, 47:
	alpha=xx*sw
	beta=xx*sqrt(1-eth**2)*cw
	bg=beta+gamma
c	DD eq 48 Delta_RomerEinstein:
	dre=alpha*(cu-er) + bg*su
c	DD eq 49 for Delta_RE':
	drep=-alpha*su + bg*cu
c	DD eq 50 for Delta_RE'':
	drepp=-alpha*cu - bg*su
c	DD eq 51 for nhat:
	anhat=an/onemecu
c  DD equations 26, 27, 57:
	sqr1me2=sqrt(1-ecc**2)
	cume=cu-ecc
	brace=onemecu-(2.d0*varsigma / varsigmasqp1) *(sw*cume+sqr1me2*cw*su)
	dlogbr=dlog(brace)
c	DD Eq 60 for Delta_Shapiro now substituted by FW2010 expression
	ds=-2.d0*h3/(varsigma**3.d0)*dlogbr
c	DD eq 27 for Delta_Aberration
	da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw)
c  Now compute d2bar, the orbital time correction in DD equation 42.
c  Delta double bar, as given by DD eq 52:
	d2bar=dre*(1-anhat*drep+(anhat**2)*(drep**2 + 0.5d0*dre*drepp -
     +    0.5d0*ecc*su*dre*drep/onemecu)) + ds + da
c	and then torb inverts it (see DD42):
	torb=-d2bar

c  Now we need the partial derivatives. Use DD equations 62a - 62k.
c  For sign convention, note that these are partials of the Deltas, not of the torbs.
	csigma=xx*(-sw*su+sqr1me2*cw*cu)/onemecu
	ce=su*csigma-xx*sw-ecc*xx*cw*su/sqr1me2
	cx=sw*cume+sqr1me2*cw*su
	comega=xx*(cw*cume-sqr1me2*sw*su)
	cgamma=su
	cdth=-ecc*ecc*xx*cw*su/sqr1me2
c   the following are the partials for h3 and varsigma
        Ch3 = (-2*dlog(onemecu - (2*varsigma*(cw*sqr1me2*su + cume*sw))/
     -        (1 + varsigma**2)))/varsigma**3
        Cvarsigma = (2*h3*((-2*varsigma*(-1 + varsigma**2)*(cw*sqr1me2* 
     -  su+cume*sw))/((1 + varsigma**2)*((1 + varsigma**2)*onemecu - 
     -  2*varsigma*(cw*sqr1me2*su + cume*sw))) + 
     -  3*dlog(onemecu - (2*varsigma*(cw*sqr1me2*su + cume*sw))/
     -  (1 + varsigma**2))))/varsigma**4
c	need to convert partials' units, usually to per spin period, by mult by f0:
	fctn(9)=cx*f0
	fctn(10)=ce*f0
	fctn(11)=-csigma*an*f0*86400.d0
	fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
	fctn(13)=comega*f0
	fctn(14)=ae*fctn(13)/(an*RAD*365.25d0*86400.d0)
	fctn(15)=cgamma*f0
	fctn(18)=0.5d-6*tt0*fctn(12)
	fctn(20)=Cvarsigma*f0
	fctn(22)=Ch3*f0
	fctn(23)=cdth*f0
	fctn(24)=tt0*fctn(9)
	fctn(25)=tt0*fctn(10)

	return
	end


