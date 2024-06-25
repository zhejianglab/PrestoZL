c      $Id$
	subroutine bnrydds(torb,fctn)

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

c  Changes by Wex, April 2014
c  - Inversion of timing model done by iteration to achieve maximum precision
c    Absolute must, if retardation becomes relevant! DD inversion not accurate enough.
 
c  Changes by Wex, July 2015 -- still experimental and need verification
c  - Include higher order corrections in the Shapiro delay
c    a) 1.5pN corrections
c    b) lensing delay (leading order)
c    c) rotational  beding delay
c    d) latitudinal bending delay

c  Changes by Wex, September 2015 -- still experimental and need verification
c  - Scaling factor for a+b+c+d (parameter SHAPHOF) that can be fitted for.

c  Variable name 'x' changed to 'xx' to avoid conflict with variable
c     'x' defined in acom.h

	implicit real*8 (a-h,o-z)
	include 'dim.h'
        include 'acom.h'
	real*8 fctn(NPAP1),k,m2
	parameter (twopi=6.28318530717958648d0,SUNMASS=4.925490947d-6)
	parameter (RAD=360.d0/twopi)
	include 'dp.h'
	include 'orbit.h'
      common /CRDS/ RCB(6),RBE(6),RCE(6),RSE(3),RCA(3),RSA(3),RBA(3),
     +             REA(3),PSID,EPSD,PC,PS,TSID,PRN(3,3),ATUT,UT1UT,ETAT,
     +             DTGR,TDIS,BCLT      
	include 'trnsfr.h'

        real*10 frb, tt0, tt, orbits

        frb = 1.d0/pb(1)
	an=twopi/pb(1)
	k=omdot/(RAD*365.25d0*86400.d0*an)
	m2=am2*SUNMASS

	tt0=(ct-t0(1))*86400.d0

	xx=a1(1)+xdot*tt0
	ecc=e(1)+edot*tt0
	er =ecc*(1.d0+dr)
	eth=ecc*(1.d0+dth)
		
c  --> Inversion of timing model by iteration: begin of loop
        epsNum = 1.0d-10
        delta  = 0.0d0
 10     delta_old = delta        
        tt  = tt0 - delta
        orbits  = tt*frb - 0.5d0*(pbdot+xpbdot)*(tt*frb)**2
        norbits = orbits
        if(orbits.lt.0.d0) norbits = norbits - 1
        phase = twopi*(orbits - norbits)

c  Compute eccentric anomaly u by iterating Keplers equation.
	u=phase+ecc*dsin(phase)*(1+ecc*dcos(phase))
100	du=(phase-(u-ecc*dsin(u)))/(1.d0-ecc*dcos(u))
	u=u+du
	if(dabs(du).gt.1.d-14) go to 100

c  DD equations 17b, 17c, 29, and 46 through 52
	su=dsin(u)
	cu=dcos(u)
	onemecu=1.d0-ecc*cu
        cume    = cu - ecc

	cae=(cu-ecc)/onemecu
	sae=dsqrt(1.d0-ecc**2)*su/onemecu
	ae1=datan2(sae,cae)
	if(ae1.lt.0.0) ae1=ae1+twopi
	ae=twopi*orbits + ae1 - phase

        if (useannorb) then
          deltai0 = -rca(1)*dsin(pra) + rca(2)*dcos(pra)
          deltaj0 = (-rca(1))*dsin(pdec)*dcos(pra) +
     +              (-rca(2))*dsin(pdec)*dsin(pra) +
     +              ( rca(3))*dcos(pdec)
          if (usefixeddist) then  ! use fixed distance in kpc
            dist = 499.00478364D0/(twopi/fixeddist/1000/3600/360)
          else                    ! use parallax distance
            dist = 499.00478364D0/(twopi*px/1000/3600/360)
          endif
          Omkopeikin = twopi/4-(PAAscNode/360*twopi)
          omegax = -(1.d0/si)/dist*
     +           (deltai0*dcos(Omkopeikin)+deltaj0*dsin(Omkopeikin))
          coti = -dsqrt(1/si**2-1)  ! hack -- hard-wired for negative cot i
          xxterm = coti/dist * 
     +           (deltai0*dsin(Omkopeikin)-deltaj0*dcos(Omkopeikin))
          xx = xx*(1.d0+xxterm)
          omega=omz(1)/RAD + k*ae + omegax
        else
          omega=omz(1)/RAD + k*ae 
        endif
	sw=dsin(omega)
	cw=dcos(omega)

        psi  = omega + ae1 ! angle w.r.t. ascending node
        spsi = DSIN(psi)
        cpsi = DCOS(psi)

c  Roemer delay (DD)
        alpha = xx*sw
        beta  = xx*DSQRT(1.d0 - eth**2)*cw
        dRoe  = alpha*(cu - er) + beta*su

c  Einstein delay (DD)
        dEin  = gamma * su

c  Shapiro delay
        sidds  = 1.d0 - DEXP(-1.d0*shapmax) ! see Kramer et al. 2006, AnP

        sqr1me2 = DSQRT(1.d0 - ecc**2)

        brace  = onemecu - sidds*(sw*cume + sqr1me2*cw*su)
        dlogbr = DLOG(brace)
        dSha   = -2.d0*m2*dlogbr

        if(nshapho.eq.1)then ! NW: higher order Shapiro corrections

           ratiompmc = 1.0714d0
       xR = xx*(1.d0 + ratiompmc) ! aR*sini/c [s]
       aR = xR/sidds             ! aR [s]

c  --> Retardation delay (Kopeikin & Schäfer 1999, eq. (130) expanded)

           epsRet1 =  an*xx/sidds*ratiompmc*ecc*su
           epsRet2 = -an*xx*sidds*ratiompmc/sqr1me2
     :            * (sw*cume + sqr1me2*cw*su)
     :            * (ecc*cw + (cw*cume - sqr1me2*sw*su)/onemecu)

       dShaRet = -2.d0*m2*(epsRet1 + epsRet2)/brace

c   --> Lensing
c       - accounts for Shapiro delay modification due to light deflection
c       - simplified version (Zschocke & Klioner)

       epsLen = 2.d0*m2/aR
       dShaLen = -2.d0*m2*epsLen/brace

c   --> Rotational bending delay (DK95) for parallel rotator (lam = i & eta = -pi/2)
c       version of DK95 is sufficient, and works for the whole orbit

           dShaBen1 = 2.d0*m2/(TWOPI*f0*xR)*cpsi/brace

c   --> Latitudinal delay due to light bending and profile changes
c       - for parallel rotator (lam = i & eta = -pi/2)
c       - simplified version of RL06, which works for the whole orbit
C       Still experimental (2016 Feb.)

       cosi = DSQRT(1.d0 - sidds**2) ! choice: > 0 => sign in cotchi0
       dShaBen2 = -2.d0*m2*cosi/(TWOPI*f0*xR)*spsi/brace * cotchi0

c   --> The sum of higher order contributions -- experimental
       dShapHo = dShaRet + dShaLen + dShaBen1 + dShaBen2
       dSha = dSha + shaphof * dShapHo

        endif

c  Aberration
        dAbe = a0*(spsi+ecc*sw)+b0*(cpsi+ecc*cw)

        delta = dRoe + dEin + dSha + dAbe

        diff  = DABS(delta - delta_old)
        if(diff.gt.epsNum) goto 10
c  --> Inversion of timing model by iteration: end of loop

        torb = -delta

c  Now we need the partial derivatives. Use DD equations 62a - 62k.
	csigma=xx*(-sw*su+sqr1me2*cw*cu)/onemecu
	ce=su*csigma-xx*sw-ecc*xx*cw*su/sqr1me2
	cx=sw*cume+sqr1me2*cw*su
	comega=xx*(cw*cume-sqr1me2*sw*su)
	cgamma=su
	cdth=-ecc*ecc*xx*cw*su/sqr1me2
	cm2=-2*dlogbr
C       csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace ! old...
        cshapmax=2*m2*(sw*cume+sqr1me2*cw*su)/brace * (1.0d0-sdds)
        cshaphof = dShapHo

	fctn(9)=cx*f0
	fctn(10)=ce*f0
	fctn(11)=-csigma*an*f0*86400.d0
	fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
	fctn(13)=comega*f0
	fctn(14)=ae*fctn(13)/(an*RAD*365.25d0*86400.d0)
	fctn(15)=cgamma*f0
	fctn(18)=0.5d-6*tt0*fctn(12)
        fctn(20)=cshapmax*f0
        fctn(39) =  cshaphof*f0
	fctn(22)=cm2*f0*SUNMASS
	fctn(23)=cdth*f0
	fctn(24)=tt0*fctn(9)
	fctn(25)=tt0*fctn(10)

	return
	end


