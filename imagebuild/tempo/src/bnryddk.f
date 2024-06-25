c      $Id$
	subroutine bnryddk(torb,fctn)

c  Damour et Deruelle modele pour le chronometrage des temps d'arrive 
c  au barycentre du systeme solaire, a la premiere approximation
c  post-newtonienne de la Relativite Generale.

C  IHS 20160126 Merge in Willen van Straten's bnryddk.f code from 2000. 
C  IHS 20160126 Maintain option to use fixed distance instead of parallax.
c  Annual-orbital parallax term (Kopeikin, 1995 ApJ 439, L5), and
c  Secular variation of x, w, i (Kopeikin, 1996 ApJ 467, L94)
c  added by Willem van Straten (August 2000)

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
c       Kopeikin's formulae are better in vector notation
        real*8 rp(3),drp_di(3),drp_dO(3),mu(3),re(3)
        real*8 I0(3),J0(3),K0(3),ip(3),jp(3),kp(3),kxrp(3),kxre(3)
	parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
        parameter (RAD=360.d0/twopi,AULTSC=499.00478364D0,ONE=1.0d0)
        parameter (OMAS=TWOPI/(360.0d0*3600.0d3),YEAR=365.25d0*86400.d0)
	include 'dp.h'
	include 'orbit.h'
        common /CRDS/ RCB(6),RBE(6),RCE(6),RSE(3),RCA(3),RSA(3),RBA(3),
     +             REA(3),PSID,EPSD,PC,PS,TSID,PRN(3,3),ATUT,UT1UT,ETAT,
     +             DTGR,TDIS,BCLT      
	include 'trnsfr.h'

        common/obsp/ site(3),pos(3),freqhz,bval,sitvel(3)

c  Calculate the secular change in i, omega, and x due to proper motion

        sin_Omega = sin(okom)
        cos_Omega = cos(okom)

c       Basis vectors rotated about the line of sight
        ip(1) = cos_Omega
        ip(2) = sin_Omega
        ip(3) = 0

        jp(1) = -sin_Omega
        jp(2) = cos_Omega
        jp(3) = 0

        kp(1) = 0
        kp(2) = 0
        kp(3) = 1

c       Find the sine and cosine of the position angles (alpha-RA delta-DEC)
c       'pos(3)' - unit vector in direction of pulsar
c       (part of common/obsp/ and calculated in timcalc)
        sin_delta = pos(3)
        cos_delta = cos(asin(sin_delta))
        sin_alpha = pos(2) / cos_delta
        cos_alpha = pos(1) / cos_delta

c       Basis vectors in the SSBC frame (Equations 10-11 in K95)
        I0(1) = -sin_alpha
        I0(2) = cos_alpha
        I0(3) = 0

        J0(1) = -cos_alpha*sin_delta
        J0(2) = -sin_alpha*sin_delta
        J0(3) = cos_delta

        K0(1) = cos_alpha*cos_delta
        K0(2) = sin_alpha*cos_delta
        K0(3) = sin_delta

c       convert proper motion from mas/year to radian/sec
        mu(1) = pmra * omas/year
        mu(2) = pmdec * omas/year
        mu(3) = 0

c       'ct' - coordinate time      (part of common/dp/ in dp.h)
c       't0' - epoch of periastron  (part of common/orbit/ in orbit.h)
        tt0=(ct-t0(1))*86400.d0

        xx = a1(1) + xdot*tt0
        oomz = omz(1)/RAD

        orbi = okin

        if (k96) then
c       Equation 10 in Kopeikin 1996
           orbi_dot = dot(mu, jp)
           orbi = orbi + orbi_dot * tt0
        endif

        si = sin (orbi)
        ci = cos (orbi)
        ti = si/ci

        if (k96) then
c       Equation 8 in Kopeikin 1996
           pm_xdot = xx * orbi_dot / ti
           xx = xx + pm_xdot * tt0
           
c       Equation 9 in Kopeikin 1996
           pm_omdot = dot(mu, ip) / si
           oomz = oomz + pm_omdot * tt0
        endif

c  Modify x and omega due to the annual-orbital parallax term
c  Willem van Straten (August 2000)

C  Calculate the SSBC earth position projected onto plane of sky:

c  'be(3)' - barycentric earth position vector
c            (part of common/trnsfr/ and produced in AU the bottom of timcalc)

c  Equations 15 and 16 in Kopeikin 1995
        d_i_not = - be(1) * sin_alpha + be(2) * cos_alpha
        d_j_not = - be(1) * sin_delta * cos_alpha
     +    - be(2) * sin_delta * sin_alpha + be(3) * cos_delta

        xpr = d_i_not * sin_Omega - d_j_not * cos_Omega
        ypr = d_i_not * cos_Omega + d_j_not * sin_Omega

        re(1) = dot(be,I0)
        re(2) = dot(be,J0)
        re(3) = dot(be,K0)
c  Equations 18 and 19 in Kopeikin 1995
c  'dpara2' - parallax in radians 
        if (usefixeddist) then  ! use fixed distance in kpc
          dpara2 = omas/(fixeddist)
        else                    ! use parallax distance
          dpara2 = dpara
        endif
        
        xx = xx + xx/ti * dpara2 * xpr
        si  = si + ci * dpara2 * xpr
        oomz = oomz - 1./si * dpara2 * ypr

c  Continue as before
	an=twopi/pb(1)
	k=omdot/(RAD*365.25d0*86400.d0*an)
	m2=am2*SUNMASS

	ecc=e(1)+edot*tt0
	er =ecc*(1.d0+dr)
	eth=ecc*(1.d0+dth)
		
	orbits=tt0/pb(1) - 0.5d0*(pbdot+xpbdot)*(tt0/pb(1))**2
	norbits=orbits
	if(orbits.lt.0.d0) norbits=norbits-1
	phase=twopi*(orbits-norbits)

c  Compute eccentric anomaly u by iterating Kepler's equation.
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
        omega=oomz + k*ae
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

c  Equations 6 and 7 in Kopeikin 1995
C IHS 20160127 These actually appear to be wrong by a factor of onemcu
C         Qu = cw * cae - sw * sae
C         Ru = sw * cae + cw * sae
        Qu = (cw * cae - sw * sae)*onemecu
        Ru = (sw * cae + cw * sae)*onemecu

c  Equation 5 in Kopeikin 1995, divided by a_p
        rp(1) = Qu*cos_Omega - Ru*sin_Omega*ci
        rp(2) = Qu*sin_Omega + Ru*cos_Omega*ci
        rp(3) = Ru*si

c  Partial derivative of rp with respect to i
        drp_di(1) = Ru*sin_Omega*si
        drp_di(2) = -Ru*cos_Omega*si
        drp_di(3) = Ru*ci

c  Partial derivative of rp with respect to Omega
        drp_dO(1) = -rp(2)
        drp_dO(2) = rp(1)
        drp_dO(3) = 0

        ap = xx/si

c  Equation 17 in Kopeikin 1995
C IHS 20160127 aopxvar is (the bracketed term in eq 17) * sini
        aopxvar = xpr*Ru*ci - ypr*Qu
C IHS 20160127 effectively xx/si*si so this is OK
        aopx = ap * dpara2 * aopxvar

        call cross (kp,rp,kxrp)
        call cross (kp,re,kxre)
C IHS 20160127 This aopx2 appears to be unused.
        aopx2 = ap * dpara2 * dot(kxrp,kxre)

C          write(*,210),aopx, aopx2
210     format ('aopx ',2e20.13)

c  partial derivatives:
        daopx_di = -xx * dpara2 * xpr*Ru
        daopx_dpx = ap * aopxvar / (rad*3600.d3*aultsc)
        daopx_dx  = aopx/xx
        daopx_dO = ap * dpara2 * ( ypr*Ru*ci + xpr*Qu )

c  Equation 1 in Kopeikin 1995 and statement at top of this file:
c  "Pulsar proper time is then TP=T+TORB."
c       torb = torb - aopx

c  px contribution from aopx
c  'dtdpx' - partial derivative with respect to parallax
c       (part of common/trnsfr/ and produced at the bottom of timcalc)
        dtdpx = dtdpx - daopx_dpx

c       contributions to proper motion on x and omega
c       dmurp_di=-xx*orbi_dot*onemecu*sin(omega+u)*tt0
C IHS 20160302 These should only be computed if using K96
        if(k96) then
          dmurp_di=ap*dot(mu, drp_di)*tt0
          dmurp_dO=ap*dot(mu, drp_dO)*tt0
          dmurp_dx=dot(mu, rp)*tt0/si
        else
          dmurp_di=0.;
          dmurp_dO=0.;
          dmurp_dx=0.;
        endif

C         write(*,211),ci*csi,dmurp_di,daopx_di
 211    format ('d/di ',3e20.13)

c  x contribution from aopx 
c IHS 20160302 and from the K96 parameters if non-zero; same for 52 and 53.
	fctn(9)=(cx+daopx_dx+dmurp_dx)*f0
	fctn(10)=ce*f0
	fctn(11)=-csigma*an*f0*86400.d0
	fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
	fctn(13)=comega*f0
	fctn(14)=ae*fctn(13)/(an*RAD*365.25d0*86400.d0)
	fctn(15)=cgamma*f0
	fctn(18)=0.5d-6*tt0*fctn(12)
c  partial derivatives with respect to i
        fctn(53)=(ci*csi+daopx_di+dmurp_di)*f0
	fctn(22)=cm2*f0*SUNMASS
	fctn(23)=cdth*f0
	fctn(24)=tt0*fctn(9)
	fctn(25)=tt0*fctn(10)
c  partial derivatives with respect to Omega
        fctn(52)=(daopx_dO+dmurp_dO)*f0

        si = 0

	return
	end


