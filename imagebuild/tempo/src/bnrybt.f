c      $Id$
	subroutine bnrybt(torb,fctn)

C  This version uses the method of Blandford and Teukolsky (ApJ 205,
C  580,1976)--including a solution for gamma.  For notation, see Taylor
C  et al, ApJ Lett 206, L53, 1976.

	implicit real*8 (a-h,o-z)
	parameter (twopi=6.28318530717958648d0)
	parameter (rad=360.d0/twopi)
	include 'dim.h'
	real*8 fctn(npap1),k
	include 'dp.h'
	include 'orbit.h'

C  Allow for more than one orbit:
C  Note: omdot, pbdot, xdot and edot are only applied to first orbit.
	torb=0.
        do 20 i=1,1+nplanets
        tt0=(ct-t0(i))*86400.d0
        if(i.eq.1)then
           orbits=tt0/pb(i) - 0.5d0*(pbdot+xpbdot)*(tt0/pb(i))**2
	   ecc=e(i)+edot*tt0
	   asini=a1(i)+xdot*tt0
	else
	   orbits=tt0/pb(i)
	   ecc=e(i)
	   asini=a1(i)
	endif

	norbits=orbits
	if(orbits.lt.0.d0) norbits=norbits-1
	phase=twopi*(orbits-norbits)
	if(i.eq.1)then
	   omega=(omz(i)+omdot*tt0/(86400.d0*365.25d0))/rad
	else
	   omega=omz(i)/rad
	endif

C  Use Pat Wallace's method of solving Kepler's equation
	ep=phase + ecc*sin(phase)*(1.d0+ecc*cos(phase))
	denom=1.d0-ecc*cos(ep)
 10	dep=(phase-(ep-ecc*sin(ep)))/denom
	ep=ep+dep
	if(abs(dep).gt.1.d-12) go to 10
	bige=ep
	if(nbin.eq.6) then
 	   u=ep
	   su=sin(u)
	   cu=cos(u)
	   onemecu=1.d0-ecc*cu
	   cae=(cu-ecc)/onemecu
	   sae=sqrt(1.d0-ecc**2)*su/onemecu
	   ae=atan2(sae,cae)
	   if(ae.lt.0.0) ae=ae+twopi
	   ae=twopi*orbits + ae - phase
	   an=twopi/pb(i)
	   k=omdot/(rad*365.25d0*86400.d0*an)
	   omega=omz(i)/rad + k*ae
	endif
	tt=1.d0-ecc*ecc
	som=sin(omega)
	com=cos(omega)
	alpha=asini*som
	beta=asini*com*sqrt(tt)
	sbe=sin(bige)
	cbe=cos(bige)
	q=alpha*(cbe-ecc) + (beta+gamma)*sbe
	r=-alpha*sbe + beta*cbe
	s=1.d0/(1.d0-ecc*cbe)
	torb=-q+(twopi/pb(i))*q*r*s + torb
	j=0
	if(i.eq.2) j=17
	if(i.eq.3) j=17+5
	fctn(9+j)=f0*(som*(cbe-ecc) + com*sqrt(tt)*sbe)
	fctn(10+j)=-f0*(alpha*(1.+sbe*sbe-ecc*cbe)*tt - 
     +     beta*(cbe-ecc)*sbe)*s/tt
	fctn(11+j)=-f0*(twopi/pb(i))*r*s*86400.d0
	fctn(12+j)=fctn(11+j)*tt0/(86400.d0*pb(i))
	fctn(13+j)=f0*asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe)
	if(i.eq.1) then
	   fctn(14)=fctn(13)*tt0/(rad*365.25d0*86400.d0)
	   fctn(15)=f0*sbe
	   fctn(18)=0.5d-6*fctn(12)*tt0
	   fctn(24)=tt0*fctn(9)
	   fctn(25)=tt0*fctn(10)
	endif
c  following commented out by DJN, 21 Nov '96
c	if(i.ge.2) fctn(18)=fctn(18) + fac*0.5d-6*fctn(12+j)*tt0
 20	continue

	return
	end
