c      $Id$
	subroutine bnrybtx(torb,fctn)

C  This version uses the method of Blandford and Teukolsky (ApJ 205,
C  580,1976)--including a solution for gamma.  For notation, see Taylor
C  et al, ApJ Lett 206, L53, 1976.

c  Modified for use with many derivatives of pb and x=asini 
c  DJN, 9 July 99, based on earlier codes by djn and za
c  DOES NOT USE PB(1) and PBDOT.  USES FB(1..NFB) instead!

	implicit real*8 (a-h,o-z)
	parameter (twopi=6.28318530717958648d0)
	parameter (rad=360.d0/twopi)
	include 'dim.h'
	real*8 fctn(npap1)
	include 'dp.h'
	include 'orbit.h'

C  Allow for more than one orbit:
C  Note: omdot, pbdot, xdot and edot are only applied to first orbit.
C IHS 20120606 Add multiple derivatives for e and om
	torb=0.
        do 20 i=1+nplanets, 1, -1
        tt0=(ct-t0(i))*86400.d0 + torb
        if(i.eq.1)then
C                         fb(1) terms of orbital frequency
           if (nfbj.eq.0) then
             orbits = (fb(1)/FBFAC)*tt0 
           else
	     orbits = (fb(1)/FBFAC)*tt0 
	     do j = 1, nfbj
	       if (tfbj(j).lt.t0(i) .and. tfbj(j).lt.ct) then
		 orbits = orbits + (ct-t0(i))*86400*fbj(j)
	       elseif (.not. (tfbj(j).gt.t0(i) .and. tfbj(j).gt.ct)) then
		 orbits = orbits + (ct-tfbj(j))*86400*fbj(j)
	       endif
	     enddo
	   endif
c                         higher order terms
           fac = 1.
	   do j = 2, NFBMAX
	     fac = fac/j
             orbits = orbits + fac*fb(j)*tt0**j / (FBFAC**j)
           enddo
	   ecc=e(i)+edot*tt0
	   fac = 1.
	   do j = 2, NEDOTMAX
	     fac = fac/j
	     ecc = ecc + fac*edot2(j)*tt0**j
           enddo
	   asini=a1(i)+xdot*tt0
	   fac = 1.
	   do j = 2, NXDOTMAX
	     fac = fac/j
	     asini = asini + fac*xdot2(j)*tt0**j
           enddo
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
	   fac = 1.
	   do j = 2, NOMDOTMAX
	     fac = fac/j
	     omega = omega + fac*omdot2(j)*tt0**j
           enddo
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
	if (i.eq.1) then
	  torb=-q+(twopi*fb(1)/FBFAC)*q*r*s + torb
        else
	  torb=-q+(twopi/pb(i))*q*r*s + torb
        endif
	ii=0
	if(i.eq.2) ii=17
	if(i.eq.3) ii=17+5
	fctn(9+ii)=f0*(som*(cbe-ecc) + com*sqrt(tt)*sbe) ! a sini
	fctn(10+ii)=-f0*(alpha*(1.+sbe*sbe-ecc*cbe)*tt - 
     +       beta*(cbe-ecc)*sbe)*s/tt ! e
	if (i.eq.1) then
				!note: we won't actually use the fctn(12+0)
				! value except to calculate other fctn(..)'s.
	  fctn(11+ii)=-f0*(twopi*fb(1)/FBFAC)*r*s*86400.d0 ! T0 (in days)
	  fctn(12+ii)=fctn(11+ii)*tt0/(86400.d0/(fb(1)/FBFAC)) ! PB
	else
	  fctn(11+ii)=-f0*(twopi/pb(i))*r*s*86400.d0 ! T0 (in days)
	  fctn(12+ii)=fctn(11+ii)*tt0/(86400.d0*pb(i)) ! PB
	endif
	fctn(13+ii)=f0*asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe) ! omega
	if(i.eq.1) then
	   fctn(14)=fctn(13)*tt0/(rad*365.25d0*86400.d0)  ! omdot
	   fctn(15)=f0*sbe                ! gamma
	   fctn(18)=0.5d-6*fctn(12)*tt0   ! pb-dot.  1e-6 cancels term in newbin
	   fctn(24)=tt0*fctn(9)           ! xdot
	   fctn(25)=tt0*fctn(10)          ! edot

	   fctn(NPAR3+1) = -fctn(12)/(fb(1)/FBFAC)**2 /FBFAC
           do j = 2, NFBMAX
             fctn(NPAR3+j) = 1.d0/j * tt0 * fctn(NPAR3+j-1) /FBFAC
           enddo
           fctn(NPAR4+1) = 0.5d0 * tt0 * fctn(24)  ! More x derivatives
           do j = 3, NXDOTMAX
             fctn(NPAR4+j-1) = (1.d0/j) * tt0 * fctn(NPAR4+j-2) 
           enddo
           fctn(NPAR7+1) = 0.5d0 * tt0 * fctn(25)  ! More e derivatives
           do j = 3, NEDOTMAX
             fctn(NPAR7+j-1) = (1.d0/j) * tt0 * fctn(NPAR7+j-2) 
           enddo
           fctn(NPAR8+1) = 0.5d0 * tt0**2 * fctn(13)  ! More om derivatives
           do j = 3, NOMDOTMAX
             fctn(NPAR8+j-1) = (1.d0/j) * tt0 * fctn(NPAR8+j-2) 
           enddo
	   do j = 1, nfbj
	     if (tfbj(j).lt.t0(i) .and. tfbj(j).lt.ct) then
               fctn(NPAR5+2*j-1) = 
     +              -fbj(j)/tt0*(fctn(NPAR3+1)*FBFAC)*86400
               fctn(NPAR5+2*j) = 
     +              ((ct-t0(i))*86400/tt0)*(fctn(NPAR3+1)*FBFAC)*1.d6
             elseif (.not. (tfbj(j).gt.t0(i) .and. tfbj(j).gt.ct)) then
	       fctn(NPAR5+2*j-1) = 
     +              -fbj(j)/tt0*(fctn(NPAR3+1)*FBFAC)*86400
               fctn(NPAR5+2*j) = 
     +              ((ct-tfbj(j))*86400/tt0)*(fctn(NPAR3+1)*FBFAC)
	     endif
	   enddo
	endif
c  following commented out by DJN, 21 Nov '96
c	if(i.ge.2) fctn(18)=fctn(18) + fac*0.5d-6*fctn(12+ii)*tt0
 20	continue

	return
	end
