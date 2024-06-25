c      $Id$
      SUBROUTINE BNRYELL1(torb,fctn)

c  Timing model for small-eccentricity binary pulsars, e<<1 (Wex 1998)
c
c  Instead of e and omega the Laplace parameters
c     epsilon1 = e*sin(omega)	
c     epsilon2 = e*cos(omega)
c  are used as new parameters. T0 is related to the ascending node (not to 
c  periastron as in BT, DD, ...)
c  
c  Time derivatives: 
c     nell1=0 -> fit for eps1dot,eps2dot 
c     nell1=1 -> fit for omdot,edot 
      
c  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch.
c  Pulsar proper time is then TP=T+TORB.
c  Units are such that c=G=1. Thus masses have units of seconds, with
c  one solar mass = 4.925490947 usec.

c  Also computes the binary orbit-related values of fctn: partial
c  derivatives of each arrival time residual with respect to the model
c  parameters.

c  Initial guesses for all parameters must be placed in common/orbit/ by the
c  calling program.  

      implicit real*8 (a-h,o-z)
      include 'dim.h'
      real*8 fctn(NPAP1),m2
      parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
      parameter (RAD=360.d0/twopi)
      include 'dp.h'
      include 'orbit.h'

      if (usefb) then
        an=twopi*fb(1)/FBFAC
      else
        an=twopi/pb(1)
      endif
      x0=a1(1)
      m2=am2*SUNMASS

c  In TEMPO we use the variable e(1) for eps1, omz(1) for eps2,
c  edot for eps1dot, and omdot for eps2dot 

      tt0=(ct-t0asc)*86400.d0
      if (usefb) then         ! orbital frequencies specified specified
        if (nfbj.eq.0) then   ! first use fb(1) and any jumps...
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
        fac = 1.              ! ...now work out higher order terms
        do j = 2, NFBMAX
          fac = fac/j
          orbits = orbits + fac*fb(j)*tt0**j / (FBFAC**j)
        enddo
      else                    ! orbital period specified
        orbits=tt0/pb(1)-0.5d0*(pbdot+xpbdot)*(tt0/pb(1))**2
      endif

      norbits=orbits
      if(orbits.lt.0.d0) norbits=norbits-1
      phase=twopi*(orbits-norbits)
	
      x=x0+xdot*tt0

      if(nell1.eq.0)then
         e1=eps1+eps1dot*tt0
         e2=eps2+eps2dot*tt0
      else
         ecc=DSQRT(eps1**2+eps2**2)
         ecc=ecc+edot*tt0
         w=DATAN2(eps1,eps2)
         w=w+omdot/(RAD*365.25d0*86400.d0)*tt0
         sw=DSIN(w)
         cw=DCOS(w)
         e1=ecc*sw
         e2=ecc*cw
      endif

      ! o(e^2) expression for Roemer delay from Norbert Wex and Weiwei Zhu
      ! This is equaiton (1) of Zhu et al (2019) but with a corrected typo:
      !    In the first line of that equation, ex->e1 and ey->e2
      !    In the other lines, ex->e2 and ey->e1
      ! See Email from Norbert and Weiwei to David on 2019-Aug-08
      ! The dre expression comes from Norbert and Weiwei; the derivatives
      ! were calculated by hand for tempo

      dre = x*(DSIN(phase)-0.5*(e1*DCOS(2.0*phase)-e2*DSIN(2.0*phase)))+
     +       (-x/8.)*(-2.*e1*e2*DCOS(phase) + 6.*e1*e2*DCOS(3.*phase) +
     +           3.*e1*e1*DSIN(phase) + 5.*e2*e2*DSIN(phase) +
     +           3.*e1*e1*DSIN(3.*phase) - 3.*e2*e2*DSIN(3.*phase))

      drep  = x*(DCOS(phase)+(e1*DSIN(2.0*phase)+e2*DCOS(2.0*phase))) +
     +       (-x/8.)*(2.*e1*e2*DSIN(phase) - 18.*e1*e2*DSIN(3.*phase) +
     +           3.*e1*e1*DCOS(phase) + 5.*e2*e2*DCOS(phase) +
     +           9.*e1*e1*DCOS(3.*phase) - 9.*e2*e2*DCOS(3.*phase))

      drepp  = x*(-DSIN(phase)+2.*(e1*DCOS(2.0*phase)-e2*DSIN(2.0*phase))) +
     +       (-x/8.)*(2.*e1*e2*DCOS(phase) - 54.*e1*e2*DCOS(3.*phase) +
     +           (-3.)*e1*e1*DSIN(phase) - 5.*e2*e2*DSIN(phase) +
     +           (-27.)*e1*e1*DSIN(3.*phase) + 27.*e2*e2*DSIN(3.*phase))
     
      ! old first-order equations 
      ! dre=x*(DSIN(phase)-0.5*(e1*DCOS(2*phase)-e2*DSIN(2*phase)))
      ! drep=x*DCOS(phase)
      ! drepp=-x*DSIN(phase)

c Shapiro delay part
      if (usefw10) then
        ds = -4.0/3.0*h3*DSIN(3.0*phase) + h4*DCOS(4.0*phase)
      else
        brace=1-si*DSIN(phase)
        dlogbr=dlog(brace)
        ds=-2*m2*dlogbr
      endif

      da=a0*sin(phase)+b0*cos(phase)

c  Now compute d2bar (cf. DD 52)

      d2bar=dre*(1-an*drep+(an*drep)**2+0.5d0*an**2*dre*drepp)+ds+da
      torb=-d2bar

c  Now we need the partial derivatives.

      Csigma   = x*DCOS(phase)
      Cx       = DSIN(phase)
      Ceps1    = -0.5*x*DCOS(2*phase)
      Ceps2    =  0.5*x*DSIN(2*phase)
      Cm2      = -2*dlogbr
      Csi      = 2*m2*DSIN(phase)/brace
c Note, these two assume nharm<=4:
      Ch3      = -4.0/3.0*DSIN(3.0*phase)
      Ch4      = DCOS(4.0*phase)

      fctn( 9) = Cx*f0
      fctn(10) = Ceps1*f0
      fctn(11) = -Csigma*an*f0*86400.d0
      if (usefb) then
        fctn(12) = fctn(11)*tt0/(86400.d0/(fb(1)/FBFAC)) ! needed to calculate other fctn()s
      else
        fctn(12) = fctn(11)*tt0/(pb(1)*86400.d0)
      endif
      fctn(13) = Ceps2*f0
      fctn(18) = 0.5d-6*tt0*fctn(12)
      if (usefw10) then
        fctn(20) = Ch4*f0
        fctn(22) = Ch3*f0
      else
        fctn(20) = Csi*f0
        fctn(22) = Cm2*f0*SUNMASS
      endif

      fctn(24) = tt0*fctn(9)

      if(nell1.eq.0)then
         Ceps1dot = tt0*Ceps1
         Ceps2dot = tt0*Ceps2
         fctn(39) = Ceps1dot*f0
         fctn(40) = Ceps2dot*f0
      else
         Comdot   = (Ceps1*e2-Ceps2*e1)*tt0
         Cedot    = (Ceps1*sw+Ceps2*cw)*tt0
         fctn(14) = Comdot*f0/(RAD*365.25d0*86400.d0)
         fctn(25) = Cedot*f0         
      endif

      if (usefb) then
        fctn(NPAR3+1) = -fctn(12)/(fb(1)/FBFAC)**2 / FBFAC
        do j = 2, NFBMAX
          fctn(NPAR3+j) = 1.d0/j * tt0 * fctn(NPAR3+j-1) / FBFAC
        enddo
        do j = 1, nfbj
          if (tfbj(j).lt.t0(i) .and. tfbj(j).lt.ct) then
            fctn(NPAR5+2*j-1) = 
     +           -fbj(j)/tt0*(fctn(NPAR3+1)*FBFAC)*86400
            fctn(NPAR5+2*j) =
     +           ((ct-t0(i))*86400/tt0)*(fctn(NPAR3+1)*FBFAC)*1.d6
          elseif (.not. (tfbj(j).gt.t0(i) .and. tfbj(j).gt.ct)) then
            fctn(NPAR5+2*j-1) =
     +           -fbj(j)/tt0*(fctn(NPAR3+1)*FBFAC)*86400
            fctn(NPAR5+2*j) =
     +           ((ct-tfbj(j))*86400/tt0)*(fctn(NPAR3+1)*FBFAC)
          endif
        enddo
      endif

      RETURN
      END
