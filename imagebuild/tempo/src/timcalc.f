c      $Id$
      subroutine TIMCALC(nmjdu,fmjdu,nmjdc,fmjdc,ctatv,etatc,wflag)
c
c     Routines taken (almost) exactly from the old Tempo BARTIM. 
c     For common block usage, see comments in ZTIM

      implicit real*8(A-H,O-Z)

	include 'dim.h'
	include 'acom.h'
        include 'dp.h'
        include 'trnsfr.h'

      common /CONST/ PI,TWOPI,SECDAY,CONVD,CONVS,AULTSC,VELC,EMRAT,OBLQ,
     +              GAUSS,RSCHW,AULTVL
      common /CRDS/ RCB(6),RBE(6),RCE(6),RSE(3),RCA(3),RSA(3),RBA(3),
     +             REA(3),PSID,EPSD,PC,PS,TSID,PRN(3,3),ATUT,UT1UT,ETAT,
     +             DTGR,TDIS,BCLT      
      common /OBSP/ SITE(3),POS(3),FREQHZ,BVAL,SITVEL(3)
	real*8 aa_star(3)
      real*8 VOBS(3), q_s(3)
      real*8 phisun
      real*8 solarn00
      integer wflag

c     get coordinates of site of observation in lt-sec
      if(ut1flag)call UT1RED(nmjdc,fmjdc,atut,ut1ut,quiet)
      fmjdu1 = fmjdu + ut1ut/secday
      nmjdu1 = nmjdu
      if(fmjdu1.ge.1.d0)then
         nmjdu1=nmjdu1+1
         fmjdu1=fmjdu1-1.d0
      endif
      call lmst(nmjdu1,fmjdu1,0.d0,tsid,sdd)  ! Greenwich mean sidereal time (turns)
      tsid=tsid*TWOPI
      CALL OBSITE(pc,prn,rea,tsid,site,sitvel)
      DO 100 I = 1, 3
        RBA(I) = RBE(I) + REA(I)
        RSA(I) = RSE(I) + REA(I)
        RCA(I) = RCE(I) + REA(I)
 100  continue

c     rotate vectors into ecliptic coordinates if necessary
      if (eclcoord) then
        call equ2ecl(rba)
        call equ2ecl(rsa)
        call equ2ecl(rca)
        call equ2ecl(rce)
        call equ2ecl(rse)
        call equ2ecl(rea)
        call equ2ecl(rce(4))
        call equ2ecl(sitvel)
      endif

         

c     compute accurate ET-AT, including diurnal & monthly terms
      ETATDM = DOT(RCE(4),REA)
      ETAT = CTATV + ETATDM + ETATC

c ###############
c ###  Irwin  ###
c ###############

      if (psrframe) then

c     start of TOD' iteration
         dt_path = 0.d0
         dt_path_old = 1.d0	!assure at least two loops
         
 101     if(abs(dt_path-dt_path_old).le.1.d-10) go to 129
         
         dt_path_old = dt_path
c     px of date for J2000 ephemeris only.
c     calculate position vector at current epoch in TDB.
         delt = (dble(nmjdu - 51545)+0.5d0+fmjdu+
     $        (etat+atut+dt_path)/secday)/36525d0
         do 10 k = 1,3
c     J2000 pulsar position (in au) at current TDB epoch.
 10      aa_star(k) = aa_q(k) + delt*aa_m(k)
c     
c     find time of flight of photon or pulse from source to solar system
c     barycentre within constant specified by epoch J2000.
         if(dpara.gt.0.d0) then
            dt_delay_num = 2.d0*delt*dot(aa_q, aa_m) +
     $           delt*delt*dot(aa_m, aa_m)
            dt_delay_denom = sqrt(dot(aa_q, aa_q) + dt_delay_num) +
     $           sqrt(dot(aa_q, aa_q))
            dt_delay = aultsc*dt_delay_num/dt_delay_denom/dpara
         else
c     if parallax is formally set to zero, it is likely to be small.
c     It is easy to show in this case that 
c     dt_delay = dt*(radial velocity + const*dt*parallax*
c     (tangential velocity)^2) + ...)/c
c     = delt*secday*36525.d0*rv/velc
c         dt_delay = delt*secday*36525.d0*rv/velc
c     RNM Since now entering pmrv = rv * px, assume rv = 0 if px = 0
            dt_delay=0.d0
         endif
         space_distance = sqrt(dot(aa_star, aa_star))
         do 20 k = 1,3
 20      pos(k) = aa_star(k)/space_distance
c     parallax of date in mas.
         px_date = 3600.d3*dpara/convd/space_distance
         
c     calculate light time from SSB to site
         dtdpxcon = convd/(3600.d3*aultsc)
         dbclt_num = -dtdpxcon*dot(rca, rca)
         bclt_num = 2.d0*dot(pos, rca) + px_date*dbclt_num
         bclt_denom = 1.d0 + sqrt(1.d0 - bclt_num*px_date*dtdpxcon)
         dbclt_denom = -0.5d0*dtdpxcon*(bclt_num + dbclt_num*px_date)/
     $        sqrt(1.d0 - bclt_num*px_date*dtdpxcon)
         bclt = bclt_num/bclt_denom
c     calculate negative derivative of bclt wrt px_date (should be good to
c     all orders, but definitely works to first order.)
         dtdpx = -(dbclt_num/bclt_denom - dbclt_denom*bclt/bclt_denom)
c     compute relativistic delay
         R = DSQRT(DOT(RSA,RSA))

c     use AA B36 formula with excellent approximation that sum of
c     distances, P+Q+E = 2Q.  We also ignore a small constant time
c     shift that depends on the distance of the star at epoch J2000.
c     in affect, this is absorbed into the TOD' definition,
c     see Irwin, McCarthy, and Marcy, 1995 (to be submitted).
         do 30 k = 1,3
c     heliocentric distance of star in light seconds.
 30      q_s(k) = aultsc*aa_star(k) + dpara*(rse(k) - rce(k))
         q_s_distance = sqrt(dot(q_s, q_s))
c     p = q - e
c     arg_ln = (q_s_distance - p_s_distance + r)/2/q_s_distance
         qmp_num = 2.d0*dot(q_s, rsa) - dpara*r*r
         qmp_denom = q_s_distance + sqrt(q_s_distance*q_s_distance -
     $        dpara*qmp_num)
         arg_ln = (qmp_num/qmp_denom + r)/q_s_distance
c     cannot apply correction if exactly centred on sun, but formula not
c     relevant closer than sun-grazing in any case.
         if(arg_ln.gt.0.d0) then
            dtgr = -2.d0*rschw*log(arg_ln)
         else
            dtgr = 0.d0
         endif
         
         DTDPPNG = 0.5D0*DTGR
         TDIS = 0D0

c        angle needed for planetary dispersion and phisun calculations
         CTH = DOT(POS,RSA)/DSQRT(DOT(RSA,RSA))
         THETH = DACOS(CTH)
c
         if (FREQHZ.LE.1D-1) then
           freqf = freqhz
           goto 115
         endif
c     correct the observing frequency to barycentric frame for dispersion
c     delay calculation
c     first add observatory velocity to EMB's (JMW)
         do 105 IAXIS = 1, 3
 105     VOBS(IAXIS) = RCE(IAXIS+3) + SITVEL(IAXIS)
         VOVERC = DOT(POS,VOBS)
C     R IS DISTANCE FROM SUN TO SITE
         FREQF=FREQHZ*(1D0-VOVERC)
C     compute interplanetary effect assuming 10 e-/cc at 1 AU
         if (SSDMFLAG.eq.0) then
           pldis = 0
         else
           solarn00 = solarn0 + solarn01*(nmjdu+fmjdu-pepoch)/365.25
           PLDIS = 2D14*THETH/R/DSQRT(1D0-CTH**2)*(solarn00/10)
           PLDIS = PLDIS/2.		
         endif
         TDIS = (BVAL+PLDIS)/FREQF**2

C     freq-dependent arrival time shifts
         do 106 i=1,NFDMAX
	   tdis = tdis + fdcof(i) * (log(freqf*1.0d-9)**i)
 106     continue
         
c     
 115     continue   
         phisun = 360*(PI-THETH)/TWOPI
         nmjdc=nmjdu
         dt_path = bclt-tdis-dtgr-dt_delay
      
         go to 101

 129     continue
c     after TOD' iteration subtract off first-order RV effect.
         dt_delay = dt_delay - delt*convd/3600.d0*pmrv/10.d0*aultsc
         fmjdc = fmjdu+(ETAT+ATUT+BCLT-TDIS-DTGR-dt_delay)/SECDAY
c         if(mod(n,100).eq.1)type *,delt,bclt,tdis,dtgr,dt_delay,fmjdc

c #############
c ###  Wex  ###
c #############

      else 

c     HERE:
c     aa_q: unit 3-vector in the direction of the pulsar (at Epoch)
c     aa_m: 3-vector in the direction of the transverse velocity of 
c           the pulsar in radians of arc per century (at Epoch)
c     pmrv: variable which is related to radial velocity * parallax 
c           in mas/yr

         pmrvrad=pmrv*CONVD/36000.d0 ! in rad/cen

c     start iteration
         dt_SSB = 0.d0
 201     dt_SSB_old = dt_SSB

         delt=(dfloat(nmjdu)-posep+fmjdu+
     :       (etat+atut+dt_SSB)/secday)/36525d0
         rr      = DOT(RCA,RCA)
         pmtrans = DSQRT(DOT(aa_m,aa_m))
         rcos1   = DOT(aa_q,RCA)
         pmtrans_rcos2 = DOT(aa_m,RCA)
         
         dt_roemer = rcos1
         dt_px     = -0.5*dpara*(rr-rcos1**2)/AULTSC
         dt_pm     = delt*pmtrans_rcos2
         dt_pmtt   = -0.5*pmtrans**2*delt**2*rcos1
         dt_pmtr   = -delt**2*pmrvrad*pmtrans_rcos2

         bclt = dt_roemer + dt_px + dt_pm + dt_pmtt + dt_pmtr

         dtdpx = 0.5*(rr-rcos1**2)*convd/(3600.d3*aultsc)
         dtdpmrv = delt**2*pmtrans_rcos2

c Compute Shapiro delay. To correct for the pulsar actual position
c just the transverse proper motion is taken into account

         DO 202 i=1,3
 202     pos(i) = aa_q(i)+delt*aa_m(i)
         pospos = pos(1)**2+pos(2)**2+pos(3)**2
         DO 203 i=1,3
 203     pos(i) = pos(i)/DSQRT(pospos)
         
         R=DSQRT(DOT(RSA,RSA))
         CTH = DOT(POS,RSA)/R
         arg_ln = (R/AULTSC)*(1.0D0+CTH)  ! R/AULTSC term added jun'00 
         dt_shapiro = -2.d0*RSCHW*DLOG(arg_ln) 


c correct the observing frequency to barycentric frame for dispersion
c delay calculation

c first add observatory velocity to EMB's (JMW)
         do 204 I=1, 3
 204     VOBS(I) = RCE(I+3) + SITVEL(I)
         VOVERC = DOT(POS,VOBS)
         if(FREQHZ.LE.1.D-1) then ! catch freq=0 (infinite frequency) case
           freqf = freqhz
         else
           FREQF=FREQHZ*(1D0-VOVERC)
         endif

c Compute interplanetary effect assuming 10 e-/cc at 1 AU
         THETH = DACOS(CTH)
         if (freqhz.le.1.d-1) then
           TDIS = 0D0
           txmx = 0.d0
         else
           if (SSDMFLAG.eq.0) then
             PLDIS = 0
           else
             solarn00 = solarn0 + solarn01*(nmjdu+fmjdu-pepoch)/365.25
             PLDIS = 2D14*THETH/R/DSQRT(1D0-CTH**2)*(solarn00/10)
             PLDIS = PLDIS/2.
           endif
           TDIS = (BVAL+PLDIS)/FREQF**2

C     freq-dependent arrival time shifts
           do 206 i=1,NFDMAX
	     tdis = tdis + fdcof(i) * (log(freqf*1.0d-9)**i)
 206       continue


C     calculate xmx delays
           txmx = 0.d0
           do i = 1, nxmx
             if (      xmxuse(i)
     +        .and.(xmxf1(i).lt.0.d0 .or. freqf.ge.xmxf1(i))
     +        .and.(xmxf2(i).lt.0.d0 .or. freqf.le.xmxf2(i))
     +        .and.(xmxr1(i).lt.0.d0 .or. nmjdu+fmjdu.ge.xmxr1(i))
     +        .and.(xmxr2(i).lt.0.d0 .or. nmjdu+fmjdu.le.xmxr2(i)))then
                txmx = txmx + xmx(i)*(freqf/1.0d6/xmxfrq0)**xmxexp(i)
c               print *,"inside",n,nmjdu,i,xmx(i),txmx
             endif
           enddo

         
         endif

         phisun = 360*(PI-THETH)/TWOPI
c           Should we worry about phisun in the case of barycenter TOAs?
c           No, because this routine is not called if nsite<0
        
         nmjdc=nmjdu

         dt_SSB = bclt-TDIS-dt_shapiro-txmx
         
         if(DABS(dt_SSB-dt_SSB_old).gt.1.0d-10) goto 201
         fmjdc = fmjdu+(ETAT+ATUT+bclt-TDIS-dt_shapiro-txmx)/SECDAY
      endif

      if(fmjdc.ge.1.d0)then
         nmjdc=nmjdc+1
         fmjdc=fmjdc-1.d0
      else if(fmjdc.lt.0.d0)then
         nmjdc=nmjdc-1
         fmjdc=fmjdc+1.d0
      endif

      if (phisun.lt.phimin) wflag = 0
      if (phisunout) write (38,fmt='(f10.6)') phisun
      if (dopplerout) write (41,fmt='(e16.9)') voverc

C     SET UP NECESSARY PARAMETERS FOR TEMPO: (JMW)

      FRQ=FREQF/1.0D6
      DO 160 ICOORD=1,3
  160 BE(ICOORD)=RCE(ICOORD)/AULTSC

      RETURN
      END

