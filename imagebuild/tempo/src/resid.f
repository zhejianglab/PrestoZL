c      $Id$
	subroutine resid(nct,fct,dnprd,dphase,dnpls,nits,jits)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*4 gasdev
	real*8 tgl(NGLT)
	integer idum

	include 'acom.h'
	include 'bcom.h'
	include 'dp.h'
	include 'orbit.h'
	include 'glitch.h'

        save

        data idum/-1/

	if(n.eq.1)then
	  ct1=ct
	  p0firs=p0
	  nf0=f0
	  ff0=f0-nf0
	  nepoch=pepoch
	  fepoch=pepoch-nepoch
	  if(ngl.gt.0)then
	    do i=1,ngl
	      nglep=glepoch(i)
	      fglep=glepoch(i)-nglep
	      ngld=nglep-nepoch
	      fgld=fglep-fepoch
	      tgl(i)=(ngld+fgld)*86400.d0
	    enddo
	  endif
	endif

	torb=0.
	if(nbin.eq.1.or.nbin.eq.6.or.nplanets.gt.0) call bnrybt(torb,x)
	if(nbin.eq.2.or.nbin.eq.5) call bnryeh(torb,n,x)
	if(nbin.eq.3) call bnrydd(torb,x)
	if(nbin.eq.4) call bnryddgr(ct,f0,n,torb,x)
	if(nbin.eq.7) call bnryddt(ct,f0,torb,x)
	if(nbin.eq.8) call bnrymss(torb,x)
	if(nbin.eq.9) call bnryell1(torb,x)
	if(nbin.eq.10) call bnrybtx(torb,x)
	if(nbin.eq.13) call bnrydds(torb,x)
        if(nbin.eq.15) call bnryddfwhiecc(torb,x)
	if(nbin.eq.14) call bnryddk(torb,x)

	ntpd=nct-nepoch
	ftpd=fct-fepoch+torb/86400.d0
	tp=(ntpd+ftpd)*86400.d0
	x(1)=1.d0
	t9=tp/1.d9
	x(2)=t9
	x(3)=t9*x(2)/2.d0
	x(4)=t9*x(3)/3.d0
	x(7)=t9*x(5)
	x(8)=t9*x(6)/cos(pdec)
	x(NPAR11+1)=t9*x(4)/4.d0
        do i = 4, nfcalc
          x(NPAR11+i-2) = t9*x(NPAR11+i-3)/(i+1.d0)
        enddo

	phase4=0.d0
	if(ngl.gt.0)then
	  do i=1,ngl
	    if(tp.ge.tgl(i))then
	      dt1=tp-tgl(i)
	      dt9=dt1/1.d9
	      td9=gltd(i)*86400.d-9
	      x(NPAR1+(i-1)*NGLP+1)=1.d0
	      x(NPAR1+(i-1)*NGLP+2)=dt9
	      x(NPAR1+(i-1)*NGLP+3)=0.5d0*dt9**2
	      if(td9.ne.0.d0)then
	        expf=exp(-dt9/td9)
	        x(NPAR1+(i-1)*NGLP+4)=td9*(1.d0-expf)
	        x(NPAR1+(i-1)*NGLP+5)=glf0d(i)*(1.d0-(1.d0+dt9/td9)*expf)
	      else
	        expf=1.d0
	      endif
	      phase4=phase4+glph(i)+glf0(i)*dt1+0.5d0*glf1(i)*dt1**2
     +          +glf0d(i)*gltd(i)*86400.d0*(1.d0-expf)
	    else
	      do j=1,NGLP
	        x(NPAR1+(i-1)*NGLP+j)=0.d0
	      enddo
	    endif
	  enddo
	endif

	phaseint=nf0*ntpd*86400.d0
	phase2=(nf0*ftpd+ntpd*ff0+ftpd*ff0)*86400.d0
	phase3=0.5d0*f1*tp**2 + (f2/6.d0)*tp**3 + (f3/24.d0)*tp**4
	if(nfcalc.ge.4) then
	  dfac=24.d0
	  do 10 i=1,nfcalc-3
	  dfac=dfac*(i+4)
10	  phase3=phase3 + (f4(i)/dfac)*tp**(i+4)
	endif

	phase5=phase2+phase3+phase4
	if(n.eq.1) phas1=mod(phase5,1.d0)
	phase5=phase5-phas1
	nphase=nint(phase5)
	dnprd=phaseint+nphase
	if (npulsein) then
	  read(35,*) dnprdin
          ddnprd = dnprd-dnprdin
          dnprd = dnprdin
        else
          ddnprd = 0
        endif
c  used to write out pulse number here;
c  moved this to arrtim, so that it would
c  come after tracking corrections
c	if (npulseout.and.jits.eq.0) then
c	  write(35,fmt='(f14.0)') dnprd
c         endif
	phasefrac=phase5-nphase
	dt=phasefrac+dphase+ddnprd
	if(nits.gt.0) then
	  iphase=dphase+dsign(0.5d0,dphase)
	  if(jits.eq.0) then
	    dnpls=dnprd-iphase
	  else
	    kpls=dnprd-iphase-dnpls
	    dt=dt+kpls
	  endif
	endif


	if(sim) dt=1.d-6*dither*gasdev(idum)*f0
	return
	end
