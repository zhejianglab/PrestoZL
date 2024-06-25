c      $Id$
        subroutine newval(chisqr,nfree,rms0,rms1,nits,jits,wmax,nboot)

	implicit real*8 (a-h,o-z)
	parameter (TWOPI=6.28318530717958648d0)
	character*1 decsgn,binflag,label*6,dmlabel*2
	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'trnsfr.h'
	include 'dp.h'
	include 'orbit.h'
	include 'clocks.h'
	include 'eph.h'
	include 'glitch.h'

	data label/'Offset'/
	data dmlabel/'DM'/

	if(gain.ne.1.d0) then
	  do i = 1, NPAP1
	    freq(i)=freq(i)*gain
	  enddo
	endif

	if(jits.lt.nits)jits=jits+1

        if(posepoch.eq.0.)then
          write(31,1038) psrname(1:12),ephname(nephem)(1:5),clklbl(nclk),
     +         pepoch
 1038     format(/'PSR ',a12,'  Ephem.: ',a,'  Clock: ',a12,
     +         '  Ref. MJD: ',f12.4)
        else
          write(31,1039) psrname(1:10),ephfile(nephem)(1:5),
     +         clklbl(nclk)(1:10),pepoch, posepoch
 1039     format(/'PSR ',a,' Ephem: ',a,' Clk: ',a,
     +         ' P Ref:',f11.4,' Pos Ref:',f11.2)
        endif
        

	if (eclcoord) then
	  write (31,1040)
 1040     format (/5x,'LAMBDA',15x,'BETA',7x,'PM LAMBDA',2x,'PM BETA',
     +         4x,'PM RV',3x,'PARALLAX'/)
          c1 = 360.d0/TWOPI     ! radians to degrees
          c = 360.*60.*60./TWOPI !radians to arcsec
          cc =c*365.25*86400000. !radians/sec to mas/yr
          write (31,1041) c1*pra, c1*pdec, pmra, pmdec, pmrv, px
 1041     format (2f18.13,4f10.5)
          freq5 = freq(5)*c1
          freq6 = freq(6)*c1
          write (31,1041) freq6, freq5, 1.d-9*freq(8)*cc, 
     +         1.d-9*freq(7)*cc, 1.d1*freq(36)*c, freq(17)
          ferr(5) = ferr(5)*c1
          ferr(6) = ferr(6)*c1
          pmrerr=1.d-9*ferr(8)*cc
          pmderr=1.d-9*ferr(7)*cc
          write(31,1041) ferr(6),ferr(5),pmrerr,pmderr,
     +       1.d1*ferr(36)*c,ferr(17)
          pra=pra+freq(6)
          pdec=pdec+freq(5)
          pmra=pmra+1.d-9*freq(8)*cc
          pmdec=pmdec+1.d-9*freq(7)*cc
          pmrv=pmrv+1.d1*freq(36)*c
          px=px+freq(17)
          write (31,1041) c1*pra, c1*pdec, pmra, pmdec, pmrv, px
	else
	  write (31,1048)
 1048     format (/6x,'RA',16x,'DEC',12x,'PM RA',4x,'PM DEC',5x,'PM RV',
     +         3x,'PARALLAX'/)

          c=360.*60./TWOPI*60.
          cc=c*365.25*86400000.

          call radian(pra,irh,irm,rsec,123,1)
          call radian(pdec,idd,idm,dsec,1,1)
          decsgn=' '
          if(pdec.lt.0.) decsgn='-'
          idd=iabs(idd)
          
          write(31,1050) irh,irm,rsec,decsgn,idd,idm,dsec,pmra,pmdec,
     +         pmrv,px
 1050     format (i2.2,i3.2,f12.8,1x,a1,i2.2,i3.2,f11.7,4f10.4)

          freq5=freq(5)*c
          freq6=freq(6)*c/15.
          write(31,1051)freq6,freq5,1.d-9*freq(8)*cc,1.d-9*freq(7)*cc,
     +         1.d1*freq(36)*c,freq(17)
 1051     format (f17.8,f18.7,4f10.4)
          
          ferr(5)=ferr(5)*c
          ferr(6)=ferr(6)*c/15.
          pmrerr=1.d-9*ferr(8)*cc
          pmderr=1.d-9*ferr(7)*cc
          write(31,1051) ferr(6),ferr(5),pmrerr,pmderr,
     +         1.d1*ferr(36)*c,ferr(17)

          pra=pra+freq(6)
          pdec=pdec+freq(5)
          call radian(pra,irh,irm,rsec,123,1)
          call radian(pdec,idd,idm,dsec,1,1)
          decsgn=' '
          if(pdec.lt.0.) decsgn='-'
          idd=iabs(idd)
          pmra=pmra+1.d-9*freq(8)*cc
          pmdec=pmdec+1.d-9*freq(7)*cc
          pmrv=pmrv+1.d1*freq(36)*c ! Convert from rad/century to mas/yr
          px=px+freq(17)
          write(31,1050) irh,irm,rsec,decsgn,idd,idm,dsec,pmra,pmdec,
     +         pmrv,px
        endif

	p0z=1.d0/f0
	p1z=-1.d15*f1/f0**2
	f0z=f0
	f1z=f1
	f2z=f2
	f3z=f3
	df0=-freq(2)*1.d-9
	df1=-freq(3)*1.d-18
	df2=-freq(4)*1.d-27
	df3=-freq(NPAR11+1)*1.d-36
	f0=f0z+df0
	f1=f1z+df1
	f2=f2z+df2
	f3=f3z+df3
	kf1=0
	kf2=0
	kf3=0
	if(abs(f1).gt.0.)kf1=int(-log10(abs(f1)))+1
	if(abs(f2).gt.0.)kf2=int(-log10(abs(f2)))+1
	if(abs(f3).gt.0.)kf3=int(-log10(abs(f3)))+1
	sf1=10**(dfloat(kf1))
	sf2=10**(dfloat(kf2))
	sf3=10**(dfloat(kf3))

	write(31,1052)kf1,kf2,kf3
 1052	format(/10x,'F0',18x,'F1(D-',i2.2,')',7x,'F2(D-',i2.2,')',5x,
     +       'F3(D-',i2.2,')'/)

	write(31,1053)f0z,sf1*f1z,sf2*f2z,sf3*f3z
1053	format(f22.17,f18.12,f14.9,f12.6)

	write(31,1053)df0,sf1*df1,sf2*df2,sf3*df3

	ef0=ferr(2)*1.d-9
	ef1=ferr(3)*1.d-18
	ef2=ferr(4)*1.d-27
	ef3=ferr(NPAR11+1)*1.d-36
	write(31,1053)ef0,sf1*ef1,sf2*ef2,sf3*ef3

	write(31,1053)f0,sf1*f1,sf2*f2,sf3*f3

	write(31,1054)
 1054	format(/10x,'P0',18x,'P1(D-15)',9x,'DM',9x,'DM1',7x,'PPN GAM'/)

	ppng=1.d0
	write(31,1055)p0z,p1z,dm,dmcof(1),ppng
 1055	format(f22.19,f18.12,3f12.6)

	p0=1.d0/f0
	p1=-1.d15*f1/f0**2
	write(31,1055)p0-p0z,p1-p1z,freq(16),freq(NPAR9+1),freq(19)

	p0e=1.d-9*ferr(2)/f0**2
	p1e=1.d-3*dsqrt((2.d9*ferr(2)*f1/f0**3)**2+(ferr(3)/f0**2)**2)
	write(31,1055)p0e,p1e,ferr(16),ferr(NPAR9+1),ferr(19)

C IHS June 3 2011: Remove the condition as fitting DM0 is now allowed even here
C        if(.not.(ndmcalc.ge.2 .and. usedmx)) dm=dm+freq(16)
        dm=dm+freq(16)
	dmcof(1)=dmcof(1)+freq(NPAR9+1)
	ppng=ppng+freq(19)
	write(31,1055)p0,p1,dm,dmcof(1),ppng

	if((nfit(7).ne.0.or.nfit(8).ne.0) .and. (.not.eclcoord))
     +       call propmo(pmra,pmdec,pmrerr,pmderr,pra,pdec)

C  Compute braking index
	if(nfit(4).ne.0) then
	  brkind=f0*f2/(f1*f1)
	  brkerr=f0*1.d-27*ferr(4)/(f1*f1)
	  write(31,1080) brkind,brkerr
1080      format(/'Braking index:',f13.4,' +/-',f13.4)
        endif

	if(nfcalc.ge.4) then
	   do j=4,nfcalc,3
	    ia=j-3
	    ib=min(ia+2,nfcalc-3)
	    write(31,1081) j,j+1,j+2
 1081	    format(/25x,'f',i2.2,14x,'f',i2.2,15x,'f',i2.2)
	    write(31,1082) (f4(i),i=ia,ib)
	    write(31,1082) (-freq(NPAR11+1+i)*(1.d-9)**(i+4),i=ia,ib)
	    write(31,1082) (ferr(NPAR11+1+i)*(1.d-9)**(i+4),i=ia,ib)
	    write(31,1082) (f4(i)-freq(NPAR11+1+i)*(1.d-9)**(i+4),i=ia,ib)
 1082	    format(5x,1p,3d22.12,0p)
	  enddo

	   do i=1,nfcalc-3
	     f4(i)=f4(i)-freq(NPAR11+1+i)*(1.d-9)**(i+4)
	   enddo
	endif

	if(ndmcalc.ge.3) then
	  do k = 2, ndmcalc-1, 8
	    kk = min(k+7,ndmcalc-1)
	    write(31,1083) (dmlabel,i,i=k,kk)
 1083	    format(/8(3x,a2,i3.3,2x))
	    write(31,1084) (dmcof(i),i=k,kk)
	    write(31,1084) (freq(NPAR9+i),i=k,kk)
	    write(31,1084) (ferr(NPAR9+i),i=k,kk)
	    do j=k,kk
	      dmcof(j)=dmcof(j)+freq(NPAR9+j)
	    enddo
	    write(31,1084) (dmcof(i),i=k,kk)
 1084	    format(8f10.6)
	  enddo
	endif

	if(ngl.gt.0)then
	  do i=1,ngl
	    glphz=glph(i)
	    glf0z=glf0(i)
	    glf1z=glf1(i)
	    glf0dz=glf0d(i)
	    gltdz=gltd(i)
	    glph(i)=glphz-freq(NPAR1+(i-1)*NGLP+1)
	    glf0(i)=glf0z-freq(NPAR1+(i-1)*NGLP+2)*1.d-9
	    glf1(i)=glf1z-freq(NPAR1+(i-1)*NGLP+3)*1.d-18
	    glf0d(i)=glf0dz-freq(NPAR1+(i-1)*NGLP+4)*1.d-9
	    gltd(i)=gltdz-freq(NPAR1+(i-1)*NGLP+5)/86400.d0
	    glf0t=glf0(i)+glf0d(i)
	    glf0e=ferr(NPAR1+(i-1)*NGLP+2)*1.d-9
	    glf0de=ferr(NPAR1+(i-1)*NGLP+4)*1.d-9
	    write(31,1065)i,glepoch(i)
	    write(31,1066)
	    write(31,1067)glphz,glf0z,glf1z,glf0dz,gltdz
	    write(31,1067)glph(i)-glphz,glf0(i)-glf0z,glf1(i)-glf1z,
     +           glf0d(i)-glf0dz,gltd(i)-gltdz
	    write(31,1067)ferr(NPAR1+(i-1)*NGLP+1),glf0e,
     +           ferr(NPAR1+(i-1)*NGLP+3)*1.d-18,glf0de,
     +           ferr(NPAR1+(i-1)*NGLP+5)/86400.d0
	    write(31,1067)glph(i),glf0(i),glf1(i),glf0d(i),gltd(i)
	    if(glf0t.ne.0.d0)then
	       glf0te=sqrt(glf0e**2 + glf0de**2)
	       qq=glf0d(i)/glf0t
	       qqe=sqrt((glf0de/glf0t)**2+(glf0d(i)*glf0te/glf0t**2)**2)
	       write(31,1068)glf0t/f0,glf0te/f0,qq,qqe
	       iph=nint(glph(i))
	       fph=glph(i)-iph
	       glep1z=glepoch(i)+dglep(i,fph)
	       glep2z=glepoch(i)+dglep(i,fph-sign(1.d0,fph))
	       glepe=ferr(NPAR1+(i-1)*NGLP+1)/
     +   	    (abs(glf0(i)+glf0d(i))*86400.d0)
	       write(31,1069)glep1z,glep2z,glepe
	    endif
	  enddo
	endif
 1065	format(/' Glitch',i2,'  MJD:',f14.6)
 1066	format('    GLPH',7x,'GLF0 (Hz)',5x,'GLF1 (Hz/s)',4x,
     :       'GLF0D (Hz)',7x,'GLTD (d)')
 1067	format(f10.6,1p,3d15.6,0p,f15.6)
 1068	format(/' DeltaF/F:',1p,d13.6,'  Err:',d13.6,
     :       0p,'  Q:',f9.6,'  Err:',f9.6)
 1069	format(' MJD for zero phase:',f14.6,' or',f14.6,'  Err:',f10.6)

	if(usedmx) then
	  koff = NPAR6
	  do j = 1, (ndmx+3)/4
	    ib = j*4
	    ia = ib-3
	    ib = min(ib,ndmx)
	    write (31,1059) ("DM Off",i,i=ia,ib)
	    write (31,1060) (dmxr1(i),dmxr2(i),i=ia,ib)
	    write (31,1061) (dmx(i),i=ia,ib)
            write (31,1061) (freq(k),k=koff+2*ia-1,koff+2*ib-1,2)
            write (31,1061) (ferr(k),k=koff+2*ia-1,koff+2*ib-1,2)
	    do i = ia, ib
	      dmx(i) = dmx(i)+freq(koff+2*i-1)
	    end do
            write (31,1061) (dmx(i),i=ia,ib)
            write (31,1059) ("DM Dot",i,i=ia,ib)
            write (31,1060) (dmxr1(i),dmxr2(i),i=ia,ib)
            write (31,1061) (dmx1(i),i=ia,ib)
            write (31,1061) (freq(k),k=koff+2*ia,koff+2*ib,2)
            write (31,1061) (ferr(k),k=koff+2*ia,koff+2*ib,2)
            do i = ia, ib
              dmx1(i) = dmx1(i)+freq(koff+2*i)
            end do
            write (31,1061) (dmx1(i),i=ia,ib)
	  end do
	end if

        do i=1, NFDMAX
          fdcof(i) = fdcof(i) + freq(NPAR10+i)
        enddo 


	koff = NPAR12
	do j = 1, (nxmx+3)/4
	  ib = j*4
	  ia = ib-3
	  ib = min(ib,nxmx)
	  write (31,1059) ("XMX",i,i=ia,ib)
	  write (31,1060) (xmxr1(i),xmxr2(i),i=ia,ib)
	  write (31,1061) (xmx(i),i=ia,ib)
          write (31,1061) (freq(k),k=koff+2*ia-1,koff+2*ib-1,2)
          write (31,1061) (ferr(k),k=koff+2*ia-1,koff+2*ib-1,2)
	  do i = ia, ib
            if (xmxuse(i)) then
   	      xmx(i) = xmx(i)+freq(koff+2*i-1)
            endif
	  end do
          write (31,1061) (xmx(i),i=ia,ib)
          write (31,1059) ("XMXEXP",i,i=ia,ib)
          write (31,1060) (xmxr1(i),xmxr2(i),i=ia,ib)
          write (31,1061) (xmxexp(i),i=ia,ib)
          write (31,1061) (freq(k),k=koff+2*ia,koff+2*ib,2)
          write (31,1061) (ferr(k),k=koff+2*ia,koff+2*ib,2)
          do i = ia, ib
            if (xmxuse(i)) then
              xmxexp(i) = xmxexp(i)+freq(koff+2*i)
            endif
          end do
          write (31,1061) (xmxexp(i),i=ia,ib)
	end do
	      

C Get rms residuals and ntoa for output
        asig=asig*p0*1000.
        rms0=1000.d0*sigma1
        rms1=1000.d0*asig
        tres=rms1
        ntoa=n

C Output new parameters
	k=index(psrname,' ')-1
	open(71,file=psrname(1:k)//'.par',status='unknown')
	call outpar(nits,irh,irm,rsec,ferr(6),decsgn,idd,idm,
     +       dsec,ferr(5))

C Output binary parameters
	if(a1(1).ne.0.0) call newbin(nits,jits)

c Output error params
	call outerrpar

	if(nxoff.gt.0) then
	  koff=NPAR2
	  do 70 k=1,(nxoff+3)/4
	  ib=k*4
	  ia=ib-3
	  ib=min(ib,nxoff)
	  write(31,1059) (label,i,i=ia,ib)
1059	  format(/1x,5(7x,a6,i3))

	  write(31,1060) ((xjdoff(j,i),j=1,2),i=ia,ib)
1060	  format(/'MJD:',5(f8.1,'-',f7.1))
	  do 63 i = ia, ib
	  do 63 j = 1, 2
63	  if(xjdoff(j,i).lt.0.001) xjdoff(j,i)=0.

	  write(31,1061) (dct(i),i=ia,ib)
1061	  format(f17.8,4f16.8)
	  write(31,1061) (freq(i)*p0,i=koff+ia,koff+ib)
	  write(31,1061) (ferr(i)*p0,i=koff+ia,koff+ib)
	  do 65 i=ia,ib
65	  dct(i)=dct(i)+freq(koff+i)*p0
	  write(31,1061) (dct(i),i=ia,ib)
70	  continue
	  if (jumpout.or.nflagjumps.gt.0) call outjumppar
	endif
	    

C Close output .par file
	close(71)

	if(nboot.gt.0) write(31,1085) nboot
1085	format('Uncertainties by bootstrap Monte Carlo:',
     +    i6,' iterations.')
	write(31,1100) rms0,rms1
1100	format(/'Weighted RMS residual: pre-fit',f10.3,
     +  ' us. Predicted post-fit',f10.3,' us.')
	if (.not.quiet) write(*,1101)  rms0,rms1
1101	format(/' Weighted RMS residual: pre-fit',f10.3,
     +  ' us. Predicted post-fit',f10.3,' us.')
	if(chisqr.ne.0.d0) then
          if (chisqr.le.999999.) then
 	    write(31,1108) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
            if (.not.quiet) 
     +        write(*,1108) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
          else if (chisqr.le.9999999.) then
 	    write(31,1109) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
            if (.not.quiet)
     +        write(*,1109) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
          else
 	    write(31,1110) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
            if (.not.quiet) 
     +        write(*,1110) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
          endif
        endif
1108	  format(' Chisqr/nfree: ',f9.2,'/',i5,' = ',f15.9,
     +    '   pre/post:',f7.2,'   Wmax:',f7.1)
1109	  format(' Chisqr/nfree: ',f10.2,'/',i5,' = ',f15.9,
     +    '   pre/post:',f8.2,'   Wmax:',f7.1)
1110	  format(' Chisqr/nfree: ',f11.2,'/',i5,' = ',f15.9,
     +    '   pre/post:',f9.2,'   Wmax:',f7.1)

	if(gro.and.(nits.eq.0.or.jits.eq.nits)) then
	  open(33,file='gro.1',status='unknown')
	  if(dt2sec.gt.0.d0) dt2sec=dt2sec-p0
	  t0geo=t0geo-dt2sec/86400.d0
	  if(abs(t0geo-pepoch).gt.1.d-3.and..not.quiet) 
     +      write(*,1112) t0geo,pepoch
1112	  format(/5x,'### Warning: t0geo=',f10.3,' and pepoch=',f10.3,
     +      ' do not match! ###')
	  dt=(t0geo-pepoch)*86400.d0
	  ff0=f0 + f1*dt + 0.5d0*f2*dt**2
	  ff1=f1 + f2*dt
	  pp0=1.d0/ff0
	  pp1=-ff1*ff0**(-2)
	  pp2=-f2*ff0**(-2) + 2.d0 * ff1**2 * ff0**(-3)
	  open(99,file='gro.99',status='unknown')
	  if(oldpar)then
	     write(99,1113) pp0,1.d15*pp1,t0geo,1.d30*pp2
 1113	     format('P',f18.16,1x,f12.5,8x,f11.5,9x,f12.1)
	  else
	     write(99,'(''F0'',f22.16)')ff0
	     write(99,'(''F1'',1p,d22.12)')ff1
	     write(99,'(''F2'',1p,d22.12)')f2
	     write(99,'(''PEPOCH'',f18.6)')t0geo
	  endif
	  close(99)
	  f0=ff0
	  f1=ff1
	  rmsmp=asig/p0
	  mjd1=amjd1
	  mjd2=amjd2+1.d0
	  binflag=' '
	  if(a1(1).ne.0.0) binflag='*'
	  write(33,1120) psrname(1:8),irh,irm,rsec,decsgn,idd,idm,dsec,
     +       mjd1,mjd2,t0geo,f0,f1,f2,rmsmp,obsflag,binflag,
     +       ephfile(nephem)(1:5),psrname
1120	  format(a8,2i3.2,f7.3,1x,a1,i2.2,i3.2,f6.2,2i6,f16.9,
     +      f18.13,1p,d13.5,d11.2,0p,f5.1,2(1x,a1),1x,a,1x,a)
	  close(33)
	endif

	return
	end

C=======================================================================

        real*8 function dglep(igl,fph)

	include 'dim.h'
	include 'glitch.h'

	real*8 fph, plim, t1, dph

	tds=gltd(igl)*86400.d0
	niter=0
	plim=1.d-6
	dph=1000.
	t1=-fph/(glf0(igl) + glf0d(igl))
	do while(abs(dph) .gt. plim)
	   dph=fph + glf0(igl)*t1 + 0.5d0*glf1(igl)*t1*t1
	   if(tds.gt.0.d0)dph=dph+glf0d(igl)*tds*(1.d0-exp(-t1/tds))
	   t1=t1-dph/(glf0(igl) + glf0d(igl))
	   niter=niter+1
	   if(niter.gt.100)then
	      write(*,*)'*** Glitch epoch convergence failed ***'
	      dglep=0.d0
	      return
	   endif
	enddo
	dglep=t1/86400.d0
	
	return
	end
