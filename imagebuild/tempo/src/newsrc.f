c      $Id$
        subroutine newsrc(nits,jits,nboot)

	implicit real*8 (A-H,O-Z)
	character DECSGN*1,path*160,str1*80,str2*80
	character*1 asite, dtmp, ntmp
        character*80 bsite
        parameter (TWOPI=6.28318530717958648d0)

	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'trnsfr.h'
	include 'clocks.h'
	include 'dp.h'
	include 'orbit.h'
	include 'eph.h'
	include 'glitch.h'

        integer sitea2n  ! external function

             ! jits = iteration number
             ! jits==0: first time through, initialize
             ! jits==-1: special case, don't initialize, then set jits=0
	if(jits.ne.0) goto 30

C  Zero out all input parameters and fit coefficients.
	call zeropar(nits)      

	if(oldpar)then
        
C --- old input format ---        

C  Columns:                 1 ... 25       27   29    33   34   40   42
	   read(50,1010) (nfit(i),i=1,25),nbin,nprnt,ntmp,iboot,ngl,nddm,
C                44    46     47      51 ... 60
     +   	nclk,nephem,ncoord,(nfit(i),i=26,35)
 1010	   format(25z1,1x,z1,1x,i1,3x,a1,z1,4x,3i2,1x,2i1,3x,10i1)

           if (ntmp.eq." ") ntmp = '1'  ! if blank, use default nits=1
           read (ntmp,*) nits

           if(nbin.ge.9)then
	      nplanets=nbin-8
	      nbin=1
	   endif

C  Set nfit(k) = 1 if the k'th parameter is to be fit for
C  If nfit(1) > 1 then no fit is done, but tz2.tmp is created
C  Set nprnt = 1, 2, ... to list every 10, 100, ... TOAs in output file
C  Set nits = max no of iterations (9 ==> iterate until convergence)

C nbin  Binary pulsar model		nclk	Reference time scale
C--------------------------------------------------------------------
C   1 	Blandford-Teukolsky		  0	Uncorrected
C   2	Epstein-Haugan			  1	UTC(NIST)
C   3	Damour-Deruelle			  2	UTC
C   4	DDGR				  3	TT(BIPM)
C   5	H88				  4	(PTB)
C   6	BT+				  5	AT1
C   7	DDT
C   8   MSS (wex main-sequence companion model)
C   9   ELL1 (wex small eccentricity model)
C  10   BTX (BT plus freq derivatives and jumps)
C  11	BT, 2 orbits
C  12	BT, 3 orbits
C  13   DDS
C  14   DDK
c  15   DDFWHE 

C  special case: nfit(3)>1 implies that frequency derivatives 1 through nfit(3) 
C  are to be fit.  (Example: if nfit(3)=5, fit F1 through F5).  Set the relevant flags.

           nfcalc = nfit(3)
           if (nfit(4).eq.1) nfcalc = max(nfit(3),4)

           if (nfit(3).ge.2) nfit(4)=1
           if (nfit(3).ge.3) then
             do i = 3, nfit(3)
               nfit(NPAR11+i-2) = 1    
             enddo
           endif
           if (nfit(3).gt.1) nfit(3)=1  

C  Get initial parameter values from lines 2-4 of data file
	   read(50,1020) psrname,pra,pdec,pmra,pmdec,rv,p0,p1,pepoch,
     +          p2,px,dm
 1020	   format(a12,8x,f20.1,f20.1,3f10.0/1x,d19.9,2d20.6,2d10.0/
     +          8x,f12.0)

	   if(nbin.ne.0) read(50,1030) a1(1),e(1),t0(1),pb(1)
 1030	   format(4f20.0)
	   if(nbin.ne.0.and.nbin.ne.7) read(50,1031) omz(1),omdot,gamma,
     +          pbdot,si,am,am2,dr,dth,a0,b0
 1031	   format(6f10.0,f6.0,3f4.3,f2.1)
	   if(nbin.eq.7) read(50,1032) omz(1),omdot,gamma,pbdot,
     +         si,am,am2,bp,bpp
 1032	   format(7f10.0,2f5.0)

	   if(nplanets.gt.0)then   ! Read params for 2nd and 3rd orbits
	      do i=2,nplanets+1
		 read(50,1030) a1(i),e(i),t0(i),pb(i)
		 read(50,1031) omz(i)
	      enddo
	   endif
           
C  Convert PB to days if greater than 3600 (ie min Pb for old style = 1h)
	   do i=1,3
	      if(pb(i).gt.3600.d0)pb(i)=pb(i)/86400.d0
	   enddo

	   if(ngl.ne.0)then
	     do i=1,ngl
	       read(50,1035)(nfit(NPAR1+(i-1)*NGLP+j),j=1,NGLP),
     +	         glepoch(i),glph(i),glf0(i),glf1(i),glf0d(i),gltd(i)
	       if((nfit(NPAR1+(i-1)*NGLP+4).ne.0.or.
     +           nfit(NPAR1+(i-1)*NGLP+5).ne.0).and.gltd(i).eq.0.d0)then
                 write(0,*)' WARNING: Exp term requested but gltd = 0;', 
     +             ' set = 1.0d'
		   gltd(i)=1.d0
		 endif
	      enddo
	   endif
 1035	   format(5i1,5x,6f12.0)

C  Calculate header lines to be skipped when reading TOAs
          nskip=4
          if(nbin.gt.0)     nskip=6
          if(nplanets.gt.0) nskip=6+2*nplanets
          if(ngl.gt.0)      nskip=nskip+ngl

	else
        
c --- free format header ---

	   call rdpar(nits)
	   if (parunit.eq.49) then      ! separate par file
               nskip = 0
               if (.not. tz) close(49)  ! leave open in tz mode for multi-day 
                                        !   calculations which re-read the par file
           endif

        endif

	if(gro)nits=1
	if(xitoa)nclk=2 ! Force correction to UTC for ITOA output

C  ELL1: convert e,omega -> eps1,eps2 (if necessary) and set nell1
C  nell1=0 -> fit for eps1dot,eps2dot 
C  nell1=1 -> fit for omdot,edot 
 
	if(nbin.eq.1.and.t0asc.ne.0) call ell12bt()

	if(nbin.eq.9)then
	   if(t0(1).ne.0.) call bt2ell1()
	   if(nfit(14).ne.0 .or. nfit(25).ne.0
     +        .or. omdot.ne.0. .or. edot .ne.0.) nell1=1
	endif

C  Convert units

	p1      = p1      *1.d-15
	p2      = p2      *1.d-30
	dr      = dr      *1.d-6
	dth     = dth     *1.d-6
	a0      = a0      *1.d-6
	b0      = b0      *1.d-6
	pbdot   = pbdot   *1.d-12
	xpbdot  = xpbdot  *1.d-12
	xdot    = xdot    *1.d-12
	edot    = edot    *1.d-12
	eps1dot = eps1dot *1.d-12
	eps2dot = eps2dot *1.d-12

	if(pepoch.gt.2400000.5d0) pepoch=pepoch-2400000.5d0
	if(posepoch.gt.0.)then
	   posep=posepoch
	else
	   posep=pepoch
	endif
	if(dmepoch.gt.0.)then
	   dmep=dmepoch
	else
	   dmep=pepoch
	endif
	ndmcalc=max(nfit(16),ndmcalc)

	if(eclcoord)then
	   pra=pra*TWOPI/360.d0
           pdec=pdec*TWOPI/360.d0
	else
	   pra=ang(3,pra)	! hhmmss.ss to radians
           pdec=ang(1,pdec)	! ddmmss.ss to radians
	endif

	if(ncoord.eq.0)then
	   pmra=pmra*1d-1
	   pmdec=pmdec*1d-1
	   px=px*1d-3
	   call fk4tofk5(pra,pdec,pmra,pmdec,px,rv)
	   pmra=pmra*1d1
	   pmdec=pmdec*1d1
	   px=px*1d3
	   if (.not.quiet) write(31,1036)
 1036	   format(' Input B1950 coords converted to J2000')
	endif

        if(f0.eq.0.)f0=1.d0/p0
        if(f1.eq.0. .and. p1.ne.0.)f1=-p1*f0**2
        if(f2.eq.0. .and. p2.ne.0.)f2=2.d0*f1*f1/f0 - p2*f0*f0

	nboot=0
	if(iboot.gt.0) nboot=2**iboot

c  Open ephemeris file
	if(nephem.lt.1.or.nephem.gt.kephem)stop 'Invalid ephemeris nr'
	kd=index(ephdir,' ')-1
	kf=index(ephfile(nephem),' ')-1
	path=ephdir(1:kd)//ephfile(nephem)(1:kf)
	call ephinit(44,path,ephbigendian(nephem))

C  Check to make sure selected parameters are consistent with model
	if(nbin.eq.0) then
 	   do 1 I = 9, 15
1	   nfit(I) =0
	   nfit(18)=0
	   do 11 i=20,35
11	   nfit(i) =0
	endif

	if(nbin.eq.1)then
	   do 2 I = 20, 23
2	   nfit(I)=0
	endif

	if(nbin.eq.2) then
	   do 21 i=21,23
21	   nfit(i)=0
	endif

	if(nbin.eq.4)then
	   nfit(15)=0
	   nfit(20)=0
	   nfit(23)=0
	endif
        
	if(nbin.ne.9)then	!Convert old style to MJD, t0(1)=0 in ELL1
	   if(t0(1).lt.39999.5d0) t0(1)=t0(1)+39999.5d0
	   if(t0(2).lt.39999.5d0) t0(2)=t0(2)+39999.5d0
	   if(t0(3).lt.39999.5d0) t0(3)=t0(3)+39999.5d0
	endif

	if(oldpar)then	       ! in old style files, xomdot and xpbdot
 	   xomdot=0.	       ! replaced omdot and pbdot in value and
	   xpbdot=0.	       ! flag fields in the header.
	   if(nbin.eq.4)then
	      xpbdot=pbdot
	      xomdot = omdot
	      omdot=0.
	      pbdot=0.
	      if(nfit(14).eq.1)then
	         nfit(14)=0
	         nfit(37)=1
	      endif
	      if(nfit(18).eq.1)then
	         nfit(18)=0
	         nfit(38)=1
	      endif
	   endif
	endif

C Read clock corrections
	if(nclk.gt.0)then
	   kc=index(clkdir,' ')-1                     ! Obs to NIST
	   path=clkdir(1:kc)//clkfile(1)
	   open(20,file=path,status='old',err=900)
	   ifile = 20
 420	   format(a1)
	   do 451 i=1,NPT-1	                      !Read the whole file
 430          continue                             !Jump here to read new card
	      read(ifile,fmt='(a80)',end=432) str1
	      goto 435
 432            continue             !end of file gets here
                close(ifile)
		if (ifile.eq.20) goto 452
	        ifile = ifile - 1
	      goto 430
 435          continue
	      call upcase(str1)
              idx = 1
	      call citem(str1,80,idx,str2,lstr2)
              if (str2(1:1).eq.'#') goto 430      !Ignore comment lines
	      if (str2(1:3).eq.'MJD') goto 430    !Ignore a commonly used...
	      if (str2(1:5).eq.'=====') goto 430  !   ...header format
	      if (str2(1:7).eq.'INCLUDE') then
		ifile = ifile + 1
		if (ifile.eq.30) then
		  write (*,*) 'Can''t nest INCLUDE''d time.dat files',
     +				' more than 10 deep.'
		  stop
	 	endif
	        call citem(str1,80,idx,str2,lstr2)
		path=clkdir(1:kc)//str2
	        open(ifile,file=path,status='old',err=900)
	        read(ifile,420)		              !Skip header lines
	        read(ifile,420)
                goto 430
	      endif

              if (str2(1:6).eq.'OFFSET') then

	        call citem(str1,80,idx,str2,lstr2)    ! Flexible format
                jsite(i) = sitea2n(str2)
	        call citem(str1,80,idx,str2,lstr2)
                read (str2,*) tdate(i)
	        call citem(str1,80,idx,str2,lstr2)
                read (str2,*) ckcorr(i)
	        call citem(str1,80,idx,str2,lstr2)
                ckflag(i) = 0
                if (lstr2.gt.0) then
                  if (str2(1:5).eq.'FIXED') ckflag(i) = 1
                endif

              else                                   ! Fixed format
              
  	        read(str1,1451) tdate(i),xlor,xjup,asite,dtmp
		
 1451	        format(f9.0,2f12.0,1x,a1,1x,a1)
	        jsite(i) = sitea2n(asite)
	        call upcase(dtmp)
	        if(xlor.gt.800.d0) xlor=xlor-818.8d0
	        ckcorr(i)=xjup-xlor
	        ckflag(i) = 0
  	        if (dtmp.eq.'F') ckflag(i) = 1  ! "fixed" -- no interpolation

              endif

 451	   continue
	   i=NPT
	   k=index(clkfile(1),' ')-1
	   write(*,'(''WARNING: '',a,'' too long, not all read'')')
     +  	clkfile(1)(1:k)
 452	   ndate=i-1
	endif

	if(nclk.ge.2)then                           ! NIST to UTC
              ! note: this conversion may be used backwards as a step in 
              !       converting ITOA-format TOAs (in UTC) to other formats
	   path=clkdir(1:kc)//clkfile(2)            
	   open(2,file=path,status='old',err=900)
	   read(2,1090)
	   read(2,1090)
	   do 23 i=1,NPT
	      read(2,1091,end=24) tutc(i),utcclk(i)
	      utcclk(i)=0.001d0*utcclk(i)
 23	   continue
	   k=index(clkfile(2),' ')-1
	   write(*,'(''WARNING: '',a,'' too long, not all read'')')
     +  	clkfile(2)(1:k)
 24	   tutc(i)=0.
	   close(2)
	endif

	if(nclk.ge.3) then                          ! NIST to other
	   path=clkdir(1:kc)//clkfile(nclk)
	   open(2,file=path,status='old',err=900) 
	   read(2,1090)
	   read(2,1090)
 1090	   format(a1)
	   do 90 i=1,NPT-1
	      read(2,1091,end=92) tbipm(i),bipm(i)
 1091	      format(f10.0,f19.0)
	      bipm(i)=0.001d0*bipm(i)
 90	   continue
	   i=NPT
	   k=index(clkfile(nclk),' ')-1
	   write(*,'(''WARNING: '',a,'' too long, not all read'')')
     +         clkfile(nclk)(1:k)
 92	   tbipm(i)=0.
	   close(2)
	endif

c  Beginning of iteration loop

 30     continue
        if (jits.eq.-1) jits=0
	if (.not.quiet) then
          write(31,'(/)')
          if(nits.gt.1)write(31,1038)jits+1,nits
 1038     format('Iteration',i3,' of',i3/)
        endif

	p0=1.d0/f0
	p1=-f1/f0**2
	if(si.gt.1.d0) si=1.d0

        if (.not.quiet) then
          write(31,1039) bmodel(nbin),nbin,nddm
 1039     format('Binary model: ',a,' nbin: ',i2,'   nddm:',i2)
          if(psrframe)write(31,'(''Parameters in pulsar frame'')')
          write (31,1040) psrname(1:index(psrname," ")-1)
 1040     format (/'Assumed parameters -- PSR ',a/)
        ENDIF

	if (eclcoord) then
          if (.not.quiet) 
     +         write (31,1042) pra*360.d0/TWOPI,pdec*360.d0/TWOPI
 1042     format ('LAMBDA:',f25.13/'BETA:',f27.13)
	else
          call radian (pra,irh,irm,rsec,123,1)
          call radian (pdec,idd,idm,dsec,1,1)
          decsgn=' '
          if(pdec.lt.0.) decsgn='-'
          idd=iabs(idd)
	  irs=rsec
	  rsec=rsec-irs
	  ids=dsec
	  dsec=dsec-ids
          if (.not.quiet) 
     +         write(31,1043) irh,irm,irs,rsec,decsgn,idd,idm,ids,dsec
 1043	  format ('RA: ',11x,i2.2,':',i2.2,':',i2.2,f9.8/
     +         'DEC:',10x,a1,i2.2,':',i2.2,':',i2.2,f9.8)
        endif

        if (.not.quiet) then
          if (eclcoord) then
            if(pmra.ne.0.)write(31,'(''PMLAMBDA (mas/yr):'',f18.4)')pmra
            if(pmdec.ne.0.)write(31,'(''PMBETA (mas/yr):'',f17.4)')pmdec
          else
            if(pmra.ne.0.0)write(31,'(''PMRA (mas/yr):'',f18.4)')pmra
            if(pmdec.ne.0.0)write(31,'(''PMDEC (mas/yr):'',f17.4)')pmdec
          endif
          if(px.ne.0.0)write(31,'(''Parallax (mas):'',f17.4)')px
          
          if(posepoch.gt.0.)then
            write (31,1100) posepoch
 1100       format('Pos Epoch (MJD):',f16.8)
          endif
          
          if(dmepoch.gt.0.)then
            write (31,1110) dmepoch
 1110       format('DM Epoch (MJD):',f16.8)
          endif
          
          write (31,1044) f0,p0,f1,p1*1.d15,f2,f3
 1044     format ('F0 (s-1): ',f22.17,6x,'(P0 (s):',f24.19,')'/
     +         'F1 (s-2): ',1p,d22.12,0p,6x,'(P1 (-15):',f22.12,')'/
     +         'F2 (s-3): ',1p,d22.9/'F3 (s-4): ',d22.6,0p)
          
          do i = 1, 5
            if (f4(i).ne.0.or.nfit(NPAR11+i+1).gt.0)
     +                       write(31,1045)i+3,-i-4,f4(i)
 1045       format ('F',i1,' (s',i2,'): ',1p,d22.9)
          enddo
          if (f4(6).ne.0.or.nfit(NPAR11+7).gt.0) 
     +                       write(31,1046)6+3,-6-4,f4(6)
 1046     format ('F',i1,' (s',i3,'): ',1p,d22.9)
          do i = 7, nfcalc
            if (f4(i).ne.0.or.nfit(NPAR11+i+1).gt.0)
     +                       write(31,1047)i+3,-i-4,f4(i)
 1047       format ('F',i2,' (s',i3,'): ',1p,d21.9)
          enddo
          write (31,1048) pepoch
 1048     format ('P Epoch (MJD):',f18.8)
          
          if(start.gt.0.)then
            write (31,1102) start
 1102       format ('Start MJD:',f22.8)
          endif
          if(finish.lt.100000.)then
            write (31,1104) finish
 1104       format ('Finish MJD:',f21.8)
          endif
          
          write (31,1106) dm
 1106     format ('DM (cm-3 pc):',f19.6)
          do i = 1, ndmcalc-1
            write (31,1120) i, dmcof(i)
          enddo
 1120     format ('DMCOF',i3.3,':',1p,d25.9)
          
          if(a1(1).ne.0.d0)then
            if(nbin.ne.9)then
	      write(31,1050) a1(1),e(1),t0(1),pb(1),omz(1)
            else
	      write(31,2050) a1(1),pb(1),t0asc,eps1,eps2
            endif
          endif
 1050     format('A1 sin(i) (s):',f18.9/'E:',f30.9/'T0 (MJD):',f23.9/
     +         'PB (d):',f25.12/'Omega0 (deg):',f19.6)
 2050     format('A1 sin(i) (s):',f18.9/'Pd (d):',f25.12/
     +         'TASC (MJD):',f21.9/'eps1:',f27.9/'eps2:',f27.9)
          
          if(omdot.ne.0.)write(31,'(''Omegadot (deg/yr):'',f14.6)')omdot
          if(xomdot.ne.0.)write(31,'(''XOMDOT (deg/yr):'',f16.3)')xomdot
          if(pbdot.ne.0.)write(31,'(''PBdot (-12):'',f20.3)')1.d12*pbdot
          if(xpbdot.ne.0.)write(31,'(''XPBDOT (-12):'',f19.3)')
     +         1.d12*xpbdot
          if(gamma.ne.0.)write(31,'(''Gamma (s):'',f22.6)')gamma
          if(si.ne.0.)write(31,'(''sin(i):'',f25.6)')si
          if(shapmax.ne.0.)write(31,'(''SHAPMAX:'',f25.6)')shapmax
          if(varsigma.ne.0.)write(31,'(''VARSIGMA:'',f25.6)')varsigma
          if(h3.ne.0.)write(31,'(''H3:'',f25.6)')h3
          if(h4.ne.0.)write(31,'(''H4:'',f25.6)')h4
          if(okom.ne.0.)write(31,'(''KOM (deg):'',f22.6)')okom*360./TWOPI
          if(okin.ne.0.)write(31,'(''KIN (deg):'',f22.6)')okin*360./TWOPI
          if(am.ne.0.)write(31,'(''M (solar):'',f22.6)')am
          if(am2.ne.0.)write(31,'(''m2 (solar):'',f21.6)')am2
          if(dr.ne.0.)write(31,'(''dr (-6):'',f24.3)')1.d6*dr
          if(dth.ne.0.)write(31,'(''dth (-6):'',f23.3)')1.d6*dth
          if(a0.ne.0.)write(31,'(''A0 (-6):'',f24.3)')1.d6*a0
          if(b0.ne.0.)write(31,'(''B0 (-6):'',f24.3)')1.d6*b0
          if(bp.ne.0.)write(31,'(''bp:'',f29.6)')bp
          if(bpp.ne.0.)write(31,'(''bpp:'',f28.6)')bpp
          if(xdot.ne.0.)write(31,'(''Xdot (-12):'',f21.6)')1.d12*xdot
          if(edot.ne.0.)
     +        write(31,'(''Edot (-12 s-1):'',f17.6)')1.d12*edot
          if(om2dot.ne.0.)write(31,'(''om2dot (rad s-2):'',5x,e10.4)')
     +         om2dot	                                     
          if(x2dot.ne.0.)write(31,'(''x2dot (s-1):'',10x,e10.4)')x2dot
          if(eps1dot.ne.0.)write(31,'(''eps1dot (-12 s-1):'',f14.6)')
     +         eps1dot*1.d12
          if(eps2dot.ne.0.)write(31,'(''eps2dot (-12 s-1):'',f14.6)')
     +         eps2dot*1.d12
          
          if(nplanets.gt.0)then
            do i=2,nplanets+1
              write(31,1051) i,a1(i),i,e(i),i,t0(i),i,pb(i),i,omz(i)
 1051         format('X(',i1,') (s):',f23.7/'E(',i1,'):',f27.9/
     +             'T0(',i1,') (MJD):',f20.9/'Pb(',i1,') (d):',f22.6/
     +             'Om(',i1,') (deg):',f20.6)
            enddo
          endif
          
          if(ngl.gt.0)then
            do i=1,ngl
              write(31,1053)i,glepoch(i),glph(i),glf0(i),
     +             glf1(i),glf0d(i),gltd(i)
            enddo
          endif
 1053     format('Glitch',i2/'  Epoch (MJD):',f18.6/
     +         '  dPHS:',f25.6/'  dF0 (s-1):',1p,d20.7/
     +         '  dF1 (s-2):',d20.7/
     +         '  dF0D (s-1):',d19.7,0p/'  TD (d):',f23.5)
        endif

	if (nbin.gt.0)then
	   do i=1,3
	      pb(i)=pb(i)*86400.d0
	   enddo
	endif
	k=0
	nfit(1)=1

	do 70 i=1,60				!Set up parameter pointers
	if(nfit(i).eq.0) go to 70
	k=k+1
	mfit(k)=i
70	continue

	if(nfit(16).ge.2) then			!Pointers to DM coeffs
	  do i=1,nfit(16)-1
	     k=k+1
	     mfit(k)=NPAR9+i
	     nfit(NPAR9+i)=1
	  enddo
	endif

	if(ngl.ne.0)then                        !Glitch parameters
	  do i=NPAR1+1,NPAR2
	    if(nfit(i).ne.0)then
	      k=k+1
	      mfit(k)=i
	    endif
	  enddo
	endif

	do i=NPAR3+1,NPAR9	!FB, XDOT, FBJ, DMX 
				!Note: NPAR9-NPAR10 handled abouve by nfit(16)
  	  if(nfit(i).ne.0) then
  	    k=k+1
	    mfit(k)=i
          endif
        enddo

	do i=NPAR10+1,NPA       ! FD terms (NPAR10+1 to NPAR11) and Fxx terms (NPAR11+1 to NPA)
	  if(nfit(i).ne.0) then
	    k=k+1
	    mfit(k)=i
	  endif
	enddo


	nparam=k        
	if (.not.quiet) write(31,1060) nparam
1060	format(/'Fit for',i3,
     +      ' parameters, including phase but excluding jumps')
c       store original values in case we need to read through TOAs twice
c       and re-do determination of new DMX and JUMP parameters
	nparam0 = nparam  
        ndmx0 = ndmx

	return

 900	write(*,'(''Failed to open clock correction file: '',a)')path
	STOP

	end

