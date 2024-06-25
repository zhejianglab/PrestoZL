c      $Id$
        program tempo

C  Analyzes pulsar time-of-arrival data to solve for period, first
C  and second derivatives, RA and DEC, proper motion, dispersion
C  measure, and orbital elements of binary pulsars.  Written by
C  J. H. Taylor, with major contributions from R. N. Manchester, 
C  J. M. Weisberg, D. J. Nice, M. F. Ryba, A. Irwin, and N. Wex.
C  Incorporates ephemeris code by E. M. Standish.

C  Recent version history:

C  Version  Date
C-----------------------------------------------------------------------------
C   7.10  ??-Oct-88  Implements H88 model for relativistic binaries
C   7.23  ??-Mar-89  Uses efficiently buffered I/O for resid1.tmp
C   7.30   2-Sep-89  When called by TZ, makes quick exit after writing
C		        computed phases to resid1.tmp
C   7.31   2-Jul-90  Fixed error interpolating nutation with JPL ephemeris
C   7.32  15-Jul-90  Changed print formats for P2, F2
C   7.33  20-Jul-90  EFAC and EMIN cmnds; DDGR prints expected DD params;
C                       date of run listed on banner
C   7.34  25-Jul-90  Better format for covariance matrix printout;
C                    very large covariances printed with full precision
C   7.35  27-Jul-90  Corrected slight errors in PMRA, PM galactic
C   7.36  29-Jul-90  Corrected constant in parallax definition (3% diff)
C   7.37  27-Aug-90  TZ shortcut saves ct2; new time.dat format uses MJDs
C   7.38  10-Sep-90  Up to 15 obsy sites; time.dat includes nsite
C   8.00  20-Sep-90  Cleaned up COMMONs, etc.
C   8.10  26-Sep-90  Incorporated functions of program tz, now obsolete
C   8.11   8-Oct-90  Help file when no arguments given; fixed GCOR printout
C   8.12  15-Nov-90  File name for NCLK=2 can be read from fit-for line
C   8.13   3-Dec-90  Overlap of +/- 5 min for poly coeffs in tzfit
C   8.14  14-Feb-91  Allows use of Parkes format TOA's.
C   8.15  10-Mar-91  Geocentric TOA's with nsite=0; MJDs listed on output.
C   8.16  12-Mar-91  New -g option gives GRO-format output 
C   8.20  13-Mar-91  All external Julian dates are now in MJD=JD-2400000.5;
C                       (old-style input parameters still OK, however).
C   8.21   6-May-91  Better UT1 extrapolation; X,Y,Z in obsys.dat; ET in ztim;
C                       half-day correction for calendar dates in tempo.lis
C   8.22   9-May-91  Projects p0,p1,p2 to t0geo, with output in fort.99
C   8.23  20-May-91  Fixed calculation of orb. phase over large no. of orbits
C   8.30  24-May-91  New relativistic binary models: 6=BT+ , 7=DDT=DD[T(lambda)]
C   8.31  11-Jun-91  Iterate "until done"; 7=DDT=DDT(bp,bpp) (bp,bpp fixed only)
C   8.32   1-Jul-91  Implemented "simulation" mode
C   8.33  15-Aug-91  PMRA fixed so that it works (and iterates) an high DEC
C   8.34   3-Sep-91  GRO option (-g) now forces writing of fort.99
C   8.35  21-Sep-91  Exits on Floating-point exception
C   8.36  20-Jul-92  Fixed incorrect constant in fk4<->fk5 routines
C   8.37  18-Aug-92  Increased no. of observatories; allow obs. X,Y,Z in TZ
C   9.00  12-Oct-92  General overhaul ... lots of additions
C   9.10   6-Nov-92  Allows non-phase-connected solutions
C   9.11  10-Nov-92  TRACK and SEARCH modes implemented
C   9.12  26-Nov-92  Better TRACK; reads ITOAF files; -x exports an ITOAF file
C   9.13  30-Nov-92  Up to 60 JUMP's allowed
C   9.14   1-Dec-92  New format for covariance matrix; Pbdot for 2 planets
C   9.15  10-Feb-93  Version control numbers integrated with SCCS
C   9.17  26-Mar-93  NCLK: 0=none; 1=time.dat; 2=bipm,ptb,or at1; 3=both 1 and 2
C   9.18  29-Mar-93  Read ut1.dat until the end of file
C   9.20   7-Apr-93  Use site "@" for barycentered input TOAs
C   9.21  13-Apr-93  Read and write ddm in itoa format
C   9.22  21-Apr-93  Up to 9-character field for terr
C   9.23  23-Apr-93  Outputs fpdot parameter when P2 is fit
C   9.24  30-Apr-93  Linearly interpolates clock offsets; for NCLK=6 subtracts
C			utccorr.tot and adds bipmfile (for itoa data)
C   9.25   5-May-93  Print warning about wmax, maximum weight assigned
C   9.26  24-May-93  Print out mass function for binaries
C   9.27  18-Jun-93  Allow number of coefficients to vary in JPL ephemeris
C   9.29  23-Jul-93  Allow MJD or d66 in "tz" input/output
C   9.30  30-Jul-93  Up to 70 characters allowed in "included" file names
C   9.31   5-Aug-93  Fixed minor bug in TRACK code in arrtim
C   9.32  11-Aug-93  JUMP ranges printed in MJD not d66.
C   9.33  16-Aug-93  Added EQUAD for systematic errors
C   9.35  17-Aug-93  New ITOA format (ktr94)
C   9.38  20-Aug-93  Fits polynomial DM=DM0 + DM1*t + DM2*t**2 + ...
C   9.38  20-Aug-93  Now you can iterate a polynomial DM fit
C   9.40  22-Aug-93  Big speedup; no resid1.tmp (unless -r); 3 orbits; bug fix
C  10.002 24-Aug-93  Rational selection of reference clock and ephemeris
C  10.004 25-Aug-93  Option for bootstrapped parameter uncertainties
C  10.008 20-Dec-93  Parkes format change, clock interp; planets in tempo -z
C  10.010 27-Dec-93  Allow universal time offsets for a given site in time.dat
C  10.011  3-Feb-94  Fix 'wait' calculation in atimfake; new '-v' option
C  10.012  7-Mar-94  clockcor.f modified slightly for Parkes data
C  10.014 12-Jun-94  Fix up itoa output
C  10.015 29-Jun-94  Glitches: fit for phase, freq, fdot and exp decay (RNM)
C  10.018  7-Aug-96  Interpol. time.dat values on days w/multiple entries (DJN)
C  10.019  6-Nov-96  Up to 9999 tztot.dat sources ("ipsr" loop in tempo) (DJN)
C  10.020 21-NOV-96  In multi-orbit fit use pbdot, omdot in 1st orb only (DJN)
C  10.020  3-FEB-97  -z: allow nbin=1 or 3 for any pulsar (DJN)
C  10.021 12-FEB-97  -z: user-specified freq.; put doppler in polyco.dat (DJN)
C  11.000 18-JUL-97  Optional free-format input, pulse frequencies as input,  
C                    Alan Irwin's improvements to lmst, geodetic to geocentric
C                    coords, tdb argument to JPL ephemeris. (RNM)
C  11.001 20-JAN-98  Ecliptic coordinates, Year 2000 fixes (DJN)
C  11.002 11-MAR-98  Separate pos'n epoch; add START/FINISH (RNM)
C                    Added "ell1" binary model (NW)
C  11.005 24-MAR-00  Add BNRYBTX model; clean up parameter indexing;
C                    add "-c" option to iterate even when seem converted;
C                    other fixes (DJN)

C Logical units			Opened by
C----------------------------------------------------------------
C   2	tempo.cfg tempo.hlp leap.sec obsys.dat time.dat ptbnist.90
C				newsrc, setup, tempo
C  11	tz parameters (tz.in)   tzinit
C  12	tz psr params           tpohdr
C  13   polyco.dat		tzinit
C  20.29 time.dat & INCLUDE'd   newsrc
C  30	resid1.tmp		tempo
C  31	tempo.lis		tempo
C  32	resid2.tmp		tempo
C  33	gro.1			tempo
C  34	gro.2			newbin
C  35   pulse number file       tempo
C  36   info.tmp                tempo/arrtim
c  37   design.tmp              tparin
C  38   phisun.tmp              tempo
C  39   <dcovfile>              glsfit
C  40   dmxn.tmp                tempo/arrtim
C  41   doppler.tmp             tempo
C  42	ut1.dat			tempo
C  43   TDB-TDT ephemeris       tempo/tdbinit
C  44   BC ephemeris            newsrc/ephinit
C  45	itoa.out		tempo
C  49   <infile>.par            tempo
C  50   <infile>		tempo
C  51,... INCLUDE files		arrtim
C  71	params.tmp		tempo
C  72	matrix.tmp		mxprt
C  99	gro.99			newval

        implicit real*8 (a-h,o-z)

	include 'dim.h'
	include 'acom.h'
        include 'array.h'
	include 'bcom.h'
	include 'clocks.h'
	include 'config.h'
	include 'dp.h'
        include 'orbit.h'
	include 'vcom.h'
	include 'eph.h'
	include 'trnsfr.h'
	include 'tz.h'
        include 'toa.h'

	logical lw, nostop
        logical memerr
	character*640 infile
        character*80 ut1file,resfile1,
     +       resfile2,listfile,fname,line,tdbfile,s,hlpfile
        character*640 obsyfile
        character*640 path
        character*160 cfgpath
	character*160 npulsefile, infofile, phisunfile, dmxnfile
        character*160 dopplerfile
	character date*9,date2*9,damoyr*9,label*12,parfile*160
	integer time
        real*8 xmean(NPA),alng(NOBSMAX)

        integer sitea2n ! external function

	common/leapsec/mjdleap(50),nleaps
	data resfile1/'resid1.tmp'/
        data infofile/'info.tmp'/
        data phisunfile/'phisun.tmp'/
        data dopplerfile/'doppler.tmp'/
        data dmxnfile/'dmxn.tmp'/
	data listfile/'tempo.lis'/
	data bmodel /'None','BT','EH','DD','DDGR','H88','BT+','DDT',
     +       'MSS','ELL1','BTX','BT1P','BT2P','DDS','DDK','DDFWHE',
     +       'ELL1H'/

        resfile2 = 'resid2.tmp'//char(0)

	version = CONFIGVERSION

	memerr = .false.
        infoout = .false.
        stflag = .false.
        tzsitedef = ' '
        tdbif99fmt = .false.

c  Get command-line arguments

        call tparin(nostop,lw,lpth,nparmax,nptsmax,version,
     +     npulsefile,infile,path,resfile1,hlpfile,parfile)

c  Parse tempo.cfg file

 	cfgpath=path(1:lpth)//'/tempo.cfg'
        call cfgin(cfgpath,ut1file,obsyfile,tdbfile)  

c  Open primary output file (tempo.lis)
        if (.not.quiet)
     +       open(31,file=listfile,status='unknown')

	call setup(version,infile,obsyfile,alng,parfile)

	nfl=index(infile,' ')-1

        nbuf = nptsmax * (nparmax+8) 
        call tmalloc(nptsmax,nbuf) ! allocate large arrays

        if (tz.and.tzsite.eq.' '.and.tzsitedef.ne.' ') then
          tzsite=tzsitedef
          nsite=sitea2n(tzsite)
        endif

c  Open leap second file, read leap second dates, close leap second file
	k=index(clkdir,' ')-1
	path=clkdir(1:k)//'leap.sec'
	lpth=index(path,' ')-1
	open(2,file=path(1:lpth),status='old',err=7)
	do 6 i=1,50
 6	read(2,*,end=8) mjdleap(i)
	print*,'Too many leap seconds'
	go to 9999
 7	print*,'Cannot open ',path(1:lpth)
	go to 9999
 8	nleaps=i-1
	close(2)

c  Open ut1 file (if present)
	path=clkdir(1:k)//ut1file
        open(42,file=path,status='old',err=10)
	ut1flag=.true.
	go to 11
 10     continue
        lpth=index(path,' ')-1
   	write(*,1005)path(1:lpth)
	write(31,1005)path(1:lpth)
 1005	format(' **** Warning - Cannot open UT1 file at ',a/
     +    ' No UT1 correction'/)
	ut1flag=.false.

 11	continue
	if (npulsein) then
	  open(35,file=npulsefile,status='old')
	else if (npulseout) then
	  open(35,file=npulsefile,status='unknown')
        endif

c  Open TDB-TDT clock offset file
	k = index(ephdir,' ')-1
	n = index(tdbfile,' ')-1
	path = ephdir(1:k)//tdbfile(1:n)
        if (tdbif99fmt) then
    	  call tdbinit2(43,path)
        else
    	  call tdbinit(43,path)
        endif

	if (tz) then  ! generate predictive ephemeris (polyco.dat),like old TZ

	  jits = 0

          if (autotz) then
            num = 1             ! always one pulsar
            parunit = 49        ! logical unit number used for .par file
            if (tzmjdstart.lt.0) then
              fmjdnow=40587+time()/86400.d0
              tzmjdstart = fmjdnow
            endif
            fmjd1 = tzmjdstart
            fmjd2 = tzmjdstart
            if (polycofile.eq.'-') then
              lupolyco = 6
            else
              lupolyco = 13
              open (lupolyco,file=polycofile,status='unknown')
            endif

          else  !.not. autotz   (use tz.in file)

            if(infile.eq.'') then
              tzfile='tz.in'
            else
              tzfile=infile
            endif
            infile='tz.tmp'
            
            if(.not.oldpar)parunit=49
            open(50,file=infile,status='unknown')
            
            call tzinit(alng,sitelng,num)
            fmjdnow=40587+time()/86400.d0
            date=damoyr(int(fmjdnow))
            
            if (.not.quiet) then
              write(*,1010) date,fmjdnow
 1010         format(/' Current date is ',a9,', or MJD',
     +             f10.3//' Enter first and last MJD,',
     +             ' or hit return to run for today: ')
            endif
            read(*,fmt='(a80)') line
            do j=80,1,-1
              if(line(j:j).ne.' ')then
                read (line,*,err=710) fmjd1,fmjd2
                go to 710
              endif
            enddo
	    fmjd1=0.d0
	    fmjd2=0.d0
 710	    continue
	    if(fmjd1.eq.0.d0) fmjd1=fmjdnow
	    if(fmjd2.eq.0.d0) fmjd2=fmjd1+0.5d0
          endif
	  date=damoyr(int(fmjd1))
	  date2=damoyr(int(fmjd2))
	  if (.not.quiet) write(*,1012) date,date2
 1012	  format(1x,a9,' through ',a9/)
 
	  do ipsr=1,num
	    start = 0.		! for ipsr>1, these will have been set by
	    finish = 0.		!   arrtime to values for the previous pulsar
            call tpohdr(oldpar,pardir,parfile,ncoord,t0,pb,p0,dm,nbin,
     +           ipsr)
	    if(start.ne.0. .and. 
     +           (fmjd1.lt.start.or.(ntzrmjd+ftzrmjd).lt.start))then
	       if(usestart)then
		  write(*,1014)
		  STOP
	       else
		  if (.not.quiet) write(*,1016)
	       endif
	    endif
	    if(finish.ne.0. .and. 
     +           (fmjd2.gt.finish.or.(ntzrmjd+ftzrmjd).gt.finish))then
	       if(usefinish)then
		  write(*,1014)
		  STOP
	       else
		  if (.not.quiet) write(*,1016)
	       endif
	    endif
 1014	    format(' ERROR: Requested MJD or TZRMJD outside parameter ',
     +         'validity range')
 1016	    format(' WARNING: Requested MJD or TZRMJD outside parameter ',
     +         'validity range')

	    if (name(ipsr).eq.'done') then
	      tsid=1.d0/1.002737909d0
	     
	      do afmjd=fmjd1,fmjd2,tsid
	  	call atimfake(afmjd,nbin,nt,sitelng,ipsr)
		rewind 50
		close(2)
                ! following moved earlier in code 2017-Aug-10
                ! left here temporarily in case we find a problem with it
		! call setup(version,infile,obsyfile,alng,parfile)
		call newsrc(nits,jits,nboot)
		call arrtim(mode,xmean,sumdt1,sumwt,dnpls(1+dnplsoff),
     +               ddmch(1+ddmchoff),ct2,alng,nz,nptsmax,
     +               nits,jits,buf(1+bufoff),npmsav(1+npmsavoff),
     +               ksav(1+ksavoff),nbuf,memerr,infofile,lw)

		call tzfit(ipsr,dnpls(2+dnplsoff),f0,dm,ct2,t0,pb)
		rewind 31
		close(44)
		rewind 50
	      enddo
	    endif
	    if(.not.oldpar)close(49)
	  enddo

          do 890 ipsr=1,num
            if (name(ipsr).ne.'done')then
	       if(parfile.eq.'def')then
		  print *,'PSR ',name(ipsr),' not in data base'
	       else
		  print *,'Parameter file not found: ',parfile
	       endif
	    endif
 890      continue
          close(21)
          if (.not.quiet) close (31)
          close(50)

	else  ! Standard TEMPO execution

c  Open parameter and residual files
	   open(50,file=infile,status='old',err=9997)
	   if (.not.oldpar) then                 ! free-form parms....
 900	     continue
	     read(50,fmt='(a)',end=910) line     ! parms in same file as toas?
	     j=1                       
	     call citem(line,80,j,s,k)
	     call upcase(s) 
	     if (s(1:4).eq.'HEAD') then          ! yes
	       parunit = 50
	     elseif ((s(1:1).eq.'C'.and.k.eq.1).or.s(1:1).eq.'#'.or.
     +		s(1:35).eq.'                                   ') then ! maybe
	       goto 900
	     else                                ! no:
		if(parfile.eq.'def')then
		   n=index(infile,'.')-1 !     open par file
		   if(n.lt.0)n=index(infile,' ')-1
		   parfile=infile(1:n)//'.par'
		endif
		open(49,file=parfile,status='old',err=9990)
		parunit = 49
	     endif
	     rewind(50)
	  endif

	  goto 911
910	  print *,"Input file is empty (except possibly for comments)"
	  goto 9999
 911	  continue

	  jits=0

	  if(oldpar.or.parunit.eq.50)then
             if (.not.quiet) write(*,1050) version,infile(1:nfl)
 1050	     format(' TEMPO v ',f6.3,
     +          ' Princeton/ATNF Pulsar Collaboration'/' Data from ',a)
	  else
             if (.not.quiet) write(*,1051) version,infile(1:nfl),parfile
 1051	     format(' TEMPO v ',f6.3,
     +       ' Princeton/ATNF Pulsar Collaboration'/
     +       ' Data from ',a, ',   Input parameters from ',a)
	  endif

          ! following moved earlier in code 2017-Aug-10
          ! left here temporarily in case we find a problem with it
	  ! call setup(version,infile,obsyfile,alng,parfile)

C         The main loop:
 60       continue

	  if (phisunout) open (38,file=phisunfile)
	  if (dmxnout) open (40,file=dmxnfile)
	  if (dopplerout) open (41,file=dopplerfile)
          call newsrc(nits,jits,nboot)

 62       continue  ! re-entry point after re-allocating arrays
          rewind 32
          call arrtim(mode,xmean,sumdt1,sumwt,dnpls(1+dnplsoff),
     +         ddmch(1+ddmchoff),ct2,alng,nz,nptsmax,nits,jits,
     +         buf(1+bufoff),npmsav(1+npmsavoff),ksav(1+ksavoff),nbuf,
     +         memerr,infofile,lw)

	  if (memerr) then
            call tfree()
            if (.not.quiet) then
	      write (*,1060)
 1060	      format ("Warning: arrays too small on first pass.")
	      write (*,1061) nptsmax, nparmax
 1061	      format ("  Current settings:        -m ",i7," -l ",i7)
	      write (*,1062) n, nparam
 1062	      format ("  Better to run:     tempo -m ",i7," -l ",i7)
	      write (*,1063)
 1063	      format ("  Re-allocating arrays and running again.")
            endif
            nptsmax = n
            nparmax = nparam+8
            memerr = .false.
            nbuf = nptsmax * (nparmax+8) 
            call tmalloc(nptsmax,nbuf) ! allocate large arrays
            nparam = nparam0    ! un-do any change to nparam from arrtim call
            ndmx = ndmx0        ! un-do any change to ndmx from arrtim call
            goto 62
          endif

	  wmax=0.0

          if (useglsfit) then
            print *,""
            print *,'Using GLS fit (this may take a while)'
            call glsfit(n,mode,chisqr,varfit,xmean,sumdt1,sumwt,nz,wmax,
     +           lw,ddmch(1+ddmchoff),
     +           buf(1+bufoff),npmsav(1+npmsavoff),ksav(1+ksavoff),
     +	         resfile2)
          else
            call fit(n,mode,chisqr,varfit,xmean,sumdt1,sumwt,nz,wmax,
     +           lw,ddmch(1+ddmchoff),
     +           buf(1+bufoff),npmsav(1+npmsavoff),ksav(1+ksavoff),
     +	         resfile2)
          endif

          asig=sqrt(varfit)
	  if(nboot.gt.0)
     +         call bootmc(n,mode,nz,nboot,nparam,mfit,freq,ferr,
     +         ddmch(1+ddmchoff),
     +         buf(1+bufoff),npmsav(1+npmsavoff),ksav(1+ksavoff),
     +         resfile2)
          if (usedmdata.and.useglsfit) then
            call newval(chisqr*(n-nz-nparam)/(2*n-nz-nparam),
     +       2*n-nz-nparam,rms0,rms1,nits,jits,wmax,nboot)
          else
            call newval(chisqr,n-nz-nparam,rms0,rms1,nits,jits,wmax,
     +                  nboot)
          endif
	  if(abs(rms1 - rms0).le.max(1.d-4*abs(rms0),0.1d0).and.
     +		.not.nostop) go to 9999

          if (phisunout) close(38)
          if (dopplerout) close(41)

          if(jits.lt.nits .or. nits.eq.9) go to 60

	endif
	go to 9999

 9990   write(*,'(/''Error opening '',a/)')
     +            parfile(1:index(parfile,' ')-1)
        goto 9999

 9997	write(*,'(/''File '',a,'' not found'')')infile(1:nfl)
        goto 9999

 9998	call system('more '//hlpfile)

 9999	continue
	end


