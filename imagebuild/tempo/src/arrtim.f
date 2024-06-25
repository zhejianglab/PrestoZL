c      $Id$
      subroutine arrtim(mode,xmean,sum,sumwt,dnpls,ddmch,ct2,
     +     alng,nz,nptsmax,nits,jits,
     +     buf,npmsav,ksav,nbuf,memerr,infofile,lw)

C  Reads pulse arrival times (at observatory, in UTC) from unit 4. Then
C  computes equivalent coordinate times at solar system barycenter,
C  and computes residuals from the assumed set of parameters.

C  RNM 27-Feb-91  Allow nsite=0 for geocentric data
C  DJN 18-Aug-92  Allow up to 36 sites

	implicit real*8 (a-h,o-z)
        save
	include 'dim.h'
	real*8 xmean(NPA),fctn(NPAP1),dnpls(*),ddmch(*),alng(36)
        real*8 buf(*)
        real*8 toff
        integer npmsav(*), ksav(*)
	real*4 gasdev
        integer ZEROWT(2000),idum
	integer ilen
        integer wflag
        integer nread
	character*640 card,card2
        character*160 infile,aitem
	character asite*1,bsite*2,comment*8,aterr*9,afmjd*15
        character*20 amjd
	logical first,offset,jdcatc,last,dithmsg,track,search
	logical memerr
	logical parsed
	integer i1, i2
        character*160 infofile
        character*80 tmp
        character*640 rawflags
	include 'acom.h'
	include 'bcom.h'
	include 'dp.h'
	include 'trnsfr.h'
	include 'orbit.h'
	include 'glitch.h'
	include 'tz.h'
        include 'toa.h'
        integer sitea2n ! external function
        character*1 siten2a ! external function
        character*2 siten2b ! external function
        character*80 getvalue ! external function

	common /CONST/ PI,TWOPI,SECDAY,CONVD,CONVS,AULTSC,VELC,EMRAT,OBLQ,
     +              GAUSS,RSCHW,AULTVL
	common/crds/ rcb(6),rbe(6),rce(6),rse(3),rca(3),rsa(3),rba(3),
     +    rea(3),psid,epsd,pc,ps,tsid,prn(3,3),atut,ut1ut,etat,dtgr,
     +    tdis,bclt
	common/obsp/ site(3),pos(3),freqhz,bval,sitvel(3)
	data first/.true./,card2/'JUMP'/,idum/-1/

	offset=.true.
	jdcatc=.false.
	do 1 j=1,NPA
1	xmean(j)=0.

        if (npulsein.or.npulseout) rewind 35
        if (dmxnout) rewind 40
        if (.not.stflag) then
          rewind 50
	  do i=1,nskip
	    read(50,1002)
 1002	    format(a1)
	  enddo
        endif

	dphase=0.d0
        dphaseflag=0.d0
        dphasetot=0.d0
	deltat=0.d0
        pha1 = 0.d0
        pha2 = 0.d0
	equad  = 0.d0
	equadsave = 0.d0
	emin = 0.d0
	efac = 1.d0
	efacsave = 1.d0
        dmefac = 1.d0
        dmefacsave = 1.d0
        dmequad = 0.d0
        dmequadsave = 0.d0
        dmjump = 0.d0
        dmjumpsave = 0.d0
	sigm = 0.d0
	amjd1=1000000.0
	amjd2=0.
	last=.false.
	dithmsg=.true.
	track=.false.
	search=.false.
	ct00=-1.d30
	fmin=0.
	fmax=1.d30
	emax=1.d30
	nblk=0
	ntrk=0
	dither=0.
        zawgt=0.
	mode=fitmode
        dmobs=0.
        dmobserr=0.
	nxoff=nflagjumps
	do i=1,NJUMP
	  nfit(NPAR2+i)=0
	enddo
	n=0
	sum=0.
	sumsq=0.
	sumwt=0.
	modscrn=10
        nw=1
        nread = 0
	if(first)then
	   call ztim(0,0.d0,nct,fct,wflag)
	   if (.not.quiet)write (31,1011)
 1011	   format(/'    N     TOBS   WORD    NJMP    WGT',
     +     '     DTIME       TIME     DPHASE      PHASE'/)
	   nz=0
	   ntzref=0
	   first = .false.
	endif
	lu = 50
        nfmt = 0

c       convert from mas/year (with cos(dec) term in) to radian/Julian century
c       in ra (without cos(dec) term), dec and rv.
	dpa = convd/3600.d0*pmra/10.d0/dcos(pdec)
	dpd = convd/3600.d0*pmdec/10.d0
	dpr = convd/3600.d0*pmrv/10.d0

c       parallax in radians from parallax in mas
	dpara = px*convd/3600.d3

	if(psrframe)then
C       A. Irwin stuff requires positions at epoch J2000. Compute space motion vectors 
C       and J2000 epoch position (Note all of this is in J2000 coordinates).
	   delt = (51545.d0 - posep)/36525.d0
	   pra2000 = pra + dpa*delt
	   pdec2000 = pdec + dpd*delt
c       radial velocity in au/Julian century
	   if(dpara.gt.0.d0)then
	      dvel = dpr/dpara
	   else
	      dvel=0.d0
	   endif
	   call space_motion(pra2000, pdec2000, dpa, dpd, dpara, dvel, 
     +   	aa_q, aa_m, 1)
c       aa_q/dpara is position vector in au at epoch J2000
c       aa_m/dpara is space velocity vector in au/Jcentury

	else
c       N. Wex transformations at pepoch
	   call pospm_conv(pra, pdec, dpa, dpd, aa_q, aa_m)
	endif

	comment='        '

C       The main loop starts here!

   10  	continue

        if (stflag)then  ! READ TOAS FROM STORAGE ARRAYS

          nread = nread + 1
          if (nread.gt.stntoa) goto 45

          nsite = stnsite(nread)
          rfrq = stfrq(nread)
          nfmjd = stnmjd(nread)
          ffmjd = stfmjd(nread)
          terr = sterr(nread)
          ddm = stddm(nread)
          asite = siten2a(nsite)

        else   ! READ TOAS (AND OTHER CARDS) DATA FROM FILE

          read(lu,1010,end=40) card
 1010     format(a640)
          if(card(1:1).lt.'A'.and.card(1:1).ne.'#') go to 50
          if((card(1:1).ge.'a').and.(card(1:1).le.'z')) go to 50
          if(card(1:4).eq.'END ') go to 45
          if(card(1:7).eq.'INCLUDE') then
            lu=lu+1
            if (lu.gt.70) then
              write (*,*) "INCLUDE statements nested too deeply;",
     +                    " maximum 21 levels nesting allowed."
              STOP
            endif
            j1 = 8
            call citem(card,640,j1,infile,ilen)
            if (.not.quiet)write(31,1012) card(1:78)
 1012       format(1x,a78)
            open(lu,file=infile(1:ilen),status='old',err=1013)
            rewind lu
            nfmt = 0
          else
            call pcard(card,mode,zawgt,deltat,fmjd,dphase,sigm,offset,
     +           jdcatc,pha1,pha2,efacsave,emin,equadsave,jits,lu,track,
     +           trkmax,search,lw,nfmt,parsed)
	    if ((.not.parsed).and.(nfmt.eq.3)) go to 50
          endif
          go to 10
          
 1013     write(*,'(''Failed to open INCLUDE file: '',a)')infile
          STOP
          
 50       continue
          if(card(1:35).eq.'                                   ') 
     +         goto 10

c  If nfmt is 3 then we've seen a FORMAT 1 line in this file
          if(nfmt.eq.3) goto 51
          
c  default: Princeton format
          nfmt=0
c  first character blank: Parkes format
          if(card(1:1).eq.' ') nfmt=1
c  first character not blank, second character not blank:
c    could be either ITOA format, in which case first field is PSR name
c    or Princeton format, but with "scan number" starting at column 2
c    (Berkeley sometimes does this).  Differentiate between the 
c    two cases by searching for the "+" or "-" indicative of a pulsar name.
          if(card(1:1).ne.' ' .and. card(2:2).ne.' ') then
            do i = 1, 80
              if(card(i:i).eq." ") then
                nfmt=0
                goto 51
              elseif (card(i:i).eq."+") then
                nfmt=2
                goto 51
              elseif (card(i:i).eq."-") then
                nfmt=2
                goto 51
              endif
            enddo
            nfmt = 2            ! shouldn't get here....but just in case
          endif
 51       continue

c Apply the default equad and efac settings.  These may be altered
c by flag-based settings for tempo2 TOAs
          efac = efacsave
          equad = equadsave
          dmefac = dmefacsave
          dmequad = dmequadsave
          dmjump = dmjumpsave
c blank out temp flag variable
          rawflags = ''

          if(nfmt.eq.0) then				! Princeton format
            read(card,10500) asite,rfrq,amjd,aterr,ddm
10500       format(a1,14x,f9.0,a20,a9,15x,f10.0)
            nsite = sitea2n(asite)

            ! Clean up MJD and separate out integer and
            ! fractional parts.  It can have either of these forms:
            ! " xxxxx.xxxxxxxxxxxxx" or "xxxxx.xxxxxxxxxxxxxx"
            ! Sometimes it has a space and then junk, such as:
            ! "xxxxx.xxxxxxxxxx xxx". anything after the space
            ! should be ignored.
            i1 = index(amjd,'.')
            i2 = index(amjd(i1+1:20),' ') 
            if (i2.eq.0) then
              i2 = 20
            else
              i2 = i1+i2
            endif
            read (amjd(1:i1-1),*) nfmjd
            read (amjd(i1:i2),*) ffmjd

            if(nfmjd.lt.30000)nfmjd=nfmjd+39126       ! Convert 1966 days to MJD

            iz = 1
            do i=1,9
              if(aterr(i:i).eq.' ') iz=i
            end do
            terr=0.
            if (iz.le.9) read(aterr(iz:9),*,err=54,end=54) terr
 54         continue

          else if(nfmt.eq.1) then ! Parkes format
            read(card,1051) rfrq,nfmjd,ffmjd,phs,terr,asite
 1051       format(25x,f9.0,i7,f14.13,f8.6,f8.1,8x,a1)
            ffmjd=ffmjd+phs*p0/86400.d0
            nsite = sitea2n(asite)
            
          else if(nfmt.eq.2) then ! ITOAF format
            read(card,1052) nfmjd,ffmjd,terr,rfrq,ddm,bsite,comment
 1052       format(9x,i5,f14.13,f6.2,f11.4,f10.6,2x,a2,2x,a8)
            nsite = sitea2n(bsite)
            asite=' '
            if (nsite.eq.-2) then
              print *,"Error: no such obervatory code as: ",bsite
              print *,"Problem TOA line:"
              print *,card
              stop
            endif

          else if(nfmt.eq.3) then ! TEMPO2 format
c First 5 fields are file, freq, TOA, err, site
c Then everything after that are flags (ignored for now)
            j1 = 1
            call citem(card,640,j1,aitem,ilen) ! File, ignore it
            call citem(card,640,j1,aitem,ilen) ! Freq
            read (aitem,*) rfrq
            call citem(card,640,j1,aitem,ilen) ! TOA
            i1 = index(aitem,'.')
            read (aitem(1:i1-1),*) nfmjd
            read (aitem(i1:ilen),*) ffmjd
            call citem(card,640,j1,aitem,ilen) ! err
            read (aitem,*) terr
            call citem(card,640,j1,aitem,ilen) ! site
            if (ilen.eq.1) then
              asite=aitem(1:1)
              nsite = sitea2n(asite)
            else 
              asite=' '
              bsite=aitem(1:2)
              nsite=sitea2n(bsite)
              if (nsite.eq.-2) then
                print *,"Error: no such obervatory code as: ",bsite
                print *,"Problem TOA line:"
                print *,card
                stop
              endif
            endif

 55         continue
            rawflags = card(j1:640)
            call getflags(card,640,j1)
            tmp = getvalue("to")
            if (tmp.ne."") then
              read(tmp,*) toff
              call mjdadd(nfmjd,ffmjd,toff)
            endif
            dphaseflag = 0.0d0
            tmp = getvalue("padd")
            if (tmp.ne."") then
              read(tmp,*) dphaseflag
            endif
            dmobs = 0.0d0
            tmp = getvalue("pp_dm")
            if (tmp.ne."") then
              read(tmp,*) dmobs
            endif
            dmobserr = 0.0d0
            tmp = getvalue("pp_dme")
            if (tmp.ne."") then
              read(tmp,*) dmobserr
            endif
            tmp = getvalue("info")
            if (tmp.ne."") then
              ! Fill in the INFO stuff -- tempo2 flag case
              infotxt = tmp
              infolen = index(tmp,' ') - 1
            endif
            if (infoflag.ne."") then
              tmp = getvalue(infoflag(2:32))
              if (tmp.ne."") then
                ! Fill in the INFO stuff
                infotxt = tmp
                infolen = index(tmp,' ') - 1
              else
                ! Flag not found, set an empty info line
                infotxt = ""
                infolen = 0
              endif
            endif
            ! Here we check for jump-related flags and do 
	    ! the right thing. If this TOA is supposed to be
	    ! "JUMPed" we need to:
	    !  - set nfit(npar2+ijump)=1 if the jump is to be fit
            !  - set x(npar2+ijump)=-1 (partial deriv array)
	    ! Also if this is the first TOA in this JUMP segment, and
            ! the jump is to be fit:
	    !  - increment nparam
	    !  - set mfit(nparam) = npar2 + ijump
	    ! These last two might be more naturally done elsewhere if
	    ! we are ok with empty JUMPs crashing things.. 
	    ! Note, nxoff should be the total number of jumps including
	    ! both flag-based and original-style.
            do i=1,nflagjumps
	      x(NPAR2+i) = 0.0d0 ! Default to not applying this jump
              tmp = getvalue(jumpflag(i)(2:32))
              if (tmp.eq.jumpflagval(i)) then
                x(NPAR2+i) = -1.0d0
                if (.not.nofitjump(i)) then
                  if (nfit(NPAR2+i).eq.0) then ! 1st TOA in JUMP
                    nparam = nparam + 1
                    mfit(nparam) = npar2 + i
                  endif
                  nfit(NPAR2+i) = 1
                endif
              endif
            enddo
	    ! Check for flag-based equad/efac and apply as
	    ! appropriate.  Currently the following behavior happens:
            !
            !   - Flag-based efac/equad values override 'old-style'
            !   efac/equad
            !
            !   - If multiple efac/equad apply to this TOA, the value
            !   that gets used depends on ordering (so avoid this!).
            !
            !   - For TOAs where no matching flags are found, any
            !   currently-set old-style EFAC/EQUAD apply.
            do i=1,nflagefac
              tmp = getvalue(efacflag(i)(2:32))
              if (tmp.eq.efacflagval(i)) then
                efac = flagefac(i)
              endif
            enddo
            do i=1,nflagequad
              tmp = getvalue(equadflag(i)(2:32))
              if (tmp.eq.equadflagval(i)) then
                equad = flagequad(i)
              endif
            enddo
            do i=1,nflagdmefac
              tmp = getvalue(dmefacflag(i)(2:32))
              if (tmp.eq.dmefacflagval(i)) then
                dmefac = flagdmefac(i)
              endif
            enddo
            do i=1,nflagdmequad
              tmp = getvalue(dmequadflag(i)(2:32))
              if (tmp.eq.dmequadflagval(i)) then
                dmequad = flagdmequad(i)
              endif
            enddo
            do i=1,nflagdmjump
              tmp = getvalue(dmjumpflag(i)(2:32))
              if (tmp.eq.dmjumpflagval(i)) then
                dmjump = flagdmjump(i)
              endif
            enddo
          endif

 56       if(ffmjd.gt.1.d0) then
            nfmjd=nfmjd+1
            ffmjd=ffmjd-1.d0
          endif
          ! nsite = sitea2n(asite)  ! this is now done earlier

        endif  ! end of read-data-from-file-or-memory if statement

	fmjd=nfmjd+ffmjd
	if(rfrq.lt.fmin .or. rfrq.gt.fmax .or. 
     +      (terr .gt. emax) .or.
     +      (usestart .and. fmjd.lt.start) .or. 
     +      (usefinish .and. fmjd.gt.finish)) go to 10

C Arrival time
	n=n+1

        if (n.gt.nptsmax) then
	  if (tz) stop ' Too many points for tempo -z'
          memerr = .true.  !too many points.  continue through the file
	  goto 10          !to figure out how many points total we will have
        endif

C Store toa for tz reference phase
C Include any TIME offset (deltat)
C It is OK if ftzrmjd is a little outside the range [0,1],
C this will be accounted for when the final value is calcualted
C and printed out.
	if(ntzref.eq.0.and..not.tz)then
	   if(fmjd.gt.pepoch) ntzref = n
           ! either:
           !   ntzref has just been set by the above if statement,
           !   in which case we want to store the tzr info
           ! or:     
           !   ntzref has not yet been set (so ntzref==0), in
           !   which case it may never be set within this
           !   part of the code, in which case we want to save
           !   the tzr info in case this is the last TOA and
           !   we set ntzref later.
	   ntzrmjd=nfmjd
	   ftzrmjd=ffmjd + deltat/86400.d0
	   ntzrsite=nsite
	   tzrfrq=rfrq
	endif

C  Get clock corrections
	if(nsite.eq.-1)then         ! no correction for sites 0 and @
	   clk2=deltat/86400.d0      
	else if (nsite.eq.0) then
                                    ! this hack tells clockcor to treat the
                                    ! TOA as UTC.  clockcor can convert it
                                    ! to TT(BIPM) if requested in the tempo header
	   call clockcor(fmjd,nsite,n,deltat,clk2,2)
	else 
	   call clockcor(fmjd,nsite,n,deltat,clk2,nfmt)
	endif
	
c Loop and apply all jumps where x(NPAR2+i).ne.0
        if(.not.jumpbarycenter.and.nxoff.gt.0) then
	  do i=1,nxoff
	    if (x(NPAR2+i).ne.0.d0) clk2=clk2+dct(i)/86400.d0
	  enddo
        endif

	if(dither.ne.0.d0 .and. (.not.sim)) then
	  clk2=clk2+dither*gasdev(idum)/86400.d6
	  if(nits.ge.2.and.dithmsg) then
	    write(*,*) 'Iterating a solution with nonzero ',
     +        'DITHER is generally a very bad idea.'
	    dithmsg=.false.
	  endif
	endif
	ffmjd=ffmjd+clk2
	fmjd=fmjd+clk2
	go to 60  ! skip over end-of-file bits

c   Get here on end-of-file
40	if(lu.gt.50)then
	   close(lu)
	   lu=lu-1
	   go to 10
	endif

c   Get here on end-of-all-files
45	continue
        if(.not.gro) go to 100
	last=.true.
	nsite=0
	rfrq=0.
	t0geo=nint((amjd1+amjd2)/2.)
	fmjd=t0geo
	nfmjd=fmjd
        ffmjd=0.d0

c   Back to processing of all TOAs
 60	yrs=(fmjd-pepoch)/365.25d0

	if(nsite.ge.1 .and. .not.siteused(nsite)) then
              write(*,1063) n,nsite,nfmt,card
1063	      format(i5,' *** Site',i3,' is not defined.'
     +                ' read in. ***    nfmt:',i2/1x,a80)
              stop
        endif
	if(nsite.gt.0)then
	  site(1)=hrd(nsite)*dcos(hlt(nsite))*499.004786d0
	  site(2)=site(1)*dtan(hlt(nsite))
	  site(3)=alng(nsite)
	else
	  site(1)=0.d0
	  site(2)=0.d0
	  site(3)=0.d0
	endif

	freqhz=rfrq*1.0d6	!NB: rfrq=0 interpreted as infinite frequency
	ddmch(n)=ddm		!Save ddm
	if(nddm.eq.0) ddm=0	!Disable DM corrections if not wanted
	dmtot=dm+ddm
        if(ndmcalc.ge.2 .and. fmjd.ge.dmvar1 .and. fmjd.le.dmvar2) then   
!        if(ndmcalc.ge.2) then
          fac = 1
	  do 61 i=1,ndmcalc-1
C IHS June 3 2011: use dmyrs now for separate dmepoch
 	    dmyrs=(fmjd-dmep)/365.25d0
            fac = fac * dmyrs / real(i)
	    dmtot=dmtot + fac*dmcof(i) 
61          continue
	endif
        if (usedmx) then
c         first see whether this point fits into an existing range
          do i = 1, ndmx
            if (nfmjd+ffmjd.ge.dmxr1(i)
     +           .and.nfmjd+ffmjd.le.dmxr2(i)) then
              if(dmxep(i).lt.10) dmxep(i)=(dmxr1(i)+dmxr2(i))/2.0
              idmx = i
              goto 80
            endif
            if(dmxep(idmx).lt.10) dmxep(idmx)=
     +           (dmxr1(idmx)+dmxr2(idmx))/2.0
          end do
c         next see whether a range can be expanded to fit it in
          do i = 1, ndmx
            if (nfmjd+ffmjd.ge.min(dmxr1(i),dmxr2(i)-dmxt)
     +           .and.nfmjd+ffmjd.le.max(dmxr2(i),dmxr1(i)+dmxt)) then
              dmxr1(i)=min(dmxr1(i),nfmjd+ffmjd)
              dmxr2(i)=max(dmxr2(i),nfmjd+ffmjd)
              if(dmxep(i).lt.10) dmxep(i)=(dmxr1(i)+dmxr2(i))/2.0
              idmx = i
              goto 80
            endif
          end do
          if (nonewdmx) then ! give this point zero weight, don't create new range
            wgt = 0
            if (jits.eq.0) then 
              nz = nz + 1
              zerowt(nz) = n
            endif
            idmx = 0
            goto 85
          endif
          if (ndmx+1.ge.NDMXMAX) then
            print *,"Error at TOA number ",n
            print *,"Too many dm offsets.  Maximum: NDMXMAX=",NDMXMAX
            print *,"  To change this, edit dim.h & recompile tempo"
            stop
          endif 
c         it doesn't fit into any existing range, so create a new range
          ndmx = ndmx + 1
          dmxr1(ndmx) = nfmjd+ffmjd
          dmxr2(ndmx) = nfmjd+ffmjd
          idmx = ndmx	
	  if (ndmx.gt.1 .or. firstdmx) then ! Skip 1st DMX if needed
            nfit(NPAR6+2*ndmx-1)=1
            nparam=nparam+1
            mfit(nparam)=NPAR6+2*ndmx-1
	    if (usedmx1.gt.0) then ! DMX and DMX1
              nfit(NPAR6+2*ndmx)=1
              nparam=nparam+1
              mfit(nparam)=NPAR6+2*ndmx
	    else ! Just DMX
              nfit(NPAR6+2*ndmx)=0
            endif
          endif
 80       continue
        endif

 85     continue
C IHS Does not enforce continuity of DM at DMX boundaries.
        if (usedmx .and. idmx>0) then
C                write(*,*) 'Adding at ',fmjd,dmx(idmx) + dmx1(idmx)*
C     +          ((nfmjd+ffmjd-dmxep(idmx))/365.25)
          dmtot = dmtot + dmx(idmx) + dmx1(idmx)*
     +          ((nfmjd+ffmjd-dmxep(idmx))/365.25)
        endif

        if (dmxnout) then
          if (usedmx .and. idmx>0) then
             write (40,*) idmx
          else
             write (40,*) 0
          endif
        endif

	bval=dmtot/2.41d-16
	  
        wflag = 1
        !   the nfmjd.gt.0 condition in the following prevents ztim
        !   from being called if 'tz' mode is set but no reference
        !   TOA has been specified, which sets nfmjd to default value
        !   of zero and nsite to 0.  This is potentially problematic
        !   because it is not possible to do a barycenter correction
        !   for MJD=0 since it is before the time of the ephemeris.
        !   (in fact, ztim just reinitializes itself if nfmjd=0, so
        !   in principle it could be run anyway, but that is ugly.)
        if (nsite.ge.0.and.nfmjd.gt.0) then
          call ztim(nfmjd,ffmjd,nct,fct,wflag)
        elseif (freqhz.lt.1) then ! barycenter, infinite frequency
          nct = nfmjd
          fct = ffmjd
          frq = freqhz * 1.e-6
        else                    ! barycenter, non-infinite frequency
          nct = nfmjd
          fct = ffmjd - bval/freqhz**2/SECDAY
          if (fct.lt.0.) then
            fct = fct + 1
            nct = nct - 1
          elseif (fct.ge.1.) then
            fct = fct - 1
            nct = nct + 1
          endif
          frq = freqhz * 1.e-6
        endif

	if(search.and.ct00.gt.0.d0.and.(abs(nct+fct-ct00).gt.trkmax)) 
     +       then
	  nblk=nblk+1
	  call pcard(card2,mode,zawgt,deltat,fmjd,dphase,sigm,offset,
     +     jdcatc,pha1,pha2,efacsave,emin,equadsave,jits,lu,track,trkmax,
     +     search,lw,nfmt,parsed)
	  if(nblk.ge.2) call pcard(card2,mode,zawgt,deltat,fmjd,
     +      dphase,sigm,offset,jdcatc,
     +      pha1,pha2,efacsave,emin,equadsave,jits,lu,track,trkmax,search,
     +      lw,nfmt,parsed)
	endif

	if(jdcatc) xjdoff(1,nxoff)=fmjd
	jdcatc=.false.
c Loop and apply all jumps where x(NPAR2+i).ne.0
        if(jumpbarycenter.and.nxoff.gt.0) then
	  do i=1,nxoff
	    if (x(NPAR2+i).ne.0.d0) fct=fct+dct(i)/86400.d0
	  enddo
        endif
	ct=nct+fct
	cp=p0+p1*(ct-pepoch)*86400.d0
	if (nsite.ge.0) then
	  call earth(cp,x(5),x(6),be,era,edc,pra,pdec,erd)
	endif
	if(ct.gt.ctmax) ctmax=ct
	wgt=1.0
	if(mode.eq.1.and.sigm.eq.0.d0) then
	  terr = sqrt(terr*terr+equad*equad)
	  terr=max(emin,efac*terr)
	  wgt=(1.d-6*terr/cp)**(-2)
	endif

	if(sigm.gt.0.d0) then
            terr=sigm
C  The following (added by JHT in version 7.03) is for PSR 1913+16 at 
C  Arecibo only.  NB: HA is abs(hour angle) in days
            if (zawgt .ne. 0.d0) then
              xlst = mod(tsid/6.283-0.1854+5.d0, 1.d0)
              ha = abs(xlst-0.8021)
              if (ha .gt. 0.0346) terr=sigm*(1.+1.8*(ha-0.0346)/0.0314)
            endif
          wgt=(1.d-6*terr/cp)**(-2)
        endif

        if (jits .eq. 0) then
            if (a1(1).ne.0.d0) then ! zero weight in case of orbital phase cut
              phi = dmod((ct-t0(1))*86400./pb(1) + 9000.,1.d0)
              if ((pha2.gt.pha1 .and. phi.gt.pha1 .and. phi.lt.pha2).or.
     +           (pha2.lt.pha1 .and.(phi.lt.pha2 .or. phi.gt.pha1)))then
                wgt = 0
                nz = nz + 1
                zerowt(nz) = n
              endif
            endif
            if (wflag.eq.0.and.wgt.ne.0) then ! zero weight in case of
                wgt = 0                       ! sun angle (phi) cut
                nz = nz + 1
                zerowt(nz) = n
              endif
          else 
            if (nw .le. nz) then
              if (n. eq. zerowt(nw)) then
                wgt = 0
                nw = nw + 1
              endif
            endif
        endif

	if(last) wgt=1.d-10*wgt

	wt=wgt

	if(xitoa) then				!Write the ITOA file
C  Write itoa file correctly, including observatory code.  (VMK, June94)
	  write(afmjd,1079) ffmjd
1079	  format(f15.13)
	  bsite=siten2b(nsite)

 70	  write(45,1080) psrname,nfmjd,afmjd(2:),terr,rfrq,
     +      ddmch(n),bsite,comment
1080	  format(a9,i5,a14,f6.2,f11.4,f10.6,2x,a2,2x,a8)
	endif

        dphasetot = dphase + dphaseflag
	call resid(nct,fct,dn,dphasetot,dnpls(n),nits,jits)

	if(track.and.n.gt.1) then
	  dt=dt+ntrk
 	  if (abs(ct-ct00).lt.trkmax) then
	    if(abs(dt+1.d0-dt00).lt.abs(dt-dt00)) then
	      dt=dt+1.d0
	      ntrk=ntrk+1
	    else if(abs(dt-1.d0-dt00).lt.abs(dt-dt00)) then
	      dt=dt-1.d0
	      ntrk=ntrk-1
	    endif
	  endif
	endif
	dt00=dt
	ct00=ct

c   write out tracking-corrected pulse number
	if (npulseout.and.jits.eq.0) then
	  write(35,fmt='(f14.0)') dn-ntrk
	endif

c   Write out "INFO" lines to info.tmp.  Write a line if:
c     -- this file has already been opened, in which case
c        there must be a line for every TOA, even if it
c        has zero length
c     -- this file has not been opened yet, but the current
c        TOA has a comment;  in this case, open the file
c        and write out lines of zero length for any previous
c        TOAs 
        if (infolen.ne.0 .or. infoout) then
          if (.not.infoout) then ! need to set up file,
            open(36,file=infofile)
            do i = 1, n-1
              write (36,fmt='()') ! write n-1 blank lines to the file
            end do        
            infoout = .true.
          endif
          write (36,fmt='(a)') infotxt(1:infolen)
        endif


C  Take a shortcut when called by TZ:
	if(tz) then
	  if(n.eq.2) ct2=ct
	  dnpls(n)=dt+dn
	  go to 10
	endif

C  DM-related partial derivatives

	x(16)=0.d0
        if (usedmx) then
          do i = 1, 2*NDMXMAX
            x(NPAR6+i) = 0
          end do
        endif

	if(frq.gt.1.d0) then
          x(16)=f0*1.0d4/(2.41d0*frq**2)
          if (usedmx.and.idmx>0) then
             x(NPAR6+2*idmx-1) = x(16)
             x(NPAR6+2*idmx) = x(16)*((nfmjd+ffmjd-dmxep(idmx))/365.25)
C Figure out min/max freq used in each DMX segment
	     if (dmxf1(idmx).eq.0) dmxf1(idmx)=rfrq
	     if (dmxf2(idmx).eq.0) dmxf2(idmx)=rfrq
	     dmxf1(idmx) = min(dmxf1(idmx),rfrq)
	     dmxf2(idmx) = max(dmxf2(idmx),rfrq)
          endif
        endif

C	if(nfit(16).ge.2) then
C IHS based on Jan 2009: change to allow dmpoly in one section
        if(nfit(16).ge.2 .and. fmjd.ge.dmvar1 .and. fmjd.le.dmvar2) then
          fac = 1.
	  do 89 i=1,nfit(16)-1
C IHS June 3 2011: use dmyrs now for separate dmepoch
            fac = fac * dmyrs / real(i)
	    x(NPAR9+i)= fac * x(16) 
89        continue
	endif

C Non-DM frequency dependent shifts (FDn)
C TODO allow arb reference freq instead of 1 GHz?
	do 898 i=1,NFDMAX
	  x(NPAR10+i) = 0.d0
898     continue
	if (frq.gt.1.d0) then
	  fac = log(frq * 1.e-3)
	  do 899 i=1,NFDMAX
	    x(NPAR10+i) = f0*fac**i
899       continue
	endif

C  XM-related partial derivatives

        do i = 1, nxmx
          if (         xmxuse(i) 
     +      .and. (xmxf1(i).lt.0.d0 .or. frq.ge.xmxf1(i))
     +      .and. (xmxf2(i).lt.0.d0 .or. frq.le.xmxf2(i))
     +      .and. (xmxr1(i).lt.0.d0 .or. nmjd+fmjd.ge.xmxr1(i))
     +      .and. (xmxr2(i).lt.0.d0 .or. nmjd+fmjd.lt.xmxr2(i)) ) then
            x(NPAR12+2*i-1) = f0 * (frq/xmxfrq0)**xmxexp(i)
            x(NPAR12+2*i  ) = x(NPAR12+2*i-1)* xmx(i) *  
     +                                          log(frq/xmxfrq0)
          else
            x(NPAR12+2*i-1) = 0.d0
            x(NPAR12+2*i  ) = 0.d0
          endif
        enddo


C Save the TOA flags for use elsewhere
C The NPTSDEF check could be removed if this array were to be 
C dynamically allocated.
	if (n.lt.NPTSDEF) stflags(n) = rawflags

C Save the DM "residual" and error (could make this part of vmemrw stuff?)
	if (usedmdata) then
          if (dmobserr.gt.0.d0) then
            dmres(n) = dmobs + dmjump - dmtot
            dmerr(n) = dmefac * sqrt(dmobserr**2 + dmequad**2)
          else
            dmres(n) = 0.0
            dmerr(n) = 0.0
          endif
	endif

	x(17)=f0*dtdpx
	x(36)=f0*dtdpmrv
	x(19)=f0*dtdppng
	do 90 j=1,nparam-1
90	fctn(j)=x(mfit(j+1))
	do 92 j=nparam,NPAP1
92	fctn(j)=0.
	do 93 j=1,nparam-1
93	xmean(j)=xmean(j)+wgt*fctn(j)
        sum=sum+wgt*dt
	sumsq=sumsq+wgt*dt*dt
	sumwt=sumwt+wgt
	call vmemw(n,fctn,ct,dt,wgt,dn,terr,frq,fmjd,rfrq,nparam-1,
     +     buf,npmsav,ksav,nbuf,memerr)

	nitsz=max0(nits,1)
	if(mod(n,20*modscrn).eq.1.and..not.quiet) 
     +        write(*,1099) max0(nitsz,1)
1099	format(/'    N       MJD       Residual (p)  Residual (us)',
     +    '  Iteration (of',i2,')'/)
	if(n.gt.200) modscrn=100
	amjd1=min(amjd1,fmjd)
	amjd2=max(amjd2,fmjd)
	if(mod(n,modscrn).eq.1.and..not.quiet) 
     +        write(*,1100)n,fmjd,dt,1d6*dt*p0,
     +    jits+1
1100	format(i7,f15.8,f11.6,f15.3,11x,i2)

        fmjdlast = fmjd

	if(.not.last) go to 10

c End of input file detected
100	continue

        if(ntzref.eq.0) ntzref=n

        if(mod(n,modscrn).ne.1.and..not.quiet) 
     +    write(*,1100)n,fmjdlast,dt,1d6*dt*p0,jits+1

	if(.not.usestart) start=amjd1-1.d-3
	if(.not.usefinish) finish=amjd2+1.d-3

	fitmode = mode ! store final mode value

	if(.not.tz) then
	  sigma1=sqrt((sumsq-sum*sum/sumwt)/sumwt)*p0*1000.d0
	  m=-1
	  call vmemw(m,fctn,ct,dt,wgt,dn,terr,frq,fmjd,rfrq,nparam-1,
     +     buf,npmsav,ksav,nbuf,memerr)

	endif

        if (infoout) then
          close (36)
          infoout = .false.
        endif
9990	continue
	return
	end

