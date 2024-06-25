c      $Id$

	logical gro,sim,xitoa,oldpar,psrframe,jumpout,eclcoord
	logical usestart, usefinish
        logical usedmx, firstdmx, nonewdmx
        logical npulsein, npulseout, infoout, quiet, polystdout
        logical phisunout, dopplerout
        logical dmxnout
        logical tz
        logical autotz
        logical useannorb   ! flag to use annual-orbital parallax
        logical k96         ! flag to use full Kopeikin 95/96 terms - Willem
        logical usefixeddist! flag to use fixed input distance in ann-orb px
        logical jumpbarycenter ! flag to apply jumps in barycenter instead of
                               ! obs frame.  Jumps used to be done in
                               ! bary frame, but now (10-jul-04) default
                               ! is obs frame
        logical tdbif99fmt  ! true if tdb file is in format used
                            ! by alan irwin (reference if99)
        logical siteused      ! flags for which observing codes are used.
        logical nofitjump     ! blocks fitting for a given JUMP
        logical useglsfit     ! use generalized-least-squares fit
        logical usedmdata     ! use input DM measurements as "data"
	logical usexmxfrq0    ! non-default xmxfrq0, write it to output

	logical xmxuse        ! whether to use a given xmx range

	integer parunit, nskip, iboot
        integer fitmode
        integer infolen
        integer ssdmflag
        integer nflagjumps  ! number of tempo2-style flag-based jumps
        integer nflagefac   ! number of tempo2-style flag-based EFAC
        integer nflagequad  ! number of tempo2-style flag-based EQUAD
        integer nflagecorr  ! number of flag-based ecorr/jitter terms
        integer nflagdmefac   ! number of tempo2-style flag-based DMEFAC
        integer nflagdmequad  ! number of tempo2-style flag-based DMEQUAD
        integer nflagdmjump   ! number of tempo2-style flag-based DMJUMP
        integer psrkeyl     ! length of PSR key in input/output .par file

	integer nxmx        ! number of xmx terms used


        real*8 phimin
        real*8 PAAscNode    ! position angle of ascending node
        real*8 solarn0      ! solar wind electron density at 1 AU (e-/cc)
        real*8 solarn01     ! solar wind electron density time derivative
        real*8 fixeddist    ! fixed input distance (instead of fit px) 
                            !    for ann-orb px
        real*8 dmepoch, dmep

        real*8 ceecl, seecl ! cosine and sine of obliquity of the ecliptic

	common pdec,pra,ba(3),bc(3),dm,dt,dt2,freq(NPAP1),
     +    ferr(NPAP1),fmin,hlt(NOBSMAX),hrd(NOBSMAX),siteused(NOBSMAX),
     +    wt,x(NPAP1),era,ec,
     +    erd,fmax,emax,tmax,phimin,start,finish,amjd1,amjd2,posepoch,
     +    posep,dither,xjdoff(2,NJUMP),dct(NJUMP),nofitjump(NJUMP),
     +    dmepoch,dmep,
     +    dmx(NDMXMAX),dmxr1(NDMXMAX),dmxr2(NDMXMAX),dmxt,ndmx,usedmx,
     +    dmx1(NDMXMAX),dmxep(NDMXMAX),dmxf1(NDMXMAX),dmxf2(NDMXMAX),
     +    flagefac(NFLAGERR),flagequad(NFLAGERR),flagecorr(NFLAGERR),
     +    flagdmefac(NFLAGERR),flagdmequad(NFLAGERR),flagdmjump(NFLAGERR),
     +    usedmx1,
     +    fdcof(NFDMAX),
     +    rnamp,rnidx,
     +    tcorr,
     +    pmra,pmdec,pmrv,dt2sec,
     +    t0geo,gain,tres,
     +    PAAscNode,solarn0,solarn01,fixeddist,
     +    ceecl, seecl,
     +    nfit(NPAP1),mfit(NPAP1),n,nscan,nparam,nxoff,nprnt,
     +    nkeep,nfq,ncoord,gro,sim,xitoa,oldpar,psrframe,jumpout,
     +	  eclcoord,usestart,usefinish,npulsein,npulseout,
     +    parunit,nskip,iboot,fitmode,ndmcalc,nflagjumps,
     +    nflagefac,nflagequad,nflagecorr,
     +    nflagdmefac,nflagdmequad,nflagdmjump,
     +    psrkeyl,
     +    nfcalc,ntoa,nparam0,ndmx0,infolen,infoout,phisunout,
     +    dopplerout,dmxnout,ssdmflag,
     +    quiet,polystdout,tz,autotz,firstdmx,nonewdmx,
     +    useannorb,usefixeddist,jumpbarycenter,useglsfit,
     +    usedmdata,
     +    tdbif99fmt

        common /xmxcom/ xmx(NXMXMAX), xmxexp(NXMXMAX),
     +                  xmxr1(NXMXMAX), xmxr2(NXMXMAX),
     +                  xmxf1(NXMXMAX), xmxf2(NXMXMAX),
     +                  xmxuse(NXMXMAX),
     +                  xmxfrq0, usexmxfrq0, nxmx

	character psrname*64,obsflag*1,pardir*80,infotxt*160
        character psrkey*32
        character obskey*5
        character jumpflag*32, jumpflagval*32
        character infoflag*32
        character efacflag*32, efacflagval*32
        character equadflag*32, equadflagval*32
        character ecorrflag*32, ecorrflagval*32
        character dmefacflag*32, dmefacflagval*32
        character dmequadflag*32, dmequadflagval*32
        character dmjumpflag*32, dmjumpflagval*32
        character eclcon*80
        character dcovfile*80

        common/acomch/psrname,pardir,obsflag,infotxt,psrkey,
     +    obskey(NOBSMAX),
     +    eclcon,dcovfile,
     +    jumpflag(NJUMP),jumpflagval(NJUMP),infoflag,
     +    efacflag(NFLAGERR),efacflagval(NFLAGERR),
     +    equadflag(NFLAGERR),equadflagval(NFLAGERR),
     +    ecorrflag(NFLAGERR),ecorrflagval(NFLAGERR),
     +    dmefacflag(NFLAGERR),dmefacflagval(NFLAGERR),
     +    dmequadflag(NFLAGERR),dmequadflagval(NFLAGERR),
     +    dmjumpflag(NFLAGERR),dmjumpflagval(NFLAGERR)


        real*8 array(NPA,NPA) ! moved here from fit.f, djn, 8-Sep-98

        common/acomcov/array
