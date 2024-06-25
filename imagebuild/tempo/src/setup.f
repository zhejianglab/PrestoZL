c       $Id$
        subroutine setup(version,infile,obsyfile,alng,parfile)

	implicit real*8 (a-h,o-z)
	parameter (TWOPI=6.28318530717958648d0)
	include 'dim.h'
	include 'acom.h'
        include 'version.h'
	real*8 alng(NOBSMAX)
	integer time
	character timstr*24,obsnam*12,damoyr*9,parfile*160
	character*640 infile
        character*(*) obsyfile
        character*640 tpocmd
        character*5 obskey0
	data ault/499.004786d0/

        character*1 siten2a      ! external function
        character*2 siten2b      ! external function
        integer sitea2n          ! external function
        character*640 getcmdline ! external function

        integer nobs             ! if non-single-character obs codes are used, this
                                 !   is set to the maximum value used; otherwise it
                                 !   is 35 (max of 1-9 and a-z).

        nobs = 35

	nsec=time()
	timstr=damoyr(40587+nsec/86400)
	nsec=mod(nsec,86400)
	nhr=nsec/3600
	nmin=(nsec-3600*nhr)/60
	nsec=mod(nsec,60)
	write(timstr(10:),1000) nhr,nmin,nsec
 1000	format(2x,i2.2,':',i2.2,':',i2.2,' UTC')

	if (.not.quiet)
     +       write(31,30) version,
     +       VERSIONID(1:10),VERSIONID(12:18),
     +       timstr(1:23)
30	format(55(1h*)/1h*,53x,1h*/
     +    '*',17x,'TEMPO version',f7.3,16x,1h*/
     +    1h*,11x,'Analysis of pulsar timing data',12x,1h*/
     +    1h*,9x'TEMPO git version date:  'a10,9x,1h*/,
     +    1h*,9x'TEMPO git version id:    'a7,12x,1h*/,
     +    1h*,9x,'Run time: 'a24,10x,1h*/
     +    1h*,53x,1h*/55(1h*)///
     +    20x,'Observatory Coordinates'//
     +    5x,'Observatory     Geodetic lat   Geodetic long   Elevation'/
     +    22x,'ddmmss.ss       ddmmss.ss         m'/)

C  Read the geodetic coordinates of up to NOBXMAX observatories.
C  Icoord = 0 if geodetic; 1 = X, Y, Z geocentric coordinates, where
C  alat = X, along = Y, and elev = Z
C
C  numerical code         one-character code    one-character code
C  internal to tempo         in obsys.dat         in TOA files
C  and used in polyco.dat                  
C     1 .. 9                    1 .. 9              1 .. 9
C    10 .. 35                   a .. z              a .. z
C    36 .. NOBSMAX                -                  (none)

	open(2,file=obsyfile,status='old',err=32)
	goto 34
 32	write(*,'(''Error opening Observatory file: '',a)')obsyfile
        stop
 34     continue
        do j = 1, NOBSMAX
          siteused(j) = .false.
        enddo
   	do 10 j=1,NOBSMAX
	  read(2,40,end=11) alat,along,elev,icoord,obsnam,obskey0
          if (obskey0(1:1).eq.' ') then
            jsite = j
          else if (obskey0(1:1).eq.'-') then
            nobs = nobs + 1
            jsite = nobs
          else
            jsite = sitea2n(obskey0)  ! only first char of obskey0 is used by sitea2n
          endif
          if (siteused(jsite)) then
            print *,"Error, site ",jsite, "listed more than once in"
            print *,obsyfile
            stop
          endif
          siteused(jsite) = .true.
          obskey(jsite) = obskey0
 40	  format(3f15.0,2x,i1,2x,a12,8x,a5)
	  if(alat.ne.0..and..not.quiet) 
     +         write(31,50) siten2a(jsite),siten2b(jsite),
     +            obsnam,alat,along,elev
 50	  format(a1,x,a2,x,a12,f15.2,f16.2,f12.1)
 
	  if(icoord.eq.0)then
	    alat=ang(1,alat)
	    alng(jsite)=ang(1,along)
           
c  old approach is IAU 1964 power series.  Tests show this is
c  good to 5.10^-10 for 1964 figure.  But series misapplied
c  below (erad series should be in alat, not hlt)
c  so actual error is 1.d-5.  Furthermore, want flexibility
c  of going to any figure to use full equation approach.
c  see AA K11.
c  IAU 1976 flattening f, equatorial radius a

	    aa_f = 1.d0/298.257d0
	    aa_a = 6378140.d0
	    aa_c = 1.d0/sqrt(1.d0+(-2.d0+aa_f)*aa_f*sin(alat)*sin(alat))
	    aa_arcf = (aa_a*aa_c+elev)*cos(alat)
	    aa_arsf = (aa_a*(1.d0-aa_f)*(1.d0-aa_f)*aa_c+elev)*sin(alat)
	    hlt(jsite) = datan2(aa_arsf,aa_arcf)
            erad = sqrt(aa_arcf*aa_arcf+aa_arsf*aa_arsf)

c	    delat=-692.7430d0*dsin(2.d0*alat)+1.1633d0*dsin(4.d0*alat)-
c     +        0.0026d0*dsin(6.d0*alat)
c	    hlt(j)=alat+delat*TWOPI/12.96d5
c	    erad=6378160.0d0*(0.998327073d0+0.001676438d0*
c     +        dcos(2.d0*hlt(jsite))-3.519d-6*dcos(4.d0*hlt(jsite))+
c     +        8.d-9*dcos(6.d0*hlt(jsite)))+elev

	    hrd(jsite)=erad/(2.99792458d8*ault)
	  else
	    erad=dsqrt(alat**2+along**2+elev**2)
	    hlt(jsite)=dasin(elev/erad)
	    alng(jsite)=datan2(-along,alat)
	    hrd(jsite)=erad/(2.99792458d8*ault)
	  endif
 10	continue
        print *, "Error: too many observatory codes in obsys.dat"
        print *, "Please increase NOBSMAX in dim.h and recompile tempo"
        stop

 11     continue
	close(2)

        if (.not.quiet) then
          len=index(infile,' ')-1
          if(oldpar.or.parunit.eq.50)then
            write (31,60) infile(1:len)
 60         format (/'Input data from ',a)
          else
            len1=index(parfile,' ')-1
            write (31,61) infile(1:len),parfile(1:len1)
 61         format (/'Input data from ',a,',  Parameters from ',a)
          endif
          write(31,62)obsyfile
 62       format('Observatory data from ',a55)
          tpocmd = getcmdline()
	  write(31,65) tpocmd(1:len_trim(tpocmd))
 65       format ('Tempo command: ',a)
        endif

	return
	end

