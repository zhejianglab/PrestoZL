c      $Id$
      SUBROUTINE ZTIM(nmjdu,fmjdu,ntdb,ftdb,wflag)

C     DERIVED FROM PROGRAM BAT ('BARYCENTRIC ARRIVAL TIME')
C ORIGINAL VERSION 1970 BY DAVID RICHARDS, BASED ON CODING IN MIT'S
C    PLANETARY EPHEMERIS PROGRAM
C REVISIONS BY J. F. CHANDLER
C    1983 GET MICROSEC. ACCURACY
C    1984 GET 0.1 MICROSEC. ACCURACY AND CONVERT TO SUBROUTINE
C FURTHER REVISIONS:
C    17 MAY 1984 JMW: ENTER AND LEAVE ROUTINE WITH TRUE JULIAN DAYS:
C	JDU=JDUTC+FRUTC+0.5
C       FRU=DMOD(FRUTC+0.5D0,1.D0)
C       JDET+FRJDET=JDC+FRC-0.5
C    17 MAY 1984 JMW: TRANSFER BARYCENTRIC OBS FREQUENCY FRQ AND
C	S.S. BC ->EARTH COORDS BE(3) TO TEMPO
C    29 MAY 1984 JMW: CALCULATE RCE(4-6)--VELOCITY OF EARTH W.R.T. SSBC
C    26 JUL 1984 JMW: (1)CHANGE COMMENTS TO SAY THAT SITE(3) IS WEST LONG.
C		      (2)PASS SITVEL THROUGH COMMON OBSP		
C		      (3)INCLUDE OBSERVATORY VELOCITY IN DISPERSION CALC
C			(FIRST POINTED OUT BY LEE FOWLER)
C		      (4)DIVIDE INTERPLANETARY DISPERSION DELAY BY 2.
C Modified for use with JPL ephemeris (D. Nice Summer'88)
C    8 Sept, 1997. MJD input and output (RNM)

C    OTHER INPUT IS PASSED VIA COMMON BLOCKS:
C              THE SOURCE, SITE, AND FREQUENCY ARE SPECIFIED IN
C              COMMON /OBSP/ AND THE OUTPUT LISTING FORTRAN UNIT
C              IS GIVEN IN /INODTA/ (SEE BELOW).
C OUTPUT:
C    ntdb,ftdb - SAME AS INPUT, BUT REPRESENTING THE TDB THE
C              SAME PULSE (AT INFINITE FREQUENCY) WOULD ARRIVE AT THE
C              SOLAR-SYSTEM BARYCENTER IN THE ABSENCE OF THE SOLAR
C              SYSTEM.
C    wflag - set to 0 in timcalc if the weight of this point should be
C            set to zero because the source-sun angle is less than PHIMIN
C    IN ADDITION, THE COORDINATES USED IN THE CALCULATIONS
C    ARE LEFT IN COMMON /CRDS/ (SEE BELOW).

c     Method:  The program reads ephemeris files containing the locations &
c              velocities of the earth-moon barycenter and the sun relative
c              to the solar system barycenter (SSB), and of the earth relative
c              to earth-moon barycenter, as well as angles of earth nutation.
c              The measured UTC arrival time of a pulse is reduced to 
c              Ephemeris time (ET, alternately called Barycentric Dynamical
c              Time, or TDB) of arrival at SSB by converting successively to
c              atomic time (AT) and ET at reception and computing the additional
c              propagation time to SSB, the general relativistic delay, long-
c              period CT-AT terms, and an estimate of interstellar and
c              interplanetary dispersions. 

c              For 0.1 usec accuracy the irregularities of the earth's rotation
c              must be included.  A table for this is supplied, but must be
c              extended continually (e.g., from B.I.H. circular D.)  Following
c              PEP, the AT used in this code is A.1 (=TAI + .0343817), and the
c              coordinate time used for interpolation is CT (=A.1 + 32.15), but
c              the epoch returned is ET (=TAI + 32.184.)
c              UT1 tables are supplied on I/O unit 42.  If 1usec accuracy is
c              adequate, the table may be omitted and replaced by two blank
c              cards

c    *** NOTE: THIS SUBROUTINE MUST BE CALLED ONCE WITH JDU=0   ***
c    ***       TO SET UP BEFORE THE 1ST REAL CALL.  WHEN JDU=0, ***
c    ***       THE OTHER ARGUMENTS ARE IGNORED.                 ***

c This program complies with the Fortran 77 Standard, except that:
c    * Double precision is assigned by "REAL*8"
c    * Double precision variables are assumed to carry at least
c      fifteen digits of accuracy

      implicit real*8(A-H,O-Z)
      save
      integer wflag


c           Common constants

      common /CONST/ PI,TWOPI,SECDAY,CONVD,CONVS,AULTSC,VELC,EMRAT,OBLQ,
     +              GAUSS,RSCHW,AULTVL

c           Coordinates

      common /CRDS/ RCB(6),RBE(6),RCE(6),RSE(3),RCA(3),RSA(3),RBA(3),
     +             REA(3),PSID,EPSD,PC,PS,TSID,PRN(3,3),ATUT,UT1UT,ETAT,
     +             DTGR,TDIS,BCLT      
C R12       - POSITION (+VELOCITY) OF BODY 2 W.R.T. BODY 1 (LS, LS/S)
C             C: SOLAR-SYSTEM BARYCENTER, S: SUN, B: EARTH-MOON BARYCENTER,
C             E: EARTH, M: MOON, A: OBSERVING SITE
C             note:  by definition, R21 = - R12
C S         - POSITION OF SITE W.R.T. EARTH CENTER OF MASS
C PSID,EPSD - NUTATIONS OF LONGITUDE AND OBLIQUITY (RAD)
C PC,PS     - EQUATORIAL AND MERIDIONAL COMPONENTS OF PSID
C PRN       - PRECESSION/NUTATION MATRIX
C ATUT      - VALUE OF A.1 - UTC IN SEC
C UT1UT     - VALUE OF UT1 - UTC IN SEC
C ETAT      - VALUE OF ET - A.1 IN SEC
C DTGR      - RELATIVISTIC TIME DELAY IN SEC
C TDIS      - DISPERSION DELAY IN SEC
C BCLT      - PROPAGATION DELAY FROM SITE TO SOLAR-SYSTEM BARYCENTER

c     observation parameters
        common /OBSP/ SITE(3),POS(3),FREQHZ,BVAL,SITVEL(3)
c         SITE - EARTH-FIXED COORDINATES OF OBSERVING SITE: R,Z IN LIGHT-SEC,
c                AND WEST(NOT EAST--JMW) LONGITUDE IN RADIANS
c         POS  - UNIT VECTOR POINTING AT TRUE SOURCE POSITION IN 1950.0
c                COORDINATES.  NOTE: THE FULL ABERRATION (INCLUDING THE E-TERM)
c                MUST BE REMOVED FROM THE APPARENT SOURCE POSITION.  OTHERWISE
c                ERRORS UP TO A MILLISECOND MAY RESULT IN THE OUTPUT TIME
c         FREQHZ - RECEIVED FREQUENCY IN HERTZ
c         BVAL - DISPERSION MEASURE SUCH THAT DELAY IN SEC. = BVAL/FREQHZ**2

c     parameters to be transferred to arrtim
        include 'trnsfr.h'

      real*8 JD(2),NUT(4),RCS(6)

C     initialize quantities
        if(nmjdu.eq.0)then
           OB=OBLQ*CONVD
           RSCHW=GAUSS**2*AULTSC**3/SECDAY**2
           AULTVL=AULTSC/SECDAY
           CTATC=32.15D0
           ETATC=32.184D0-0.0343817D0
           if(ut1flag)CALL UT1RED(0,DUM,DUM,DUM)
           RETURN
        endif

c     get approx. CT epoch in sec since 1950 Jan 1
        ATUT=A1UTCF(nmjdu,fmjdu)
        nmjdc=nmjdu
        fmjdc=fmjdu + (ctatc+atut)/SECDAY
        if(fmjdc.ge.1.d0)then
           nmjdc=nmjdc+1
           fmjdc=fmjdc-1.d0
        endif

c     real JED for JPL ephemeris and BdL CTATVs (Ephemeris Time!!!)
        nmjde = nmjdc
        jd(2) = fmjdc + (ETATC-CTATC)/SECDAY +0.5d0
        if(jd(2).ge.1.d0)then
           nmjde=nmjde+1
           jd(2)=jd(2)-1.d0
        endif
        jd(1) = nmjde + 2400000
c     get TDB - TDT
        call tdbread(jd(1),jd(2),ctatv)

c     adjust jd to TDB for ephemeris read
        jd(2) = jd(2) + ctatv/secday
c     call EPHREAD to get interpolated ephemeris data
        call ephread(JD,RCB,RBE,RCS,NUT)

c	mean obliquity from AA (AWI change).
c	jd is TDB.
	toblq=((jd(1)-2451545.0d0) + jd(2))/36525.d0
c	mean obliquity in deg AA 1984 S21.
	oblq = (((1.813d-3*toblq-5.9d-4)*toblq-4.6815d1)*toblq +
     +    84381.448d0)/3600.d0
          OB=OBLQ*CONVD

c     Calculate distances not found explicitly in the ephemeris
        do 80 I = 1, 3
   80     RSE(I)  = -RCS(I)+RCB(I)+RBE(I)
        do 81 I = 1, 6
   81     RCE(I)  =  RCB(I)+RBE(I)

c     Get variations from calculated Precession & Nutation
        PSID = NUT(1)
        EPSD = NUT(2)
c	this first one is important equation of equinoxes.
c	this agrees with tables, but not completely clear that AA tables are correct.
        PC = cos(ob+epsd)*PSID
        PS = sin(ob+epsd)*PSID

c     Compute precession-nutation matrix
        fmjd = dfloat(nmjde) + jd(2) - 0.5d0
        call PRCNUT(fmjd, psid, epsd, ob ,prn)

c     Do the rest of the arrival time calculation in subroutine TIMCALC
c     Enter with approx CT, return with TDB.
        call TIMCALC(nmjdu,fmjdu,nmjdc,fmjdc,ctatv,etatc,wflag)
        ntdb=nmjdc
        ftdb=fmjdc

        return
        
        end
