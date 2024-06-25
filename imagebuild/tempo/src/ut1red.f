c      $Id$
      SUBROUTINE UT1RED(nmjd,fmjd,ATUT,UT1,quiet)
C***********************************************************************
C         R.KING   OCT 1980    SUBROUTINE UT1RED
C         READ UT1 VALUES FROM AN EXTERNAL DATA SET

C  INPUT:   NMJD,FMJD - INT MJD AND FRACTION FROM 0H ET
C                      IF JD.EQ.0, THEN INITIALIZE ONLY
C           ATUT     - VALUE OF A.1 - UTC AT NMJD+FMJD  (IN SEC)

C  OUTPUT:  UT1      - VALUE OF UT1 - UTC AT NMJD+FMJD  (IN SEC)

      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dim.h'
      logical quiet

      DIMENSION TAB(4),Y1(2),Y2(2),IUT(NUTMAX)
      CHARACTER VARFMT*32,flag*3

      DATA JUT1/42/,nmsg/2/

C         THE EXTERNAL DATA SET (JUT1) MUST HAVE THE FOLLOWING FORM
C     RECORD 1:  TITLE  (10A8)
C     RECORD 2:  INFORMATION DESCRIBING THE TABLE
C         COLUMNS
C          1-32  FORMAT OF DATA ENTRIES  (4A8)
C         33-34  KIND OF DATA IN TABLE  (I2)
C                  =1 UT1-UTC
C                  =2 TAI-UT1
C                  =3  A1-UT1
C         52-53  MAXIMUM NUMBER OF VALUES PER RECORD  (I2)
C         55-56  INTERVAL IN DAYS BETWEEEN TABULAR VALUES  (I2)
C         58-72  UNITS OF VALUES IN SECONDS  (E15.7)
C     RECORDS 3-END:  DATA ENTRIES - FORMAT IS READ IN RECORD 2, BUT
C                     MUST BE OF THE FORM INTEGER MJD FOLOWED BY UP
C                     TO 12 INTEGER VALUES OF UT1.  IF A RECORD IS
C                     SHORT, THE NUMBER OF VALUES IN THAT RECORD IS
C                     GIVEN IN COLUMNS 79-80 (I2).
      save

      IF(nmjd.eq.0)then

C     READ AND WRITE HEADER RECORDS, SET UP ARRAYS
         read(jut1,20)
         READ(JUT1,20) VARFMT,KIND,NPR,INT,UNITS
 20      FORMAT(A32,I2,17X,I2,1X,I2,1X,E15.7)
         
C     READ  DATA RECORDS INTO STORAGE
         nvtot=0
 100     READ(JUT1,VARFMT,END=110) flag, MJD,(IUT(NVTOT+I),I=1,NPR),NVR
         if(flag.eq.'END')go to 110
         
         if(nvtot.eq.0)mjd1=mjd+int
         if(nvr.eq.0)nvr=npr
         nvrold=nvr
         nvtot=nvtot+nvr
         if(nvtot+npr.gt.NUTMAX)then
            print *,'UT1 table overflow'
            stop
         endif
         mjdlast=mjd
         go to 100
         
 110     mjd2=mjdlast+(nvrold-2)*int
         itsav=0
         imsg=0
         return
      endif
C     MJD1 AND MJD2 ARE THE LIMITS OF USABLE VALUES IN ARRAY
C     AND DEPEND ON THE INTERPOLATION SCHEME

C     IS MJD WITHIN RANGE OF TABLE?
      if(nmjd.gt.mjd1 .and. nmjd.lt.mjd2)then

C     CALCULATE INTERPOLATION TIMES AND VALUE OF TABULAR POINTS
         t=((nmjd-mjd1)+fmjd)/int
         it=t
         t=t-it
         s=1.d0-t
         if(it.ne.itsav)then
            do i=1,4
               j=it+i
               tab(i)=iut(j)
            enddo
            do i=1,2
               nr=i+1
               f2=1.d0/6.d0 * (tab(nr+1)+tab(nr-1))
               y1(i)=4.d0/3.d0 * tab(nr) - f2
               y2(i)=-1.d0/3.d0 * tab(nr) + f2
            enddo
            itsav=it
         endif

C        SECOND DIFFERENCE INTERPOLATION
         UT1=(T*(Y1(2)+T*T*Y2(2)) + S*(Y1(1)+S*S*Y2(1)))*UNITS

      else
         imsg=imsg+1
         if(imsg.lt.nmsg.and..not.quiet)then
            write(*,500)nmjd,mjd1,mjd2
            write(31,500)nmjd,mjd1,mjd2
 500        format(' *** Warning - MJD =',i6,
     +            ' outside UT1 table range (',i5,'-',i5,')') 
         else if(imsg.eq.nmsg.and..not.quiet)then
            write(*,510)
            write(31,510)
 510        format(' *** Warning - Further UT1 messages suppressed')     
         endif
C     Extrapolate using last value in table
         if (nmjd.le.mjd1) then
           ut1 = iut(1)*units
         else
           ut1=iut(nvtot)*units   
         endif
      endif

C     CONVERT TABLE VALUES TO UT1-UTC
      if(kind.eq.1)then
         return                 ! Table is already UT1 - UTC
      else if(kind.eq.2)then
         ut1=atut -0.0343817d0 - ut1 ! Table is TAI-UT1
C     UT1-UTC = A1-UTC +  (TAI-A1) - (TAI-UT1)
      else if(kind.eq.3)then
         ut1=atut-ut1           ! Table is A1-UT1
C     UT1-UTC = A1-UTC - (A1-UT1)
      endif

      RETURN

      END



