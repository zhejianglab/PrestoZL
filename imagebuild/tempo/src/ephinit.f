c      $Id$
      subroutine ephinit (iunitin, ephfile, ephbigendian)

c     set up a JPL ephemeris file for reading

c     --open the file with a bogus record length
c     --read some header info
c     --figure out the real record length
c     --close the file
c     --re-open it using the correct record length

c     based on JPL routine "fsizer2.f"

c     DJN 29 July 1997

      include 'ephcom.h'

c     inputs:
      integer iunitin           ! fortran unit number for the file
      character*256 ephfile     ! ephemeris file name
      logical ephbigendian
      logical bigendian
      integer nrecl
      logical openflag


c     stuff read from the data file header;  not all is used, but I've
c     left in the JPL variable names for possible future use;  some of
c     these are defined in ephcom.h, but I've listed them here as well
c     for clarity


      character*6 ttl(14,3)
      character*6 cnam(400)     ! names of ephemeris constants
c      real*8 ss(3)              ! start date, end date, dates per block
      integer ncon              ! number of ephemeris constants
      real*8 au                 ! value of 1 AU in km
c      real*8 emrat              ! ratio of earth mass to moon mass
c      integer ipt(3,13)         ! indecies within records of 13 data sets
      integer numde            

c     other variables
      integer i, j
      integer khi, kmx          ! mysterious names used in JPL software
      integer nd, nfl

c     the definition of "recl" in file open statements in system-dependent
c     nrecl=4 if recl is in bytes (sun compiler; gcc compiler)
c     nrecl=1 if recl is in 4-byte words (vax compiler; intel compiler)
c     the following (fortran90) should work on all modern systems
      real*4 tmp
      INQUIRE (IOLENGTH=nrecl) tmp

      clight = 299792.458d0     ! defined to be constant

      nrl = -1                  ! initialize last-record-read variable

      iunit = iunitin
      nfl = index(ephfile,' ')-1
c     

      inquire(unit=iunit,opened=openflag)
      if (openflag) close(iunit)

      open (iunit, file=ephfile(1:nfl), access='DIRECT', 
     +     form='UNFORMATTED',recl=nrecl*1000, status='OLD',err=999)

      read (iunit, rec=1) ttl, cnam, ss, ncon, au, emrat,
     +     ((ipt(i,j),i=1,3),j=1,12), numde, (ipt(i,13),i=1,3)

      close (iunit)
      if(ephbigendian.and..not.bigendian() .or.
     +  .not.ephbigendian.and.bigendian()) then
        call dbyterev(ss,3)
        call byterev(ncon,1)
        call dbyterev(au,1)
        call dbyterev(emrat,1)
        call byterev(ipt,39)
        call byterev(numde,1)
      endif

c     ipt(1,i) gives location within a record of body i, so finding
c     the highest-valued ipt(1,i) gives the location of the start of
c     data for the last body;  ipt(2,i)*ipt(3,i) is the number of blocks
c     times the number of coefficients for that body, so for ND dimensions
c     [ND=3 for planets, ND=2 for nutations], the number of real*8 words
c     in a record is ipt(1,i) - 1 + 3*ipt(2,i) 

      khi = 0
      kmx = 0
      do i = 1, 13
        if (ipt(1,i) .gt.kmx) then
          kmx = ipt(1,i)
          khi = i
        endif
      enddo

      nd = 3
      if (khi.eq.12) nd = 2

      ncoeff = (ipt(1,khi)-1+nd*ipt(2,khi)*ipt(3,khi))

      open (iunit, file=ephfile(1:nfl), access='DIRECT', 
     +     form='UNFORMATTED', recl=2*nrecl*ncoeff, status='OLD')

      return

 999  write(*,'(''Failed to open ephemeris file: '',a)')ephfile(1:nfl)
      STOP

      end
