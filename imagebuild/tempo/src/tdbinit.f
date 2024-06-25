c      $Id$
      subroutine tdbinit (tdbunitin, tdbfile)

c     set up the TDB-TDT ephemeris file for reading
c     modelled after the JPL ephemeris file init/read routines

c     --open the file with a bogus record length
c     --read some header info
c     --figure out the real record length
c     --close the file
c     --re-open it using the correct record length

c     DJN 19 August 1997

      include 'tdbcom.h'

c     inputs:
      integer tdbunitin         ! fortran unit number for the file
      character*640 tdbfile     ! ephemeris file name

      integer nfl

      integer nrecl

      logical bigendian         ! external routine

c     the definition of "recl" in file open statements in system-dependent
c     nrecl=4 if recl is in bytes (sun compiler; gcc compiler)
c     nrecl=1 if recl is in 4-byte words (vax compiler; intel compiler)
c     the following (fortran90) should work on all modern systems
      real*4 tmp
      INQUIRE (IOLENGTH=nrecl) tmp

      tdbnrl = -1               ! initialize last-record-read variable

      tdbunit = tdbunitin
      nfl = index(tdbfile,' ')-1

c     open file with a bogus record length, get header info, close file

      open (tdbunit, file=tdbfile(1:nfl), access='DIRECT', 
     +     form='UNFORMATTED',recl=6*nrecl, status='OLD',err=999)

      read (tdbunit, rec=1) tdbd1, tdbd2, tdbdt, tdbncf

      close (tdbunit)
      if(.not.bigendian()) then
        call dbyterev(tdbd1,1)
        call dbyterev(tdbd2,1)
        call byterev(tdbdt,1)
        call byterev(tdbncf,1)
      endif
      ! following defined for compatability with irwin-type files
      tdbrecsize = tdbncf      
      ipt(1,1) = 1
      ipt(2,1) = tdbncf
      ipt(3,1) = 1

c     now re-open file with correct record length

      open (tdbunit, file=tdbfile(1:nfl), access='DIRECT', 
     +     form='UNFORMATTED', recl=2*nrecl*tdbncf, status='OLD')

      return

 999  write(*,'(''Error opening TDB file: '',a)')tdbfile(1:nfl)
      STOP

      end
