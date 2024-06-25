c      $Id$
      subroutine tdbinit2(tdbunitin, tdbfile)

c     set up the TDB-TDT ephemeris file for reading
c     modelled after the JPL ephemeris file init/read routines

c     --open the file with a bogus record length
c     --read some header info
c     --figure out the real record length
c     --close the file
c     --re-open it using the correct record length

c     This version of the code is adapted to use Alan Irwin's
c     tdb ephemeris files, which hopefully have the a very
c     similar in format to the JPL ephemeris files.
c     It makes use of Irwin's code time_fsizer2.f, which
c     presumably is based on the JPL routine "fsizer.f"

c     DJN 30 July 2004; based on tdbinit

      include 'tdbcom.h'

c     inputs:
      integer tdbunitin         ! fortran unit number for the file
      character*640 tdbfile     ! ephemeris file name

      integer nfl

      integer nrecl

      logical bigendian         ! external routine

c     stuff read from the data file header;  not all is used, but I've
c     left in the variable names for possible future use;  


      character*6 ttl(14,3)
      character*6 cnam(2)       ! names of ephemeris constants
      integer ncon              ! number of ephemeris constants
c     note: tdbt1, tdbt2 are defined real*8 in tdbcom.h
c     note: tdbdt is defined int in tdbcom.h
c     note: ipt is defined int(3,2) in tdbcom.h
      real*8 tmptdbdt

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
     +     form='UNFORMATTED',recl=79*nrecl, status='OLD',err=999)

      read (tdbunit, rec=1) ttl, cnam, tdbd1, tdbd2, tmptdbdt, ncon, ipt

      close (tdbunit)
      if(bigendian()) then  ! irwin ephem make on linux/x86 (little endian)
        call dbyterev(tdbd1,1)
        call dbyterev(tdbd2,1)
        call dbyterev(tmptdbdt,1)
        call byterev(ncon,1)
        call byterev(ipt,6)
      endif
      tdbdt = int(tmptdbdt)

      tdbncf = ipt(2,1)   ! number of coefficients
      tdbrecsize = (ipt(1,2)-1 + 3*ipt(2,2)*ipt(3,2)) !'ksize/2' in irwin's code

c     now re-open file with correct record length

      open (tdbunit, file=tdbfile(1:nfl), access='DIRECT', 
     +     form='UNFORMATTED', recl=2*nrecl*tdbrecsize, status='OLD')

      return

 999  write(*,'(''Error opening TDB file: '',a)')tdbfile(1:nfl)
      STOP

      end
