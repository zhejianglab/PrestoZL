c      $Id$
c     ephcom.h

c     DJN 29 Jul 1997

c     ephemris common stuff:

      real*8 ss(3)              ! start, stop, record length in days
      real*8 emrat              ! earth/moon mass ratio
      real*8 clight             ! speed of light
      integer ipt(3,13)         ! indecies of ephemeris records
      integer ncoeff            ! total number of coefs per record
      integer nrl               ! last record number
      integer iunit             ! fortran unit number for data file

      common /ephcom/ ss, emrat, clight, ipt, ncoeff, iunit, nrl

