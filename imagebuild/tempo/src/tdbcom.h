c      $Id$
c     common block for tdbinit.f and tdbread.f 

      real*8 tdbd1, tdbd2       ! start, end JD of the file
      integer tdbdt             ! number of days per block of coefficients
      integer tdbncf            ! number of coefficients per block
      integer tdbrecsize        ! number of coefficients per block of all
                                !   types (time and velocity records for
                                !   irwin-type files)

      integer tdbunit           ! fortran i/o unit number
      integer tdbnrl            ! last record read
      integer ipt(3,2)

      common /tdbcom/ tdbd1, tdbd2, tdbdt, tdbncf, tdbunit, tdbnrl,
     +                  tdbrecsize, ipt

