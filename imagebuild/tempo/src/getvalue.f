c      $Id$
      character*80 function getvalue(f)

c     return the value of a flag set on a TOA line

      implicit none

      character*(*) f
      integer i

      include 'flagcom.h'

      getvalue = ""
      do i = 1, nflag
        if (f.eq.flagname(i)) then
          getvalue = flagvalue(i)
          exit
        endif
      enddo

      return
      end
