      character*2 function siten2b(nsite)

C     convert integer n into an ASCII observatory code

      implicit real*8 (a-h,o-z)
      include 'dim.h'
      include 'acom.h'

      integer nsite

      siten2b = "--"   ! default for no assigned two-character code

      if (nsite.eq.-1) then
        siten2b = '@@'
      else if (nsite.eq.0) then
        siten2b = '00'
      else
        siten2b = obskey(nsite)(4:5)
      endif

      return
      end
