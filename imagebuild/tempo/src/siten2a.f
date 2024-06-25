c      $Id$
c
      character*1 function siten2a(n)

C     convert integer n into an ASCII observatory code


      implicit none

      integer n

      siten2a = "-"   ! default for no assigned single-character code

      if (n.eq.-1) then
        siten2a = '@'
      else if (n.ge.0.and.n.le.9) then
        siten2a = char(48+n)
      else if (n.ge.10.and.n.le.35) then
        siten2a = char(55+n)
      endif
      return

      end
