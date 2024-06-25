c     $Id$
      subroutine mjdadd(nmjd,fmjd,t)

c     add time t (in seconds) to MJD nmjd+fmjd
c     while retaining 0 <= fmjd < 1

      implicit none
      integer nmjd
      real*8 fmjd
      real*8 t

      integer nt
      real*8 ft
  
      nt = int(t/86400.)
      ft = t/86400.-nt

      nmjd = nmjd + nt
      fmjd = fmjd + ft
      if (fmjd.ge.1) then
        nmjd = nmjd + 1
        fmjd = fmjd - 1.
      endif
      if (fmjd.lt.0) then
        nmjd = nmjd - 1
        fmjd = fmjd + 1.
      endif
      return
      end

