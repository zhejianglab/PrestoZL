      subroutine getecliptic

      implicit real*8 (a-h,o-z)

      include 'dim.h'
      include 'acom.h'
      include 'eph.h'

c     factor to convert arcseconds to radians
c     2*pi/(360*3600)
      real*8 FAC
      parameter (FAC=4.848136811095359935899141d-6)

      integer i, j, k, l, n
      character*640 s, s1, s2
      real*8 ecl

      l = index(eclcon,' ')-1
      call upcase(eclcon)

      k = index(ephdir,' ')-1
      open (2,action='read',status='old',
     +              file=ephdir(1:k)//'ecliptic.dat')

      do while (.TRUE.)
        read (2,fmt='(a640)',end=2001) s

        n = index(s,"#") 
        if (n.gt.0) then
          do i = n, 640
            s(i:i) = " "
          end do
        end if
 
        j = 1
        call citem(s,640,j,s1,k)
        call upcase(s1)
        if (s1(1:l).eq.eclcon(1:l))  then 
          call citem(s,640,j,s2,k)
          if (k.eq.0) then
            write (*,*) "Error in ecliptic.dat, no angle in this line:"
            write (*,*) s
            stop
          endif

          read (s2,*) ecl
          ceecl = cos(FAC*ecl)
          seecl = sin(FAC*ecl)

          goto 9000

       endif
      enddo
 
 2001 continue

      write (*,2011) eclcon(1:l)
 2011 format ("Error: cannot find obliquity of ecliptic for ",A)
      stop


 9000 continue
      return
      end

