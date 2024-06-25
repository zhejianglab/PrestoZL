c      $Id$
c
      subroutine covar(l,m,c,err)

c     Sets 'c' to the element of the covariance matrix array()
c     (generated in fit.f) pointed to by mfit(l+1) and mfit(m+1).  

      implicit real*8 (a-h,o-z)

      integer l
      integer m
      real*8 c
      logical err

      include 'dim.h'
      include 'acom.h'


      c = 0.
      err = .true.

      do i = 1, nparam - 1
        if (mfit(i+1).eq.l) then
          ll = i
          goto 100
        endif
      end do
      goto 990   ! no match found

 100  continue
      do i = 1, nparam-1
        if (mfit(i+1).eq.m) then
          mm = i
          goto 110
        endif
      end do
      goto 990   ! no match found

      ! both indecies match
 110  continue

      err = .false.
      c = array(ll,mm)


 990  continue

      return
      end


      
