c      $Id$
      subroutine propmo(fa,fd,ea,ed,ara,adc)

      implicit none

      real*8 twopi
      parameter (twopi=6.28318530717958648d0)

      real*8 fa, fd, ea, ed, ara, adc

      real*8 pm, phi, pmerr, phierr
      real*8 c

      real*8 q, galang

      logical err


c  pm  (pmerr) is space proper motion
c  phi (phierr) is position angle of space vector
c  galang is position angle in galactic coordinate

      pm=dsqrt(fd**2+fa**2)
      phi=datan2(fa,fd)
      if (phi.lt.0.) phi=phi+twopi
      phi=phi*180./3.1415926
      call covar(7,8,c,err)
      if (.not. err) then
        pmerr = sqrt(fa*fa*ea*ea+fd*fd*ed*ed+2*c*fa*fd)/pm
        phierr = sqrt(ea*ea*fd*fd + fa*fa*ed*ed + 2*c*fa*fd) / pm**2
        phierr = phierr *180./3.1415926
        call galco (ara,adc,q,phi,galang)
        write(31,10) pm,pmerr,phi,phierr,galang
 10     format(/'Composite proper motion =',f9.4,' +/-',f9.4,' mas/yr'/
     +       'Position angle:  Celestial:',F8.2,' +/-',F7.2,' deg',
     +       '  Galactic:',F8.2,' deg')
      else
        write(31,20) pm, phi, galang
 20     format(/'Composite proper motion =',f9.4,' mas/yr'/
     +       'Position angle:  Celestial:',F8.2,' deg',
     +       '  Galactic:',F8.2,' deg')
      endif


      return
      end
