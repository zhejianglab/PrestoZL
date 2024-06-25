c      $Id$
      subroutine ephread(JD,RCB,RBE,RCS,NUT)

c     ephread -- read and interpolate JPL planetary ephemeris

c     based on JPL code state.f

c     DJN 29 July 1997

      implicit real*8 (a-h,o-z) ! added 12-Nov-13, needed for eph.h
      include 'dim.h'           ! added 12-Nov-13, needed for eph.h
      include 'eph.h'           ! added 12-Nov-13, needed for bigendian-mess

      include 'ephcom.h' 

      save buf

c     input:
      real*8 JD(2)              ! input JD, make JD(1) int part, JD(2) fraction
                                !  (yes, JD(1) is still real*8!)
c     outputs:
                                ! Following  are x,y,z posn then x,y,z velocity
      real*8 RCB(6)             ! position of Earth-Moon Barycenter wrt SSB
      real*8 RBE(6)             ! position of Earth wrt Earth-Moon Barycetner
      real*8 RCS(6)             ! position of Sun wrt SSB
      real*8 NUT(4)             ! nutations: 
                                ! nut(1) = d psi (nutation in longitude)
                                ! nut(2) = d epsilon (notation in obliquity)
                                ! nut(3) = d psi dot
                                ! nut(4) = d epsilon dot

c     data buffer:
      real*8 buf(1500)          ! max size of jpl ephemeris record

c     work variables:
      real*8 jda, jdb
      real*8 t(2)
      integer nr
      integer i
      logical bigendian

c     Ephemeris always starts on an MJD that is some integer+0.5, and record
c     length is always an integer number of days.  So set jda to the nearest
c     value of int+0.5 below jd(1)+jd(2), and set jdb to the remaining fraction
c     of a day.  Note:  assumes jd(1) an integer and 0<=jd(2)<1.  


      if (jd(2).ge.0.5d0) then
        jda = jd(1) + 0.5d0
        jdb = jd(2) - 0.5d0
      else
        jda = jd(1) - 0.5d0
        jdb = jd(2) + 0.5d0
      endif

      nr = int((jda-ss(1))/ss(3))+3 ! record number in file

      if (nr.lt.3 .or. jda+jdb.gt.ss(2)) then
        write (*,*) 'Ephemeris date ',jda+jdb,' out of range (',
     +       ss(1),'-',ss(2),')'
        stop
      endif

      if (nr.ne.nrl) then
        nrl = nr
        read (iunit, rec=nr) (buf(i),i=1,ncoeff)
        if(ephbigendian(nephem).and..not.bigendian() .or.
     +     .not.ephbigendian(nephem).and.bigendian()  )
     +        call dbyterev(buf,ncoeff) 
      endif

      t(1) = ((jda-((nr-3)*ss(3)+ss(1)))+jdb)/ss(3) ! fraction within record

      t(2) = ss(3)*86400.d0     ! length of interval in seconds

c     get ephmeris values (EMB=earth-moon barycenter, SSB=solar sys barycenter)
      call interp(buf(ipt(1,11)),t,ipt(2,11),3,ipt(3,11),2,rcs) !sun wrt SSB
      call interp(buf(ipt(1, 3)),t,ipt(2, 3),3,ipt(3, 3),2,rcb) !EMB wrt SSB
      call interp(buf(ipt(1,10)),t,ipt(2,10),3,ipt(3,10),2,rbe) !Moon wrt Earth
      call interp(buf(ipt(1,12)),t,ipt(2,12),2,ipt(3,12),2,nut) !nutations

      do i = 1, 6 ! convert Moon wrt Earth to Earth wrt EMB
        rbe(i) = -rbe(i)/(1.d0+emrat) ! emrat is earth/moon mass ratio
      enddo

      do i = 1, 6               ! convert from kilometers to light-seconds
        rcb(i) = rcb(i) / clight
        rbe(i) = rbe(i) / clight
        rcs(i) = rcs(i) / clight
      enddo

      return
      end
      

      
