c      $Id$
      subroutine tdbread(JD1,JD2,ctatv)

c     tdbread -- read and interpolate TDB-TDT file

c     file should have been written by tdbgen, which uses the Fairhead
c     et al. approximate analytical formula to generate Chebyshev polynomial
c     coefficients, which are read and interpolated here.

      implicit real*8 (a-h,o-z)
      include 'dim.h'
      include 'acom.h'
      include 'tdbcom.h'

c     input variables: 
      real*8 jd1, jd2           ! jd1+jd2 is TDT:
                                ! jd1 is integer part,
                                ! jd2 is fractional part

      real*8 ctatv              ! output TDB-TDT

      real*8 buf(226)            ! max size of data buffer

      real*8 jda, jdb, t(2)

      integer i
      integer nr

      logical bigendian

      save buf                  ! For Linux

c following copied from JPL ephemeris routines -- is it also valid here?

c     Ephemeris always starts on an MJD that is some integer+0.5, and record
c     length is always an integer number of days.  So set jda to the nearest
c     value of int+0.5 below jd1+jd2, and set jdb to the remaining fraction
c     of a day.  Note:  assumes jd1 an integer and 0<=jd2<1.  


      if (jd2.ge.0.5d0) then
        jda = jd1 + 0.5d0
        jdb = jd2 - 0.5d0
      else
        jda = jd1 - 0.5d0
        jdb = jd2 + 0.5d0
      endif

c     Not sure what's going oni with the +2 for old format but +3 for new
c     format.  The +3 is taken from Allan Irwin's time_state.f code; 
c     Not sure where the +2 is from.

      if (tdbif99fmt) then  ! copied from time_state.f in irwin code
        nr = int((jda-tdbd1)/tdbdt)+3 
      else
        nr = int((jda-tdbd1)/tdbdt)+2 ! record number in file;  "+2" skips
                                ! skips over the hdr rec 
      endif

      if (nr.lt.1 .or. jd1+jd2.gt.tdbd2) then
        write (*,*) "Date ",jda+jdb," out of range of TDB-TDT table (",
     +       tdbd1,"-",tdbd2,")"
        stop
      endif

      if (nr.ne.tdbnrl) then
        tdbnrl = nr
c        read (tdbunit, rec=nr) (buf(i),i=1,tdbncf)
        read (tdbunit, rec=nr) (buf(i),i=1,tdbrecsize)
	if (tdbif99fmt.eqv.bigendian()) call dbyterev(buf,tdbncf)
      endif

      if (tdbif99fmt) then  ! copied from time_state.f in irwin code
        t(1) = ((jda-((nr-3)*tdbdt+tdbd1))+jdb)/tdbdt ! fraction within record
        t(2) = tdbdt          
      else
        t(1) = ((jda-((nr-2)*tdbdt+tdbd1))+jdb)/tdbdt ! fraction within record
        t(2) = 1.                 ! unused 
      endif

      call interp(buf(ipt(1,1)),t,ipt(2,1),1,ipt(3,1),1,ctatv)
      if (tdbif99fmt)  then
        ctatv=ctatv*86400         ! Convert days to seconds
        ctatv=ctatv-65.564518e-6  ! Constant offset; see eqn 17 of
                                  ! Irwin & Fukushima 1999, A&A 348: 642.
      endif

      return
      end

