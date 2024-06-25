      program avtime

c     reduces tempo files to daily averages
c     revision history
c     feb 90  djn  ????
c     11 may 92  djn  now follows "INCLUDE" statements up to 10 deep
c     28 oct 92  djn  variable precision of "error" quantity on output
c     29 nov 92  djn  fixed bug--didn't do nbin>1 header properly
c        ???     djn  modified to use MJD's as well as fut's
c     19 feb 93  djn  modified to allow slightly shifted TOA column
c     10 aug 95  djn  shifted human-readable date output by one column
c     21 oct 97  djn  skip old and new style header
c     02 nov 00  djn  read ITOAs;  not extensively tested

c     known bugs/features:
c     confused by "PHAS" statements--can't handle zero-weight points


      parameter (NPTSMAX=1000)

      character nsite(NPTSMAX)
      integer*4 nsc(NPTSMAX),nfut(NPTSMAX)
      logical isend,rotated,incflag
      real*8 dname,fut(NPTSMAX),w(NPTSMAX),dt2(NPTSMAX)
      real*8 freq(NPTSMAX),ffut(NPTSMAX)
      real*8 dnamez
      real*8 parm,ffadd,ffut0
      real*8 dtave, sum, sq, sumwt, ave, sdmean, sdev
      real*8 err
      real*8 ct,dtp,phase,period,fdop
      character*9 date(NPTSMAX)
      character*80 infile,outfile
      character*80 card
      character*80 fmtstring
      integer i, ii, j, idx, nargs, nmid
      integer nav,imax
      integer navmin, navdef
      integer ilen
      real*8 tmin, tmindef
      real*8 resmax, resdef
      integer nhdr
      integer nfile
      logical newstyle
      logical checkcard    ! external subroutine
      logical skip
      character*2 bsite
      character*9 aterr
      real*8 emax, terr
      character*15 ffuts
      character*8 ffutsfmt
      integer nffuts

      idx=0
      isend=.false.
      incflag = .false.
      tmindef=500.0
      resdef=0.1
      navdef=3

      emax = 1.e20

      nargs=iargc()
      if (nargs.eq.2) then
        call getarg(1,infile)
        call getarg(2,outfile)
      else
        infile = 'raw.dat'
        outfile = 'avtime.dat'
      endif
      nfile = 20
      open(nfile,file=infile)
      open(30,file=outfile)
      open(32,file='resid2.tmp',form='unformatted')
      rewind 20
      rewind 30
      rewind 32

      write (*,*) 'Daily TOA/Tempo averager: Version 2.X'
      write (*,*)
      write (*,*) 'tmin (<ret>=',tmindef,')  ?'
      read '(F10.0)',tmin
      write (*,*) 'resmax (<ret>=',resdef,')  ?'
      read '(F10.0)',resmax
      write (*,*)  'navmin (<ret>=',navdef,')  ?'
      read '(I10)',navmin
      if (tmin.eq.0.) tmin = tmindef
      if (resmax.eq.0) resmax = resdef
      if (navmin.eq.0) navmin = navdef
      write(6,100) tmin,resmax,navmin
 100  format(' tmin=',f6.2,'  resmax=',f10.7,'   navmin=',i3//)

c     skip over Tempo header
      newstyle = .false.
      read (nfile,fmt='(A80)',end=999) card
      if (checkcard(card,'HEAD').or.checkcard(card,'#')) 
     +                      newstyle = .true.
      write (30,fmt='(A80)') card
      nhdr = 4
      if (card(27:27).ge.'1') nhdr = 6
      do 200 i = 2, nhdr
        read (nfile,fmt='(A80)',end=999) card
        if (checkcard(card,'HEAD')) newstyle = .true.
        if (card.eq.'HEAD') newstyle = .true.
        write (30,fmt='(A80)') card
        if (checkcard(card,'TOA')) goto 400 ! new style header, done reading it
 200  continue

      if (newstyle) then   ! still more header to read
 300    continue
        read (nfile,fmt='(A80)',end=999) card
        write (30,fmt='(A80)') card
        if (.not.checkcard(card,'TOA')) goto 300
      endif

 400  continue

c     used to monitor Tempo input "TIME" corrections
      ffadd = 0.

c     read through file to get one day's residuals
 500  continue

      skip = .false.

      do 1020 i=1,NPTSMAX  

 1005   read(nfile,fmt='(A80)',end=999) card

        if (card(1:4).eq.'SKIP') then
          skip = .true.
          goto 1005
        endif
        if (card(1:6).eq.'NOSKIP') then
          skip = .false.
          goto 1005
        endif
        if (skip) goto 1005
        if (card(1:4).eq.'TIME') then
          read (card,1007) parm
 1007     format (5x,f20.0)
          ffadd = ffadd + parm/86400.
          goto 1005
        elseif (card(1:4).eq.'EMAX') then
          read (card,1007) emax
          goto 1005
        elseif (card(1:7).eq.'INCLUDE') then
          incflag = .true.
          nav = idx
          goto 1051
        elseif (card(1:3).eq.'END') then
          goto 999
        elseif (card(1:1).eq.'C') then
          goto 1005
        elseif (card(1:1).eq.'#') then
          goto 1005
        else
                                ! Figure out whether card is a TOA
                                !   (Side benefit: we also figure out
                                !    here what format if it is a TOA)
          
          nfmt = -1		! default == not a TOA
          if ((card(1:1).ge.'1'.and.card(1:1).le.'9')
     +         .or.(card(1:1).ge.'a'.and.card(1:1).le.'z') ) then
            if (card(2:2).eq.' ') then
              nfmt = 0          ! Traditional Princeton Format
            else                ! Could be P'ton or ITOA
              do j = 1, 80
                if (card(j:j).eq." ") then
                  nfmt = 0
                  goto 1009
                elseif (card(j:j).eq."+") then
                  nfmt = 2
                  goto 1009
                elseif (card(j:j).eq."-") then
                  nfmt = 2
                  goto 1009
                elseif (card(j:j).lt."0".or.card(j:j).gt."9") then
                  goto 1009
                endif
              enddo
 1009         continue
            endif
          endif
          if ((card(80:80).ge.'1'.and.card(80:80).le.'9')
     +         .or.(card(80:80).ge.'a'.and.card(80:80).le.'z') ) 
     +         nfmt = 1         ! Parkes/Jodrell Format
          if (nfmt.eq.-1) then
            nav = idx
            goto 1051
          endif
        endif


        idx=idx+1
        if (nfmt.eq.0) then     ! Princeton Format
          if (card(30:30).eq.'.') then
            read(card,1010) nsite(idx),nsc(idx),dname,
     +           freq(idx),nfut(idx),ffuts,aterr,date(idx)
 1010       format(a1,i5,1x,a8,f9.0,i5,a15,a9,1x,a9)
            ! following avoids arcane problem in some TOA lines
            ! where there's extra data within the columns allocated
            ! for the fut...the data follows the fut and a space...
            ! problem seen in 1913 mark1 data (djn-8 August 2007)
            n = index(ffuts,' ')
            if (n.eq.0) then
              read(ffuts,fmt='(f15.14)') ffut(idx)
            else
              write(ffutsfmt,1011) n-1, n-2
 1011         format("(f",i2.2,".",i2.2,")")
              read(ffuts(1:index(ffuts,' ')-1),fmt=ffutsfmt)ffut(idx)
            endif
          else
c     (assume card(31;31).eq.'.')
            read(card,1013) nsite(idx),nsc(idx),dname,
     +           freq(idx),nfut(idx),ffut(idx),aterr,date(idx)
 1013       format(a1,i5,1x,a8,f9.0,i6,f14.13,a9,1x,a9)
          endif
c         mysterious code for interpreting aterr, straight out
c         of tempo arrtim.f
          iz = 1
          do ii = 1, 9         
            if (aterr(ii:ii).eq.' ') iz=ii
          end do
          terr=0.
          if (iz.le.9) read(aterr(iz:9),*,err=1014,end=1014) terr
 1014     continue

        else if (nfmt.eq.1) then !Parkes/Jodrell
          read(card,1015) dname,freq(idx),
     +         nfut(idx),ffut(idx),terr,nsite(idx)
          nsc(idx) = 0  ! don't try to translate scan numbers
          ! sometimes they are in columns 8-12
 1015     format(7x,5x,5x,a8,f9.0,i7,f14.13,8x,f8.1,8x,a1)
          call modayrs(nfut(idx),date(idx))

        else                    ! Must be ITOA
          read(card,1017) nfut(idx),ffut(idx),terr,freq(idx),bsite
 1017     format(9x,i5,f14.13,f6.2,f11.4,10x,2x,a2)
          if (bsite.eq."AO") then
            nsite(idx) = "3"
          else
            print *,"Can't handle site ",bsite," in ITOA format"
          endif
        endif

        if (terr.gt.emax) print *, "YES ",terr,emax
        if (terr.gt.emax) goto 1005

        if (nfut(idx).lt.20000) nfut(idx)=nfut(idx)+39126

c     make TIME adjustments to ffut, nfut:
        ffut(idx) = ffut(idx) + ffadd
        if (ffut(idx).ge.0.d0) then
          nfut(idx) = nfut(idx) + int(ffut(idx))
          ffut(idx) = ffut(idx) - int(ffut(idx))
        else
          nfut(idx) = nfut(idx) + int(ffut(idx)) - 1
          ffut(idx) = ffut(idx) - (int(ffut(idx)) - 1)
        endif
c     create fut() array for approximate calculations:
        fut(idx) = ffut(idx) + nfut(idx)
        if (dabs(freq(idx)-freq(1)).gt.100.d0 .or. 
     +       dabs(fut(idx)-fut(1)).gt.tmin/1440.) goto 1030
        dnamez=dname
        read(32) ct,dtp,dt2(idx),phase,fdop,w(idx)
c     throw out points with zero weight
        if (w(idx).eq.0.) idx = idx - 1
        if (dtp.ne.0.) period = dt2(idx)/dtp
 1020 continue
      write (*,*) "Avtime can only average ",NPTSMAX," points. To "
      write (*,*) "analyze more points, change NPTSMAX in avtime.f"
      write (*,*) "or modify the input parameters with this data set."
 1030 continue
      backspace nfile
 1040 nav = idx - 1
 1051 continue
      idx = 0
      if (nav.lt.navmin) goto 99
      nmid = (nav+1)/2


c     calculate statistics on day's residuals
      rotated = .false.
 1052 sum = 0.
      sq = 0.
      sumwt = 0.
      do 1055 i = 1, nav 
        sumwt = sumwt + w(i)
        sum = sum + w(i)*dt2(i)
        sq = sq + w(i)*dt2(i)**2
 1055 continue
      ave = sum/sumwt
      imax=nav
      sdev = dsqrt((sq-sum*sum/sumwt)/sumwt)
      if (sdev.ge.0.25*period) then
        if (rotated) then
          write (*,1081) nfut(nmid)
 1081     format (' Day ',i5,': std dev > 1/4*period')
          goto 23024
        else
          rotated = .true.
          do 1090 i = 1, imax
            if (dt2(i).gt.0.5*period) 
     +           dt2(i) = dt2(i) - period
 1090     continue
          goto 1052
        endif
      endif


c     eliminate outriders
      do 1080 i = 1, imax
        if (dabs(dt2(i)-ave).lt.resmax) goto 1080
        write(6,1070) fut(i), dt2(i)-ave
 1070   format(' bad point -- ',f20.11,'   resid=',f12.9)
        sumwt = sumwt - w(i)
        sum = sum - w(i)*dt2(i)
        sq = sq - w(i)*dt2(i)**2
        nav = nav - 1
 1080 continue

      if(nav.lt.navmin) goto 23024
      dtave = sum/sumwt
      sdev = dsqrt((sq-sum*sum/sumwt)/sumwt)

      sdmean = sdev/sqrt(nav-1.0)
      err = sdmean*1.d6
      ffut0 = ffut(nmid) + (dtave-dt2(nmid))/86400.
      if (ffut0.gt.1.0d0) then
        nfut(nmid)=nfut(nmid)+1
        ffut0=ffut0-1.d0
      else if (ffut0.lt.0.d0) then
        nfut(nmid)=nfut(nmid)-1
        ffut0=ffut0+1.d0
      endif
      fmtstring = '(a1,i5,1x,a8,f9.3,i5,f15.13,f9.2,x,a9,x,i4,f11.7)'
      if (freq(nmid).lt.1.d3) fmtstring(17:17)='4'
      if (err.ge.1.d3) fmtstring(32:32)='1'
      if (err.ge.1.d4) fmtstring(32:32)='0'
      write (6,fmtstring) nsite(nmid),nsc(nmid),dnamez,freq(nmid),
     +     nfut(nmid)/10,mod(nfut(nmid),10)+ffut0,err,
     +     date(nmid),nav,dtave
      write (30,fmtstring) nsite(nmid),nav,dnamez,freq(nmid),
     +     nfut(nmid)/10,mod(nfut(nmid),10)+ffut0,err,
     +     date(nmid)
23024 continue
 99   continue
      if (incflag) then
        ilen = index(card(9:80),' ')
        if (ilen.eq.0) ilen = 72
        infile = card(9:ilen+9)
        nfile = nfile + 1
        open(nfile,file=infile)
        incflag = .false.
      endif
      if (isend.and.nfile.lt.20) stop
      isend = .false.
      if (card(1:1).gt.'9'  .and.  card(1:4).ne.'TIME'
     +     .and.  card(1:7).ne.'INCLUDE'
     +     .and.  card(1:1).ne.'C'    .and.
     +     (card(1:1).lt.'a'.or.card(1:1).gt.'z') ) call wtrim(card)
      goto 500

c     end of file
 999  isend = .true.
      idx = idx + 1
      close (nfile)
      nfile = nfile - 1
      goto 1040
      end

C----------------------------------------------

      subroutine wtrim(line)

      character*80 line
      integer i, len

      do 3010 i=80,1,-1
        if(line(i:i).ne.' ') go to 3020
 3010 continue

 3020 len=i
      if(len.gt.0) write(30,3030) (line(i:i),i=1,len)
      if(len.gt.0) write(6,3030) (line(i:i),i=1,len)
 3030 format(80a1)
      return
      end

C--------------------------------------------------------------------------

      SUBROUTINE MODAYRS(IMJD,S)

      IMPLICIT REAL*8 (A-H,O-Z)
      character*9 s
      INTEGER*4 MODAY(12)
      CHARACTER*3 MONTH(12),MO
      common /MDY/ month,moday
      DATA MONTH/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     +  'Oct','Nov','Dec'/
      DATA MODAY/31,28,31,30,31,30,31,31,30,31,30,31/

C  XJD IS THE JULIAN DAY.
C  IDAY IS THE NUMBER OF DAYS SINCE JAN 1.0, 1966.


      IDAY=IMJD-39126
      IYR=1966
      DO 10 I=1,33
      ND=365
      IF(MOD(IYR,4).EQ.0) ND=366
      IF(IDAY.LE.ND) GO TO 20
      IDAY=IDAY-ND
10    IYR=IYR+1

20    MODAY(2)=28
      IF(ND.EQ.366) MODAY(2)=29
      DO 30 J=1,12
      IF(IDAY.LE.MODAY(J)) GO TO 40
30    IDAY=IDAY-MODAY(J)
      J=12
40    MO=MONTH(J)

      iyr = iyr - 100*int(iyr/100)
      write (s,fmt='(i2,"-",a3,"-",i2.2)') iday, mo, iyr

      RETURN
      END


C*************************************************************************** 

      subroutine upcase(w)

C  Converts string (up to first blank) to upper case.

      character*(*) w
      do 10 i=1,len(w)
         if(w(i:i).eq.' ') go to 20
         j=ichar(w(i:i))
         if(j.ge.97.and.j.le.122) w(i:i)=char(j-32)
 10   continue
 20   return
      end


C***************************************************************************

      logical function checkcard(a,c)

c     check to see whether a(1:80) has c as key (first non-whitespace string)

      implicit none

      character*80 a
      character*(*) c

      character*80 b
      integer j
      integer k

      j = 1
      call citem(a,80,j,b,k)
      call upcase(b)
      checkcard = .false.
      if (b(1:len(c)).eq.c) checkcard = .true.

      return
      end
      


C***************************************************************************

	subroutine citem(line,ll,j1,item,li)

c  Searches line(j1:ll) for next .not.(space.or.tab)
c  Returns item of non-blank length li, and j1 = index of following character.
c  RNM March, 1995

        implicit none
        integer ll,j,jn,j1,li
        character line*(*),item*(*)

	item=' '
	do j=j1,ll
	  if(line(j:j).ne.' '.and.ichar(line(j:j)).ne.9)go to 10
	enddo
	li=0         ! No non-space character found
	return

10	jn=j

	do j=jn,ll
	  if(line(j:j).eq.' '.or.ichar(line(j:j)).eq.9)go to 20
	enddo
        j1=ll
        go to 30

 20	j1=j-1

 30	li=j1-jn+1
        item(1:li)=line(jn:j1)
        j1=j1+1

        return
	end
