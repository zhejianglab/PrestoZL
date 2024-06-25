c      $Id$
	subroutine clockcor(fmjd,nsite,nn,deltat,clk,nfmt)

c       modified DJN 22 Jul '06: if nfmt=2 (itoa), allow
c       UTC->TT(BIPM) correction

	implicit real*8 (a-h,o-z)
	save
	character*9 date,damoyr,datez2,datez3,csite*35
	include 'dim.h'
	include 'clocks.h'
	include 'acom.h'
	data lsite/0/,ii/1/,td1/0.d0/,nmsg/0/,maxmsg/20/
        data nutcmsg/0/
	data csite/'123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

c SECTION 1: Observatory to NIST correction
c File entries at UT 00h interpolated to fmjd. Works for different sites 
c and out-of-order toas. Clock entries for a given site must be in increasing 
c time order, but entries for different sites can be in any order.
	clk1=0.
	if(nclk.gt.0.and.nfmt.ne.2)then
	   if(nsite.ne.lsite)then
	      lsite=nsite
	      ii=1
	   endif
	   if(fmjd.lt.td1)ii=1
	   if (ii.eq.1) td1=0.d0
	   
 10	   td2=tdate(ii)
	   if(jsite(ii).ne.nsite)then
	      ii=ii+1
	      if(ii.gt.ndate)then
		 nmsg=nmsg+1
		 if(nmsg.le.maxmsg.and..not.quiet) then
		    write(*,1005) fmjd,nsite
		    write(31,1005) fmjd,nsite
 1005		    format(' Clk corr error: TOA at',f12.5,
     +             ' greater than last clk file entry for nsite ',i3)
		 endif
		 if(nmsg.eq.maxmsg+1.and..not.quiet) print*,
     +      '*** Additional clock-correction messages suppressed. ***'
		 ii=1
		 clk1=0.d0
		 go to 56
	      endif
	      go to 10
	   else
	      if(fmjd.gt.td2)then
		 td1=td2
		 ckc1=ckcorr(ii)
		 flag1=ckflag(ii)
		 ii=ii+1
		 if(ii.gt.ndate)then
		    nmsg=nmsg+1
		    if(nmsg.le.maxmsg) then
		       write(*,1005) fmjd,csite(nsite:nsite)
		       write(31,1005) fmjd,csite(nsite:nsite)
		    endif
		    if(nmsg.eq.maxmsg+1) print*,
     +    '*** Additional clock-correction messages suppressed. ***'

		    ii=1
		    clk1=0.d0
		    go to 56
		 endif
		 go to 10
	      else
		 if(ii.eq.1)then
		    write(*,1010)fmjd,td2,nsite
 1010		    format(' Clk corr error: TOA at',f12.5,' less ',
     +             ' than first clk entry',f12.5,' for nsite ',i3)
		    stop
		 end if
		 if (ckflag(ii).eq.0 .and. flag1.eq.0) then
		    slope=(ckcorr(ii)-ckc1)/(td2-td1)
		    clk1=(ckc1+slope*(fmjd-td1))*1.d-6
		 else
		    off1 = fmjd-td1
		    off2 = td2-fmjd
		    if (off1.gt.1 .and. off2.gt.1) then
		       nmsg=nmsg+1
		       if(nmsg.le.maxmsg.and..not.quiet) then
			  write (*,1020) fmjd,csite(nsite:nsite),
     +		             td1,td2
	                write (31,1020) fmjd,csite(nsite:nsite),
     +                       td1,td2
 1020                     format('No clock corr for TOA at',f9.2,
     +                       ' site ',
     +                       a1,'; nearest at ',f9.2,' and ',f9.2)
		       endif
		       if(nmsg.eq.maxmsg+1.and..not.quiet) print*,
     +    '*** Additional clock-correction messages suppressed. ***'
		       clk1 = 0
		    elseif (off1.lt.off2) then
		       clk1 = ckc1*1.d-6
		    else
		       clk1 = ckcorr(ii)*1.d-6
		    endif
		 endif
		 goto 56
	      endif
	   endif
	endif

56	continue 

c SECTION 2: This has two uses which give the opposite sign:
c  1.  For "normal" observatory TOAs which were converted to UTC(NIST) in 
c      Section 1, this does a further conversion from UTC(NIST) to UTC
c  2.  For TOAs already in UTC (i.e., TOAs in ITOA format), this section
c      can do a back conversion from UTC back to UTC(NIST);  this is useful
c      only if section 3 is invoked (which does UTC(NIST) to "other")
c      and hence is only run if nclk==3 [UTC(NIST)] AND nfmt==2 [ITOA]

	clk2=0.d0				
	clkoff=0.d0      
	if(nclk.eq.2.and.nfmt.ne.2 .or. nclk.ge.3.and.nfmt.eq.2) then
	   do i=1,NPT
	      if(abs(fmjd-tutc(i)).le.5.d0) then
		 if (fmjd.lt.tutc(i)) then
		    if (i.eq.1) then
		       clkoff=utcclk(i)
		    else
		       slope=(utcclk(i)-utcclk(i-1))/(tutc(i)-tutc(i-1))
		       const=utcclk(i)-slope*tutc(i)
		       clkoff=slope*fmjd+const
		    endif
		 else
		    if (tutc(i+1).lt.1.) then
		       clkoff=utcclk(i)
		    else
		       slope=(utcclk(i+1)-utcclk(i))/(tutc(i+1)-tutc(i))
		       const=utcclk(i)-slope*tutc(i)
		       clkoff=slope*fmjd+const
		    endif
		 endif

		 clk2=clkoff*1.d-6
		 go to 49
	      endif
	      if(tutc(i).eq.0.d0) go to 48
	   enddo

 48	   date=damoyr(int(fmjd))
           if(date.ne.datez2.and..not.quiet) then
             nutcmsg = nutcmsg + 1
             if (nutcmsg.le.maxmsg) then
               write(*,1055) nn,date,clkfile(2)(1:40)
 1055          format(i5,2x,a9,' No entry in clock file ',a40)
             else if (nutcmsg.eq.maxmsg+1) then
               write(*,*) " Additional messages suppressed for ",
     +                         clkfile(2)(1:40)
             endif
           endif
	   datez2=date
	endif
49      continue
        if (nfmt.eq.2) clk2=-clk2     !  backconvert UTC->UTC(NIST)

c  SECTION 3: UTC to TT(BIPM)
	clk3=0.d0				
	if(nclk.ge.3) then
	   do i=1,NPT
	      if(abs(fmjd-tbipm(i)).le.5.d0) then
		 if (fmjd.lt.tbipm(i)) then
		    if (i.eq.1) then
		       clkoff=bipm(i)
		    else
		       slope=(bipm(i)-bipm(i-1))/(tbipm(i)-tbipm(i-1))
		       const=bipm(i)-slope*tbipm(i)
		       clkoff=slope*fmjd+const
		    endif
		 else
		    if (tbipm(i+1).lt.1.) then
		       clkoff=bipm(i)
		    else
		       slope=(bipm(i+1)-bipm(i))/(tbipm(i+1)-tbipm(i))
		       const=bipm(i)-slope*tbipm(i)
		       clkoff=slope*fmjd+const
		    endif
		 endif
		 clk3 = clkoff/1.d6
		 go to 59
	      endif
	      if(tbipm(i).eq.0.d0) go to 58
	   enddo

 58	   date=damoyr(int(fmjd))
	   if(date.ne.datez3) write(*,1055) nn,date,clkfile(nclk)(1:40)
	   datez3=date
	endif

59	clk=(clk1+clk2+clk3+deltat)/86400.d0

	return
	end
