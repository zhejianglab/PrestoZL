	subroutine atimfake(afmjd,nbin,nt,alng,ipsr)
	implicit real*8 (a-h,o-z)
	parameter (nn=31,pi=3.141592653589793d0)
	real*8 xx(nn)
        real*8 maxha

	include 'dim.h'
	include 'tz.h'
        include 'toa.h'
        include 'acom.h'
        integer sitea2n ! external function

        character*20 amjd
        integer i1, i2

	if(oldpar)then
	   write(50,1012) nbin,ncoord
 1012	   format('2',25x,i1,18x,'1',i1)
	   if(nbin.ge.1) then
	      np=6
	   else
	      np=4
	   end if
	   write(50,1000) (params(i),i=1,np)
 1000	   format(a80)
	   read(params(np),1013) tzrsite, tzrfrq, amjd
 1013	   format(a1,14x,f9.0,a20)
c          parse TOA into int+fracion (code from arrtim.f; see notes there)
           i1 = index(amjd,'.')
           i2 = index(amjd(i1+1:20),' ')
           if (i2.eq.0) then
             i2 = 20
           else
             i2 = i1+i2
           endif
           read (amjd(1:i1-1),*) ntzrmjd
           read (amjd(i1:i2),*) ftzrmjd
           if (ntzrmjd.lt.30000) ntzrmjd=ntzrmjd+39126   !convert to MJD

	endif

c       store reference TOA line
        stflag = .true.
        stntoa = 1
        ntzrsite = sitea2n(tzrsite)
        stnsite(stntoa) = ntzrsite
        stfrq(stntoa) = tzrfrq
        stnmjd(stntoa) = ntzrmjd
        stfmjd(stntoa) = ftzrmjd
        sterr(stntoa) = 0.
        stddm(stntoa) = 0.

	nspan=nsp(ipsr)
	maxha=mxha(ipsr)
	tzfreq=tzof(ipsr)
	if(tzfreq.lt.0.)tzfreq=tzrfrq

        if (autotz) then
          fmjd1 = afmjd + nspan/2880.
          nmjd1 = int(fmjd1)
          fmjd1 = fmjd1 - nmjd1
          ! Round fmjd1 to nearest second.  Not really necessary, but it will bring
          ! hhmmss and mjd dates in polyco.dat into closer agreement if tz start
          ! time is not an integer second.
          fmjd1 = int(fmjd1*86400+0.5)/86400.d0
          nsets = (maxha*60+(nspan-1))/nspan  ! the "nspan-1" forces rounding up
        else 
          hlst=24.d0*dmod(1.002737909d0*afmjd+0.154374d0 -
     +         alng/6.2831853071795864d0,1.d0)	
          read(pname,1030) irah,iram		
 1030     format(2i2)
          rax=irah+(iram+2)/60.			
          wait=(rax-hlst)*0.99727 !Solar to sidereal
          if(wait.lt.-maxha) wait=wait + 23.9345
          fmjd1=afmjd+(wait-maxha)/24.+nspan/2880. !Compute start time
          nmjd1 = int(fmjd1)    ! break into integer, fractional parts
          fmjd1 = fmjd1 - nmjd1
          fmjd1=nint(48*fmjd1)/48.d0 ! round to nearest half-hour?!
          nsets=(120*maxha+nspan-1)/nspan
        endif
 
	b=nspan/2 + 5
	a=-b
	bma=0.5*(b-a)
	bpa=0.5*(b+a)
	do 30 k=1,nn
          xx(k)=cos(pi*(k-0.5)/nn)*bma+bpa
 30     continue

	i=0
	do 50 j=1,nsets				
          fmjd2 = fmjd1 + (j-1)*nspan/1440.d0
          ntmp1 = int(fmjd2*1.d8) ! round fmjd2 to 1.e-10 (messy since ints
          ntmp2 = int((fmjd2 *1.d8-ntmp1)*1.d2) ! don't have 10 significant
          fmjd2 = ntmp2*1.d-10 + ntmp1*1.d-8 !    digits)
          do 40 k=1,nn
            i = i + 1
            if(i.gt.NTZARR) stop ' Nspan too small'
            ftmjd(i) = fmjd2 + xx(k)/1440.d0
            ntmjd(i) = nmjd1
            if (ftmjd(i).ge.2.) then
              ftmjd(i) = ftmjd(i) - 2.
              ntmjd(i) = ntmjd(i) + 2
            elseif (ftmjd(i).ge.1.) then
              ftmjd(i) = ftmjd(i) - 1.
              ntmjd(i) = ntmjd(i) + 1
            elseif (ftmjd(i).lt.0.) then
              ftmjd(i) = ftmjd(i) + 1.
              ntmjd(i) = ntmjd(i) - 1
            endif
            tmin(i) = 1440.d0*((ntmjd(i)-ntmjd(1))+ftmjd(i)-ftmjd(1))
            stntoa = stntoa + 1
            stnsite(stntoa) = sitea2n(tzsite)
            stfrq(stntoa) = tzfreq
            stnmjd(stntoa) = ntmjd(i)
            stfmjd(stntoa) = ftmjd(i)
            sterr(stntoa) = 0.
            stddm(stntoa) = 0.            
 40       continue
 50	continue

	nt=i
	return
	end
