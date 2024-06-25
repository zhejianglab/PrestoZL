c      $Id$
      subroutine fit(npts,mode,chisqr,varfit,xmean,ymean,sum,nz,
     +     wmax,lw,ddmch,
     +     buf,npmsav,ksav,
     +     resfile2)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'vcom.h'
	include 'orbit.h'
	include 'tz.h'

c       DJN changes Sep '98
c       --reformat code (indent loops, end all loops with 'continue')
c       --get rid of unneeded 'rmul' calculation
c       --get rid of separate sigmax() calculation since sigmax(i)=array(i,i)
c       --don't calculate gcor()=array(i,i) just before matinv() since 
c         array(i,i)=1.0 always here
c       --add 120/121 loops to un-normalize array() rather than doing so
c         in sigma0 calculation (this way un-normalized array is eventually
c         saved in common block, as desired for covar().)

c       moved declaration of real*8 array(NPA,NPA) to acom.h, djn, 8 Sep 98
	real*8 xmean(NPA),sigmax(NPA),fctn(NPAP1)
	real*8 r(NPA),a(NPA),sigmaa(NPA),gcor(NPA)
	logical lw
	real*8 ddmch(*)
        real*8 buf(*)
        integer npmsav(*), ksav(*)
	character mark*1,adn*14,date*9,damoyr*9
	character*80 resfile2

        integer isetflag

cc       Under linux/g77, output is always flushed.  this
cc       causes huge problems when writing out 80 bytes at
cc       a time to resid2.tmp (unit 32), especially when 
cc       running with a cross-mounted nfs disk. 
cc       A kludgey solution is to gather the resid2.tmp
cc       data into a single memory block and to dump it
cc       at once.  
cc       The data format of fortran output of 72-byte records
cc       is emulated:  in particular, there is a 4-byte integer
cc       with the record length ("72"), followed by the 72 bytes
cc       of actual data, followed by the record length ("72") again,
cc       so that the whole thing ends up being 80 bytes of output
cc       per record.
  
        integer resn(20)  ! small buffer to hold one toa's worth of information
	real*8 resr(9)
        equivalence (resn(2),resr(1)) ! real*8 will not be on 8-byte boundaries
c	integer ipointer, resboff
c	integer resb(1)   ! the big buffer dumped to resid2.tmp; malloc'd below

	integer fd
	integer nwrt
	integer flags, filemode  
 	integer open, close, write

c	rewind 32

	mprt=10**nprnt
	sigma=0.
	chisq=0.

	nterms=nparam-1
	do 29 j=1,nterms
          sigmax(j)=0.
          r(j)=0.
          a(j)=0.
          sigmaa(j)=0.
          do 28 k=1,nterms
            array(j,k)=0.
 28       continue
 29     continue

        ymean=ymean/sum
	do 53 j=1,nterms
          xmean(j)=xmean(j)/sum
 53     continue
	fnpts=npts-nz
	wmean=sum/fnpts

	if(ldesign) rewind(37)
	if(ldesign) write(37) npts, nterms

	do 67 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          weight=weight/wmean
	  if(ldesign) write(37) ct, weight, (fctn(j)-xmean(j), j=1,nterms)
          sigma=sigma+weight*(y-ymean)**2
          do 66 j=1,nterms
            fx=fctn(j)-xmean(j)
            wfx=weight*fx
            sigmax(j)=sigmax(j)+wfx*fx
            r(j)=r(j)+wfx*(y-ymean)
            do 65 k=1,j
              array(j,k)=array(j,k) + wfx * (fctn(k)-xmean(k))
 65         continue
 66       continue
 67     continue

	do 68 i = 1, nterms
	  sigmax(i) = array(i,i)
c If any of the sigmas are zero that probably means something like an 
c empty JUMP or DMX range.  We might as well catch this here instead
c of spewing out garbage later.. PBD 2013/11/20
          if(sigmax(i).eq.0.d0) stop 'fit 68 (unconstrained parameter)'
 68	continue

	free1=fnpts-1
	sigma=dsqrt(sigma/free1)
	do 78 j=1,nterms
          sigmax(j)=dsqrt(sigmax(j)/free1)
          fs=free1*sigmax(j)
          r(j)=r(j)/(fs*sigma)
          do 77 k=1,j
            array(j,k)=array(j,k)/(fs*sigmax(k))
            array(k,j)=array(j,k)
 77       continue
 78     continue
	call matinv(array,nterms,det)
	if(det.eq.0.d0) stop 'fit 78'
	do 80 i=1,nterms
          gcor(i)=sqrt(abs(1.d0 - 1.d0/array(i,i)))
 80     continue
	aa0=ymean

	if(lw) then
	  write(31,1170)
1170	  format(//'Normalized covariance matrix in "number of 9''s"',
     +      ' format:'/)
	  call mxprt(array,gcor,nterms,mfit,nbin,eclcoord)
	  unit=1.d-6
	  if(p0firs.gt.0.1) unit=1.d-3
	  nunit=dlog10(unit)-0.0001
	  write(31,1050) nunit
1050	  format(//'   N     TOBS         Date',
     +     '      FREQ    WGT      NPULSE     R1(P)  R2(',i2,') PH/DT'/)
	endif

	do 106 j=1,nterms
          do 104 k=1,nterms
            a(j)=a(j)+r(k)*array(j,k)
 104      continue
          a(j)=a(j)*sigma/sigmax(j)
          aa0=aa0-a(j)*xmean(j)
 106    continue

c        open(32,file=resfile2,form='unformatted',status='unknown',
c     +		recl = 80*npts)
c	ipointer = mallocxi(resb(1),20*npts,4,resboff)
        flags = isetflag()
	filemode  = 6*64 + 6*8 + 2  ! octal 662 = rw-rw-r--
	if (lw) fd = open(resfile2,flags,filemode)
	resn(1) = 72
	resn(20) = 72

	do 108 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          dt2=y-aa0
          weight=weight/wmean
          do 107 j=1,nterms
            dt2=dt2-a(j)*fctn(j)
 107      continue
          nct=ct
          mark=char(mod(nct,26)+65)
          date=damoyr(int(fmjd))
          dt2sec=dt2*p0firs
          phase=0.
          if(a1(1).ne.0.) then
            if (nbin.ne.9) then
               phase=dmod((ct-t0(1))*86400.d0/pb(1)+1000000.d0,1.d0)
            else
               phase=dmod((ct-t0asc)*86400.d0/pb(1)+1000000.d0,1.d0)
            endif
          endif

          if(nprnt.eq.0.or.mod(i,mprt).eq.1) then
            write(adn,1106) dn
 1106       format(f14.0)

            if(lw) write(31,1107) i,fmjd,mark,date,
     +           rfrq,weight,adn(1:13),y,dt2sec/unit,phase
 1107       format(i4,f12.5,1x,a1,1x,a9,f9.3,f6.3,a13,
     +           f8.3,f8.2,f7.3)
          endif

C Correct tz ref TOA
          if(i.eq.ntzref) then
            ftzrmjd=ftzrmjd-dt2sec/8.64d4
            if (ftzrmjd.lt.0.) then
              ftzrmjd = ftzrmjd + 1
              ntzrmjd = ntzrmjd - 1
            elseif (ftzrmjd.ge.1.) then
              ftzrmjd = ftzrmjd - 1
              ntzrmjd = ntzrmjd + 1
            endif
          endif

c          if((i.lt.npts.or.(.not.gro)).and.lw) write(32,rec=i) 
c     +         ct,dt2,dt2sec,phase,frq,weight,terr,y,ddmch(i)
	  resr(1) = ct
	  resr(2) = dt2
	  resr(3) = dt2sec
	  resr(4) = phase
	  resr(5) = frq
	  resr(6) = weight
	  resr(7) = terr
	  resr(8) = y
	  resr(9) = ddmch(i)
c	  do j = 1, 20
c            resb(20*(i-1)+j+resboff) = resn(j)
c          enddo
          if (lw) nwrt = write(fd,resn,80)
          wmax=max(wmax,weight)
          chisq=chisq+weight*dt2**2
 108    continue
c	write (32,rec=1) (resb(resboff+i),i=1,npts)
c	call freexi(ipointer)
c	close (32)
	if (lw) nwrt = close(fd)

	freen=fnpts-nterms-1
	chisqr=chisq*wmean/freen
	if(mode.eq.0) varnce=chisqr
	if(mode.ge.1) varnce=1./wmean
	varfit=chisq/fnpts
	if(mode.eq.0) chisqr=0.

	do 120 j = 1, nterms
	  do 121 k = 1, nterms
	    array(j,k) = array(j,k)*varnce/(free1*sigmax(j)*sigmax(k))
 121	  continue
 120	continue

	do 133 j=1,nterms
          sigmaa(j)=dsqrt(array(j,j))
 133    continue
	sigma0=varnce/fnpts
	do 145 j=1,nterms
          do 144 k=1,nterms
            sigma0=sigma0+xmean(j)*xmean(k)*array(j,k) 
 144      continue
 145    continue

	sigma0=dsqrt(sigma0)

	do 146 j=1,NPAP1
          freq(j)=0.
          ferr(j)=0.
 146      continue

	freq(1)=aa0
	ferr(1)=sigma0

	do 150 j=1,nterms
          freq(mfit(j+1))=a(j)
          ferr(mfit(j+1))=sigmaa(j)
 150    continue


	return
	end
