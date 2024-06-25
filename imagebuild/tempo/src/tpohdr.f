c      $Id$
	subroutine tpohdr(oldpar,pardir,parfile,ncoord,t0,pb,p0,dm,
     :      nbin,ipsr)


	implicit real*8 (a-h,o-z)
	character pardir*(*),parfile*(*),path*160
	character*12 zname
	logical oldpar

	include 'dim.h'
	include 'tz.h'
 
	pname=name(ipsr)
	if(oldpar)then
	   kd=index(tzdir,' ')-1
	   path=tztfile
	   open(12,file=path,status='old',err=998)
	   
 10	   read(12,1000,end=999) params(1)		
 1000	   format(a)
	   read(12,1000,end=999) (params(i),i=1,3)	
	   zname=params(1)(1:8)
	   if (params(2)(1:1).eq.'B') then
	      nbin=1
	      ncoord=0
	   else if (params(2)(1:1).eq.'P') then
	      nbin=0
	      ncoord=0
	   else if (params(2)(1:1).eq.'Q') then
	      nbin=0
	      ncoord=1
	   else if (params(2)(1:1).eq.'C') then
	      nbin=1
	      ncoord=1
	   else if (params(2)(1:1).eq.'D') then
	      nbin=9
	      ncoord=1
	   else if (params(2)(1:1).eq.'E') then
	      nbin=10
	      ncoord=1
	   else if (params(2)(1:1).eq.'F') then
	      nbin=3
	      ncoord=0
	   else if (params(2)(1:1).eq.'G') then
	      nbin=3
	      ncoord=1
	   else
	      print *,' Pulsar type ',params(2)(1:1),
     +   	   ' illegal, check tztot.dat'
	      stop
	   end if
	   np = 4
	   if (nbin.eq.0) then
	      np = 4
	   else if (nbin.eq.1 .or. nbin.eq.3) then
	      np = 6
	   else if (nbin.eq.9) then
	      np = 8
	   else if (nbin.eq.10) then
	      np = 10
	   endif
	   read(12,1000)(params(i),i=4,np)

	   if(zname(1:8).ne.pname(1:8))go to 10
 
	   close(12)
	   read(params(2),1011) p0			
 1011	   format(1x,f19.0)
	   read(params(3),1012) dm			
 1012	   format(8x,f20.0)
	   pb=0
	   if(nbin.gt.0) read(params(4),1013) t0,pb	
 1013	   format(40x,2f20.0)

	else                                       ! New-style parameters

	   if(parfile.eq.'def')then
	      kd=index(pardir,' ')-1
	      kn=index(pname,' ')-1
	      path=pardir(1:kd)//pname(1:kn)//'.par'
	   else
	      path=parfile
	   endif
	   open(49,file=path,status='old',err=998)
	   call zeropar(nits)	! To get TZ ref params
	   call rdpar(nits)	! To get TZ ref params
	   rewind 49
	endif

	name(ipsr)='done'
	return

 998	write(*,'(''Failed to open '',a)')path
 
 999	return

	end




