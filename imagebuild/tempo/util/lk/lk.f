	program lk

	implicit real*8 (a-h,o-z)
	parameter(maxpts=40000)
	real*8 rct(maxpts),rphase(maxpts),resid1(maxpts),resid2(maxpts),
     +    xindex(maxpts),dayofyr(maxpts)
	character*1 an
	character*32 infile
	data infile/'resid2.tmp'/

	if(iargc().eq.1) call getarg(1,infile)
        if (infile.eq."-h") goto 9998
	open(32,file=infile,form='unformatted',status='old',err=9997)

1	rewind 32
	ctmin=1.d30
	ctmax=0.

	do 3 i=1,maxpts
	read(32,end=4) ct,dt2,dt2sec,phase,freq,w,terr,y
	rct(i)=ct
	if(ct.lt.ctmin) ctmin=ct
	if(ct.gt.ctmax) ctmax=ct
	resid1(i)=y
	p0=dt2sec/dt2
	rphase(i)=phase
	resid2(i)=dt2
	xindex(i)=i
3	dayofyr(i)=mod(ct,365.25d0)

4	npts=i-1

10	xnpts=npts

	r1=0.
	r2=0.
	do 11 i=1,npts
	r1=dmax1(r1,dabs(resid1(i)))
11	r2=dmax1(r2,dabs(resid2(i)))

12	print*,'Enter command: '
	read(*,1014) an
1014	format(a1)
	if(an.eq.'x') stop
	if(an.eq.'0') stop

	if(an.ne.'1') go to 20
	print*,'Pre-fit residuals vs. date'
	call plt(rct,ctmin,ctmax,resid1,-r1,r1,npts,rct,p0)

20	if(an.ne.'3') go to 30
	print*,'Pre-fit residuals vs. orbit phase'
	call plt(rphase,0.d0,1.d0,resid1,-r1,r1,npts,rct,p0)

30	if(an.ne.'2') go to 40
	print*,'Post-fit residuals vs. date'
	call plt(rct,ctmin,ctmax,resid2,-r2,r2,npts,rct,p0)

40	if(an.ne.'4') go to 50
	print*,'Post-fit residuals vs. orbit phase'
	call plt(rphase,0.d0,1.d0,resid2,-r2,r2,npts,rct,p0)

50	if(an.ne.'5') go to 60
	print*,'Pre-fit residuals serially'
	call plt(xindex,1.0d0,xnpts,resid1,-r1,r1,npts,rct,p0)

60	if(an.ne.'6') go to 70
	print*,'Post-fit residuals serially'
	call plt(xindex,1.0d0,xnpts,resid2,-r2,r2,npts,rct,p0)

70	if(an.ne.'7') go to 80
	print*,'Pre-fit residuals vs. day of year'
	call plt(dayofyr,0.d0,365.25d0,resid1,-r1,r1,npts,rct,p0)

80	if(an.ne.'8') go to 100
	print*,'Post-fit residuals vs. day of year'
	call plt(dayofyr,0.d0,365.25d0,resid2,-r2,r2,npts,rct,p0)

100	go to 12

9997    continue
        print *, "Error: no residual file"
        print *, "Type 'lk -h' if you need instructions"
        goto 9999

9998    continue
        print *, "Usage:"
        print *, "  lk [residual file | -h]"
        print *, "  residual file: must be output by tempo"
        print *, "                 default is resid2.tmp"  
        print *, "  -h: help message"
        print *, ""
        print *, "Typically one runs tempo and then runs lk in the"
        print *, "  same directory to see the residuals."
        print *, ""
        print *, "At the lk 'Enter command' prompt, enter:"
        print *, "  0: quit"
        print *, "  1: pre-fit residuals vs date"
        print *, "  2: post-fit residuals vs date"
        print *, "  3: pre-fit residuals vs orbital phase"
        print *, "  4: post-fit residuals vs orbital phase"
        print *, "  5: pre-fit residuals vs serial number"
        print *, "  6: post-fit residuals vs serial number"
        print *, "  7: pre-fit residuals vs day of year"
        print *, "  8: post-fit residuals vs day of year"
     
9999    continue
	end

	subroutine plt(x,xmin,xmax,y,ymin,ymax,npts,ct,p0)

	implicit real*8 (a-h,o-z)
	real*8 x(npts),y(npts),ct(npts)
	character*1 iplt(70,17)
	character*26 letter
	data letter/'abcdefghijklmnopqrstuvwxyz'/

	ys=16.999/(ymax-ymin)
	jz=-ymin*ys+1.5
	jz=18-jz
	if(jz.lt.1) jz=1
	if(jz.gt.17) jz=17
	do 20 k=1,70
	do 10 j=1,17
10	iplt(k,j)=' '
20	iplt(k,jz)='-'

	sum=0.
	sq=0.
	xs=69.999/(xmax-xmin)
	do 50 i=1,npts
	sum=sum+y(i)
	sq=sq+y(i)**2
	k=xs*(x(i)-xmin)+1.0001
	j=ys*(y(i)-ymin)+1.5
	j=18-j
	if(j.lt.1) j=1
	if(j.gt.17) j=17
	mark=mod(ifix(sngl(ct(i))),26)+1
50	iplt(k,j)=letter(mark:mark)

	rms=dsqrt(sq/npts - (sum/npts)**2)
	write(6,1030)
1030	format(1x,'+0------1------2------3------4------5------6',
     +    '------7------8------9------+')
	write(6,1040) iplt
1040	format(' |',70a1,'|')
	write(6,1030)
	write(6,1050) ymax,ymax*p0*1.e6,rms,rms*p0*1.e6
1050	format(' Max:',f15.9,' P0,',f9.3,' us.',5x,'RMS:',f15.9,
     +     ' P0,',f9.3,' us.')

	return

	end
