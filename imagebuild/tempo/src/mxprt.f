c      $Id$
        subroutine mxprt(a,gcor,nn,mfit,nbin,eclcoord)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*8 a(NPA,NPA),sig(NPA),gcor(NPA)
	character*5 param(47),paramj
	character*1 agcor,line(NPA)
	character*11 mark
	integer*4 mfit(NPAP1)
	logical eclcoord

        data param/'     ','   f0','   f1','   f2','  Dec','   RA',
     +       ' pmdc',' pmra','    x','    e','   T0','   Pb','   Om',
     +       ' Omdt','gamma','   DM','   px',' Pbdt',' PPNg','    s',
     +       '    M','   m2','  dth',' xdot',' edot','   x2','   e2',
     +       '  T02','  Pb2','  Om2','   x3','   e3','  T03','  Pb3',
     +       '  Om3',' PMRV',' XOmd',' Xpbd',' om2d','  x2d',
     +       '   f3','   f4','   f5','   f6','   f7','   f8','   f9'/
	data mark/' 123456789*'/

c       Note: indexes in param() now have same numbering as
c             other tempo arrays, such as x(), freq(), and ferr()
c                                                  --DJN 19 Dec 2001 



	if(nbin.eq.9)then
	   param(10)='eps1'
	   param(13)='eps2'
	   param(39)='e1dt'
	   param(40)='e2dt'
	endif
	if (eclcoord) then
	  param(5)='beta'
	  param(6)='lambda'
	  param(7)='pmbeta'
	  param(8)='pmlambda'
	endif
	  

	close(72)
	open(72,file='matrix.tmp',status='unknown',form='unformatted')

	do 5 j=1,nn
5	sig(j)=sqrt(a(j,j))

	write(31,1005) (mod(i,10),i=1,nn)
1005	format(12x,'G  ',500i1)
	write(31,1006) ('-',i=1,nn)
1006	format(12x,'-  ',500a1)

	do 20 j=1,nn
	do 10 k=1,nn
	x=abs(a(j,k)/(sig(j)*sig(k)))
	i=1
	if(x.gt.0.9d0) i=2
	if(x.gt.0.99d0) i=3
	if(x.gt.0.999d0) i=4
	if(x.gt.0.9999d0) i=5
	if(x.gt.0.99999d0) i=6
	if(x.gt.0.999999d0) i=7
	if(x.gt.0.9999999d0) i=8
	if(x.gt.0.99999999d0) i=9
	if(x.gt.0.999999999d0) i=10
	if(x.gt.0.9999999999d0) i=11
10	line(k)=mark(i:i)
	y=abs(gcor(j))
	i=1
	if(y.gt.0.9d0) i=2
	if(y.gt.0.99d0) i=3
	if(y.gt.0.999d0) i=4
	if(y.gt.0.9999d0) i=5
	if(y.gt.0.99999d0) i=6
	if(y.gt.0.999999d0) i=7
	if(y.gt.0.9999999d0) i=8
	if(y.gt.0.99999999d0) i=9
	if(y.gt.0.999999999d0) i=10
	if(y.gt.0.9999999999d0) i=11
	agcor=mark(i:i)
	jj=mfit(j+1)

	if(jj.le.40) then
	  paramj=param(jj)			!One of the listed params
!#### 
!PARAMETERS 50 THROUGH 59 HAVE BEEN SUPERCEDED ... NO LONGER DM derivatitves
!NOT CURRENTLY IN USE
	else if(jj.le.50) then     
          write(*,1017) jj
 1017	  format ('Internal tempo error.  Parameter number ',i5/
     +         'encountered in mxprt.f, but it is not currently'/
     +         'assigned to any parameter.')

!	  write(paramj,1017) jj-40		!DM polynomial coeffs
!1017	  format('  DM',z1)
!#### END OF SUPERCEDED SECTION
	else if(jj.le.NPAR1) then
	  write(paramj,1018) jj-50+2		!Freq derivatives
 1018     format('  f',i2.2)
	else if(jj.lt.NPAR2)then		!Gltiches
	  jj1=jj-(NPAR1+1)
	  jj2=mod(jj1,NGLP)+1
	  jj1=jj1/NGLP+1
	  write(paramj,1030)jj1,jj2
1030	  format(' GL',2i1)
	elseif (jj.lt.NPAR3)then
	  write(paramj,1019) jj-(NPAR2)		!Offsets
1019	  format(' O',i3.3)
	  if(paramj(2:3).eq.'O ') paramj(2:3)=' O'
	elseif (jj.lt.NPAR4) then
	  write(paramj,1050) jj-NPAR3-1
 1050	  format(' FB',i2.2)
	elseif (jj.lt.NPAR5) then
	  write(paramj,1060) jj-NPAR4+1
 1060	  format(' XD',i2.2)
	elseif (jj.lt.NPAR6) then
	  if (2*int(jj/2).ne.jj) then
	    write(paramj,1070) (jj-NPAR5+1)/2
 1070	    format('FJ',i3.3)
	  else
	    write(paramj,1080) (jj-NPAR5)/2
 1080	    format('TJ',i3.3)
	  endif
	else if(jj.lt.NPAR7) then
          if (mod(jj-NPAR6,2).eq.1) then
	    write (paramj,1090) (jj-NPAR6+1)/2
 1090	    format('DX',i3.3)
          else 
	    write (paramj,1091) (jj-NPAR6)/2
 1091	    format('D1',i3.3)
          endif
	elseif (jj.lt.NPAR8) then
	  write(paramj,1093) jj-NPAR7+1
 1093	  format(' ED',i2.2)
	elseif (jj.lt.NPAR9) then
	  write(paramj,1095) jj-NPAR8+1
 1095	  format('OMD',i2.2)
	elseif (jj.lt.NPAR10) then
	  write(paramj,1100) jj-NPAR9		!DM polynomial coeffs
 1100	  format('DM',i3.3)
	elseif (jj.lt.NPAR11) then              ! FD values
	  write(paramj,1101) jj-NPAR10
 1101     format('FD',i1)
	elseif (jj.lt.NPAR12) then              ! More frequency derivatives
	  write(paramj,1102) jj-NPAR11+2
 1102     format(' f',i1)
	else					! XMX parameters
          if (mod(jj-NPAR12,2).eq.1) then
	    write (paramj,1111) (jj-NPAR12+1)/2
 1111	    format('XX',i3.3)
          else 
	    write (paramj,1112) (jj-NPAR12)/2
 1112	    format('XE',i3.3)
          endif
	endif

	write(72) nn,j,paramj,gcor(j),sig(j),
     +    (a(j,k)/(sig(j)*sig(k)),k=1,nn)
        write(31,1020,advance="no") float(j),paramj,agcor
1020	format(f4.0,2x,a5,'|',a1,2x)
        do i = 1, nn
          write (31,fmt="(a1)",advance="no") line(i)
        end do
        write (31,*)  ! end line
20      continue

	write(31,1006) ('-',i=1,nn)
	write(31,1005) (mod(i,10),i=1,nn)

	return
	end
