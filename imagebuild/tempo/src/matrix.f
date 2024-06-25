c      $Id$
	program matrix
	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*8 a(NPA,NPA),sig(NPA),gcor(NPA),fac(NPA)
	character*5 param(NPA)
	character fmt*32,fmt2*32,fmt3*32,fmt4*32
        character arg*32

	narg=iargc()
	arg=' '
	if(narg.ge.1) then
	  call getarg(1,arg)
          if(narg.gt.1.or.(arg(1:2).lt.'-0'.or.arg(1:2).gt.'-9'))
     +      goto 997
	endif
	ndec=2
	if(arg(1:1).eq.'-') read(arg(2:3),*,err=997) ndec

	open(72,file='matrix.tmp',status='old',form='unformatted',err=998)
	read(72) nn
	rewind 72

c          "XXXXXXXXX" in the following line gets replaced by 
c          something like 037f10.04, where 037 is the number of
c          parameters and 10.04 is the field width and number of
c          to be printed
	fmt='(f4.0,2x,a5,1h|,1x,XXXXXXXXX)'
	write(fmt(20:28),1001) nn,4+ndec,ndec
1001    format(i3.3,'f',i2.2,'.',i2.2)

	do 10 n=1,nn
10	read(72) mm,j,param(j),gcor(j),sig(j),(a(j,k),k=1,mm)

	fmt2='(12x,'//fmt(20:22)//'a'//fmt(24:25)//')'
	write(*,fmt2) (param(j),j=1,nn)
        fmt4='(12x,XXXX(1h-))'
        write(fmt4(6:9),fmt='(i4.4)') nn*(ndec+4)+2
        write(*,fmt4)

	do 20 j=1,nn
20	write(*,fmt) float(j),param(j),(a(j,k),k=1,j)

	fmt3='(/7x,6hgcor:,'//fmt(20:)
	write(*,fmt3) (gcor(j),j=1,nn)
        write(*,fmt4)
	write(*,fmt2) (param(j),j=1,nn)
        write(*,fmt4)

	read(72,end=999) nboot
	write(*,1020) nboot
1020	format(/' Uncertainty factors from bootstrap method (',
     +    i4,' iterations):')

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,6h rms:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,6h  -2:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,6h  -1:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,6h  +1:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,6h  +2:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)
	go to 999

997     continue
	    print*,'Usage: matrix [-n]'
	    print*,
     +       '       n = number of decimal places [default 2]'
	    go to 999

998     continue
   	print*,'Cannot open matrix.tmp'

999	end
