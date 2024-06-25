c      $Id$
	subroutine matinv (array,norder,det)

C  Inverts a square symmetric matrix array(norder,norder).
C  Result is returned in the same array, determinant in det.

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	dimension array(NPA,NPA)
	integer ik(NPA),jk(NPA)

	det=1.
	do 100 k=1,norder
	amax=0.
21	do 30 i=k,norder
	do 30 j=k,norder
	if(dabs(amax).gt.dabs(array(i,j)))go to 30
	amax=array(i,j)
	ik(k)=i
	jk(k)=j
30	continue

C  Interchange rows and columns to put amax in array(k,k)

	if(amax.eq.0.) go to 900
	i=ik(k)
	if(i-k) 21,51,43
43	do 50 j=1,norder
	save=array(k,j)
	array(k,j)=array(i,j)
50	array(i,j)=-save
51	j=jk(k)
	if(j-k)21,61,53
53	do 60 i=1,norder
	save=array(i,k)
	array(i,k)=array(i,j)
60	array(i,j)=-save

C  Accumulate elements of inverse matrix

61	do 70 i=1,norder
	if(i.eq.k) go to 70
	array(i,k)=-array(i,k)/amax
70	continue
	do 80 i=1,norder
	if(i.eq.k) go to 80
	do 79 j = 1, norder
	if (j .eq. k) go to 79
	array(i,j)=array(i,j)+array(i,k)*array(k,j)
79	continue
80	continue
	do 90 j=1,norder
	if(j.eq.k) go to 90
	array(k,j)=array(k,j)/amax
90	continue
	array(k,k)=1./amax
100	det=det*amax

C  Restore normal ordering of matrix.

	do 130 l=1,norder
	k=norder-l+1
	j=ik(k)
	if(j.le.k) go to 111
	do 110 i=1,norder
	save=array(i,k)
	array(i,k)=-array(i,j)
110	array(i,j)=save
111	i=jk(k)
	if(i.le.k) go to 130
	do 120 j=1,norder
	save=array(k,j)
	array(k,j)=-array(i,j)
120	array(i,j)=save
130	continue
	return

900	det=0.
	return
	end
