c      $Id$
      subroutine bootmc(npts,mode,nz,nboot,nparam,mfit,freq,ferr,ddmch,
     +     buf,npmsav,ksav,resfile2)
C  Implements the "bootstrap method" of Monte Carlo error estimation.
C  See Numerical Recipes, Second Edition, pp 686 ff.

	implicit real*8 (a-h,o-z)
	save
	include 'dim.h'
	real*8 ddmch(*)
        real*8 buf(*)
        integer npmsav(*), ksav(*)
	real*8 xmean(NPA),freq(NPAP1),ferr(NPAP1),a(NBOOTMAX,NPAP1)
	real*8 p0(NPAP1),e0(NPAP1),psq(NPAP1),fctn(NPAP1)
	real*8 fl1(NPAP1),fl2(NPAP1),fh1(NPAP1),fh2(NPAP1)
	integer mfit(NPAP1)
	real*4 ran1
        character*80 resfile2
	logical lw,first
	data lw/.false./,first/.true./,idum/-999/

	first=.false.
	write(*,1000) nboot
1000	format(' Uncertainties by bootstrap Monte Carlo:',
     +    i6,' iterations.')

	nterms=nparam-1
	do 1 j=1,nterms				!Save the fitted params
	k=mfit(j+1)				! and uncertainties
	p0(j)=freq(k)
	e0(j)=ferr(k)
1	psq(j)=0.				!Clear variance accumulator
	fac=sqrt(dfloat(npts))
	x1=0.342*nboot
	x2=0.477*nboot
	xmid=0.5*(nboot+1)
	il1= nint(xmid - x1)
	il2= nint(xmid - x2)
	ih1= nint(xmid + x1)
	ih2= nint(xmid + x2)

	do 50 iter=1,nboot			!Do the bootstrap Monte-Carlo
	sum=0.					!Clear accumulators
	sumwt=0.
	do 10 j=1,nterms
10	xmean(j)=0.

	do 30 i=1,npts
	  ii=npts*ran1(idum) + 1.0		!Randomize the data index
	  call vmemr(ii,fctn,ct,dt,wgt,dn,terr,frq,fut,rfrq,nterms,
     +     buf,npmsav,ksav)
	  do 20 j=1,nterms
20	  xmean(j)=xmean(j)+wgt*fctn(j)
	  sum=sum+wgt*dt
	  sumwt=sumwt+wgt
30	continue

        call fit(npts,mode,chisqr,varfit,xmean,sum,sumwt,nz,wmax,
     +       lw,ddmch,buf,npmsav,ksav,resfile2)
	do 40 j=1,nterms
	k=mfit(j+1)
	x=fac*(freq(k)-p0(j))/e0(j)
	a(iter,j)=x
40	psq(j)=psq(j) + x*x
50	continue

	do 60 j=1,nterms
	k=mfit(j+1)
	freq(k)=p0(j)				!Restore fitted params
	ferr(k)=e0(j)*sqrt(psq(j)/nboot)	!Bootstrapped error estimate
	call sort(nboot,a(1,j))
	fl1(k)=a(il1,j)
	fl2(k)=a(il2,j)
	fh1(k)=a(ih1,j)
60	fh2(k)=a(ih2,j)

	write(72) nboot
	write(72) nterms,(ferr(mfit(j+1))/e0(j),j=1,nterms)
	write(72) nterms,(fl2(mfit(j+1)),j=1,nterms)
	write(72) nterms,(fl1(mfit(j+1)),j=1,nterms)
	write(72) nterms,(fh1(mfit(j+1)),j=1,nterms)
	write(72) nterms,(fh2(mfit(j+1)),j=1,nterms)

	return
	end
