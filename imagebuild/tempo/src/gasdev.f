c      $Id$
	function gasdev(idum)

C  Returns a normally distributed deviate with zero mean and
C  unit variance, using ran1(idum) as the source of uniform 
C  deviates.

        save iset, gset
	data iset/0/

	if(iset.eq.0) then
1	  v1=2.*ran1(idum)-1.
	  v2=2.*ran1(idum)-1.
	  r=v1**2 + v2**2
	  if(r.ge.1.0) go to 1
	  fac=sqrt(-2.*log(r)/r)
	  gset=v1*fac
	  gasdev=v2*fac
	  iset=1
	else
	  gasdev=gset
	  iset=0
	endif
	return
	end
