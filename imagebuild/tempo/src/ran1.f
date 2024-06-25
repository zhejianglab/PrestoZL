c      $Id$
	real*4 function ran1(idum)

C  Returns a random deviate between 0.0 and 1.0.  Set IDUM to
C  any negative value to initialize or reinitialize the sequence.

	save r
        save iff, ix1, ix2, ix3
	dimension r(97)
	parameter (m1=259100,ia1=7141,ic1=54773)
	parameter (m2=134456,ia2=8121,ic2=28411)
	parameter (m3=243000,ia3=4561,ic3=51349)
	data iff/0/
	data ix1/0/
	data ix2/0/
	data ix3/0/

	rm1=1./m1
	rm2=1./m2
	
	if(idum.lt.0.or.iff.eq.0) then
	  iff=1
	  ix1=mod(ic1-idum,m1)
	  ix1=mod(ia1*ix1+ic1,m1)
	  ix2=mod(ix1,m2)
	  ix1=mod(ia1*ix1+ic1,m1)
	  ix3=mod(ix1,m3)

	  do 11 j=1,97
	    ix1=mod(ia1*ix1+ic1,m1)
	    ix2=mod(ia2*ix2+ic2,m2)
11	    r(j)=(ix1+ix2*rm2)*rm1
	  idum=1

	endif

	ix1=mod(ia1*ix1+ic1,m1)
	ix2=mod(ia2*ix2+ic2,m2)
	ix3=mod(ia3*ix3+ic3,m3)
	j=1+(97*ix3)/m3
	if(j.gt.97.or.j.lt.1) stop 'ran1'
	ran1=r(j)
	r(j)=(ix1+ix2*rm2)*rm1
	return
	end
