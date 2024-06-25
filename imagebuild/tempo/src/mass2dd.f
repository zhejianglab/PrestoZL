c      $Id$
c
	subroutine mass2dd(am,am2,x,ecc,an,arr,ar,xk,si,gamma,pbdot)

c       given system masses m, m2 and keplerian parameters x, ecc, an
c       calculate values of arr, ar, si, gamma, pbdot under GR

	implicit none
	real*8 TWOPI,SUNMASS,ARRTOL
	parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
	parameter (ARRTOL=1.d-10)

	real*8 am,am2,x,ecc,an,arr,ar,xk,si,gamma,pbdot
	real*8 m, m2, m1, arr0, arrold
	m=am*SUNMASS
	m2=am2*SUNMASS
	m1=m-m2
	arr0=(m/(an**2)) ** (1.d0/3)
	arr = arr0
 5	continue
	arrold = arr
	arr=arr0*(1+(m1*m2/m**2 - 9)*0.5d0*m/arr) ** (2.d0/3)
	if (abs((arr-arrold)/arr).gt.ARRTOL) goto 5
	arr=arr0*(1+(m1*m2/m**2 - 9)*0.5d0*m/arr) ** (2.d0/3)
	ar=arr*m2/m
	si=x/ar
C        xk=3*m/(arr*(1-ecc**2))                 
C IHS based on 060327 changing to non-tw89 defn for omdot.
        xk=3*m/(arr0*(1-ecc**2))                 
C IHS 100420 changing gamma usage as well as per Norbert's advice
C	gamma=ecc*m2*(m1+2*m2)/(an*arr*m)
	gamma=ecc*m2*(m1+2*m2)/(an*arr0*m)
	pbdot=-(96*twopi/5) * an**(5.d0/3) * (1-ecc**2)**(-3.5d0) * 
     +    (1+ (73.d0/24)*ecc**2 + (37.d0/96)*ecc**4) * 
     +    m1*m2*m**(-1.d0/3) 
	return
	end
