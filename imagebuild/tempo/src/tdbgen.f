c      $Id$
      program tdbgen

c     This program uses the "tdb1ns" routine by Fairhead et al. to generate
c     sets of Chebyshev polynomial coefficients to be used for calculating
c     TDB-TDT in Tempo.

c     Based on convert.f and dconvert.f, older programs which generated
c     tempo support files.

c     DJN  19 August 1997

c     Modifications:
c     DJN  12 August 1998   added byte-swapping for little endian systems

      logical bigendian

      integer NN, NRECL
      parameter (NN=20)         ! #coeffs calculated (not all are retained)


      real*8 c(NN), f(NN), jd
      real*8 pi
      data pi/3.141592653589793d0/

      integer dt, ncf
      real*8 bma, bpa
      real*8 fac
      real*8 y
      real*8 sum


      integer i, j, k
      real*8 d1, d2
      character*160 s
      character*160 outfile

      real*8 tmp1, tmp2
      integer tmp3, tmp4

      integer iargc              ! built-in function

c     the definition of "recl" in file open statements in system-dependent
c     nrecl=4 if recl is in bytes (sun compiler; gcc compiler)
c     nrecl=1 if recl is in 4-byte words (vax compiler; intel compiler)
c     the following (fortran90) should work on all modern systems
      real*4 tmp
      INQUIRE (IOLENGTH=nrecl) tmp



      if (iargc().eq.3) then
        call getarg(1,s)
        read (s,*) d1
        call getarg(2,s)
        read (s,*) d2
        call getarg(3,outfile)
      else
        write (*,*) 'Use:  dconvert d1 d2 outputfile'
        write (*,*) 'd1 and d2 are start, end date in years'
        stop
      endif

c
c	step size (days), number of coefficients per interval
c	used to calculate CTATV (aka TDT-TDB).
c	this should really be a user input.
c	16 days, 8 coefficients gives < 5 ns accuracy
c
      dt = 16                   ! days per interval
      ncf = 8                   ! coefficients retained within an interval

c       convert to JD if Years were entered
      if (d1 .LT. 3000D0) d1 = 365.25*d1 + 1721045
      if (d2 .LT. 3000D0) d2 = 365.25*d2 + 1721045
c       neaten them up a bit
      d1 = dint(d1)
      d2 = d1 + dint((d2-d1)/dt)*dt

      open (8,file=outfile,access='DIRECT',
     +     status='NEW',recl=2*NCF*nrecl)

      if (bigendian()) then
        write (8,rec=1) d1, d2, dt, ncf, (0,i=7,2*NCF) ! pad out to full recl
      else
        tmp1 = d1
        tmp2 = d2
        tmp3 = dt
        tmp4 = ncf
        call dbyterev(tmp1,1)
        call dbyterev(tmp2,1)
        call byterev(tmp3,1)
        call byterev(tmp4,1)
        write (8,rec=1) tmp1,tmp2,tmp3,tmp4, (0,i=7,2*NCF) ! pad out 
      endif
        

c     calculate coefficients from Fairhead et al routine.
c     for information on numerical routines & Chebyshev polynomials, see 
c     Numerical Recipes, sec. 5.6;  note that our "c(1)" is half their "c(1)".
c     Note that NN=20 coefficients are calculated but only NCF are retained 
c     in the file.


      bma = 0.5d0 * dt
      fac = 2.d0 / nn
      do 300 i = 0, (d2-d1)/dt-1
        jd = d1 + i*dt
        bpa = jd + bma
        do 310 j = 1, nn
          y = cos(pi*(j-0.5d0)/nn)*bma
          call tdbtrans (bpa,y,f(j))
 310    continue
        do 320 j = 1, nn
          sum = 0.d0
          do 325 k = 1, nn
            sum = sum + f(k)*cos((pi*(j-1))*((k-0.5d0)/nn))
 325      continue
          c(j) = fac * sum
 320    continue
        c(1) = 0.5d0 * c(1) 
        if (.not. bigendian()) call dbyterev(c,ncf)
        write (8,rec=i+2) (c(j),j=1,ncf) ! i+2 because(a)i=0..(b)rec1 is hdr
          
 300  continue

      close (8)

      end
