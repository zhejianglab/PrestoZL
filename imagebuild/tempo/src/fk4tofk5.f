c      $Id$
      subroutine FK4TOFK5(RA,DEC,PMRA,PMDEC,PX,V)
C
C     Converts stellar positions from FK4 (B1950.0) coordinates to
C     FK5 (J2000.0) coordinates, using method from (e.g.) Astronomical
C     Almanac, 1987 edition, p. B42
C
C     WARNING:  Many parts of the program below are normally stored in
C               "Comment" lines, and are not executed, even though it
C               looks like they should at first glance.
C
C     input/output parameters:
C     RA,DEC     -- Position in radians
C     PMRA,PMDEC -- Proper Motion in seconds of arc/tropical century
C     PX         -- Paralax in seconds of arc
C     V          -- Radial velocity in km/sec
C
      implicit real*8 (A-H,O-Z)
	parameter (TWOPI=6.28318530717958648d0)
      real*8 M(6,6)
C	real*8 A(3),ADOT(3)
      real*8 R0(3),R0DOT(3),R1(3),R1DOT(3),R(3),RDOT(3)
      equivalence (R0,R1)
      equivalence (R0DOT,R1DOT)

c     astronomical almanac matrix M is listed below
c
      data M/ 0.9999256782,      0.0111820610,      0.0048579479,
     +          -0.000551,          0.238514,         -0.435623, 
     +       -0.0111820611,      0.9999374784,     -0.0000271474,
     +          -0.238565,         -0.002667,          0.012254,
     +       -0.0048579477,     -0.0000271765,      0.9999881997,
     +           0.435739,         -0.008541,          0.002117,
     +        0.00000242395018,  0.00000002710663,  0.00000001177656,
     +           0.99994704,        0.01118251,        0.0048567,
     +       -0.00000002710663,  0.00000242397878, -0.00000000006582,
     +          -0.01118251,        0.99995883,       -0.00002718,
     +       -0.00000001177656, -0.00000000006587,  0.00000242410173,
     +          -0.004857567,      -0.00002718,        1.00000956/
c
c     standish's M matrix is listed below
c
c      data M/0.9999256782190734, 0.0111820581287209, 0.0048579435553111      
c     +, -0.0005266390034373, 0.2371579690136167,-0.4373690764163668 
c     +,     -0.0111820580744435, 0.9999374784651157,-0.0000271730844707
c     +, -0.2371432945996690,-0.0026518108275197,-0.0011512525743752
c     +,     -0.0048579436802253,-0.0000271507426487, 0.9999881997539576
c     +,  0.4373353030376622, 0.0048893075298284, 0.0021251718008283
c     +,      0.0000024238984065, 0.0000000271047938, 0.0000000117786069
c     +,  0.9999256794956877, 0.0111814832391717, 0.0048590037723143
c     +,     -0.0000000271047938, 0.0000024239270237,-0.0000000000658629
c     +, -0.0111814832204662, 0.9999374848933135,-0.0000271702937440
c     +,     -0.0000000117786070,-0.0000000000658443, 0.0000024240499480
c     +, -0.0048590038153592,-0.0000271625947142, 0.9999881946023742/
c
c     standish's M4 matrix is listed below
c
c      data M/ .9999257079523629, .0111789381264276, .0048590038414544,
c     +0,0,0, -.0111789381377700, .9999375133499888,-.0000271579262585,
c     +0,0,0, -.0048590038153592,-.0000271625947142, .9999881946023742,
c     +0,0,0,                                                    0,0,0,
c     + .9999257079523629, .0111789381264276, .0048590038414544, 0,0,0,
c     +-.0111789381377700, .9999375133499888,-.0000271579262585, 0,0,0,
c     +-.0048590038153592,-.0000271625947142, .9999881946023742/       
c
C      data A/-1.62557E-6,-0.31919E-6,-0.13843E-6/
C      data ADOT/1.244E-3,-1.579E-3,-0.660E-3/
C
C     calculate rectangular coordinates in 1950 epoch
C
      R0(1) = dcos(ra)*dcos(dec)
      R0(2) = dsin(ra)*dcos(dec)
      R0(3) = dsin(dec)
      R0DOT(1) = -PMRA*R0(2)-PMDEC*dcos(RA)*dsin(DEC)+21.095*V*PX*R0(1)
      R0DOT(2) =  PMRA*R0(1)-PMDEC*dsin(RA)*dsin(DEC)+21.095*V*PX*R0(2)
      R0DOT(3) =  PMDEC*dcos(DEC)+21.095*V*PX*R0(3)
C
C     Pxy = dot product of x & y
C
C     Remove the effects of the E-terms of aberration
C     (Take the "C" off the front of the next five lines and remove the
C      "equivalence" statements above to include this routine)
C
C      PR0A = DOT(R0,A)
C      PR0ADOT = DOT(R0,ADOT)
C      DO 10 I = 1, 3
C        R1(I) = R0(I) - A(I) + PR0A * R0(I)
C 10     R1DOT(I) = R0DOT(I) - ADOT(I) + PR0ADOT * R0(I)
c
c     Find the rectangular coordinates in the 2000 epoch
c
      DO 20 I = 1, 3
        R(I) = 0
        RDOT(I) = 0
        DO 20 J = 1, 3
          R(I) = R(I) + M(I,J)*R1(J) + M(I,J+3)*R1DOT(J)
20        RDOT(I) = RDOT(I) + M(I+3,J)*R1(J) + M(I+3,J+3)*R1DOT(J)
C
C     Convert back to spherical coordinates
C
      XY = R(1)*R(1) + R(2)*R(2)
      ABSR = DSQRT(XY+R(3)*R(3))
	cd=sqrt(xy)
	dec=atan2(r(3),cd)
	ra=atan2(r(2),r(1))
      if (RA .LT. 0D0) RA = RA + TWOPI
      PMRA = (R(1)*RDOT(2)-R(2)*RDOT(1))/XY
      PMDEC = (RDOT(3)*XY-R(3)*(R(1)*RDOT(1)+R(2)*RDOT(2)))/
     +                                       (ABSR*ABSR*DSQRT(XY))
      IF (PX .NE. 0) V = DOT(R,RDOT)/(21.095*PX*ABSR)
      PX = PX/ABSR
C
      RETURN
      END
