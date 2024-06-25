c      $Id$
c
c  Adjustable arrays:
	integer*8 dnplsoff
	real*8 dnpls(1)
        integer*8 ddmchoff
        real*8 ddmch(1)
        integer*8 bufoff
        real*8 buf(1)
        integer*8 npmsavoff
        integer npmsav(1)
        integer*8 ksavoff
        integer ksav(1)

        common /array/
     +    dnpls,   ddmch,   buf,   npmsav,   ksav,
     +    dnplsoff,ddmchoff,bufoff,npmsavoff,ksavoff

