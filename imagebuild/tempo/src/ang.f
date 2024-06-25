c      $Id$
	real*8 function ang (i,f)

c  converts input real*8 no. to radians or to fraction of 2pi
c   for i=1,2, input=ddmmss.s? for i=3,4, input=hhmmss.s
c   for i=1,3, output in radians? i=2,4, output in fraction of 2pi

	implicit real*8 (a-h,o-z)
	parameter (TWOPI=6.28318530717958648D0)
	ia=f/1.d4
	ib=(f-ia*1.d4)/1.d2
	g=f-ia*1.d4-ib*1.d2
	ang=(ia+(ib+g/6.d1)/6.d1)/36.d1
	if (i.gt.2) ang=ang*15.d0
	if (i.eq.1.or.i.eq.3) ang=ang*TWOPI

	return
	end
