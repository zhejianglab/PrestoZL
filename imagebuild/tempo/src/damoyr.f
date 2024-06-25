c      $Id$
	character*9 function damoyr(mjd)

C  Converts integer*4 MJD to dd-mmm-yy

	implicit real*8 (a-h,o-z)
	integer*2 moday(12)
	character*3 mo,month(12),date*9
	data month/3HJAN,3HFEB,3HMAR,3HAPR,3HMAY,3HJUN,3HJUL,3HAUG,3HSEP,
     +  3HOCT,3HNOV,3HDEC/
	data moday/31,28,31,30,31,30,31,31,30,31,30,31/

	iday=mjd-39126+1
	iyr=1966
	do 10 i=1,133
	nd=365
	if((mod(iyr,4).eq.0 .and. mod(iyr,100).ne.0).or.
     +       mod(iyr,400).eq.0) nd=366
	if(iday.le.nd) go to 20
	iday=iday-nd
10	iyr=iyr+1
20	moday(2)=28
	if(nd.eq.366) moday(2)=29
	do 30 j=1,12
	if(iday.le.moday(j)) go to 40
30	iday=iday-moday(j)
	j=12
40	mo=month(j)
        iyr = mod(iyr,100)
	write(date,1040) iday,mo,iyr
1040	format(i2,'-',a3,'-',i2.2)
	damoyr=date

	return
	end
