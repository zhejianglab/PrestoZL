c      $Id$
	subroutine dskrd(jrec,djj,emb,ctemb,embc,xmoon,dpsi,deps)

	implicit real*8 (a-h,o-z)
	real*8 emb(60),embc(60),xmoon(80,6),dpsi(8),deps(8)
	real*8 temb(6),tembc(6),txmoon(48)

	read(44,rec=jrec) djj,temb,ctemb,tembc,txmoon,dpsi,deps
	do 10 i=1,6
	k=10*i-9
	emb(k)=temb(i)
10	embc(k)=tembc(i)

	k=0
	do 20 i=1,8
	do 20 j=1,6
	k=k+1
20	xmoon(i,j)=txmoon(k)

	return
	end
