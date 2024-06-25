c      $Id$
	subroutine tzinit(alng,sitelng,num)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'tz.h'
        include 'acom.h'

        real*8 alng(NOBSMAX)
	character path*96
	character line*80,item*40
	character*5 obskey0
	real*8 maxhadef

	integer sitea2n ! external function


	kd=index(tzdir,' ')-1
	path=tzdir(1:kd)//tzfile
	open(11,file=path,status='old',err=90)
	go to 92
 90	write(*,'(''Failed to open '',a)')path
	stop

 92	continue
	lupolyco=13
	open(lupolyco,file=polycofile,status='unknown')

c Read default values of asite,maxha,nspan,ncoeff & freq from 1st line of tz.in
c Free format, but must be in order
	read(11,'(a)')line
	jn=1
	ll=80
	call citem(line,ll,jn,item,jl)                  ! asite
	tzsite=item
	nsite = sitea2n(tzsite)
        if (nsite.eq.-2) then
	   write(*,'(''Invalid site in tz.in: '',a)')item(1:jl)
	   stop
	endif

	call citem(line,ll,jn,item,jl) ! maxha
	read(item,*)maxhadef
	
	call citem(line,ll,jn,item,jl) ! nspan
	read(item,*)nspdef
	
	call citem(line,ll,jn,item,jl) ! ncoeff
	read(item,*)ncodef

	call citem(line,ll,jn,item,jl) ! freq
	if(jl.eq.0)then
	   freqdef=-1.		! No freq entered, assume tzref freq
	else 
	   read(item,*)freqdef
	endif
	
C Read pulsar list and nspan, ncoeff, maxha and freq overrides
	read(11,*)
	read(11,*)
        if (.not.quiet) then
  	  write(*,'(''TZ source list for site = '',a/)')tzsite
	  write(*,'(''    PSR        Nspan  Ncoeffs  Maxha    Freq'')')
	  write(*,'(''----------------------------------------------'')')
        endif
	
	do 310 i=1,ntzmax

	  read(11,'(a)',end=320)line

	  jn=1
	  name(i)=' '
	  call citem(line,ll,jn,item,jl)
	  if(jl.gt.12)jl=12
	  name(i)(1:jl)=item(1:jl)

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     nsp(i)=nspdef
	     nco(i)=ncodef
	     mxha(i)=maxhadef
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)nsp(i)
	  endif

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     nco(i)=ncodef
	     mxha(i)=maxhadef
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)nco(i)
	  endif

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     mxha(i)=maxhadef
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)mxha(i)
	  endif

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)tzof(i)
	  endif

 312      continue
          if (.not.quiet) then
    	    if(tzof(i).gt.0.)then
              write(*,'(1x,a,2i8,x,f7.2,1x,f12.5)') name(i),nsp(i),
     +           nco(i), mxha(i),tzof(i)
            else
              write(*,'(1x,a,3i8,2x,''from tzref'')') 
     +           name(i),nsp(i),nco(i)
            endif
          endif
 310	continue
	print*,' Sorry, ',ntzmax,' pulsar limit'
 320	num=i-1
	close(11)

        if(nsite.le.0) then    ! geocentric, barycentric
          sitelng = 0.
        else	
          sitelng=alng(nsite)
        endif

	return
	end

