c      $Id$
c
      subroutine tparin(nostop,lw,lpth,nparmax,nptsmax,
     +     version,npulsefile,infile,path,resfile1,hlpfile,parfile)

c     tparin -- Tempo PARamater INput
c     parses tempo command-line parameters and sets variables appropriately

      implicit real*8 (a-h,o-z)

      include 'dim.h'
      include 'acom.h'
      include 'vcom.h'
      include 'tz.h'
      include 'version.h'

c     following variables are set in this routine:

      logical nostop
      logical lw
      integer lpth
      integer nparmax, nptsmax
      real*8 version
      character*160 npulsefile
      character*640 infile
      character*640 path
      character*80 resfile1
      character*80 hlpfile
      character*160 parfile
      character*160 key
      character*160 val
      integer err

c     following variables are also set, but are in common blocks
c     and hence not formally declared here      
c      logical xitoa 
c      character*1 obsflag
c      logical gro
c      logical jumpout
c      logical lresid1
c      logical npulsein
c      logical npulseout
c      logical oldpar
c      logical psrframe
c      logical useglsfit
 

c     variables used internally in this routine:

      integer narg, iarg, ipar
      integer iargc
      integer i, ii
      character*160 s, s2

      integer vall
      character*1 vallast
      real tmp


      character*160 getparm  ! external function


c     default values of parameters


      gro = .false.
      jumpout = .false.
      lresid1 = .false.
      ldesign = .false.
      ldcov = .false.
      nostop = .false.
      npulsein = .false.
      npulseout = .false.
      phisunout = .false.
      dopplerout = .false.
      dmxnout = .false.
      nonewdmx = .false.
      oldpar = .false.
      psrframe = .false.
      lw = .true.
      ssdmflag = 1
      tz = .false.
      autotz = .false.
      xitoa = .false.
      quiet = .false.
      polystdout = .false.
      jumpbarycenter = .false.
      useglsfit = .false.
      nparmax = NPARDEF
      nptsmax = NPTSDEF
      parfile = 'def'
      polycofile = 'polyco.dat'

      infile = ''	

      call getenv ('TEMPO',path)
      lpth = index(path,' ')-1
      if (lpth.eq.0) then
        print *, "You must set the $TEMPO environment variable in order to run tempo."
        stop
      endif
      hlpfile = path(1:lpth)//'/tempo.hlp'	


      narg = iargc()
      iarg = 0
      ipar = 0
      if (narg.lt.1) goto 9999  ! error

      do while (iarg.lt.narg)

        iarg = iarg + 1
        call getarg (iarg,s)
        ii = index(s,' ') - 1

        if (s(1:1).eq.'-') then

          i = 2
          do while (i.le.ii)
            if (s(i:i).eq.'a') then
              phisunout = .true.
            else if (s(i:i).eq.'b') then
              nonewdmx = .true.
            else if (s(i:i).eq.'c') then
              nostop = .true.
            else if (s(i:i).eq.'C') then
              ldcov = .true.
            else if (s(i:i).eq.'d') then
              path = getparm(s,i,ii,iarg,narg)
              if (path.eq.'') goto 9999 ! error
              lpth = index(path,' ')-1
              hlpfile = path(1:lpth)//'/tempo.hlp'
	    else if (s(i:i).eq.'D') then
	      ldesign = .true.
	      open(37,file='design.tmp',form='unformatted',
     +                                       status='unknown')
            else if (s(i:i).eq.'e') then
              dopplerout = .true.
            else if (s(i:i).eq.'f') then
              parfile = getparm(s,i,ii,iarg,narg)
              if (parfile.eq.'') goto 9999 ! error
            else if (s(i:i).eq.'g') then
              gro = .true.
              obsflag = 'P'
              if (s(i+1:i+1).ne.' ') obsflag = s(i+1:i+1)
              i = 9999
            else if (s(i:i).eq.'G') then
              useglsfit = .true.
            else if (s(i:i).eq.'h') then
              goto 9999         ! error/help
            else if (s(i:i).eq.'j') then
              jumpout = .true.
            else if (s(i:i).eq.'J') then
              jumpbarycenter = .true.
            else if (s(i:i).eq.'k') then
              dmxnout = .true.
            else if (s(i:i).eq.'l') then
              s2 = getparm(s,i,ii,iarg,narg)
              read (s2,*) nparmax
            else if (s(i:i).eq.'m') then
              s2 = getparm(s,i,ii,iarg,narg)
              read (s2,*) nptsmax
            else if (s(i:i).eq.'n') then
              i = i + 1
              if (s(i:i).eq.'i') then
                npulsein = .true.
                npulsefile = getparm(s,i,ii,iarg,narg)
              elseif (s(i:i).eq.'o') then
                npulseout = .true.
                npulsefile = getparm(s,i,ii,iarg,narg)
              else
                goto 9999       ! error
              endif
            else if (s(i:i).eq.'o') then
              oldpar = .true.
            else if (s(i:i).eq.'p') then
              psrframe = .true.
            else if (s(i:i).eq.'r') then
              lresid1 = .true.
              open(30,file=resfile1,form='unformatted',status='unknown')
            else if (s(i:i).eq.'s') then
              s2 = getparm(s,i,ii,iarg,narg)
              read (s2,*) ssdmflag
            else if (s(i:i).eq.'v') then
              write(*,1001) version,trim(VERSIONID)
 1001         format (' Tempo v ',f6.3,' (',a,')')
              stop
            else if (s(i:i).eq.'w') then
              lw = .false.
            else if (s(i:i).eq.'x') then
              xitoa = .true.
              open(45,file='itoa.out',status='unknown')
            else if (s(i:i).eq.'z') then
              tz = .true.
            else if (s(i:i).eq.'Z') then
              if (.not.autotz) then
                tz = .true.
                quiet = .true.
                autotz = .true.
                ! set default autotz variables
                name(1) = ' '
                nco(1) = 15
                nsp(1) = 60
                mxha(1) = 1
                tzof(1) = 1420 ! arbitrary default frequency
                tzsite = ''
                tzmjdstart = -1
                ! nsite = sitea2n(tzsite)
              endif
              s2 = getparm(s,i,ii,iarg,narg)
              call spliteq(s2,key,val,vall,err)
              if (err.eq.1) then
                write (*,*) "Error: -Z parameter, no value: ",s(i+1:ii)
                write (*,*)
                goto 9999
              endif
              if (err.eq.0) then
                if (key(1:3).eq."PSR".or.key(1:6).eq."PULSAR") then
                  name(1) = val(1:12)
                else if (key(1:3).eq."OUT") then
                  polycofile = val
                else if (key(1:6).eq."NCOEFF") then
                  read(val,*) nco(1)
                else if (key(1:4).eq."SPAN") then
                  vallast = val(vall:vall)
                  call upcase(vallast)
                  if (vallast.eq."H") then
                       read(val(1:vall-1),*) tmp    
                       nsp(1) = tmp*60              ! convert hr to min
                  else if (vallast.eq."M") then
                       read(val(1:vall-1),*) tmp    
                       nsp(1) = tmp                 ! minutes: no conversion
                  else if (vallast.eq."S") then
                       read(val(1:vall-1),*) tmp    ! conver min to sec
                       nsp(1) = tmp/60              ! convert sec to min
                  else       !    no trailing letter, interpret as seconds
                       read(val,*) tmp   
                       nsp(1) = tmp/60              ! convert sec to min
                  endif
                else if (key(1:4).eq."TOBS") then
                  vallast = val(vall:vall)
                  call upcase(vallast)
                  if (vallast.eq."M") then
                       read(val(1:vall-1),*) tmp
                       mxha(1) = tmp/60.              ! convert min to hr
                  else if (vallast.eq."S") then
                       read(val(1:vall-1),*) tmp     
                       mxha(1) = tmp/3600.            ! convert sec to hr
                  else if (vallast.eq."H") then
                       read(val(1:vall-1),*) mxha(1)  ! hr: no conversion
                  else       !    no trailing letter, interpret as hr      
                       read(val,*) mxha(1)            ! hr: no conversion
                  endif
                else if (key(1:4).eq."FREQ") then
                  read (val,*) tzof(1)
                else if (key(1:4).eq."SITE".or.key(1:3).eq."OBS") then
                  tzsite = val
                  ! nsite = sitea2n(tzsite)
                else if (key(1:5).eq."START") then
                  read (val,*) tzmjdstart
                else
                  write (*,*) "Error: incorrect -Z key: ",key
                  write (*,*)
                  goto 9999
                endif
              endif
            else
              goto 9999         ! error
            endif
              
            i = i + 1
          enddo

        else

          ipar = ipar + 1  
          if (ipar.eq.1) then
            infile = s
          else
            goto 9999           ! error
          endif
          
        endif
        
      enddo

      if (ipar.lt.1 .and. .not. tz) goto 9999   ! error -- not enough parameter

      return

 9999 call system ('more '//hlpfile)
      stop
      end

C----------------------------------------------------------------

        character*160 function getparm(s,i,ii,iarg,narg)

c       Extract the parameter associated with a command-line switch.
c
c       Many command-line switches have one and only one parameter.
c       These can be entered either as, e.g., "-xAAAA" or as "-x AAAA",
c       where "x" is the switch and "AAAA" is the parameter.  The switch
c       may be bundled with other switches, as long as those switches
c       do not require parameters, e.g., "-abcxAAAA".
c
c       s:  input to this routine is command-ilne parameter string s
c       i:  index within s pionting to the location of the flag
c       ii: length of s
c       iarg: argument number on command line
c       narg: number of arguments on command line
c
c       strategy:
c       if there is anything in s(i+1:ii), it becomes the return value,
c       otherwise use the next command-line-parameter as return value.

        implicit none

        character*160 s
        integer i, ii, iarg, narg

        integer i0

        getparm = "" ! default (indicates error if this is read upon return)

        i0 = i+1
        if (i.eq.ii) then
          iarg = iarg + 1
          if (iarg.gt.narg) goto 9000
          call getarg(iarg,s)
          i0 = 1
          ii = index(s,' ')
        endif
        
        i = 9999


        getparm = s(i0:ii)

 9000   continue
        return
	end

C----------------------------------------------------------------

	subroutine spliteq(s,key,val,vall,err)

	character*160 s
        character*160 key
        character*160 val
        integer vall
        integer err
        integer idx 
        integer j

c       Examine string s(i+1:ii) for a structure of the form "xxx=yyy".
c       If such a structure is present, return xxx as key and yyy as val.
c       Return length of val as integer vall (length is defined as
c         position of leftmost non-space character)
c       Return values for err:
c         0 = key and val returned
c         1 = no equals sign detected, return s(i+1:ii) as key; val empty
c         2 = empty string
c       Key and Val are explicitly space-padded.  This isn't technically
c       necessary, at least in sun fortran, but it seems worth doing.


        idx = index(s,'=')
        if (idx.gt.0) then
          err = 0
          key = s(1:idx-1)
          do j = idx, len(key)
            key(j:j) = " "
          enddo
          val = s(idx+1:len(s))
          do j = len(s)-idx+1, len(val)
            val(j:j) = " "
          enddo
        else
          err = 1
          key = s(1:len(s))
          if (len(key).gt.len(s)) then   ! silly but may be helpful some day
	    l = len(s)                   ! silence compiler warning if len(s)>len(key)
            do j = l+1, len(key)
              key(j:j) = " "
            enddo
            val = ""
          endif
        endif

c       find length of vall (position of last non-space character)
c       robust (if not efficient) to do this independently of previous code
        vall=160  ! default: maximum length
        do i = 1, len(val)
          if (val(i:i).eq.' ') then
            vall=i-1
            goto 500
          endif
        enddo
  500   continue
        
        call upcase(key)
        
        return
        end            
