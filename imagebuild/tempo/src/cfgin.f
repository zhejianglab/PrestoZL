      subroutine cfgin(cfgfile,ut1file,obsyfile,tdbfile)

c     parse tempo configuration file, tyipcally $TEMPO/tempo.cfg
c     djn 30 August 2013

c     most of this was formerly part of the main tempo.f routine

      implicit real*8 (a-h,o-z)

      include 'dim.h'
      include 'acom.h'
      include 'clocks.h'
      include 'eph.h'
      include 'tz.h'

c     input parameters
      character*(*) cfgfile
c     output parameters
      character*(*) ut1file, obsyfile, tdbfile

c     local variables
      character*640 line, label, fname, s
      character*640 key, val
      integer vall, err
      integer i, j, k
      integer iclkbipm, ix

      open(2,file=cfgfile,status='old',err=990)

      clklbl(0)='UNCORR'
      iclkbipm = 6 ! index number of next-added tt(BIPM) entry

      kephem = 0

      do 510 i=1,999
        read(2,1000,end=511)line
 1000   format(a)
        j=1
        call citem(line,640,j,label,k)
        call upcase(label)
        call citem(line,640,j,fname,k)

        if(label(1:8).eq.'OBS_NIST')then
          clkfile(1)=fname(1:k) ! Obs to UTC(NIST)
          clklbl(1)='UTC(NIST)'
        else if(label(1:8).eq.'NIST_UTC')then
          clkfile(2)=fname(1:k) ! UTC(NIST) to UTC(BIPM)
          clklbl(2)='UTC(BIPM)'
        else if(label(1:8).eq.'NIST_PTB')then
          clkfile(3)=fname(1:k) ! UTC(NIST) to UTC(PTB)
          clklbl(3)='PTB'
        else if(label(1:8).eq.'NIST_AT1')then
          clkfile(4)=fname(1:k) ! UTC(NIST) to AT1
          clklbl(4)='AT1'
        else if(label(1:9).eq.'NIST_BIPM')then
          if (label(10:10).eq.' ') then
            clkfile(5)=fname(1:k) ! UTC(NIST) to TT(BIPM)
            clklbl(5)='TT(BIPM)'
          else
            if (iclkbipm>NCLKMAX) then
              print *,"Error, too many clocks listed in tempo.cfg"
              print *,"Maximum number of clocks is NCLKMAX=",NCLKMAX
              stop
            endif
            ix = index(label," ")
            clkfile(iclkbipm)=fname(1:k) ! UTC(NIST) to TT(BIPMxx)
            clklbl(iclkbipm) = 'TT(BIPM' // label(10:ix-1) //')'
            iclkbipm = iclkbipm + 1
          endif
        else if(label(1:3).eq.'UT1')then
          ut1file=fname(1:k)  ! UT1 - UT
        else if(label(1:6).eq.'EPHDIR')then
          ephdir=fname(1:k)   ! BC ephem directory
        else if(label(1:6).eq.'CLKDIR')then
          clkdir=fname(1:k)   ! Clock files directory
        else if(label(1:6).eq.'PARDIR')then
          pardir=fname(1:k)   ! .par files for -z
        else if(label(1:6).eq.'TZDIR')then
          tzdir=fname(1:k)    ! Dir for tz.in
        else if(label(1:5).eq.'OBSYS')then
          obsyfile=fname(1:k) ! Observatory parameters
        else if(label(1:5).eq.'TZTOT')then
          tztfile=fname(1:k)  ! param file for -z
        else if(label(1:7).eq.'TDBFILE')then
          tdbfile=fname(1:k)
        else if(label(1:7).eq.'EPHFILE')then
           kephem=kephem+1     ! BC ephemeris names
           ephnamel(kephem) = 0
           ephbigendian(kephem) = .true.
           do while (k.gt.0)
             if (fname(1:2).eq.'--') then
               call spliteq(fname(3:158),key,val,vall,err)
               call upcase(key)
               if (key(1:3).eq.'EPH') then
                 ephname(kephem) = val(1:vall)
                 ephnamel(kephem) = vall
               else if (key(1:6).eq.'ENDIAN') then
                 call upcase(val)
                 if (val(1:1).eq.'L') ephbigendian(kephem) = .false.
               else
                 print *,"Warning: don't understand tempo.cfg entry:"
                 print *,line
               endif
             else
               ephfile(kephem) = fname(1:k)
             endif
             call citem(line,640,j,fname,k)
           enddo
           if (ephnamel(kephem).eq.0) then
             ! historical default: ephemeris name (e.g. "DE405") is first five chars of file name
             ! historical default: ephemeris in bigendian format
             ephname(kephem) = ephfile(kephem)(1:5)
             ephnamel(kephem) = 5
           endif
        else if (label(1:6).eq.'TZSITE') then
          tzsitedef = fname(1:1)
        else if (label(1:9).eq.'TDBFMT') then
          call upcase(fname)
          if (fname(1:4).eq.'IF99') tdbif99fmt = .true.
        else if (label(1:1).ne."#" .and. label(1:1).ne." ") then
          goto 991
        endif
 510  continue ! end of main loop

 511  continue ! jump to here when loop is done

      close(2)

      goto 999

 990  continue
      print *,"Cannot open configuration file at ", cfgfile
      stop

 991  continue
      print *,"Error reading configuration file ",cfgfile
      print *,"unrecognised label ", label
      stop

 992  continue   
      print *,"Error reading configuration file ",cfgfile
      print *,"Too many ephemeris files"
      stop

 999  return
      end     


