c      $Id$
      subroutine getflags(s,ll,j1)

      implicit none

      character*640 s
      integer ll, j1
      character*640 s1, s2
      integer s1len, s2len
      integer j1x

      include 'flagcom.h'

      nflag = 0

      do while (j1.lt.ll)
        j1x = j1
        call citem(s,ll,j1,s1,s1len)
        if (s1len.eq.0) exit
        if (s1(1:1).ne.'-') goto 900
        call citem(s,ll,j1,s2,s2len)
        if (s2len.eq.0) goto 900
        nflag = nflag + 1
        flagname(nflag) = s1(2:80)
        flagvalue(nflag)= s2
      end do

      goto 999

 900  continue
      print *,'Tempo getflags: Error parsing flags in the following'
      print *,'  TOA line starting at position ',j1x
      print *,'Full line:'
      print *,s
      print *,'Trouble parsing starting here:'
      print *,s(j1x:)
      stop

 999  continue
      return
      end
        
        
  
               
