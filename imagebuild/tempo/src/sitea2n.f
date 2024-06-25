      integer function sitea2n(s)

C     convert first character of ascii string s into a numerical
C     observatory code

      implicit real*8 (a-h,o-z)
      include 'dim.h'
      include 'acom.h'

      character*(*) s
      character*2 ss

      if (len(s).eq.1) then
        ss(1:1)=s(1:1)
        ss(2:2)= " "
      else
        ss(1:2) = s(1:2)
      endif

      call upcase(ss)

      if (ss(2:2).eq.' ') then       ! one-character code

        if(ss(1:1).ge.'0'.and.ss(1:1).le.'9')then
          sitea2n=ichar(ss(1:1))-48
        else if(ss(1:1).ge.'A'.and.ss(1:1).le.'Z')then
          sitea2n = ichar(ss(1:1))-55
        else if (ss(1:1).eq.'@') then
          sitea2n = -1
        else
          sitea2n = -2
        endif

      else                           ! two-character code

        sitea2n = -2
        do i = 1, NOBSMAX
          if (ss.eq.obskey(i)(4:5)) then
            sitea2n = i
            goto 100
          endif
        enddo
100     continue

      endif

      return
      end
