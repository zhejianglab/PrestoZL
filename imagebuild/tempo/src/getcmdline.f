	character*640 function getcmdline()

	character*640 a
        character*640 s

	a = ""
        la = 0
	do i = 0, iargc()
          call getarg(i,s)
          ls = len_trim(s)
          if (i.gt.0 .and. la.lt.640) then
            a = a(1:la) // " "
            la = la + 1
          endif
          if (la+ls.gt.640)  then
            a = "echo Command line too long to save"
            goto 100
          endif
          a = a(1:la) // s(1:ls)
          la = la + ls
        enddo

 100	continue
	getcmdline = a
	end
