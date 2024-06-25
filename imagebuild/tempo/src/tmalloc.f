c      $Id$
c
	subroutine tmalloc(nptsmax,nbuf)

        implicit none

	integer nptsmax, nbuf

	include 'array.h'

        call mallocxd(dnpls,nptsmax,8,dnplsoff)
        call mallocxd(ddmch,nptsmax,8,ddmchoff)
        call mallocxi(ksav,nptsmax,4,ksavoff)
        call mallocxi(npmsav,nptsmax,4,npmsavoff)
        call mallocxd(buf,nbuf,8,bufoff)

        return
        end
	

C----------------------------------------------------------------

	subroutine tfree()

        implicit none
        include 'array.h'

        call freexd(dnpls)
        call freexd(ddmch)
        call freexi(ksav)
        call freexi(npmsav)
        call freexd(buf)

        return
        end
