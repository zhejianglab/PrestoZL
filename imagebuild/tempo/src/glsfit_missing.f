c      $Id$
      subroutine glsfit(npts,mode,chisqr,varfit,xmean,ymean,sum,nz,
     +     wmax,lw,ddmch,
     +     buf,npmsav,ksav,
     +     resfile2)

c This is a dummy version of glsfit() to be linked in when lapack
c is not available.  It will print a error and exit
c PBD 2014/04

        stop "Tempo was not compiled with GLS fit support " 
     +    // "(requires LAPACK)."

	return
	end
