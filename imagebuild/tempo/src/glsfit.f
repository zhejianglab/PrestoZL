c      $Id$
      subroutine glsfit(npts,mode,chisqr,varfit,xmean,ymean,sum,nz,
     +     wmax,lw,ddmch,
     +     buf,npmsav,ksav,
     +     resfile2)

	implicit real*8 (a-h,o-z)
	implicit integer*8 (i-k)
	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'vcom.h'
	include 'orbit.h'
	include 'tz.h'
	include 'dp.h'
        include 'toa.h'

c This is a version of the original tempo fit() function, rewritten
c to use a generalized least squares (GLS) approach, which can use
c a non-diagonal 'weighting' matrx.  This uses the LAPACK routines
c for various things.  As much as possible this is consistent with
c the original fit() in terms of input/outputs.
c PBD 2014/04

c       moved declaration of real*8 array(NPA,NPA) to acom.h, djn, 8 Sep 98
	real*8 xmean(NPA),fctn(NPAP1)
	real*8 a(NPAP1),atmp(NPAP1),sigmaa(NPAP1),gcor(NPA)
        real*8 mcov(NPAP1), mscal(NPAP1)
	logical lw
	real*8 ddmch(*)
        real*8 buf(*)
        integer npmsav(*), ksav(*)
	character mark*1,adn*14,date*9,damoyr*9
        character*80 resfile2

        integer isetflag

        integer resn(20)  ! small buffer to hold one toa's worth of information
	real*8 resr(9)
        equivalence (resn(2),resr(1)) ! real*8 will not be on 8-byte boundaries

c Just define some huge matrices here for now.  Once it's working, try
c replacing using the "malloc" routines?
c        real*8 Adm(2*NPTSDEF,NPAP1) ! weighted fit design matrix
c                                    ! factor of 2 for TOA+DM
        real*8 VTsvd(NPAP1,NPAP1) ! SVD "V-transpose" result
        real*8 sv(NPAP1) ! The singular values
        real*8 r(2*NPTSDEF) ! The weighted prefit TOA+DM resids
        real*8 rr(NPTSDEF)  ! red residuals
        real*8 cts(NPTSDEF) ! copy of the times
        real*8 ecorr(NPTSDEF) ! ECORR values per-TOA
        real*8 work(10*NPAP1*NPAP1)
        integer lwork
        integer iwork(10*NPAP1)
        integer inforv

c Design matrix to be allocated. Sizes:
c   Ntoa * Nparam normally
c   2 * Ntoa * Nparam if DMDATA
        real*8,allocatable :: Adm(:)
        integer Admrows
C        integer*8 Admoff

c DM-data-related stuff
c Currently requires:
c   1. DM cov matrix is diagonal.
c   2. DM and TOAs are not covariant.
c   3. Only DMX (no DMX1 or DMnn) is fit.
        real*8 dmwt
        integer ndmparam       ! number of DM fit params

c packed cov matrix, to be malloced
        logical havecov
        data havecov/.false./
        logical diagcov
        data diagcov/.true./
        integer*8 ncovpts  
        integer*8 idx
        real*8,allocatable :: dcov(:) ! data cov matrix
        real*8,allocatable :: rcov(:) ! red noise portion
        real*8 detcov,cmax,cmin

        real*8 rmean, r2mean

        character*80 getvalue ! need to declare ext function..
        character*80 flagtmp

	integer fd
	integer nwrt
	integer flags, filemode  
 	integer open, close, write

c save the cov matrix stuff so we can iterate faster
        save havecov, ncovpts, dcov, detcov, Adm
	save rcov

        lwork = 10*NPAP1*NPAP1
	mprt=10**nprnt
	sigma=0.
	chisq=0.
        dm_chisq=0.
        toa_chisq=0.
        rmean=0.
        r2mean=0.
        ndmparam=0

c GLS mode does not work if MODE 0 is specified
        if (mode.eq.0) stop "glsfit: MODE 0 is not allowed"

c nparam is total number of fit params (including mean)
c Zero out various matrices
	nterms=nparam-1
	do 29 j=1,nterms
          xmean(j) = xmean(j)/sum
          do 28 k=1,nterms
            array(j,k)=0.
 28       continue
 29     continue

c Allocate design matrix
c TODO what about iterations...
c index using Adm(i,j) -> Adm(Admoff + i + (j-1)*Admrows)
        if (usedmdata) then
          Admrows = 2*npts
        else 
          Admrows = npts
        endif
        if (.not. havecov) then 
          allocate(Adm(Admrows*nparam),stat=istat)
          if (istat.ne.0) stop "glsfit: can't allocate Adm"
        endif
        do i=1,Admrows*nparam
          Adm(i) = 0.
        enddo

        do 32 j=1,NPAP1
          a(j)=0.
          atmp(j)=0.
          sigmaa(j)=0.
          sv(j)=0.
          mcov(j)=0.
          mscal(j)=0.
          do 30 k=1,NPAP1
            VTsvd(k,j)=0.
 30       continue
 32     continue 

        write (*,'(''  glsfit Ntoa='',i6)') npts
        if (nz.gt.0) print *,"Warning: nz>0 in glsfit"

c Alloc space for COV matrix.  
c Packed upper triangular storage means element (i,j) is accessed
c using index i+j*(j-1)/2, with i<=j.
c If red noise is enabled we save a separate copy of the red noise
c portion of the cov matrix so that 'white' residuals can be calculated
c at the end.
        if (.not. havecov) then
          print *,'  ... allocate matrices'
          ncovpts = int(npts,8)*(npts-1)/2 + npts
          allocate(dcov(ncovpts),stat=istat)
          if (istat.ne.0)
     +      stop "glsfit: can't allocate memory for dcov"
          if (rnamp.gt.0) then
            allocate(rcov(ncovpts),stat=istat)
            if (istat.ne.0)
     +        stop "glsfit: can't allocate memory for rcov"
          endif
          do i=1,ncovpts 
            dcov(i)=0.0
            if (rnamp.gt.0) rcov(i)=0.0
          enddo
        endif

c This is the part where we read in the data
c vmemr call reads info for TOA number i from memory
c Inputs:
c   y is the (pre-fit) residual, in pulse phase (turns)
c   fctn are the fit basis funcs (param derivs) evaluated for this TOA
c   weight is the TOA's weight (units phase^-2)
c   terr is the TOA's uncertainty (us, including efac/equad)
c   ct is the infinite freq barycentric time (MJD, double)
c   fmjd is the site time (MJD, double)
c Computed here:
c   Adm, design matrix
c   r, pre-fit resids
c   dcov, cov matrix (diagonal part only)
        print *,'  ... fill design matrix'
	if(ldesign) rewind(37)
	if(ldesign) write(37) npts, nterms
	do 67 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          if(ldesign) write(37) ct, weight, (fctn(j)-xmean(j), j=1,nterms)
          !cts(i) = ct
          cts(i) = fmjd
          r(i) = y
          rr(i) = 0.0
          Adm(i) = 1d0 ! Constant phase term
          do 66 j=2,nparam
            Adm(i+(j-1)*Admrows) = (fctn(j-1)-xmean(j-1))
            if (usedmdata .and. (mfit(j).eq.16 .or. 
     +        (mfit(j).gt.NPAR6 .and. mfit(j).le.NPAR7))) then
              ! This is a DM/DMX param, fill in the "extra" DM part of
              ! Adm.  
              !  TODO:
              !    - This assumes(!) only DMXs (no DMX1s) are set.
              !    - Extend to deal with DM polynomials.
              ! The test for fctn!=0 makes sure this TOA is affected by
              ! this param (ie is in the relevant DMX bin).
              if (fctn(j-1).ne.0d0) Adm(i+npts+(j-1)*Admrows) = 1.0
              ndmparam = ndmparam + 1
            endif
 66       continue
          if (.not.havecov) dcov(i+i*(i-1)/2) = 1d0/weight
 67     continue

        if (.not.havecov) then

C compute red noise here so that we can save it separately
c Power-law red noise:
c in "timing power" units (us/sqrt(yr^-1))
          !rnamp = 0.028997
          !rnidx = -13d0/3d0
          if (rnamp.gt.0) then
C call once to initialize plnoise values:
            print *,'  ... compute covariance (rn)'
            z = plnoise_interp(0d0,rnidx,1d0,rnamp,.true.)
            do i=1,npts
              cp = p0+p1*(cts(i)-peopch)*86400.d0
              do j=1,i
                idx=j+i*(i-1)/2
                rcov(idx) =
     +            plnoise_interp(abs(cts(i)-cts(j)),rnidx,1d0,rnamp,.false.)
     +            / cp**2 / 1d12
              enddo
            enddo
          endif

c read inverted cov matrix from disk here if option selected
c TODO some kind of check that the file works would be nice... 
          if (dcovfile.ne."") then

            print *,'  ... read cov matrix from disk'
	    open(39,file=dcovfile,form='unformatted',
     +            status='unknown')
	    read(39) (dcov(i), i=1,ncovpts)
	    close(39)

          else

c Add in extra cov matrix terms
c Use flag-based ECORR ("jitter") values
C This should probably be a subroutine?
            if (nflagecorr.gt.0) then
              if (tcorr.eq.0.0) tcorr = 2d0*p0/86400d0
              print *,'  ... compute covariance (ecorr)'
              diagcov = .false.
              do i=1,npts
                cp = p0+p1*(cts(i)-peopch)*86400.d0
                ecorr(i) = getecorr(stflags(i)) * 1d-6 / cp
                do j=1,i
                  idx=j+i*(i-1)/2
                  if (abs(cts(i)-cts(j)).lt.tcorr) then
                    dcov(idx) = dcov(idx) + ecorr(i)*ecorr(j)
                  endif
                enddo
              enddo
            endif

C Add in red noise component if it exists
            if (rnamp.gt.0) then
              cmax = dcov(1)
              cmin = dcov(1)
              diagcov = .false.
              do i=1,npts
                do j=1,i
                  idx=j+i*(i-1)/2
                  dcov(idx) = dcov(idx)
     +              + rcov(idx)
                  if (dcov(idx).gt.cmax) cmax = dcov(idx)
                  if (dcov(idx).lt.cmin) cmin = dcov(idx)
                enddo
              enddo
c            print *,"dcov(max)=",cmax
c            print *,"dcov(min)=",cmin
c            do i=1,10
c              print *,i,dcov(i)
c            enddo
            endif

            if (diagcov) then
C if cov matrix is diagonal just do it the easy way 

              print *,'  ... invert diagonal matrix'
              do i=1,npts
                idx=i+i*(i-1)/2
                dcov(idx) = 1.0/sqrt(dcov(idx))
              enddo

            else
C Non-diagonal cov, do the full cholesky thing

c Notes for GLS:
c Cholesky factorize: dpotrf() or dpptrf() for packed
              print *,'  ... Cholesky decomp'
              call dpptrf('U',npts,dcov,inforv)
              !print *,"dpptrf info=",inforv
              if (inforv.ne.0) stop "glsfit: Cholesky decomp failed"
              detcov=0d0
              cmax = dcov(1)
              cmin = dcov(1)
              do i=1,npts
                idx = i+i*(i-1)/2
                detcov = detcov + log(dcov(idx))
                if (dcov(idx).gt.cmax) cmax = dcov(idx)
                if (dcov(idx).lt.cmin) cmin = dcov(idx)
              enddo
              print *,"      log(det(cov))=",detcov
              !print *,"cmax=",cmax
              !print *,"cmin=",cmin

c invert triangular matrix: dtrtri() or dtptri() for packed
              print *,'  ... matrix inverse'
              call dtptri('U','N',npts,dcov,inforv)
              if (inforv.ne.0) stop "glsfit: invert cov matrix failed"

            endif

c write inverted/cholesky'd matrix to disk here 
c TODO make this optional
            if (ldcov) then
              print *,'  ... write matrix to disk'
              open(39,file='datacov.tmp',form='unformatted',
     +              status='unknown')
              rewind(39)
              write(39) (dcov(i), i=1,ncovpts)
              close(39)
            endif

          endif ! readcov

          havecov = .true.
        else
          print *,'  ... using cached cov values'
        endif ! not havecov

c        do i=1,10
c          print *,i,dcov(i)
c        enddo

c scale columns of design matrix
        do j=1,nparam
          mscal(j) = dnrm2(Admrows,Adm(1+(j-1)*Admrows),1)
        enddo
	do i=1,Admrows
	  do j=1,nparam
            Adm(i+(j-1)*Admrows)=Adm(i+(j-1)*Admrows)/mscal(j)
	  enddo
	enddo

c Multiply by inv cov matrix
        print *,'  ... matrix multiply'
        call dtpmv('U','T','N',npts,dcov,r,1)
        do 68 i=1,nparam
          call dtpmv('U','T','N',npts,dcov,
     +      Adm(1+(i-1)*Admrows),1)
 68     continue

c Scale DM and DM part of Adm by inv DM uncertainties (only
c allow diagonal DM cov matrix for now)
        if (usedmdata) then
          print *,'  ... using DM data points'
          do i=1,npts
            dmwt = 0.0
            if (dmerr(i).gt.0.0) dmwt=1.0/dmerr(i)
            r(i+npts) = dmres(i)*dmwt
            do j=1,nparam
              Adm(i+npts+(j-1)*Admrows) = 
     +          Adm(i+npts+(j-1)*Admrows)*dmwt
            enddo
          enddo
        endif

c remove cov with mean from other params
c Not sure how necessary this is since functions were already
c mean-subtracted (unweighted)..
        do j=1,nparam
          mcov(j) = ddot(npts,Adm(1),1,
     +                     Adm(1+(j-1)*Admrows),1)
          if (j.gt.1) then
            mcov(j) = mcov(j)/mcov(1)
            call daxpy(npts,-mcov(j),Adm(1),1,
     +                    Adm(1+(j-1)*Admrows),1)
          endif
        enddo

c Call SVD routine.  On output, Adm will be replaced by "U".
        print *,'  ... SVD'
        call dgesdd('O',Admrows,nparam,Adm(1),Admrows,sv,
     +    Adm(1),Admrows,VTsvd,NPAP1,work,lwork,iwork,inforv)

c TODO could test for low SV's here
        write(*,'(''      log10(cond number) = '',f5.2)')
     +    log10(sv(1)) - log10(sv(nparam))

c Fill in "array" param cov matrix values
        do j=1,nterms
          do k=1,nterms
            array(j,k) = 0.0
            do i=1,nparam
              array(j,k)=array(j,k)+VTsvd(i,j+1)*VTsvd(i,k+1)/sv(i)/sv(i)
          enddo
        enddo
       enddo
       do j=1,nterms
         do k=1,nterms
           array(j,k) = array(j,k)/mscal(j+1)/mscal(k+1)
         enddo
       enddo

c uncertainty on mean
       mcov(1) = 0.0
       do i=1,nparam
         mcov(1) = mcov(1) + VTsvd(i,1)*VTsvd(i,1)/sv(i)/sv(i)
       enddo

c Need to keep "gcor" results for output?
c what is the reference for this formula??
	do 80 i=1,nterms
c          gcor(i)=sqrt(abs(1.d0 - zzz/array(i,i)))
          gcor(i)=0.0
 80     continue

	if(lw) then
	  write(31,1170)
1170	  format(//'Normalized covariance matrix in "number of 9''s"',
     +      ' format:'/)
	  call mxprt(array,gcor,nterms,mfit,nbin,eclcoord)
	  unit=1.d-6
	  if(p0firs.gt.0.1) unit=1.d-3
	  nunit=dlog10(unit)-0.0001
	  write(31,1050) nunit
1050	  format(//'   N     TOBS         Date',
     +     '      FREQ    WGT      NPULSE     R1(P)  R2(',i2,') PH/DT'/)
	endif

c Compute post-fit resids in "whitened" basis
c r_post = r_pre - U U^t r_pre
c TODO make tempo report this as "the" chi2?
        call dgemv('T',Admrows,nparam,1d0,Adm(1),
     +              Admrows,r,1,0d0,atmp,1)
        call dgemv('N',Admrows,nparam,-1d0,Adm(1),
     +              Admrows,atmp,1,1d0,r,1)
        if (usedmdata) then
          toa_chisq = ddot(npts,r,1,r,1)
          dm_chisq = ddot(npts,r(npts+1),1,r(npts+1),1)
          print *,"TOA chisq=", toa_chisq
          print *,"DM chisq=", dm_chisq
          print *,"Ndof=",2*npts-nparam
          chisq = ddot(2*npts,r,1,r,1)
        else
          chisq = ddot(npts,r,1,r,1)
        endif
        print *,"GLS chisq=", chisq

c Here is where we compute red residuals
        if (rnamp.gt.0) then
          print *,'  ... computing red residuals'
          call dtpmv('U','N','N',npts,dcov,r,1)
          call dspmv('U',npts,1d0,rcov,r,1,0d0,rr,1)
        else
          do i=1,npts
            rr(i) = 0d0
          enddo
        endif

c Computes the post-fit param values in a
	do i=1,nparam
          atmp(i) = atmp(i)/sv(i)
        enddo
        call dgemv('T',nparam,nparam,1d0,VTsvd,NPAP1,atmp,1,0d0,a,1)
        do i=1,nparam
          a(i) = a(i)/mscal(i)
        enddo

c Note a(1) is the const phase term, aa0 is var name from orig fit.f
        aa0 = a(1)

        flags = isetflag()
	filemode  = 6*64 + 6*8 + 2  ! octal 662 = rw-rw-r--
	if (lw) fd = open(resfile2,flags,filemode)
	resn(1) = 72
	resn(20) = 72

c re-read data, compute post-fit residuals in dt2, output 
c them to resid2.tmp and tempo.lis
        fnpts = npts-nz
        wmean = sum/fnpts
	do 108 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          dt2=y-aa0
          weight=weight/wmean
          do 107 j=1,nterms
c j+1 here due to a(1)==aa0
            dt2=dt2-a(j+1)*(fctn(j)-xmean(j)-mcov(j+1)*mscal(j+1)/mscal(1))
 107      continue
          nct=ct
          mark=char(mod(nct,26)+65)
          date=damoyr(int(fmjd))
          dt2sec=dt2*p0firs
          phase=0.
          if(a1(1).ne.0.) then
            if (nbin.ne.9) then
               phase=dmod((ct-t0(1))*86400.d0/pb(1)+1000000.d0,1.d0)
            else
               phase=dmod((ct-t0asc)*86400.d0/pb(1)+1000000.d0,1.d0)
            endif
          endif

          if(nprnt.eq.0.or.mod(i,mprt).eq.1) then
            write(adn,1106) dn
 1106       format(f14.0)

            if(lw) write(31,1107) i,fmjd,mark,date,
     +           rfrq,weight,adn(1:13),y,dt2sec/unit,phase
 1107       format(i4,f12.5,1x,a1,1x,a9,f9.3,f6.3,a13,
     +           f8.3,f8.2,f7.3)
          endif

C Correct tz ref TOA
          if(i.eq.ntzref) then
            ftzrmjd=ftzrmjd-dt2sec/8.64d4
            if (ftzrmjd.lt.0.) then
              ftzrmjd = ftzrmjd + 1
              ntzrmjd = ntzrmjd - 1
            elseif (ftzrmjd.ge.1.) then
              ftzrmjd = ftzrmjd - 1
              ntzrmjd = ntzrmjd + 1
            endif
          endif

	  resr(1) = ct
	  resr(2) = dt2
	  resr(3) = dt2sec
	  resr(4) = phase
	  resr(5) = frq
	  resr(6) = weight
	  resr(7) = terr
	  resr(8) = y
	  resr(9) = rr(i)*p0firs
          if (lw) nwrt = write(fd,resn,80)
          wmax=max(wmax,weight)
          r2mean=r2mean+weight*dt2**2
          rmean=rmean+weight*dt2
 108    continue
	if (lw) nwrt = close(fd)

	freen=fnpts-nterms-1 ! TODO fix for DMs
	!chisqr=r2mean*wmean/freen ! "old" chi2
        chisqr=chisq/freen ! GLS chi2
        rmean=rmean/fnpts
        r2mean=r2mean/fnpts
	varfit=r2mean - rmean*rmean
        !if(mode.eq.0) chisqr=0. ! don't think mode=0 makes sense here

	do 133 j=1,nterms
          sigmaa(j+1)=dsqrt(array(j,j))
 133    continue

	sigmaa(1)=dsqrt(mcov(1))

	do 146 j=1,NPAP1
          freq(j)=0.
          ferr(j)=0.
 146      continue

	do 150 j=1,nparam
          freq(mfit(j))=a(j)
          ferr(mfit(j))=sigmaa(j)
 150    continue

	return
	end

c ===========================================================

        real*8 function getecorr(rawflags)

	implicit real*8 (a-h,o-z)
        include 'dim.h'
        include 'acom.h'

        character*80 rawflags
        character*80 tmp
        character*80 getvalue

        getecorr = 0d0

        j1 = 1
        call getflags(rawflags,640,j1)
        do i=1,nflagecorr
          tmp = getvalue(ecorrflag(i)(2:32))
          if (tmp.eq.ecorrflagval(i)) then
            getecorr = flagecorr(i)
          endif
        enddo

        return
        end

