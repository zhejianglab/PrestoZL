c      $Id$
      subroutine tzfit(ipsr,phase,f0,dmpsr,ct2,t0,pb)

      implicit real*8 (a-h,o-z)
      parameter (nn=31,pi=3.141592653589793d0)
      character date*9,damoyr*9,binpha*16

      include 'dim.h'
      include 'acom.h'
      include 'tz.h'
      include 'trnsfr.h'

      real*8 t(nn),ph(nn),c(nn),d(nn),phase(NTZARR)
      common /obsp/ site(3), pos(3), freqhz, bval, sitvel(3)

      nspan=nsp(ipsr)
      ncoeff=nco(ipsr)
      tzfreq=tzof(ipsr)

      nh=nn/2
      do 100 nb=1,nsets
        iref=nb*nn-nh
        nrmjd=ntmjd(iref)
        frmjd=ftmjd(iref)				
        rtime=tmin(iref)			
        rphase=phase(iref)	
        
        i=iref-nh-1
        do 10 j=1,nn	
          i=i+1
          t(j)=tmin(i)-rtime
          ph(j)=phase(i)-rphase-t(j)*f0*60.d0	
 10     continue
        fac=2.d0/nn
        do 13 j=1,nn
          sum=0.d0
          do 12 k=1,nn
            sum=sum+ph(k)*cos((pi*(j-1))*((k-0.5d0)/nn))
 12       continue
          c(j)=fac*sum
 13     continue
        
        b=nspan/2 + 5
        a=-b
        call chebpc(c,d,ncoeff)
        call pcshft(a,b,d,ncoeff)
        sq=0.
        do 15 j=1,nn
          phtst=d(1)
          do 14 k=2,ncoeff
            phtst=phtst+d(k)*t(j)**(k-1)
 14       continue
          sq=sq+(phtst-ph(j))**2
 15     continue
        rms=1.d6*sqrt(sq/nn)/f0
        
        nutsec=mod(86400.d0*(nrmjd+frmjd)+0.005d0,86400.d0)
        nuthrs=nutsec/3600
        nutmin=(nutsec-3600*nuthrs)/60
        uts=nutsec-3600*nuthrs-60*nutmin
        utprint=10000.d0*nuthrs + 100.d0*nutmin + uts
        date=damoyr(int(nrmjd))
        if (.not.quiet) 
     +       write(*,1051) pname,date,nuthrs,nutmin,rms
 1051   format(' PSR ',a12,1x,a9,'  UTC: ',i4.2,i2.2,'   rms:',f9.3,
     +       ' us.')
        
        binpha=' '
        if(pb.gt.0.d0) then
          ct=ct2+(nrmjd-ntmjd(1))+(frmjd-ftmjd(1))
          phifac=8.64d4/pb
          phi0=dmod((ct-t0)*phifac+90000.d0,1.d0)
          write(binpha,1055) phi0,phifac
 1055     format(f7.4,f9.4)
        endif

c	get doppler shift
c       ztim will shift freqhz (in Hz) to frqf (in MHz);
c       these are passed in common blocks defined above
        if (nsite.ge.0) then
          freqhz = 1.d9
          call ztim(nrmjd,frmjd,nct,fct)
          freqhz = freqhz * 1.d-6
          z = (frq-freqhz)/freqhz
          z = z * 1.e4
        else
          z = 0.d0
        endif

c  log_10 of rms fit error in periods
        rms = dlog10(1.d-6 * rms * f0)
        
        nx = nrmjd/10
        fx = (nrmjd-10*nx) + frmjd
        write(lupolyco,1070) pname(1:10),date,utprint,nx,fx,dmpsr,z,rms,
     +       rphase,f0,nsite,nspan,ncoeff,tzfreq,binpha	
 1070   format(a10,1x,a9,f11.2,i7,f13.11,f21.6,1x,f6.3,f7.3/
     +       f20.6,f18.12,3i5,f10.3,a16)
        write(lupolyco,1080) (d(i),i=1,ncoeff)		
 1080   format(3d25.17)
        
 100  continue
      return
      end



