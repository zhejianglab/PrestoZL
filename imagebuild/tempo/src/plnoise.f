c      $Id$

c Functions for computing covariance for power-law red noise
c spectra.

	real*8 function plnoise(dt,alpha,f0,amp,init)

c Compute cov for lag dt (in years).  Basically follows
c formula from van Haasteren et al (2011).  Set init to true
c the first time to set up coeffs, etc.  Input lag and freq in 
c years, output in in us^2.

        real*8 dt,alpha,f0,amp
        logical init

        real*8 fl
        data fl/0.01/

        real*8 pi
        data pi/3.141592653589793d0/

        parameter(NC=32)
        real*8 scale0, scale1
        real*8 log_coeffs(NC), signs(NC)
        real*8 ctmp,xtmp,lxtmp,ssum
        real*8 exval

        save pi, fl, scale0, scale1, signs, log_coeffs

c exval is the value that appears in many exponents in the
c formula.  If input alpha is the strain spectral index, then:
        !exval = 2.0-2.0*alpha
c OR, if alpha is instead the timing power index, then:
        exval = -1.0-alpha

c This output scaling assumes input amp is in strain units:
c          scale0 = f0**(exval-2.0)/(12.0*pi**2)/(fl**(exval))
c          scale0 = scale0 * (365.24*86400.0)**2 * 1d12 ! us^2 output
c If instead amp is in "observer units" of us / sqrt(yr^-1) then:
          scale0 = f0 * (f0/fl)**exval

c First-time initialization, only depends on spec idx
        if (init) then
          scale1 = gamma(-exval) * cos(pi*(1.0-exval/2.0))
          do i=0,NC-1
            if (mod(i,2).eq.1) then
              signs(i+1) = -1.0
            else
              signs(i+1) = 1.0
            endif
            ctmp = 2.0d0*i - exval
            if (ctmp.lt.0.0) then
              signs(i+1) = -1.0 * signs(i+1)
              ctmp = abs(ctmp)
            endif
            log_coeffs(i+1) = -log_gamma(2.0d0*i+1.0) - log(ctmp)
            !print *,"LC",i,log_coeffs(i+1)
          enddo
        endif

        if (dt.eq.0.0) then
          plnoise = (amp**2)*scale0*(-signs(1)*exp(log_coeffs(1)))
          return
        endif

        xtmp = abs(2.0d0*pi*fl*dt)
        lxtmp = log(xtmp)
        ssum = 0.0
        do i=0,NC-1
          ssum=ssum+signs(i+1)*exp(2.0d0*i*lxtmp + log_coeffs(i+1))
        enddo
        plnoise = (amp**2)*scale0*(-scale1*(xtmp**(exval))-ssum)

        return
        end

c ===================================================================

	real*8 function plnoise_interp(dt,alpha,f0,amp,init)

        real*8 plnoise

        real*8 dt,alpha,f0,amp
        logical init

        integer idx0
        real*8 frac

        parameter(NPTS=10000) ! max lag in days
        real*8 lu(NPTS)

        save lu

        if (init) then
          lu(1) = plnoise(0d0,alpha,f0,amp,.true.)
c          print *,"PL0 =", lu(1)
          do i=1,NPTS-1
            lu(i+1) = plnoise(i/365.24d0,alpha,f0,amp,.false.)
          enddo
        endif

c Assume dt is in days.. should really make these all consistent...
        idx0 = floor(abs(dt))
        frac = abs(dt) - idx0
        plnoise_interp = (1d0 - frac)*lu(idx0+1) + frac*lu(idx0+2)

        return
        end


