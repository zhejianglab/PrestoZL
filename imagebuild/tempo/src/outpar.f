c      $Id$
      subroutine outpar(nits,irh,irm,rsec,ers,decsgn,idd,idm,dsec,eds)

      implicit real*8 (A-H,O-Z)
      character decsgn*1, fit1*3

      parameter (TWOPI=6.28318530717958648d0)

      include 'dim.h'     
      include 'acom.h'
      include 'bcom.h'
      include 'trnsfr.h'
      include 'clocks.h'
      include 'dp.h'
      include 'glitch.h'
      include 'eph.h'
      include 'tz.h'

      character*1 siten2a  ! external function
      character*2 siten2b  ! external function

      c1=360.d0/TWOPI
      c=360.*3600./TWOPI
      cc=1.d-9*c*365.25*8.64d7
      fit1='  1'

      write(71,'(a,''              '',a)')
     +      psrkey(1:psrkeyl),psrname(1:index(psrname," ")-1)

      if(eclcoord)then
         if(nfit(6).gt.0)then
            write(71,1011)c1*pra,fit1,ferr(6)
         else
            write(71,1011)c1*pra
         endif
 1011    format('LAMBDA',f20.13,a,f20.13)
         
         if(nfit(5).gt.0)then
            write(71,1012)c1*pdec,fit1,ferr(5)
         else
            write(71,1012)c1*pdec
         endif
 1012    format('BETA',f22.13,a,f20.13)
         
         if(pmra.ne.0.)then
            if(nfit(8).gt.0)then
               write(71,1013)pmra,fit1,ferr(8)*cc
            else
               write(71,1013)pmra
            endif
         endif
 1013    format('PMLAMBDA',f18.4,a,f20.4)
         
        if(pmdec.ne.0.)then
           if(nfit(7).gt.0)then
              write(71,1014)pmdec,fit1,ferr(7)*cc
           else
              write(71,1014)pmdec
           endif
        endif
 1014   format('PMBETA',f20.4,a,f20.4)
        
      else
        irs=rsec
        rs=rsec-irs
        ids=dsec
        ds=dsec-ids
        if(nfit(6).gt.0)then
           write(71,1021)irh,irm,irs,rs,fit1,ers
        else
           write(71,1021)irh,irm,irs,rs
        endif
 1021   format('RAJ',i8.2,':',i2.2,':',i2.2,f9.8,a,f20.8)
        
        if(nfit(5).gt.0)then
           write(71,1022)decsgn,idd,idm,ids,ds,fit1,eds
        else
           write(71,1022)decsgn,idd,idm,ids,ds
        endif
 1022   format('DECJ',5x,a,i2.2,':',i2.2,':',i2.2,f8.7,a,f20.7)
        
        if(pmra.ne.0.)then
           if(nfit(8).gt.0)then
              write(71,1023)pmra,fit1,ferr(8)*cc
           else
              write(71,1023)pmra
           endif
        endif
 1023   format('PMRA',f22.4,a,f20.4)
        
        if(pmdec.ne.0.)then
           if(nfit(7).gt.0)then
              write(71,1024)pmdec,fit1,ferr(7)*cc
           else
              write(71,1024)pmdec
           endif
        endif
 1024   format('PMDEC',f21.4,a,f20.4)
        
      endif

      if(pmrv.ne.0.)then
         if(nfit(36).gt.0)then
            write(71,1025)pmrv,fit1,10.*c*ferr(36)
         else
            write(71,1025)pmrv
         endif
      endif
 1025 format('PMRV',f22.4,a,f20.4)

      if(px.ne.0.)then
         if(nfit(17).gt.0)then
            write(71,1026)px,fit1,ferr(17)
         else
            write(71,1026)px
         endif
      endif
 1026 format('PX',f24.4,a,f20.4)

      if(posepoch.gt.0.)write(71,1027)posepoch
 1027 format('POSEPOCH',f18.4)

        
      if(nfit(2).gt.0)then
         write(71,1031)f0,fit1,ferr(2)*1.d-9
      else
         write(71,1031)f0
      endif
 1031 format('F0',f24.16,a,f20.16)

      if(nfit(3).gt.0)then
         write(71,1032)f1,fit1,ferr(3)*1.d-18
      else
         write(71,1032)f1
      endif
 1032 format('F1',1p,d24.12,a,d20.12)

      if(nfit(4).gt.0)then
         write(71,1033)f2,fit1,ferr(4)*1.d-27
      else if (f2.ne.0.) then
         write(71,1033)f2
      endif
 1033 format('F2',1p,d24.12,a,d20.12)

      if(nfit(NPAR11+1).gt.0)then
         write(71,1034)f3,fit1,ferr(NPAR11+1)*1.d-36
      else if (f3.ne.0)then
         write(71,1034)f3
      endif
 1034 format('F3',1p,d24.12,a,d20.12)

      do i = 1, nfcalc-3
         if(nfit(NPAR11+1+i).gt.0)then   ! if fit, print 
	    if (i+3.lt.10) then
            write(71,1035)i+3,f4(i),fit1,ferr(NPAR11+1+i)*(1.d-9)**(i+4)
            else
            write(71,1036)i+3,f4(i),fit1,ferr(NPAR11+1+i)*(1.d-9)**(i+4)
            endif	
         else 
           if (i+3.lt.10) then
             write(71,1035)i+3,f4(i)
           else
             write(71,1036)i+3,f4(i)
           endif
         endif
      enddo
 1035 format('F',i1,1p,d24.12,a,d20.12)
 1036 format('F',i2,1p,d23.12,a,d20.12)

      write(71,'(''PEPOCH'',f20.6)')pepoch
      if (usestart) then
        write(71,'(''START'',f21.3,a)')start,fit1
      else
        write(71,'(''START'',f21.3)')start
      endif
      if (usefinish) then
        write(71,'(''FINISH'',f20.3,a)')finish,fit1
      else
        write(71,'(''FINISH'',f20.3)')finish
      endif

      if(nfit(16).gt.0)then
         write(71,1050)dm,fit1,ferr(16)
      else
         write(71,1050)dm
      endif
 1050 format('DM',f24.6,a,f20.6)

      do i=1,ndmcalc-1
         if(dmcof(i).ne.0)then
            if(nfit(NPAR9+i).gt.0)then
               write(71,1051)i,dmcof(i),fit1,ferr(NPAR9+i)
            else
               write(71,1051)i,dmcof(i)
            endif
         endif
      enddo
 1051 format('DM',i3.3,1p,d21.12,a,d20.12)

      if(dmepoch.gt.0.)write(71,1052)dmepoch
 1052 format('DMEPOCH',f18.4)

      if(ngl.gt.0)then
         do i=1,ngl
            ii=NPAR1+(i-1)*NGLP
            write(71,'(''GLEP_'',i1,f18.6)')i,glepoch(i)
            if(nfit(ii+1).gt.0)then
               write(71,1061)i,glph(i),fit1,ferr(ii+1)
            else
               write(71,1061)i,glph(i)
            endif
 1061       format('GLPH_',i1,f20.6,a,f20.6)
            if(nfit(ii+2).gt.0)then
               write(71,1062)i,glf0(i),fit1,ferr(ii+2)*1.d-9
            else
               write(71,1062)i,glf0(i)
            endif
 1062       format('GLF0_',i1,1p,d20.8,a,d20.8)
            if(nfit(ii+3).gt.0)then
               write(71,1063)i,glf1(i),fit1,ferr(ii+3)*1.d-18
            else
               write(71,1063)i,glf1(i)
            endif
 1063       format('GLF1_',i1,1p,d20.8,a,d20.8)
            if(nfit(ii+4).gt.0)then
               write(71,1064)i,glf0d(i),fit1,ferr(ii+4)*1.d-9
            else
               write(71,1064)i,glf0d(i)
            endif
 1064       format('GLF0D_',i1,1p,d19.8,a,d20.8)
            if(nfit(ii+5).gt.0)then
               write(71,1065)i,gltd(i),fit1,ferr(ii+5)/86400.d0
            else
               write(71,1065)i,gltd(i)
            endif
 1065       format('GLTD_',i1,f20.4,a,f20.4)
         enddo
      endif

      if (usedmx) then
        write (71,1079),dmxt
 1079   format ('DMX',f23.6)
        do i = 1, ndmx
          write (71,1080) i,dmx(i),nfit(NPAR6+2*i-1),ferr(NPAR6+2*i-1)
 1080     format('DMX_',i4.4,1p,d18.8,i3,d20.8)
          if (dmx1(i).ne.0 .or. nfit(NPAR6+2*i).ne.0) 
     +      write (71,1081) i,dmx1(i),nfit(NPAR6+2*i),ferr(NPAR6+2*i)
 1081     format('DMX1_',i4.4,1p,d17.8,i3,d20.8)
          write (71,1082) i,dmxep(i)
 1082     format('DMXEP_',i4.4,f16.5)
c The following lines round DMXR1 down and DMXR2 up at the 
c output precision (5 frac digits) so that the ranges will
c work correctly in tempo2:
          write (71,1083) i,dint(dmxr1(i)*1.d5)/1.d5
 1083     format('DMXR1_',i4.4,f16.5)
          write (71,1084) i,dint(dmxr2(i)*1.d5+1.0)/1.d5
 1084     format('DMXR2_',i4.4,f16.5)
          write (71,1085) i,dmxf1(i)
 1085     format('DMXF1_',i4.4,f16.3)
          write (71,1086) i,dmxf2(i)
 1086     format('DMXF2_',i4.4,f16.3)
        enddo
      endif
      if (dmvar1.ne.0.) write (71,1095) 1,dmvar1
      if (dmvar2.ne.999999.) write (71,1095) 2,dmvar2
 1095 format ('DMVAR',i1,f25.8)
      if (.not.firstdmx) write (71,1096)
 1096 format ('DMXFIX 1')

      do i=1,NFDMAX
        if (fdcof(i).ne.0) then
          if(nfit(NPAR10+i).gt.0) then
            write(71,1100)i,fdcof(i),fit1,ferr(NPAR10+i)
          else
            write(71,1100)i,fdcof(i)
          endif
 1100     format('FD',i1,1p,d16.8,a,d20.8)
	endif
      enddo

      if (usexmxfrq0) then
         write (71, 1201) xmxfrq0
 1201    format('XMXFRQ0',12x,f16.5)
      endif
      do i= 1, nxmx
        if (xmxuse(i)) then
          write (71,1202) i,xmx(i),nfit(NPAR12+2*i-1),ferr(NPAR12+2*i-1)
 1202     format('XMX_',i4.4,1p,d18.8,i3,d20.8)
          write (71,1203) i,xmxexp(i),nfit(NPAR12+2*i),ferr(NPAR12+2*i)
 1203     format('XMXEXP_',i4.4,1p,d18.8,i3,d20.8)
          if (xmxr1(i).gt.0.) write(71,1204) "R1",i,xmxr1(i)
          if (xmxr2(i).gt.0.) write(71,1204) "R2",i,xmxr2(i)
          if (xmxf1(i).gt.0.) write(71,1204) "F1",i,xmxf1(i)
          if (xmxf2(i).gt.0.) write(71,1204) "F2",i,xmxf2(i)
 1204     format('XMX',a2,'_',i4.4,f16.5)
        endif
      enddo

      write(71,'(''SOLARN0'',f19.2)')solarn0
      if (solarn01.ne.0.) then
        write(71,'(''SOLARN0_1'',f19.2)')solarn01
      endif
      write(71,'(''EPHEM'',15x,a)')ephfile(nephem)(1:5)
      if (eclcon.ne."DEFAULT") then
        write(71,'(''ECL'',17x,a)')eclcon(1:index(eclcon," ")-1)
      endif
      write(71,'(''CLK'',17x,a)')clklbl(nclk)
c tempo2-compatibility:
      write(71,'(''UNITS'',15x,''TDB'')')
      write(71,'(''TIMEEPH'',13x,''FB90'')')
      write(71,'(''T2CMETHOD'',11x,''TEMPO'')')
      write(71,'(''CORRECT_TROPOSPHERE'',1x,''N'')')
      write(71,'(''PLANET_SHAPIRO'',6x,''N'')')
      write(71,'(''DILATEFREQ'',10x,''N'')')
      write(71,'(''NTOA'',i22)')ntoa
      write(71,'(''TRES'',f22.2)')tres
      nx = ntzrmjd/10
      fx = (ntzrmjd-10*nx)+ftzrmjd
      if (fx.lt.0.d0) then
        fx = fx + 10.d0
        nx = nx - 10
      endif
      if (fx.ge.10.d0) then
        fx = fx - 10.d0
        nx = nx + 10
      endif
      write(71,'(''TZRMJD '',i5,f16.14)')nx,fx
      write(71,'(''TZRFRQ '',f19.3)')tzrfrq
      if (nsite.le.0 .or. siten2b(ntzrsite).eq.'--') then
        write(71,'(''TZRSITE '',17x,a)')siten2a(ntzrsite)
      else
        write(71,'(''TZRSITE '',17x,a)')siten2b(ntzrsite)
      endif
      write(71,'(''MODE'',i22)')fitmode
      if(nprnt.gt.0)write(71,'(''NPRNT'',i21)')nprnt
      if(nits.gt.0)write(71,'(''NITS'',i22)')nits
      if(iboot.gt.0)write(71,'(''IBOOT'',i21)')iboot
      if(nddm.gt.0)write(71,'(''NDDM'',i22)')nddm
      if(usedmdata)write(71,'(''DMDATA'',i20)')1
      if(infoflag.ne."")write(71,'(''INFO'',1x,a)')infoflag

      return
      end

c=======================================================================

      subroutine outbinpar

      implicit real*8 (A-H,O-Z)

      parameter (TWOPI=6.28318530717958648d0)

      include 'dim.h'     
      include 'acom.h'
      include 'bcom.h'
      include 'trnsfr.h'
      include 'dp.h'
      include 'orbit.h'

      character fit1*3

      fit1='  1'

      nbinout = nbin
      if (nbin.eq.9.and.usefw10) nbinout=16

      write(71,'(''BINARY'',12x,a)')bmodel(nbinout)
      
      if(nbin.eq.10 .and. nplanets.gt.0)then
        write(71,'(''PLAN'',1x,i2)')nplanets
      endif

      if(nfit(9).gt.0)then
         write(71,1009)a1(1),fit1,ferr(9)
      else
         write(71,1009)a1(1)
      endif
 1009 format('A1',f24.9,a,f20.9)

      if(nbin.ne.9)then
         if(nfit(10).gt.0)then
            write(71,1010)e(1),fit1,ferr(10)
         else
            write(71,1010)e(1)
         endif
      endif
 1010 format('E',f25.10,a,f20.10)

      if(nbin.ne.9)then
         if(nfit(11).gt.0)then
            write(71,1011)t0(1),fit1,ferr(11)
         else
            write(71,1011)t0(1)
         endif
      endif
 1011 format('T0',f24.10,a,f20.10)

      if (nbin.ne.10 .and. ((.not.usefb).or.nbin.ne.9)) then
        if(nfit(12).gt.0)then
          write(71,1012)pb(1),fit1,ferr(12)
        else
          write(71,1012)pb(1)
        endif
      endif
 1012 format('PB',f24.14,a,f20.14)

      if(nbin.ne.9)then
         if(nfit(13).gt.0)then
            write(71,1013)omz(1),fit1,ferr(13)*57.295
         else
            write(71,1013)omz(1)
         endif
      endif
 1013 format('OM',f24.12,a,f20.12)

      if(nbin.eq.9)then
         if(nfit(11).gt.0)then
            write(71,2011)t0asc,fit1,ferr(11)
         else
            write(71,2011)t0asc
         endif
 2011    format('TASC',f22.9,a,f20.9)

         if(nfit(10).gt.0)then
            write(71,2010)eps1,fit1,ferr(10)
         else
            write(71,2010)eps1
         endif
 2010    format('EPS1',f22.10,a,f20.10)

         if(nfit(13).gt.0)then
            write(71,2013)eps2,fit1,ferr(13)
         else
            write(71,2013)eps2
         endif
 2013    format('EPS2',f22.10,a,f20.10)
      endif

      if(omdot.ne.0.)then
         if(nfit(14).gt.0)then
            write(71,1014)omdot,fit1,ferr(14)
         else
            write(71,1014)omdot
         endif
      endif
 1014 format('OMDOT',f21.7,a,f20.7)

      if(xomdot.ne.0.)then
         if(nfit(37).gt.0)then
            write(71,1037)xomdot,fit1,ferr(37)
         else
            write(71,1037)xomdot
         endif
      endif
 1037 format('XOMDOT',f20.7,a,f20.7)

      if(gamma.ne.0.)then
         if(nfit(15).gt.0)then
            write(71,1015)gamma,fit1,ferr(15)
         else
            write(71,1015)gamma
         endif
      endif
 1015 format('GAMMA',f21.9,a,f20.9)

      if(pbdot.ne.0.)then
         if(nfit(18).gt.0)then
            write(71,1018)pbdot*1.d12,fit1,ferr(18)*1.d6
         else
            write(71,1018)pbdot*1.d12
         endif
      endif
 1018 format('PBDOT',f21.7,a,f20.7)

      do i = 1, NFBMAX
        if (fb(i).ne.0.or.ferr(NPAR3+i).ne.0) then
         if(nfit(NPAR3+i).gt.0)then
           if (i-1.lt.10) then
             write(71,1070)i-1,fb(i)/(FBFAC**i),fit1,
     +                               ferr(NPAR3+i)/(FBFAC**i)
           else
             write(71,1071)i-1,fb(i)/(FBFAC**i),fit1,
     +                               ferr(NPAR3+i)/(FBFAC**i)
           endif
         else
           if (i-1.lt.10) then
             write(71,1070)i-1,fb(i)/(FBFAC**i)
           else
             write(71,1071)i-1,fb(i)/(FBFAC**i)
           endif
         endif
        endif
      enddo
 1070 format('FB',i1,1p,d23.12,a,d20.12)
 1071 format('FB',i2,1p,d22.12,a,d20.12)

      do i = 1, nfbj
        if (tfbj(i).ne.0.or.ferr(NPAR5+2*i-1).ne.0) then
          if (ferr(NPAR5+2*i-1).gt.0)then
            write (71,1045)i,tfbj(i),fit1,ferr(NPAR5+2*i-1)
          else
            write (71,1045)i,tfbj(i)
          endif
        endif
        if (fbj(i).ne.0.or.ferr(NPAR5+2*i).ne.0) then
          if (ferr(NPAR5+2*i).gt.0)then
            write (71,1046)i,fbj(i),fit1,ferr(NPAR5+2*i)
          else
            write (71,1046)i,fbj(i)
          endif
        endif
      enddo
 1045 format('TFBJ_',z1,f20.6,a,f20.6)
 1046 format('FBJ_',z1,1p,d21.12,a,d20.12)
      
      if(xpbdot.ne.0.)then
         if(nfit(38).gt.0)then
            write(71,1038)xpbdot*1.d12,fit1,ferr(38)*1.d6
         else
            write(71,1038)xpbdot*1.d12
         endif
      endif
 1038 format('XPBDOT',f20.7,a,f20.7)

c --- params 20 and 22 have different meaning depending on Shapiro
c --- parameterization (standard vs FW10):

      if(usefw10)then
          if(varsigma.ne.0.)then
              if(nfit(20).gt.0)then
                  write(71,1120)varsigma,fit1,ferr(20)
              else
                  write(71,1120)varsigma
              endif
          endif
1120      format('VARSIGMA',f22.6,a,f20.6)
          if(h3.ne.0.)then
              if(nfit(22).gt.0)then
                  write(71,1122)h3,fit1,ferr(22)
              else
                  write(71,1122)h3
              endif
          endif
1122      format('H3',f24.12,a,f20.12)
          if(h4.ne.0.)then
              if(nfit(20).gt.0)then
                  write(71,1124)h4,fit1,ferr(20)
              else
                  write(71,1124)h4
              endif
          endif
1124      format('H4',f24.12,a,f20.12)
      else
          if((si.ne.0.).and.(shapmax.eq.0))then
             if(nfit(20).gt.0)then
                write(71,1020)si,fit1,ferr(20)
             else
                write(71,1020)si
             endif
          endif
1020      format('SINI',f22.6,a,f20.6)

          if(shapmax.ne.0.)then
             if(nfit(20).gt.0)then
                write(71,1099) shapmax,fit1,ferr(20)
             else
                write(71,1099) shapmax
             endif
          endif
1099      format('SHAPMAX',f19.6,a,f20.6)

c         --> NW: higher order Shapiro
          if(nshapho.eq.1.)then
             if(nfit(39).gt.0)then
                    write(71,1199) shaphof,fit1,ferr(39)
                 else
                    write(71,1199) shaphof
                 endif

                 if(cotchi0.ne.0.)then
                write(71,1399) cotchi0
             endif
          endif
1199      format('SHAPHOF',f19.6,a,f20.6)
1399      format('COTCHI0',f19.6)

          if(am.ne.0.)then
             if(nfit(21).gt.0)then
                write(71,1021)am,fit1,ferr(21)
             else
                write(71,1021)am
             endif
          endif
1021      format('MTOT',f22.6,a,f20.6)

          if((am2.ne.0.).and.(h3.eq.0))then
             if(nfit(22).gt.0)then
                write(71,1022)am2,fit1,ferr(22)
             else
                write(71,1022)am2
             endif
          endif
1022      format('M2',f24.6,a,f20.6)
      endif

      if(dth.ne.0.)then
         if(nfit(23).gt.0)then
            write(71,1023)dth*1.d6,fit1,ferr(23)*1.d6
         else
            write(71,1023)dth*1.d6
         endif
      endif
 1023 format('DTHETA',f20.6,a,f20.6)

      if(xdot.ne.0.)then
         if(nfit(24).gt.0)then
            write(71,1024)xdot*1.d12,fit1,ferr(24)*1.d12
         else
            write(71,1024)xdot*1.d12
         endif
      endif
 1024 format('XDOT',f22.6,a,f20.6)

      do i = 2, NXDOTMAX
        if (xdot2(i).ne.0.or.ferr(NPAR4+(i-1)).ne.0) then
          if(nfit(NPAR4+(i-1)).gt.0)then
            write(71,1025)i,xdot2(i),fit1,ferr(NPAR4+(i-1))
          else
            write(71,1025)i,xdot2(i)
          endif
        endif
      enddo
 1025 format('XDOT',z1,1p,d21.12,a,d20.12)
      
      if(edot.ne.0.)then
         if(nfit(25).gt.0)then
            write(71,1030)edot*1.d12,fit1,ferr(25)*1.d12
         else
            write(71,1030)edot*1.d12
         endif
      endif
 1030 format('EDOT',f22.6,a,f20.6)

      do i = 2, NEDOTMAX
        if (edot2(i).ne.0.or.ferr(NPAR7+(i-1)).ne.0) then
          if(nfit(NPAR7+(i-1)).gt.0)then
            write(71,1042)i,edot2(i),fit1,ferr(NPAR7+(i-1))
          else
            write(71,1042)i,edot2(i)
          endif
        endif
      enddo
 1042 format('EDOT',z1,1p,d21.12,a,d20.12)
      
      do i = 2, NOMDOTMAX
        if (omdot2(i).ne.0.or.ferr(NPAR8+(i-1)).ne.0) then
          if(nfit(NPAR8+(i-1)).gt.0)then
            write(71,1041)i,omdot2(i),fit1,ferr(NPAR8+(i-1))
          else
            write(71,1041)i,omdot2(i)
          endif
        endif
      enddo
 1041 format('OMDOT',z1,1p,d20.12,a,d20.12)
      
      if(nbin.eq.8)then
         if(om2dot.ne.0.)then
            if(nfit(39).gt.0)then
               write(71,1039)om2dot,fit1,ferr(39)
            else
               write(71,1039)om2dot
            endif
         endif
 1039    format('OM2DOT',10x,e10.4,a,10x,e10.4)
         
         if(x2dot.ne.0)then
            if(nfit(40).gt.0)then
               write(71,1040)x2dot,fit1,ferr(40)
            else
               write(71,1040)x2dot
            endif
         endif   
 1040    format('X2DOT',11x,e10.4,a,10x,e10.4)
      endif

      if(nbin.eq.9)then
         if(eps1dot.ne.0.)then
            if(nfit(39).gt.0)then
               write(71,2039)eps1dot*1.d12,fit1,ferr(39)*1.d12
            else
               write(71,2039)eps1dot*1.d12
            endif
         endif   
 2039    format('EPS1DOT',f19.6,a,f20.6)

         if(eps2dot.ne.0.)then
            if(nfit(40).gt.0)then
               write(71,2040)eps2dot*1.d12,fit1,ferr(40)*1.d12
            else
               write(71,2040)eps2dot*1.d12
            endif
         endif   
 2040    format('EPS2DOT',f19.6,a,f20.6)
      endif

      if(nbin.eq.14)then
         okom_deg = okom * 360.0d0 / twopi
         okom_err_deg = ferr(52) * 360.0d0 / twopi
         if(nfit(52).gt.0)then
            write(71,2100)90.d0-okom_deg,fit1,okom_err_deg
         else
            write(71,2100)90.d0-okom_deg
         endif
 2100    format('KOM',f23.3,a,f20.3)

         okin_deg = okin * 360.0d0 / twopi
         okin_err_deg = ferr(53) * 360.0d0 / twopi
         if(nfit(53).gt.0)then
            write(71,2150)180.d0-okin_deg,fit1,okin_err_deg
         else
            write(71,2150)180.d0-okin_deg
         endif
 2150    format('KIN',f23.3,a,f20.3)

         if (k96) then
            write(71,'("K96",i26)')1
         endif
      endif

      if(afac.ne.0.)write (71,2200) 'AFAC  ',afac
      if(dr.ne.0.)  write (71,2200) 'DR    ',dr*1.d6
      if(a0.ne.0.)  write (71,2200) 'A0    ',a0*1.d6
      if(b0.ne.0.)  write (71,2200) 'B0    ',b0*1.d6
      if(bp.ne.0.)  write (71,2200) 'BP    ',bp
      if(bpp.ne.0.) write (71,2200) 'BPP   ',bpp
 2200 format(a6,f20.7)

      if(useannorb) write(71,'(''PAASCNODE'',f19.2)')PAAscNode

      if(nplanets.gt.0)then
         do j=2,nplanets+1
            jj=16+j*5
            if(nfit(jj).gt.0)then
               write(71,3026)j,a1(j),fit1,ferr(jj)
            else
               write(71,3026)j,a1(j)
            endif
 3026       format('A1_',i1,f22.9,a,f20.9)
            if(nfit(jj+1).gt.0)then
               write(71,3027)j,e(j),fit1,ferr(jj+1)
            else
               write(71,3027)j,e(j)
            endif
 3027       format('E_',i1,f23.9,a,f20.9)
            if(nfit(jj+2).gt.0)then
               write(71,3028)j,t0(j),fit1,ferr(jj+2)
            else
               write(71,3028)j,t0(j)
            endif
 3028       format('T0_',i1,f22.9,a,f20.9)
            if(nfit(jj+3).gt.0)then
               write(71,3029)j,pb(j),fit1,ferr(jj+3)
            else
               write(71,3029)j,pb(j)
            endif
 3029       format('PB_',i1,f22.12,a,f20.12)
            if(nfit(jj+4).gt.0)then
               write(71,3030)j,omz(j),fit1,ferr(jj+4)*57.295
            else
               write(71,3030)j,omz(j)
            endif
 3030       format('OM_',i1,f22.12,a,f20.12)
         enddo
      endif

      return
      end

c=======================================================================

      subroutine outjumppar

      implicit real*8 (A-H,O-Z)

      parameter (TWOPI=6.28318530717958648d0)

      include 'dim.h'     
      include 'acom.h'
      include 'bcom.h'
      include 'trnsfr.h'
      include 'dp.h'
      include 'orbit.h'

      character fit1*3

      fit1='  1'

      do i = 1, nxoff
        if (nfit(NPAR2+i).gt.0) then 
          fit1='  1'
	else
	  fit1='  0'
        endif
        if (i.le.nflagjumps) then
          write (71,1093) trim(jumpflag(i)),trim(jumpflagval(i)),
     +      dct(i),fit1,ferr(NPAR2+i)/f0
        else if (i.lt.10) then
          write (71,1090) i,dct(i),fit1,ferr(NPAR2+i)/f0
	else if (i.lt.100) then
          write (71,1091) i,dct(i),fit1,ferr(NPAR2+i)/f0
        else
          write (71,1092) i,dct(i),fit1,ferr(NPAR2+i)/f0
        endif
 1090   format('JUMP_',i1,f20.9,a,f20.9)
 1091   format('JUMP_',i2,f19.9,a,f20.9)
 1092   format('JUMP_',i3,f18.9,a,f20.9)
 1093   format('JUMP ',a,' ',a,f18.9,a,f18.9)
      enddo

      return
      end

c=======================================================================
     
      subroutine outerrpar

      implicit real*8 (A-H,O-Z)

      include 'dim.h'     
      include 'acom.h'

      if (useglsfit.and.dcovfile.ne."") then
        write(71,'(''DCOVFILE'',12x,a)')
     +          dcovfile(1:index(dcovfile," ")-1)
        write(71,'(a)') "# NOTE: TOA uncertainties and covariances were"
        write(71,'(a)') "# read from the DCOVFILE entry listed above.  The"
        write(71,'(a)') "# following noise parameters were present in"
        write(71,'(a)') "# the input par file, but were ignored in this"
        write(71,'(a)') "# fit.  They are given here for reference but"
        write(71,'(a)') "# may not be consistent with the DCOVFILE"
        write(71,'(a)') "# contents!"
        write(71,'(a)') "# ------- begin noise params -------"
      endif

      if (useglsfit) then
        if (rnidx.ne.0) then
          write(71,'(''RNAMP'',d20.5)') rnamp
          write(71,'(''RNIDX'',f20.5)') rnidx
        endif
      endif

      do i=1,nflagefac
        write(71,1094) trim(efacflag(i)),trim(efacflagval(i)),
     +    flagefac(i)
1094    format('T2EFAC ',a,' ',a,' ',f7.3)
      enddo

      do i=1,nflagequad
        write(71,1095) trim(equadflag(i)),trim(equadflagval(i)),
     +    flagequad(i)
1095    format('T2EQUAD ',a,' ',a,' ',f9.5)
      enddo

      if (useglsfit) then
        do i=1,nflagecorr
          write(71,1096) trim(ecorrflag(i)),trim(ecorrflagval(i)),
     +      flagecorr(i)
1096      format('ECORR ',a,' ',a,' ',f9.5)
        enddo
      endif

      if (useglsfit) then
        do i=1,nflagdmefac
          write(71,1097) trim(dmefacflag(i)),trim(dmefacflagval(i)),
     +      flagdmefac(i)
1097      format('DMEFAC ',a,' ',a,' ',f9.5)
        enddo
      endif

      if (useglsfit) then
        do i=1,nflagdmequad
          write(71,1098) trim(dmequadflag(i)),trim(dmequadflagval(i)),
     +      flagdmequad(i)
1098      format('DMEQUAD ',a,' ',a,' ',f9.5)
        enddo
      endif

      if (useglsfit) then
        do i=1,nflagdmjump
          write(71,1099) trim(dmjumpflag(i)),trim(dmjumpflagval(i)),
     +      flagdmjump(i)
1099      format('DMJUMP ',a,' ',a,' ',f9.5)
        enddo
      endif

      if (useglsfit.and.dcovfile.ne."") then
        write(71,'(a)') "# -------- end noise params -------"
      endif

      return
      end

