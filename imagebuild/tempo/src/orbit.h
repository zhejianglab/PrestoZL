c      $Id$

        parameter(NMODELS=16)
        parameter(FBFAC=1.d10)
        ! fb(i) quantities  are scaled by FBFAC**i
        !    input values are multipiled by FBFAC**i
        !    output values arde divided by FBFAC**i
        !    associated elements of fctn() are divided
        !       by FBFAC**i since they are derivatives
        !       with dfb(i) in the denominator

        character bmodel(0:NMODELS)*8

        logical usefb,usefw10

	common/orbitp/ a1(4),e(4),t0(4),pb(4),omz(4),
     +       eps1,eps2,eps1dot,eps2dot,t0asc,okom,okin,
     +       omdot,gamma,pbdot,si,am,am2,dth,xomdot,xpbdot,dr,a0,b0,xk,
     +       bp,bpp,xdot,edot,a0aligned,afac,om2dot,x2dot, shapmax,
     +       varsigma,h3,h4,shaphof, cotchi0,
     +       fb(NFBMAX), xdot2(NXDOTMAX), fbj(NFBJMAX), tfbj(NFBJMAX),
     +       edot2(NEDOTMAX),omdot2(NOMDOTMAX), 
     +       usefb,usefw10

        common/orbitn/ nbin,nplanets,nell1,nfbj,nshapho,k96
	common/bmodel/ bmodel
