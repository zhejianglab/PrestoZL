c      $Id$
        character clkfile(NCLKMAX)*80, clklbl(0:NCLKMAX)*32,clkdir*80
	integer ckflag(NPT)

	common/clocks/ tdate(NPT),ckcorr(NPT),jsite(NPT),ndate,nclk,
     +    tbipm(NPT),bipm(NPT),utcclk(NPT),tutc(NPT),ckflag
        common/clkfile/clkdir,clkfile,clklbl
