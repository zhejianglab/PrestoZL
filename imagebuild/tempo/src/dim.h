c      $Id$

c       NPTSDEF  default maximum number of TOAs
c       NPT      maximum number of clock corrections

c       NFMAX    maximum number of frequency derivatives
c       NGLT     maximum number of glitches
c       NGLP     number of parameters per glitch
c       NJUMP    maximum number of jumps
c       NFBMAX   maximum order of orbital frequency taylor expansion, orbit 1
c       NXDOTMAX maximum number of orbital size (x) derivatives, orbit 1
c       NEDOTMAX maximum number of orbital size (x) derivatives, orbit 1
c       NOMDOTMAX maximum number of orbital size (x) derivatives, orbit 1
c       NFBJMAX  paximum number of orbital period jumps, orbit 1
c       NOBSMAX  maximum number of observatories

c                parameters are in following order:
c       1       to NPAR1  sixty basic parameters
c       NPAR1+1 to NPAR2  glitch parameters (NGLT*NGLP)
c       NPAR2+1 to NPAR3  jump parameters (NJUMP)
c       NPAR3+1 to NPAR4  Orbital frequency & derivatives (1 to NFBMAX)
c       NPAR4+1 to NPAR5  Orbital size (ASINI) derivatives (2 to NXDOTMAX)
c       NPAR5+1 to NPAR6  Orbital period jump parameters
c       NPAR6+1 to NPAR7  DM offsets and first derivatives (1 to 2*NDMXMAX)
c       NPAR7+1 to NPAR8  Orbital E derivatives (2 to NEDOTMAX)
c       NPAR8+1 to NPAR9  Orbital OM derivatives (2 to NOMDOTMAX)
c       NPAR9+1 to NPAR10 DM derivatives (1 to NDMCOFMAX)
c       NPAR10+1 to NPAR11 Non-DM freq-dependent terms (1 to NFDMAX)
c       NPAR11+1 to NPA12 Non-DM freq-dependent terms above F2 (3 to NFMAX) 
c       NPAR12+1 to NPA   XMX values and spectral indexes
c       NPA      total number of parameters (basic+glitch+jump)
c       NPAP1    total number of parameters plus one
c       NBUFDEF  default size of virtual memory buffer (why 35*NPTSMAX???)
c       NBOOTMAX max size of bootstrap Monte Carlo integraions
c       NCLKMAX  max number of clock correction files (obs to UTC, UTC to xx)
c       NEPHMAX  max number of ephemerides
c       NUTMAX   max number of ut1 corrections
c       NFLAGERR max number of flag-based EFAC/EQUAD

	parameter (NPTSDEF=100000)
	parameter (NPT=200000)
        parameter (NFMAX=20)
	parameter (NGLT=9,NGLP=5)
	parameter (NJUMP=250)
	parameter (NFBMAX=20)
	parameter (NXDOTMAX=10)
	parameter (NEDOTMAX=10)
	parameter (NOMDOTMAX=10)
        parameter (NFBJMAX=120)
        parameter (NDMXMAX=500)
        parameter (NXMXMAX=500)
        parameter (NDMCOFMAX=30)
        parameter (NFDMAX=10)
        parameter (NPAR1=60)
        parameter (NPAR2=NPAR1+NGLT*NGLP)
	parameter (NPAR3=NPAR2+NJUMP)
	parameter (NPAR4=NPAR3+(NFBMAX))
        parameter (NPAR5=NPAR4+(NXDOTMAX-1))
	parameter (NPAR6=NPAR5+2*NFBJMAX)
        parameter (NPAR7=NPAR6+2*NDMXMAX)
        parameter (NPAR8=NPAR7+(NEDOTMAX-1))
        parameter (NPAR9=NPAR8+(NOMDOTMAX-1))
        parameter (NPAR10=NPAR9+NDMCOFMAX)
        parameter (NPAR11=NPAR10+NFDMAX)
        parameter (NPAR12=NPAR11+NFMAX-2)
        parameter (NPA=NPAR12+2*NXMXMAX)
        parameter (NPAP1=NPA+1)
	parameter (NPARDEF=28)
        parameter (NBOOTMAX=1024)
        parameter (NTZMAX=1000,NCLKMAX=40,NEPHMAX=10,NUTMAX=6000)
        parameter (NTZARR=1500)
        parameter (NFLAGERR=64)
        parameter (NOBSMAX=500)
