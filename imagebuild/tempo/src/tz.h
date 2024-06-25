c      $Id$
       character tzdir*80,tzfile*80,tztfile*80,tzsite*80
       character tzrsite*80,tzrsitedefault*80
       character tzsitedef*1
       character polycofile*160
       character name(NTZMAX)*12,pname*12,params(6)*80
       integer ntzsite,ntzrsite
       integer nsp(NTZMAX),nco(NTZMAX)
       real*8 mxha(NTZMAX)
       integer ntmjd(NTZARR)
       real*8 ftmjd(NTZARR),tmin(NTZARR),tzof(NTZMAX)
       real*8 ftzrmjd,ftzrmjddefault
       integer ntzrmjd,ntzrmjddefault
       real*8 tzmjdstart
       integer lupolyco

c      non-character variables:
       common/tz/tzrfrq,tzrfrqdefault,ftzrmjd,ftzrmjddefault,ftmjd,tmin,tzof,
     +    tzmjdstart,
     +    mxha,ntzref,nsp,nco,nsets,params,nsite,
     +    ntzrmjd,ntzrmjddefault,ntmjd,
     +    lupolyco,ntzsite,ntzrsite
c      character variables:
       common/tzch/tzdir,tzfile,tztfile,tzrsite,tzrsitedefault,
     +    tzsite,name,pname,
     +    tzsitedef,polycofile

