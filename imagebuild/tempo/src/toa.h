c      $Id$
c
c   toa.h: 
c
c   Variables for storing raw TOAs and associated information, 
c   as read in from a data file.
c
c   Initial use is to allow storage in memory of fake TOAs in
c   tz mode, an application that could probably be implemented
c   using existing variables.  However, the stored raw TOA
c   mechanism could also be used to store TOAs read from files,
c   reducing the need to re-read files when iterating.  To allow
c   for this future expansion, arrays are defined for all
c   quantities allowed on a TOA line, even those (like ddm) which
c   are not needed in tz mode

        logical stflag          ! flag, .TRUE. means there are stored TOAs
        integer stntoa          ! number of stored TOAs

        integer stnsite(NPTSDEF) ! observing site, integer
        integer stnmjd(NPTSDEF) ! integer part of TOA
        real*8 stfmjd(NPTSDEF)  ! real part of TOA
        real*8 stfrq(NPTSDEF)   ! observing frequency, MHz
        real*8 sterr(NPTSDEF)   ! error in microseconds
        real*8 stddm(NPTSDEF)   ! dm correction

        character*640 stflags(NPTSDEF) ! TOA flags

c For "DM data" mode:
        real*8 dmres(NPTSDEF)  ! Current DM "residuals"
        real*8 dmerr(NPTSDEF)  ! input DM uncertainties per-TOA
        
        common /sttoa/ stfmjd, stfrq, sterr, stddm, stnmjd,
     +       stnsite, stntoa, stflag, stflags, dmres, dmerr
