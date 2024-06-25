c      $Id$
	character*80 ephdir,ephfile(NEPHMAX),ephname(NEPHMAX)
        integer ephnamel(NEPHMAX)
        logical ephbigendian(NEPHMAX)
        integer kephem, nephem

        common/ephfile/ ephdir,ephfile,ephname,ephnamel,ephbigendian,
     +                  kephem,nephem

