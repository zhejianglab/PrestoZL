c      $Id$
        integer NFLAGMAX 
        parameter(NFLAGMAX=32)
        integer nflag
        character flagname(NFLAGMAX)*80
        character flagvalue(NFLAGMAX)*80

	common/flagcom/ nflag, flagname, flagvalue

