bin_PROGRAMS = tempo

EXTRA_PROGRAMS = matrix tdbgen
CLEANFILES = $(EXTRA_PROGRAMS) version.h

tempo_SOURCES = \
        a1ut1f.f a1utcf.f ang.f arrtim.f atimfake.f bigendian.f \
        blkdbat.f bnrybt.f bnrybtx.f bnrydd.f bnryddgr.f bnryddk.f bnryddt.f bnryddfwhiecc.f\
        bnryeh.f bnryell1.f bnrymss.f bnrydds.f bootmc.f bt2ell1.f \
        cfgin.f chebpc.f citem.f clockcor.f covar.f cross.f damoyr.f dot.f \
        earth.f  ephinit.f ephread.f equ2ecl.f fit.f fk4tofk5.f galco.f \
        gasdev.f getcmdline.f getecliptic.f getflags.f getvalue.f \
        inpar.f interp.f lmst.f mass2dd.f matinv.f mjdadd.f mxprt.f  \
        newbin.f newsrc.f newval.f nutation.f obsite.f outpar.f pcard.f pcshft.f \
        pospm_conv.f prcnut.f precession.f propmo.f radian.f ran1.f resid.f \
        setup.f sitea2n.f siten2a.f siten2b.f \
	sort.f space_motion.f ell12bt.f tdbinit.f tdbinit2.f tdbread.f tempo.f \
        timcalc.f tmalloc.f tparin.f tpohdr.f tzfit.f tzinit.f ut1red.f \
        vmemw.f  ztim.f \
        cutil.c \
        acom.h array.h dim.h bcom.h clocks.h dp.h eph.h ephcom.h flagcom.h glitch.h \
	orbit.h tdbcom.h trnsfr.h tz.h vcom.h toa.h 

if HAVE_LAPACK
tempo_SOURCES += glsfit.f plnoise.f
tempo_LDADD = $(LDADD) @LAPACK_LIBS@ @BLAS_LIBS@
glsfit.$(OBJEXT): glsfit.f acom.h bcom.h orbit.h dim.h tz.h vcom.h
else
tempo_SOURCES += glsfit_missing.f
glsfit_missing.$(OBJEXT): glsfit_missing.f
endif

matrix_SOURCES = matrix.f

tdbgen_SOURCES = tdbgen.f tdb1ns.f bigendian.f cutil.c


version.h: get_version_id.sh FORCE
	./get_version_id.sh
FORCE:

EXTRA_DIST = version.h
dist-hook:
	echo "#! /bin/sh" > $(distdir)/get_version_id.sh
	echo "cat <<HERE >version.h" >> $(distdir)/get_version_id.sh
	cat version.h >> $(distdir)/get_version_id.sh
	echo "HERE" >> $(distdir)/get_version_id.sh
	chmod a+x $(distdir)/get_version_id.sh

# dependencies of Tempo subroutines

arrtim.$(OBJEXT): arrtim.f acom.h orbit.h dp.h bcom.h trnsfr.h clocks.h dim.h glitch.h tz.h
atimfake.$(OBJEXT): atimfake.f dim.h acom.h tz.h toa.h
bigendian.$(OBJECTX): bigendiam.f config.h
bnrybt.$(OBJEXT): bnrybt.f dp.h orbit.h dim.h
bnrybtx.$(OBJEXT): bnrybt.f dp.h orbit.h dim.h
bnrydd.$(OBJEXT): bnrydd.f dp.h orbit.h acom.h dim.h trnsfr.h
bnryddfwhiecc.$(OBJEXT): bnryddfwhiecc.f dp.h orbit.h acom.h dim.h trnsfr.h
bnryddgr.$(OBJEXT): bnryddgr.f orbit.h dim.h
bnryddk.$(OBJEXT): bnrydd.f dp.h orbit.h acom.h dim.h trnsfr.h
bnrydds.$(OBJEXT): bnrydds.f dim.h acom.h dp.h orbit.h trnsfr.h
bnryddt.$(OBJEXT): bnryddt.f orbit.h dim.h
bnryeh.$(OBJEXT): bnryeh.f dp.h orbit.h dim.h
bnryell1.$(OBJEXT): bnryell1.f dp.h orbit.h dim.h
bnrymss.$(OBJEXT): bnrymss.f dp.h orbit.h dim.h
bootmc.$(OBJEXT): bootmc.f dim.h
bt2ell1.$(OBJEXT): bt2ell1.f orbit.h dim.h
cfgin.$(OBJEXT): cfgin.f acom.h clocks.h dim.h eph.h tz.h 
clockcor.$(OBJEXT): clockcor.f dim.h acom.h clocks.h
covar.$(OBJEXT): acom.h dim.h
cutil.$(OBJEXT): cutil.c
ell12bt.$(OBJEXT): ell12bt.f dim.h orbit.h
ephinit.$(OBJEXT): ephinit.f ephcom.h
ephread.$(OBJEXT): dim.h ephread.f eph.h ephcom.h 
equ2ecl.$(OBJEXT): equ2ecl.f dim.h acom.h
fit.$(OBJEXT): fit.f acom.h bcom.h orbit.h dim.h tz.h vcom.h
getcmdline.$(OBJEXT): getcmdline.f
getecliptic.$(OBJEXT): getecliptic.f dim.h acom.h eph.h
getflags.$(OBJEXT): getflags.f flagcom.h
getvalue.$(OBJEXT): getvalue.f flagcom.h
inpar.$(OBJEXT): inpar.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h 
matinv.$(OBJEXT): matinv.f dim.h
mjdadd.$(OBJEXT): mjdadd.f
mxprt.$(OBJEXT): mxprt.f dim.h
newbin.$(OBJEXT): newbin.f acom.h orbit.h dim.h
newsrc.$(OBJEXT): newsrc.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h
newval.$(OBJEXT): newval.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h
outpar.$(OBJEXT): outpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h eph.h
pcard.$(OBJEXT): pcard.f acom.h dim.h dp.h
resid.$(OBJEXT): resid.f acom.h bcom.h orbit.h dp.h dim.h glitch.h
setup.$(OBJEXT): setup.f acom.h dim.h version.h
sitea2n.$(OBJEXT): sitea2n.f acom.h dim.h
siten2b.$(OBJEXT): siten2b.f acom.h dim.h
tdbinit.$(OBJEXT): tdbinit.f tdbcom.h
tdbinit2.$(OBJEXT): tdbinit2.f tdbcom.h
tdbread.$(OBJEXT): tdbread.f acom.h tdbcom.h dim.h
tempo.$(OBJEXT): tempo.f acom.h bcom.h dp.h array.h clocks.h config.h eph.h vcom.h orbit.h dim.h trnsfr.h tz.h toa.h
timcalc.$(OBJEXT): timcalc.f acom.h dim.h dp.h trnsfr.h
tmalloc.$(OBJEXT): array.h
tparin.$(OBJEXT): tparin.f acom.h dim.h vcom.h tz.h version.h
tpohdr.$(OBJEXT): tpohdr.f dim.h tz.h
tzfit.$(OBJEXT): tzfit.f dim.h acom.h trnsfr.h tz.h
tzinit.$(OBJEXT): tzinit.f dim.h acom.h tz.h
ut1red.$(OBJEXT): ut1red.f dim.h
vmemw.$(OBJEXT): vmemw.f dim.h vcom.h
ztim.$(OBJEXT): ztim.f trnsfr.h

# dependency of matrix:
matrix.$(OBJEXT): matrix.f dim.h


