
INSTALLDIR = /usr/local/apache/htdocs/tempo


REFINSTALLDIR = $(INSTALLDIR)/ref_man_sections


REF_MAN_SECTIONS = description.txt \
			options.txt \
			inputfile.txt \
			freeformat.txt \
			control.txt \
			parameter.txt \
			fixedformat.txt \
			binary.txt \
			toa.txt \
			toacommand.txt \
			output.txt \
			config.txt \
			clock.txt \
			ephem.txt \
			tdb-tt.txt \
			ecliptic.txt \
			obsys.txt \
			tz.txt \
			tz-in.txt \
			tz-param.txt \
			tz-polyco.txt \
			notes.txt 

TPO_PAGES = tempo11.html tempo11.000_install.html tempo_install.html \
		tempo_idx.html tempo_update.html reference_manual.html

REFDIR = ref_man_sections

EXTRA_DIST = $(TPO_PAGES)
SUBDIRS = $(REFDIR)

install:
	if [ ! -d $(INSTALLDIR) ] ; then mkdir $(INSTALLDIR) ; fi
	if [ ! -d $(REFINSTALLDIR) ] ; then mkdir $(REFINSTALLDIR) ; fi
	echo "" > $(REFINSTALLDIR)/usage.txt
	cd $(REFDIR); cat $(REF_MAN_SECTIONS) >> $(REFINSTALLDIR)/usage.txt
	cd $(REFDIR); cp $(REF_MAN_SECTIONS) $(REFINSTALLDIR)
	cp $(TPO_PAGES) $(INSTALLDIR)
	chmod ug+w,o+rx $(INSTALLDIR)
	chmod ug+w,o+rx $(REFINSTALLDIR)
	cd $(REFINSTALLDIR); chmod ug+w,o+r *
	cd $(INSTALLDIR); chmod ug+w,o+r $(TPO_PAGES)
	if [ ! -r $(INSTALLDIR)/index.html ] ; \
	  then ln -s $(INSTALLDIR)/tempo_idx.html $(INSTALLDIR)/index.html ; \
          fi

archive:
	tar -cvf tpowww.tar Makefile $(TPO_PAGES) $(REFDIR)/*.txt
	gzip tpowww.tar
