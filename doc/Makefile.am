# Copyright (c) 2011 FFLAS-FFPACK
# This file is part of FFLAS-FFPACK
# written by BB <bboyer@imag.fr>
# see COPYING for licence

if FFLASFFPACK_BUILD_DOC
USE_TARGETS = docs
INSTALL_TARGETS = install-doc
endif
docdir=$(FFLASFFPACK_DOC_PATH)

#man1_MANS = fflas-ffpack-config.1

all all-local: $(USE_TARGETS)

install install-data-local: $(USE_TARGETS) $(INSTALL_TARGETS)

docs:
	sed -i 's/^\\version.*/\\version\ $(VERSION)/' mainpage.doxy
	if test -d fflas-ffpack-html ; then echo exists; else mkdir fflas-ffpack-html ; fi
	if test -d fflas-ffpack-dev-html ; then echo exists; else mkdir fflas-ffpack-dev-html ; fi
	cp ../INSTALL fflas-ffpack-html/
	cp ../COPYING fflas-ffpack-html/
	cp ../AUTHORS fflas-ffpack-html/
	doxygen Doxyfile

docs_dev:
	make docs
	cp ../INSTALL fflas-ffpack-dev-html/
	cp ../COPYING fflas-ffpack-dev-html/
	cp ../AUTHORS fflas-ffpack-dev-html/
	doxygen DoxyfileDev

#  cp tutorial.html fflas-ffpack-dev-html/
#  cp install-dist.html fflas-ffpack-dev-html/
#  cp install-dev.html fflas-ffpack-dev-html/



install-doc:
	mkdir -p $(DESTDIR)/$(docdir)
	cp -rp fflas-ffpack-html     $(DESTDIR)/$(docdir)/fflas-ffpack-html
	cp -rp fflas-ffpack-dev-html $(DESTDIR)/$(docdir)/fflas-ffpack-dev-html
	cp -p fflas-ffpack.html      $(DESTDIR)/$(docdir)/fflas-ffpack.html

EXTRA_DIST= \
			Doxyfile            \
			DoxyfileDev         \
			mainpage.doxy       \
			fflas-ffpack.html
# \
#			doc.doxy            \
#			tutorial.doxy       \
#			fflas-ffpack.html         \
#			fflas-ffpack-config.1 \
#tutorial.html       \
#	install-dev.html    \
#	index-dev.html      \
#	install-dist.html


clean-local :
	rm -rf fflas-ffpack-html
	rm -rf fflas-ffpack-dev-html


