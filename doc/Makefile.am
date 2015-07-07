# Copyright (c) 2011 FFLAS-FFPACK
# written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
# adapted from LinBox documentation
#
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/


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


