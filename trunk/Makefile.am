# Copyright (c) 2011 FFLAS-FFPACK
# written by BB <bboyer@imag.fr>
# adapted from LinBox configuration
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


ACLOCAL_AMFLAGS = -I macros

SUBDIRS=fflas-ffpack tests doc benchmarks macros optimiser


#  include_HEADERS=fflas-ffpack-config.h
#!@todo add examples dir XXX

docs:doc/fflas-ffpack-html/index.html

doc/fflas-ffpack-html/index.html:
	(cd doc; make docs)

docs_dev:doc/fflas-ffpack-dev-html/index.html

doc/fflas-ffpack-dev-html/index.html:
	(cd doc; make docs_dev)

#  examples:
#  (cd examples; make examples)

benchmarks:
	(cd benchmarks; make benchmarks)

uninstall-hook:
	(test -d "$(includedir)/fflas-ffpack" && rm -rf \
	       	"$(mandir)/man1" \
		"$(mandir)" \
		"$(includedir)/fflas-ffpack/fflas" \
		"$(includedir)/fflas-ffpack/fflas/fflas_fgemm" \
		"$(includedir)/fflas-ffpack/ffpack" \
		"$(includedir)/fflas-ffpack/field" \
		"$(includedir)/fflas-ffpack/utils" \
		"$(includedir)/fflas-ffpack/"  \
		"$(datarootdir)/" \
		) || true

.PHONY:examples benchmarks

bin_SCRIPTS=fflas-ffpack-config

VERSION=2.0.0

#  EXTRA_DIST=incremente-versions
