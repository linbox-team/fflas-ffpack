# Copyright (c) 2011 FFLAS-FFPACK
# This file is part of FFLAS-FFPACK
# written by BB <bboyer@imag.fr>
# see COPYING for licence

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=fflas-ffpack tests doc benchmark macros utils


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

.PHONY:examples benchmarks

bin_SCRIPTS=fflas-ffpack-config

VERSION=1.4.3

#  EXTRA_DIST=incremente-versions
