# Copyright (c) 2011 Fflas-Ffpack
# This file is part of Fflas-Ffpack
# written by BB <bboyer@imag.fr>
# see COPYING for licence

ACLOCAL_AMFLAGS = -I macros

SUBDIRS=fflas-ffpack tests doc benchmark macros utils

include_HEADERS=fflasffpack-config.h

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

bin_SCRIPTS=fflasffpack-config

VERSION=1.4.2

#  EXTRA_DIST=incremente-versions
