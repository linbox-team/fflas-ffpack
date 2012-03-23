# Copyright (c) 2011 FFLAS-FFPACK
# written by BB <bboyer@imag.fr>
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

SUBDIRS =

BENCHMARK=check-dgemm  check-dtrsm  check-dtrtri check-dgetri check-dgetrf

bench: $(BENCHMARK)

AM_CPPFLAGS=-I$(top_srcdir)
AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS) -I$(top_srcdir)/fflas-ffpack/utils/ -I$(top_srcdir)/fflas-ffpack/fflas/  -I$(top_srcdir)/fflas-ffpack/ffpack  -I$(top_srcdir)/fflas-ffpack/field  $(CBLAS_FLAG) $(GMP_CFLAGS)  $(GIVARO_CFLAGS)
LDADD = $(BLAS_LIBS) $(GIVARO_LIBS) $(GMP_LIBS)
AM_LDFLAGS=-static



#  EXTRA_DIST=

CLEANFILES=           \
		$(BENCHMARK)

EXTRA_PROGRAMS =       \
		$(BENCHMARK)

check_dgemm_SOURCES = check-dgemm.C
check_dtrsm_SOURCES = check-dtrsm.C
check_dtrtri_SOURCES = check-dtrtri.C
check_dgetri_SOURCES = check-dgetri.C
check_dgetrf_SOURCES = check-dgetrf.C

# for compilation of new benches
FFLASFFPACK_BIN=@bindir@
%:%.C
	    $(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)
