# Copyright (c) 2014 FFLAS-FFPACK
# written by JGD <Jean-Guillaume.Dumas@imag.fr>
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
examples:  $(EXAMPLES)

AM_CXXFLAGS = $(FFLASFFPACK_CXXFLAGS) $(GIVARO_CFLAGS) $(BLAS_CFLAGS) $(PARFLAGS) $(PRECOMPILE_FLAGS)
AM_CPPFLAGS = -I$(top_srcdir)
LDADD = $(GIVARO_LIBS) $(BLAS_LIBS) $(PARLIBS)
if FFLASFFPACK_PRECOMPILED
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libfflas.la
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libffpack.la
endif

EXAMPLES =  \
	det		\
	pluq	\
	rank	\
	charpoly\
	matmul 	\
	solve

det_SOURCES = det.C
pluq_SOURCES = pluq.C
rank_SOURCES = rank.C
charpoly_SOURCES = charpoly.C
matmul_SOURCES = matmul.C
solve_SOURCES = solve.C

CLEANFILES = $(EXAMPLES)

EXTRA_PROGRAMS = $(EXAMPLES)

define other_compilation
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)
endef

%:%.C
	$(other_compilation)

%:%.cpp
	$(other_compilation)
