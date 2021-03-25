# Copyright (c) 2017 FFLAS-FFPACK
# written by JGD <Jean-Guillaume.Dumas@imag.fr>
# Philippe LEDENT <philippe.ledent@etu.univ-grenoble-alpes.fr>
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

tutorials:  $(TUTORIALS)

AM_CXXFLAGS = $(FFLASFFPACK_CXXFLAGS) $(GIVARO_CFLAGS) $(BLAS_CFLAGS) $(PARFLAGS) $(PRECOMPILE_FLAGS)
AM_CPPFLAGS = -I$(top_srcdir)
LDADD = $(GIVARO_LIBS) $(BLAS_LIBS) $(PARLIBS)
if FFLASFFPACK_PRECOMPILED
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libfflas.la
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libffpack.la
endif

TUTORIALS = fflas_101 fflas_101_lvl1 fflas-101_1 fflas-101_3 2x2-fgemm 101-fgemm 2x2-pluq 2x2-ftrsv ffpack-solve ffpack-fgesv
2x2_fgemm_SOURCES = 2x2-fgemm.C
2x2_ftrsv_SOURCES = 2x2-ftrsv.C
101_fgemm_SOURCES = 101-fgemm.C
2x2_pluq_SOURCES = 2x2-pluq.C
fflas_101_SOURCES = fflas_101.C
fflas_101_lvl1_SOURCES = fflas_101_lvl1.C
fflas_101_1_SOURCES = fflas-101_1.C
fflas_101_3_SOURCES = fflas-101_3.C
ffpack_solve_SOURCES = ffpack-solve.C
ffpack_fgesv_SOURCES = ffpack-fgesv.C

CLEANFILES = $(TUTORIALS) 

EXTRA_PROGRAMS = $(TUTORIALS)

define other_compilation
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) ${DEFAULT_INCLUDES} $(PARFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)
endef

%:%.C
	$(other_compilation)

%:%.cpp
	$(other_compilation)
