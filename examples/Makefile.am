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

AM_CPPFLAGS=-I$(top_srcdir) -g
AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS +=  $(CBLAS_FLAG) $(GIVARO_CFLAGS) $(OPTFLAGS) -I$(top_srcdir)/fflas-ffpack/utils/ -I$(top_srcdir)/fflas-ffpack/fflas/  -I$(top_srcdir)/fflas-ffpack/ffpack  -I$(top_srcdir)/fflas-ffpack/field $(CUDA_CFLAGS) $(PARFLAGS)
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(CUDA_LIBS)
AM_LDFLAGS=-static $(PARLIBS)

FFLA_EXAMP =    2x2-fgemm 
2x2_fgemm_SOURCES = 2x2-fgemm.C

if FFLASFFPACK_HAVE_LAPACK
USE_LAPACK_EXAMP = $(LAPA_EXAMP)
endif

EXAMPLES =			\
        $(FFLA_EXAMP)   \
        $(BLAS_EXAMP)   \
        $(USE_LAPACK_EXAMP) \
        $(USE_OMP_EXAMP)

CLEANFILES = $(EXAMPLES) 

EXTRA_PROGRAMS = $(EXAMPLES)


define other_compilation
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) $(PARFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)
endef

%:%.C
	$(other_compilation)

%:%.cpp
	$(other_compilation)
