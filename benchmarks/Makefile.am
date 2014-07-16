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
benchmarks:  $(BENCHMARKS)

AM_CPPFLAGS=-I$(top_srcdir)
AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS)  -I$(top_srcdir)/fflas-ffpack/utils/ -I$(top_srcdir)/fflas-ffpack/fflas/  -I$(top_srcdir)/fflas-ffpack/ffpack  -I$(top_srcdir)/fflas-ffpack/field  $(CBLAS_FLAG) $(GMP_CFLAGS)  $(GIVARO_CFLAGS) $(CUDA_CFLAGS) $(PARFLAGS)
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(GMP_LIBS) $(CUDA_LIBS)
AM_LDFLAGS=-static $(PARLFLAGS)

FFLA_BENCH =    benchmark-fgemm benchmark-wino benchmark-pfgemm  benchmark-ftrsm  benchmark-ftrtri  benchmark-inverse  benchmark-lqup 
BLAS_BENCH =    benchmark-dgemm  benchmark-dtrsm
LAPA_BENCH =    benchmark-dtrtri benchmark-dgetri benchmark-dgetrf


OMP_BENCH = benchmark-pfgemm benchmark-pftrsm

if FFLASFFPACK_HAVE_OMP
USE_OMP_BENCH = $(OMP_BENCH)
endif

if FFLASFFPACK_HAVE_LAPACK
USE_LAPACK_TESTS = $(LAPA_BENCH)
endif



BENCHMARKS =			\
        $(FFLA_BENCH)   \
        $(BLAS_BENCH)   \
        $(USE_LAPACK_TESTS) \
        $(USE_OMP_BENCH)

CLEANFILES = $(BENCHMARKS)

EXTRA_PROGRAMS = $(BENCHMARKS)

benchmark_dgemm_SOURCES = benchmark-dgemm.C
benchmark_dtrsm_SOURCES = benchmark-dtrsm.C

benchmark_dtrtri_SOURCES = benchmark-dtrtri.C
benchmark_dgetri_SOURCES = benchmark-dgetri.C
benchmark_dgetrf_SOURCES = benchmark-dgetrf.C

benchmark_pfgemm_SOURCES = benchmark-pfgemm.C
benchmark_pftrsm_SOURCES = benchmark-pftrsm.C

benchmark_fgemm_SOURCES = benchmark-fgemm.C
benchmark_wino_SOURCES = benchmark-wino.C
benchmark_ftrsm_SOURCES = benchmark-ftrsm.C
benchmark_ftrtri_SOURCES = benchmark-ftrtri.C
benchmark_inverse_SOURCES = benchmark-inverse.C
benchmark_lqup_SOURCES = benchmark-lqup.C



# for compilation of new benchmarks
FFLASFFPACK_BIN=@bindir@

define other_compilation
	    $(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) $(PARFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)
endef

%:%.C
	$(other_compilation)

%:%.cpp
	$(other_compilation)
