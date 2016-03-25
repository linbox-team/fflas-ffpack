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

AM_CPPFLAGS=-I$(top_srcdir) -g
AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS +=  $(CBLAS_FLAG) $(GIVARO_CFLAGS) $(OPTFLAGS) -I$(top_srcdir)/fflas-ffpack/utils/ -I$(top_srcdir)/fflas-ffpack/fflas/  -I$(top_srcdir)/fflas-ffpack/ffpack  -I$(top_srcdir)/fflas-ffpack/field $(CUDA_CFLAGS) $(PARFLAGS)
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(CUDA_LIBS)
AM_LDFLAGS=-static $(PARLIBS)

PERFPUBLISHERFILE=benchmarks-report.xml

FFLA_BENCH =    benchmark-fgemm benchmark-wino benchmark-ftrsm  benchmark-ftrtri  benchmark-inverse  benchmark-lqup benchmark-pluq benchmark-charpoly benchmark-fgemm-mp benchmark-fgemv-mp benchmark-ftrsm-mp benchmark-lqup-mp
BLAS_BENCH =    benchmark-sgemm$(EXEEXT) benchmark-dgemm benchmark-dtrsm
LAPA_BENCH =    benchmark-dtrtri benchmark-dgetri benchmark-dgetrf


if FFLASFFPACK_HAVE_LAPACK
USE_LAPACK_BENCH = $(LAPA_BENCH)
benchmark_dtrtri_SOURCES = benchmark-dtrtri.C
benchmark_dgetri_SOURCES = benchmark-dgetri.C
benchmark_dgetrf_SOURCES = benchmark-dgetrf.C
endif

BENCHMARKS =			\
        $(FFLA_BENCH)   \
        $(BLAS_BENCH)   \
        $(USE_LAPACK_BENCH) \
        $(USE_OMP_BENCH)

CLEANFILES = $(BENCHMARKS) $(PERFPUBLISHERFILE)

EXTRA_PROGRAMS = $(BENCHMARKS)

benchmark_sgemm_SOURCES = benchmark-dgemm.C
benchmark_dgemm_SOURCES = benchmark-dgemm.C
benchmark_dtrsm_SOURCES = benchmark-dtrsm.C

benchmark_fgemm_SOURCES = benchmark-fgemm.C
benchmark_fgemm_mp_SOURCES = benchmark-fgemm-mp.C
benchmark_fgemv_mp_SOURCES = benchmark-fgemv-mp.C
benchmark_wino_SOURCES = benchmark-wino.C
benchmark_ftrsm_SOURCES = benchmark-ftrsm.C
benchmark_ftrsm_mp_SOURCES = benchmark-ftrsm-mp.C
benchmark_ftrtri_SOURCES = benchmark-ftrtri.C
benchmark_inverse_SOURCES = benchmark-inverse.C
benchmark_charpoly_SOURCES = benchmark-charpoly.C
benchmark_lqup_SOURCES = benchmark-lqup.C
benchmark_lqup_mp_SOURCES = benchmark-lqup-mp.C
benchmark_pluq_SOURCES = benchmark-pluq.C

benchmark_sgemm_CXXFLAGS = $(AM_CXXFLAGS) -D__SGEMM__

# Perfpublisher script interaction - AB 2014/11/17
perfpublisher:
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(BENCHMARKS)" "$(CXX)"

# for compilation of new benchmarks
FFLASFFPACK_BIN=@bindir@

define other_compilation
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) $(PARFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)
endef

%:%.C
	$(other_compilation)

%:%.cpp
	$(other_compilation)
