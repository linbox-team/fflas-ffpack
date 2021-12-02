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

AM_CXXFLAGS = $(FFLASFFPACK_CXXFLAGS) $(GIVARO_CFLAGS) $(BLAS_CFLAGS) $(PARFLAGS) $(PRECOMPILE_FLAGS)
AM_CPPFLAGS = -I$(top_srcdir)
LDADD = $(GIVARO_LIBS) $(BLAS_LIBS) $(PARLIBS)
if FFLASFFPACK_PRECOMPILED
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libfflas.la
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libffpack.la
endif

PERFPUBLISHERFILE=benchmarks-report.xml

FFLA_BENCH =    benchmark-fgemm benchmark-fgemm-rns benchmark-wino benchmark-ftrsm  benchmark-fgesv benchmark-ftrsv benchmark-ftrtri  benchmark-inverse benchmark-fsytrf benchmark-fsyrk  benchmark-lqup benchmark-fsyr2k benchmark-pluq benchmark-charpoly benchmark-charpoly-mp benchmark-fgemm-mp benchmark-fgemv-mp benchmark-ftrsm-mp benchmark-lqup-mp benchmark-checkers benchmark-fadd-lvl2 benchmark-fdot benchmark-fgemv benchmark-quasisep benchmark-storage-transpose

BLAS_BENCH =    benchmark-sgemm$(EXEEXT) benchmark-dgemm benchmark-dtrsm
LAPA_BENCH =    benchmark-dtrtri benchmark-dgetri benchmark-dgetrf benchmark-dsytrf


if FFLASFFPACK_HAVE_LAPACK
USE_LAPACK_BENCH = $(LAPA_BENCH)
benchmark_dtrtri_SOURCES = benchmark-dtrtri.C
benchmark_dgetri_SOURCES = benchmark-dgetri.C
benchmark_dgetrf_SOURCES = benchmark-dgetrf.C
benchmark_dsytrf_SOURCES = benchmark-dsytrf.C
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
benchmark_fgemm_rns_SOURCES = benchmark-fgemm-rns.C
benchmark_fgemv_mp_SOURCES = benchmark-fgemv-mp.C
benchmark_wino_SOURCES = benchmark-wino.C
benchmark_ftrsm_SOURCES = benchmark-ftrsm.C
benchmark_fgesv_SOURCES = benchmark-fgesv.C
benchmark_ftrsm_mp_SOURCES = benchmark-ftrsm-mp.C
benchmark_ftrsv_SOURCES = benchmark-ftrsv.C
benchmark_ftrtri_SOURCES = benchmark-ftrtri.C
benchmark_inverse_SOURCES = benchmark-inverse.C
benchmark_fsytrf_SOURCES = benchmark-fsytrf.C
benchmark_fsyr2k_SOURCES = benchmark-fsyr2k.C
benchmark_fsyrk_SOURCES = benchmark-fsyrk.C
benchmark_quasisep_SOURCES = benchmark-quasisep.C
benchmark_charpoly_SOURCES = benchmark-charpoly.C
benchmark_charpoly_mp_SOURCES = benchmark-charpoly-mp.C
benchmark_lqup_SOURCES = benchmark-lqup.C
benchmark_lqup_mp_SOURCES = benchmark-lqup-mp.C
benchmark_pluq_SOURCES = benchmark-pluq.C
benchmark_checkers_SOURCES = benchmark-checkers.C
benchmark_fadd_lvl2_SOURCES = benchmark-fadd-lvl2.C
benchmark_fdot_SOURCES = benchmark-fdot.C
benchmark_fgemv_SOURCES = benchmark-fgemv.C
benchmark_storage_transpose_SOURCES = benchmark-storage-transpose.C

benchmark_sgemm_CXXFLAGS = $(AM_CXXFLAGS) -D__SGEMM__

# Perfpublisher script interaction - AB 2014/11/17
perfpublisher:
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(BENCHMARKS)" "$(CXX)"
	mv $(PERFPUBLISHERFILE) ..
	make clean
	mv ../$(PERFPUBLISHERFILE) .

# for compilation of new benchmarks
FFLASFFPACK_BIN=@bindir@

define other_compilation
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(LIBS)
endef

%:%.C
	$(other_compilation)

%:%.cpp
	$(other_compilation)
