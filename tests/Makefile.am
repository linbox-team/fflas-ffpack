# Copyright (c) 2011 FFLAS-FFPACK
# written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

SUBDIRS = data
check:
	$(BASE_TESTS)

AM_CPPFLAGS=-I$(top_srcdir)
AM_CXXFLAGS = @TESTS_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS)  -I$(top_srcdir)/fflas-ffpack/ -I$(top_srcdir)/fflas-ffpack/utils/ -I$(top_srcdir)/fflas-ffpack/fflas/  -I$(top_srcdir)/fflas-ffpack/ffpack  -I$(top_srcdir)/fflas-ffpack/field $(GIVARO_CFLAGS) $(CBLAS_FLAG) $(CUDA_CFLAGS) $(PARFLAGS) $(PRECOMPILE_FLAGS)

AM_LDFLAGS=-static  #-L$(prefix)/lib   -lfflas -lffpack -lfflas_c -lffpack_c

EXTRA_DIST= test-utils.h

PERFPUBLISHERFILE=tests-report.xml

BASIC_TESTS =               \
		test-lu             \
		test-det            \
		test-echelon        \
		test-rankprofiles   \
		test-compressQ      \
		test-fadd           \
		test-finit          \
		test-fscal          \
		test-fgemm          \
		test-fger           \
		test-ftrsm          \
		test-multifile      \
		regression-check

if FFLASFFPACK_PRECOMPILED
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(CUDA_LIBS) $(PARFLAGS) \
	$(top_builddir)/fflas-ffpack/interfaces/libs/libfflas.la \
	$(top_builddir)/fflas-ffpack/interfaces/libs/libffpack.la
INTERFACE_TESTS= test-interfaces-c 
test_interfaces_c_LDFLAGS = $(LDADD) \
	$(top_builddir)/fflas-ffpack/interfaces/libs/libfflas_c.la \
	$(top_builddir)/fflas-ffpack/interfaces/libs/libffpack_c.la
else
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(CUDA_LIBS) $(PARFLAGS)
endif
NOT_A_TEST =  \
		test-lqup2             \
		test-charpoly          \
		benchlqup              \
		test-fsquare           \
		test-redcolechelon     \
		benchfgemm             \
		test-rank              \
		test-krylov-elim       \
		test-rowechelon        \
		test-fgemv             \
		test-colechelon        \
		test-fullranksubmatrix \
		test-redrowechelon     \
		test-ftrtri            \
		test-redechelon        \
		test-frobenius         \
		test-fgesv             \
		test-invert            \
		test-nullspace

INFINITE_TEST= \
		testeur_fgemm          \
		testeur_ftrsm          \
		testeur_lqup


EXTRA_PROGRAMS =       \
		$(BASIC_TESTS)  $(USE_OMP_TESTS) $(INTERFACE_TESTS)

CLEANFILES =           \
		$(NOT_A_TEST) $(EXTRA_PROGRAMS) $(PERFPUBLISHERFILE)

TESTS =     $(EXTRA_PROGRAMS)

test_compressQ_SOURCES         = test-compressQ.C
test_lu_SOURCES              = test-lu.C
#test_lqup2_SOURCES              = test-lqup2.C
test_det_SOURCES               = test-det.C
test_echelon_SOURCES           = test-echelon.C
test_rankprofiles_SOURCES           = test-rankprofiles.C
test_fgemm_SOURCES             = test-fgemm.C
test_fger_SOURCES             = test-fger.C
test_multifile_SOURCES             = test-multifile1.C test-multifile2.C
#  test_fgemm_SOURCES             = test-fgemm.C
#  test_charpoly_SOURCES          = test-charpoly.C
#  benchfgemm_SOURCES             = benchfgemm.C
#  test_fsquare_SOURCES           = test-fsquare.C
#  test_rank_SOURCES              = test-rank.C
#  benchlqup_SOURCES              = benchlqup.C
#  test_ftrmm_SOURCES             = test-ftrmm.C
#  test_redcolechelon_SOURCES     = test-redcolechelon.C
#  testeur_fgemm_SOURCES          = testeur_fgemm.C
test_ftrsm_SOURCES             = test-ftrsm.C
#  test_redechelon_SOURCES        = test-redechelon.C
#  testeur_ftrsm_SOURCES          = testeur_ftrsm.C
#  test_ftrtri_SOURCES            = test-ftrtri.C
#  test_redrowechelon_SOURCES     = test-redrowechelon.C
#  testeur_lqup_SOURCES           = testeur_lqup.C
#  test_fullranksubmatrix_SOURCES = test-fullranksubmatrix.C
#  test_rowechelon_SOURCES        = test-rowechelon.C
#  test_invert_SOURCES            = test-invert.C
#  test_fgemv_SOURCES             = test-fgemv.C
#  test_krylov_elim_SOURCES       = test-krylov-elim.C
#  test_colechelon_SOURCES        = test-colechelon.C
#  test_fgesv_SOURCES             = test-fgesv.C
#  test_frobenius_SOURCES         = test-frobenius.C
#  test_nullspace_SOURCES         = test-nullspace.C
test_fadd_SOURCES = test-fadd.C
test_fscal_SOURCES = test-fscal.C
test_finit_SOURCES = test-finit.C
test_interfaces_c_SOURCES = test-interfaces-c.c
#test_interfaces_c_CFLAGS= -std=c11 -I/$(prefix)/include $(AM_CPPFLAGS) $(AM_CXXFLAGS) $(PARFLAGS)
#test_interfaces_c_LDFLAGS= $(LDFLAGS) $(LDADD) $(AM_LDFLAGS) -L/$(prefix)/lib/ -lfflas_c -lffpack_c -lstdc++
#  test_fspmv_SOURCES = test-fspmv.C

regression_check_SOURCES = regression-check.C


dense_generator: dense_generator.C
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) dense_generator.C -o dense_generator

#  dense_generator_SOURCES        = dense_generator.C

# Perfpublisher script interaction - AB 2014/11/17
perfpublisher:
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(EXTRA_PROGRAMS)" "$(CXX)"

# for compilation of new tests
FFLASFFPACK_BIN=@bindir@

new_examp_comp = $(CXX) $(CXXFLAGS) $(AM_CXXFLAGS)  ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)

%:%.C
	$(new_examp_comp)

%:%.cpp
	$(new_examp_comp)
