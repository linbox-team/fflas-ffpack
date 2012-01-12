# Coypright (c) 2011 FFLAS-FFPACK
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
check:
	$(BASE_TESTS)

AM_CPPFLAGS=-I$(top_srcdir)
AM_CXXFLAGS = @TESTS_CFLAGS@
AM_CPPFLAGS += $(OPTFLAGS) -I$(top_srcdir)/fflas-ffpack/utils/ -I$(top_srcdir)/fflas-ffpack/fflas/  -I$(top_srcdir)/fflas-ffpack/ffpack  -I$(top_srcdir)/fflas-ffpack/field  $(CBLAS_FLAG) $(GMP_CFLAGS)  $(GIVARO_CFLAGS)
LDADD = $(BLAS_LIBS) $(GIVARO_LIBS) $(GMP_LIBS)
AM_LDFLAGS=-static

BASIC_TESTS =              \
		test-lqup              \
		test-compressQ

NOT_A_TEST =  \
		test-charpoly          \
		benchlqup              \
		test-fsquare           \
		test-det               \
		test-redcolechelon     \
		benchfgemm             \
		test-rank              \
		test-echelon           \
		test-krylov-elim       \
		test-rowechelon        \
		test-fgemv             \
		test-colechelon        \
		test-fullranksubmatrix \
		test-redrowechelon     \
		test-ftrtri            \
		test-redechelon        \
		test-ftrsm             \
		test-fgemm             \
		test-frobenius         \
		test-fgesv             \
		test-invert            \
		test-nullspace

INFINITE_TEST= \
		testeur_fgemm          \
		testeur_ftrsm          \
		testeur_lqup

EXTRA_DIST= timer.C timer.h test-utils.h

CLEANFILES =           \
		$(BASIC_TESTS) $(NOT_A_TEST)

EXTRA_PROGRAMS =       \
		$(BASIC_TESTS)

TESTS =                \
	    $(BASIC_TESTS)

test_compressQ_SOURCES         = test-compressQ.C
test_lqup_SOURCES              = test-lqup.C
#  test_charpoly_SOURCES          = test-charpoly.C
#  dense_generator_SOURCES        = dense_generator.C
#  test_det_SOURCES               = test-det.C
#  benchfgemm_SOURCES             = benchfgemm.C
#  test_fsquare_SOURCES           = test-fsquare.C
#  test_rank_SOURCES              = test-rank.C
#  benchlqup_SOURCES              = benchlqup.C
#  test_echelon_SOURCES           = test-echelon.C
#  test_ftrmm_SOURCES             = test-ftrmm.C
#  test_redcolechelon_SOURCES     = test-redcolechelon.C
#  testeur_fgemm_SOURCES          = testeur_fgemm.C
#  test_ftrsm_SOURCES             = test-ftrsm.C
#  test_redechelon_SOURCES        = test-redechelon.C
#  testeur_ftrsm_SOURCES          = testeur_ftrsm.C
#  test_ftrtri_SOURCES            = test-ftrtri.C
#  test_redrowechelon_SOURCES     = test-redrowechelon.C
#  testeur_lqup_SOURCES           = testeur_lqup.C
#  test_fullranksubmatrix_SOURCES = test-fullranksubmatrix.C
#  test_rowechelon_SOURCES        = test-rowechelon.C
#  test_fgemm_SOURCES             = test-fgemm.C
#  test_invert_SOURCES            = test-invert.C
#  test_fgemv_SOURCES             = test-fgemv.C
#  test_krylov_elim_SOURCES       = test-krylov-elim.C
#  test_colechelon_SOURCES        = test-colechelon.C
#  test_fgesv_SOURCES             = test-fgesv.C
#  test_frobenius_SOURCES         = test-frobenius.C
#  test_nullspace_SOURCES         = test-nullspace.C

dense_generator:
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) dense_generator.C -o dense_generator


# for compilation of new tests
FFLASFFPACK_BIN=@bindir@
%:%.C
	    $(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES)
