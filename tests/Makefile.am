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

AM_CXXFLAGS = $(FFLASFFPACK_CXXFLAGS) $(GIVARO_CFLAGS) $(BLAS_CFLAGS) $(PARFLAGS) $(PRECOMPILE_FLAGS) -g
AM_CPPFLAGS = -I$(top_srcdir)
LDADD = $(GIVARO_LIBS) $(BLAS_LIBS) $(PARLIBS)
if FFLASFFPACK_PRECOMPILED
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libfflas.la
LDADD += $(top_builddir)/fflas-ffpack/interfaces/libs/libffpack.la
endif

PERFPUBLISHERFILE=tests-report.xml

BASIC_TESTS =               \
		test-fdot           \
		test-finit          \
		test-fscal          \
		test-fadd           \
		test-fger           \
		test-ftrsv          \
		test-ftrtri         \
		test-ftrmv          \
		test-ftrsm          \
		test-ftrstr          \
		test-ftrssyr2k      \
		test-ftrsm-check    \
		test-ftrmm          \
		test-fgemm          \
		test-fgemm-check    \
		test-permutations   \
		test-rpm   \
		test-compressQ      \
		test-lu             \
		test-quasisep       \
		test-sss            \
		test-pluq-check     \
		test-fsyrk          \
		test-fsyr2k          \
		test-fsytrf          \
		test-invert-check   \
		test-rankprofiles   \
		test-det            \
		test-det-check     \
		test-echelon        \
		test-charpoly       \
		test-charpoly-check \
		test-minpoly		\
		test-multifile      \
		test-io      \
		test-maxdelayeddim \
		test-solve \
		test-fgesv             \
		test-simd \
		test-fgemv \
		test-nullspace \
		test-storage-transpose	\
		regression-check

if FFLASFFPACK_PRECOMPILED
INTERFACE_TESTS= test-interfaces-c
#test_interfaces_c_CFLAGS = $(AM_CXXFLAGS)
test_interfaces_c_LDADD = \
	$(top_builddir)/fflas-ffpack/interfaces/libs/libfflas_c.la \
	$(top_builddir)/fflas-ffpack/interfaces/libs/libffpack_c.la \
	-lm -lstdc++
endif

NOT_A_TEST =  \
		test-krylov-elim       \
		test-fullranksubmatrix \
		test-frobenius         \
		test-invert            \
		test-sparse

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
test_permutations_SOURCES         = test-permutations.C
test_rpm_SOURCES         = test-rpm.C

test_lu_SOURCES              = test-lu.C
test_quasisep_SOURCES              = test-quasisep.C
test_sss_SOURCES               = test-sss.C
test_det_SOURCES               = test-det.C
test_pluq_check_SOURCES = test-pluq-check.C
test_det_check_SOURCES = test-det-check.C
test_fgemm_check_SOURCES = test-fgemm-check.C
test_ftrsm_check_SOURCES = test-ftrsm-check.C
test_invert_check_SOURCES = test-invert-check.C
test_charpoly_check_SOURCES = test-charpoly-check.C
test_minpoly_SOURCES = test-minpoly.C
test_echelon_SOURCES           = test-echelon.C
test_rankprofiles_SOURCES           = test-rankprofiles.C
test_fgemm_SOURCES             = test-fgemm.C
test_fger_SOURCES             = test-fger.C
test_multifile_SOURCES             = test-multifile1.C test-multifile2.C
test_io_SOURCES             = test-io.C
#  test_fgemm_SOURCES             = test-fgemm.C
test_charpoly_SOURCES          = test-charpoly.C
#  test_ftrmm_SOURCES             = test-ftrmm.C
#  testeur_fgemm_SOURCES          = testeur_fgemm.C
test_ftrsm_SOURCES             = test-ftrsm.C
test_ftrstr_SOURCES             = test-ftrstr.C
test_ftrssyr2k_SOURCES             = test-ftrssyr2k.C
test_ftrmm_SOURCES             = test-ftrmm.C
test_ftrsv_SOURCES             = test-ftrsv.C
test_ftrmv_SOURCES             = test-ftrmv.C
test_fsyrk_SOURCES             = test-fsyrk.C
test_fsyr2k_SOURCES             = test-fsyr2k.C
test_fsytrf_SOURCES            = test-fsytrf.C
#  test_redechelon_SOURCES        = test-redechelon.C
#  testeur_ftrsm_SOURCES          = testeur_ftrsm.C
test_ftrtri_SOURCES            = test-ftrtri.C
#  testeur_lqup_SOURCES           = testeur_lqup.C
#  test_fullranksubmatrix_SOURCES = test-fullranksubmatrix.C
#  test_invert_SOURCES            = test-invert.C
#  test_krylov_elim_SOURCES       = test-krylov-elim.C
test_fgesv_SOURCES             = test-fgesv.C
#  test_frobenius_SOURCES         = test-frobenius.C
test_nullspace_SOURCES = test-nullspace.C
test_fdot_SOURCES = test-fdot.C
test_fadd_SOURCES = test-fadd.C
test_fscal_SOURCES = test-fscal.C
test_finit_SOURCES = test-finit.C
#test_sparse_SOURCES = test-sparse.C
test_interfaces_c_SOURCES = test-interfaces-c.c
test_maxdelayeddim_SOURCES = test-maxdelayeddim.C
#  test_fspmv_SOURCES = test-fspmv.C
test_storage_transpose_SOURCES = test-storage-transpose.C
regression_check_SOURCES = regression-check.C

test_solve_SOURCES = test-solve.C
test_simd_SOURCES = test-simd.C
test_fgemv_SOURCES = test-fgemv.C
#test_pfgemm_DSL_SOURCES = test-pfgemm-DSL.C


# Perfpublisher script interaction - AB 2014/11/17
perfpublisher:
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(EXTRA_PROGRAMS)" "$(CXX)"
	mv $(PERFPUBLISHERFILE) ..
	make clean
	mv ../$(PERFPUBLISHERFILE) .

# for compilation of new tests
FFLASFFPACK_BIN=@bindir@

new_examp_comp = $(CXX) $(AM_CXXFLAGS) ${CXXFLAGS}  ${INCLUDES} $(AM_CPPFLAGS) $^ -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(LIBS)

%:%.C
	$(new_examp_comp)

%:%.cpp
	$(new_examp_comp)
