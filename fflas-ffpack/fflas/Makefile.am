# Copyright (c) 2011 FFLAS-FFPACK
# written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
# adapted from LinBox configuration
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


pkgincludesubdir=$(pkgincludedir)/fflas

SUBDIRS=fflas_fgemm fflas_igemm fflas_simd fflas_sparse

sparse=fflas_sparse.h \
       fflas_sparse.inl

multiprecision= fflas_ftrsm_mp.inl \
		fflas_fscal_mp.inl \
		fflas_freduce_mp.inl \
		fflas_fger_mp.inl  \
		fflas_fgemv_mp.inl

pkgincludesub_HEADERS=       \
	   fflas_bounds.inl      \
	   fflas_fassign.h       \
	   fflas_fassign.inl       \
	   fflas_ftrmm.inl       \
	   fflas.h               \
	   fflas_level1.inl \
	   fflas_level2.inl \
	   fflas_level3.inl \
	   fflas_fadd.h        \
	   fflas_fadd.inl        \
	   fflas_fdot.inl        \
	   fflas_ftrmm_src.inl   \
	   fflas_fgemm.inl       \
	   fflas_pfgemm.inl      \
	   fflas_pftrsm.inl      \
	   fflas_ftrsm.inl       \
	   fflas_fsyrk.inl       \
	   fflas_fsyrk_strassen.inl       \
	   fflas_fsyr2k.inl       \
	   fflas_fgemv.inl       \
	   fflas_freivalds.inl       \
	   fflas_fscal.h       \
	   fflas_fscal.inl       \
	   fflas_transpose.h	\
	   fflas_ftrsm_src.inl   \
	   fflas_faxpy.inl       \
	   fflas_fger.inl        \
	   fflas_ftrsv.inl       \
	   fflas_freduce.h         \
	   fflas_freduce.inl       \
	   fflas_helpers.inl     \
	   fflas_simd.h          \
	   fflas_enum.h          \
	   ${sparse}		 \
	   ${multiprecision}

EXTRA_DIST=fflas.doxy
