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


pkgincludesubdir=$(pkgincludedir)/ffpack

multiprecision= ffpack_ludivine_mp.inl \
		ffpack_pluq_mp.inl     \
		ffpack_charpoly_mp.inl \
		ffpack_det_mp.inl


pkgincludesub_HEADERS=        \
		ffpack_charpoly_danilevski.inl        \
		ffpack_charpoly.inl                   \
		ffpack_charpoly_kgfastgeneralized.inl \
		ffpack.h                              \
		ffpack_charpoly_kgfast.inl            \
		ffpack_krylovelim.inl                 \
		ffpack_charpoly_kglu.inl              \
		ffpack_echelonforms.inl               \
		ffpack_ludivine.inl                   \
		ffpack_pluq.inl                       \
		ffpack_ppluq.inl \
		ffpack_frobenius.inl                  \
		ffpack_minpoly.inl \
		ffpack.inl\
		ffpack_invert.inl\
		ffpack_fgesv.inl\
		ffpack_fgetrs.inl\
		ffpack_fsytrf.inl\
		ffpack_permutation.inl\
		ffpack_ftrtr.inl\
		ffpack_rankprofiles.inl\
		ffpack_ftrstr.inl\
		ffpack_ftrssyr2k.inl\
		ffpack_bruhatgen.inl\
		ffpack_sss.inl                        \
		$(multiprecision)


EXTRA_DIST=ffpack.doxy
