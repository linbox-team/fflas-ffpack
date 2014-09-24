# Copyright (c) 2011 FFLAS-FFPACK
# written by BB <bboyer@imag.fr>
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


pkgincludesubdir=$(pkgincludedir)/utils

SIMD128= simd128.inl \
	 simd128_double.inl \
	 simd128_float.inl   \
	 simd128_int16.inl \
	 simd128_int32.inl \
	 simd128_int64.inl

SIMD256= simd256.inl \
	 simd256_double.inl \
	 simd256_float.inl   \
	 simd256_int16.inl \
	 simd256_int32.inl \
	 simd256_int64.inl


pkgincludesub_HEADERS=       \
	     align-allocator.h  \
	     args-parser.h  \
	     debug.h  \
	     fflas_memory.h \
	     flimits.h \
	     Makefile.am \
	     Matio.h  \
	     print-utils.h \
	     simd.h \
	     $(SIMD128) \
	     $(SIMD256)


### is this really extra dist ?
### this is wrong in a header only lib
EXTRA_DIST=timer.h


