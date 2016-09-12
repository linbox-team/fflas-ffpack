# Copyright (c) 2016 FFLAS-FFPACK
# written by Clément Pernet <Clement.Pernet@imag.fr>
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

AUTOTUNE_FGEMM = winograd-modular-float winograd-modular-double winograd-modularbalanced-float winograd-modularbalanced-double

AUTOTUNE = $(AUTOTUNE_FGEMM)

EXTRA_PROGRAMS=$(AUTOTUNE)

autotune: autotune_fgemm
	cat fgemm-thresholds.h > ../fflas-ffpack/fflas-ffpack-thresholds.h

autotune_fgemm: $(AUTOTUNE_FGEMM)
	./tune_fgemm.sh

winograd_modular_float_SOURCES=winograd.C
winograd_modular_float_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::Modular<float>"
winograd_modular_double_SOURCES=winograd.C
winograd_modular_double_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::Modular<double>"
winograd_modularbalanced_float_SOURCES=winograd.C
winograd_modularbalanced_float_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::ModularBalanced<float>"
winograd_modularbalanced_double_SOURCES=winograd.C
winograd_modularbalanced_double_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::ModularBalanced<double>"

CLEANFILES = $(AUTOTUNE)

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CXXFLAGS +=  $(CBLAS_FLAG) $(GIVARO_CFLAGS) $(OPTFLAGS)  $(PARFLAGS)
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(CUDA_LIBS)
