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

AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CXXFLAGS +=  $(CBLAS_FLAG) $(GIVARO_CFLAGS) $(OPTFLAGS)  $(PARFLAGS)
LDADD = $(CBLAS_LIBS) $(GIVARO_LIBS) $(CUDA_LIBS)

AUTOTUNE_FGEMM = winograd-modular-float winograd-modular-double winograd-modularbalanced-float winograd-modularbalanced-double
AUTOTUNE_PLUQ = pluq

AUTOTUNE = $(AUTOTUNE_FGEMM) $(AUTOTUNE_PLUQ)

EXTRA_PROGRAMS=$(AUTOTUNE)
EXTRA_DIST=tune_fgemm.sh tune_pluq.sh

define run_tune_fgemm
	./tune_fgemm.sh
endef	
define run_tune_pluq
	./tune_pluq.sh
endef	
define merge_thresholds
	cat fgemm-thresholds.h pluq-threshold.h > ../fflas-ffpack/fflas-ffpack-thresholds.h
endef	

# This forces the autotune runs to be sequential
autotune: autotune_fgemm autotune_pluq
	${run_tune_fgemm}
	${run_tune_pluq}
	${merge_thresholds}

# Tune fgemm only
tune_fgemm: autotune_fgemm
	${run_tune_fgemm}
	${merge_thresholds}

# Tune pluq only
tune_pluq: autotune_pluq
	${run_tune_pluq}
	${merge_thresholds}

autotune_fgemm: $(AUTOTUNE_FGEMM)

autotune_pluq: $(AUTOTUNE_PLUQ)



pluq_SOURCES=pluq.C
winograd_modular_float_SOURCES=winograd.C
winograd_modular_float_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::Modular<float>"
winograd_modular_double_SOURCES=winograd.C
winograd_modular_double_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::Modular<double>"
winograd_modularbalanced_float_SOURCES=winograd.C
winograd_modularbalanced_float_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::ModularBalanced<float>"
winograd_modularbalanced_double_SOURCES=winograd.C
winograd_modularbalanced_double_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::ModularBalanced<double>"

CLEANFILES = $(AUTOTUNE)

