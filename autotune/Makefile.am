# Copyright (c) 2016 FFLAS-FFPACK
# written by Clément Pernet <Clement.Pernet@imag.fr>
# Philippe LEDENT <philippe.ledent@etu.univ-grenoble-alpes.fr>
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

AM_CXXFLAGS = $(FFLASFFPACK_CXXFLAGS) $(GIVARO_CFLAGS) $(BLAS_CFLAGS) $(PARFLAGS)
AM_CPPFLAGS = -I$(top_srcdir)
LDADD = $(GIVARO_LIBS) $(BLAS_LIBS) $(PARLIBS)

AUTOTUNE_FGEMM = winograd-modular-float winograd-modular-double winograd-modularbalanced-float winograd-modularbalanced-double
AUTOTUNE_PLUQ = pluq
AUTOTUNE_CHARPOLY = charpoly-LUK-ArithProg charpoly-Danilevskii-LUK arithprog
AUTOTUNE_FTRTRI = ftrtri
AUTOTUNE_FSYTRF = fsytrf
AUTOTUNE_FSYRK = fsyrk

AUTOTUNE = $(AUTOTUNE_FGEMM) $(AUTOTUNE_PLUQ) $(AUTOTUNE_CHARPOLY) $(AUTOTUNE_FTRTRI) $(AUTOTUNE_FSYTRF) $(AUTOTUNE_FSYRK)

EXTRA_PROGRAMS=$(AUTOTUNE)
EXTRA_DIST=tune_fgemm.sh tune_pluq.sh tune_charpoly.sh tune_ftrtri.sh tune_fsytrf.sh tune_fsyrk.sh

define run_tune_fgemm
	${top_srcdir}/autotune/tune_fgemm.sh
endef	
define run_tune_pluq
	${top_srcdir}/autotune/tune_pluq.sh
endef	
define run_tune_charpoly
	${top_srcdir}/autotune/tune_charpoly.sh
endef
define run_tune_ftrtri
	${top_srcdir}/autotune/tune_ftrtri.sh
endef
define run_tune_fsytrf
	${top_srcdir}/autotune/tune_fsytrf.sh
endef
define run_tune_fsyrk
	${top_srcdir}/autotune/tune_fsyrk.sh
endef
define merge_thresholds
	cat fgemm-thresholds.h pluq-threshold.h charpoly-LUK-ArithProg-threshold.h charpoly-Danilevskii-LUK-threshold.h arithprog-blocksize.h ftrtri-threshold.h fsytrf-threshold.h fsyrk-threshold.h> ${top_srcdir}/fflas-ffpack/fflas-ffpack-thresholds.h
endef	

# This forces the autotune runs to be sequential
autotune: autotune_fgemm autotune_pluq autotune_charpoly autotune_ftrtri autotune_fsytrf autotune_fsyrk
	${run_tune_fgemm}
	${run_tune_pluq}
	${run_tune_charpoly}
	${run_tune_ftrtri}
	${run_tune_fsyrk}
	${run_tune_fsytrf}
	${merge_thresholds}

# Tune fgemm only
tune_fgemm: autotune_fgemm
	${run_tune_fgemm}
	${merge_thresholds}

# Tune pluq only
tune_pluq: autotune_pluq
	${run_tune_pluq}
	${merge_thresholds}

# Tune ftrtri only
tune_ftrtri: autotune_ftrtri
	${run_tune_ftrtri}
	${merge_thresholds}

# Tune fsytrf only
tune_fsytrf: autotune_fsytrf
	${run_tune_fsytrf}
	${merge_thresholds}

# Tune fsyrk only
tune_fsyrk: autotune_fsyrk
	${run_tune_fsyrk}
	${merge_thresholds}

# Tune charpoly only
tune_charpoly: autotune_charpoly
	${run_tune_charpoly}
	${merge_thresholds}

autotune_fgemm: $(AUTOTUNE_FGEMM:%=%$(EXEEXT))

autotune_pluq: $(AUTOTUNE_PLUQ:%=%$(EXEEXT))

autotune_ftrtri: $(AUTOTUNE_FTRTRI:%=%$(EXEEXT))

autotune_fsytrf: $(AUTOTUNE_FSYTRF:%=%$(EXEEXT))

autotune_fsyrk: $(AUTOTUNE_FSYRK:%=%$(EXEEXT))

autotune_charpoly: $(AUTOTUNE_CHARPOLY:%=%$(EXEEXT))

pluq_SOURCES=pluq.C
ftrtri_SOURCES=ftrtri.C
fsytrf_SOURCES=fsytrf.C
fsyrk_SOURCES=fsyrk.C
charpoly_Danilevskii_LUK_SOURCES=charpoly.C
charpoly_Danilevskii_LUK_CXXFLAGS=$(AM_CXXFLAGS) -DVARIANT1="FfpackDanilevski" -DVARIANT2="FfpackLUK" -DNSTART=32 -DNFIRSTSTEP=32 -DNMAX=1000 -DNPREC=1 -DITER=500
charpoly_LUK_ArithProg_SOURCES=charpoly.C
charpoly_LUK_ArithProg_CXXFLAGS=$(AM_CXXFLAGS) -DVARIANT1="FfpackLUK" -DVARIANT2="FfpackArithProg" -DNSTART=1024 -DNFIRSTSTEP=128 -DNMAX=3000 -DNPREC=8 -DITER=1
arithprog_SOURCES=arithprog.C
#arithprog_CXXFLAGS=$(#AM_CXXFLAGS)
winograd_modular_float_SOURCES=winograd.C
winograd_modular_float_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::Modular<float>"
winograd_modular_double_SOURCES=winograd.C
winograd_modular_double_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::Modular<double>"
winograd_modularbalanced_float_SOURCES=winograd.C
winograd_modularbalanced_float_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::ModularBalanced<float>"
winograd_modularbalanced_double_SOURCES=winograd.C
winograd_modularbalanced_double_CXXFLAGS=$(AM_CXXFLAGS) -DFIELD="Givaro::ModularBalanced<double>"

CLEANFILES = $(AUTOTUNE)

