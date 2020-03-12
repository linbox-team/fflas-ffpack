/* fflas.h
 * Copyright (C) 2005,2013,2014 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file fflas.h
 * @author Clément Pernet.
 * @brief <b>F</b>inite <b>F</b>ield <b>L</b>inear <b>A</b>lgebra <b>S</b>ubroutines
 */

#ifndef __FFLASFFPACK_fflas_H
#define __FFLASFFPACK_fflas_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/config.h"
#include "fflas-ffpack/config-blas.h"
#include <cmath>
#include <cstring>

#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

// namespace FFLAS {
#ifndef WINOTHRESHOLD
#define WINOTHRESHOLD __FFLASFFPACK_WINOTHRESHOLD
#endif
// }

/** Thresholds determining which floating point representation to use, depending
 * on the cardinality of the finite field. This is only used when the element
 * representation is not a floating point type.
 * @bug to be benchmarked.
 */
#ifndef DOUBLE_TO_FLOAT_CROSSOVER
#define DOUBLE_TO_FLOAT_CROSSOVER 800
#endif

#include <float.h>

/// @brief FFLAS: <b>F</b>inite <b>F</b>ield <b>L</b>inear <b>A</b>lgebra <b>S</b>ubroutines.

#include <algorithm>

#include "fflas_enum.h"

#include "fflas-ffpack/utils/fflas_memory.h"
#include "fflas-ffpack/paladin/parallel.h"

//---------------------------------------------------------------------
// Level 1 routines
//---------------------------------------------------------------------
#include "fflas_level1.inl"

//---------------------------------------------------------------------
// Level 2 routines
//---------------------------------------------------------------------
#include "fflas_level2.inl"

//---------------------------------------------------------------------
// Level 3 routines
//---------------------------------------------------------------------
#include "fflas_level3.inl"

#ifdef FFLAS_COMPILED
#include "fflas-ffpack/interfaces/libs/fflas_L1_inst.h"
#include "fflas-ffpack/interfaces/libs/fflas_L2_inst.h"
#include "fflas-ffpack/interfaces/libs/fflas_L3_inst.h"
#endif

//---------------------------------------------------------------------
// Checkers
#include "fflas-ffpack/checkers/checkers_fflas.h"
//---------------------------------------------------------------------

//---------------------------------------------------------------------
// specialisations and implementation
//---------------------------------------------------------------------

#include "fflas_freduce.h"
#include "fflas_fadd.h"
#include "fflas_fscal.h"
#include "fflas_fassign.h"

#include "fflas_fgemm.inl"
#include "fflas_pfgemm.inl"
// fgemm must be before fgemv according to ScalAndReduce function declaration ?!? PG
#include "fflas_fgemv.inl"
#include "fflas-ffpack/paladin/pfgemv.inl"
#include "fflas_freivalds.inl"
#include "fflas_fger.inl"
#include "fflas_fsyrk.inl"
#include "fflas_fsyrk_strassen.inl"
#include "fflas_fsyr2k.inl"
#include "fflas_ftrsm.inl"
#include "fflas_pftrsm.inl"
#include "fflas_ftrmm.inl"
#include "fflas_ftrsv.inl"
#include "fflas_faxpy.inl"
#include "fflas_fdot.inl"

//---------------------------------------------------------------------
// MultiPrecision routines
//---------------------------------------------------------------------

// include multiprecision fields for specialisation

#include "fflas-ffpack/field/rns.h" //forward declaration of the multiprecision field
#include "fflas_fscal_mp.inl"
#include "fflas_freduce_mp.inl"
#include "fflas-ffpack/fflas/fflas_fger_mp.inl"
#include "fflas_fgemm/fgemm_classical_mp.inl"
#include "fflas_ftrsm_mp.inl"
#include "fflas_fgemv_mp.inl"
#include "fflas-ffpack/field/rns.inl" // real implementation of the multiprecision field



#include "fflas-ffpack/paladin/fflas_plevel1.h"

//---------------------------------------------------------------------
// Sparse routines
//---------------------------------------------------------------------

#include "fflas_sparse.h"

//---------------------------------------------------------------------
// Checkers
//---------------------------------------------------------------------
#include "fflas-ffpack/checkers/checkers_fflas.inl"

#endif // __FFLASFFPACK_fflas_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
