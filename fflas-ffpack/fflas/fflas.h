/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#include "fflas_enum.h"

namespace FFLAS{ namespace Protected {

		template <class X, class Y> class AreEqual {
		public:
			static const bool value = false;
		};

		template <class X> class AreEqual<X, X> {
		public:
			static const bool value = true;
		};
	} // Protected
} // class FFLAS
#include <algorithm>

namespace FFLAS {

template <class T> const T &min3(const T &m, const T &n, const T &k) { return std::min(m, std::min(n, k)); }

template <class T> const T &max3(const T &m, const T &n, const T &k) { return std::max(m, std::min(n, k)); }

template <class T> const T &min4(const T &m, const T &n, const T &k, const T &l) {
    return std::min(std::min(m, n), std::min(k, l));
}

template <class T> const T &max4(const T &m, const T &n, const T &k, const T &l) {
    return std::max(std::max(m, n), std::max(k, l));
}

} // FFLAS


#include "fflas-ffpack/utils/fflas_memory.h"
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

//---------------------------------------------------------------------
// specialisations
//---------------------------------------------------------------------

#include "fflas_freduce.h"
#include "fflas_fadd.h"
#include "fflas_fscal.h"
#include "fflas_fassign.h"

#include "fflas_fgemm.inl"
#include "fflas_pfgemm.inl"
// fgemm must be before fgemv according to ScalAndReduce function declaration ?!? PG
#include "fflas_fgemv.inl"
#include "fflas_freivalds.inl"
#include "fflas_fger.inl"
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

//---------------------------------------------------------------------
// Sparse routines
//---------------------------------------------------------------------

#include "fflas_sparse.h"

#endif // __FFLASFFPACK_fflas_H