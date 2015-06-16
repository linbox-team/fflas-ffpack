/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas_enum.h
 * Copyright (C) The FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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
#ifndef __FFLASFFPACK_enum_INL
#define __FFLASFFPACK_enum_INL

namespace FFLAS {

	/// Storage by row or col ?
	enum FFLAS_ORDER {
		FflasRowMajor=101, /**< row major */
		FflasColMajor=102  /**< col major */
	};
	// public:
	/// Is matrix transposed ?
	enum FFLAS_TRANSPOSE {
		FflasNoTrans = 111, /**< Matrix is not transposed */
		FflasTrans   = 112  /**< Matrix is transposed */
	};
	/// Is triangular matrix's shape upper ?
	enum FFLAS_UPLO {
		FflasUpper = 121, /**< Triangular matrix is Upper triangular (if \f$i>j\f$ then \f$T_{i,j} = 0\f$)*/
		FflasLower = 122  /**< Triangular matrix is Lower triangular (if \f$i<j\f$ then \f$T_{i,j} = 0\f$)*/
	};

	/// Is the triangular matrix implicitly unit diagonal ?
	enum FFLAS_DIAG {
		FflasNonUnit = 131, /**< Triangular matrix has an explicit arbitrary diagonal */
		FflasUnit    = 132 /**< Triangular matrix has an implicit unit diagonal (\f$T_{i,i} = 1\f$)*/ /**< */
	};

	/// On what side ?
	enum FFLAS_SIDE {
		FflasLeft  = 141,/**< Operator applied on the left */
		FflasRight = 142 /**< Operator applied on the rigth*/
	};

	/** \p FFLAS_BASE  determines the type of the element representation for Matrix Mult kernel. (deprecated, should not be used) */
	enum FFLAS_BASE {
		FflasDouble  = 151,  /**<  to use the double precision BLAS */
		FflasFloat   = 152,  /**<  to use the single precison BLAS */
		FflasGeneric = 153   /**< for any other domain, that can not be converted to floating point integers */
	};
}
#endif // __FFLASFFPACK_enum_INL
