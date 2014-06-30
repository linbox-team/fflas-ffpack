/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas.h
 * Copyright (C) 2005,2013,2014 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Written by BB <bbboyer@ncsu.edu>
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

#include <cmath>
#include <cstring>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/field/unparametric.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-positive.h"

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

//#define LB_TRTR
/// @brief FFLAS: <b>F</b>inite <b>F</b>ield <b>L</b>inear <b>A</b>lgebra <b>S</b>ubroutines.
namespace FFLAS {

	// public:
	/// Is matrix transposed ?
	enum FFLAS_TRANSPOSE
	{
		FflasNoTrans=111, /**< Matrix is not transposed */
		FflasTrans  =112  /**< Matrix is transposed */
	};
	/// Is triangular matrix's shape upper ?
	enum FFLAS_UPLO
	{
		FflasUpper=121,  /**< Triangular matrix is Upper triangular (if \f$i>j\f$ then \f$T_{i,j} = 0\f$)*/
		FflasLower=122   /**< Triangular matrix is Lower triangular (if \f$i<j\f$ then \f$T_{i,j} = 0\f$)*/
	};

	/// Is the triangular matrix implicitly unit diagonal ?
	enum FFLAS_DIAG
	{
		FflasNonUnit=131 ,  /**< Triangular matrix has an explicit general diagonal */
		FflasUnit   =132    /**< Triangular matrix has an implicit unit diagonal (\f$T_{i,i} = 1\f$)*//**< */
	};

	/// On what side ?
	enum FFLAS_SIDE
	{
		FflasLeft  = 141, /**< Operator applied on the left */
		FflasRight = 142  /**< Operator applied on the rigth*/
	};

	/** \p FFLAS_BASE  determines the type of the element representation for Matrix Mult kernel.  */
	enum FFLAS_BASE
	{
		FflasDouble  = 151,  /**<  to use the double precision BLAS */
		FflasFloat   = 152,  /**<  to use the single precison BLAS */
		FflasGeneric = 153   /**< for any other domain, that can not be converted to floating point integers */
	};





	/* Representations of Z with floating point elements*/

	typedef FFPACK::UnparametricField<float> FloatDomain;
	typedef FFPACK::UnparametricField<double> DoubleDomain;

	namespace Protected {

		template <class X,class Y>
		class AreEqual {
		public:
			static const bool value = false;
		};

		template <class X>
		class AreEqual<X,X> {
		public:
			static const bool value = true;
		};
	} // Protected
} // class FFLAS

namespace FFLAS {

	template <class T>
	const T& min3(const T & m, const T & n , const T & k)	{return std::min(m,std::min(n,k));}

	template <class T>
	const T& max3(const T & m, const T & n , const T & k)	{return std::max(m,std::min(n,k));}

	template <class T>
	const T& min4(const T & m, const T & n , const T & k, const T & l)
	{return std::min(std::min(m,n),std::min(k,l));}

	template <class T>
	const T& max4(const T & m, const T & n , const T & k, const T & l)
	{return std::max(std::max(m,n),std::max(k,l));}


} // FFLAS


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
// specialisatons
//---------------------------------------------------------------------

#include "fflas_fadd.inl"
#include "fflas_fscal.inl"
#include "fflas_fcopy.inl"
#include "fflas_finit.inl"

//#ifdef __FFLASFFPACK_USE_OPENMP
#include "fflas_blockcuts.inl"
#include "fflas_pfgemm.inl"
#include "fflas_pftrsm.inl"
//#endif

#include "fflas_fgemv.inl"
#include "fflas_fger.inl"
#include "fflas_ftrsm.inl"
#include "fflas_ftrmm.inl"
#include "fflas_ftrsv.inl"
#include "fflas_faxpy.inl"
#include "fflas_fdot.inl"

//---------------------------------------------------------------------
// Sparse routines
//---------------------------------------------------------------------


//#include "fflas_sparse_fgemv.h"

#if 0
//BB
#ifdef LB_TRTR
#include "fflas_ftrtr.inl"
#endif
#endif


#undef LB_TRTR

#endif // __FFLASFFPACK_fflas_H


