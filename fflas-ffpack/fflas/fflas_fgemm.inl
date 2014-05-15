/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fgemm.inl
 * Copyright (C) 2005 Clement Pernet
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet  < Clement.Pernet@imag.fr >
 *            BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fgemm_INL
#define __FFLASFFPACK_fgemm_INL

namespace FFLAS {

	// Traits and categories will need to be placed in a proper file later
	namespace FieldCategories {
		//! generic ring.
		struct GenericTag{};
		//! If it can init/convert elements to/from floating point types: float, double
		struct FloatingPointConvertibleTag : public  GenericTag{};
		//! If it is a Modular or ModularBalanced templated by float or double
		struct ModularFloatingPointTag : public GenericTag{};
		//! If it is Modular<float> or ModularBalanced<float>
		struct ModularFloatTag : public ModularFloatingPointTag{};
		//! If it is Modular<double> or ModularBalanced<double>
		struct ModularDoubleTag : public ModularFloatingPointTag{};
		//! If it is a multiprecision field
		struct MultiPrecisionTag : public  GenericTag{};
		//! If it is DoubleDomain or a FloatDomain
		struct FloatingPointTag : public GenericTag{};

	}

	/*! FieldTrait
	*/

	template <class Field>
	struct FieldTraits {typedef typename FieldCategories::GenericTag value;};
	template<>
	struct FieldTraits<FFPACK::Modular<double> > {typedef  FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<FFPACK::Modular<float> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<FFPACK::ModularBalanced<double> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<FFPACK::ModularBalanced<float> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<DoubleDomain> {typedef FieldCategories::FloatingPointTag value;};
	template<>
	struct FieldTraits<FloatDomain> {typedef FieldCategories::FloatingPointTag value;};
	template<typename  Element>
	struct FieldTraits<FFPACK::Modular<Element> > {typedef FieldCategories::FloatingPointConvertibleTag value;};


	//	template<> struct FieldTraits<Modular<integer> > {typedef FieldCategories::MultiPrecisionTag value;};


	template <typename FieldT>
	struct ClassicHelper;

	template <typename FieldT>
	struct Winograd2Helper;

} // FFLAS

namespace FFLAS {


	DoubleDomain associatedDomain (const FFPACK::Modular<double> & ){return DoubleDomain();}
	DoubleDomain associatedDomain (const FFPACK::ModularBalanced<double> & ){return DoubleDomain();}
	FloatDomain associatedDomain (const FFPACK::Modular<float> & ){return FloatDomain();}
	FloatDomain associatedDomain (const FFPACK::ModularBalanced<float> & ){return FloatDomain();}
	// This last defintion is useless, but necessary for compilation
	//(some template combinations are compiled but never run  in practice)
	template<class Field>
	const Field& associatedDomain (const Field& F){return F;}
} // FFLAS

namespace FFLAS { namespace Protected {

	template  < typename FloatElement, class Field >
	inline void fgemm_convert (const Field& F,
				   const FFLAS_TRANSPOSE ta,
				   const FFLAS_TRANSPOSE tb,
				   const size_t m, const size_t n,const size_t k,
				   const typename Field::Element alpha,
				   const typename Field::Element * A, const size_t lda,
				   const typename Field::Element * B, const size_t ldb,
				   const typename Field::Element beta,
				   typename Field::Element* C, const size_t ldc,
				   const ClassicHelper<FieldCategories::FloatingPointConvertibleTag> & H
				  );

	template  < typename FloatElement, class Field >
	inline void fgemm_convert (const Field& F,
				   const FFLAS_TRANSPOSE ta,
				   const FFLAS_TRANSPOSE tb,
				   const size_t m, const size_t n,const size_t k,
				   const typename Field::Element alpha,
				   const typename Field::Element * A, const size_t lda,
				   const typename Field::Element * B, const size_t ldb,
				   const typename Field::Element beta,
				   typename Field::Element* C, const size_t ldc,
				   const Winograd2Helper<FieldCategories::FloatingPointConvertibleTag> & H
				  );
} // Protected
} // FFLAS

// #include "fflas_fgemm/matmul_algos.inl"
#include "fflas_fgemm/fgemm_classical.inl"
#include "fflas_fgemm/fgemm_winograd.inl"
// #include "fflas_fgemm/gemm_bini.inl"

// fgemm
namespace FFLAS {

	template<class Field>
	typename Field::Element*
	fgemm( const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       const typename Field::Element* A, const size_t lda,
	       const typename Field::Element* B, const size_t ldb,
	       const typename Field::Element beta,
	       typename Field::Element* C, const size_t ldc
	       , const int w = (int) -1)
	{
		// typedef typename Field::Element Element ;
		if (!k) {
			fscalin(F, m, n, beta, C, ldc);
			return C;
		}

		if (!m || !n) {
			return C;
		}

		if (F.isZero (alpha)){
			fscalin(F, m, n, beta, C, ldc);
			return C;
		}

#ifndef NDEBUG
		/*  check if alpha is invertible.
		 *  XXX do it in F.isInvertible(Element&) ?
		 *  XXX do it in return status of F.inv(Element&,Element&)
		 */
		typename Field::Element e ;
		F.init(e,beta);
		// F.init(e,F.one);
		F.divin(e,alpha);
		F.mulin(e,alpha);
		// FFLASFFPACK_check(F.isOne(e));
		FFLASFFPACK_check(F.areEqual(e,beta));
#endif

#if 0
		// detect fgemv
		if (n == 1 and ...) {}
		// detect fger
		if (k==1 and ...) {}
#endif

		Winograd2Helper<typename FieldTraits<Field>::value> H (w);
		H.computeParameters(F,m,n,k,alpha,beta);
		fgemm2 (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, H);
		return C;
	}

} // FFLAS


// fsquare
namespace FFLAS {
	template < class Field >
	inline typename Field::Element*
	fsquare (const Field& F,
		 const FFLAS_TRANSPOSE ta,
		 const size_t n, const typename Field::Element alpha,
		 const typename Field::Element* A, const size_t lda,
		 const typename Field::Element beta,
		 typename Field::Element* C, const size_t ldc)
	{

		double alphad, betad;
		F.convert (alphad, alpha);
		if (F.isMOne (beta))
			betad = -1.0;
		else
			F.convert (betad, beta);

		//! @bug why double ?
		// Double  matrices initialisation
		DoubleDomain::Element * Ad = new DoubleDomain::Element[n*n];
		DoubleDomain::Element * Cd = new DoubleDomain::Element[n*n];
		// Conversion finite Field = >  double
		fconvert (F, n, n, Ad, n, A, lda);
		if (!F.isZero(beta)) fconvert(F, n, n, Cd, n, C, ldc);

		// Call to the blas Multiplication
		FFLASFFPACK_check(n);
		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE)ta,
			     (CBLAS_TRANSPOSE)ta, (int)n, (int)n, (int)n,
			     (DoubleDomain::Element) alphad, Ad, (int)n, Ad, (int)n,
			     (DoubleDomain::Element) betad, Cd, (int)n);
		// Conversion double = >  Finite Field
		delete[] Ad;
		finit (F,n,n, Cd, n, C, ldc);
		delete[] Cd;
		return C;
	}

	namespace Protected {

		// F is Modular(Balanced)<float/double>
		template < class Field >
		inline typename Field::Element*
		fsquareCommon (const Field& F,
			       const FFLAS_TRANSPOSE ta,
			       const size_t n, const typename Field::Element alpha,
			       const typename Field::Element* A, const size_t lda,
			       const typename Field::Element beta,
			       typename Field::Element* C, const size_t ldc)
		{
			typedef typename Field::Element Element ; // double or float
			if (C==A) {
				Element * Ad = new Element[n*n];
				fcopy(F,n,n,A,lda,Ad,n);
				fgemm (F, ta, ta, n, n, n, alpha, Ad, n, Ad, n, beta, C, ldc);
				delete[] Ad;
			}
			else
				fgemm (F, ta, ta, n, n, n, alpha, A, lda, A, lda, beta, C, ldc);
			// Conversion double/float = >  Finite Field
			finit(F,n,n,C,ldc);
			return C;

		}

	} // Protected

	template <>
	inline double* fsquare (const  FFPACK:: ModularBalanced<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  FFPACK:: ModularBalanced<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline double* fsquare (const  FFPACK:: Modular<double> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const double alpha,
				const double* A, const size_t lda,
				const double beta,
				double* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

	template <>
	inline float * fsquare (const  FFPACK:: Modular<float> & F,
				const FFLAS_TRANSPOSE ta,
				const size_t n, const float alpha,
				const float* A, const size_t lda,
				const float beta,
				float* C, const size_t ldc)
	{
		return Protected::fsquareCommon(F,ta,n,alpha,A,lda,beta,C,ldc);
	}

} // FFLAS

#endif // __FFLASFFPACK_fgemm_INL
