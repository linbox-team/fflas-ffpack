/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2008, 2014 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

/** @file fflas_fgemm/fgemm_classical.inl
 * @brief Classical \f$2n^3\$f matrix multiplication.
 * @warning The domain is supposed to be a field since some divisions are required for efficiency purposes
 * An alternative has to be written for finite rings if necessary
 */

#ifndef __FFLASFFPACK_fflas_fflas_fgemm_classical_INL
#define __FFLASFFPACK_fflas_fflas_fgemm_classical_INL

#include "fflas_bounds_classic.inl"

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
		//! If it a multiprecision field
		struct MultiPrecisionTag : public  GenericTag{};
		//! If it DoubleDomain or a FloatDomain
		struct FloatingPointTag : public GenericTag{};

	}

	/*! FieldTrait
	 */

	template <class Field>
	struct FieldTraits { typedef typename FieldCategories::GenericTag value;};

	template<> struct FieldTraits<FFPACK::Modular<double> > {typedef  FieldCategories::ModularFloatingPointTag value;};
	template<> struct FieldTraits<FFPACK::Modular<float> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<> struct FieldTraits<FFPACK::ModularBalanced<double> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<> struct FieldTraits<FFPACK::ModularBalanced<float> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<> struct FieldTraits<DoubleDomain> {typedef FieldCategories::FloatingPointTag value;};
	template<> struct FieldTraits<FloatDomain> {typedef FieldCategories::FloatingPointTag value;};
	template<typename  Element> struct FieldTraits<FFPACK::Modular<Element> > {typedef FieldCategories::FloatingPointConvertibleTag value;};


//	template<> struct FieldTraits<Modular<integer> > {typedef FieldCategories::MultiPrecisionTag value;};

	template <typename FieldT>
	struct ClassicHelper : public MMParameters {
		size_t     kmax ;
		FFLAS_BASE base ;
		// ijk order
		ClassicHelper() :
			kmax(0), base(FflasDouble)
		{
		}
		ClassicHelper(size_t k, FFLAS_BASE b):
			kmax(k),base(b)
		{
		}
		template<class Field>
		void computeParameters(const Field & F
				       , const size_t & m
				       , const size_t & n
				       , const size_t & k
				       , const typename Field::Element &alpha
				       , const typename Field::Element &beta)
		{
			typename Field::Element gamma;
			F.div(gamma,beta,alpha);
			Protected::MatMulParametersClassic (F, m, n, k, gamma, kmax, base);
		}

	};

} // FFLAS

namespace FFLAS { namespace BLAS3 {

	// G is (Float/Double)Domain
	template  < class Field, class FloatField >
	inline void ClassicMatmulFloat (const Field& F,
					const FFLAS_TRANSPOSE ta,
					const FFLAS_TRANSPOSE tb,
					const size_t m, const size_t n,const size_t k,
					const typename Field::Element alpha,
					const typename Field::Element * A, const size_t lda,
					const typename Field::Element * B, const size_t ldb,
					const typename Field::Element beta,
					typename Field::Element* C, const size_t ldc,
					// const size_t kmax, const FFLAS_BASE base
					const FloatField & G, const size_t k2,
					ClassicHelper<FieldCategories::FloatingPointConvertibleTag> & H
				       )
	{
		typedef typename FloatField::Element FloatElement ;
		FloatElement alphad, betad;
		FloatElement * Add = new FloatElement[m*k2];
		FloatElement * Bdd = new FloatElement[k2*n];
		FloatElement * Cd = new FloatElement[m*n];

		size_t nblock = k / H.kmax;
		size_t remblock = k % H.kmax;
		if (!remblock) {
			remblock = H.kmax;
			--nblock;
		}
		if (F.isMOne( beta)) betad = -1.0;
		else F.convert (betad, beta);

		if (F.isMOne( alpha)) alphad = -1.0;
		else {
			alphad = 1.0;
			if (! F.isOne( alpha)) {
				// Compute y = A*x + beta/alpha.y
				// and after y *= alpha
				typename Field::Element tmp;
				FFLASFFPACK_check(!F.isZero(alpha));
				F.div (tmp, beta, alpha);
				F.convert (betad, tmp);
			}
		}

		size_t dlda, dldb;
		if (!F.isZero(beta))
			fconvert (F, m, n, Cd, n, C, ldc);

		if (ta == FflasTrans) {
			dlda = m;
			fconvert(F, remblock, m, Add, dlda, A+k2*nblock*lda, lda);
		}
		else {
			dlda = k2;
			fconvert(F, m, remblock, Add, dlda, A+k2*nblock, lda);
		}
		if (tb == FflasTrans) {
			dldb = k2;
			fconvert (F, n, remblock, Bdd, k2, B+k2*nblock, ldb);
		}
		else {
			dldb = n;
			fconvert(F, remblock, n, Bdd, dldb, B+k2*nblock*ldb, ldb);
		}

		fgemm2 (G, ta, tb, m, n, remblock, alphad, Add, dlda,
			Bdd, dldb, betad, Cd, n, ClassicHelper<FieldCategories::FloatingPointTag>() );

		finit (F, m, n, Cd, n, C, ldc);
		fconvert(F, m, n, Cd, n, C, ldc);

		for (size_t i = 0; i < nblock; ++i) {
			if (ta == FflasTrans) fconvert(F, k2, m, Add, dlda, A+k2*i*lda, lda);
			else fconvert(F, m, k2, Add, dlda,  A+k2*i, lda);
			if (tb == FflasTrans) fconvert(F, n, k2, Bdd, dldb, B+k2*i, ldb);
			else fconvert(F, k2, n, Bdd, dldb, B+k2*i*ldb, ldb);

			fgemm2 (G, ta, tb, m, n, k2, alphad, Add, dlda,
				Bdd, dldb, 1.0, Cd, n, ClassicHelper<FieldCategories::FloatingPointTag>());
			finit(F, m, n, Cd, n, C, ldc);
			fconvert(F, m, n, Cd, n, C, ldc);
		}
		if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
			fscalin(F,m,n,alpha,C,ldc);
		}
		delete[] Add;
		delete[] Bdd;
		delete[] Cd;

	}

	// F is Modular(Balanced)<float/double>
	template<class Field>
	inline void ClassicMatmulCommon (const Field & F,
					 const FFLAS_TRANSPOSE ta,
					 const FFLAS_TRANSPOSE tb,
					 const size_t m, const size_t n,const size_t k,
					 const typename Field::Element alpha,
					 const typename Field::Element * A, const size_t lda,
					 const typename Field::Element * B, const size_t ldb,
					 const typename Field::Element beta,
					 typename Field::Element* C, const size_t ldc,
					 // const size_t kmax, const FFLAS_BASE base
					 ClassicHelper<FieldCategories::ModularFloatingPointTag> & H
					)
	{
		typename Field::Element _alpha, _beta;
		// To ensure the initial computation with beta
		size_t k2 = std::min(k,H.kmax);
		size_t nblock = k / H.kmax;
		size_t remblock = k % H.kmax;
		if (!remblock) {
			remblock = H.kmax;
			--nblock;
		}
		if (F.isMOne( beta)) _beta = -1.0;
		else _beta = beta;
		if (F.isMOne( alpha)) _alpha = -1.0;
		else{
			_alpha = 1.0;
			if (! F.isOne( alpha)) {
				// Compute y = A*x + beta/alpha.y
				// and after y *= alpha
				FFLASFFPACK_check(!F.isZero(alpha));
				F.divin (_beta, alpha);
			}
		}
		size_t shiftA, shiftB;
		if (ta == FflasTrans) shiftA = k2*lda;
		else shiftA = k2;
		if (tb == FflasTrans) shiftB = k2;
		else shiftB = k2*ldb;
		
		ClassicHelper<typename FieldCategories::FloatingPointTag> associatedH;
		fgemm2 (associatedDomain(F), ta, tb, m, n, remblock, _alpha, A+nblock*shiftA, lda,
			B+nblock*shiftB, ldb, _beta, C, ldc, associatedH);
		finit(F,m,n,C,ldc);
		for (size_t i = 0; i < nblock; ++i) {
			fgemm2 (associatedDomain(F), ta, tb, m, n, k2, _alpha, A+i*shiftA, lda,
				B+i*shiftB, ldb, F.one, C, ldc, associatedH);
			finit(F,m,n,C,ldc);
		}
		if ((!F.isOne( alpha)) && (!F.isMOne( alpha))) {
			fscalin(F,m,n,alpha,C,ldc);
		}
	}
} // BLAS3
} // FFLAS

namespace FFLAS {

	// Classic Multiplication over double
	// Classic multiplication over a finite field
	template  < class Field >
	inline void fgemm2 (const Field& F,
			    const FFLAS_TRANSPOSE ta,
			    const FFLAS_TRANSPOSE tb,
			    const size_t m, const size_t n,const size_t k,
			    const typename Field::Element alpha,
			    const typename Field::Element * A, const size_t lda,
			    const typename Field::Element * B, const size_t ldb,
			    const typename Field::Element beta,
			    typename Field::Element* C, const size_t ldc,
			    // const size_t kmax, const FFLAS_BASE base
			    ClassicHelper<typename FieldTraits<Field>::categoryTag> & H
			   )
	{

		size_t k2 = std::min(k,H.kmax); // Size of the blocks

		if (k2 > 1) {
			if (H.base == FflasFloat) {
				FloatDomain G ;
				BLAS3::ClassicMatmulFloat(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H,G,k2);
			}
			else { // base == FflasDouble
				DoubleDomain G ;
				BLAS3::ClassicMatmulFloat(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H,G,k2);
			}
		}
		else { // k2 == 1
			FFLASFFPACK_check(k2 == 1);
			// Standard algorithm is performed over the Field, without conversion
			if (F.isZero (beta))
				for (size_t i = 0; i < m; ++i)
					for (size_t j = 0; j < n; ++j)
						F.assign (*(C+i*ldc+j), F.zero);
			else {
				typename Field::Element betadivalpha;
				FFLASFFPACK_check(!F.isZero(alpha));
				F.div (betadivalpha, beta, alpha);
				fscalin(F,m,n,betadivalpha,C,ldc);
			}
			if (ta == FflasNoTrans)
				if (tb == FflasNoTrans)
					for (size_t i = 0; i < m; ++i)
						for (size_t l = 0; l < k; ++l)
							for (size_t j = 0; j < n; ++j)
								F.axpyin (*(C+i*ldc+j), *(A+i*lda+l), *(B+l*ldb+j));
				else
					for (size_t i = 0; i < m; ++i)
						for (size_t j = 0; j < n; ++j)
							for (size_t l = 0; l < k; ++l)
								F.axpyin (*(C+i*ldc+j), *(A+i*lda+l), *(B+j*ldb+l));
			else
				if (tb == FflasNoTrans)
					for (size_t i = 0; i < m; ++i)
						for (size_t l = 0; l < k; ++l)
							for (size_t j = 0; j < n; ++j)
								F.axpyin (*(C+i*ldc+j), *(A+l*lda+i), *(B+l*ldb+j));
				else
					for (size_t i = 0; i < m; ++i)
						for (size_t j = 0; j < n; ++j)
							for (size_t l = 0; l < k; ++l)
								F.axpyin (*(C+i*ldc+j), *(A+l*lda+i), *(B+j*ldb+l));
			if (! F.isOne(alpha))
				fscalin(F,m,n,alpha,C,ldc);
		}
	}

	template<>
	inline void fgemm2(const DoubleDomain& ,
			   const FFLAS_TRANSPOSE ta,
			   const FFLAS_TRANSPOSE tb,
			   const size_t m, const size_t n,const size_t k,
			   const DoubleDomain::Element alpha,
			   const DoubleDomain::Element * Ad, const size_t lda,
			   const DoubleDomain::Element * Bd, const size_t ldb,
			   const DoubleDomain::Element beta,
			   DoubleDomain::Element * Cd, const size_t ldc,
			   // const size_t kmax, const FFLAS_BASE base
			   ClassicHelper<FieldCategories::FloatingPointTag> &
			  )
	{

		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);
		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
			     (int)m, (int)n, (int)k, (DoubleDomain::Element) alpha,
			     Ad, (int)lda, Bd, (int)ldb, (DoubleDomain::Element) beta,Cd, (int)ldc);
	}

	template  <>
	inline void fgemm2 (const FloatDomain& F,
			    const FFLAS_TRANSPOSE ta,
			    const FFLAS_TRANSPOSE tb,
			    const size_t m, const size_t n,const size_t k,
			    const FloatDomain::Element alpha,
			    const FloatDomain::Element * Ad, const size_t lda,
			    const FloatDomain::Element * Bd, const size_t ldb,
			    const FloatDomain::Element beta,
			    FloatDomain::Element * Cd, const size_t ldc,
			    // const size_t kmax, const FFLAS_BASE base
			    ClassicHelper<FieldCategories::FloatingPointTag> &
			   )
	{
		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);

		cblas_sgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
			     (int)m, (int)n, (int)k, (FloatDomain::Element) alpha,
			     Ad, (int)lda, Bd, (int)ldb, (FloatDomain::Element) beta,Cd, (int)ldc);
	}

	template <>
	inline void fgemm2 (const FFPACK:: ModularBalanced<double> & F,
			    const FFLAS_TRANSPOSE ta,
			    const FFLAS_TRANSPOSE tb,
			    const size_t m, const size_t n,const size_t k,
			    const double alpha,
			    const double * A, const size_t lda,
			    const double * B, const size_t ldb,
			    const double beta,
			    double* C, const size_t ldc,
			    // const size_t kmax, const FFLAS_BASE base
			    ClassicHelper<FieldCategories::ModularFloatingPointTag> & H
			   )

	{
		return BLAS3::ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
	}


	template <>
	inline void fgemm2 (const  FFPACK:: ModularBalanced<float> & F,
			    const FFLAS_TRANSPOSE ta,
			    const FFLAS_TRANSPOSE tb,
			    const size_t m, const size_t n,const size_t k,
			    const float alpha,
			    const float * A, const size_t lda,
			    const float * B, const size_t ldb,
			    const float beta,
			    float* C, const size_t ldc,
			    // const size_t kmax, const FFLAS_BASE base
			    ClassicHelper<FieldCategories::ModularFloatingPointTag> & H
			   )
	{
		return BLAS3::ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
	}


	template <>
	inline void fgemm2 (const  FFPACK:: Modular<double> & F,
			    const FFLAS_TRANSPOSE ta,
			    const FFLAS_TRANSPOSE tb,
			    const size_t m, const size_t n,const size_t k,
			    const double alpha,
			    const double * A, const size_t lda,
			    const double * B, const size_t ldb,
			    const double beta,
			    double* C, const size_t ldc,
			    // const size_t kmax, const FFLAS_BASE base
			    ClassicHelper<FieldCategories::ModularFloatingPointTag> & H
			   )
	{
		return BLAS3::ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
	}

	template <>
	inline void fgemm2 (const  FFPACK:: Modular<float> & F,
			    const FFLAS_TRANSPOSE ta,
			    const FFLAS_TRANSPOSE tb,
			    const size_t m, const size_t n,const size_t k,
			    const float alpha,
			    const float * A, const size_t lda,
			    const float * B, const size_t ldb,
			    const float beta,
			    float* C, const size_t ldc,
			    // const size_t kmax, const FFLAS_BASE base
			    ClassicHelper<FieldCategories::ModularFloatingPointTag> & H
			   )
	{
		return BLAS3::ClassicMatmulCommon(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H);
	}

} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_fgemm_classical_INL
