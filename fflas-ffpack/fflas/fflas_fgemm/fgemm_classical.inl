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

//#include "fflas_bounds_classic.inl"
#include "fflas-ffpack/field/field-general.h"

namespace FFLAS {

	// F is Modular(Balanced)<float/double>
	template<class Field>
	inline void fgemm (const Field & F,
                           const FFLAS_TRANSPOSE ta,
                           const FFLAS_TRANSPOSE tb,
                           const size_t m, const size_t n,const size_t k,
                           const typename Field::Element alpha,
                           typename Field::Element_ptr A, const size_t lda,
                           typename Field::Element_ptr B, const size_t ldb,
                           const typename Field::Element beta,
                           typename Field::Element_ptr C, const size_t ldc,
                           MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::DelayedModularFloatingPointTag> & H)
	{
		
                // Input matrices are unreduced: need to figure out the best option between:
                // - reducing them
                // - making possibly more blocks (smaller kmax)
		
		typename MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::DelayedModularFloatingPointTag>::DelayedField_t::Element alphadf, betadf;
		F.convert (betadf, beta);
		if (F.isMOne (alpha)) {
			alphadf = -1.0;
		} else {
			alphadf = 1.0;
			if (! F.isOne( alpha)) {
				    // Compute y = A*x + beta/alpha.y
				    // and after y *= alpha
				FFLASFFPACK_check(!F.isZero(alpha));
				F.div (betadf, beta, alpha);
			}
		}
			
		if (F.isMOne(betadf)) betadf = -1.0;
		
		size_t kmax = H.MaxDelayedDim (betadf);

		if (kmax <=  k/2 ){
                        // Might as well reduce inputs
                        if (H.Amin < H.FieldMin || H.Amax>H.FieldMax){
				H.initA();
				finit(F, (ta==FflasNoTrans)?m:k, (ta==FflasNoTrans)?k:m, A, lda);
			}
			if (H.Bmin < H.FieldMin || H.Bmax>H.FieldMax){
				H.initB();
				finit(F, (tb==FflasNoTrans)?k:n, (tb==FflasNoTrans)?n:k, B, ldb);
			}
			if (H.Cmin < H.FieldMin || H.Cmax>H.FieldMax){
				H.initC();
				finit(F, m, n, C, ldc);
			}
			kmax = H.MaxDelayedDim (betadf);
		}
		
		if (!kmax){
			MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::GenericTag> HG(H);
			H.initOut();
			return fgemm (F, ta, tb, m,n,k,alpha, A, lda, B, ldb, beta, C, ldc, HG);
		}
		
		size_t k2 = std::min(k,kmax);
		size_t nblock = k / kmax;
		size_t remblock = k % kmax;
		if (!remblock) {
			remblock = kmax;
			--nblock;
		}
                
		size_t shiftA, shiftB;
		if (ta == FflasTrans) shiftA = k2*lda;
		else shiftA = k2;
		if (tb == FflasTrans) shiftB = k2;
		else shiftB = k2*ldb;

		MMHelper<typename associatedDelayedField<Field>::value, 
			 MMHelperAlgo::Classic,
			 typename FieldCategories::FloatingPointTag > Hfp(H);

		fgemm (H.delayedField, ta, tb, m, n, remblock, alphadf, A+nblock*shiftA, lda,
		       B+nblock*shiftB, ldb, betadf, C, ldc, Hfp);

		for (size_t i = 0; i < nblock; ++i) {
			finit(F,m,n,C,ldc);
			Hfp.initC();
			fgemm (H.delayedField, ta, tb, m, n, k2, alphadf, A+i*shiftA, lda,
			       B+i*shiftB, ldb, F.one, C, ldc, Hfp);
		}

                if (!F.isOne(alpha) && !F.isMOne(alpha)){
                    double al; F.convert(al, alpha);
                        if (fabs(al)*std::max(-Hfp.Outmin, Hfp.Outmax)>Hfp.MaxStorableValue){
				finit(F,m,n,C,ldc);
				Hfp.initOut();
			}

			fscalin(H.delayedField, m,n,alpha,C,ldc);

			if (alpha>0){
				H.Outmin = alpha*Hfp.Outmin;
				H.Outmax = alpha*Hfp.Outmax;
			} else {
				H.Outmin = alpha*Hfp.Outmax;
				H.Outmax = alpha*Hfp.Outmin;
			}
		}else {
			H.Outmin = Hfp.Outmin;
			H.Outmax = Hfp.Outmax;
		}
	}
} // FFLAS

namespace FFLAS {

	// Classic multiplication over a generic finite field

	template  < class Field>
	inline void fgemm (const Field& F,
			   const FFLAS_TRANSPOSE ta,
			   const FFLAS_TRANSPOSE tb,
			   const size_t m, const size_t n,const size_t k,
			   const typename Field::Element alpha,
			   typename Field::Element_ptr A, const size_t lda,
			   typename Field::Element_ptr B, const size_t ldb,
			   const typename Field::Element beta,
			   typename Field::Element_ptr C, const size_t ldc,
			   MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::GenericTag> & H)
	{
                // Standard algorithm is performed over the Field, without conversion
                if (F.isZero (beta))
                        fzero (F, m, n, C, ldc);
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
		fscalin(F,m,n,alpha,C,ldc);
        }

	inline void fgemm (const DoubleDomain& ,
			   const FFLAS_TRANSPOSE ta,
			   const FFLAS_TRANSPOSE tb,
			   const size_t m, const size_t n,const size_t k,
			   const DoubleDomain::Element alpha,
			   DoubleDomain::Element_ptr Ad, const size_t lda,
			   DoubleDomain::Element_ptr Bd, const size_t ldb,
			   const DoubleDomain::Element beta,
			   DoubleDomain::Element_ptr Cd, const size_t ldc,
			   MMHelper<DoubleDomain, MMHelperAlgo::Classic, FieldCategories::FloatingPointTag> &H)
	{
		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);
		
                H.setOutBounds(k, alpha, beta);

		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
			     (int)m, (int)n, (int)k, (DoubleDomain::Element) alpha,
			     Ad, (int)lda, Bd, (int)ldb, (DoubleDomain::Element) beta, Cd, (int)ldc);
	}

	inline void fgemm (const FloatDomain& F,
			   const FFLAS_TRANSPOSE ta,
			   const FFLAS_TRANSPOSE tb,
			   const size_t m, const size_t n,const size_t k,
			   const FloatDomain::Element alpha,
			   FloatDomain::Element_ptr Ad, const size_t lda,
			   FloatDomain::Element_ptr Bd, const size_t ldb,
			   const FloatDomain::Element beta,
			   FloatDomain::Element_ptr Cd, const size_t ldc,
			   MMHelper<FloatDomain, MMHelperAlgo::Classic, FieldCategories::FloatingPointTag> & H)
	{
		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);

                H.setOutBounds(k, alpha, beta);
 
		cblas_sgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
			     (int)m, (int)n, (int)k, (FloatDomain::Element) alpha,
			     Ad, (int)lda, Bd, (int)ldb, (FloatDomain::Element) beta,Cd, (int)ldc);
	}
} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_fgemm_classical_INL
