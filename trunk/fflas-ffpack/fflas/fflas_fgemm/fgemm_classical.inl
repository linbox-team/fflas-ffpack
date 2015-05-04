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

#include "fflas-ffpack/field/field-traits.h"

namespace FFLAS {

	// F is a field supporting delayed reductions
	template<class Field>
	inline void fgemm (const Field & F,
                           const FFLAS_TRANSPOSE ta,
                           const FFLAS_TRANSPOSE tb,
                           const size_t m, const size_t n,const size_t k,
                           const typename Field::Element alpha,
                           typename Field::ConstElement_ptr A, const size_t lda,
                           typename Field::ConstElement_ptr B, const size_t ldb,
                           const typename Field::Element beta,
                           typename Field::Element_ptr C, const size_t ldc,
                           MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> & H)
	{		
                // Input matrices are unreduced: need to figure out the best option between:
                // - reducing them
                // - making possibly more blocks (smaller kmax)
		typedef MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::LazyTag> HelperType;
		typename HelperType::DelayedField::Element alphadf, betadf;
		betadf = beta;
		if (F.isMOne (alpha)) {
			alphadf = -H.delayedField.one;
		} else {
			alphadf = F.one;
			if (! F.isOne( alpha)) {
				    // Compute y = A*x + beta/alpha.y
				    // and after y *= alpha
				FFLASFFPACK_check(!F.isZero(alpha));
				typename Field::Element betadalpha;
				F.init(betadalpha);
				F.div (betadalpha, beta, alpha);
				betadf = betadalpha;
			}
		}

		if (F.isMOne(betadf)) betadf = -F.one;

		size_t kmax = H.MaxDelayedDim (betadf);

		if (kmax <=  k/2 ){
                        // Might as well reduce inputs
                        if (H.Amin < H.FieldMin || H.Amax>H.FieldMax){
				H.initA();
				freduce_constoverride (F, (ta==FflasNoTrans)?m:k, (ta==FflasNoTrans)?k:m, A, lda);
			}
			if (H.Bmin < H.FieldMin || H.Bmax>H.FieldMax){
				H.initB();
				freduce_constoverride (F, (tb==FflasNoTrans)?k:n, (tb==FflasNoTrans)?n:k, B, ldb);
			}
			if (H.Cmin < H.FieldMin || H.Cmax>H.FieldMax){
				H.initC();
				freduce (F, m, n, C, ldc);
			}
			kmax = H.MaxDelayedDim (betadf);
		}
		
		if (!kmax){
			MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> HG(H);
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

		typedef MMHelper<typename HelperType::DelayedField, MMHelperAlgo::Classic, ModeCategories::DefaultTag> DelayedHelper_t;
		DelayedHelper_t Hfp(H);
//		write_field(H.delayedField,std::cout<<"classical fgemm"<<std::endl,(DFCElt_ptr)B, k,n,ldb);
		typedef typename HelperType::DelayedField::Element DFElt;
		typedef typename HelperType::DelayedField::Element_ptr DFElt_ptr;
		typedef typename HelperType::DelayedField::ConstElement_ptr DFCElt_ptr;
		
		fgemm (H.delayedField, ta, tb, m, n, remblock, alphadf, 
		       (DFCElt_ptr)A +nblock*shiftA, lda,
		       (DFCElt_ptr)B +nblock*shiftB, ldb, betadf, 
		       (DFElt_ptr)C, ldc, Hfp);

		for (size_t i = 0; i < nblock; ++i) {
			freduce (F, m, n, C, ldc);
			Hfp.initC();
			fgemm (H.delayedField, ta, tb, m, n, k2, alphadf, 
			       (DFCElt_ptr)A +i*shiftA, lda,
			       (DFCElt_ptr)B +i*shiftB, ldb, F.one, 
			       (DFElt_ptr)C, ldc, Hfp);
		}

                if (!F.isOne(alpha) && !F.isMOne(alpha)){
			DFElt al; F.convert(al, alpha);
			if (al<0) al = -al;
			if (al*std::max(-Hfp.Outmin, Hfp.Outmax)>Hfp.MaxStorableValue){
				freduce (F, m, n, C, ldc);
				Hfp.initOut();
			}

			fscalin(H.delayedField, m,n,alpha,(typename DelayedHelper_t::DelayedField_t::Element_ptr)C,ldc);

			if (alpha>0){
				H.Outmin = (const DFElt)(alpha) * Hfp.Outmin;
				H.Outmax = (const DFElt)alpha * Hfp.Outmax;
			} else {
				H.Outmin = (const DFElt)alpha * Hfp.Outmax;
				H.Outmax = (const DFElt)alpha * Hfp.Outmin;
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
			   typename Field::ConstElement_ptr A, const size_t lda,
			   typename Field::ConstElement_ptr B, const size_t ldb,
			   const typename Field::Element beta,
			   typename Field::Element_ptr C, const size_t ldc,
			   MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag> & H)
	{
                // Standard algorithm is performed over the Field, without conversion
                if (F.isZero (beta))
                        fzero (F, m, n, C, ldc);
		else {
			typename Field::Element betadivalpha;
			F.init(betadivalpha);
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
		H.setOutBounds(k,alpha,beta);
		// write_field(F,std::cerr<<"classical \n A = "<<std::endl,A,m,k,lda);
		// write_field(F,std::cerr<<"B = "<<std::endl,B,k,n,ldb);
		// write_field(F,std::cerr<<"C = "<<std::endl,C,m,n,ldc);

        }

	inline void fgemm (const Givaro::DoubleDomain& F,
			   const FFLAS_TRANSPOSE ta,
			   const FFLAS_TRANSPOSE tb,
			   const size_t m, const size_t n,const size_t k,
			   const Givaro::DoubleDomain::Element alpha,
			   Givaro::DoubleDomain::ConstElement_ptr Ad, const size_t lda,
			   Givaro::DoubleDomain::ConstElement_ptr Bd, const size_t ldb,
			   const Givaro::DoubleDomain::Element beta,
			   Givaro::DoubleDomain::Element_ptr Cd, const size_t ldc,
			   MMHelper<Givaro::DoubleDomain, MMHelperAlgo::Classic> &H)
	{
		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);

                H.setOutBounds(k, alpha, beta);

		cblas_dgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
			     (int)m, (int)n, (int)k, (Givaro::DoubleDomain::Element) alpha,
			     Ad, (int)lda, Bd, (int)ldb, (Givaro::DoubleDomain::Element) beta, Cd, (int)ldc);
	}

	inline void fgemm (const Givaro::FloatDomain& F,
			   const FFLAS_TRANSPOSE ta,
			   const FFLAS_TRANSPOSE tb,
			   const size_t m, const size_t n,const size_t k,
			   const Givaro::FloatDomain::Element alpha,
			   Givaro::FloatDomain::ConstElement_ptr Ad, const size_t lda,
			   Givaro::FloatDomain::ConstElement_ptr Bd, const size_t ldb,
			   const Givaro::FloatDomain::Element beta,
			   Givaro::FloatDomain::Element_ptr Cd, const size_t ldc,
			   MMHelper<Givaro::FloatDomain, MMHelperAlgo::Classic> & H)
	{
		FFLASFFPACK_check(lda);
		FFLASFFPACK_check(ldb);
		FFLASFFPACK_check(ldc);
                H.setOutBounds(k, alpha, beta);
		cblas_sgemm (CblasRowMajor, (CBLAS_TRANSPOSE) ta, (CBLAS_TRANSPOSE) tb,
			     (int)m, (int)n, (int)k, (Givaro::FloatDomain::Element) alpha,
			     Ad, (int)lda, Bd, (int)ldb, (Givaro::FloatDomain::Element) beta,Cd, (int)ldc);
	}
} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_fgemm_classical_INL