/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fgemv.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *            Brice Boyer <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_fgemv_INL
#define __FFLASFFPACK_fgemv_INL

namespace FFLAS{ namespace Protected {
	template <typename FloatElement, class Field>
	inline typename Field::Element_ptr
	fgemv_convert (const Field& F,
		       const FFLAS_TRANSPOSE ta,
		       const size_t M, const size_t N,
		       const typename Field::Element alpha,
		       typename Field::ConstElement_ptr A,const size_t lda,
		       typename Field::ConstElement_ptr X,const size_t incX,
		       const typename Field::Element beta,
		       typename Field::Element_ptr Y, const size_t incY)
	{
		FFLASFFPACK_check(lda);

		FFPACK::ModularBalanced<FloatElement> G((FloatElement) F.characteristic());
		FloatElement tmp,alphaf, betaf;
		F.convert (tmp, beta);
		G.init(betaf,tmp);
		F.convert (tmp, alpha);
		G.init(alphaf,tmp);
		size_t ma, na;
		if (ta == FflasTrans) { ma = N; na = M; }
		else { ma = M; na = N; }
		// sizet ldaf = na;
		FloatElement* Af = new FloatElement[M*N];
		FloatElement* Xf = new FloatElement[na];
		FloatElement* Yf = new FloatElement[ma];

		fconvert(F, M, N, Af, N, A, lda);
		finit (G, M, N, Af, N);
		fconvert(F, na, Xf, 1, X, incX);
		finit (G, na, Xf, 1);

		if (!F.isZero(beta)){
			fconvert (F, ma, Yf, 1, Y, incY);
			finit (G, ma, Yf, 1);
		}

		fgemv (G, ta, M, N, alphaf, Af, N, Xf, 1, betaf, Yf, 1);

		finit(F, ma, Yf, 1, Y, incY);
		fflas_delete (Af);
		fflas_delete (Xf);
		fflas_delete (Yf);
		return Y;
	}
	}// Protected
}// FFLAS

namespace FFLAS {
	template<class Field>
	inline  typename Field::Element_ptr
	fgemv (const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       typename Field::Element_ptr A, const size_t lda,
	       typename Field::Element_ptr X, const size_t incX,
	       const typename Field::Element beta,
	       typename Field::Element_ptr Y, const size_t incY,
	       MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::FloatingPointConvertibleTag> & H)
	{
		if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
			return Protected::fgemv_convert<float,Field>(F,ta,M,N,alpha,A,lda,X, incX, beta,Y,incY);
		else
			return Protected::fgemv_convert<double,Field>(F,ta,M,N,alpha,A,lda,X, incX, beta,Y,incY);
	}
}// FFLAS

namespace FFLAS {

	//---------------------------------------------------------------------
	// fgemv: GEneral Matrix Vector Multiplication
	// Computes  Y <- alpha.op(A).X + beta.Y
	// A is M*N,
	//---------------------------------------------------------------------

	template<class Field>
	inline typename Field::Element_ptr
	fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr X, const size_t incX,
	       const typename Field::Element beta,
	       typename Field::Element_ptr Y, const size_t incY)
	{
		MMHelper<Field, MMHelperAlgo::Classic > HW (F, 0);
		return 	fgemv (F, ta, M, N, alpha,
			       FFPACK::fflas_const_cast<typename Field::Element_ptr>(A), lda,
			       FFPACK::fflas_const_cast<typename Field::Element_ptr>(X), incX,
			       beta, Y, incY, HW);
	}

	template<class Field>
	inline typename Field::Element_ptr
	fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr X, const size_t incX,
	       const typename Field::Element beta,
	       typename Field::Element_ptr Y, const size_t incY,
	       MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::ModularFloatingPointTag> & H)
	{

		if (!M) {return Y;}
		size_t Ydim = (ta == FflasNoTrans)?M:N;
		if (!N || F.isZero (alpha)){
			fscalin(F, Ydim, beta, Y, incY);
			return Y;
		}

		typename Field::Element alpha_,beta_;
		F.assign (alpha_,alpha);
		F.assign (beta_,beta);
		if (Protected::AreEqual<Field, FFPACK::Modular<double> >::value ||
		    Protected::AreEqual<Field, FFPACK::ModularBalanced<double> >::value){
			    // Modular<double> need to switch to float if p too small
			if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
				return Protected::fgemv_convert<float,Field>(F,ta,M,N,alpha,A,lda,X,incX,beta,Y,incY);
		}
		if ( !F.isOne(alpha) && !F.isMOne(alpha)){
			F.assign (alpha_, F.one);
			F.div (beta_, beta, alpha);
		}
		MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::DelayedModularFloatingPointTag > HD(F,0);

		fgemv (F, ta, M, N, alpha_,
		       FFPACK::fflas_const_cast<typename Field::Element_ptr>(A), lda,
		       FFPACK::fflas_const_cast<typename Field::Element_ptr>(X), incX,
		       beta_, Y, incY, HD);

		Protected::ScalAndInit (F, Ydim, alpha, Y, incY, HD);
		H.initOut();

		return Y;
	}



}

namespace FFLAS{
	template<class Field>
	inline typename Field::Element_ptr
	fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::ConstElement_ptr X, const size_t incX,
	       const typename Field::Element beta,
	       typename Field::Element_ptr Y, const size_t incY,
	       MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::GenericTag> & H)

	{
		size_t Ydim = (ta==FflasNoTrans)?M:N;

		if (F.isZero (beta))
			fzero (F, Ydim, Y, incY);
		else {
			typename Field::Element betadivalpha;
			FFLASFFPACK_check(!F.isZero(alpha));
			F.div (betadivalpha, beta, alpha);
			fscalin (F, Ydim, betadivalpha, Y, incY);
		}
		if (ta == FflasNoTrans)
			for (size_t i = 0; i < M; ++i)
				Y[i*incY] = fdot(F, N, A+i*lda, 1, X, incX);
		else
			for (size_t i = 0; i < M; ++i)
				Y[i*incY] = fdot(F, M, A+i, lda, X, incX);
		fscalin (F, Ydim, alpha, Y, incY);
		return Y;
	}
}

namespace FFLAS{
	template<class Field>
	inline typename Field::Element_ptr
	fgemv (const Field& F, const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       typename Field::Element_ptr A, const size_t lda,
	       typename Field::Element_ptr X, const size_t incX,
	       const typename Field::Element beta,
	       typename Field::Element_ptr Y, const size_t incY,
	       MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::DelayedModularFloatingPointTag> & H)
	{
		typename MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::DelayedModularFloatingPointTag>::DelayedField_t::Element alphadf=alpha, betadf=beta;

		 size_t Ydim = (ta==FflasNoTrans)?M:N;
		 size_t Xdim = (ta==FflasNoTrans)?N:M;

		if (F.isMOne (alpha)) alphadf = -1.0;
		else {
			alphadf = 1.0;
			if (! F.isOne( alpha)) {
				// Compute y = A*x + beta/alpha.y, then y *= alpha
				FFLASFFPACK_check(!F.isZero(alpha));
				F.div (betadf, beta, alpha);
			}
		}
		if (F.isMOne(betadf)) betadf = -1.0;

		size_t kmax = H.MaxDelayedDim (betadf);

		if (kmax <=  Xdim/2 ){
                        // Might as well reduce inputs
                        if (H.Amin < H.FieldMin || H.Amax>H.FieldMax){
				H.initA();
				finit(F, M, N, A, lda);
			}
			if (H.Bmin < H.FieldMin || H.Bmax>H.FieldMax){
				H.initB();
				finit(F, Xdim, X, incX);
			}
			if (H.Cmin < H.FieldMin || H.Cmax>H.FieldMax){
				H.initC();
				finit(F, Ydim, Y, incY);
			}
			kmax = H.MaxDelayedDim (betadf);
		}

		if (!kmax){
			MMHelper<Field, MMHelperAlgo::Classic, FieldCategories::GenericTag> HG(H);
			H.initOut();
			return fgemv (F, ta, M, N, alpha, A, lda, X, incX, beta, Y, incY, HG);
		}
		size_t k2 = std::min (Xdim, kmax);
		size_t nblock = Xdim / kmax;
		size_t remblock = Xdim % kmax;
		if (!remblock) {
			remblock = kmax;
			--nblock;
		}
		size_t shiftA, M1, N1, Mi, Ni;
		if (ta == FflasTrans) {
			shiftA = k2*lda;
			M1 = remblock;
			Mi = k2;
			Ni = N1 = N;
		}else {
			shiftA = k2;
			Mi = M1 = M;
			N1 = remblock;
			Ni = k2;
		}
		MMHelper<typename associatedDelayedField<Field>::value,
			 MMHelperAlgo::Classic,
			 typename FieldCategories::FloatingPointTag > Hfp(H);

		fgemv (H.delayedField, ta, M1, N1, alphadf, A+nblock*shiftA, lda,
		       X+nblock*k2*incX, incX, betadf, Y, incY, Hfp);

		for (size_t i = 0; i < nblock; ++i) {
			finit(F, Ydim ,Y, incY);
			Hfp.initC();
			fgemv (H.delayedField, ta, Mi, Ni, alphadf, A+i*shiftA, lda,
			       X+i*k2*incX, incX, F.one, Y, incY, Hfp);
		}

                if (!F.isOne(alpha) && !F.isMOne(alpha)){
                        if (abs(alpha)*std::max(-Hfp.Outmin, Hfp.Outmax)>Hfp.MaxStorableValue){
				finit (F, Ydim, Y, incY);
				Hfp.initOut();
			}
			fscalin(H.delayedField, Ydim, alpha, Y, incY);
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
		return Y;
	}
}

namespace FFLAS{
	inline DoubleDomain::Element_ptr
	fgemv (const DoubleDomain& F, const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const DoubleDomain::Element alpha,
	       const DoubleDomain::Element_ptr A, const size_t lda,
	       const DoubleDomain::Element_ptr X, const size_t incX,
	       const DoubleDomain::Element beta,
	       DoubleDomain::Element_ptr Y, const size_t incY,
	       MMHelper<DoubleDomain, MMHelperAlgo::Classic, FieldCategories::FloatingPointTag> & H)
	{
		FFLASFFPACK_check(lda);
		// FFLASFFPACK_check(ldb);
		// FFLASFFPACK_check(ldc);

                H.setOutBounds((ta ==FflasNoTrans)?N:M, alpha, beta);

		cblas_dgemv (CblasRowMajor, (CBLAS_TRANSPOSE) ta,
			     (int)M, (int)N, (DoubleDomain::Element) alpha,
			     A, (int)lda, X, (int)incX, (DoubleDomain::Element) beta, Y, (int)incY);
		return Y;
	}

	inline FloatDomain::Element_ptr
	fgemv (const FloatDomain& F, const FFLAS_TRANSPOSE ta,
	       const size_t M, const size_t N,
	       const FloatDomain::Element alpha,
	       const FloatDomain::Element_ptr A, const size_t lda,
	       const FloatDomain::Element_ptr X, const size_t incX,
	       const FloatDomain::Element beta,
	       FloatDomain::Element_ptr Y, const size_t incY,
	       MMHelper<FloatDomain, MMHelperAlgo::Classic, FieldCategories::FloatingPointTag> & H)
	{
		FFLASFFPACK_check(lda);
		// FFLASFFPACK_check(ldb);
		// FFLASFFPACK_check(ldc);

		H.setOutBounds((ta ==FflasNoTrans)?N:M, alpha, beta);

		cblas_sgemv (CblasRowMajor, (CBLAS_TRANSPOSE) ta,
			     (int)M, (int)N, (FloatDomain::Element) alpha,
			     A, (int)lda, X, (int)incX, (FloatDomain::Element) beta, Y, (int)incY);
		return Y;
	}

}

#endif //  __FFLASFFPACK_fgemv_INL
