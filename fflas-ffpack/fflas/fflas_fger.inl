/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fger.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __FFLASFFPACK_fger_INL
#define __FFLASFFPACK_fger_INL
namespace FFLAS {

	template<class Field>
	inline void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      const typename Field::Element * x, const size_t incx,
	      const typename Field::Element * y, const size_t incy,
	      typename Field::Element * A, const size_t lda)
	{
		MMHelper<MMHelperAlgo::Classic, typename FieldTraits<Field>::value, Field > H(F,0);
		fger (F, M, N, alpha, const_cast<typename Field::Element*>(x), incx, const_cast<typename Field::Element*>(y), incy, A, lda, H);
		finit (F, M, N, A, lda);
	}
} //FFLAS

namespace FFLAS { namespace Protected {
	template<class FloatElement, class Field>
	inline void
	fger_convert (const Field& F, const size_t M, const size_t N,
		      const typename Field::Element alpha,
		      const typename Field::Element * x, const size_t incx,
		      const typename Field::Element * y, const size_t incy,
		      typename Field::Element * A, const size_t lda)
	{
		FFPACK::Modular<FloatElement> G((FloatElement) F.characteristic());
		FloatElement alphaf;
		F.convert (alphaf, alpha);

		FloatElement * Af = new FloatElement[M*N];
		FloatElement * Xf = new FloatElement[M];
		FloatElement * Yf = new FloatElement[N];

		fconvert(F, M, N, Af, N, A, lda);
		fconvert(F, M, Xf, 1, x, incx);
		fconvert(F, M, Yf, 1, y, incy);

		fger (G, M, N, alphaf, Xf, 1, Yf, 1, Af, N);

		finit(F, M, N, Af, N, A, lda);

		delete[] Af;
		delete[] Xf;
		delete[] Yf;
	}
}// Protected
}// FFLAS

namespace FFLAS{
	template<class Field>
	inline void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      typename Field::Element * x, const size_t incx,
	      typename Field::Element * y, const size_t incy,
	      typename Field::Element * A, const size_t lda,
	      MMHelper<MMHelperAlgo::Classic, FieldCategories::FloatingPointConvertibleTag, Field> & H)
	{
		if (F.characteristic() < DOUBLE_TO_FLOAT_CROSSOVER)
			return Protected::fger_convert<float,Field>(F,M,N,alpha,x, incx, y,incy, A, lda);
		else
			return Protected::fger_convert<double,Field>(F,M,N,alpha,x, incx, y,incy, A, lda);
	}

	template<class Field>
	inline void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      typename Field::Element * x, const size_t incx,
	      typename Field::Element * y, const size_t incy,
	      typename Field::Element * A, const size_t lda,
	      MMHelper<MMHelperAlgo::Classic, FieldCategories::DelayedModularFloatingPointTag, Field> & H)
	{
		typename MMHelper<MMHelperAlgo::Classic, FieldCategories::DelayedModularFloatingPointTag, Field>::DelayedField_t::Element alphadf;
		if (F.isMOne( alpha)) alphadf = -1.0;
		else alphadf = 1.0;

		MMHelper<MMHelperAlgo::Classic,
			 typename FieldCategories::FloatingPointTag, 
			 typename associatedDelayedField<Field>::value > Hfp(H);

		if (Hfp.MaxDelayedDim(1.0) < 1){
			
			if (Hfp.Amin < H.FieldMin || Hfp.Amax>H.FieldMax){
				Hfp.initA();
				finit(F, M, x, incx);
			}
			if (Hfp.Bmin < H.FieldMin || Hfp.Bmax>H.FieldMax){
				Hfp.initB();
				finit(F, N, y, incy);
			}
			if (Hfp.Cmin < H.FieldMin || Hfp.Cmax>H.FieldMax){
				Hfp.initC();
				finit(F, M, N, A, lda);
			}
		}
		Hfp.Outmin = Hfp.FieldMin;
		Hfp.Outmax = Hfp.FieldMax;

		fger (H.delayedField, M, N, alphadf, x, incx, y, incy, A, lda, Hfp);

		if (!F.isOne(alpha) && !F.isMOne(alpha)){
                        if (abs(alpha)*std::max(-Hfp.Outmin, Hfp.Outmax)>Hfp.MaxStorableValue){
				finit (F, M, N, A, lda);
				Hfp.initOut();
			}
			fscalin(H.delayedField, M, N, alpha, A, lda);
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

	template<class Field>
	inline void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      const typename Field::Element * x, const size_t incx,
	      const typename Field::Element * y, const size_t incy,
	      typename Field::Element * A, const size_t lda,
	      MMHelper<MMHelperAlgo::Classic, FieldCategories::GenericTag, Field> & H)
	{

		typename Field::Element tmp;
		const typename Field::Element* xi=x, *yj=y;
		typename Field::Element* Ai=A;

		if (F.isZero(alpha)) return ;


		if ( M < N ){
			if ( F.isOne( alpha ) )
				for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
					yj = y;
					for (size_t j = 0; j < N; ++j, yj+=incy )
						F.axpyin( *(Ai+j), *xi, *yj );
				}
			else if ( F.isMOne( alpha ) )
				for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
					F.neg( tmp, *xi );
					yj = y;
					for (size_t j = 0; j < N; ++j, yj+=incy )
						F.axpyin( *(Ai+j), tmp, *yj );
				}
			else
				for ( ; Ai < A+M*lda; Ai+=lda, xi+=incx ){
					F.mul( tmp, alpha, *xi );
					yj = y;
					for (size_t j = 0; j < N; ++j, yj+=incy )
						F.axpyin( *(Ai+j), tmp, *yj );
				}
		} else {
			if ( F.isOne( alpha ) ){
				for ( ; Ai < A+N; ++Ai, yj+=incy ){
					xi = x;
					for (size_t i = 0; i < M; ++i, xi+=incx )
						F.axpyin( *(Ai+i*lda), *xi, *yj );
				}
			}
			else if ( F.isMOne( alpha ) )
				for ( ; Ai < A+N; ++Ai, yj+=incy ){
					F.neg( tmp, *yj );
					xi = x;
					for (size_t i = 0; i < M; ++i, xi+=incx )
						F.axpyin( *(Ai+i*lda), *xi, tmp );
				}
			else
				for ( ; Ai < A+N; ++Ai, yj+=incy ){
					F.mul( tmp, alpha, *yj );
					xi = x;
					for (size_t i = 0; i < M; ++i, xi+=incx )
						F.axpyin( *(Ai+i*lda), *xi, tmp );
				}
		}

	}

	inline void
	fger( const DoubleDomain& F, const size_t M, const size_t N,
	      const DoubleDomain::Element alpha,
	      const DoubleDomain::Element * x, const size_t incx,
	      const DoubleDomain::Element * y, const size_t incy,
	      DoubleDomain::Element * A, const size_t lda,
	      MMHelper<MMHelperAlgo::Classic, FieldCategories::FloatingPointTag, DoubleDomain> & H)
	{
		if (F.isZero(alpha)) return ;

		FFLASFFPACK_check(lda);
		cblas_dger( CblasRowMajor, (int)M, (int)N, alpha, x, (int)incx, y, (int)incy, A, (int)lda );
		H.setOutBounds (1, alpha, 1.0);
	}

	inline void
	fger( const FloatDomain& F, const size_t M, const size_t N,
	      const FloatDomain::Element alpha,
	      const FloatDomain::Element * x, const size_t incx,
	      const FloatDomain::Element * y, const size_t incy,
	      FloatDomain::Element * A, const size_t lda,
	      MMHelper<MMHelperAlgo::Classic, FieldCategories::FloatingPointTag, FloatDomain> & H)
	{
		if (F.isZero(alpha)) return ;

		FFLASFFPACK_check(lda);
		cblas_sger( CblasRowMajor, (int)M, (int)N, alpha, x, (int)incx, y, (int)incy, A, (int)lda );
		H.setOutBounds (1, alpha, 1.0);
	}

} // FFLAS
#endif // __FFLASFFPACK_fger_INL
