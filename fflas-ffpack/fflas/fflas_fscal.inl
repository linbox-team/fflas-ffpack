/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_faxpy.inl
 * Copyright (C) 2014 FFLAS-FFPACK group
 *
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

#ifndef __FFLASFFPACK_fscal_INL
#define __FFLASFFPACK_fscal_INL

namespace FFLAS { namespace vectorised {

#ifdef __FFLASFFPACK_USE_SIMD

	template<class SimdT, class Element>
	inline typename std::enable_if<is_simd<SimdT>::value, void>::type
	VEC_SCAL(SimdT & C, SimdT & ALPHA, SimdT & Q, SimdT & T, SimdT & P, SimdT & NEGP, SimdT & INVP, SimdT & MIN, SimdT & MAX)
	{
		using simd = Simd<Element>;
		Q = simd::mul(C,INVP);
		C = simd::mul(C,ALPHA);
		Q = simd::floor(Q);
		C = simd::fnmadd(C,Q,P);
		Q = simd::greater(C,MAX);
		T = simd::lesser(C,MIN);
		Q = simd::vand(Q,NEGP);
		T = simd::vand(T,P);
		Q = simd::vor(Q,T);
		C = simd::add(C,Q);
	}

	template<class Element, class T1, class T2>
	inline typename std::enable_if<std::is_floating_point<Element>::value, void>::type
	scalp(Element *T, const Element alpha, const Element * U, const size_t n, const Element p, const Element invp, const T1 min_, const T2 max_)
	{
		Element min = (Element)min_, max=(Element)max_;
		using simd = Simd<Element>;
		using vect_t = typename simd::vect_t;

		size_t i = 0;

		if (n < simd::vect_size)
		{
			for (; i < n ; i++)
			{
				T[i]=monfmod(alpha*U[i], p);
				T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
			return;

		}
		vect_t C,Q,P,NEGP,INVP,TMP,MIN,MAX,ALPHA;
		ALPHA = simd::set1(alpha);
		P   = simd::set1(p);
		NEGP = simd::set1(-p);
		INVP = simd::set1(invp);
		MIN = simd::set1(min);
		MAX = simd::set1(max);
		long st = long(T) % simd::alignment;
		if (st)
		{ // the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)
			for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j+=sizeof(Element), i++)
			{
				T[i] = monfmod(alpha*U[i], p);
				T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
		}
		FFLASFFPACK_check((long(T+i) % simd::alignment == 0));
		if ((long(U+i)%simd::alignment==0))
		{
			// perform the loop using 256 bits SIMD
			for (;i <= n - simd::vect_size ; i += simd::vect_size)
			{
				C = simd::load(U+i);
				VEC_SCAL<vect_t,Element>(C, ALPHA, Q, TMP, P, NEGP, INVP, MIN, MAX);
				simd::store(T+i,C);
			}
		}
		// perform the last elt from T without SIMD
		for (; i < n ; i++)
		{
			T[i] = monfmod(alpha*U[i],p);
			T[i] -= (T[i] > max) ? p : 0;
			T[i] += (T[i] < min) ? p : 0;
		}
	}

#else

	template<class Element, class T1, class T2>
	void
	scalp(Element *T, const Element alpha, const Element * U, const size_t n, const Element p, const Element invp, const T1 min_, const T2 max_)
	{
		Element min = (Element)min_, max=(Element)max_;

		size_t i = 0;

		{
			for (; i < n ; i++)
			{
				T[i]=monfmod(alpha*U[i], p);
				T[i] -= (T[i] > max) ? p : 0;
				T[i] += (T[i] < min) ? p : 0;
			}
			return;

		}

	}

#endif // __FFLASFFPACK_USE_SIMD
} // vectorised
} // FFLAS

namespace FFLAS {

	/***************************/
	/*         LEVEL 1         */
	/***************************/


	template<class Field>
	inline void
	fscal( const Field& F, const size_t N,
	       const typename Field::Element a,
	       typename Field::ConstElement_ptr X, const size_t incX,
	       typename Field::Element_ptr Y, const size_t incY )
	{
		// details::fscal(F,N,a,X,incX,Y,incY, typename FieldTraits<Field>::value() );
		if (F.isOne(a)) {
			fassign(F,N,X,incX,Y,incY);
			return ;
		}

		typename Field::ConstElement_ptr Xi = X;
		typename Field::Element_ptr Yi = Y;
		if (F.areEqual(a,F.mOne)){
			fneg(F,N,X,incX,Y,incY);
			return;
		}

		if (F.isZero(a)){
			fzero(F,N,Y,incY);
			return;
		}

		if (incX == 1 && incY == 1)
			for (size_t i = 0 ; i < N ; ++i)
				F.mul( Y[i], a, X[i] );
		else
		for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
			F.mul( *Yi, a, *Xi );
	}

	template<class Field>
	inline void
	fscalin (const Field& F, const size_t n, const typename Field::Element a,
	       typename Field::Element_ptr X, const size_t incX)
	{
		if (F.isOne(a))
			return ;

		if (F.isMOne(a)){
			fnegin(F,n,X,incX);
			return;
		}

		if (F.isZero(a)){
			fzero(F,n,X,incX);
			return;
		}

		typename Field::Element_ptr Xi = X ;

		if ( incX == 1)
			for (size_t i = 0 ; i < n ; ++i)
				F.mulin( X[i], a);
		else

			for (; Xi < X+n*incX; Xi+=incX )
				F.mulin( *Xi, a);
	}

	template<>
	inline void
	fscal( const Givaro::DoubleDomain& , const size_t N,
	       const Givaro::DoubleDomain::Element a,
	       Givaro::DoubleDomain::ConstElement_ptr x, const size_t incx,
	       Givaro::DoubleDomain::Element_ptr y, const size_t incy )
	{
		cblas_dcopy( (int)N, x, (int)incy, y, (int)incy);
		cblas_dscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscal( const Givaro::FloatDomain& , const size_t N,
	       const Givaro::FloatDomain::Element a,
	       Givaro::FloatDomain::ConstElement_ptr x, const size_t incx,
	       Givaro::FloatDomain::Element_ptr y, const size_t incy )
	{
		cblas_scopy( (int)N, x, (int)incy, y, (int)incy);
		cblas_sscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscalin( const Givaro::DoubleDomain& , const size_t N,
	       const Givaro::DoubleDomain::Element a,
	       Givaro::DoubleDomain::Element_ptr y, const size_t incy )
	{

		cblas_dscal( (int)N, a, y, (int)incy);
	}

	template<>
	inline void
	fscalin( const Givaro::FloatDomain& , const size_t N,
	       const Givaro::FloatDomain::Element a,
	       Givaro::FloatDomain::Element_ptr y, const size_t incy )
	{

		cblas_sscal( (int)N, a, y, (int)incy);
	}


	template<>
	inline void
	fscalin( const Givaro::Modular<float>& F , const size_t N,
	       const float a,
	       float * X, const size_t incX )
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=a/p;
			vectorised::scalp(X,a,X,N,p,invp,0,p-1);
		}
		else {
			float * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}

	template<>
	inline void
	fscalin( const Givaro::ModularBalanced<float>& F , const size_t N,
	       const float a,
	       float * X, const size_t incX )
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=a/p;
			vectorised::scalp(X,a,X,N,p,invp,F.minElement(),F.maxElement());
		}
		else {
			float * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}

	template<>
	inline void
	fscalin( const Givaro::Modular<double>& F , const size_t N,
	       const double a,
	       double * X, const size_t incX )
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=a/p;
			vectorised::scalp(X,a,X,N,p,invp,0,p-1);
		}
		else {
			double * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}

	template<>
	inline void
	fscal( const Givaro::Modular<double>& F , const size_t N,
	       const double a,
		 const double * X, const size_t incX,
		 double * Y, const size_t incY )
	{
		if(incX == 1 && incY==1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=a/p;
			vectorised::scalp(Y,a,X,N,p,invp,0,p-1);
		}
		else {
			const double * Xi = X ;
			double * Yi = Y ;
			for (; Xi < X+N*incX; Xi+=incX,Yi+=incY )
				F.mul(*Yi, *Xi , a);

		}
	}

	template<>
	inline void
	fscalin( const Givaro::ModularBalanced<double>& F , const size_t N,
	       const double a,
	       double * X, const size_t incX )
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=a/p;
			vectorised::scalp(X,a,X,N,p,invp,F.minElement(),F.maxElement());
		}
		else {
			double * Xi = X ;
			for (; Xi < X+N*incX; Xi+=incX )
				F.mulin( *Xi , a);

		}
	}


	/***************************/
	/*         LEVEL 2         */
	/***************************/



	template<class Field>
	void
	fscalin (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element a,
	       typename Field::Element_ptr A, const size_t lda)
	{
		if (F.isOne(a)) {
			return ;
		}
		else if (F.isZero(a)) {
			fzero(F,m,n,A,lda);
		}
		else if (F.isMOne(a)) {
			fnegin(F,m,n,A,lda);
		}
		else {
			if (lda == n) {
				fscalin(F,n*m,a,A,1);
			}
			else {
				for (size_t i = 0 ; i < m ; ++i)
					fscalin(F,n,a,A+i*lda,1);
			}

			return;
		}
	}

	template<class Field>
	void
	fscal (const Field& F, const size_t m , const size_t n,
	       const typename Field::Element a,
	       typename Field::ConstElement_ptr A, const size_t lda,
	       typename Field::Element_ptr B, const size_t ldb)
	{
		if (F.isOne(a)) {
			fassign(F,m,n,A,lda,B,ldb) ;
		}
		else if (F.isZero(a)) {
			fzero(F,m,n,B,ldb);
		}
		else if (F.isMOne(a)) {
			fneg(F,m,n,A,lda,B,ldb);
		}
		else {
			if (n == lda && m == lda)
				fscal(F,m*n,a,A,lda,B,ldb);
			else {
				for (size_t i = 0; i < m ; ++i)
					fscal(F,n,a,A+i*lda,1,B+i*ldb,1);
			}
		}

		return;
	}

} // FFLAS

#endif // __FFLASFFPACK_fscal_INL
