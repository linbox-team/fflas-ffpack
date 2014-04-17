/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_finit.inl
 * Copyright (C) 2014 Pascal Giorgi
 *
 * Written by Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
 * BB<bboyer@ncsu.edu>
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
#ifndef __FFLASFFPACK_fflas_init_INL
#define __FFLASFFPACK_fflas_init_INL

#include "fflas-ffpack/fflas/fflas_avx_functions.h"
#include "fflas-ffpack/field/unparametric.h"

#define FFLASFFPACK_COPY_INIT 100

namespace FFLAS {

	/***************************/
	/*         LEVEL 1         */
	/***************************/



	template<>
	void finit (const FFPACK:: Modular<double> & F, const size_t m,
		    double * A, const size_t incX)
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=1./p;
			vectorised::modp<true,false>(A,A,m,p,invp,0,p-1);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				double * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				double * Ac = new double[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}
		}
	}


	template<>
	void finit (const FFPACK:: Modular<float> & F, const size_t m,
		    float * A, const size_t incX)
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=1.f/p;
			vectorised::modp<true,false>(A,A,m,p,invp,0,p-1);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				float * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				float * Ac = new float[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}

		}

	}


	template<>
	void finit (const FFPACK:: ModularBalanced<double> & F, const size_t m,
		    double * A, const size_t incX)
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=1./p;
			double pmax = (p-1)/2 ;
			double pmin = pmax-p+1;
			vectorised::modp<false,false>(A,A,m,p,invp,pmin,pmax);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				double * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				double * Ac = new double[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}
		}

	}

	template<>
	void finit (const FFPACK:: ModularBalanced<float> & F, const size_t m,
		    float * A, const size_t incX)
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=1.f/p;
			float pmax = (p-1)/2 ;
			float pmin = pmax-p+1;
			vectorised::modp<false,false>(A,A,m,p,invp,pmin,pmax);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_INIT) {
				float * Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.init( *Xi , *Xi);
			}
			else {
				float * Ac = new float[m] ;
				fcopy(F,m,A,incX,Ac,1);
				finit(F,m,Ac,1);
				fcopy(F,m,Ac,1,A,incX);
				delete[] Ac;
			}
		}

	}



	template<>
	void finit (const FFPACK:: Modular<double> & F, const size_t m,
		    const double * B, const size_t incY,
		    double * A, const size_t incX)
	{
		if(incX == 1 && incY == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=1./p;
			vectorised::modp<true,false>(A,B,m,p,invp,0,p-1);
		}
		else {
			double * Xi = A ;
			const double * Yi = B ;
			for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
				F.init( *Xi , *Yi);
		}
	}


	template<>
	void finit (const FFPACK:: Modular<float> & F, const size_t m,
		    const float * B, const size_t incY,
		    float * A, const size_t incX)
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=1.f/p;
			vectorised::modp<true,false>(A,B,m,p,invp,0,p-1);
		}
		else {
			float * Xi = A ;
			const float * Yi = B ;
			for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
				F.init( *Xi , *Yi);
		}

	}


	template<>
	void finit (const FFPACK:: ModularBalanced<double> & F, const size_t m,
		    const double * B, const size_t incY,
		    double * A, const size_t incX)
	{
		if(incX == 1) {
			double p, invp;
			p=(double)F.cardinality();
			invp=1./p;
			double pmax = (p-1)/2 ;
			double pmin = pmax-p+1;
			vectorised::modp<false,false>(A,B,m,p,invp,pmin,pmax);
		}
		else {
			double * Xi = A ;
			const double * Yi = B ;
			for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
				F.init( *Xi , *Yi);
		}
	}

	template<>
	void finit (const FFPACK:: ModularBalanced<float> & F, const size_t m,
		    const float * B, const size_t incY,
		    float * A, const size_t incX)
	{
		if(incX == 1) {
			float p, invp;
			p=(float)F.cardinality();
			invp=1.f/p;
			float pmax = (p-1)/2 ;
			float pmin = pmax-p+1;
			vectorised::modp<false,false>(A,B,m,p,invp,pmin,pmax);
		}
		else {
			float * Xi = A ;
			const float * Yi = B ;
			for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
				F.init( *Xi , *Yi);
		}

	}


	// just to make sure
#if 0
	template<class T>
	void finit (const FFPACK:: UnparametricField<T> & , const size_t ,
		    T * , const size_t)
	{
		return ;
	}
#endif

	template<class Field>
	void
	finit (const Field& F, const size_t n,
	       typename Field::Element * X, const size_t incX)
	{
		typename Field::Element * Xi = X ;

		if (incX == 1 )
			for (; Xi < X + n ; ++Xi)
				F.init( *Xi, *Xi);
		else
			for (; Xi < X+n*incX; Xi+=incX )
				F.init( *Xi, *Xi );
	}

	template<class Field, class OtherElement>
	void
	finit (const Field& F, const size_t n,
	       const OtherElement * Y, const size_t incY,
	       typename Field::Element * X, const size_t incX)
	{
		typename Field::Element * Xi = X ;
		const OtherElement * Yi = Y ;

		if (incX == 1 && incY == 1)
			for (; Yi < Y + n ; ++Xi, ++Yi)
				F.init( *Xi , *Yi);
		else
			for (; Yi < Y+n*incY; Xi+=incX, Yi += incX )
				F.init( *Xi , *Yi);
	}


	/***************************/
	/*         LEVEL 2         */
	/***************************/



	template<class Field>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       typename Field::Element * A, const size_t lda)
	{
		if (n == lda)
			finit(F,n*m,A,1);
		else
			for (size_t i = 0 ; i < m ; ++i)
				finit(F,n,A+i*lda,1);
		return;
	}

	template<class Field, class OtherElement>
	void
	finit (const Field& F, const size_t m , const size_t n,
	       const OtherElement * B, const size_t ldb,
	       typename Field::Element * A, const size_t lda)
	{
		if (n == lda && n == ldb)
			finit(F,n*m,B,1,A,1);
		else
			for (size_t i = 0 ; i < m ; ++i)
				finit(F,n,B+i*ldb,1,A+i*lda,1);
		return;
	}

} // end of namespace FFLAS
#endif // __FFLASFFPACK_fflas_init_INL

