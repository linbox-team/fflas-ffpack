/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas-ffpack/ffpack/ffpack_minpoly_construct.inl
 * Copyright (C) 2003 Clement Pernet
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

#ifndef __FFLASFFPACK_ffpack_minpoly_construct_INL
#define __FFLASFFPACK_ffpack_minpoly_construct_INL

#error "not included anywhere"

namespace FFPACK {
	//#define LB_DEBUG
	//---------------------------------------------------------------------
	// MinPoly: Compute the minimal polynomial of (A,v) using an LUP
	// factorization of the Krylov Base (v, Av, .., A^kv)
	// U must be (n+1)*n
	//---------------------------------------------------------------------
	template <class Field, class Polynomial>
	Polynomial&
	MinPoly( const Field& F, Polynomial& minP, const size_t N,
		 typename Field::ConstElement_ptr A, const size_t lda,
		 typename Field::Element_ptr U, size_t ldu, typename Field::Element_ptr X, size_t ldx,
		 size_t* P)
	{

		// nRow is the number of row in the krylov base already computed
		size_t j, k, nRow = 2;
		typename Field::Element_ptr B = FFLAS::fflas_new (F, N, N);
		typename Polynomial::iterator it;
		typename Field::Element_ptr Xi, *Ui;
		typename Field::RandIter g (F);
		bool KeepOn=true;
		// Creating the Krylov Base copy matrix X where to factorize
		//typename Field::Element_ptr X = new elt[(N+1)*N];
#ifdef LB_DEBUG
		for (j=0;j<(N+1)*N;j++)
			X[j] = zero;
#endif
		// Creating the copy of A, where to compute A^2^i
		// Try memcopy here
		for (size_t i=0; i<N; ++i)
			for ( size_t j=0; j<N; ++j)
				*(B+i*N+j) = *(A+i*N+j);

		// Picking a non zero vector
		do{
			for (Ui=U, Xi = X; Ui<U+N; ++Ui, ++Xi){
				g.random (*Ui);
				*Xi = *Ui;
				if (!F.isZero(*Ui))
					KeepOn = false;
			}
		}while(KeepOn);


		nRow = 1;
		size_t nUsedRow = 0;
		// LUP factorization of the Krylov Base Matrix

		k = Protected::LUdivine_construct(F, FflasUnit, N+1, N, B, N, U, ldu, X, N, P,
				       &nRow, N+1, &nUsedRow );
		fflas_delete (B);
		minP.resize(k+1);
		minP[k] = one;
		if (k==1 && F.isZero(*(X+N))){ // minpoly is X
			return minP;
		}
		// m contains the k first coefs of the minpoly
		typename Field::Element_ptr m= FFLAS::fflas_new (F,k,1);
		fcopy( F, k, X+k*N, 1, m, 1);
		ftrsv( F, FflasLower, FflasTrans, FflasNonUnit, k, X, N, m, 1);
		//delete[] X;
		it = minP.begin();
		for (j=0; j<k; ++j, it++){
			F.neg(*it, m[j]);
		}
		fflas_delete (m);
		return minP;
	}
} // FFPACK

#undef LB_DEBUG

#endif //__FFLASFFPACK_ffpack_minpoly_construct_INL
