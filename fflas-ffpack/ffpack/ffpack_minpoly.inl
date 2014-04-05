/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* ffpack/ffpack_minpoly.inl
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
#ifndef __FFLASFFPACK_ffpack_minpoly_INL
#define __FFLASFFPACK_ffpack_minpoly_INL
namespace FFPACK {

	template <class Field, class Polynomial>
	Polynomial&
	MinPoly( const Field& F, Polynomial& minP, const size_t N
		 ,const typename Field::Element *A, const size_t lda
		 ,typename Field::Element* X, const size_t ldx
		 ,size_t* P
		 ,const FFPACK_MINPOLY_TAG MinTag// = FfpackDense
		 ,const size_t kg_mc// =0
		 ,const size_t kg_mb//=0
		 ,const size_t kg_j //=0
		 )
	{

		typedef typename Field::Element elt;
		// nRow is the number of row in the krylov base already computed
		size_t j, k ;
		//size_t	nRow = 2;
		typename Polynomial::iterator it;
		elt* Xi, *Ui;
		typename Field::RandIter g (F);
		bool KeepOn=true;
		elt* U = new elt[N];
		// Picking a non zero vector
		do{
			for (Ui=U, Xi = X; Ui<U+N; ++Ui, ++Xi){
				g.random (*Ui);
				*Xi = *Ui;
				if (!F.isZero(*Ui))
					KeepOn = false;
			}
		}while(KeepOn);

		// nRow = 1;
		// LUP factorization of the Krylov Base Matrix
		k = Protected::LUdivine_construct (F, FFLAS::FflasUnit, N+1, N, A, lda, X, ldx, U, P, true,
					MinTag, kg_mc, kg_mb, kg_j);
		//delete[] U;
		minP.resize(k+1);
		minP[k] = F.one;
		if ( (k==1) && F.isZero(*(X+ldx))){ // minpoly is X
			delete[] U;
			for (size_t i=0; i<k; ++i)
				minP[i] = F.zero;
			return minP;
		}
		// U contains the k first coefs of the minpoly
		//elt* m= new elt[k];
		FFLAS::fcopy( F, k, X+k*ldx, 1, U, 1);
		ftrsv( F, FFLAS::FflasLower, FFLAS::FflasTrans, FFLAS::FflasNonUnit, k, X, ldx, U, 1);
		it = minP.begin();
		for (j=0; j<k; ++j, it++){
			F.neg(*it, U[j]);
		}
		delete[] U;
		return minP;
	}

} // FFPACK
#endif // __FFLASFFPACK_ffpack_minpoly_INL
