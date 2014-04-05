/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack/ffpack_ludivine.inl
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

#ifndef __FFLASFFPACK_ffpack_ludivine_INL
#define __FFLASFFPACK_ffpack_ludivine_INL


//#define LB_DEBUG
namespace FFPACK {
	template<class Field>
	inline size_t
	LUdivine_gauss( const Field& F, const FFLAS::FFLAS_DIAG Diag,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda, size_t*P,
			size_t *Q, const FFPACK::FFPACK_LUDIVINE_TAG LuTag)
	{
		size_t MN = std::min(M,N);
		typename Field::Element * Acurr = A;
		size_t r = 0;

		for (size_t k = 0; k < MN; ++k){
			size_t p = r;
			Acurr = A+k*lda+r;
			while ((p < N) && F.isZero (*(Acurr++)))
				p++;
			if (p < N){
				P[r] = p;
				if (r < k){
					FFLAS::fcopy (F, N-r, (A+k*lda+r),1, (A + r*(lda+1)), 1);
					Acurr = A+r+k*lda;
					for (size_t i=r; i<N; ++i)
						F.assign(*(Acurr++),F.zero);
				}

				FFLAS::fswap (F, M, A+r, lda, A+p, lda);
				Q[r] = k;
				r++;
			}
			if (k+1<M){
				ftrsv (F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit, r, A, lda, A+(k+1)*lda, 1);
				fgemv (F, FFLAS::FflasTrans, r, N-r, F.mOne, A+r, lda, A+(k+1)*lda, 1, F.one, A+(k+1)*lda+r, 1);
			}
			else
				return r;
		}

		return r;
	}

	template<class Element>
	class callLUdivine_small;



	template<class Field>
	inline size_t
	LUdivine_small( const Field& F, const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda, size_t*P,
			size_t *Q, const FFPACK::FFPACK_LUDIVINE_TAG LuTag)
	{
		return callLUdivine_small <typename Field::Element> ()
		(F, Diag, trans, M, N, A, lda, P, Q, LuTag);
	}

	template<class Element>
	class callLUdivine_small {
	public:
		template <class Field>
		inline size_t
		operator()( const Field& F, const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
			    const size_t M, const size_t N,
			    typename Field::Element * A, const size_t lda, size_t*P,
			    size_t *Q, const FFPACK::FFPACK_LUDIVINE_TAG LuTag)
		{

			if ( !(M && N) ) return 0;
			typedef typename Field::Element elt;
			elt * Aini = A;
			elt * Acurr;
			size_t rowp = 0;
			size_t R = 0;
			size_t k = 0;
			//size_t kmax = FFLAS::Protected::DotProdBound (F, 0, one) -1; // the max number of delayed operations
			while ((rowp<M) && (k<N)){
				size_t colp;

				//Find non zero pivot
				colp = k;
				Acurr = Aini;
				while ((F.isZero(*Acurr)) || (F.isZero (F.init (*Acurr, *Acurr))))
					if (++colp == N){
						if (rowp==M-1)
							break;
						colp=k; ++rowp;
						Acurr = Aini += lda;
					}
					else
						++Acurr;

				if ((rowp == M-1)&&(colp == N))
					break;
				R++;
				P[k] = colp;
				Q[k] = rowp;

				// Permutation of the pivot column
				FFLAS::fswap (F, M, A+k, lda, A + colp , lda);

				//Normalization
				elt invpiv;
				F.init(*Aini,*Aini);
				F.inv (invpiv,*Aini);

				for (size_t j=1; j<N-k; ++j)
					if (!F.isZero(*(Aini+j)))
						F.init(*(Aini+j), *(Aini+j));
				for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
					if (!F.isZero(*(Aini+i)))
						F.init(*(Aini+i), *(Aini+i));


				if (Diag == FFLAS::FflasUnit) {
					for (size_t j=1; j<N-k; ++j)
						if (!F.isZero(*(Aini+j)))
							F.mulin (*(Aini+j),invpiv);
				}
				else
					for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
						if (!F.isZero(*(Aini+i)))
							F.mulin (*(Aini+i),invpiv);

				//Elimination
				//Or equivalently, but without delayed ops :
				FFLAS::fger (F, M-rowp-1, N-k-1, F.mOne, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);

				Aini += lda+1; ++rowp; ++k;
			}

			// Compression the U matrix
			size_t l;
			if (Diag == FFLAS::FflasNonUnit){
				Aini = A;
				l = N;
			}
			else {
				Aini = A+1;
				l=N-1;
			}
			for (size_t i=0; i<R; ++i, Aini += lda+1) {
				if (Q[i] > i){
					FFLAS::fcopy (F, l-i, Aini+(Q[i]-i)*lda, 1, Aini, 1);
					for (size_t j=0; j<l-i; ++j)
						F.assign (*(Aini+(Q[i]-i)*lda+j), F.zero);
				}
			}
			return R;
		}
	};

	template<>
	class callLUdivine_small<double> {
	public:
		template <class Field>
		inline size_t
		operator()( const Field& F,
			    const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
			    const size_t M, const size_t N,
			    typename Field::Element * A, const size_t lda, size_t*P,
			    size_t *Q, const FFPACK::FFPACK_LUDIVINE_TAG LuTag)
		{

			if ( !(M && N) ) return 0;
			typedef typename Field::Element elt;
			elt * Aini = A;
			elt * Acurr;
			size_t rowp = 0;
			size_t R = 0;
			size_t k = 0;
			size_t delay =0;
			size_t kmax = FFLAS::Protected::DotProdBound (F, 0, F.one, FFLAS::FflasDouble) -1; // the max number of delayed operations
			while ((rowp<M) && (k<N)){
				size_t colp;

				//Find non zero pivot
				colp = k;
				Acurr = Aini;
				while ((F.isZero(*Acurr)) || (F.isZero (F.init (*Acurr, *Acurr))))
					if (++colp == N){
						if (rowp==M-1)
							break;
						colp=k; ++rowp;
						Acurr = Aini += lda;
					}
					else
						++Acurr;

				if ((rowp == M-1)&&(colp == N))
					break;
				R++;
				P[k] = colp;
				Q[k] = rowp;

				// Permutation of the pivot column
				FFLAS::fswap (F, M, A+k, lda, A + colp , lda);

				//Normalization
				elt invpiv;
				F.init(*Aini,*Aini);
				F.inv (invpiv,*Aini);

				for (size_t j=1; j<N-k; ++j)
					if (!F.isZero(*(Aini+j)))
						F.init(*(Aini+j), *(Aini+j));
				for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
					if (!F.isZero(*(Aini+i)))
						F.init(*(Aini+i), *(Aini+i));


				if (Diag == FFLAS::FflasUnit) {
					for (size_t j=1; j<N-k; ++j)
						if (!F.isZero(*(Aini+j)))
							F.mulin (*(Aini+j),invpiv);
				}
				else
					for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
						if (!F.isZero(*(Aini+i)))
							F.mulin (*(Aini+i),invpiv);

				if (delay++ >= kmax){ // Reduction has to be done
					delay = 0;
					for (size_t i=1; i<M-rowp; ++i)
						for (size_t j=1; j<N-k; ++j)
							F.init(	*(Aini+i*lda+j),*(Aini+i*lda+j));
				}
				//Elimination
				for (size_t i=1; i<M-rowp; ++i)
					for (size_t j=1; j<N-k; ++j)
						*(Aini+i*lda+j) -= *(Aini+i*lda) * *(Aini+j);
				//Or equivalently, but without delayed ops :
				//FFLAS::fger (F, M-rowp-1, N-k-1, F.mOne, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);

				Aini += lda+1; ++rowp; ++k;
			}

			// Compression the U matrix
			size_t l;
			if (Diag == FFLAS::FflasNonUnit){
				Aini = A;
				l = N;
			}
			else {
				Aini = A+1;
				l=N-1;
			}
			for (size_t i=0; i<R; ++i, Aini += lda+1) {
				if (Q[i] > i){
					FFLAS::fcopy (F, l-i, Aini+(Q[i]-i)*lda, 1, Aini, 1);
					for (size_t j=0; j<l-i; ++j)
						F.assign (*(Aini+(Q[i]-i)*lda+j), F.zero);
				}
			}
			return R;
		}
	};

	template<>
	class callLUdivine_small<float> {
	public:
		template <class Field>
		inline size_t
		operator()( const Field& F,
			    const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
			    const size_t M, const size_t N,
			    typename Field::Element * A, const size_t lda, size_t*P,
			    size_t *Q, const FFPACK::FFPACK_LUDIVINE_TAG LuTag)
		{

			if ( !(M && N) ) return 0;
			typedef typename Field::Element elt;
			elt * Aini = A;
			elt * Acurr;
			size_t rowp = 0;
			size_t R = 0;
			size_t k = 0;
			size_t delay =0;
			size_t kmax = FFLAS::Protected::DotProdBound (F, 0, F.one, FFLAS::FflasFloat) -1; // the max number of delayed operations
			while ((rowp<M) && (k<N)){
				size_t colp;

				//Find non zero pivot
				colp = k;
				Acurr = Aini;
				while ((F.isZero(*Acurr)) || (F.isZero (F.init (*Acurr, *Acurr))))
					if (++colp == N){
						if (rowp==M-1)
							break;
						colp=k; ++rowp;
						Acurr = Aini += lda;
					}
					else
						++Acurr;

				if ((rowp == M-1)&&(colp == N))
					break;
				R++;
				P[k] = colp;
				Q[k] = rowp;

				// Permutation of the pivot column
				FFLAS::fswap (F, M, A+k, lda, A + colp , lda);

				//Normalization
				elt invpiv;
				F.init(*Aini,*Aini);
				F.inv (invpiv,*Aini);

				for (size_t j=1; j<N-k; ++j)
					if (!F.isZero(*(Aini+j)))
						F.init(*(Aini+j), *(Aini+j));
				for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
					if (!F.isZero(*(Aini+i)))
						F.init(*(Aini+i), *(Aini+i));


				if (Diag == FFLAS::FflasUnit) {
					for (size_t j=1; j<N-k; ++j)
						if (!F.isZero(*(Aini+j)))
							F.mulin (*(Aini+j),invpiv);
				}
				else
					for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
						if (!F.isZero(*(Aini+i)))
							F.mulin (*(Aini+i),invpiv);

				if (delay++ >= kmax){ // Reduction has to be done
					delay = 0;
					for (size_t i=1; i<M-rowp; ++i)
						for (size_t j=1; j<N-k; ++j)
							F.init(	*(Aini+i*lda+j),*(Aini+i*lda+j));
				}
				//Elimination
				for (size_t i=1; i<M-rowp; ++i)
					for (size_t j=1; j<N-k; ++j)
						*(Aini+i*lda+j) -= *(Aini+i*lda) * *(Aini+j);
				//Or equivalently, but without delayed ops :
				//FFLAS::fger (F, M-rowp-1, N-k-1, F.mOne, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);

				Aini += lda+1; ++rowp; ++k;
			}

			// Compression the U matrix
			size_t l;
			if (Diag == FFLAS::FflasNonUnit){
				Aini = A;
				l = N;
			}
			else {
				Aini = A+1;
				l=N-1;
			}
			for (size_t i=0; i<R; ++i, Aini += lda+1) {
				if (Q[i] > i){
					FFLAS::fcopy (F, l-i, Aini+(Q[i]-i)*lda, 1, Aini, 1);
					for (size_t j=0; j<l-i; ++j)
						F.assign (*(Aini+(Q[i]-i)*lda+j), F.zero);
				}
			}
			return R;
		}
	};

	template <class Field>
	inline size_t
	LUdivine (const Field& F,
		  const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
		  const size_t M, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  size_t*P, size_t *Q
		  , const FFPACK::FFPACK_LUDIVINE_TAG LuTag // =FFPACK::FfpackLQUP
		  , const size_t cutoff // =__FFPACK_LUDIVINE_CUTOFF
		 )
	{

		if ( !(M && N) ) return 0;
		typedef typename Field::Element elt;
		size_t MN = std::min(M,N);

		size_t incRow, incCol, rowDim, colDim;
		if (trans == FFLAS::FflasTrans){
			incRow = 1;
			incCol = lda;
			colDim = M;
			rowDim = N;
		}
		else {
			incRow = lda;
			incCol = 1;
			colDim = N;
			rowDim = M;
		}

		if ((rowDim < cutoff) && (colDim < 2*cutoff)) { // the coeff 2 is experimentally determined!
			return LUdivine_small (F, Diag, trans, M, N, A, lda, P, Q, LuTag);
		}
		else { // recursively :
			if (MN == 1){
				size_t ip=0;
				//while (ip<N && !F.isUnit(*(A+ip)))ip++;
				while (F.isZero (*(A+ip*incCol)))
					if (++ip == colDim)
						break;
				*Q=0;
				if (ip == colDim){ // current row is zero
					*P=0;
					if (colDim == 1){
						//while (ip<M && !F.isUnit(*(A+ip*lda)))
						while (ip<rowDim && F.isZero(*(A + ip*incRow))){
							Q[ip]=ip;
							ip++;
						}
						if (ip == rowDim) {
							return 0;
						}
						else {
							size_t oldip = ip;
							if ( Diag == FFLAS::FflasNonUnit ){
								elt invpiv;
								F.inv(invpiv,*(A+ip*incRow));
								while(++ip<rowDim)
									F.mulin(*(A + ip*incRow), invpiv);
								elt tmp;
								F.assign(tmp, *(A+oldip*incRow));
								F.assign( *(A+oldip*incRow), *A);
								F.assign( *A, tmp);
							}
							*Q=oldip;

							return 1;
						}
					}
					else{
					       	*Q=0; return 0;
					}
				}
				*P=ip;
				if (ip!=0){
					// swap the pivot
					typename Field::Element tmp=*A;
					*A = *(A + ip*incCol);
					*(A + ip*incCol) = tmp;
				}
				elt invpiv;
				F.inv(invpiv, *A);
				if ( Diag == FFLAS::FflasUnit ){
					// Normalisation of the row
					for (size_t k=1; k<colDim; k++)
						F.mulin(*(A+k*incCol), invpiv);
				}
				else  {
					if ( colDim==1 )
					while(++ip<rowDim)
						F.mulin(*(A + ip*incRow), invpiv);
				}
				return 1;
			}
			else { // MN>1
				size_t Nup = rowDim >> 1;
				size_t Ndown =  rowDim - Nup;
				// FFLASFFPACK_check(Ndown < rowDim);
				// Recursive call on NW
				size_t R, R2;
				if (trans == FFLAS::FflasTrans){
					R = LUdivine (F, Diag, trans, colDim, Nup, A, lda, P, Q,
						      LuTag, cutoff);

					typename Field::Element *Ar = A + Nup*incRow;   // SW
					typename Field::Element *Ac = A + R*incCol;     // NE
					typename Field::Element *An = Ar+ R*incCol;     // SE

					if (!R){
						if (LuTag == FFPACK::FfpackSingular )
							return 0;
					}
					else {
						FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
								Ndown, 0,(int) R, Ar, lda, P);
						// Ar <- L1^-1 Ar
						ftrsm( F, FFLAS::FflasLeft, FFLAS::FflasLower,
						       FFLAS::FflasNoTrans, Diag, R, Ndown,
						       F.one, A, lda, Ar, lda);
						// An <- An - Ac*Ar
						fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, colDim-R, Ndown, R,
						       F.mOne, Ac, lda, Ar, lda, F.one, An, lda);
					}
					// Recursive call on SE
					R2 = LUdivine (F, Diag, trans, colDim-R, Ndown, An, lda, P + R, Q + Nup, LuTag, cutoff);
					for (size_t i = R; i < R + R2; ++i)
						P[i] += R;
					if (R2) {
						// An <- An.P2
						FFPACK::applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
								Nup,(int) R, (int)(R+R2), A, lda, P);
					}
					else {
						if (LuTag == FFPACK::FfpackSingular)
							return 0;
					}

				}
				else { // trans == FFLAS::FflasNoTrans
					R = LUdivine (F, Diag, trans, Nup, colDim, A, lda, P, Q, LuTag, cutoff);
					typename Field::Element *Ar = A + Nup*incRow;   // SW
					typename Field::Element *Ac = A + R*incCol;     // NE
					typename Field::Element *An = Ar+ R*incCol;     // SE


					if (!R){
						if (LuTag == FFPACK::FfpackSingular )
							return 0;
					}
					else { /*  R>0 */
						// Ar <- Ar.P
						FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
								Ndown, 0,(int) R, Ar, lda, P);
						// Ar <- Ar.U1^-1
						ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper,
						       FFLAS::FflasNoTrans, Diag, Ndown, R,
						       F.one, A, lda, Ar, lda);
						// An <- An - Ar*Ac
						fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ndown, colDim-R, R,
						       F.mOne, Ar, lda, Ac, lda, F.one, An, lda );

					}
					// Recursive call on SE
					R2=LUdivine (F, Diag, trans, Ndown, N-R, An, lda,P+R, Q+Nup, LuTag, cutoff);
					for (size_t i = R; i < R + R2; ++i)
						P[i] += R;
					if (R2)
						// An <- An.P2
						FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
								Nup,(int) R, (int)(R+R2), A, lda, P);
					else if (LuTag == FFPACK::FfpackSingular)
						return 0;

				}
				// Non zero row permutations
				for (size_t i = Nup; i < Nup + R2; i++)
					Q[i] += Nup;
				if (R < Nup){
					// Permutation of the 0 rows
					if (Diag == FFLAS::FflasNonUnit){
						for ( size_t i = Nup, j = R ; i < Nup + R2; ++i, ++j){
							FFLAS::fcopy( F, colDim - j, A + i*incRow + j*incCol, incCol, A + j * (lda + 1), incCol);
							for (typename Field::Element *Ai = A + i*incRow + j*incCol;
							     Ai != A + i*incRow + colDim*incCol; Ai+=incCol)
								F.assign (*Ai, F.zero);
							///@todo std::swap ?
							size_t t = Q[j];
							Q[j]=Q[i];
							Q[i] = t;
						}
					}
					else { // Diag == FFLAS::FflasUnit
						for ( size_t i = Nup, j = R+1 ; i < Nup + R2; ++i, ++j){
							FFLAS::fcopy( F, colDim - j,
								      A + i*incRow + j*incCol, incCol,
								      A + (j-1)*incRow + j*incCol, incCol);
							for (typename Field::Element *Ai = A + i*incRow + j*incCol;
							     Ai != A + i*incRow + colDim*incCol; Ai+=incCol)
								F.assign (*Ai, F.zero);
							size_t t = Q[j-1];
							Q[j-1]=Q[i];
							Q[i] = t;
						}
					}
				}
				return R + R2;
			}
		}
	}

	namespace Protected {

		//---------------------------------------------------------------------
		// LUdivine_construct: (Specialisation of LUdivine)
		// LUP factorisation of the Krylov base matrix of A^t and v.
		// When all rows have been factorized in A, and rank is full,
		// then new krylov vectors are computed and then triangularized
		// P is the permutation matrix stored in the lapack style
		// nRowX is the number of Krylov vectors already computed,
		// nUsedRowX is the number of Krylov vectors already triangularized
		//---------------------------------------------------------------------

		template <class Field>
		size_t
		LUdivine_construct( const Field& F, const FFLAS::FFLAS_DIAG Diag,
				    const size_t M, const size_t N,
				    const typename Field::Element * A, const size_t lda,
				    typename Field::Element * X, const size_t ldx,
				    typename Field::Element * u, size_t* P,
				    bool computeX
				    , const FFPACK::FFPACK_MINPOLY_TAG MinTag //= FFPACK::FfpackDense
				    , const size_t kg_mc// =0
				    , const size_t kg_mb// =0
				    , const size_t kg_j // =0
				  )
		{

			size_t MN = std::min(M,N);

			if (MN == 1){
				size_t ip=0;
				while (ip<N && F.isZero(*(X+ip))){ip++;}
				if (ip==N){ // current row is zero
					*P=0;
					return 0;
				}
				*P=ip;
				if (ip!=0){
					// swap the pivot
					typename Field::Element tmp=*X;
					*X = *(X+ip);
					*(X+ip) = tmp;
				}
				if ( Diag == FFLAS::FflasUnit ){
					typename Field::Element invpiv;
					F.inv(invpiv, *X);

					// Normalisation of the row
					for (size_t k=1; k<N; k++)
						F.mulin(*(X+k), invpiv);
				}
				if (N==1 && M>1 && computeX)// Only appends when A is 1 by 1
					F.mul(*(X+ldx),*X, *A);

				return 1;
			}
			else{ // MN>1
				size_t Nup = MN>>1;
				size_t Ndown =  M - Nup;

				// Recursive call on NW
				size_t R = LUdivine_construct(F, Diag, Nup, N, A, lda, X, ldx, u,
							      P, computeX, MinTag, kg_mc, kg_mb, kg_j );
				if (R==Nup){
					typename Field::Element * Xr = X + Nup*ldx; //  SW
					typename Field::Element * Xc = X + Nup;     //  NE
					typename Field::Element * Xn = Xr + Nup;    //  SE
					typename Field::Element * Xi = Xr;
					if ( computeX ){
						if (MinTag == FFPACK::FfpackDense)
							for (size_t i=0; i< Ndown; ++i, Xi+=ldx){
								fgemv(F, FFLAS::FflasNoTrans, N, N, F.one,
								      A, lda, u, 1, F.zero, Xi,1);
								FFLAS::fcopy(F, N,Xi, 1, u,1);
							}
						else // Keller-Gehrig Fast algorithm's matrix
							for (size_t i=0; i< Ndown; ++i, Xi+=ldx){
								FFPACK::Protected::fgemv_kgf( F, N, A, lda, u, 1, Xi, 1,
											      kg_mc, kg_mb, kg_j );
								FFLAS::fcopy(F, N,Xi, 1, u,1);
							}
					}
					// Apply the permutation on SW
					FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
							Ndown, 0,(int) R, Xr, ldx, P);
					// Triangular block inversion of NW and apply to SW
					// Xr <- Xr.U1^-1
					ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag,
					       Ndown, R, F.one, X, ldx, Xr, ldx);

					// Update of SE
					// Xn <- Xn - Xr*Xc
					fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ndown, N-Nup, Nup,
					       F.mOne, Xr, ldx, Xc, ldx, F.one, Xn, ldx);

					// Recursive call on SE

					size_t R2 = LUdivine_construct(F, Diag, Ndown, N-Nup, A, lda,
								       Xn, ldx, u, P + Nup,
								       false, MinTag, kg_mc, kg_mb, kg_j);
					for ( size_t i=R;i<R+R2;++i) P[i] += R;

					FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
							Nup, (int)R, (int)(R+R2), X, ldx, P);

					return R+=R2;
				}
				else
					return R;
				// Rank deficient matrices can only be factorized
				// under the condition: the first R rows are linearly independent
				// If not, the lower block is never factorized as soon as the
				// upper block is rank defficient
			}
		}

	} // Protected

	//---------------------------------------------------------------------
	// TURBO: rank computation algorithm
	//---------------------------------------------------------------------

	template <class Field>
	inline size_t
	TURBO (const Field& F, const size_t M, const size_t N,
	       typename Field::Element* A, const size_t lda, size_t * P, size_t * Q, const size_t cutoff)
	{

		size_t mo2 = (M>>1);
		size_t no2 = (N>>1);

		typename Field::Element * NW = A;
		typename Field::Element * NE = A + no2;
		typename Field::Element * SW = A + mo2*lda;
		typename Field::Element * SE = SW + no2;

		size_t ld1, ld2, ld3, ld4;
		ld1 = ld2 = ld3 = ld4 = lda;

		if ( !(M && N) ) return 0;

		// Column permutation
		size_t * P1 = new size_t[no2];
		size_t * P2 = new size_t[N-no2];
		// Row Permutation
		size_t * Q1 = new size_t[mo2];
		size_t * Q2 = new size_t[M-mo2];
		for (size_t i=0; i<mo2; ++i)
			Q1[i] = 0;
		for (size_t i=0; i<M-mo2; ++i)
			Q2[i] = 0;
		size_t q1,q2,q3,q3b,q4;
		/*q1=q2=q3=  */q3b=q4=0;


		// Step 1: NW = L1.Q1.U1.P1
		size_t mloc = mo2;
		size_t nloc ;
#if 0
		Timer tim;
		tim.clear();
		tim.start();
#endif
		q1 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, no2, NW, ld1, P1, Q1, FFPACK::FfpackLQUP, cutoff);

#if 0
		tim.stop();
		cerr<<"LQUP1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif
#if LB_DEBUG
		std::cerr<<"NW= L1.Q1.U1.P1"<<std::endl;
		write_field(F,std::cerr,NW,M,N,lda);
#endif
		// B1 = L^-1.NE
#ifdef LB_DEBUG
		std::cerr<<"avant B1 = L^-1.NE"<<std::endl;
		write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif
		solveLB( F, FFLAS::FflasLeft, mo2, N-no2, q1, NW, ld1, Q1, NE, ld2);
#ifdef LB_DEBUG
		std::cerr<<"B1 = L^-1.NE"<<std::endl;
		write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif

		// NE = Q^-1.NE

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
				N-no2, 0,(int) mo2, NE, ld2, Q1);
#ifdef LB_DEBUG
		std::cerr<<"NE=Q^-1.NE"<<std::endl;
		write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif

		// SW = SW.P1
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
				M-mo2, 0,(int) q1, SW, ld3, P1 );
#ifdef LB_DEBUG
		std::cerr<<"SW = SW.P1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif

#if 0
		tim.stop();
		std::cerr<<"L^-1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// N1 = SW_{1,q1} . U1^-1
		ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, M-mo2, q1, F.one, NW, ld1 , SW, ld3 );
#ifdef LB_DEBUG
		std::cerr<<" N1 = SW_{1,q1} . U1^-1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif

#if 0
		tim.stop();
		std::cerr<<"trsm:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// I1 = SW_{q1+1,n} - N1.G1
		fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-mo2,  no2-q1, q1, F.mOne, SW, ld3, NW+q1, ld1, F.one, SW+q1, ld3);
#ifdef LB_DEBUG
		std::cerr<<" I1 = SW_{q1+1,n} - N1.G1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif

#if 0
		tim.stop();
		std::cerr<<"fgemm1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// E1 = SE - N1.B1_{1,q1}
		fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-mo2, N-no2, q1, F.mOne, SW, ld3, NE, ld2, F.one, SE, ld4);
#ifdef LB_DEBUG
		std::cerr<<"  E1 = SE - N1.B1_{1,q1}"<<std::endl;
		write_field(F,std::cerr,SE,M-mo2,N-no2,ld4);
#endif

#if 0
		tim.stop();
		std::cerr<<"fgemm2:"<<tim.realtime()<<std::endl;
		tim.start();
#endif


		//Step 2: E1 = L2.Q2.U2.P2
		mloc = M-mo2;
		nloc = N-no2;
		q2 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, SE, ld4, P2, Q2, FFPACK::FfpackLQUP, cutoff);
#ifdef LB_DEBUG
		std::cerr<<"  E1 = L2.Q2.U2.P2"<<std::endl;
		write_field(F,std::cerr,SE,M-mo2,N-no2,ld4);
#endif

#if 0
		tim.stop();
		std::cerr<<"LQUP2:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// [I2;F2] = L2^-1.I1
		solveLB( F, FFLAS::FflasLeft, mloc, no2-q1, q2, SE, ld4, Q2, SW+q1, ld3);
#ifdef LB_DEBUG
		std::cerr<<"  [I2;F2] = L2^-1.I1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif
		// I1 = Q2^-1.I1
		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
				no2-q1, 0,(int) mloc, SW+q1, ld3, Q2 );
#ifdef LB_DEBUG
		std::cerr<<"I1 = Q2^-1.I1"<<std::endl;
		write_field(F,std::cerr,SW,mloc,no2,ld3);
#endif

		// B1 = B1.P2
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
				mo2, 0,(int) q2, NE, ld2, P2 );
#ifdef LB_DEBUG
		std::cerr<<"B1 = B1.P2"<<std::endl;
		write_field(F,std::cerr,NE,mo2,N-no2,ld2);
#endif
		// Updating P
#if 0
		for (size_t i=no2;i<N;++i)
			P[i] += no2;
		tim.stop();
		std::cerr<<"L2^-1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		//alternative: de 0 a q2 avant
		// N2 = B1_{q1+1,mo2} . V2^-1
		ftrsm(F, FFLAS::FflasRight, FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit, mo2-q1, q2, F.one, SE, ld4, NE+q1*ld2,ld2);
#if 0
		tim.stop();
		std::cerr<<"trsm2:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// H2 = B1_{q1+1,mo2;q2,N-no2} - N2.E2
		fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mo2-q1, N-no2-q2, q2, F.mOne, NE+q1*ld2, ld2, SE+q2, ld4, F.one, NE+q1*ld2+q2, ld2);

#if 0
		tim.stop();
		std::cerr<<"fgemm12:"<<tim.realtime()<<std::endl;
		tim.start();
		O2 = NW_{q1+1,mo2;q1+1,N-no2} = - N2.I2
		write_field (F,cerr<<"avant O2"<<endl, A, M, N, lda);
#endif

		fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mo2-q1, no2-q1, q2, F.mOne, NE+q1*ld2, ld2, SW+q1, ld3, F.zero,
		      NW+q1*(ld1+1), ld1);
		//	write_field (F,cerr<<"apres O2"<<endl, A, M, N, lda);
#if 0
		tim.stop();
		std::cerr<<"fgemm22:"<<tim.realtime()<<std::endl;
		tim.start();
#endif


		//Step 3: F2 = L3.Q3.U3.P3
		mloc = M-mo2-q2;
		nloc = no2-q1;
		q3 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, SW+q2*ld3+q1, ld3, P1+q1, Q2+q2, FFPACK::FfpackLQUP, cutoff);

		// Updating P1,Q2
		for (size_t i=q1;i<no2;++i)
			P1[i] += q1;
		for (size_t i=q2;i<q2+q3;++i)
			Q2[i] += q2;

		//Step 3bis: H2 = L3b.Q3b.U3b.P3b
		mloc = mo2-q1;
		nloc = N-no2-q2;

		q3b = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, NE+q1*ld2+q2, ld2, P2+q2, Q1+q1, FFPACK::FfpackLQUP, cutoff);

		// Updating P2, Q1
		for (size_t i = q2; i < q2+q3b; ++i)
			P2[i] += q2;

#if 0
		tim.stop();
		std::cerr<<"LQUP3et3bis:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		if (( q3 < no2-q1) && (q3b<mo2-q1)){

			// [O3;_] = L3b^-1.O2
			if (q3b>0){
#if 0
				if ( mo2-q1 < N-no2-q2+q1)
					// L is expanded to a Lower triangular matrix
					solveLB( F, FFLAS::FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, rP3b, NW+q1*(ld1+1), ld1);
				else
#endif
					//std::cerr<<"USING SOLVELB2"<<std::endl;
					//no modification of L
					solveLB2( F, FFLAS::FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, Q1+q1, NW+q1*(ld1+1), ld1);
#ifdef LB_DEBUG
				std::cerr<<"O2 avant="<<std::endl;
				write_field(F,std::cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif

				// O2 = Q3b^-1.O2
				FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
						no2-q1, 0,(int) mloc, NW+q1*(ld1+1), ld1, Q1+q1 );
#ifdef LB_DEBUG
				std::cerr<<"O2 apres="<<std::endl;
				write_field(F,std::cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif

				//updating Q
#if 0
				size_t tmp;
				for (size_t j=0;j<mo2-q1;++j)
					if (rP3b[j]!=j){
						//	std::cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<std::endl;
						tmp = Q[j+q1];
						Q[j+q1] = Q[rP3b[j]+q1];
						Q[rP3b[j]+q1] = tmp;
					}
#endif

				// X2 = X2.P3
				// Si plusieurs niveaux rec, remplacer X2 par [NW;I2]
				FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
						mo2-q1-q3b,(int) q1, (int)(q1+q3),
						NW/*+(q1+q3b)*ld1*/, ld1, P1);
				FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
						q2,(int) q1, (int)(q1+q3),
						SW/*+(q1+q3b)*ld1*/, ld3, P1);


				// A faire si plusieurs niveaux recursifs
				// B2 = B2.P3b
				FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
						q1,(int) q2, (int)(q2+q3b),
						NW, ld2, P2);
				//flaswp(F,q1,NE,lda,no2+q2,no2+q2+q3b,P,1);
				// E2 = E2.P3b
				FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
						q2,(int) q2, (int)(q2+q3b),
						SE, ld4, P2);
				//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1);
			}

			// N3 = X2 . D3^-1
			ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, mo2-q1-q3b, q3, F.one, SW+q2*ld3+q1, ld3 ,NW+(q1+q3b)*ld1+q1,ld1);

			// T2 = T2 - N3.F3
			fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mo2-q1-q3b, no2-q1-q3,q3, F.mOne, NW+(q1+q3b)*ld1+q1, ld1, SW+q2*ld3+q3+q1, ld3, F.one, NW+(q1+q3b)*ld1+q1+q3, ld1 );


			//Step 4: T2 = L4.Q4.U4.P4
			mloc = mo2-q1-q3b;
			nloc = no2-q1-q3;

#if 0
			size_t * rP4 = new size_t[mloc];
			for (size_t j=0;j<mo2-q1;++j)
				rP4[j]=0;
#endif
			q4 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, NW+(q1+q3b)*ld1+q1+q3, ld1, P1+q1+q3, Q1+q1+q3b, FFPACK::FfpackLQUP, cutoff);

			// Updating P
			for (size_t i=q1+q3;i<q1+q3+q4;++i)
				P1[i] += q3;

#if 0
			size_t tmp;
			if (rP4[j]!=j){
				//	std::cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<std::endl;
				tmp = Q[j+q1+q3b];
				Q[j+q1+q3b] = Q[rP3b[j]+q1+q3b];
				Q[rP3b[j]+q1+q3b] = tmp;
			}
#endif

			// A faire si plusieurs niveaux recursifs
			// [G1;O3] = [G1;O3].P4
			FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					q1+q3b, (int)(q1+q3), (int)(q1+q3+q4),
					NW, ld1, P1);
			//flaswp(F,q1+q3b,NE,lda,no2+q2,no2+q2+q3b,P,1);
			// [I2;F3] = [I2;F3].P4
			FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					q2+q3, (int)(q1+q3),(int) (q1+q3+q4),
					SW, ld3, P1);
			//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1);
		}
		//!!!!!! Attention a appliquer Q4, Q2, Q3, Q3b a gauche !!!!!!!

		//updating Q1
		for (size_t i = q1; i < q1+q3b; ++i)
			Q1[i] += q1;
		for (size_t i=q1+q3b;i<q1+q3b+q4;++i)
			Q1[i] += q1 + q3b;

		for (size_t i=0; i<q1; ++i)
			P[i] = P1[i];
		for (size_t i=q1; i<q1+q2; ++i)
			P[i] = P2[i-q1] + no2;
		for (size_t i=q1+q2; i<q1+q2+q3; ++i)
			P[i] = P1[i-q2];
		for (size_t i=q1+q2+q3; i<q1+q2+q3+q3b; ++i)
			P[i] = P2[i-q1-q3]+no2;
		for (size_t i=q1+q2+q3+q3b; i<q1+q2+q3+q3b+q4; ++i)
			P[i] = P1[i-q2-q3b];
		delete[] P1;
		delete[] P2;

		for (size_t i=0; i<q1; ++i)
			Q[i] = Q1[i];
		for (size_t i=q1; i<q1+q2; ++i)
			Q[i] = Q2[i-q1] + mo2;
		for (size_t i=q1+q2; i<q1+q2+q3; ++i)
			Q[i] = Q2[i-q1] + mo2;
		for (size_t i=q1+q2+q3; i<q1+q2+q3+q3b; ++i)
			Q[i] = Q1[i-q2-q3];
		for (size_t i=q1+q2+q3+q3b; i<q1+q2+q3+q3b+q4; ++i)
			P[i] = Q1[i-q2-q3];
		delete[] Q1;
		delete[] Q2;


		//write_field (F, cerr<<"avant reordonnancement"<<endl, A, M,N, lda)<<endl;
		typename Field::Element * R = new typename Field::Element[M*N];
		size_t ldr = N;
		// Copying first q1 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q1, NW+i*ld1,1, R+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q1, SW+(i-q1)*ld3,1, R+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q1, NW+(i-q2-q3)*ld1,1, R+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q1, SW+(i-mo2)*ld3,1, R+i*ldr, 1);
		// Copying q1..q2 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q2, NE+i*ld2,1, R+q1+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q2, SE+(i-q1)*ld4,1, R+q1+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q2, NE+(i-q2-q3)*ld2,1, R+q1+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q2, SE+(i-mo2)*ld4,1, R+q1+i*ldr, 1);
		// Copying q2..q3 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q3, NW+q1+i*ld1,1, R+q1+q2+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q3, SW+q1+(i-q1)*ld3,1, R+q1+q2+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q3, NW+q1+(i-q2-q3)*ld1,1, R+q1+q2+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q3, SW+q1+(i-mo2)*ld3,1, R+q1+q2+i*ldr, 1);
		// Copying q3..q3b cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q3b, NE+q2+i*ld2,1, R+q1+q2+q3+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q3b, SE+q2+(i-q1)*ld4,1, R+q1+q2+q3+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q3b, NE+q2+(i-q2-q3)*ld2,1, R+q1+q2+q3+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q3b, SE+q2+(i-mo2)*ld4,1, R+q1+q2+q3+i*ldr, 1);
		// Copying q3b..q4 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q4, NW+q1+q3+i*ld1,1, R+q1+q2+q3+q3b+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q4, SW+q1+q3+(i-q1)*ld3,1, R+q1+q2+q3+q3b+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q4, NW+q1+q3+(i-q2-q3)*ld1,1, R+q1+q2+q3+q3b+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q4, SW+q1+q3+(i-mo2)*ld3,1, R+q1+q2+q3+q3b+i*ldr, 1);
		// Copying the last cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, NW+q1+q3+q4+i*ld1,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, SW+q1+q3+q4+(i-q1)*ld3,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, NW+q1+q3+q4+(i-q2-q3)*ld1,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, SW+q1+q3+q4+(i-mo2)*ld3,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		// Copying the last cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, NE+q2+q3b+i*ld2,1, R+no2+q2+q3b+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, SE+q2+q3b+(i-q1)*ld4,1, R+no2+q2+q3b+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, NE+q2+q3b+(i-q2-q3)*ld2,1, R+no2+q2+q3b+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, SE+q2+q3b+(i-mo2)*ld4,1, R+no2+q2+q3b+i*ldr, 1);

		// A=R : to be improved (avoid allocation of R). To be changed if rec data structure are used
		for (size_t i=0; i<M; ++i)
			FFLAS::fcopy (F, N, R+i*ldr,1, A+i*lda, 1);

		delete[] R;
		//delete[] Q;
		// Necessaire:
		// 1 traiter les flaswp manquants
		// Facultatif:
		// 2 permutations de lignes doivent etre coherentes
		// 3 effectuer les dernieres permutations lignes et colonnes
		//std::cerr<<q1<<" "<<q2<<" "<<q3<<" "<<q3b<<" "<<q4<<std::endl;
		return q1+q2+q3+q3b+q4;
	}


	/*!
	 * @brief Updates an existing LU factorisation with more rows.
	 *
	 * @param F Field on which arithmetic is done
	 * @param Diag Is \p L unit ? (if so, \c FFLAS::FflasUnit)
	 * @param trans Not used yet, should be \c FFLAS::FflasNoTrans
	 * @param M rows in \p A
	 * @param N cols in \p A
	 * @param A \p A is already in \c LU factorisation
	 * @param lda leading dimension of \p A
	 * @param R rank of \p A
	 * @param K rows in \p B
	 * @param B more rows to append to \p A
	 * @param ldb leading dimension of \p B (not tested if != lda)
	 * @param P permutation for \c LU in \p A. Should be big enough so it can store the permutation for \c LU of \p A and \p B
	 * @param Q same as \p P
	 * @param LuTag see \c LUdivine.
	 * @param cutoff see \c LUdivine.
	 *
	 * @return rank of <code>A.append(B)</code>
	 * @bug may be bogus.
	 */
	template <class Field>
	size_t LUpdate (const Field& F,
			const FFLAS::FFLAS_DIAG Diag, const FFLAS::FFLAS_TRANSPOSE trans,
			const size_t M, const size_t N,
			typename Field::Element * A, const size_t lda,
			const size_t R,
			const size_t K,
			typename Field::Element * B, const size_t ldb,
			size_t*P, size_t *Q , const FFPACK::FFPACK_LUDIVINE_TAG LuTag , const size_t cutoff)
	{
		if (trans == FFLAS::FflasTrans)
			throw Failure(__func__,__FILE__,__LINE__,"Transposed version is not implemented yet");
		// size_t MN = std::min(M,N);

		size_t incRow, incCol, rowDim, colDim;
#if 0 /*  not working */
		if (trans == FFLAS::FflasTrans){
			incRow = 1;
			incCol = lda;
			colDim = M;
			rowDim = N;
		}
		else
#endif
		{ // trans == FFLAS::FflasNoTrans
			incRow = lda;
			incCol = 1;
			colDim = N;
			rowDim = M;
		}

		size_t Nup   = rowDim;
		size_t Ndown = K;


		if ( !K ) {  // no line to append
			std::cout << "no row to append" << std::endl;
			return R;
		}
		if ( !R ) { // A was null
			// std::cout << "A was 0" << std::endl;
			// LU on B
			size_t R2 = LUdivine(F,Diag,trans,K,N,B,lda,P,Q+Nup,LuTag,cutoff);
			if (!R2)
				return 0 ;
			// Non zero row permutations
			for (size_t i = Nup; i < Nup + R2; ++i)
				Q[i] += Nup;
			{ // Permutation of the 0 rows Could probably be improved !
				if (Diag == FFLAS::FflasNonUnit){
					for ( size_t i = 0, j = R ; i < R2; ++i, ++j){
						FFLAS::fcopy( F, colDim - j, B + i*incRow + j*incCol, incCol, A + j * (lda + 1),
							      incCol);
						typename Field::Element *Ai = B + i*incRow + j*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						///@todo std::swap ?
						size_t t = Q[j];
						Q[j]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
				else { // Diag == FFLAS::FflasUnit
					for ( size_t i = 0, ii = R+1 ; i < R2; ++i, ++ii){
						if (ii < M)
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      A + (ii-1)*incRow + ii*incCol, incCol);
						else {
							std::cout << "dangerous zone" << std::endl;
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      B + (ii-M-1)*incRow + ii*incCol, incCol);
						}

						typename Field::Element *Ai   = B + i*incRow + ii*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						size_t t = Q[ii-1];
						Q[ii-1]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
			}

			return R2 ;
		}


		{
			// Recursive call on NW
			size_t R2;
			{ // trans == FFLAS::FflasNoTrans
				typename Field::Element *Ac = A  + R*incCol;    // NE
				typename Field::Element *Ar = B;                // SW
				typename Field::Element *An = Ar + R*incCol;    // SE


				// Ar <- Ar.P
				FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
						Ndown, 0,(int) R, Ar, lda, P);
				// Ar <- Ar.U1^-1
				ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper,
				       FFLAS::FflasNoTrans, Diag, Ndown, R,
				       F.one, A, lda, Ar, lda);
				// An <- An - Ar*Ac
				fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ndown, colDim-R, R,
				       F.mOne, Ar, lda, Ac, lda, F.one, An, lda);

				// LU call on SE
				R2=LUdivine (F, Diag, trans, Ndown, N-R, An, lda,P+R, Q+Nup,
					     LuTag, cutoff);
				if (R2) {
					for (size_t i = R; i < R + R2; ++i)
						P[i] += R;
					// An <- An.P2
					FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
							Nup,(int) R, (int)(R+R2), A, lda, P);

				}
				else {
					if (LuTag == FFPACK::FfpackSingular)
						return 0;
				}

			}
			// Non zero row permutations
			for (size_t i = Nup; i < Nup + R2; ++i)
				Q[i] += Nup;
			if (R < Nup){ // Permutation of the 0 rows
				if (Diag == FFLAS::FflasNonUnit){
					for ( size_t i = 0, j = R ; i < R2; ++i, ++j){
						FFLAS::fcopy( F, colDim - j, B + i*incRow + j*incCol, incCol, A + j * (lda + 1),
							      incCol);
						typename Field::Element *Ai = B + i*incRow + j*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						///@todo std::swap ?
						size_t t = Q[j];
						Q[j]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
				else { // Diag == FFLAS::FflasUnit
					for ( size_t i = 0, ii = R+1 ; i < R2; ++i, ++ii){
						if (ii < M)
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      A + (ii-1)*incRow + ii*incCol, incCol);
						else {
							// std::cout << "dangerous zone" << std::endl;
							FFLAS::fcopy( F, colDim - ii,
								      B + i*incRow + ii*incCol, incCol,
								      B + (ii-M-1)*incRow + ii*incCol, incCol);
						}

						typename Field::Element *Ai   = B + i*incRow + ii*incCol ;
						typename Field::Element *Aend = B + colDim*incCol ;
						for (; Ai != Aend + i*incRow ; Ai+=incCol)
							F.assign (*Ai, F.zero);
						size_t t = Q[ii-1];
						Q[ii-1]=Q[Nup+i];
						Q[Nup+i] = t;
					}
				}
			}
			return R + R2;
		}
	}



} // FFPACK

#undef LB_DEBUG
#endif //__FFLASFFPACK_ffpack_ludivine_INL

