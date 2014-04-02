/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* ffpack/ffpack_charpoly.inl
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
#ifndef __FFLASFFPACK_charpoly_INL
#define __FFLASFFPACK_charpoly_INL

namespace FFPACK {


	template <class Field, class Polynomial>
	std::list<Polynomial>&
	CharPoly (const Field& F, std::list<Polynomial>& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  const FFPACK_CHARPOLY_TAG CharpTag)
	{
		switch (CharpTag) {
		case FfpackLUK:
			{
				typename Field::Element * X = new typename Field::Element[N*(N+1)];
				Protected::LUKrylov (F, charp, N, A, lda, X, N);
				delete[] X;
				return charp;
			}
		case FfpackKG:
			{
				return Protected::KellerGehrig (F, charp, N, A, lda);
				// break;
			}
		case FfpackDanilevski:
			{
				return Danilevski (F, charp, N, A, lda);
				// break;
			}
		case FfpackKGFast:
			{
				size_t mc, mb, j;
				if (Protected::KGFast (F, charp, N, A, lda, &mc, &mb, &j)){
					std::cerr<<"NON GENERIC MATRIX PROVIDED TO KELLER-GEHRIG-FAST"<<std::endl;
				}
				return charp;
				// break;
			}
		case FfpackKGFastG:
			{
				return Protected::KGFast_generalized (F, charp, N, A, lda);
			}
		case FfpackHybrid:
			{
				typename Field::Element * X = new typename Field::Element[N*(N+1)];
				Protected::LUKrylov_KGFast (F, charp, N, A, lda, X, N);
				delete[] X;
				return charp;
			}
		case FfpackArithProg:
			{
				size_t attempts=0;
				bool cont = false;
				FFLAS_INT_TYPE p;
				F.characteristic(p);
				// Heuristic condition (the pessimistic theoretical one being p<2n^2.
				if ((unsigned long) (p) < N)
					return CharPoly (F, charp, N, A, lda, FfpackLUK);

				do{
					try {
						CharpolyArithProg (F, charp, N, A, lda, __FFPACK_CHARPOLY_THRESHOLD);
					}
					catch (CharpolyFailed){
						if (attempts++ < 2)
							cont = true;
						else
							return CharPoly(F, charp, N, A, lda, FfpackLUK);

					}
				} while (cont);
				return charp;
			}
		default:
			{
				typename Field::Element * X = new typename Field::Element[N*(N+1)];
				Protected::LUKrylov (F, charp, N, A, lda, X, N);
				delete[] X;
				return charp;
			}
		}
	}

	template<class Polynomial, class Field>
	Polynomial & mulpoly(const Field& F, Polynomial &res, const Polynomial & P1, const Polynomial & P2)
	{
		size_t i,j;
		// Warning: assumes that res is allocated to the size of the product
		res.resize(P1.size()+P2.size()-1);
		FFLAS::fzero(F,res.size(),&res[0],1);
		for ( i=0;i<P1.size();i++)
			for ( j=0;j<P2.size();j++)
				F.axpyin(res[i+j],P1[i],P2[j]);
		return res;
	}

	template <class Field, class Polynomial>
	Polynomial&
	CharPoly( const Field& F, Polynomial& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  const FFPACK_CHARPOLY_TAG CharpTag/*= FfpackArithProg*/)
	{

		std::list<Polynomial> factor_list;
		CharPoly (F, factor_list, N, A, lda, CharpTag);
		typename std::list<Polynomial >::const_iterator it;
		it = factor_list.begin();

		charp.resize(N+1);

		Polynomial P = charp = *(it++);

		while( it!=factor_list.end() ){
			mulpoly (F,charp, P, *it);
			P = charp;
			++it;
		}

		return charp;
	}


	namespace Protected {
		template <class Field, class Polynomial>
		std::list<Polynomial>&
		LUKrylov (const Field& F, std::list<Polynomial>& charp, const size_t N,
			  typename Field::Element * A, const size_t lda,
			  typename Field::Element * X, const size_t ldx)
		{

			typedef typename Field::Element elt;
			elt* Ai, *Xi, *X2=X;
			int Ncurr=int(N);
			charp.clear();
			int nbfac = 0;
			while (Ncurr > 0){
				size_t *P = new size_t[(size_t)Ncurr];
				Polynomial minP;//=new Polynomial();
				FFPACK::MinPoly (F, minP, (size_t)Ncurr, A, lda, X2, ldx, P);
				int k = int(minP.size()-1); // degre of minpoly
				if ((k==1) && F.isZero ((minP)[0])){ // minpoly is X
					Ai = A;
					int j = Ncurr*Ncurr;
					while (j-- && F.isZero(*(Ai++))) ;
					if (!j){ // A is 0, CharPoly=X^n
						minP.resize((size_t)Ncurr+1);
						(minP)[1] = F.zero;
						(minP)[(size_t)Ncurr] = F.one;
						k=Ncurr;
					}
				}
				nbfac++;
				charp.push_front (minP);
				if (k==Ncurr){
					delete[] P;
					return charp;
				}
				size_t Nrest = (size_t)(Ncurr-k);
				elt * X21 = X2 + k*(int)ldx;
				elt * X22 = X21 + k;
				// Compute the n-k last rows of A' = PA^tP^t in X2_
				// A = A . P^t
				applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					(size_t)Ncurr, 0, (int)k, A, lda, P);
				// Copy X2_ = (A'_2)^t
				for (Xi = X21, Ai = A+k; Xi != X21 + Nrest*ldx; Ai++, Xi+=ldx-(size_t)Ncurr)
					for (size_t jj=0; jj<(size_t)Ncurr*lda; jj+=lda)
						*(Xi++) = *(Ai+jj);
				// A = A . P : Undo the permutation on A
				applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					(size_t)Ncurr, 0, (int)k, A, lda, P);
				// X2_ = X2_ . P^t (=  (P A^t P^t)2_)
				applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					Nrest, 0, (int)k, X21, ldx, P);
				delete[] P ;
				// X21 = X21 . S1^-1
				ftrsm(F, FFLAS::FflasRight, FFLAS::FflasUpper,
				      FFLAS::FflasNoTrans, FFLAS::FflasUnit, Nrest, (size_t)k,
				      F.one, X2, ldx, X21, ldx);
				// Creation of the matrix A2 for recurise call
				for (Xi = X22, Ai = A;
				     Xi != X22 + Nrest*ldx;
				     Xi += (ldx-Nrest), Ai += (lda-Nrest))
					for (size_t jj=0; jj<Nrest; ++jj)
						*(Ai++) = *(Xi++);
				fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Nrest, Nrest, (size_t)k, F.mOne,
				       X21, ldx, X2+k, ldx, F.one, A, lda);
				X2 = X22;
				Ncurr = int(Nrest);
			}
			return charp;
		}


		template <class Field, class Polynomial>
		std::list<Polynomial>&
		LUKrylov_KGFast (const Field& F, std::list<Polynomial>& charp, const size_t N,
				 typename Field::Element * A, const size_t lda,
				 typename Field::Element * X, const size_t ldx)
		{

			typedef typename Field::Element elt;

			size_t kg_mc, kg_mb, kg_j;

			if (!KGFast (F, charp, N, A, lda, &kg_mc, &kg_mb, &kg_j))
				return charp;
			else{// Matrix A is not generic
				Polynomial *minP = new Polynomial();
				const elt* Ai;
				elt* A2i, *Xi;
				size_t *P = new size_t[N];

				FFPACK::MinPoly (F, *minP, N, A, lda, X, ldx, P, FfpackKGF, kg_mc, kg_mb, kg_j);
				size_t k = minP->size()-1; // degre of minpoly
				if ((k==1) && F.isZero ((*minP)[0])){ // minpoly is X
					Ai = A;
					int j = int(N*N);
					while (j-- && F.isZero(*(Ai++))) ;
					if (!j){ // A is 0, CharPoly=X^n
						minP->resize((size_t)N+1);
						(*minP)[1] = F.zero;
						(*minP)[N] = F.one;
						k=N;
					}
				}

				if (k==N){
					charp.clear();
					charp.push_front(*minP); // CharPoly = MinPoly
					delete[] P;
					return charp;
				}

				size_t Nrest = N-k;
				elt * X21 = X + k*ldx;
				elt * X22 = X21 + k;

				// Creates the matrix A
				//size_t lambda = std::max(0,N - kg_mc*(kg_j+1) - kg_mb);  // uint >= 0 !!!
				size_t lambda =   kg_mc*(kg_j+1) + kg_mb;
				if (lambda > N)
					lambda = 0 ;
				else
					lambda = N - lambda ;

				size_t imax = kg_mc+kg_mb;
				// First Id
				for (size_t j = 0; j < lambda; ++j){
					for (size_t i=0; i<imax; ++i)
						F.assign (*(A+j+i*lda), F.zero);
					F.assign (*(A+j+imax*lda), F.one);
					for (size_t i=imax+1; i<N; ++i)
						F.assign (*(A+j+i*lda), F.zero);
					++imax;
				}
				// Column block B
				for (typename Field::Element* Aj=A; Aj<A+N*lda; Aj+=lda)
					FFLAS::fcopy (F, kg_mb, Aj+lambda, 1, Aj+N-kg_mc-kg_mb, 1);

				// Second Id block
				imax = N- kg_j*kg_mc;
				for (size_t j = 0; j< kg_j*kg_mc; ++j){
					for (size_t i = 0; i<imax; ++i)
						F.assign (*(A+lambda+kg_mb+j+i*lda), F.zero);
					F.assign (*(A+lambda+kg_mb+j+imax*lda), F.one);
					for (size_t i = imax+1; i<N; ++i)
						F.assign (*(A+lambda+kg_mb+j+i*lda), F.zero);
					++imax;
				}

				// Compute the n-k last rows of A' = PA^tP^t in X2_

				// A = P . A
				applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
					N, 0,(int) k,
					const_cast<typename Field::Element* &>(A), lda, P);

				// Copy X2_ = (A'2_)
				for (Xi = X21, Ai = A+k*lda; Xi != X21 + Nrest*ldx; Ai+=lda-N, Xi+=ldx-N){
					for (size_t jj=0; jj<N; ++jj){
						*(Xi++) = *(Ai++);
					}
				}

				// A = P^t . A : Undo the permutation on A
				applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
					N, 0,(int) k,
					const_cast<typename Field::Element* &>(A), lda, P);

				// X2_ = X2_ . P^t (=  (P A P^t)2_)
				applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					Nrest, 0,(int) k, X21, ldx, P);

				// X21 = X21 . S1^-1
				ftrsm(F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Nrest, k,
				      F.one, X, ldx, X21, ldx);

				// Creation of the matrix A2 for recurise call
				elt * A2 = new elt[Nrest*Nrest];

				for (Xi = X22, A2i = A2;
				     Xi != X22 + Nrest*ldx;
				     Xi += (ldx-Nrest)){
					for (size_t jj=0; jj<Nrest; ++jj){
						*(A2i++) = *(Xi++);
					}
				}
				fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Nrest, Nrest, k, F.mOne,
				       X21, ldx, X+k, ldx, F.one, A2, Nrest);

				// Recursive call on X22
				LUKrylov_KGFast (F, charp, Nrest, A2, Nrest, X22, ldx);
				charp.push_front (*minP);
				delete[] P;
				delete[] A2;
				return charp;
			}
		}

	} // Protected

} // FFPACK

#endif // __FFLASFFPACK_charpoly_INL
