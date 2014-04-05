/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas-ffpack/ffpack/ffpack_frobenius.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <cpernet@uwaterloo.ca>
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

#include <fflas-ffpack/field/nonzero-randiter.h>

//---------------------------------------------------------------------
// CharpolyArithProg: Las Vegas algorithm to compute the Charpoly
// over a large field (Z/pZ, s.t.  p > 2n^2)
//---------------------------------------------------------------------
//
//
namespace FFPACK { namespace Protected {
	template <class Field>
	void CompressRows (Field& F, const size_t M,
			   typename Field::Element * A, const size_t lda,
			   typename Field::Element * tmp, const size_t ldtmp,
			   const size_t * d, const size_t nb_blocs);

	template <class Field>
	void CompressRowsQK (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d,const size_t deg, const size_t nb_blocs);

	template <class Field>
	void DeCompressRows (Field& F, const size_t M, const size_t N,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs);
	template <class Field>
	void DeCompressRowsQK (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t deg, const size_t nb_blocs);

	template <class Field>
	void CompressRowsQA (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs);
	template <class Field>
	void DeCompressRowsQA (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t nb_blocs);
	} // Protected
} // FFPACK


template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::CharpolyArithProg (const Field& F, std::list<Polynomial>& frobeniusForm,
			   const size_t N, typename Field::Element * A, const size_t lda,
			   const size_t c)
{

	FFLASFFPACK_check(c);

	size_t * rp = new size_t[2*N];
	size_t noc = static_cast<size_t>(ceil(double(N)/double(c)));
	size_t Nnoc = N*noc;

	// Building the workplace matrix
	typename Field::Element *K  = new typename Field::Element[Nnoc*c];
	typename Field::Element *K2 = new typename Field::Element[Nnoc*c];
	// for (size_t i = 0 ; i < Nnoc*c ; ++i)
		// K[i] = F.zero;
	size_t ldk = N;

	size_t *dA = new size_t[N]; //PA
	size_t *dK = new size_t[noc*c];
	for (size_t i=0; i<noc; ++i)
		dK[i]=0;

	// Picking a random noc x N block vector U^T
	typename Field::RandIter g (F);
	NonzeroRandIter<Field> nzg (F,g);
	for (size_t i = 0; i < noc; ++i)
 		for (size_t j = 0; j < N; ++j)
 			g.random( *(K + i*ldk +j) );
	for (size_t i = 0; i < noc; ++i)
		nzg.random (*(K + i*ldk +i));

	// Computing the bloc Krylov matrix [U AU .. A^(c-1) U]^T
	for (size_t i = 1; i<c; ++i){
// #warning "leaks here"
		fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasTrans,  noc, N, N,F.one,
		       K+(i-1)*Nnoc, ldk, A, lda, F.zero, K+i*Nnoc, ldk);
	}
	// K2 <- K (re-ordering)
	//! @todo swap to save space ??
	size_t w_idx = 0;
	for (size_t i=0; i<noc; ++i)
		for (size_t j=0; j<c; ++j, w_idx++)
			FFLAS::fcopy(F, N, (K+(i+j*noc)*ldk), 1, (K2+(w_idx)*ldk), 1);

	// Copying K <- K2
	for (size_t i=0; i<noc*c; ++i)
		FFLAS::fcopy (F, N, K2+i*ldk, 1, (K+i*ldk), 1);

	size_t * Pk = new size_t[N];
	size_t * Qk = new size_t[N];
	for (size_t i=0; i<N; ++i)
		Qk[i] = 0;
	for (size_t i=0; i<N; ++i)
		Pk[i] = 0;

	size_t R = LUdivine(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, N, N, K, ldk, Pk, Qk, FfpackLQUP);

	size_t row_idx = 0;
	size_t ii=0;
	size_t dold = c;
	size_t nb_full_blocks = 0;
	size_t Mk = 0;
	// Determining the degree sequence dK
	for (size_t k = 0; k<noc; ++k){
		size_t d = 0;
		while ( (d<c) && (row_idx<R) && (Qk[row_idx] == ii)) {ii++; row_idx++; d++;}
		if (d > dold){
			// std::cerr << "FAIL in preconditionning phase:"
			//           << " degree sequence is not monotonically not increasing"
			// 	     << std::endl;
			delete[] rp; delete[] K;
			delete[] Pk; delete[] Qk; delete[] dA; delete[] dK;
			throw CharpolyFailed();
		}
		dK[k] = dold = d;
		Mk++;
		if (d == c)
			nb_full_blocks++;
		if (row_idx < N)
			ii = Qk[row_idx];
	}

	// Selection of the last iterate of each block

	typename Field::Element * K3 = new typename Field::Element[Mk*N];
	typename Field::Element * K4 = new typename Field::Element[Mk*N];
	size_t bk_idx = 0;
	for (size_t i = 0; i < Mk; ++i){
		FFLAS::fcopy (F, N, (K2 + (bk_idx + dK[i]-1)*ldk), 1, (K3+i*ldk), 1);
		bk_idx += c;
	}
	delete[] K2;

	// K <- K A^T
	fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasTrans, Mk, N, N,F.one,  K3, ldk, A, lda, F.zero, K4, ldk);

	// K <- K P^T
	applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
		Mk, 0,(int) R, K4, ldk, Pk);

	// K <- K U^-1
	ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, Mk, R,F.one, K, ldk, K4, ldk);

	// L <-  Q^T L
	applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
	       N, 0,(int) R, K, ldk, Qk);

	// K <- K L^-1
	ftrsm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Mk, R,F.one, K, ldk, K4, ldk);

	//undoing permutation on L
	applyP(F, FFLAS::FflasLeft, FFLAS::FflasTrans,
	       N, 0,(int) R, K, ldk, Qk);

	// Recovery of the completed invariant factors
	size_t Ma = Mk;
	size_t Ncurr = R;
	size_t offset = Ncurr-1;
	for (size_t i=Mk-1; i>=nb_full_blocks+1;  --i){
		if (dK[i] >= 1){
			for (size_t j = offset+1; j<R; ++j)
				if (!F.isZero(*(K4 + i*ldk + j))){
					//std::cerr<<"FAIL C != 0 in preconditionning"<<std::endl;
					delete[] K3; delete[] K4; delete[] K;
					delete[] Pk; delete[] Qk; delete[] rp;
					delete[] dA; delete[] dK;
					throw CharpolyFailed();
				}
			Polynomial P (dK [i]+1);
			F.assign(P[dK[i]],F.one);
			for (size_t j=0; j < dK [i]; ++j)
				F.neg (P [dK [i]-j-1], *(K4 + i*ldk + (offset-j)));
			frobeniusForm.push_front(P);
			offset -= dK [i];
			Ncurr -= dK [i];
			Ma--;
		}
	}
	Mk = Ma;

	if (R<N){
		for (size_t i=0; i<nb_full_blocks + 1; ++i)
			for (size_t j=R; j<N; ++j){
				if (!F.isZero( *(K4+i*ldk+j) )){
					delete[] K3; delete[] K4; delete[] K;
					delete[] Pk; delete[] Qk; delete[] rp;
					delete[] dA; delete[] dK;
					throw CharpolyFailed();
				}
			}

		//std::cerr<<"Preconditionning failed; missing rank = "<<N-R
		//	 <<" completing the Krylov matrix"
		//	 <<std::endl;
		size_t Nrest = N-R;
		typename Field::Element * K21 = K + R*ldk;
		typename Field::Element * K22 = K21 + R;
		typename Field::Element * Ki, *Ai;

		//  Compute the n-k last rows of A' = P A^T P^T in K2_
		// A = A . P^t
		applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
			N, 0,(int) R, A, lda, Pk);

		// Copy K2_ = (A'_2)^t
		for (Ki = K21, Ai = A+R; Ki != K21 + Nrest*ldk; Ai++, Ki+=ldk-N)
			for ( size_t j=0; j<N*lda; j+=lda )
				*(Ki++) = *(Ai+j);

		// A = A . P : Undo the permutation on A
		applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
			N, 0,(int) R, A, lda, Pk);

		// K2_ = K2_ . P^t (=  ( P A^t P^t )2_ )
		applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
			Nrest, 0,(int) R, K21, ldk, Pk);

		// K21 = K21 . S1^-1
		ftrsm (F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, Nrest, R,
		      F.one, K, ldk, K21, ldk);

		typename Field::Element * Arec = new typename Field::Element[Nrest*Nrest];
		size_t ldarec = Nrest;

		// Creation of the matrix A2 for recursive call
		for (Ki = K22,  Ai = Arec;
		     Ki != K22 + Nrest*ldk;
		     Ki += (ldk-Nrest) )
			for ( size_t j=0; j<Nrest; ++j )
				*(Ai++) = *(Ki++);
		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Nrest, Nrest, R,F.mOne,
		       K21, ldk, K+R, ldk,F.one, Arec, ldarec);

		std::list<Polynomial> polyList;
		polyList.clear();

		// Recursive call on the complementary subspace
		CharPoly(F, polyList, Nrest, Arec, ldarec);
		delete[] Arec;
		frobeniusForm.merge(polyList);
	}

	delete[] Pk;
	delete[] Qk;
	size_t deg = c+1;
	for (size_t i=0; i<Mk; ++i)
 		dA[i] = dK[i];
	bk_idx = 0;

	typename Field::Element *Arp = new typename Field::Element[Ncurr*Ma];
	typename Field::Element *Ac = new typename Field::Element[Ncurr*Ma];
	size_t ldac = Ma;
	size_t ldarp = Ncurr;

	for (size_t i=0; i < Ncurr; ++i)
 		for (size_t j=0; j<Ma; ++j)
			*(K+i*ldk+j) = *(Ac + i*Ma +j) = *(K4 + i + (j)*ldk);
	delete[] K4;


	// Main loop of the arithmetic progession
	while ((nb_full_blocks >= 1) && (Mk > 1)) {
		size_t block_idx, it_idx, rp_val;
		delete[] K;
		delete[] K3;
		K = new typename Field::Element[Ncurr*Ma];
		K3 = new typename Field::Element[Ncurr*Ma];
		ldk = Ma;

		// Computation of the rank profile
		for (size_t i=0; i < Ncurr; ++i)
			for (size_t j=0; j < Ma; ++j)
				*(Arp + j*ldarp + Ncurr-i-1) = *(Ac + i*ldac + j);
		for (size_t i=0; i<2*Ncurr; ++i)
			rp[i] = 0;
		size_t RR;
		try{
			RR = SpecRankProfile (F, Ma, Ncurr, Arp, ldarp, deg-1, rp);
		} catch (CharpolyFailed){
			delete[] Arp; delete[] Ac; delete[] K; delete[] K3;
			delete[] rp; delete[] dA; delete[] dK;
			throw CharpolyFailed();
		}
		if (RR < Ncurr){
			//std::cerr<<"FAIL RR<Ncurr"<<std::endl;
			delete[] Arp; delete[] Ac; delete[] K; delete[] K3;
			delete[] rp; delete[] dA; delete[] dK;
			throw CharpolyFailed();
		}

		// Computation of the degree vector dK
		it_idx = 0;
		rp_val = 0;
		size_t gg = 0;
		size_t dtot=0;
		block_idx = 0;
		nb_full_blocks = 0;
		while (dtot<Ncurr){
			do {gg++; rp_val++; it_idx++;}
			while ( /*(gg<Ncurr ) &&*/ (rp[gg] == rp_val) && (it_idx < deg ));
			if ((block_idx)&&(it_idx > dK[block_idx-1])){
				delete[] Arp; delete[] Ac;delete[] K; delete[] K3;
				delete[] rp; delete[] dA; delete[] dK;
				throw CharpolyFailed();
				//std::cerr<<"FAIL d non decroissant"<<std::endl;
				//exit(-1);
			}
			dK[block_idx++] = it_idx;
			dtot += it_idx;
			if (it_idx == deg)
				nb_full_blocks ++;
			it_idx=0;
			rp_val = rp[gg];
		}

		Mk = block_idx;

		// Selection of dense colums of K
		for (size_t i=0; i < nb_full_blocks; ++i){
			FFLAS::fcopy (F, Ncurr, Ac+i, ldac, K+i, ldk);
		}

		// K <- QK K
		size_t pos = nb_full_blocks*(deg-1);
		for (size_t i = nb_full_blocks; i < Mk; ++i){
			for (size_t j=0; j<Ncurr; ++j)
				F.assign (*(K + i + j*ldk), F.zero);
			F.assign (*(K + i + (pos + dK[i]-1)*ldk),F.one);
			pos += dA[i];
		}

		// Copying K3 <- K
		for (size_t i=0; i<Mk; ++i)
			FFLAS::fcopy (F, Ncurr, K+i, ldk, K3+i, ldk);
		Protected::CompressRowsQK (F, Mk, K3 + nb_full_blocks*(deg-1)*ldk, ldk,
				Arp, ldarp, dK+nb_full_blocks, deg, Mk-nb_full_blocks);

		// K <- PA K
		Protected::CompressRows (F, nb_full_blocks, K, ldk, Arp, ldarp, dA, Ma);

		// A <- newQA^T K (compress)
		Protected::CompressRowsQA (F, Ma, Ac, ldac, Arp, ldarp, dA, Ma);

		// K <- A K
		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ncurr-Ma, nb_full_blocks, Ma,F.one,
		       Ac, ldac, K+(Ncurr-Ma)*ldk, ldk,F.one, K, ldk);
		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ma, nb_full_blocks, Ma,F.one,
		       Ac+(Ncurr-Ma)*ldac, ldac, K+(Ncurr-Ma)*ldk, ldk, F.zero, Arp, ldarp);
		for (size_t i=0; i< Ma; ++i)
			FFLAS::fcopy(F, nb_full_blocks, Arp+i*ldarp, 1, K+(Ncurr-Ma+i)*ldk, 1);

		// Copying the last rows of A times K
		offset = (deg-2)*nb_full_blocks;
		for (size_t i = nb_full_blocks; i < Mk; ++i) {
			for (size_t j=0; j<Ncurr; ++j)
				F.assign(*(K+i+j*ldk), F.zero);
			if (dK[i] == dA[i]) // copy the column of A
				FFLAS::fcopy (F, Ncurr, Ac+i, ldac, K+i, ldk);
			else{
				F.assign (*(K + i + (offset+dK[i]-1)*ldk),F.one);
			}
			offset += dA[i]-1;
		}

		// K <- QA K
		Protected::DeCompressRowsQA (F, Mk, Ncurr, K, ldk, Arp, ldarp, dA, Ma);

		// K <- QK^T K
		Protected::CompressRowsQK (F, Mk, K + nb_full_blocks*(deg-1)*ldk, ldk, Arp, ldarp,
				dK+nb_full_blocks, deg, Mk-nb_full_blocks);

		// K <- K^-1 K
		size_t *P=new size_t[Mk];
		size_t *Q=new size_t[Mk];
		if (LUdivine (F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, Mk, Mk , K3 + (Ncurr-Mk)*ldk, ldk, P, Q, FfpackLQUP) < Mk){
			// should never happen (not a LAS VEGAS check)
			//std::cerr<<"FAIL R2 < MK"<<std::endl;
			//			exit(-1);
		}
		ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, Mk, Mk,F.one,
		       K3 + (Ncurr-Mk)*ldk, ldk, K+(Ncurr-Mk)*ldk, ldk);
		ftrsm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, Mk, Mk,F.one,
		       K3+(Ncurr-Mk)*ldk, ldk, K+(Ncurr-Mk)*ldk, ldk);
		applyP (F, FFLAS::FflasLeft, FFLAS::FflasTrans,
			Mk, 0,(int) Mk, K+(Ncurr-Mk)*ldk,ldk, P);
		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, Ncurr-Mk, Mk, Mk,F.mOne,
		       K3, ldk, K+(Ncurr-Mk)*ldk,ldk,F.one, K, ldk);
		delete[] P;
		delete[] Q;

		// K <- PK^T K
		Protected::DeCompressRows (F, Mk, Ncurr, K, ldk, Arp, ldarp, dK, Mk);

		// K <- K PK (dA <- dK)
		if (nb_full_blocks*deg < Ncurr)
			Ma = nb_full_blocks+1;
		else
			Ma = nb_full_blocks;

		for (size_t i=0; i< Ma; ++i)
			dA[i] = dK[i];

		// Recovery of the completed invariant factors
		offset = Ncurr-1;
		size_t oldNcurr = Ncurr;
		for (size_t i=Mk-1; i>=nb_full_blocks+1;  --i)
			if (dK[i] >= 1){
				Polynomial  PP (dK [i]+1);
				F.assign(PP[dK[i]],F.one);
				for (size_t j=0; j < dK[i]; ++j)
					F.neg( PP[dK[i]-j-1], *(K + i + (offset-j)*ldk));
				frobeniusForm.push_front(PP);
				offset -= dK[i];
				Ncurr -= dK[i];
			}
		for (size_t i= offset+1; i<oldNcurr; ++i)
			for (size_t j=0; j<nb_full_blocks+1; ++j){
				if (!F.isZero( *(K+i*ldk+j) )){
					//std::cerr<<"FAIL C != 0"<<std::endl;
					delete[] rp; delete[] Arp; delete[] Ac;
					delete[] K; delete[] K3;
					delete[] dA; delete[] dK;
					throw CharpolyFailed();
				}
			}

		// A <- K
		delete[] Ac; delete[] Arp;
		Ac = new typename Field::Element[Ncurr*Mk];
		ldac = Mk;
		Arp = new typename Field::Element[Ncurr*Mk];
		ldarp=Ncurr;
		for (size_t i=0; i < Ncurr; ++i )
			FFLAS::fcopy (F, Mk, K + i*ldk, 1, Ac + i*ldac, 1);

		deg++;

	}

	// Recovery of the first invariant factor
	Polynomial Pl(dK [0]+1);
	F.assign(Pl[dK[0]],F.one);
	for (size_t j=0; j < dK[0]; ++j)
		F.neg( Pl[j], *(K  + j*ldk));
	frobeniusForm.push_front(Pl);
	delete[] rp; delete[] Arp; delete[] Ac; delete[] K; delete[] K3;
	delete[] dA; delete[] dK;
	return frobeniusForm;
}

namespace FFPACK { namespace Protected {
template <class Field>
void CompressRowsQK (Field& F, const size_t M,
			   typename Field::Element * A, const size_t lda,
			   typename Field::Element * tmp, const size_t ldtmp,
			   const size_t * d, const size_t deg,const size_t nb_blocs)
{

	int currtmp = 0;
	size_t currw = d[0]-1;
	size_t currr = d[0]-1;
	for (int i = 0; i< int(nb_blocs)-1; ++i){
		// FFLAS::fcopy(F,deg-d[i],M,A+currr*lda,lda,tmp+(size_t)currtmp*ldtmp);
		for (int j = int(d[i]-1); j<int(deg)-1; ++j, ++currr, ++currtmp)
			FFLAS::fcopy(F, M,  A + currr*lda, 1, tmp + (size_t)currtmp*ldtmp, 1);
		// currr += (deg - d[i]);
		for (int j=0; j < int(d[i+1]) -1; ++j, ++currr, ++currw){
			FFLAS::fcopy(F, M, A+(currr)*lda, 1, A + (currw)*lda, 1);
		}
	}
	for (int i=0; i < currtmp; ++i, ++currw){
		FFLAS::fcopy (F, M, tmp + (size_t)i*ldtmp, 1, A + (currw)*lda, 1);
	}
}

template <class Field>
void CompressRows (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs)
{

	size_t currd = d[0]-1;
	size_t curri = d[0]-1;
	for (int i = 0; i< int(nb_blocs)-1; ++i){
		FFLAS::fcopy(F, M,  A + currd*lda, 1, tmp + i*(int)ldtmp, 1);
		for (int j=0; j < int(d[i+1]) -1; ++j){
			FFLAS::fcopy(F, M, A+(currd+(size_t)j+1)*lda, 1, A + (curri++)*lda, 1);
		}
		currd += d[i+1];
	}
	for (int i=0; i < int(nb_blocs)-1; ++i){
		FFLAS::fcopy (F, M, tmp + i*(int)ldtmp, 1, A + (curri++)*lda, 1);
	}
}

template <class Field>
void DeCompressRows (Field& F, const size_t M, const size_t N,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs)
{

	for (int i=0; i<int(nb_blocs)-1; ++i)
		FFLAS::fcopy(F, M, A + (N-nb_blocs+(size_t)i)*lda, 1, tmp + i*(int)ldtmp, 1);

	size_t w_idx = N - 2;
	size_t r_idx = N - nb_blocs - 1;
	int i = int(nb_blocs)-1 ;
	for (; i--; ){
		for (size_t j = 0; j<d[i+1]-1; ++j)
			FFLAS::fcopy (F, M, A + (r_idx--)*lda, 1, A + (w_idx--)*lda, 1);
		FFLAS::fcopy (F, M, tmp + i*(int)ldtmp, 1, A + (w_idx--)*lda, 1);
	}
}

template <class Field>
void DeCompressRowsQK (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t deg,const size_t nb_blocs)
{

	size_t zeroblockdim = 1; // the last block contributes with 1
	size_t currtmp = 0;
	for (int i=0; i<int(nb_blocs)-1; ++i)
		zeroblockdim += deg - d[i];
	for (size_t i=0; i < zeroblockdim - 1; ++i, ++currtmp)
		FFLAS::fcopy(F, M,  A + (N - zeroblockdim +i)*lda, 1, tmp + currtmp*ldtmp, 1);
	currtmp--;
	size_t w_idx = N - 2;
	size_t r_idx = N - zeroblockdim - 1;

	int i = int(nb_blocs)-1 ;
	for (; i--;){
		for (size_t j = 0; j < d [i+1] - 1; ++j)
			FFLAS::fcopy (F, M, A + (r_idx--)*lda, 1, A + (w_idx--)*lda, 1);
		for (size_t j = 0; j < deg - d[i]; ++j)
			FFLAS::fcopy (F, M, tmp + (currtmp--)*ldtmp, 1, A + (w_idx--)*lda, 1);
	}
}

template <class Field>
void CompressRowsQA (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs)
{

	size_t currd = 0;
	size_t curri = 0;
	for (size_t i = 0; i< nb_blocs; ++i){
		FFLAS::fcopy(F, M,  A + currd*lda, 1, tmp + i*ldtmp, 1);
		for (size_t j=0; j < d[i] -1; ++j)
			FFLAS::fcopy(F, M, A+(currd+j+1)*lda, 1, A + (curri++)*lda, 1);
		currd += d[i];
	}
	for (size_t i=0; i < nb_blocs; ++i)
		FFLAS::fcopy (F, M, tmp + i*ldtmp, 1, A + (curri++)*lda, 1);
}

template <class Field>
void DeCompressRowsQA (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t nb_blocs)
{

	for (size_t i=0; i<nb_blocs; ++i)
		FFLAS::fcopy(F, M, A + (N-nb_blocs+i)*lda, 1, tmp + i*ldtmp, 1);

	size_t w_idx = N - 1;
	size_t r_idx = N - nb_blocs - 1;
	int i = int(nb_blocs) ;
	for (; i--; ){
		for (size_t j = 0; j<d[i]-1; ++j)
			FFLAS::fcopy (F, M, A + (r_idx--)*lda, 1, A + (w_idx--)*lda, 1);
		FFLAS::fcopy (F, M, tmp + i*(int)ldtmp, 1, A + (w_idx--)*lda, 1);
	}
}

} // Protected
} //FFPACK
