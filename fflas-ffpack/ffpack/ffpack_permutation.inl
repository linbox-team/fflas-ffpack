/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 FFLAS-FFACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * BB <bbboyer@ncsu.edu>
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

#ifndef __FFLASFFPACK_ffpack_permutation_INL
#define __FFLASFFPACK_ffpack_permutation_INL


namespace FFPACK {

	template<class Field>
	void
	applyP( const Field& F,
		const FFLAS::FFLAS_SIDE Side,
		const FFLAS::FFLAS_TRANSPOSE Trans,
		const size_t M, const int ibeg, const int iend,
		typename Field::Element * A, const size_t lda, const size_t * P )
	{

		if ( Side == FFLAS::FflasRight ) {
			typename Field::Element tmp;
			if ( Trans == FFLAS::FflasTrans )
				for (size_t j = 0 ; j < M ; ++j){
					for ( size_t i=(size_t)ibeg; i<(size_t) iend; ++i)
						if ( P[i]!= i ) {
							F.assign(tmp,A[j*lda+P[i]]);
							F.assign(A[j*lda+P[i]],A[j*lda+i]);
							F.assign(A[j*lda+i],tmp);
							// std::swap(A[j*lda+P[i]],A[j*lda+i]);
						}
					//FFLAS::fswap( F, M, A + P[i]*1, lda, A + i*1, lda );
				}
			else // Trans == FFLAS::FflasNoTrans
				for (size_t j = 0 ; j < M ; ++j){
					for (int i=iend; i-->ibeg; )
						if ( P[i]!=(size_t)i ) {
							F.assign(tmp,A[j*lda+P[i]]);
							F.assign(A[j*lda+P[i]],A[j*lda+(size_t)i]);
							F.assign(A[j*lda+(size_t)i],tmp);
							// std::swap(A[j*lda+P[i]],A[j*lda+(size_t)i]);
						}
					//FFLAS::fswap( F, M, A + P[i]*1, lda, A + i*1, lda );
				}
		}
		else { // Side == FFLAS::FflasLeft
			if ( Trans == FFLAS::FflasNoTrans )
				for (size_t i=(size_t)ibeg; i<(size_t)iend; ++i){
					if ( P[i]!= (size_t) i )
						FFLAS::fswap( F, M,
							      A + P[i]*lda, 1,
							      A + i*lda, 1 );
				}
			else // Trans == FFLAS::FflasTrans
				for (int i=iend; i-->ibeg; ){
					if ( P[i]!= (size_t) i ){
						FFLAS::fswap( F, M,
							      A + P[i]*lda, 1,
							      A + (size_t)i*lda, 1 );
					}
				}
		}

	}

	template <class Element>
	inline void applyS (Element* A, const size_t lda, const size_t width,
			    const size_t M2,
			    const size_t R1, const size_t R2,
			    const size_t R3, const size_t R4)
	{
		Element * tmp = new Element [(M2-R1-R2)*width];
		// std::cerr<<"ici"<<std::endl;
		for (size_t i = 0, j = R1+R2; j < M2; ++i, ++j)
			for (size_t k = 0; k<width; ++k)
				tmp [i*width + k] = A [j*lda + k];
		for (size_t i = M2, j = R1+R2; i < M2+R3+R4; ++i, ++j)
			for (size_t k = 0; k<width; ++k)
				A [j*lda + k] = A [i*lda +k];
		for (size_t i = 0, j = R1+R2+R3+R4; i < M2-R1-R2; ++i, ++j)
			for (size_t k = 0; k<width; ++k)
				A [j*lda + k] = tmp [i*width + k];
		delete[] tmp;
	}


	template <class Element>
	inline void applyT (Element* A, const size_t lda, const size_t width,
			    const size_t N2,
			    const size_t R1, const size_t R2,
			    const size_t R3, const size_t R4)
	{
		Element * tmp = new Element[(N2-R1)*width];
		for (size_t k = 0; k < width; ++k){
			for (size_t i = 0, j = R1; j < N2; ++i, ++j){
				tmp [i + k*(N2-R1)] = A [k*lda + j];
			}

			for (size_t i = N2, j = R1; i < N2+R2; ++i, ++j)
				A [k*lda + j] = A [k*lda + i];

			for (size_t i = 0, j = R1+R2; i < R3; ++i, ++j)
				A [k*lda + j] = tmp [k*(N2-R1) + i];

			for (size_t i = N2+R2, j = R1+R2+R3; i < N2+R2+R4; ++i,++j)
				A [k*lda + j] = A [k*lda + i];
			for (size_t i = R3, j = R1+R2+R3+R4; i < N2-R1; ++i,++j)
				A [k*lda + j] = tmp [k*(N2-R1) + i];
		}
		delete[] tmp;
	}


            /**
	     * Conversion of a permutation from LAPACK format to Math format
	     */
	inline void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP,
				  const size_t N)
	{
		for (size_t i=0; i<N; i++)
			MathP[i] = i;
		for (size_t i=0; i<N; i++){
			if (LapackP[i] != i){
				size_t tmp = MathP[i];
				MathP[i] = MathP[LapackP[i]];
				MathP[LapackP[i]] = tmp;
			}
		}
	}

	    /**
	     * Conversion of a permutation from Maths format to LAPACK format
	     */
	inline void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP,
					 const size_t N)
	{
		size_t * T = new size_t[N];
		size_t * Tinv = new size_t[N];
		for (size_t i=0; i<N; i++){
			T[i] =i;
			Tinv[i] = i;
		}
		for (size_t i=0; i<N; i++){
			size_t j = Tinv [MathP [i]];
			LapackP [i] = j;
			size_t tmp = T[j];
			T[j] = T[i];
			Tinv[T[i]] = j;
			T[i] = tmp;
			Tinv[tmp] = i;
		}
		delete[] T;
		delete[] Tinv;
	}

	    /**
	     * Computes P1 [ I_R     ] stored in MathPermutation format
	     *             [     P_2 ]
	     */
	inline void composePermutationsP (size_t * MathP,
				  const size_t * P1,
				  const size_t * P2,
				  const size_t R, const size_t N)
	{
		for (size_t i=0; i<N; ++i)
			MathP[i] = i;
		LAPACKPerm2MathPerm (MathP, P1, N);

		for (size_t i=R; i<N; i++){
			if (P2[i-R] != i-R){
				size_t tmp = MathP[i];
				MathP[i] = MathP[P2[i-R]+R];
				MathP[P2[i-R]+R] = tmp;
			}
		}
	}

	inline void composePermutationsQ (size_t * MathP,
				  const size_t * Q1,
				  const size_t * Q2,
				  const size_t R, const size_t N)
	{
		for (size_t i=0; i<N; ++i)
			MathP[i] = i;
		LAPACKPerm2MathPerm (MathP, Q1, N);

		for (size_t i=R; i<N; i++){
			if (Q2[i-R] != i-R){
				size_t tmp = MathP[i];
				MathP[i] = MathP[Q2[i-R]+R];
				MathP[Q2[i-R]+R] = tmp;
			}
		}
	}

	inline void
	cyclic_shift_mathPerm (size_t * P,  const size_t s)
	{
                size_t tmp;
                tmp = *(P+s-1);
		    //memmove(P+1, P, (s)*sizeof(size_t));
		size_t * Pi = P;
		std::copy(Pi, Pi+s-1, Pi+1);

                *(P)=tmp;
	}

	template<typename Base_t>
	inline void cyclic_shift_row_col(Base_t * A, size_t m, size_t n, size_t lda)
	{

#ifdef MEMCOPY
//		std::cerr << "BEF m: " << m << ", n: " << n << std::endl;

		if (m > 1) {
			const size_t mun(m-1);
			if (n > 1) {
//     std::cerr << "m: " << m << ", n: " << n << std::endl;
				const size_t nun(n-1);
				const size_t blo(sizeof(Base_t));
				// const size_t bmu(blo*mun);
				const size_t bnu(blo*nun);
				Base_t * b = new Base_t[mun];
				for(size_t i=0; i<mun; ++i) b[i] = A[i*lda+nun];
				Base_t * dc = new Base_t[n];
				memcpy(dc+1,A+mun*lda,bnu);
				*dc = A[mun*lda+nun]; // this is d
				    // dc = [ d c ]

				for(size_t i=mun; i>0; --i)
					memcpy(A+1+i*lda, A+(i-1)*lda, bnu);

				memcpy(A, dc, bnu+blo);
				for(size_t i=0; i<mun; ++i) A[i*lda+lda] = b[i];
				delete [] dc;
				delete [] b;

			} else if (n != 0) {
				Base_t d = A[mun*lda];
				for(size_t i=mun; i>0; --i) A[i*lda]=A[(i-1)*lda];
				*A=d;
			}
		} else {
			if ((m!=0) && (n > 1)) {
				const size_t nun(n-1);
				const size_t blo(sizeof(Base_t));
				const size_t bnu(blo*nun);
				Base_t d = A[nun];
//  std::cerr << "d: " << d << std::endl;
				Base_t * tmp = new Base_t[nun];
				memcpy(tmp,A,bnu);
				memcpy(A+1,tmp,bnu);
//				std::copy(A,A+nun,A+1);
				*A=d;
				delete [] tmp;

			}
		}
//		std::cerr << "AFT m: " << m << ", n: " << n << std::endl;

#else

		    //	std::cerr << "BEF m: " << m << ", n: " << n << std::endl;
		if (m > 1) {
			const size_t mun(m-1);
			if (n > 1) {
				const size_t nun(n-1);

				Base_t * b = new Base_t[mun];
				Base_t * Ainun = A+nun;
				for(size_t i=0; i<mun; ++i, Ainun+=lda) b[i] = *Ainun;

				    // dc = [ d c ]
				Base_t * dc = new Base_t[n];
				std::copy(Ainun-nun, Ainun, dc+1);

				    // this is d
				*dc = *Ainun;

				Base_t * Ai=A+(mun-1)*lda;
				for(size_t i=mun; i>0; --i, Ai-=lda)
					std::copy(Ai, Ai+nun, Ai+1+lda);

				std::copy(dc, dc+n, A);

				Base_t * Aipo = A+lda;
				for(size_t i=0; i<mun; ++i, Aipo+=lda) *Aipo = b[i];

				delete [] dc;
				delete [] b;
			} else if (n != 0) {
				Base_t * Ai=A+mun*lda;
				Base_t d = *Ai;
				for(; Ai != A; Ai-=lda) *Ai= *(Ai-lda);
				*A=d;
			}
		} else {
			if ((m!=0) && (n > 1)) {
				const size_t nun(n-1);
				Base_t d = A[nun];
				std::copy(A,A+nun,A+1);
				*A=d;
			}
		}

#endif
	}

	template<typename Base_t>
	inline void cyclic_shift_row(Base_t * A, size_t m, size_t n, size_t lda)
	{

#ifdef MEMCOPY
		if (m > 1) {
			const size_t mun(m-1);

			Base_t * b = new Base_t[n];
			Base_t * Ai = A+mun*lda;
			memcpy (b,Ai,n*sizeof(Base_t));

			for(Base_t * Ac = A+mun*lda; Ac!=A;Ac-=lda)
				memcpy (Ac, Ac-lda, n*sizeof(Base_t));

			memcpy ( A, b, n*sizeof(Base_t));
			delete [] b;
		}

#else
		if (m > 1) {
			const size_t mun(m-1);

			Base_t * b = new Base_t[n];
			Base_t * Ai = A+mun*lda;
			for(size_t i=0; i<n; ++i, Ai+=1) b[i] = *Ai;

			for(Base_t * Ac = A+mun*lda; Ac!=A;Ac-=lda)
                std::copy(Ac-lda,Ac-lda+n, Ac);

			Base_t * Aii = A;
			for(size_t i=0; i<n; ++i, Aii+=1) *Aii = b[i];

			delete [] b;
		}

#endif
	}

	template<typename Element>
	inline void cyclic_shift_col(Element * A, size_t m, size_t n, size_t lda)
	{
		if (n > 1) {
			const size_t nun(n-1);
			for(Element*Ai=A; Ai!= A+m*lda; Ai+=lda)
			{
				Element tmp = Ai[nun];
				std::copy_backward(Ai, Ai+nun, Ai+n);
				*Ai=tmp;
			}
		}
	}



#ifdef __FFLASFFPACK_USE_OPENMP
	template<class Field>
	void
	papplyP( const Field& F,
		 const FFLAS::FFLAS_SIDE Side,
		 const FFLAS::FFLAS_TRANSPOSE Trans,
		 const size_t m, const int ibeg, const int iend,
		 typename Field::Element * A, const size_t lda, const size_t * P )
	{
		int numthreads = omp_get_max_threads();
		size_t BLOCKSIZE=std::max(2*m/numthreads,(size_t)1); // Assume that there is at least 2 ApplyP taking place in parallel
		size_t NBlocks = m/BLOCKSIZE;
		size_t LastBlockSize = m % BLOCKSIZE;
		if (LastBlockSize)
			NBlocks++;
		else
			LastBlockSize=BLOCKSIZE;
		for (size_t t = 0; t < NBlocks; ++t)
		{
			size_t BlockDim = BLOCKSIZE;
			if (t == NBlocks-1)
				BlockDim = LastBlockSize;
#pragma omp task shared (A, P, F) firstprivate(BlockDim)
			applyP(F, Side, Trans, BlockDim, ibeg, iend, A+BLOCKSIZE*t*((Side == FFLAS::FflasRight)?lda:1), lda, P);
		}
#pragma omp taskwait
	}

	template <class Element>
	void papplyT (Element* A, const size_t lda, const size_t width,
		      const size_t N2,
		      const size_t R1, const size_t R2,
		      const size_t R3, const size_t R4)
	{
		int numthreads = omp_get_max_threads();
		size_t BLOCKSIZE=std::max(width/numthreads,(size_t)1);
		size_t NBlocks = width/BLOCKSIZE;
		size_t LastBlockSize = width % BLOCKSIZE;
		if (LastBlockSize)
			NBlocks++;
		else
			LastBlockSize=BLOCKSIZE;
		for (size_t t = 0; t < NBlocks; ++t)
		{
			size_t BlockDim = BLOCKSIZE;
			if (t == NBlocks-1)
				BlockDim = LastBlockSize;
#pragma omp task shared (A) firstprivate(BlockDim,R1,R2,R3,R4,N2)
			applyT (A+BLOCKSIZE*t*lda, lda, BlockDim, N2, R1, R2, R3, R4);
		}
#pragma omp taskwait
	}


	template <class Element>
	void papplyS (Element* A, const size_t lda, const size_t width,
		      const size_t M2,
		      const size_t R1, const size_t R2,
		      const size_t R3, const size_t R4)
	{
		int numthreads = omp_get_max_threads();
		size_t BLOCKSIZE=std::max(width/numthreads,(size_t)1);
		size_t NBlocks = width/BLOCKSIZE;
		size_t LastBlockSize = width % BLOCKSIZE;
		if (LastBlockSize)
			NBlocks++;
		else
			LastBlockSize=BLOCKSIZE;
		for (size_t t = 0; t < NBlocks; ++t)
		{
			size_t BlockDim = BLOCKSIZE;
			if (t == NBlocks-1)
				BlockDim = LastBlockSize;
#pragma omp task shared (A) firstprivate(BlockDim,R1,R2,R3,R4,M2)
			applyS(A+BLOCKSIZE*t, lda, BlockDim, M2, R1, R2, R3, R4);
		}
#pragma omp taskwait
	}

#endif // __FFLASFFPACK_USE_OPENMP

} // FFPACK

#endif // __FFLASFFPACK_ffpack_permutation_INL
