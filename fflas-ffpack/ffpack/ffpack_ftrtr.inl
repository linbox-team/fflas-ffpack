/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 FFLAS-FFACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_ffpack_ftrtr_INL
#define __FFLASFFPACK_ffpack_ftrtr_INL


namespace FFPACK {


template<class Field>
	void
	ftrtri (const Field& F, const FFLAS::FFLAS_UPLO Uplo, const FFLAS::FFLAS_DIAG Diag,
		const size_t N, typename Field::Element_ptr A, const size_t lda)
	{
		if (N == 1){
			if (Diag == FFLAS::FflasNonUnit)
				F.invin (*A);
		}
		else {
			size_t N1 = N/2;
			size_t N2 = N - N1;
			ftrtri (F, Uplo, Diag, N1, A, lda);
			ftrtri (F, Uplo, Diag, N2, A + N1*(lda+1), lda);
			if (Uplo == FFLAS::FflasUpper){
				ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, N1, N2,
				       F.one, A, lda, A + N1, lda);
				ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, N1, N2,
				       F.mOne, A + N1*(lda+1), lda, A + N1, lda);
			}
			else {
				ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, N2, N1,
				       F.one, A + N1*(lda+1), lda, A + N1*lda, lda);
				ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, N2, N1,
				       F.mOne, A, lda, A + N1*lda, lda);
			}
		}
	}


template<class Field, class ParSeq>
	void
	ftrtri (const Field& F, const FFLAS::FFLAS_UPLO Uplo, const FFLAS::FFLAS_DIAG Diag,
		const size_t N, typename Field::Element_ptr A, const size_t lda, ParSeq& H)
	{
        if ( (N == 1) || (H.numthreads()<=1)) return ftrtri (F,Uplo,Diag,N,A,lda);
        SYNCH_GROUP(			
            size_t N1 = N/2;
            size_t N2 = N - N1;
            
                // [ A1 | A3 ]
                // [  0 | A2 ]
            typename Field::Element * A1 = A;
            typename Field::Element * A2 = A + N1*(lda+1);
            ParSeq Hleft(std::max(H.numthreads()/2,size_t(1)));
            ParSeq Hright(std::max(H.numthreads()-Hleft.numthreads(),size_t(1)));

            TASK(MODE(CONSTREFERENCE(F, A1) READWRITE(A1[0])),
                 ftrtri (F, Uplo, Diag, N1, A1, lda, Hleft); );

            TASK(MODE(CONSTREFERENCE(F, A2) READWRITE(A2[0])),
                 ftrtri (F, Uplo, Diag, N2, A2, lda, Hright); );

            CHECK_DEPENDENCIES;

            if (Uplo == FFLAS::FflasUpper){
                typename Field::Element * A3 = A + N1;

                TASK(MODE(READ(A1[0]) READWRITE(A3[0]) CONSTREFERENCE(F, A1, A3)),
                ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, 
                       N1, N2, F.one, A1, lda, A3, lda, H); );

                TASK(MODE(READ(A2[0]) READWRITE(A3[0]) CONSTREFERENCE(F, A2, A3)),
                ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, 
                       N1, N2, F.mOne, A2, lda, A3, lda, H); );
                
                CHECK_DEPENDENCIES;
            }
            else {
                typename Field::Element * A3 = A + N1*lda;

                TASK(MODE(READ(A2[0]) READWRITE(A3[0]) CONSTREFERENCE(F, A2, A3)),
                ftrmm (F, FFLAS::FflasLeft, Uplo, FFLAS::FflasNoTrans, Diag, 
                       N2, N1, F.one, A2, lda, A3, lda, H); );

                TASK(MODE(READ(A1[0]) READWRITE(A3[0]) CONSTREFERENCE(F, A1, A3)),
                ftrmm (F, FFLAS::FflasRight, Uplo, FFLAS::FflasNoTrans, Diag, 
                       N2, N1, F.mOne, A1, lda, A3, lda, H); );

                CHECK_DEPENDENCIES;
            }
        );   
	}


	template<class Field>
	void
	ftrtrm (const Field& F, const FFLAS::FFLAS_DIAG diag, const size_t N,
		typename Field::Element_ptr A, const size_t lda)
	{

		if (N == 1)
			return;
		size_t N1 = N/2;
		size_t N2 = N-N1;

		ftrtrm (F, diag, N1, A, lda);

		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N1, N1, N2, F.one,
		       A+N1, lda, A+N1*lda, lda, F.one, A, lda);

		ftrmm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans,
		       (diag == FFLAS::FflasUnit) ? FFLAS::FflasNonUnit : FFLAS::FflasUnit,
		       N1, N2, F.one, A + N1*(lda+1), lda, A + N1, lda);

		ftrmm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, diag, N2, N1,
		       F.one, A + N1*(lda+1), lda, A + N1*lda, lda);

		ftrtrm (F, diag, N2, A + N1*(lda+1), lda);

	}

	template<class Field, class ParSeq>
	void
	ftrtrm (const Field& F, const FFLAS::FFLAS_DIAG diag, const size_t N,
		typename Field::Element_ptr A, const size_t lda, ParSeq& H)
	{

		if (N == 1) return;
        if (H.numthreads()<=1) return ftrtrm(F,diag,N,A,lda);

		size_t N1 = N/2;
		size_t N2 = N-N1;

		ftrtrm (F, diag, N1, A, lda, H);

		fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N1, N1, N2, F.one,
		       A+N1, lda, A+N1*lda, lda, F.one, A, lda, H);

		ftrmm (F, FFLAS::FflasRight, FFLAS::FflasLower, FFLAS::FflasNoTrans,
		       (diag == FFLAS::FflasUnit) ? FFLAS::FflasNonUnit : FFLAS::FflasUnit,
		       N1, N2, F.one, A + N1*(lda+1), lda, A + N1, lda, H);

		ftrmm (F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, diag, N2, N1,
		       F.one, A + N1*(lda+1), lda, A + N1*lda, lda, H);

		ftrtrm (F, diag, N2, A + N1*(lda+1), lda, H);

	}

	template<class Field>
	void trinv_left( const Field& F, const size_t N, typename Field::ConstElement_ptr L, const size_t ldl,
			 typename Field::Element_ptr X, const size_t ldx )
	{
		FFLAS::fassign(F,N,N,L,ldl,X,ldx);
		ftrtri (F, FFLAS::FflasLower, FFLAS::FflasUnit, N, X, ldx);
		//invL(F,N,L,ldl,X,ldx);
	}

} // FFPACK

#endif // __FFLASFFPACK_ffpack_ftrtr_INL
