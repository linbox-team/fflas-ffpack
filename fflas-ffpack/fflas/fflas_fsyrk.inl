/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2016 the FFLAS-FFPACK group
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

#ifndef __FFLASFFPACK_fflas_fsyrk_INL
#define __FFLASFFPACK_fflas_fsyrk_INL

namespace FFLAS {
    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::ConstElement_ptr A, const size_t lda,
           const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc){
        
            //@TODO: write an optimized iterative basecase
        if (N==1){ // Base case
            F.mulin (*C, beta);
			size_t incA = (trans==FFLAS::FflasNoTrans)?1:lda;
            F.axpyin (*C, alpha, fdot (F, K, A, incA, A, incA));
            return C;
        } else if (K==1){
            if (!F.isOne(beta))
                fscalin (F, N, N, beta, C, ldc);
			size_t incA = (trans==FFLAS::FflasNoTrans)?lda:1;
            fger (F, N, N, alpha, A, incA, A, incA, C, ldc);
            return C;
        } else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            size_t K1 = K>>1;
            size_t K2 = K - K1;
                // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans
			size_t incRow,incCol;
			FFLAS_TRANSPOSE oppTrans;
			if (trans==FflasNoTrans) {incRow=lda,incCol=1;oppTrans=FflasTrans;}
			else {incRow = 1; incCol = lda;oppTrans=FflasNoTrans;}

            typename Field::ConstElement_ptr A21 = A + N1*incRow;
            typename Field::ConstElement_ptr A22 = A21 + K1*incCol;
            typename Field::ConstElement_ptr A12 = A + K1*incCol;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;
                // C11 <- alpha A11 x A11^T + beta C11
            fsyrk (F, UpLo, trans, N1, K1, alpha, A, lda, beta, C, ldc);
                // C11 <- alpha A12 x A12^T + C11
            fsyrk (F, UpLo, trans, N1, K2, alpha, A12, lda, F.one, C, ldc);
                // C22 <- alpha A21 x A21^T + beta C22
            fsyrk (F, UpLo, trans, N2, K1, alpha, A21, lda, beta, C22, ldc);
                // C22 <- alpha A22 x A22^T + C22
            fsyrk (F, UpLo, trans, N2, K2, alpha, A22, lda, F.one, C22, ldc);

            if (UpLo == FflasUpper) {
						// C12 <- alpha A11 * A21^T + beta C12
					fgemm (F, trans, oppTrans, N1, N2, K1, alpha, A, lda, A21, lda,
						   beta, C12, ldc);
						// C12 <- alpha A12 * A22^T + C12
					fgemm (F, trans, oppTrans, N1, N2, K2, alpha, A12, lda, A22, lda,
						   F.one, C12, ldc);
			} else {
					// C21 <- alpha A21 * A11^T + beta C21
				fgemm (F, trans, oppTrans, N2, N1, K1, alpha, A21, lda, A, lda,
					   beta, C21, ldc);
					// C21 <- alpha A12 * A22^T + C21
				fgemm (F, trans, oppTrans, N2, N1, K2, alpha, A22, lda, A12, lda,
					   F.one, (UpLo==FflasUpper)?C12:C21, ldc);

			}
			return C;
        }
    }
    template<class Field>
    inline typename Field::Element_ptr
    fsyrk (const Field& F,
           const FFLAS_UPLO UpLo,
           const FFLAS_TRANSPOSE trans,
           const size_t N,
           const size_t K,
           const typename Field::Element alpha,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::ConstElement_ptr D, const size_t incD,
		   const typename Field::Element beta,
           typename Field::Element_ptr C, const size_t ldc){

            //@TODO: write an optimized iterative basecase
        if (N==1){ // Base case
            F.mulin (*C, beta);
			size_t incA = (trans==FFLAS::FflasNoTrans)?1:lda;
			typename Field::Element acc;
			F.init(acc,F.zero);
			for (typename Field::ConstElement_ptr Ai=A, Di=D; Ai != A+k*incA; k++, Ai+=incA,Di+=incD){
				F.assign (tmp, *Ai);
				F.mulin (*Ai, *Di);
				F.axpyin (acc, tmp, *Ai);}
			F.axpyin (*C, alpha, acc));
            return C;
        } else if (K==1){
            if (!F.isOne(beta))
                fscalin (F, N, N, beta, C, ldc);
			size_t incA = (trans==FFLAS::FflasNoTrans)?lda:1;
			typename Field::Element tmp;
			F.init(tmp, alpha);
			F.mulin (tmp, *D);
			fger (F, N, N, tmp, A, incA, A, incA, C, ldc);
			fscalin (F, N, *D, A, incA);
            return C;
        } else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
            size_t K1 = K>>1;
            size_t K2 = K - K1;
                // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans
			size_t incRow,incCol;
			FFLAS_TRANSPOSE oppTrans;
			if (trans==FflasNoTrans) {incRow=lda;incCol=1;oppTrans=FflasTrans;}
			else {incRowA = 1; incColA = lda;oppTrans=FflasNoTrans;}

            typename Field::ConstElement_ptr A21 = A + N1*incRow;
            typename Field::ConstElement_ptr A22 = A21 + K1*incCol;
            typename Field::ConstElement_ptr A12 = A + K1*incCol;
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;
            typename Field::Element_ptr D2 = D + K1*incD;
                // C11 <- alpha A11 x D1 x A11^T + beta C11
            fsyrk (F, UpLo, trans, N1, K1, alpha, A, lda, D, incD, beta, C, ldc);
                // C11 <- alpha A12 x D2 x A12^T + C11
            fsyrk (F, UpLo, trans, N1, K2, alpha, A12, lda, D2, incD, F.one, C, ldc);
                // C22 <- alpha A21 x D1 x A21^T + beta C22
            fsyrk (F, UpLo, trans, N2, K1, alpha, A21, lda, D, incD, beta, C22, ldc);
                // C22 <- alpha A22 x D2 x A22^T + C22
            fsyrk (F, UpLo, trans, N2, K2, alpha, A22, lda, D2, incD, F.one, C22, ldc);

            if (UpLo == FflasUpper) {
					// A21 <- A21 x D1
				typename Field::Element_ptr temp = fflas_new (F, N2*K1);
				if (trans==FflasNoTrans) {ldt=K1; incRowT=ldt; incColT=1;} else {ldt = N2; incRowT=1; incColT=ldt;}
				for (typename Field::Element_ptr Ai = A21, Ti = temp, Di = D;
					 Ai != A21 + K1*incCol;
					 Ai += incCol, Ti += incColT, Di+=incD){
					typename Field::Element inv; F.init(inv);
					F.inv(inv, *Di);
					fscal (F, N2, inv, Ai, incRow, Ti, incRowT);
				}
					// C12 <- alpha A11 x A21^T + beta C12
				fgemm (F, trans, oppTrans, N1, N2, K1, alpha, A, lda, temp, ldt, beta, C12, ldc);
				fflas_delete (temp);
					// A12 <- A12 x D2
				temp = fflas_new (F, K2*N1);
				if (trans==FflasNoTrans) {ldt=K2; incRowT=ldt; incColT=1;} else {ldt = N1; incRowT=1; incColT=ldt;}
				for (typename Field::Element_ptr Ai = A12, Ti = temp, Di = D2;
					 Ai != A12 + K2*incCol;
					 Ai += incCol, Ti += incColT, Di+=incD){
					typename Field::Element inv; F.init(inv);
					F.inv(inv, *Di);
					fscal (F, N1, inv, Ai, incRow, temp, ldt);
				}
					// C12 <- alpha A12 x A22^T + C12
				fgemm (F, trans, oppTrans, N1, N2, K2, alpha, temp, ldt, A22, lda, F.one, C12, ldc);
				fflas_delete (temp);
			} else {
					// C21 <- alpha A21 x D1 x A11^T + beta C21
				fgemm (F, trans, oppTrans, N2, N1, K1, alpha, A21, lda, A, lda, beta, C21, ldc);
					// C21 <- alpha A12 x D2 x A22^T + C21
				fgemm (F, trans, oppTrans, N2, N1, K2, alpha, A22, lda, A12, lda, F.one, (UpLo==FflasUpper)?C12:C21, ldc);

			}
			return C;
        }
    }
}

#endif //__FFLASFFPACK_fflas_fsyrk_INL
