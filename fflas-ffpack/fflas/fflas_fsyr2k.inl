/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2018 the FFLAS-FFPACK group
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

#ifndef __FFLASFFPACK_fflas_fsyr2k_INL
#define __FFLASFFPACK_fflas_fsyr2k_INL

namespace FFLAS {
    template<class Field>
    inline typename Field::Element_ptr
    fsyr2k (const Field& F,
			const FFLAS_UPLO UpLo,
			const FFLAS_TRANSPOSE trans,
			const size_t N,
			const size_t K,
			const typename Field::Element alpha,
			typename Field::ConstElement_ptr A, const size_t lda,
			typename Field::ConstElement_ptr B, const size_t ldb,
			const typename Field::Element beta,
			typename Field::Element_ptr C, const size_t ldc){
        
            //@TODO: write an optimized iterative basecase
        if (N==1){ // Base case
            F.mulin (*C, beta);
            size_t incA = (trans==FFLAS::FflasNoTrans)?1:lda;
            size_t incB = (trans==FFLAS::FflasNoTrans)?1:ldb;
            typename Field::Element two;
            F.init(two, 2);
            F.mulin (two,fdot (F, K, A, incA, B, incB));
            F.axpyin(*C, alpha, two);
            return C;
        } else {
            size_t N1 = N>>1;
            size_t N2 = N - N1;
                // Comments written for the case UpLo==FflasUpper, trans==FflasNoTrans
            FFLAS_TRANSPOSE oppTrans;
            if (trans==FflasNoTrans) {oppTrans=FflasTrans;}
            else {oppTrans=FflasNoTrans;}

            typename Field::ConstElement_ptr A2 = A + N1*(trans==FflasNoTrans?lda:1);
            typename Field::ConstElement_ptr B2 = B + N1*(trans==FflasNoTrans?ldb:1);
            typename Field::Element_ptr C12 = C + N1;
            typename Field::Element_ptr C21 = C + N1*ldc;
            typename Field::Element_ptr C22 = C12 + N1*ldc;
                // C11 <- alpha (A1 x B1^T + B1 x A1^T) + beta C11
            fsyr2k (F, UpLo, trans, N1, K, alpha, A, lda, B, ldb, beta, C, ldc);
                // C22 <- alpha (A2 x B2^T +B2 x A2^T) + beta C22
            fsyr2k (F, UpLo, trans, N2, K, alpha, A2, lda, B2, ldb, beta, C22, ldc);

            if (UpLo == FflasUpper) {
					// C12 <- alpha A1 * B2^T + beta C12
				fgemm (F, trans, oppTrans, N1, N2, K, alpha, A, lda, B2, ldb, beta, C12, ldc);
					// C12 <- alpha B1 * A2^T + C12
				fgemm (F, trans, oppTrans, N1, N2, K, alpha, B, ldb, A2, lda, F.one, C12, ldc);
            } else {
					// C21 <- alpha A2 * B1^T + beta C12
				fgemm (F, trans, oppTrans, N2, N1, K, alpha, A2, lda, B, ldb, beta, C21, ldc);
					// C21 <- alpha B2 * A1^T + C12
				fgemm (F, trans, oppTrans, N2, N1, K, alpha, B2, ldb, A, lda, F.one, C21, ldc);
            }
            return C;
        }
    }
} // namespace FFLAS
#endif //__FFLASFFPACK_fflas_fsyr2k_INL
