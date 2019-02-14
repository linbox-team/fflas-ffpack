/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Clément Pernet
 * This file is Free Software and part of FFLAS-FFPACK.
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


#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

Givaro::Timer tperm, tgemm, tBC, ttrsm,trest,timtot;

#include "fflas-ffpack/ffpack/ffpack.h"

using namespace std;
using namespace FFLAS;
using namespace FFPACK;
using Givaro::Modular;

bool checkMonotonicApplyP(FFLAS_SIDE Side, FFLAS_TRANSPOSE trans, size_t * P, size_t N, size_t R){
    bool ok = true;

    typedef Modular<double> Field;
    Field F(101);
    size_t M = 2;
    size_t lda = (Side == FflasLeft)? M : N;
    size_t ldb = lda;
    Field::Element_ptr A = fflas_new(F, M, N);
    Field::Element_ptr B = fflas_new(F, M, N);
    if (Side == FflasLeft) {
        for (size_t i = 0; i<N; ++i) {
            for (size_t j = 0; j<M; ++j) {
                F.init(A[i*lda+j], i*10+j);
            }
        }
    } else {
        for (size_t i = 0; i<N; ++i) {
            for (size_t j = 0; j<M; ++j) {
                F.init(A[i+j*lda], i*10+j);
            }
        }
    }

    fassign(F, N,M, A, lda, B, ldb);

    MonotonicApplyP (F, FflasLeft, FflasNoTrans, M, 0, N, A, lda, P, R);

    // checking that cols have not been permuted
    typename Field::Element x;
    F.init(x);
    for (size_t i=0; i<N; i++){
        F.sub(x,A[lda*i+1],A[lda*i]);
        if (!F.isOne(x)){
            std::cerr<<"ERROR: A["<<i<<", 1] = "<<A[i*lda+1]<<" != "<<A[i*lda]+1<<" = A["<<i<<", 0]+1"<<std::endl;
            ok = false;
        }
    }

    // Checking that the non pivot rows are monotonically increasing
    for (size_t i=R; i<N-1; i++){
        if (A[i*lda] >= A[(i+1)*lda]){
            std::cerr<<"ERROR: A["<<i<<", 0] = "<<A[i*lda]<<" >= "<<A[(i+1)*lda]<<" = A["<<i+1<<", 0]"<<std::endl;
            ok = false;
        }
    }

    // Checking that the first R rows have been permuted correctly
    applyP(F, Side, trans, M, 0, R, B, ldb, P);
    if (!fequal(F, R, M,A, lda, B, ldb)){
        std::cerr<<"ERROR: first R rows are not permuted correctly"<<std::endl;
        ok =false;
    }
    fflas_delete(A);
    fflas_delete(B);

    return ok;
}

int main(){


    bool ok = true;

    size_t  P1[10] = {0,5,6,6,7,9,6,7,8,9};
    ok = ok && checkMonotonicApplyP(FflasLeft, FflasNoTrans, P1, 10, 6);
    size_t  P2[10] = {0,3,3,6,6,5,6,7,8,9};
    ok = ok && checkMonotonicApplyP(FflasLeft, FflasNoTrans, P2, 10, 5);
    size_t  P3[10] = {0,4,2,4,5,5,6,7,8,9};
    ok = ok && checkMonotonicApplyP(FflasLeft, FflasNoTrans, P3, 10, 6);

    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
