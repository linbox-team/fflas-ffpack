/* Copyright (c) FFLAS-FFPACK
 * Written by Philippe LEDENT <philippe.ledent@etu.univ-grenoble-alpes.fr>
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
 */

#include <fflas-ffpack/fflas/fflas.h>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/fflas_io.h"
#include <iostream>

using namespace FFLAS;

int main(int argc, char** argv) {
    std::cout << "" << std::endl;
    typedef Givaro::Modular<float> Float_Field;
    Float_Field F(101);
    uint64_t seed = time(NULL);
    typename Float_Field::RandIter G(F,seed);

    // Let p and q be naturals
    const size_t p = 7, q = 11, lda = q;
    // Let A be a p times q random matrix
    Float_Field::Element_ptr A;
    A = fflas_new(F,p,q);
    frand(F,G,p,q,A,lda);

    // Let m and n be naturals such that (0 < m <= p && 0 < n <= q)
    const size_t m = 5, n = 10;
    // Let B be a m times n sub matrix of A where B[0][0] = A[1][2]
    const size_t ldb = lda;
    Float_Field::Element_ptr B = A + 1*lda + 2;

    // Let C be a m times n random matrix
    const size_t ldc = m;
    Float_Field::Element_ptr C;
    C = fflas_new(F,m,n);
    frand(F,G,m,n,C,ldc);

    // Let D be a m times n matrix
    // where D := B - C
    const size_t ldd = m;
    Float_Field::Element_ptr D;
    D = fflas_new(F,m,n);
    fsub(F,m,n,B,ldb,C,ldc,D,ldd);

    // Similar routines :
    // add(F,m,n,B,ldb,C,ldc,D,ldd); // D := B + C
    // subin(F,m,n,B,ldb,C,ldc);     // C := B - C
    // addin(F,m,n,B,ldc,C,ldc);     // C := B + C

    // if B == C then D == zeroMatrix
    bool res = fiszero(F,m,n,D,ldd);


    //Output
    WriteMatrix(std::cout<<"A:="<<std::endl,F,p,q,A,lda)<<std::endl;
    WriteMatrix(std::cout<<"B:="<<std::endl,F,m,n,B,ldb)<<std::endl;
    WriteMatrix(std::cout<<"C:="<<std::endl,F,m,n,C,ldc)<<std::endl;
    WriteMatrix(std::cout<<"D:="<<std::endl,F,m,n,D,ldd)<<std::endl;
    std::cout<<"Is D the zero matrix ?"<< std::endl;
    if(res)
        std::cout<<"TRUE"<< std::endl;
    else
        std::cout<<"FALSE"<< std::endl;

    // Clearing up the memory
    fflas_delete(A);
    fflas_delete(C);
    fflas_delete(D);

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
