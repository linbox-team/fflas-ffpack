/* Copyright (c) FFLAS-FFPACK
 * Written by ZHU Hongguang <zhuhongguang2014@gmail.com>
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
#include <fflas-ffpack/ffpack/ffpack.h>
#include <iostream>

using namespace FFLAS;
using namespace FFPACK;

int main(int argc, char** argv) {

    typedef Givaro::Modular<float> Field;
    Field F(17);



    // Let A be a M times M square matrix
    const size_t M = 4, lda = M;


    // Let A be a M times M random square matrix
    Field::Element_ptr A;
    A = fflas_new(F,M,M);
    Field::Element_ptr A2;
    A2 = fflas_new(F,M,M);

    // Fulfill the square matrix A so that A is invertible
    F.assign(A[0], F.one);
    F.assign(A[1],F.zero);
    F.assign(A[2],F.zero);
    F.assign(A[3],F.zero);
    F.init(A[4],15);
    F.assign(A[5],F.one);
    F.assign(A[6],F.zero);
    F.init(A[7],3);
    F.assign(A[8],F.zero);
    F.assign(A[9],F.zero);
    F.assign(A[10],F.one);
    F.init(A[11],2);
    F.init(A[12],10);
    F.assign(A[13],F.zero);
    F.assign(A[14],(F.zero));
    F.assign(A[15],F.one);

    FFLAS::fassign (F, M, M, A, lda, A2, lda);


    // Print out matrix A to verify
    WriteMatrix(std::cout<<"A:="<<std::endl,F,M,M,A,lda)<<std::endl;

    // Let x be a M dimensional vector
    const size_t incx = 1;
    Field::Element_ptr x;
    x = fflas_new(F,M,1);
    fzero(F,M,x,incx);//initialize all elements to zero

    // Let b be a M dimensional vector
    const size_t incb = 1;
    Field::Element_ptr b;
    b = fflas_new(F,M,1);

    // Fulfill the vector b with desired values
    F.init(b[0],1);
    F.init(b[1],3);
    F.init(b[2],6);
    F.init(b[3],5);

    // Print out matrix A to verify
    WriteMatrix(std::cout<<"b:="<<std::endl,F,M, 1, b, incb)<<std::endl;


    //Solve the linear system Ax=b for x
    /**
PS: the function Solve will modify the matrix A so here we used a duplicate matrix A2 otherwise A*x will not be equal to b for the later verification stage
*/
    Solve( F, M, A2, lda, x, incx, b, incb );

    // Print out x to verify
    WriteMatrix(std::cout<<"x:="<<std::endl,F,M, 1, x, incx)<<std::endl;




    // Let res be a M times 1 vector
    const size_t incres = 1;
    Field::Element_ptr res;
    res = fflas_new(F,M,1);
    fzero(F,M,res,incres);//initialize all elements to zero

    // Verify if A*x == b to confirm the found the solution
    std::cout<<"Verification:"<<std::endl;
    fgemv(F, FflasNoTrans, M, M, F.one, A, lda, x, incx, F.zero, res, incres);
    WriteMatrix(std::cout<<"A*x:="<<std::endl,F,M,1,res,incres)<<std::endl;
    if( !fequal (F, M, res, incres, b, incb)  )
    {
        std::cout<<"Results are incorrect!"<<std::endl;
    }
    else
    {
        std::cout<<"Results are correct!"<<std::endl;
    }



    // Clearing up the memory
    fflas_delete(A);
    fflas_delete(x);
    fflas_delete(b);
    fflas_delete(res);
    fflas_delete(A2);
    return 0;

}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
