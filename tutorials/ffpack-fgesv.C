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


    // Let A be a M times M square matrix of coefficients
    const size_t M = 4, lda = M;


    Field::Element_ptr A;
    A = fflas_new(F,M,M);


    // Fulfill the square matrix A so that A is invertible
    F.assign(A[0], F.one);
    F.assign (A[1],F.zero);
    F.assign(A[2],F.one);
    F.assign (A[3],F.zero);
    F.assign(A[4],F.zero);
    F.assign (A[5],F.one);
    F.assign(A[6],F.zero);
    F.assign (A[7],F.one);
    F.assign(A[8],F.zero);
    F.assign (A[9],F.zero);
    F.assign(A[10],F.one);
    F.assign (A[11],F.zero);
    F.assign(A[12],F.zero);
    F.assign (A[13],F.zero);
    F.assign(A[14],(F.zero));
    F.assign (A[15],F.one);



    WriteMatrix(std::cout<<"A:="<<std::endl,F,M,M,A,lda)<<std::endl;

    // Let X be a M times 2 matrix of variables
    const size_t ldx = 2;
    Field::Element_ptr x;
    x = fflas_new(F,M,2);
    fiszero (F, M, 2, x, ldx); //initialize all elements to zero

    WriteMatrix(std::cout<<"x:="<<std::endl,F,M, 2, x, ldx)<<std::endl;

    // Let b be a M times 2 matrix of solutions
    const size_t ldb = 2;
    Field::Element_ptr b;
    b = fflas_new(F,M,2);


    // Fulfill the matrix b with desired values
    F.init(b[0],4);
    F.init(b[1],4);
    F.init(b[2],6);
    F.init(b[3],3);
    F.init(b[4],3);
    F.init(b[5],6);
    F.init(b[6],4);
    F.init(b[7],4);

    WriteMatrix(std::cout<<"b:="<<std::endl,F,M, 2, b, ldb)<<std::endl;

    // make a copy of b into x
    fassign(F,M,2,b,ldb,x,2);
    WriteMatrix(std::cout<<"copied b:="<<std::endl,F,M, 2, x, ldx)<<std::endl;

    //Solve the system
    int state;
    size_t rank = fgesv(F, FflasLeft, M, 2, A, lda, x, ldx, &state);
    if(rank!=M)std::cout<<"Results are incorrect after the fgesv()!"<<std::endl;
    WriteMatrix(std::cout<<"x:="<<std::endl,F,M, 2, x, ldx)<<std::endl;



    // Let res be a M times 2 matrix
    const size_t ldres = 2;
    Field::Element_ptr res;
    res = fflas_new(F,M,2);
    fiszero (F, M, 2, res, ldres); //initialize all elements to zero

    // Verify if A*x == b to confirm the found the solution
    std::cout<<"Verification:"<<std::endl;
    fgemm(F, FflasNoTrans, FflasNoTrans, M, 2, M, F.one, A, lda, x, ldx, F.zero, res, ldres);
    WriteMatrix(std::cout<<"A*x:="<<std::endl,F,M,2,res,ldres)<<std::endl;

    if( !fequal (F, M, 2, res, ldres, b, ldb) ) {
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

    return 0;

}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
