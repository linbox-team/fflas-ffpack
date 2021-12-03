/* Copyright (c) FFLAS-FFPACK
 * Written by Philippe LEDENT
 * philippe.ledent@etu.univ-grenoble-alpes.fr
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
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include <iostream>

using namespace FFLAS;

int main(int argc, char** argv) {
    std::cout << "" << std::endl;
    typedef Givaro::Modular<float> Float_Field;
    Float_Field F(101);

    // Let m be a natural
    const size_t m = 11, inca = 1;
    // Let a be a m by 1 random vector
    Float_Field::Element_ptr a;
    a = fflas_new(F,m,1);
    uint64_t seed = time(NULL);
    typename Float_Field::RandIter G(F,seed);
    frand(F,G,m,a,inca);

    // Let n be natural
    // Let b be a 1 by n random vector
    const size_t n = 13, incb = 1;
    Float_Field::Element_ptr b;
    b = fflas_new(F,1,n);

    frand(F,G,n,b,incb);

    // Let A be an m by n matrix obtained by the outer product between a and b
    const size_t lda = n;
    Float_Field::Element_ptr A;
    A = fflas_new(F,m,n);
    fger(F,m,n,F.one,a,inca,b,incb,A,lda);



    // Let C be an n by 1 vector such that C
    const size_t incc = 1;
    Float_Field::Element_ptr c;
    c = fflas_new(F,m,1);

    frand(F,G,m,c,incc);

    // Let d be a scalar where d = b dot c
    Float_Field::Element d;
    d = fdot(F,n,b,incb,c,incc);

    // Let e be a copy of a
    Float_Field::Element_ptr e;
    e = fflas_new(F,m,1);
    fassign(F,m,1,a,inca,e,inca);

    // Compute   e :=  (d scalar e) - (A times c)
    // Therefore e := -(A times c)  + (d scalar e)
    fgemv(F,FFLAS::FflasNoTrans,m,n,F.mOne,A,lda,c,incc,d,e,inca);

    // If e is the zero vector then

    // Is a the  zero vector ?
    bool res = fiszero(F,m,e,inca);


    //Output
    WriteMatrix(std::cout<<"a:=\n",F,m,1,a,inca)<<std::endl;
    WriteMatrix(std::cout<<"b:= ",F,1,n,b,incb)<<std::endl;
    WriteMatrix(std::cout<<"A:=\n",F,m,n,A,lda)<<std::endl;
    WriteMatrix(std::cout<<"c:=\n",F,m,1,c,incc)<<std::endl;
    //WriteMatrix(std::cout<<"d:=",F,1,1,d,1)<<std::endl;
    WriteMatrix(std::cout<<"e:=\n",F,m,1,e,inca)<<std::endl;
    std::cout<<"Is e the zero vector ?"<< std::endl;
    if(res)
        std::cout<<"TRUE"<< std::endl;
    else
        std::cout<<"FALSE"<< std::endl;

    // note :
    // There are many routines that create random matrices.
    // They can be found in fflas-ffpack/utils/fflas_randommatrix.





    // Clearing up the memory
    fflas_delete(a);
    fflas_delete(b);
    fflas_delete(A);
    fflas_delete(c);
    fflas_delete(e);

}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
