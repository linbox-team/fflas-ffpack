/*
 * Copyright (C) FFLAS-FFPACK
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


//--------------------------------------------------------------------------
//          Test for the krylov-elimination
//--------------------------------------------------------------------------
// usage: test-krylov-elim p A, to compute the rank profile of the (n+m)xn matrix B
// formed by the n identity vectors and the mxn matrix A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <list>
#include <vector>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace std;
typedef Givaro::Modular<double> Field;

//! @bug does not belong here
template<class T>
std::ostream& printvect(std::ostream& o, vector<T>& vect){
    for(size_t i=0; i < vect.size(); ++i)
        o << vect[i] << " " ;
    return o << std::endl;
}

int main(int argc, char** argv)
{


    static Argument as[] = {
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    // int m,n;

    Field F(65521);
    size_t N = 17;
    double * A = FFLAS::fflas_new<double>(N*N);
    double * tmp = FFLAS::fflas_new<double>(N*N);
    size_t * deg = FFLAS::fflas_new<size_t>(N);

    for (size_t i=0; i<(size_t)N*N; ++i)
        A[i] = 0;
    for (size_t i=0; i<3; ++i)
        A[i+i*N] = 1;

    for (size_t i=3; i<6; ++i)
        A[i+1+i*N] = 1;
    for (size_t i=6; i<9; ++i)
        A[i+2+i*N] = 1;

    A[12+9*N] = 1;
    A[13+10*N] = 1;
    A[14+12*N] = 1;
    A[15+15*N] = 1;
    A[16+16*N] = 1;
    deg[0]  = 4; deg[1] = 4; deg[2] = 4;deg[3] = 2; deg[4] = 1; deg[5] =2;
    for (size_t i=0; i<size_t(N); ++i)
        A[11+i*N] = A[7+i*N] = A[3+i*N] = double(i % 10);

    double * B = FFLAS::fflas_new<double>(N*N) ;
    FFLAS::fassign(F,N*N,A,1,B,1);

    FFPACK::Protected::CompressRowsQK (F, N, A+9*N, N, tmp, N, deg+3, 4, 3 );

    FFPACK::Protected::DeCompressRowsQK (F, N, N-9, A+9*N, N, tmp, N, deg+3, 4, 3 );

    int ok = 0 ;
    for (size_t i = 0 ; i < (size_t)N * (size_t)N ; ++i)
        if (A[i] != B[i])
        {
            ok = 1 ;
            break ;
        }

    FFLAS::fflas_delete( A );
    FFLAS::fflas_delete( tmp) ;
    FFLAS::fflas_delete(deg) ;
    FFLAS::fflas_delete(  B );

    return ok ;

}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
