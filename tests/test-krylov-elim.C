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

#include <iostream>
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/timer.h"
using namespace std;
#include "givaro/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"


using namespace FFPACK;
typedef Givaro::Modular<double> Field;

template<class T>
std::ostream& printvect(std::ostream& o, T* vect, size_t dim)
{
    for(size_t i=0; i<dim; ++i)
        o << vect[i] << " " ;
    return o << std::endl;
}

int main(int argc, char** argv){

    size_t m,n;


    if (argc!=3){
        cerr<<"usage : test-lqup <p> <A>"<<endl
        <<"         to compute the rank profile of the (n+m)xn matrix B formed by the n identity vectors and the mxn matrix A over Z/pZ"
        <<endl;
        exit(-1);
    }
    Field F(atoi(argv[1]));
    Field::Element* A;
    FFLAS::ReadMatrix (argv[2],F,m,n,A);

    Field::Element * B = FFLAS::fflas_new<Field::Element>((m+n)*n);
    for (size_t i=0; i<(n+m)*n;++i) *(B+i)=0;

    size_t deg = (n-1)/m+1;
    size_t curr_row = 0;
    size_t it_idx = 0;
    size_t bk_idx = 0;
    for (size_t i=0; i<m; ++i){
        for (size_t j=0; j<deg; ++j){
            if (curr_row < n+m -1){
                F.assign( *(B + curr_row*n + n-1 - it_idx), F.one);
                curr_row++;
                it_idx++;
            }
        }
        for (size_t j=0; j<n; ++j)
            *(B + curr_row*n + j) = *(A + bk_idx*n + j);
        bk_idx++;
        curr_row++;
    }
    FFLAS::WriteMatrix (cout<<"A = "<<endl, F, m, n, A, n);
    FFLAS::WriteMatrix (cout<<"B = "<<endl, F, (m+n), n,B, n);

    size_t *rp = FFLAS::fflas_new<size_t>(n);

    FFPACK::SpecRankProfile(F, m, n, A, n, deg,rp);

    size_t * P = FFLAS::fflas_new<size_t>(n);
    size_t * Q = FFLAS::fflas_new<size_t>(n+m);
    FFPACK::LUdivine(F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,(int)m+n, n, B, n, P, Q);

    printvect (cout<<"RankProfile (A) = "<<endl, rp, n)<<endl;

    printvect (cout<<"RankProfile (B) = "<<endl, Q, n)<<endl;

    FFLAS::fflas_delete( rp );
    FFLAS::fflas_delete( A );
    FFLAS::fflas_delete( B );

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
