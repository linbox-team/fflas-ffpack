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
//                        Test for rank
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/ffpack/ffpack.h"



using namespace std;
using namespace FFPACK;

typedef Givaro::Modular<double> Field;

int main(int argc, char** argv){

    int n,m;
    cerr<<setprecision(10);
    if (argc !=  3)	{
        cerr<<"Usage : test-fullranksubmatrix <p> <A> <<i>"
        <<endl;
        exit(-1);
    }
    Field F (atoi(argv[1]));
    Field::Element * A;
    Field::Element * X;

    FFLAS::ReadMatrix (argv[2],F,m,n,A);
    FFMAS::WriteMatrix (cerr<<"A = "<<endl, F, m, n, A, n);

    FFLAS::Timer tim,t; t.clear();tim.clear();
    size_t R;

    FFPACK::ColRankProfileSubmatrix (F, m, n, A, n, X, R);

    FFLAS::WriteMatrix (cerr<<"X = "<<endl, F, R, R, X, R);

    size_t r2 = FFPACK::Rank(F, R,R, X, R);
    if (r2 != R)
        std::cerr<<"Fail : Rank (X) != Rank (A)"<<std::endl;



    FFLAS::fflas_delete(X);
    FFLAS::fflas_delete(A);
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
