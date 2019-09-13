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


#include <iostream>
#include "givaro/modular.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/ffpack/ffpack.h"

using namespace std;
using namespace FFLAS;
using namespace FFPACK;

bool checkRPM (size_t M, size_t N, size_t R){

    size_t * rows = fflas_new<size_t>(R);
    size_t * cols = fflas_new<size_t>(R);

    RandomRankProfileMatrix (M, N, R, rows, cols);

    std::sort(rows,rows+R);
    std::sort(cols,cols+R);
    bool ok = true;
    // checking that all row and col indices are distinct
    for (size_t i=0; i<R-1; i++){
        if (rows[i] == rows[i+1] || cols[i] == cols[i+1]){
            ok =false;
            cerr<<"Error in RPM: rows["<<i<<"] = "<<rows[i]<<endl
            <<"              rows["<<i+1<<"] = "<<rows[i+1]<<endl
            <<"              cols["<<i<<"] = "<<cols[i]<<endl
            <<"              cols["<<i+1<<"] = "<<cols[i+1]<<endl;
        }
    }
    fflas_delete(rows,cols);
    return ok;
}
bool checkSymmetricRPM (size_t N, size_t R){

    size_t * rows = fflas_new<size_t>(R);
    size_t * cols = fflas_new<size_t>(R);

    bool ok = true;

    RandomSymmetricRankProfileMatrix (N, R, rows, cols);

    // checking that the matrix is symmetric
    Givaro::Modular<float> F(11);
    float * A=fflas_new(F,N*N);
    fzero(F,N,N,A,N);
    for (size_t i=0; i<R; i++)
        F.assign(A[rows[i]*N+cols[i]],F.one);
    for (size_t i=0; i<N; i++)
        for (size_t j=i+1; j<N; j++)
            if (!F.areEqual(A[i*N+j],A[j*N+i])){
                ok=false;
                cerr<<"A["<<i<<", "<<j<<"] = "<<A[i*N+j]<<" != A["<<j<<", "<<i<<"]"<<std::endl;
            }
    fflas_delete(A);

    if (!ok){
        std::cerr<<"RPM FAILED"<<std::endl;
        return false;
    }

    // checking that all row and col indices are distinct
    std::sort(rows,rows+R);
    std::sort(cols,cols+R);

    for (size_t i=0; i<R-1; i++){
        if (rows[i] == rows[i+1] || cols[i] == cols[i+1]){
            ok =false;
            cerr<<"Error in SymRPM: rows["<<i<<"] = "<<rows[i]<<endl
            <<"                 rows["<<i+1<<"] = "<<rows[i+1]<<endl
            <<"                 cols["<<i<<"] = "<<cols[i]<<endl
            <<"                 cols["<<i+1<<"] = "<<cols[i+1]<<endl;
        }
    }
    fflas_delete(rows,cols);

    if (!ok) std::cerr<<"Symmetric RPM FAILED"<<std::endl;

    return ok;
}

int main(int argc, char** argv){
    size_t m=120;
    size_t n=120;
    size_t r=70;
    size_t iters=3;
    bool loop=false;
    uint64_t seed = getSeed();

    Argument as[] = {
        { 'm', "-m M", "Set the row dimension of the matrix.",      TYPE_INT , &m },
        { 'n', "-n N", "Set the column dimension of the matrix.", TYPE_INT , &n },
        { 'r', "-r R", "Set the rank.", TYPE_INT , &r },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    if (r > std::min (m,n))
        r = std::min (m, n);

    srand(seed);

    bool ok = true;
    size_t it = iters;
    do{
        ok = ok && checkRPM(m,n,r);
        ok = ok && checkSymmetricRPM(n,r);
        it--;
    } while ((loop || it)  && ok);

    if (!ok)
        std::cerr<<"with seed = "<<seed<<std::endl;
    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
