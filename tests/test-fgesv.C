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

#include <iostream>
#include <iomanip>
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"

#include <givaro/modular.h>

using namespace std;
using namespace FFLAS;
using namespace FFPACK;


template <class Field,  class RandIter>
bool test_square_fgesv (Field& F, FFLAS_SIDE side, string fileA, string fileB, size_t m, size_t k, size_t r, RandIter& G){

    typename Field::Element_ptr A, B, X, Acop, R=NULL;
    size_t lda, ldb, ldx, ldr, brows, bcols;
    if (side == FflasLeft){brows = m; bcols = k;}
    else {brows = k; bcols = m;}
    
    if (!fileA.empty()){
        ReadMatrix (fileA, F, m, m, A);
        lda = m;
    } else {
        lda = m+(rand() % 13);
        A = fflas_new (F, m,lda);
        RandomMatrixWithRankandRandomRPM (F, m, m, r, A, lda, G);
    }
    if (!fileB.empty()){
        ReadMatrix (fileB, F, brows, bcols, B);
        ldb = ldx = ldr = bcols;
    } else {
        ldb = ldx = ldr = bcols+(rand() % 13);
        B = fflas_new (F, brows,ldb);
        if (r==m)
            RandomMatrix(F, brows, bcols, B, ldb, G);
        else{
            if (side == FflasLeft){
                typename Field::Element_ptr S=fflas_new(F,m,k);
                RandomMatrix(F, m, k, S, k, G);
                fgemm (F, FflasNoTrans, FflasNoTrans, m, k, m, F.one, A, lda, S, k, F.zero, B, ldb);
                fflas_delete(S);
            } else {
                typename Field::Element_ptr S=fflas_new(F,k,m);
                RandomMatrix(F, k, m, S, m, G);
                fgemm (F, FflasNoTrans, FflasNoTrans, k, m, m, F.one, S, m, A, lda, F.zero, B, ldb);
                fflas_delete(S);
            }
        }
    }
    Acop = fflas_new(F, m, lda);
    X = fflas_new(F, brows, ldx);
    R = fflas_new(F, brows, ldr);
    fassign (F, brows, bcols, B,ldb, X, ldx);
    fassign (F, m, m, A, lda, Acop, lda);

    int info=0;
    bool ok = true;
    fgesv(F, side, brows, bcols, A, lda, X, ldx, &info);

    if (side == FflasLeft)
        fgemm (F, FflasNoTrans, FflasNoTrans, m, k, m, F.one, Acop, lda, X, ldx, F.zero, R, ldr);
    else
        fgemm (F, FflasNoTrans, FflasNoTrans, k, m, m, F.one, X, ldx, Acop, lda, F.zero, R, ldr);

    ok = ok && fequal (F, brows, bcols, R, ldr, B, ldb) && (info == 0);

    if (!ok){
        if (side == FflasLeft) std::cerr<<"ERROR A X != B"<<std::endl;
        else std::cerr<<"ERROR X A != B"<<std::endl;
        WriteMatrix(std::cerr<<"A ="<<std::endl,F,m,m,Acop,lda);
        WriteMatrix(std::cerr<<"X ="<<std::endl,F,brows,bcols,X,ldx);
        WriteMatrix(std::cerr<<"B ="<<std::endl,F,brows,bcols,B,ldb);
        WriteMatrix(std::cerr<<"AX ="<<std::endl,F,brows,bcols,R,ldr);
    }
    cout<<".";

    fflas_delete(A,B,X,Acop,R);
    return ok;
}

template <class Field,  class RandIter>
bool test_rect_fgesv (Field& F, FFLAS_SIDE side, string fileA, string fileB, size_t m, size_t n, size_t k, size_t r, RandIter& G){

    typename Field::Element_ptr A, B, X, Acop, R=NULL;
    size_t lda, ldb, ldx, ldr, brows, bcols, xrows, xcols, nbeq;
       
    if (!fileA.empty()){
        ReadMatrix (fileA, F, m, n, A);
        lda = n;
    } else {
        lda = n+(rand() % 13);
        A = fflas_new (F, m,lda);
        RandomMatrixWithRankandRandomRPM (F, m, n, r, A, lda, G);
    }
    if (!fileB.empty()){
        ReadMatrix (fileB, F, brows, bcols, B);
        ldb = ldr = bcols;
        if (side == FflasLeft) {ldx = xcols = bcols; xrows = n; nbeq = m;}
        else {ldx = xcols = m; xrows = brows; nbeq = n;}
    } else {
        if (side == FflasLeft) {brows = m; bcols = k;ldx = xcols = k; xrows = n; nbeq = m;}
        else {brows = k; bcols = n;ldx = xcols = m; xrows = k; nbeq = n;}
        ldb = ldr = bcols+(rand() % 13);
        B = fflas_new (F, brows,ldb);
        if (r == nbeq) // any B matrix makes a consistent system
            RandomMatrix(F, brows, bcols, B, ldb, G);
        else{ // need to sample B from the row/col-space of A
            if (side == FflasLeft){
                typename Field::Element_ptr S = fflas_new(F, n, k);
                RandomMatrix(F, n, k, S, k, G);
                fgemm (F, FflasNoTrans, FflasNoTrans, m, k, n, F.one, A, lda, S, k, F.zero, B, ldb);
                fflas_delete(S);
            } else {
                typename Field::Element_ptr S = fflas_new(F, k, m);
                RandomMatrix(F, k, m, S, m, G);
                fgemm (F, FflasNoTrans, FflasNoTrans, k, n, m, F.one, S, m, A, lda, F.zero, B, ldb);
                fflas_delete(S);
            }
        }
    }
  
    Acop = fflas_new(F, m, lda);
    X = fflas_new(F, xrows, ldx);
    R = fflas_new(F, brows, ldr);
    fassign (F, m, n, A, lda, Acop, lda);

    int info=0;
    bool ok = true;
    fgesv(F, side, m, n, k, A, lda, X, ldx, B, ldb, &info);

    if (side == FflasLeft)
        fgemm (F, FflasNoTrans, FflasNoTrans, m, k, n, F.one, Acop, lda, X, ldx, F.zero, R, ldr);
    else
        fgemm (F, FflasNoTrans, FflasNoTrans, k, n, m, F.one, X, ldx, Acop, lda, F.zero, R, ldr);

    ok = ok && fequal(F, brows, bcols, R, ldr, B, ldb) && (info == 0);

    if (!ok){
        if (side == FflasLeft) std::cerr<<"ERROR A X != B"<<std::endl;
        else std::cerr<<"ERROR X A != B"<<std::endl;
        WriteMatrix(std::cerr<<"A ="<<std::endl,F,m,n,Acop,lda);
        WriteMatrix(std::cerr<<"X ="<<std::endl,F,xrows,xcols,X,ldx);
        WriteMatrix(std::cerr<<"B ="<<std::endl,F,brows,bcols,B,ldb);
        WriteMatrix(std::cerr<<"AX ="<<std::endl,F,brows,bcols,R,ldr);
    }

    cout<<".";

    fflas_delete(A,B,X,Acop,R);
    return ok;
}
template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t m, size_t n, size_t k, size_t r, size_t iters, string fileA, string fileB,  uint64_t& seed){
    bool ok = true ;
    int nbit=(int)iters;

    while (ok &&  nbit){
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
        if (F==nullptr)
            return true;
        typename Field::RandIter G(*F,seed++);

        ostringstream oss;
        F->write(oss);

        cout.fill('.');
        cout<<"Checking ";
        cout.width(40);
        cout<<oss.str();
        cout<<" ...";

            // Full rank instances
        ok = ok && test_square_fgesv (*F, FflasLeft, fileA, fileB, m, k, m, G);
        ok = ok && test_square_fgesv (*F, FflasRight, fileA, fileB, m, k, m, G);
        ok = ok && test_rect_fgesv (*F, FflasLeft, fileA, fileB, m, n, k, std::min(m,n), G);
        ok = ok && test_rect_fgesv (*F, FflasRight, fileA, fileB, m, n, k, std::min(m,n), G);
            // Rank defficient instances
        ok = ok && test_square_fgesv (*F, FflasLeft, fileA, fileB, m, k, r, G);
        ok = ok && test_square_fgesv (*F, FflasRight, fileA, fileB, m, k, r, G);
        ok = ok && test_rect_fgesv (*F, FflasLeft, fileA, fileB, m, n, k, r, G);
        ok = ok && test_rect_fgesv (*F, FflasRight, fileA, fileB, m, n, k, r, G);


            // Randomized dimensions
        size_t MM = (rand() % m)+50;
        size_t NN = (rand() % n)+50;
        size_t KK = (rand() % k)+50;
        size_t RR = (rand() % std::min(NN,MM));
            // Full rank instances
        ok = ok && test_square_fgesv (*F, FflasLeft, fileA, fileB, MM, KK, MM, G);
        ok = ok && test_square_fgesv (*F, FflasRight, fileA, fileB, MM, KK, MM, G);
        ok = ok && test_rect_fgesv (*F, FflasLeft, fileA, fileB, MM, NN, KK, std::min(MM,NN), G);
        ok = ok && test_rect_fgesv (*F, FflasRight, fileA, fileB, MM, NN, KK, std::min(MM,NN), G);
            // Rank defficient instances
        ok = ok && test_square_fgesv (*F, FflasLeft, fileA, fileB, MM, KK, RR, G);
        ok = ok && test_square_fgesv (*F, FflasRight, fileA, fileB, MM, KK, RR, G);
        ok = ok && test_rect_fgesv (*F, FflasLeft, fileA, fileB, MM, NN, KK, RR, G);
        ok = ok && test_rect_fgesv (*F, FflasRight, fileA, fileB, MM, NN, KK, RR, G);

        delete F;

        nbit--;
        if (!ok)
            cout << "FAILED "<<endl;
        else
            cout << "PASSED "<<endl;
    }
    return ok;
}


int main(int argc, char** argv){
    cerr<<setprecision(20);

    Givaro::Integer q = -1;
    size_t b = 0;
    size_t m = 257;
    size_t n = 210;
    size_t k = 112;
    size_t r = 127;
    size_t iters = 4 ;
    bool loop=false;
    uint64_t seed = getSeed();
    string fileA = "";
    string fileB = "";
    Argument as[] = {
        { 'q', "-q Q", "Set the field cardinality.",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
        { 'm', "-m M", "Set the row dimension the system matrix A.", TYPE_INT , &m },
        { 'n', "-n N N", "Set the column dimension of the system matrix.", TYPE_INT , &n },
        { 'k', "-k K", "Set the dimension of the right/left-handside of the system.", TYPE_INT , &k },
        { 'r', "-r R", "Set the rank of the matrix.", TYPE_INT , &r },
        { 'i', "-i I", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        { 'f', "-f fileA", "Set input file for system matrix A", TYPE_STR, &fileA },
        { 'g', "-g fileB", "Set input file for right/left hand side matrix B", TYPE_STR, &fileB },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);

    srand(seed);
    cerr << "with seed = " << seed << endl;
    if (r>std::min(m,n)) r=std::min(m,n);
    bool ok=true;
    do{
        ok = ok && run_with_field<Givaro::Modular<float> >           (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::Modular<double> >          (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<float> >   (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<double> >  (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::Modular<int32_t> >         (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::Modular<int64_t> >         (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,m,n,k,r,iters,fileA, fileB,seed);
        ok = ok && run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:128),m/4+1,n/4+1,k/4+1,r/4+1,iters,fileA, fileB,seed);
    } while (loop && ok);

    if (!ok) cerr << "with seed = " << seed << endl;

    return !ok;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
