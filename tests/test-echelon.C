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
//          Test for the echelon factorisation
//--------------------------------------------------------------------------

//#define __LUDIVINE_CUTOFF 1

#define  __FFLASFFPACK_SEQUENTIAL
#define __FFLASFFPACK_GAUSSJORDAN_BASECASE 25
#define __FFLASFFPACK_PLUQ_THRESHOLD 25
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <iomanip>
#include <givaro/modular-balanced.h>
#include <givaro/udl.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_io.h"

#include <random>
#include <chrono>

using namespace std;
using namespace FFPACK;
using namespace FFLAS;
using Givaro::Modular;
using Givaro::ModularBalanced;

template<class Field, class RandIter>
bool
test_colechelon(Field &F, size_t m, size_t n, size_t r, size_t iters, FFPACK::FFPACK_LU_TAG LuTag, RandIter& G, bool par)
{
    typedef typename Field::Element Element ;
    Element * A = FFLAS::fflas_new (F,m,n);
    Element * B = FFLAS::fflas_new (F,m,n);
    Element * L = FFLAS::fflas_new (F,m,n);
    Element * U = FFLAS::fflas_new (F,n,n);
    Element * X = FFLAS::fflas_new (F,m,n);
    size_t lda = n; //!@todo check lda

    size_t *P = FFLAS::fflas_new<size_t>(n);
    size_t *Q = FFLAS::fflas_new<size_t>(m);
    size_t R = (size_t)-1;

    bool pass=true;

    for (size_t  l=0;l<iters;l++){
        R = (size_t)-1;
        RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);
        FFLAS::fassign(F,m,n,A,lda,B,lda);
        
        for (size_t j=0;j<n;j++) P[j]=0;
        for (size_t j=0;j<m;j++) Q[j]=0;
        if(!par)
            R = FFPACK::ColumnEchelonForm (F, m, n, A, n, P, Q, true, LuTag);
        else
            R = FFPACK::pColumnEchelonForm (F, m, n, A, n, P, Q, true, 0, LuTag);
        if (R != r) {pass = false; break;}

        FFPACK::getEchelonTransform (F, FFLAS::FflasLower, FFLAS::FflasUnit, m,n,R,P,Q,A,lda,U,n, LuTag);

        FFPACK::getEchelonForm (F, FFLAS::FflasLower, FFLAS::FflasUnit, m,n,R,Q,A,n,L,n,false, LuTag);

        // Testing if C is in col echelon form
        size_t nextpiv = 0;
        for (size_t j=0; j<R; ++j){
            size_t i=0;
            while ((i < m) && F.isZero (L[i*n+j])) i++;
            if (i==m) // zero column in the first R columns
                pass = false;
            if (i < nextpiv)  // not in echelon form
                pass = false;
            nextpiv = i+1;
        }
        pass = pass && FFLAS::fiszero (F, m, n-R, L+R, n);

        // Testing A U = L
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,n, 1.0, B, n, U, n, 0.0, X,n);

        pass = pass && FFLAS::fequal(F, m, n, L, n, X, n);

        if (!pass) {
            std::cerr<<"FAIL (column echelon LuTag="<<LuTag<<")"<<std::endl;
            FFLAS::WriteMatrix (std::cerr<<"A = "<<std::endl,F,m,n,B,lda);
            FFLAS::WriteMatrix (std::cerr<<"InplaceEchelon = "<<std::endl,F,m,n,A,lda);
            FFLAS::WritePermutation(std::cerr<<"P = [",P,n)<<std::endl;
            FFLAS::WritePermutation(std::cerr<<"Q = [",Q,m)<<std::endl;
            FFLAS::WriteMatrix(std::cerr<<"ColEchelon = "<<std::endl,F,m,n,L,n);
            FFLAS::WriteMatrix(std::cerr<<"Transform = "<<std::endl,F,n,n,U,n);
            FFLAS::WriteMatrix(std::cerr<<"B x X  = "<<std::endl,F,m,n,X,n);
            break;
        }
    }
    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( X);
    FFLAS::fflas_delete( B);
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);
    return pass;
}

template<class Field, class RandIter>
bool
test_rowechelon(Field &F, size_t m, size_t n, size_t r, size_t iters, FFPACK::FFPACK_LU_TAG LuTag, RandIter& G, bool par)
{
    typedef typename Field::Element Element ;
    Element * A = FFLAS::fflas_new (F,m,n);
    Element * B = FFLAS::fflas_new (F,m,n);
    Element * L = FFLAS::fflas_new (F,m,m);
    Element * U = FFLAS::fflas_new (F,m,n);
    Element * X = FFLAS::fflas_new (F,m,n);
    size_t lda = n; //!@todo check lda

    size_t *P = FFLAS::fflas_new<size_t>(m);
    size_t *Q = FFLAS::fflas_new<size_t>(n);
    size_t R = (size_t)-1;

    bool pass=true;

    for (size_t  l=0;l<iters;l++){
        R = (size_t)-1;
        RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);
        FFLAS::fassign(F,m,n,A,lda,B,lda);

        for (size_t j=0;j<m;j++) P[j]=0;
        for (size_t j=0;j<n;j++) Q[j]=0;

        if(!par)
            R = FFPACK::RowEchelonForm (F, m, n, A, n, P, Q, true, LuTag);
        else
            R = FFPACK::pRowEchelonForm (F, m, n, A, n, P, Q, true, 0, LuTag);

        if (R != r) {pass = false; break;}

        FFPACK::getEchelonTransform (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m,n,R,P,Q,A,lda,L,m, LuTag);

        FFPACK::getEchelonForm (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m,n,R,Q,A,n,U,n, false, LuTag);

        // Testing if U is in row echelon form
        size_t nextpiv = 0;
        for (size_t j=0; j<R; ++j){
            size_t i=0;
            while ((i < n) && F.isZero (U[i+j*n])) i++;
            if (i==n) // zero row in the first R columns
                pass = false;
            if (i < nextpiv)  // not in echelon form
                pass = false;
            nextpiv = i+1;
        }
        pass = pass && FFLAS::fiszero (F, m-R, n, U+R*n, n);

        // Testing A U = L
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,m, 1.0, L, m, B, n, 0.0, X,n);

        pass = pass && FFLAS::fequal(F, m, n, U, n, X, n);

        if (!pass) {
            std::cerr<<"FAIL (row echelon LuTag="<<LuTag<<")"<<std::endl;
            FFLAS::WriteMatrix (std::cerr<<"A = "<<std::endl,F,m,n,B,lda);
            FFLAS::WriteMatrix (std::cerr<<"InplaceEchelon = "<<std::endl,F,m,n,A,lda);
            FFLAS::WritePermutation(std::cerr<<"P = [",P,m)<<std::endl;
            FFLAS::WritePermutation(std::cerr<<"Q = [",Q,n)<<std::endl;
            FFLAS::WriteMatrix (std::cerr<<"RowEchelon = "<<std::endl,F,m,n,U,n);
            FFLAS::WriteMatrix (std::cerr<<"Transform = "<<std::endl,F,m,m,L,m);
            break;
        }
    }
    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( X);
    FFLAS::fflas_delete( B);
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);
    return pass;
}

template<class Field, class RandIter>
bool
test_redcolechelon(Field &F, size_t m, size_t n, size_t r, size_t iters, FFPACK::FFPACK_LU_TAG LuTag, RandIter& G, bool par)
{
    typedef typename Field::Element Element ;
    Element * A = FFLAS::fflas_new (F,m,n);
    Element * B = FFLAS::fflas_new (F,m,n);
    Element * L = FFLAS::fflas_new (F,m,n);
    Element * U = FFLAS::fflas_new (F,n,n);
    Element * X = FFLAS::fflas_new (F,m,n);
    size_t lda = n; //!@todo check lda

    size_t *P = FFLAS::fflas_new<size_t>(n);
    size_t *Q = FFLAS::fflas_new<size_t>(m);
    size_t R = (size_t)-1;
    bool pass=true;

    for (size_t  l=0;l<iters;l++){
        R = (size_t)-1;
        RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);
        FFLAS::fassign(F,m,n,A,lda,B,lda);

        for (size_t j=0;j<n;j++) P[j]=0;
        for (size_t j=0;j<m;j++) Q[j]=0;

        if(!par)
            R = FFPACK::ReducedColumnEchelonForm (F, m, n, A, n, P, Q, true, LuTag);
        else
            R = FFPACK::pReducedColumnEchelonForm (F, m, n, A, n, P, Q, true, 0, LuTag);

        if (R != r) {pass = false; break;}

        FFPACK::getReducedEchelonTransform (F, FFLAS::FflasLower, m,n,R,P,Q,A,lda,U,n, LuTag);

        FFPACK::getReducedEchelonForm (F, FFLAS::FflasLower, m,n,R,Q,A,n,L,n, false, LuTag);

        // Testing if C is in reduced col echelon form
        size_t nextpiv = 0;
        for (size_t j=0; j<R; ++j){
            size_t i=0;
            while ((i < m) && F.isZero (L[i*n+j])) i++;
            if (i==m) // zero column in the first R columns
                pass = false;
            if (i < nextpiv)  // not in echelon form
                pass = false;
            if (j) // is pivot row reduced
                pass = pass && FFLAS::fiszero(F, j-1, L + i*n, 1);
            pass = pass && F.isOne(L[j+i*n]);
            nextpiv = i+1;
        }
        pass = pass && FFLAS::fiszero (F, m, n-R, L+R, n);
        // Testing A U = L
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,n, 1.0, B, n, U, n, 0.0, X,n);

        pass = pass && FFLAS::fequal(F, m, n, L, n, X, n);

        if (!pass) {
            std::cerr<<"FAIL (reduced column echelon LuTag="<<LuTag<<")"<<std::endl;
            break;
        }
    }
    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( X);
    FFLAS::fflas_delete( B);
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);
    return pass;
}
template<class Field, class RandIter>
bool
test_redrowechelon(Field &F, size_t m, size_t n, size_t r, size_t iters, FFPACK::FFPACK_LU_TAG LuTag, RandIter& G, bool par)
{
    typedef typename Field::Element Element ;
    Element * A = FFLAS::fflas_new (F,m,n);
    Element * B = FFLAS::fflas_new (F,m,n);
    Element * L = FFLAS::fflas_new (F,m,m);
    Element * U = FFLAS::fflas_new (F,m,n);
    Element * X = FFLAS::fflas_new (F,m,n);
    size_t lda = n; //!@todo check lda

    size_t *P = FFLAS::fflas_new<size_t>(m);
    size_t *Q = FFLAS::fflas_new<size_t>(n);
    size_t R = (size_t)-1;


    bool pass=true;

    for (size_t  l=0;l<iters;l++){
        R = (size_t)-1;

        RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);

        FFLAS::fassign(F,m,n,A,lda,B,lda);

        for (size_t j=0;j<m;j++) P[j]=0;
        for (size_t j=0;j<n;j++) Q[j]=0;

        if(!par)
            R = FFPACK::ReducedRowEchelonForm (F, m, n, A, n, P, Q, true, LuTag);
        else
            R = FFPACK::pReducedRowEchelonForm (F, m, n, A, n, P, Q, true, 0, LuTag);


        if (R != r) {pass = false; break;}

        FFPACK::getReducedEchelonTransform (F, FFLAS::FflasUpper, m,n,R,P,Q,A,lda,L,m, LuTag);

        FFPACK::getReducedEchelonForm (F, FFLAS::FflasUpper, m,n,R,Q,A,n,U,n, false, LuTag);

        // Testing if U is in row echelon form
        size_t nextpiv = 0;
        for (size_t j=0; j<R; ++j){
            size_t i=0;
            while ((i < n) && F.isZero (U[i+j*n])) i++;
            if (i==n) // zero row in the first R rows
                pass = false;
            if (i < nextpiv)  // not in echelon form
                pass = false;
            if (j) // is pivot row reduced
                pass = pass && FFLAS::fiszero(F, j-1, U + i, n);
            pass = pass && F.isOne(U[j*n+i]);
            nextpiv = i+1;
        }
        pass = pass && FFLAS::fiszero (F, m-R, n, U+R*n, n);

        // Testing A U = L
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,m, 1.0, L, m, B, n, 0.0, X,n);

        pass = pass && FFLAS::fequal(F, m, n, U, n, X, n);

        if (!pass) {
            std::cerr<<"FAIL (reduced row echelon LuTag="<<LuTag<<")"<<std::endl;
            FFLAS::WriteMatrix (std::cerr<<"B = "<<std::endl,F,m,n,B,lda);
            FFLAS::WriteMatrix (std::cerr<<"RedRowEchelon = "<<std::endl,F,m,n,U,n);
            FFLAS::WriteMatrix (std::cerr<<"X x B  = "<<std::endl,F,m,n,X,n);
            FFLAS::WriteMatrix (std::cerr<<"Transform = "<<std::endl,F,m,m,L,m);
            FFLAS::WriteMatrix (std::cerr<<"InplaceEchelon = "<<std::endl,F,m,n,A,lda);
            std::cerr<<"R = "<<R<<std::endl;
            std::cerr<<"P = [";	for (size_t i=0; i<m; ++i) std::cerr<<P[i]<<", ";std::cerr<<"]\n";
            std::cerr<<"Q = [";	for (size_t i=0; i<n; ++i) std::cerr<<Q[i]<<", ";std::cerr<<"]\n";
            break;
        }
    }
    FFLAS::fflas_delete( U);
    FFLAS::fflas_delete( L);
    FFLAS::fflas_delete( X);
    FFLAS::fflas_delete( B);
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( P);
    FFLAS::fflas_delete( Q);
    return pass;
}

template <class Field>
bool run_with_field (Givaro::Integer q, uint64_t b, size_t m, size_t n, size_t r, size_t iters, uint64_t seed){
    bool ok = true ;

    int nbit=(int)iters;

    while (ok &&  nbit){
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
        if (F==nullptr)
            return true;
        typename Field::RandIter G(*F,seed++);
        std::ostringstream oss;
        F->write(oss);
        std::cout.fill('.');
        std::cout<<"Checking ";
        std::cout.width(40);
        std::cout<<oss.str();
        std::cout<<" .";

        ok = ok && test_colechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, false);
        std::cout<<".";
        ok = ok && test_colechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, false);
        std::cout<<".";
        ok = ok && test_redcolechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, false);
        std::cout<<".";
        ok = ok && test_redcolechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, false);
        std::cout<<".";
        ok = ok && test_rowechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, false);
        std::cout<<".";
        ok = ok && test_rowechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, false);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, false);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, false);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackGaussJordanSlab, G, false);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackGaussJordanTile, G, false);
        std::cout<<".";

        ok = ok && test_colechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, true);
        std::cout<<".";
        ok = ok && test_colechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, true);
        std::cout<<".";
        ok = ok && test_redcolechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, true);
        std::cout<<".";
        ok = ok && test_redcolechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, true);
        std::cout<<".";
        ok = ok && test_rowechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, true);
        std::cout<<".";
        ok = ok && test_rowechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, true);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackSlabRecursive, G, true);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackTileRecursive, G, true);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackGaussJordanSlab, G, true);
        std::cout<<".";
        ok = ok && test_redrowechelon(*F,m,n,r,iters, FFPACK::FfpackGaussJordanTile, G, true);
        std::cout<<".";

        nbit--;
        if ( !ok )
            std::cout << "FAILED with seed = "<<seed<<std::endl;
        else
            std::cout << "PASSED "<<std::endl;
        delete F;
    }
    return ok;
}

int main(int argc, char** argv){
    std::cerr<<std::setprecision(20);

    Givaro::Integer q = -1;
    size_t b = 0;
    size_t m = 124;
    size_t n = 64;
    size_t r = 34;
    size_t iters = 3 ;
    bool loop = false;
    uint64_t seed =  getSeed();
    static Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic.",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.", TYPE_INT , &b },
        { 'n', "-n N", "Set the number of cols in the matrix.", TYPE_INT , &n },
        { 'm', "-m N", "Set the number of rows in the matrix.", TYPE_INT , &m },
        { 'r', "-r r", "Set the rank of the matrix."          , TYPE_INT , &r },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-l Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        // { 'f', "-f file", "Set input file", TYPE_STR, &file },
        END_OF_ARGUMENTS
    };
    FFLAS::parseArguments(argc,argv,as);
    r = std::min(r, std::min(m,n));

    srand(seed);

    bool ok = true;
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,m,n,r,iters,seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,m,n,r,iters,seed);
        ok = ok && run_with_field<Modular<float> >(q,b,m,n,r,iters,seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,m,n,r,iters,seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,m,n,r,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,m,n,r,iters,seed);
        ok = ok && run_with_field<Modular<int64_t> >(q,b,m,n,r,iters,seed);
        //ok = ok && run_with_field<Modular<RecInt::rint<7> > >(q,b,m,n,r,iters,seed); // BUG: not available yet (missing division in the field
        ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,m,n,r,iters,seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:128_ui64),m/8+1,n/8+1,r/8+1,iters,seed);

    } while (loop && ok);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
