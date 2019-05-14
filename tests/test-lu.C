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


//-------------------------------------------------------------------------
//      Test suite for the Gaussian elimination routines: LUdivine and PLUQ
//-------------------------------------------------------------------------

// #define MONOTONIC_CYLCES
// #define MONOTONIC_MOREPIVOTS
// #define MONOTONIC_FEWPIVOTS
#ifdef MONOTONIC_CYLCES
#define MONOTONIC_APPLYP
#endif
#ifdef MONOTONIC_MOREPIVOTS
#define MONOTONIC_APPLYP
#endif
#ifdef MONOTONIC_FEWPIVOTS
#define MONOTONIC_APPLYP
#endif

#define BASECASE_K 37 // Forcing a lower base case to be able to test a few recursive steps with smallish dimensions


#define  __FFLASFFPACK_SEQUENTIAL
#define __LUDIVINE_CUTOFF 1
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-balanced.h>
#include <iostream>
#include <iomanip>
Givaro::Timer tperm, tgemm, tBC, ttrsm,trest,timtot;
size_t mvcnt = 0;
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/test-utils.h"

#include "fflas-ffpack/utils/args-parser.h"

#include <random>

using namespace std;
using namespace FFPACK;
using namespace FFLAS;

/*! Tests the LUdivine routine.
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in U
 * @tparam Trans
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param m rows
 * @param n cols
 * @param lda leading dim of A
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS::FFLAS_DIAG diag, FFLAS_TRANSPOSE trans>
bool test_LUdivine(const Field & F,
                   typename Field::ConstElement_ptr A, size_t lda,
                   size_t r, size_t m, size_t n)
{
    bool fail = false;
    typename Field::Element_ptr B = fflas_new(F,m,lda) ;
    fassign(F,m,n,A,lda,B,lda);

    size_t maxP, maxQ ;

    if (trans == FflasTrans){
        maxP = m;
        maxQ = n;
    }
    else{ // trans == FflasNoTrans
        maxP = n;
        maxQ = m;
    }

    size_t * P = fflas_new<size_t>(maxP) ;
    size_t * Q = fflas_new<size_t>(maxQ) ;

    size_t R = LUdivine (F, diag, trans, m, n, B, lda, P, Q);

    if (R != r) {
        std::cout << "rank is wrong (expecting " << r << " but got " << R << ")" << std::endl;
        fflas_delete( B );
        fflas_delete( P );
        fflas_delete( Q );
        return fail = true;
    }

    /*  Build L,U */
    typename Field::Element_ptr L = fflas_new(F, m, R);
    typename Field::Element_ptr U = fflas_new(F, R, n);

    if (trans == FflasNoTrans){
        getTriangular (F, FflasUpper, diag, m, n, R, B, lda, U, n, true);
        getEchelonForm (F, FflasLower, (diag==FflasNonUnit)?FflasUnit:FflasNonUnit, m, n, R, Q, B, lda, L, R, true);
        applyP (F, FflasRight, FflasNoTrans, R, 0, R, U, n, P);
    } else {
        getTriangular (F, FflasLower, diag, m, n, R, B, lda, L, R, true);
        getEchelonForm (F, FflasUpper, (diag==FflasNonUnit)?FflasUnit:FflasNonUnit, m, n, R, Q, B, lda, U, n, true);
        applyP (F, FflasLeft, FflasTrans, R, 0, R, L, R, P);
    }
    fgemm (F, FflasNoTrans, FflasNoTrans, m, n, R, 1.0, L, R, U, n, 0.0, B, lda);

    fail |= !fequal(F, m, n, A, lda, B, lda);

    if (fail){
        FFLAS::WriteMatrix(cerr<<"A = "<<endl,F,m,n,A,lda);
        FFLAS::WriteMatrix(cerr<<"LU = "<<endl,F,m,n,B,lda);
        FFLAS::WriteMatrix(cerr<<"L = "<<endl,F,m,R,L,R);
        FFLAS::WriteMatrix(cerr<<"U = "<<endl,F,R,n,U,n);
    }

    fflas_delete( P);
    fflas_delete( L);
    fflas_delete( U);
    fflas_delete( Q);
    fflas_delete( B);
    return fail;
}

/*! Verifies that B = PLUQ where A stores [L\U]
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in U
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param m rows
 * @param n cols
 * @param lda leading dim of A
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS_DIAG diag>
bool verifPLUQ (const Field & F, typename Field::ConstElement_ptr A, size_t lda,
                typename Field::Element_ptr PLUQ, size_t ldpluq,
                size_t * P, size_t * Q, size_t m, size_t n, size_t R)
{


    typename Field::Element_ptr X = fflas_new (F, m, n);
    typename Field::Element_ptr L = fflas_new (F, m, R);
    typename Field::Element_ptr	U = fflas_new (F, R, n);
    fzero(F, m, R, L, R);
    fzero(F, R, n, U, n);

    // FFLAS::WriteMatrix(std::cerr<<"PLUQ = "<<std::endl,F,m,n,PLUQ,ldpluq);
    getTriangular(F, FflasUpper, diag, m,n,R, PLUQ, ldpluq, U, n, true);
    getTriangular(F, FflasLower, (diag==FflasNonUnit)?FflasUnit:FflasNonUnit,
                  m,n,R, PLUQ, ldpluq, L, R, true);
    applyP( F, FflasLeft, FflasTrans, R,0,m, L, R, P);
    applyP (F, FflasRight, FflasNoTrans, R,0,n, U, n, Q);
    fgemm (F, FflasNoTrans, FflasNoTrans, m,n,R, F.one, L,R, U,n, F.zero, X,n);

    // FFLAS::WritePermutation(std::cerr<<"P = ",P,m);
    // FFLAS::WritePermutation(std::cerr<<"Q = ",Q,n);
    // FFLAS::WriteMatrix(std::cerr<<"L = "<<std::endl,F,m,R,L,R);
    // FFLAS::WriteMatrix(std::cerr<<"U = "<<std::endl,F,R,n,U,n);


    bool fail = false;
    for(size_t i=0; i<m; ++i)
        for (size_t j=0; j<n; ++j)
            if (!F.areEqual (*(A+i*lda+j), *(X+i*n+j))){
                std::cerr << std::endl<<" A ["<<i<<","<<j<<"] = " << (*(A+i*lda+j))
                << " PLUQ ["<<i<<","<<j<<"] = " << (*(X+i*n+j));
                fail=true;
            }
    //FFLAS::WriteMatrix(std::cerr<<"X = "<<std::endl,F, m,n,X,n);
    if (fail)
        std::cerr << std::endl;

    fflas_delete( U);
    fflas_delete( L);
    fflas_delete( X);
    return fail;
}
/*! Tests the LUdivine routine.
 * @tparam Field Field
 * @tparam Diag  Unit diagonal in U
 * @tparam Trans
 * @param F field
 * @param A Matrix (preallocated)
 * @param r rank of A
 * @param m rows
 * @param n cols
 * @param lda leading dim of A
 * @return 0 iff correct, 1 otherwise
 */
template<class Field, FFLAS_DIAG diag, class RandIter>
bool test_pluq (const Field & F,
                typename Field::ConstElement_ptr A,
                size_t r, size_t m, size_t n, size_t lda, RandIter& G)
{
    bool fail = false;
    typedef typename Field::Element_ptr Element_ptr ;
    Element_ptr B = fflas_new(F,m,lda) ;
    fassign(F,m,n,A,lda,B,lda);

    size_t * P = fflas_new<size_t> (m);
    size_t * Q = fflas_new<size_t> (n);

    ForceCheck_PLUQ<Field> checker (G,m,n,B,lda);

    size_t R = PLUQ (F, diag, m, n, B, lda, P, Q);

    try {
        checker.check(B,lda,diag,R,P,Q);
    } catch(FailurePLUQCheck &e) {
        std::cout << m << 'x' << n << " pluq verification failed!\n";
    }

    if (R != r) {
        std::cout << "rank is wrong (expected " << r << " but got " << R << ")" << std::endl;
        fflas_delete (B);
        fflas_delete (P);
        fflas_delete (Q);
        return fail = true;
    }
    fail |=  verifPLUQ<Field,diag> (F,A, lda, B, lda, P, Q, m, n, r);
    fflas_delete (B);
    fflas_delete(P);
    fflas_delete(Q);
    return fail;
}

template<class Field, FFLAS_DIAG diag, FFLAS_TRANSPOSE trans, class RandIter>
bool launch_test(const Field & F,
                 size_t r,
                 size_t m, size_t n, RandIter& G)
{
    //typedef typename Field::Element Element ;
    typedef typename Field::Element_ptr Element_ptr ;
    bool fail = false ;
    { /*  user given and lda bigger */
        size_t lda = n+10 ;
        Element_ptr A = fflas_new (F, m, lda);
        RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);
        fail |= test_LUdivine<Field,diag,trans>(F,A,lda,r,m,n);
        RandomMatrixWithRankandRandomRPM(F,m,n,r,A,lda,G);
        fail |= test_pluq<Field,diag>(F,A,r,m,n,lda,G);
        if (fail) std::cout << "failed at big lda" << std::endl;
        fflas_delete( A );
    }
    { /*  user given and lda bigger. Rank is max */
        size_t lda = n+10 ;
        size_t R = std::min(m,n);
        Element_ptr A = fflas_new (F, m, lda);
        RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
        fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,m,n);
        RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
        fail |= test_pluq<Field,diag>(F,A,R,m,n,lda,G);
        if (fail) std::cout << "failed at big lda max rank" << std::endl;
        fflas_delete( A );
    }
    { /*  user given and lda bigger. Rank is min */
        size_t lda = n+10 ;
        size_t R = 0;
        Element_ptr A = fflas_new (F, m, lda);
        RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
        fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,m,n);
        RandomMatrixWithRankandRandomRPM(F,m,n,R,A,lda,G);
        fail |= test_pluq<Field,diag>(F,A,R,m,n,lda,G);
        if (fail) std::cout << "failed at big lda, rank 0" << std::endl;
        fflas_delete( A );
    }
    { /*  square  */
        size_t M = std::max(m,n);
        size_t N = M ;
        size_t R = M/2 ;
        size_t lda = N+10 ;
        Element_ptr A = fflas_new (F, M, lda);
        RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
        fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
        RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
        fail |= test_pluq<Field,diag>(F,A,R,M,N,lda,G);
        if (fail) std::cout << "failed at square" << std::endl;
        fflas_delete( A );
    }
    { /*  wide  */
        size_t M = std::max(m,n);
        size_t N = 2*M ;
        size_t R = 3*M/4 ;
        size_t lda = N+5 ;
        Element_ptr A = fflas_new (F, M, lda);
        RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
        fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
        RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
        fail |= test_pluq<Field,diag>(F,A,R,M,N,lda,G);
        if (fail) std::cout << "failed at wide" << std::endl;
        fflas_delete( A );
    }
    { /*  narrow  */
        size_t M = std::max(m,n);
        size_t N = M/2 ;
        size_t R = 3*M/8 ;
        size_t lda = N+5 ;
        Element_ptr A = fflas_new (F, M, lda);
        RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
        fail |= test_LUdivine<Field,diag,trans>(F,A,lda,R,M,N);
        RandomMatrixWithRankandRandomRPM(F,M,N,R,A,lda,G);
        fail |= test_pluq<Field,diag>(F,A,R,M,N,lda,G);
        if (fail) std::cout << "failed at narrow" << std::endl;
        fflas_delete( A );
    }
    return !fail;
}

template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t m, size_t n, size_t r, size_t iters, uint64_t seed){
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
        std::cout<<" ... ";


        ok = ok && launch_test<Field,FflasUnit,FflasNoTrans>    (*F,r,m,n,G);
        ok = ok && launch_test<Field,FflasUnit,FflasTrans>      (*F,r,m,n,G);
        ok = ok && launch_test<Field,FflasNonUnit,FflasNoTrans> (*F,r,m,n,G);
        ok = ok && launch_test<Field,FflasNonUnit,FflasTrans>   (*F,r,m,n,G);

#if 0 /*  may be bogus */
        ok = ok && launch_test_append<Field,FflasUnit,FflasNoTrans>   (*F,r,m,n,G);
        ok = ok && launch_test_append<Field,FflasNonUnit,FflasNoTrans>(*F,r,m,n,G);
        ok = ok && launch_test_append<Field,FflasUnit,FflasTrans>     (*F,r,m,n,G);
        ok = ok && launch_test_append<Field,FflasNonUnit,FflasTrans>  (*F,r,m,n,G);
#endif
        nbit--;
        if ( !ok )
            //std::cout << "\033[1;31mFAILED\033[0m "<<std::endl;
            std::cout << "FAILED "<<std::endl;
        else
            //std::cout << "\033[1;32mPASSED\033[0m "<<std::endl;
            std::cout << "PASSED "<<std::endl;
        delete F;
    }
    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(20);
    Givaro::Integer q=-1;
    size_t b=0;
    size_t m=90;
    size_t n=93;
    size_t r=50;
    size_t iters=3;
    bool loop=false;
    uint64_t seed = getSeed();

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.",  TYPE_INT , &b },
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

    bool ok=true;
    do{
        ok = ok &&run_with_field<Givaro::Modular<float> >           (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<double> >          (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<float> >   (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<double> >  (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int32_t> >         (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int64_t> >         (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,m,n,r,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,5,m/6,n/6,r/6,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:512),m/6,n/6,r/6,iters,seed);
    } while (loop && ok);

    return !ok;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
