/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
 *
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
 *
 */

//--------------------------------------------------------------------------
//          Test for Checker_fgemm
//--------------------------------------------------------------------------

#define ENABLE_ALL_CHECKINGS 1

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"

using namespace Givaro;
using namespace FFLAS;
using namespace FFPACK;
template<class Field, class RandIter>
bool launch_MM_dispatch(const Field &F, const int mm, const int nn, const int kk, const typename Field::Element alpha,
                        const typename Field::Element beta, const size_t iters, RandIter& G)
{
    size_t m,n,k;
    size_t lda,ldb,ldc;
    //!@bug test for ldX equal
    //!@bug test for transpo
    //!@todo does nbw actually do nbw recursive calls and then call blas (check ?) ?
    size_t ld = 13 ;

    FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans ;
    FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans ;
    if (random()%2) ta = FFLAS::FflasTrans ;
    if (random()%2) tb = FFLAS::FflasTrans ;

    if (mm<0)
        m = 1+(size_t)random() % -mm;
    else m = mm;
    if (nn<0)
        n = 1+(size_t)random() % -nn;
    else n = nn;
    if (kk<0)
        k = 1+(size_t)random() % -kk;
    else k = kk;

    lda = std::max(k,m)+(size_t)random()%ld;
    ldb = std::max(n,k)+(size_t)random()%ld;
    ldc = n+(size_t)random()%ld;
#ifdef __FFLASFFPACK_DEBUG
    std::cerr <<"q = "<<F.characteristic()<<" m,k,n = "<<m<<", "<<k<<", "<<n<<" C := "
    <<alpha<<".A"<<((ta==FFLAS::FflasTrans)?"^T":"")
    <<" * B"<<((tb==FFLAS::FflasTrans)?"^T":"");
    if (!F.isZero(beta))
        std::cerr<<" + "<<beta<<" C";
#endif

    typename Field::Element_ptr A, B, C;
    C = FFLAS::fflas_new (F,m,ldc);
    FFLASFFPACK_check(ldc >= n);
    size_t Arows,Acols, Brows,Bcols;
    if (ta == FFLAS::FflasNoTrans){
        FFLASFFPACK_check(lda >= k);
        Arows = m; Acols = k;
    } else {
        FFLASFFPACK_check(lda >= m);
        Arows = k; Acols = m;
    }
    if (tb == FFLAS::FflasNoTrans){
        FFLASFFPACK_check(ldb >= n);
        Brows = k; Bcols = n;
    } else {
        FFLASFFPACK_check(ldb >= k);
        Brows = n; Bcols = k;
    }
    A = FFLAS::fflas_new (F, Arows, lda);
    FFLAS::fzero(F,Arows,lda,A,lda);
    B = FFLAS::fflas_new (F, Brows, ldb);
    FFLAS::fzero(F,Brows,ldb,B,ldb);
    for(size_t i = 0; i<iters;++i){
        RandomMatrix(F, Arows, Acols, A, lda, G);
        RandomMatrix(F, Brows, Bcols, B, ldb, G);
        RandomMatrix(F, m, n, C, ldc, G);

        FFLAS::Checker_fgemm<Field> checker(F,m,n,k,beta,C,ldc);
        try {
            FFLAS::fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
            checker.check(ta,tb,alpha,A,lda,B,ldb,C);
            //std::cout << "PASSED\n";
        } catch (FailureFgemmCheck &e) {
            std::cout << "FAILED\n";
            FFLAS::fflas_delete(A,B,C);
            return false;
        }
    }
    FFLAS::fflas_delete(A,B,C);
    return true;
}

template <class Field>
bool run_with_field (Givaro::Integer q, uint64_t b, int m, int n, int k, size_t iters, uint64_t seed){
    bool ok = true ;
    uint64_t local_seed = seed;
    int nbit=(int)iters;
    while (ok &&  nbit){
        typedef typename Field::Element Element ;
        // choose Field
        srand(local_seed);
        Field* F= chooseField<Field>(q,b,local_seed);
        if (F==nullptr)
            return true;

        std::ostringstream oss;
        F->write(oss);
        std::cout.fill('.');
        std::cout<<"Checking ";
        std::cout.width(50);
        std::cout<<oss.str();
        std::cout<<" ... ";

#ifdef __FFLASFFPACK_DEBUG
        F->write(std::cerr) << std::endl;
#endif
        typedef typename Field::Element  Element ;
        typename Field::RandIter R(*F,local_seed++);
        typename Field::NonZeroRandIter NZR(R);

        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->one,F->zero,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->zero,F->zero,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,F->zero,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->one ,F->one,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->zero,F->one,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,F->one,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->one ,F->mOne,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->zero,F->mOne,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,F->mOne,iters, R);

        Element alpha,beta ;
        NZR.random(alpha);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->one ,alpha,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->zero,alpha,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,F->mOne,alpha,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,alpha,F->one ,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,alpha,F->zero,iters, R);
        ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,alpha,F->mOne,iters, R);

        for (size_t j = 0 ; j < 3 ; ++j) {
            R.random(alpha);
            R.random(beta);
            ok = ok && launch_MM_dispatch<Field>(*F,m,n,k,alpha,beta,iters, R);
        }
        nbit--;
        if ( !ok )
            std::cout << "FAILED with seed = "<<local_seed-1<<std::endl;
        else
            std::cout << "PASSED with seed = "<<local_seed-1<<std::endl;
        delete F;
    }
    return ok;
}

int main(int argc, char** argv)
{
    std::cout<<std::setprecision(17);
    std::cerr<<std::setprecision(17);
    uint64_t seed = getSeed();
    size_t iters = 3 ;
    Givaro::Integer q = -1 ;
    uint64_t b = 0 ;
    int m = 50 ;
    int n = 50 ;
    int k = 50 ;
    bool loop = false;
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
        { 'm', "-m M", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |n|).",      TYPE_INT , &m },
        { 'n', "-n N", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |n|).",      TYPE_INT , &n },
        { 'k', "-k K", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |k|).",      TYPE_INT , &k },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-l Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments(argc,argv,as);
    bool ok = true;
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,m,n,k,iters, seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,m,n,k,iters, seed);
        ok = ok && run_with_field<Modular<float> >(q,b,m,n,k,iters, seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,m,n,k,iters, seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,m,n,k,iters, seed);
        ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,m,n,k,iters, seed);
        seed++;
    } while (loop && ok);

    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
