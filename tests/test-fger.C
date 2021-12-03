/*
 * Copyright (C) FFLAS-FFPACK
 * Written by JG Dumas
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
//                        Test for fger : 1 computation
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

// #define __FFLASFFPACK_DEBUG
#define TIME 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iomanip>
#include <iostream>
#include <givaro/modular-integral.h>
#include <givaro/modular-balanced.h>
#include <givaro/givintprime.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"

#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_io.h"

using namespace std;
using namespace FFPACK;
using namespace FFLAS;
using Givaro::Modular;
using Givaro::ModularBalanced;

// checks that D = alpha . x . y^T + C

// WARNING
template<class Field>
bool check_fger(const Field                   & F,
                const typename Field::Element_ptr  Cd, // c0
                const size_t                    m,
                const size_t                    n,
                const typename Field::Element & alpha,
                const typename Field::Element_ptr  x,
                const size_t                    incx,
                const typename Field::Element_ptr  y,
                const size_t                    incy,
                const typename Field::Element_ptr  C, // res
                const size_t                    ldc
               ) {
    bool wrong = false;

    typedef typename Field::Element Element;
    typedef typename Field::Element_ptr Element_ptr;

    // 	std::cerr << "with(LinearAlgebra):" << std::endl;
    //         WriteMatrix(std::cerr <<"X:=",F, m, 1, x,incx, FflasMaple) << ';' << std::endl;
    //         WriteMatrix(std::cerr <<"Y:=Transpose(", F, n, 1, y, incy, FflasMaple) << ");" << std::endl;
    //         WriteMatrix(std::cerr <<"A:=",F, m, n, Cd, ldc, FflasMaple) << ';' << std::endl;
    //         F.write(std::cerr << "a:=", alpha) << ';' << std::endl;
    // 	std::cerr << "q:=" << F.characteristic() << ';' << std::endl;
    Element_ptr D  = fflas_new (F,m,n);
    fassign(F,m,n,Cd,n,D,n);
    for(size_t i=0; i<m; ++i) {
        Element tmp; F.init(tmp);
        F.mul(tmp, alpha, *(x+i*incx) );
        for(size_t j=0; j<n; j+=incy) {
            F.axpyin(*(D+i*n+j), tmp, *(y+j) );
            if ( !F.areEqual( *(D+i*n+j), *(C+i*ldc+j) ) ) {
                wrong = true;
            }
        }
    }
    //     WriteMatrix(std::cerr <<"d:=",F, m, n, D, n, FflasMaple) << ';' << std::endl;
    // 	F.write(std::cerr, alpha) << "*X.Y+A,d;";
    // 	F.write(std::cerr, alpha) << "*X.Y+A-d mod q;" << std::endl;
    if ( wrong ){
        size_t ici = 20 ;
        std::cout<<"FAIL"<<std::endl;
        std::cout << "a   :" << alpha<<std::endl;
        std::cout << "m   :" << m   << ", n   : " <<  n  << std::endl;
        std::cout << "incx :" << incx << ", incy : " << incy << ", ldC : " << ldc << std::endl;
        for (size_t i=0; i<m && ici; ++i){
            for (size_t j =0; j<n && ici; ++j)
                if (!F.areEqual( *(C+i*ldc+j), *(D+i*n+j) ) ) {
                    std::cout<<"Error C["<<i<<","<<j<<"]="
                    <<(*(C+i*ldc+j))<<" D["<<i<<","<<j<<"]="
                    <<(*(D+i*n+j))<<std::endl;
                    ici--;
                }
        }
        if (m<80 && n<80) {
            for (size_t i=0; i<m ; ++i){
                for (size_t j =0; j<n ; ++j) {
                    if ( !F.areEqual( *(C+i*ldc+j), *(D+i*n+j) ) )
                        std::cout << 'X' ;
                    else
                        std::cout << '.' ;
                }
                std::cout << std::endl;
            }
        }
    }
    fflas_delete (D);

    return !wrong ;
}


template<class Field, class RandIter>
bool launch_fger(const Field & F,
                 const size_t   m,
                 const size_t   n,
                 const typename Field::Element alpha,
                 const size_t ldc,
                 const size_t inca,
                 const size_t incb,
                 size_t iters,
                 RandIter& G)
{
    bool ok = true;

    typedef typename Field::Element_ptr Element_ptr;
    Element_ptr A ;
    FFLASFFPACK_check(inca >= 1);
    Element_ptr B ;
    FFLASFFPACK_check(incb >= 1);
    Element_ptr C = fflas_new (F,m,ldc);
    FFLASFFPACK_check(ldc >= n);
    fzero(F,m,n,C,ldc);
    Element_ptr D = fflas_new (F, m, n);
    for(size_t i = 0;i<iters;++i){
        A = fflas_new (F, m, inca);
        RandomMatrix(F, m, inca, A, inca, G);
        B = fflas_new (F, n, incb);
        RandomMatrix(F, n, incb, B, incb, G);
        RandomMatrix(F, m, n, C, ldc, G);
        fassign(F,m,n,C,ldc,D,n);
        fger (F,m,n,alpha, A, inca, B, incb, C,ldc);
        ok = ok && check_fger(F, D, m,n,alpha, A, inca, B, incb, C,ldc);

        fflas_delete(A);
        fflas_delete(B);

        if (!ok)
            break;


    }
    fflas_delete (C);
    fflas_delete (D);

    return ok ;
}


template<class Field, class RandIter>
bool launch_fger_dispatch(const Field &F,
                          const size_t nn,
                          const typename Field::Element alpha,
                          const size_t iters,
                          RandIter& G)
{
    bool ok = true;
    size_t m,n;
    size_t inca,incb,ldc;
    //!@bug test for incx equal
    //!@bug test for transpo
    //!@todo does nbw actually do nbw recursive calls and then call blas (check ?) ?
    // size_t ld = 13 ;
    {
        m = 1+(size_t)random()%nn;
        n = 1+(size_t)random()%nn;


        // 		lda = m+(size_t)random()%ld;
        // 		ldb = 1+(size_t)random()%ld;

        inca = 1;
        incb = 1;

        // 		ldc = n+(size_t)random()%ld;
        ldc = n;

#ifdef __FFLASFFPACK_DEBUG
        std::cout <<"q = "<<F.characteristic()<<" m,n = "<<m<<", "<<n<<" C := "
        <<alpha<<".x * y^T + C";
#endif
        ok = ok && launch_fger<Field>(F,m,n, alpha, ldc, inca, incb, iters, G);
#ifdef __FFLASFFPACK_DEBUG
        std::cout<<(ok?" -> ok ":" -> KO")<<std::endl;
#endif
    }
    return ok ;
}
template <class Field>
bool run_with_field (int64_t q, uint64_t b, size_t n, size_t iters, uint64_t seed){
    bool ok = true ;
    int nbit=(int)iters;
    while (ok &&  nbit){
        typedef typename  Field::Element Element ;
        typedef typename Field::Element  Element ;

        Field* F= chooseField<Field>(q,b,seed);
        if (F==NULL) return true;
        std::ostringstream oss;
        F->write(oss);

        std::cout.fill('.');
        std::cout<<"Checking ";
        std::cout.width(45);
        std::cout<<oss.str();
        std::cout<<"... ";

        typename Field::RandIter R(*F,seed++);
        typename Field::NonZeroRandIter NZR(R);

        //size_t k = 0 ;
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_fger_dispatch<Field>(*F,n,F->one,iters, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_fger_dispatch<Field>(*F,n,F->zero,iters, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_fger_dispatch<Field>(*F,n,F->mOne,iters, R);
        //std::cout << k << "/24" << std::endl; ++k;

        Element alpha ;
        R.random(alpha);

        ok = ok && launch_fger_dispatch<Field>(*F,n,alpha,iters, R);
        if (!ok)
            //std::cout << "\033[1;31mFAILED\033[0m "<<std::endl;
            std::cout << "FAILED "<<std::endl;
        else
            //std::cout << "\033[1;32mPASSED\033[0m "<<std::endl;
            std::cout << "PASSED "<<std::endl;
        //std::cout<<std::endl;
        nbit--;
        delete F;
    }
    return ok;
}

int main(int argc, char** argv)
{
    std::cout<<setprecision(17);
    std::cerr<<setprecision(17);

    size_t iters = 3 ;
    long long q = -1 ;
    uint64_t b = 0 ;
    size_t n = 50 ;
    bool loop = false;
    uint64_t seed = getSeed();
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_LONGLONG , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
        { 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-loop Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
        { 's', "-s N", "Set the seed.",                         TYPE_UINT64 , &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    bool ok = true;
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<float> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,n,iters,seed);
        ok = ok && run_with_field<Modular<int64_t> >(q,b,n,iters,seed);
        ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,n,iters,seed);
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;

    return !ok ;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
