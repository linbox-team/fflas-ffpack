/*
 * Copyright (C) the FFLAS-FFPACK group
 * Written by Clément Pernet
 *            Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

// #ifndef NEWINO
// #define NEWWINO
// #endif

// #define WINOTHRESHOLD 100
// #define OLD_DYNAMIC_PEELING



#include "fflas-ffpack/fflas-ffpack-config.h"

#include <iomanip>
#include <iostream>

#include <givaro/modular.h>

#include <recint/rint.h>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"

#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"

using namespace std;
using namespace FFPACK;
using namespace FFLAS;

using Givaro::Modular;
using Givaro::ModularBalanced;


// checks that D = beta . Y + alpha . A ^ta * X
template<class Field>
bool check_MV(const Field                   & F,
              const typename Field::Element_ptr  Cd, // c0
              enum FFLAS_TRANSPOSE   & ta,
              const size_t                    m,
              const size_t                    k,
              const typename Field::Element & alpha,
              const typename Field::Element_ptr  A, size_t lda,
              const typename Field::Element_ptr  X, size_t incX,
              const typename Field::Element & beta,
              const typename Field::Element_ptr Y, size_t incY)
{
    bool wrong = false;
    typename Field::Element_ptr D;
    if (ta == FflasNoTrans){
        D = fflas_new(F,m);
        fassign (F, m, Cd, 1, D, 1);
        for (size_t i=0; i<m; i++){
            F.mulin (D[i], beta);
            typename Field::Element tmp;
            F.init(tmp);
            for (size_t j=0; j<k; j++){
                F.axpyin (tmp, A[i*lda+j], X[j*incX]);
            }
            F.axpyin(D[i],alpha,tmp);
        }
        wrong = !fequal(F, m, D, 1, Y, incY);
    } else {
        D = fflas_new(F,k);
        fassign (F, k, Cd, 1, D, 1);
        for (size_t i=0; i<k; i++){
            F.mulin (D[i], beta);
            typename Field::Element tmp;
            F.init(tmp);
            for (size_t j=0; j<m; j++){
                F.axpyin (tmp, A[i+j*lda], X[j*incX]);
            }
            F.axpyin(D[i],alpha,tmp);
        }
        wrong = !fequal(F, k, D, 1, Y, incY);
    }
    size_t Ydim = (ta==FflasNoTrans)? m : k;
    if ( wrong ){
        size_t canprint = 20 ;
        std::cerr<<"FAIL"<<std::endl;
        std::cerr << "alpha   :" << alpha<<", beta   : " << beta << std::endl;
        std::cerr << "m   :" << m   << ", k   : " << k << std::endl;
        std::cerr << "ldA :" << lda << ", incX : " << incX << ", incY : " << incY << std::endl;
        for (size_t i=0; i<Ydim && canprint; ++i){
            if (!F.areEqual( Y[i*incY], D[i] ) ) {
                std::cerr<<"Error Y["<<i<<"]="<<Y[i*incY]<<" D["<<i<<"]="<<D[i]<<std::endl;
                canprint--;
            }
        }
        if (Ydim<80) {
            for (size_t i=0; i<Ydim ; ++i){
                if (!F.areEqual( Y[i*incY], D[i] ) )
                    std::cout << 'X' ;
                else
                    std::cout << '.' ;
            }
            std::cout << std::endl;
        }
    }

    fflas_delete (D);

    return !wrong ;

}


template<class Field, class RandIter>
bool launch_MV(const Field & F,
               const size_t   m,
               const size_t   k,
               const typename Field::Element alpha,
               const typename Field::Element beta,
               const size_t lda,
               enum FFLAS_TRANSPOSE    ta,
               const size_t incX,
               const size_t incY,
               size_t iters,
               bool par,
               RandIter& G)
{

    bool ok = true;

    typedef typename Field::Element_ptr Element_ptr;
    Element_ptr A ;
    for(size_t i = 0;i<iters;++i){
        FFLASFFPACK_check(lda >= k);
        A = fflas_new (F, m, lda);
        fzero(F,m,lda,A,lda);
        RandomMatrix(F, m, k, A, lda, G);
        size_t Xdim = (ta == FflasNoTrans)? k : m;
        size_t Ydim = (ta == FflasNoTrans)? m : k;
        Element_ptr X = fflas_new (F, Xdim, incX);
        Element_ptr Y = fflas_new (F, Ydim, incY);
        fzero (F, Xdim, incX, X, incX);
        fzero (F, Ydim, incY, Y, incY);
        Element_ptr D = fflas_new (F, Ydim);

        RandomMatrix (F, Xdim, 1, X, incX, G);
        RandomMatrix (F, Ydim, 1, Y, incY, G);
        fassign (F, Ydim, Y, incY, D, 1);

        //Y will be modified so keep a copy of Y as Y2
        Element_ptr Y2 =  fflas_new (F, Ydim, incY);
        fassign (F, Ydim, Y, incY, Y2, incY);


        if (par){
            {
                MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Recursive, StrategyParameter::Threads> >  WH;

                PAR_BLOCK{
                    fgemv(F, ta, m,k,alpha, A,lda, X, incX, beta, Y, incY, WH);
                }
            }
        }else{
            //MMHelper<Field,MMHelperAlgo::Auto,typename ModeTraits<Field>::value> WH(F,nbw,ParSeqHelper::Sequential());
            fgemv(F, ta, m, k,alpha, A,lda, X, incX, beta, Y, incY);
        }

        ok = ok && check_MV(F, D, ta, m, k,alpha, A, lda, X, incX, beta, Y, incY);



        if (!ok){
            fflas_delete (A, X, Y, Y2, D);
            break;
        }



        fassign (F, Ydim, Y2, incY, D, 1);
        fassign (F, Ydim, Y2, incY, Y, incY);
        if (par){
            {
                MMHelper<Field, MMHelperAlgo::Classic, ModeTraits<Field>, ParSeqHelper::Parallel<CuttingStrategy::Row, StrategyParameter::Threads> >  WH;

                PAR_BLOCK{
                    fgemv(F, ta, m,k,alpha, A,lda, X, incX, beta, Y, incY, WH);
                }
            }
        }else{
            //MMHelper<Field,MMHelperAlgo::Auto,typename ModeTraits<Field>::value> WH(F,nbw,ParSeqHelper::Sequential());
            fgemv(F, ta, m, k,alpha, A,lda, X, incX, beta, Y, incY);
        }

        ok = ok && check_MV(F, D, ta, m, k,alpha, A, lda, X, incX, beta, Y, incY);


        if (!ok){
            fflas_delete (A, X, Y, Y2, D);
            break;
        }




        fassign (F, Ydim, Y2, incY, D, 1);
        fassign (F, Ydim, Y2, incY, Y, incY);
        if (par){
            {
                ParSeqHelper::Parallel<CuttingStrategy::Row,StrategyParameter::Grain>  WH(4);

                PAR_BLOCK{
                    fgemv(F, ta, m,k,alpha, A,lda, X, incX, beta, Y, incY, WH);
                }
            }
        }else{
            //MMHelper<Field,MMHelperAlgo::Auto,typename ModeTraits<Field>::value> WH(F,nbw,ParSeqHelper::Sequential());
            fgemv(F, ta, m, k,alpha, A,lda, X, incX, beta, Y, incY);
        }

        ok = ok && check_MV(F, D, ta, m, k,alpha, A, lda, X, incX, beta, Y, incY);

        fflas_delete (A, X, Y, Y2, D);
        if (!ok){

            break;
        }



    }
    return ok ;
}


template<class Field, class RandIter>
bool launch_MV_dispatch(const Field &F,
                        const int mm,
                        const int kk,
                        const typename Field::Element alpha,
                        const typename Field::Element beta,
                        const size_t iters,
                        const bool par,
                        RandIter& G)
{
    bool ok = true;
    size_t m,k;
    size_t lda,incX, incY;
    size_t ld = 13 ;
    {
        //FFLAS_TRANSPOSE ta = FflasNoTrans ;
        //if (! par) {
        //if (random()%2) ta = FflasTrans ;
        //}

        if (mm<0)
            m = 1+(size_t)random() % -mm;
        else m = mm;
        if (kk<0)
            k = 1+(size_t)random() % -kk;
        else k = kk;

        lda = k+(size_t)random()%ld;
        incX = 1+(size_t)random()%ld;
        incY = 1+(size_t)random()%ld;

        ok = ok && launch_MV (F, m, k, alpha,beta, lda, FflasNoTrans, incX, incY, iters, par, G);
        ok = ok && launch_MV (F, m, k, alpha,beta, lda, FflasTrans, incX, incY, iters, par, G);
    }
    return ok ;
}

template <class Field>
bool run_with_field (Givaro::Integer q, uint64_t b, int m, int k, size_t iters, bool par, uint64_t seed){
    bool ok = true ;

    int nbit=(int)iters;

    while (ok &&  nbit){
        typedef typename Field::Element Element ;
        // choose Field
        Field* F= chooseField<Field>(q,b,seed);
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
        typename Field::RandIter R(*F,seed++);
        typename Field::NonZeroRandIter NZR(R);

        //size_t k = 0 ;
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->one,F->zero,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->zero,F->zero,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->mOne,F->zero,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->one ,F->one,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->zero,F->one,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->mOne,F->one,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->one ,F->mOne,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->zero,F->mOne,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->mOne,F->mOne,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;

        Element alpha,beta ;
        NZR.random(alpha);
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->one ,alpha,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->zero,alpha,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,F->mOne,alpha,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,alpha,F->one ,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,alpha,F->zero,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;
        ok = ok && launch_MV_dispatch<Field>(*F,m,k,alpha,F->mOne,iters, par, R);
        //std::cout << k << "/24" << std::endl; ++k;

        for (size_t j = 0 ; j < 3 ; ++j) {
            R.random(alpha);
            R.random(beta);
            ok = ok && launch_MV_dispatch<Field>(*F,m,k,alpha,beta,iters, par, R);
            //std::cout << k << "/24" << std::endl; ++k;
        }
        //std::cout<<std::endl;
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
    std::cout<<setprecision(17);
    std::cerr<<setprecision(17);

    uint64_t seed = getSeed();
    size_t iters = 3 ;
    Givaro::Integer q = -1 ;
    uint64_t b = 0 ;
    int m = -50 ;
    int k = -50 ;
    int nbw = -1 ;
    bool loop = false;
    bool p = false;
    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'b', "-b B", "Set the bitsize of the random characteristic.",         TYPE_INT , &b },
        { 'm', "-m M", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |n|).",      TYPE_INT , &m },
        { 'k', "-k K", "Set the dimension of the matrix (negative values, mean, any random value between 0 and |k|).",      TYPE_INT , &k },
        { 'w', "-w N", "Set the number of winograd levels (-1 for random).",    TYPE_INT , &nbw },
        { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
        { 'l', "-l Y/N", "run the test in an infinte loop.", TYPE_BOOL , &loop },
        { 'p', "-p Y/N", "run the parallel fgemv.", TYPE_BOOL , &p },
        { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    bool ok = true;
    srand(seed);
    do{
        ok = ok && run_with_field<Modular<double> >(q,b,m,k,iters,p, seed);
        ok = ok && run_with_field<ModularBalanced<double> >(q,b,m,k,iters,p, seed);
        ok = ok && run_with_field<Modular<float> >(q,b,m,k,iters,p, seed);
        ok = ok && run_with_field<ModularBalanced<float> >(q,b,m,k,iters,p, seed);
        ok = ok && run_with_field<Modular<int32_t> >(q,b,m,k,iters,p, seed);
        ok = ok && run_with_field<ModularBalanced<int32_t> >(q,b,m,k,iters,p, seed);
        ok = ok && run_with_field<Modular<int64_t> >(q,b,m,k,iters, p, seed);
        ok = ok && run_with_field<ModularBalanced<int64_t> >(q,b,m,k,iters, p, seed);

        ok = ok && run_with_field<Modular<RecInt::rint<8> > >(q,b?b:127_ui64,m,k,iters, p, seed);
        ok = ok && run_with_field<Modular<RecInt::ruint<7>,RecInt::ruint<8> > >(q,b?b:127_ui64,m,k,iters, p, seed);
        ok = ok && run_with_field<Modular<Givaro::Integer> >(q,(b?b:512_ui64),m,k,iters,p, seed);
        ok = ok && run_with_field<Givaro::ZRing<Givaro::Integer> >(0,(b?b:512_ui64),m,k,iters,p, seed);
    } while (loop && ok);




    return !ok ;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
