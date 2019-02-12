/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Ziad Sultan
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
//                        DSL test for pfgemm
//
//--------------------------------------------------------------------------
// Ziad Sultan
//-------------------------------------------------------------------------
/*
*/
#define NEWWINO
#ifndef TIME
#define TIME 1
#endif

#define __FFLASFFPACK_DEBUG 1
#include <iomanip>
#include <iostream>
using namespace std;

#define  __FFLASFFPACK_USE_OPENMP
//#define  __FFLASFFPACK_USE_KAAPI

//#define __FFLASFFPACK_FORCE_SEQ


#include "givaro/modular.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "time.h"

/*
#ifdef __FFLASFFPACK_USE_KAAPI
#include <kaapi++>
#endif

#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif
*/


using namespace FFPACK;

typedef Givaro::Modular<double> Field;
//typedef Givaro::Modular<float> Field;
//typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;
//typedef Givaro::Modular<int> Field;

BEGIN_PARALLEL_MAIN(int argc, char** argv)
{

    if (argc != 8)  {
        cerr<<"Testing pfgemm with : test-fgemm-DSL <p> <file-matrixA> <File-matrixB> <w> <i> <alpha> <beta>"
        <<endl;
        exit(-1);
    }
    srand48( FFLAS::BaseTimer::seed());
    size_t m,n, k;

    Field F(atoi(argv[1]));

    typename Field::Element *A;
    FFLAS::ReadMatrix (argv[2],F,m,k,A);
    typename Field::Element *B;
    FFLAS::ReadMatrix (argv[3],F,k,n,B);


    int nbw=atoi(argv[4]); // number of winograd levels
    int nbit=atoi(argv[5]); // number of times the product is performed
    cerr<<setprecision(10);
    Field::Element alpha,beta;

    F.init( alpha);
    F.assign( alpha, Field::Element(atoi(argv[6])));
    F.init( beta);
    F.assign( beta, Field::Element(atoi(argv[7])));

    size_t lda=m;
    size_t ldb=n;


    enum FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans;
    enum FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans;

    Field::Element * C=NULL;
    struct timespec t0,t1;
    double delay, avrg;
    double t_total=0;

    //const FFLAS::CuttingStrategy meth; //= FFLAS::BLOCK
    //const FFLAS::StrategyParameter strat;// = FFLAS::THREADS;
    FFLAS::ParSeqHelper::Parallel PSH(MAX_THREADS, FFLAS::CuttingStrategy::Block,
                                      FFLAS::StrategyParameter::Threads);
    FFLAS::MMHelper<Field,
    FFLAS::MMHelperAlgo::Winograd,
    FFLAS::FieldTraits<Field>::category,
    FFLAS::ParSeqHelper::Parallel> pWH (F, nbw, PSH);
    for(int i = 0;i<nbit;++i){
        C = FFLAS::fflas_new<Field::Element>(m*n);
        clock_gettime(CLOCK_REALTIME, &t0);

        //PAR_INSTR{

        FFLAS::fgemm(F, ta, tb,m,n,k,alpha, A,lda, B,ldb,
                     beta,C,n, pWH);
        //}
        BARRIER;
        clock_gettime(CLOCK_REALTIME, &t1);
        delay = (double)(t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)/1000000000;

        if (i)
            t_total+=delay;

    }
    avrg = t_total/(nbit-1);

#if TIME

    double mflops = (2.0*(m*k-((!F.isZero(beta))?m:0))/1000000.0)*n/avrg;

    cerr<<m<<" "<<n<<" "<<k<<" "<<nbw/*<<" "<<RBLOCKSIZE<<" "<<CBLOCKSIZE*/<<" "<<alpha<<" "<<beta<<" "
    <<mflops<<" "<<avrg<<endl;
#endif


#if __FFLASFFPACK_DEBUG
    bool wrong = false;
    Field::Element * Cd;
    Cd  = FFLAS::fflas_new<Field::Element>(m*n);
    for (int i=0; i<m*n; ++i)
        F.assign (*(Cd+i), F.zero);

    Field::Element aij, bij,  tmp;
    // Field::Element beta_alpha;
    //F.div (beta_alpha, beta, alpha);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j){
            F.mulin(*(Cd+i*n+j),beta);
            F.assign (tmp, F.zero);
            for ( int l = 0; l < k ; ++l ){
                if ( ta == FFLAS::FflasNoTrans )
                    aij = *(A+i*lda+l);
                else
                    aij = *(A+l*lda+i);
                if ( tb == FFLAS::FflasNoTrans )
                    bij = *(B+l*ldb+j);
                else
                    bij = *(B+j*ldb+l);
                //F.mul (tmp, aij, bij);
                //F.axpyin( *(Cd+i*n+j), alpha, tmp );
                F.axpyin (tmp, aij, bij);
            }
            F.axpyin (*(Cd+i*n+j), alpha, tmp);
            //F.mulin( *(Cd+i*n+j),alpha );
            if ( !F.areEqual( *(Cd+i*n+j), *(C+i*n+j) ) ) {
                wrong = true;
            }
        }
    if ( wrong ){
        cerr<<"FAIL"<<endl;
        for (int i=0; i<m; ++i){
            for (int j =0; j<n; ++j)
                if (!F.areEqual( *(C+i*n+j), *(Cd+i*n+j) ) )
                    cerr<<"Erreur C["<<i<<","<<j<<"]="
                    <<(*(C+i*n+j))<<" C[d"<<i<<","<<j<<"]="
                    <<(*(Cd+i*n+j))<<endl;
        }
    }
    else{
        cerr<<"PASS"<<endl;
    }
    FFLAS::fflas_delete( Cd);
#endif

    FFLAS::fflas_delete( C);
    FFLAS::fflas_delete( A);
    FFLAS::fflas_delete( B);


}
END_PARALLEL_MAIN()
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
