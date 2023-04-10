/* Copyright (c) FFLAS-FFPACK
 * Written by Clément Pernet <clement.pernet@imag.fr>
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
 */

//#include "goto-def.h"

// Please do not commit with any of these defines on - AB 2015-01-12
//#define __FFLASFFPACK_USE_TBB
//#define __FFLASFFPACK_USE_OPENMP
//#define __FFLASFFPACK_USE_DATAFLOW
//#define WINO_PARALLEL_TMPS
//#define __FFLASFFPACK_FORCE_SEQ
//#define PFGEMM_WINO_SEQ 32
//#define CLASSIC_SEQ
#define CLASSIC_HYBRID
//#define WINO_SEQ
//#define FFT_PROFILER
//#define PROFILE_FGEMM_MP
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular-balanced.h>

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"

#ifdef __FFLASFFPACK_USE_KAAPI
#include "libkomp.h"
#endif


using namespace std;
using namespace FFLAS;


int main(int argc, char** argv) {

    size_t iter = 3 ;
    Givaro::Integer q = 131071 ;
    size_t m = 2000 ;
    size_t k = 2000 ;
    size_t n = 2000 ;
    int nbw = -1 ;
    int p=0;
    int t=MAX_THREADS;
    int NBK = -1;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",         TYPE_INTEGER , &q },
        { 'm', "-m M", "Set the row dimension of A.",      TYPE_INT , &m },
        { 'k', "-k K", "Set the col dimension of A.",      TYPE_INT , &k },
        { 'n', "-n N", "Set the col dimension of B.",      TYPE_INT , &n },
        { 'w', "-w N", "Set the number of winograd levels (-1 for random).",    TYPE_INT , &nbw },
        { 'i', "-i R", "Set number of repetitions.",       TYPE_INT , &iter },
        { 'p', "-p P", "0 for sequential, 1 for 2D iterative, 2 for 2D rec, 3 for 2D rec adaptive, 4 for 3D rc in-place, 5 for 3D rec, 6 for 3D rec adaptive.", TYPE_INT , &p },
        { 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
        { 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
        END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);

    if (NBK==-1) NBK = t;
    //  typedef Givaro::Modular<Givaro::Integer> Field;
    // typedef Givaro::Modular<int64_t> Field;
    typedef Givaro::Modular<double> Field;
    //  typedef Givaro::Modular<float> Field;
    //  typedef Givaro::ModularBalanced<float> Field;
    //  typedef Givaro::ModularBalanced<double> Field;
    // typedef Givaro::ModularBalanced<int64_t> Field;
    //  typedef Givaro::Modular<Givaro::Integer> Field;
    typedef Field::Element Element;

    Field F(q);
    if (q > F.maxCardinality()) return 1;

    Timer chrono, TimFreivalds;
    double timev=0.0;
    double *time=new double[iter];

    Element * A, * B, * C;

    Field::RandIter G(F);
    A = fflas_new(F,m,k,Alignment::CACHE_PAGESIZE);
    //#pragma omp parallel for collapse(2) schedule(runtime)
    PAR_BLOCK { pfrand(F,G, m,k,A,m/size_t(NBK)); }

    B = fflas_new(F,k,n,Alignment::CACHE_PAGESIZE);
    //#pragma omp parallel for collapse(2) schedule(runtime)
    PAR_BLOCK { pfrand(F,G, k,n,B,k/NBK); }

    C = fflas_new(F,m,n,Alignment::CACHE_PAGESIZE);
    //#pragma omp parallel for collapse(2) schedule(runtime)
    PAR_BLOCK { pfzero(F, m,n,C,m/NBK); }


    for (size_t i=0;i<=iter;++i){

        chrono.clear();
        if (p && p!=7){
            // CuttingStrategy meth = RECURSIVE;
            // StrategyParameter strat = THREADS;

            typedef CuttingStrategy::Block block;
            typedef CuttingStrategy::Recursive rec;
            typedef StrategyParameter::Threads threads;
            typedef StrategyParameter::TwoD  twod;
            typedef StrategyParameter::TwoDAdaptive  twoda;
            typedef StrategyParameter::ThreeD  threed;
            typedef StrategyParameter::ThreeDAdaptive  threeda;
            typedef StrategyParameter::ThreeDInPlace  threedip;
            PAR_BLOCK{
                if (i) { chrono.start(); }

                switch (p){
                case 1:{
                           MMHelper<Field, MMHelperAlgo::Winograd, typename ModeTraits<Field>::value, ParSeqHelper::Parallel<block,threads> > WH(F,nbw, SPLITTER(t,block,threads));
                           fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n, WH);
                           break;}
                case 2:{
                           MMHelper<Field, MMHelperAlgo::Winograd, typename ModeTraits<Field>::value, ParSeqHelper::Parallel<rec,twod> > WH(F,nbw, SPLITTER(t,rec,twod));
                           fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n, WH);
                           break;
                       }
                case 3:{
                           MMHelper<Field, MMHelperAlgo::Winograd, typename ModeTraits<Field>::value, ParSeqHelper::Parallel<rec,twoda> > WH(F,nbw, SPLITTER(t,rec,twoda));
                           fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n, WH);
                           break;
                       }
                case 4:{
                           MMHelper<Field, MMHelperAlgo::Winograd, typename ModeTraits<Field>::value, ParSeqHelper::Parallel<rec,threedip> > WH(F,nbw, SPLITTER(t,rec,threedip));
                           fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n, WH);
                           break;
                       }
                case 5:{
                           MMHelper<Field, MMHelperAlgo::Winograd, typename ModeTraits<Field>::value, ParSeqHelper::Parallel<rec,threed> > WH(F,nbw, SPLITTER(t,rec,threed));
                           fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n, WH);
                           break;
                       }
                case 6:{
                           MMHelper<Field, MMHelperAlgo::Winograd, typename ModeTraits<Field>::value, ParSeqHelper::Parallel<rec,threeda> > WH(F,nbw, SPLITTER(t,rec,threeda));
                           fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n, WH);
                           break;
                       }
                default:{
                            MMHelper<Field, MMHelperAlgo::Winograd, typename ModeTraits<Field>::value, ParSeqHelper::Parallel<block,threads> > WH(F,nbw, SPLITTER(t,block,threads));
                            fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n, WH);
                            break;
                        }
                }
            }
            if (i) {chrono.stop(); time[i-1]=chrono.realtime();}
        }else{
            if(p==7){

                int nrec = 0;
                int dim = m;
                //                      if(dim < 19000)
                nrec--;
                while(dim >= __FFLASFFPACK_WINOTHRESHOLD*2){
                    dim=dim/2;
                    nrec++;
                }
                nrec=std::max(1,nrec);
                //              std::cout<<" WINO_THREShold"<<__FFLASFFPACK_WINOTHRESHOLD<<" nrec = "<<nrec<<" dim = "<<dim<<std::endl;
                if(nbw != -1)
                    nrec=nbw;
                nbw=nrec;
                if (i) chrono.start();
                PAR_BLOCK
                {
                    MMHelper<Field, MMHelperAlgo::WinogradPar,ModeTraits<Field>::value,ParSeqHelper::Parallel<> >  WH (F, nrec, ParSeqHelper::Parallel<>(t));
                    fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n,WH);
                }
                if (i) {chrono.stop(); time[i-1]=chrono.realtime();}


                // MMHelper<Field, MMHelperAlgo::WinogradPar>
                //          WH (F, nbw, ParSeqHelper::Sequential());
                //      //              cout<<"wino parallel"<<endl;
                // if (i) chrono.start();
                // PAR_BLOCK
                // {
                //          fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n,WH);
                // }
                // if (i) {chrono.stop(); time+=chrono.realtime();}
            }
            else{

                MMHelper<Field,MMHelperAlgo::Winograd>//,
                //typename FieldTraits<Field>::value,
                //ParSeqHelper::Sequential>
                WH (F, nbw, ParSeqHelper::Sequential());
                if (i) chrono.start();
                fgemm (F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, F.zero, C,n,WH);
                if (i) {chrono.stop(); time[i-1]=chrono.realtime();}
            }
        }

        TimFreivalds.clear();
        TimFreivalds.start();

        bool pass = freivalds(F, FflasNoTrans, FflasNoTrans, m,n,k, F.one, A, k, B, n, C,n);
        TimFreivalds.stop();
        timev+=TimFreivalds.usertime();
        if (!pass)
            std::cout<<"FAILED"<<std::endl;
    }
    fflas_delete( A);
    fflas_delete( B);
    fflas_delete( C);
    std::sort(time, time+iter);
    double mediantime = time[iter/2];
    delete[] time;

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/11/14
    std::cout << "Time: " << mediantime
    << " Gfops: " << (2.*double(m)/1000.*double(n)/1000.*double(k)/1000.0) / mediantime;
    writeCommandString(std::cout, as) << std::endl;

#if __FFLASFFPACK_DEBUG
    std::cout<<"Freivalds vtime: "<<timev/(double)iter<<std::endl;
#endif

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
