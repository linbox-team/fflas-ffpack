/* Copyright (c) FFLAS-FFPACK
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


// declare that the call to openblas_set_numthread will be made here, hence don't do it
// everywhere in the call stack
#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas/fflas.h"

#include <iostream>

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/args-parser.h"

#ifdef __FFLASFFPACK_USE_KAAPI
#include "libkomp.h"
#endif

using namespace FFLAS;

typedef FFPACK::rns_double RNS;
typedef FFPACK::RNSInteger<RNS> Field;
typedef Field::Element_ptr Element_ptr;
typedef Field::ConstElement_ptr ConstElement_ptr;

typedef StrategyParameter::Threads THREADS;
typedef StrategyParameter::Grain GRAIN;
typedef StrategyParameter::TwoD TWOD;
typedef StrategyParameter::TwoDAdaptive TWODA;
typedef StrategyParameter::ThreeD THREED;
typedef StrategyParameter::ThreeDAdaptive THREEDA;
typedef StrategyParameter::ThreeDInPlace THREEDIP;

typedef ParSeqHelper::Sequential PSeq;

template<typename ModParSeqTrait, typename FgemmParSeqTrait>
static inline void
bench_fgemm (const Field ZZ, const size_t m, const size_t n, const size_t k,
             ConstElement_ptr A, ConstElement_ptr B, Element_ptr C,
             int moduli_th, int fgemm_th, int nbw)
{
    PAR_BLOCK
    {
        typedef ParSeqHelper::Compose<ModParSeqTrait, FgemmParSeqTrait> RNSParSeqHelper;
        RNSParSeqHelper RNSParSeq (moduli_th, fgemm_th);

        typedef MMHelper<Field, MMHelperAlgo::Classic, ModeCategories::DefaultTag, RNSParSeqHelper> RNSHelper;
        RNSHelper WH (ZZ, nbw, RNSParSeq);

        fgemm (ZZ, FflasNoTrans, FflasNoTrans, m, n, k, ZZ.one, A, k, B, n, ZZ.zero,
               C, n, WH);
    }
}

template<typename T>
static inline void
bench_do_it (const Field ZZ, const size_t m, const size_t n, const size_t k,
             ConstElement_ptr A, ConstElement_ptr B, Element_ptr C,
             int moduli_th, int p, int fgemm_th, int nbw)
{
    if (p == 0)
        bench_fgemm<T, PSeq> (ZZ, m, n, k, A, B, C, moduli_th, fgemm_th, nbw);
    else if (p == 1)
    {
        typedef ParSeqHelper::Parallel<CuttingStrategy::Block, THREADS> PPar;
        bench_fgemm<T, PPar> (ZZ, m, n, k, A, B, C, moduli_th, fgemm_th, nbw);
    }
    else if (p == 2)
    {
        typedef ParSeqHelper::Parallel<CuttingStrategy::Recursive, TWOD> PPar;
        bench_fgemm<T, PPar> (ZZ, m, n, k, A, B, C, moduli_th, fgemm_th, nbw);
    }
    else if (p == 3)
    {
        typedef ParSeqHelper::Parallel<CuttingStrategy::Recursive, TWODA> PPar;
        bench_fgemm<T, PPar> (ZZ, m, n, k, A, B, C, moduli_th, fgemm_th, nbw);
    }
    else if (p == 4)
    {
        typedef ParSeqHelper::Parallel<CuttingStrategy::Recursive, THREEDIP> PPar;
        bench_fgemm<T, PPar> (ZZ, m, n, k, A, B, C, moduli_th, fgemm_th, nbw);
    }
    else if (p == 5)
    {
        typedef ParSeqHelper::Parallel<CuttingStrategy::Recursive, THREED> PPar;
        bench_fgemm<T, PPar> (ZZ, m, n, k, A, B, C, moduli_th, fgemm_th, nbw);
    }
    else /* p == 6 */
    {
        typedef ParSeqHelper::Parallel<CuttingStrategy::Recursive, THREEDA> PPar;
        bench_fgemm<T, PPar> (ZZ, m, n, k, A, B, C, moduli_th, fgemm_th, nbw);
    }
}

int
main(int argc, char *argv[])
{
 
#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t pbits = 20;
    size_t r = 8;
    size_t m = 1000;
    size_t k = 1000;
    size_t n = 1000;
    int nbw = -1;
    size_t iter = 3;
    int fgemm_th = MAX_THREADS;
    int p = 0;
    int moduli_th = 1;
    int q = 0;

    Argument as[] = {
        { 'r', "-r R", "Number of RNS moduli", TYPE_INT , &r},
        { 'b', "-b B", "Number of bits of the moduli (in [10, 26])", TYPE_INT,
            &pbits},
        { 'm', "-m M", "Row dimension of A", TYPE_INT, &m},
        { 'k', "-k K", "Col dimension of A", TYPE_INT, &k},
        { 'n', "-n N", "Col dimension of B", TYPE_INT, &n},
        { 'w', "-w N", "Number of Winograd levels (-1 means computed by the lib)",
            TYPE_INT, &nbw},
        { 'i', "-i R", "Number of repetitions", TYPE_INT, &iter},
        { 'p', "-p P", "0 for sequential, 1 for 2D iterative, 2 for 2D rec, "
            "3 for 2D rec adaptive, 4 for 3D rc in-place, "
            "5 for 3D rec, 6 for 3D rec adaptive.", TYPE_INT , &p },
        { 't', "-t T", "Number of threads for the fgemm computations.", TYPE_INT,
            &fgemm_th},
        { 'q', "-q Q", "Strategy parameter for the moduli: 0 for sequential, "
            "1 for threads", TYPE_INT , &q},
        { 'u', "-u U", "Number of threads for handling the moduli "
            "(depends on the value of q) .", TYPE_INT, &moduli_th},
        END_OF_ARGUMENTS
    };

    parseArguments (argc, argv, as);
    std::ostringstream oss;
    /* Check some command line parameters */
    if (fgemm_th <= 0 || moduli_th <= 0)
    {
        std::cerr << "Error, the number of threads must be positive" << std::endl;
        return 1;
    }
    if (pbits < 10 || pbits > 26)
    {
        std::cerr << "Error, the number of bits of the moduli must be in [10, 26]"
        << std::endl;
        return 1;
    }

    /* Some warnings */
    if (q == 0 && moduli_th > 1)
    {
        oss << "Warning, -u is set to 1 because -q is 0 (sequential)"
        << std::endl;
        moduli_th = 1;
    }
    if (p == 0 && fgemm_th > 1)
    {
        oss << "Warning, -t is set to 1 because -p is 0 (sequential)"
        << std::endl;
        fgemm_th = 1;
    }

    RNS rns (pbits, r);
    Field ZZ(rns);

    Timer chrono, TimFreivalds;
    double time=0.0, timev=0.0;

    Element_ptr A, B, C;

    Field::RandIter G(ZZ);

    A = fflas_new (ZZ, m, k, Alignment::CACHE_PAGESIZE);
    frand (ZZ, G, m*k, A, 0);
    B = fflas_new (ZZ, k, n, Alignment::CACHE_PAGESIZE);
    frand (ZZ, G, k*n, B, 0);
    C = fflas_new (ZZ, m, n, Alignment::CACHE_PAGESIZE);
    fzero (ZZ, m*n, C, 0);

    for (size_t i=0; i<=iter; ++i)
    {
        chrono.clear();
        if (i)
            chrono.start();

        if (q == 0) /* moduli are done sequentially */
            bench_do_it<PSeq> (ZZ, m, n, k, A, B, C, moduli_th, p, fgemm_th, nbw);
        else /* (q == 1) */
        {
            typedef ParSeqHelper::Parallel<CuttingStrategy::RNSModulus, THREADS> PPar;
            bench_do_it<PPar> (ZZ, m, n, k, A, B, C, moduli_th, p, fgemm_th, nbw);
        }

        if (i)
        {
            chrono.stop();
            time+=chrono.realtime();
        }

        TimFreivalds.clear();
        TimFreivalds.start();

        bool pass = freivalds (ZZ, FflasNoTrans, FflasNoTrans, m, n, k, ZZ.one, A, k, B, n, C,n);
        TimFreivalds.stop();
        timev += TimFreivalds.usertime();
        if (!pass)
            std::cout << "FAILED" << std::endl;
    }

    fflas_delete (A);
    fflas_delete (B);
    fflas_delete (C);

    // -----------
    // Standard output for benchmark
    double time_per_iter = time / double(iter);
    double Gflops = double(iter) * double (r) * (2. * double(m)/1000. * double(n)/1000. * double(k)/1000.0) / time;
    std::cout << "Time: " << time_per_iter << " Gflops: " << Gflops;
    writeCommandString (std::cout, as) << std::endl;
#if __FFLASFFPACK_DEBUG
    std::cout << "Freivalds vtime: " << timev / (double)iter << std::endl;
#endif
    std::cerr<<oss.str();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
