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

#define ENABLE_ALL_CHECKINGS 1 // DO NOT CHANGE
#define _NR_TESTS 5
#define _MAX_SIZE_MATRICES 1000

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/checkers/checkers_fflas.h"
#include "fflas-ffpack/checkers/checkers_ffpack.h"
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
    size_t NR_TESTS = _NR_TESTS;
    int    q    = 131071;
    size_t    MAX_SIZE_MATRICES    = _MAX_SIZE_MATRICES;
    size_t Range = 500;
    size_t seed( (int) time(NULL) );
    std::string file("checkers_report.txt");

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
        { 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &MAX_SIZE_MATRICES },
        { 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &NR_TESTS },
        { 'r', "-r R", "Set the range of matrix sizes.",                     TYPE_INT , &Range },
        { 's', "-s N", "Set the seed.", TYPE_INT , &seed },
        { 'f', "-f FILE", "Set the output file.",  TYPE_STR , &file },
        END_OF_ARGUMENTS
    };

    std::ofstream stats_f(file.c_str());

    FFLAS::parseArguments(argc,argv,as);

    srand (seed);

    typedef Givaro::Modular<double> Field;
    typedef Givaro::Poly1Dom<Field> PolRing;
    typedef PolRing::Element Polynomial;

    Field F(q);
    Field::RandIter Rand(F,seed);
    Field::NonZeroRandIter NZRand(Rand);

    size_t pass;
    FFLAS::Timer chrono,global;
    double gffop(0.);
    global.start();
    double time1, time2;

    Field::Element_ptr A = FFLAS::fflas_new(F,MAX_SIZE_MATRICES+Range,MAX_SIZE_MATRICES+Range);
    Field::Element_ptr B = FFLAS::fflas_new(F,MAX_SIZE_MATRICES+Range,MAX_SIZE_MATRICES+Range);
    Field::Element_ptr C = FFLAS::fflas_new(F,MAX_SIZE_MATRICES+Range,MAX_SIZE_MATRICES+Range);
    typename Field::Element alpha,beta,tmp;
    F.init(alpha, rand()%1000+1);
    F.init(beta,  rand()%1000+1);
    size_t m,n,k,lda,ldb,ldc;
    FFLAS::FFLAS_TRANSPOSE ta,tb;

    stats_f << "     Matrix size\tSuccess rate\t\tTime comput.\t\tTime checker\n\n";

    // #####   FGEMM   #####
    stats_f << "FGEMM:\n";
    for (size_t i=0; i<MAX_SIZE_MATRICES; i+=Range) {
        pass = 0; time1 = 0.0; time2 = 0.0;
        for (size_t j=0; j<NR_TESTS; ++j) {
            m = rand() % Range + i;
            n = rand() % Range + i;
            k = rand() % Range + i;
            gffop += (2.*double(m)/1000.*double(n)/1000.*double(k)/1000.0);

            ta = FFLAS::FflasNoTrans;//rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans,
            tb = FFLAS::FflasNoTrans;//rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans;
            lda = ta == FFLAS::FflasNoTrans ? k : m,
                ldb = tb == FFLAS::FflasNoTrans ? n : k,
                ldc = n;

            PAR_BLOCK { FFLAS::pfrand(F,Rand, m,k,A,m/MAX_THREADS); }
            PAR_BLOCK { FFLAS::pfrand(F,Rand, k,n,B,k/MAX_THREADS); }
            PAR_BLOCK { FFLAS::pfrand(F,Rand, m,n,C,n/MAX_THREADS); }

            chrono.clear(); chrono.start();
            FFLAS::ForceCheck_fgemm<Field> checker1(Rand,m,n,k,beta,C,ldc);
            chrono.stop(); time1 += chrono.usertime();

            chrono.clear(); chrono.start();
            FFLAS::fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
            chrono.stop(); time2 += chrono.usertime();

            chrono.clear(); chrono.start();
            pass += checker1.check(ta,tb,alpha,A,lda,B,ldb,C) ? 1 : 0;
            chrono.stop(); time1 += chrono.usertime();
        }
        time1 /= NR_TESTS;
        time2 /= NR_TESTS;
        stats_f << "     " << i << "-" << i+Range << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2
        << "\t\t" << time1 << endl;
    }
    stats_f << endl;



    // #####   FTRSM   #####
    stats_f << "FTRSM:\n";
    for (size_t i=0; i<MAX_SIZE_MATRICES; i+=Range) {
        pass = 0; time1 = 0.0; time2 = 0.0;
        for (size_t j=0; j<NR_TESTS; ++j) {
            m = rand() % Range + i;
            n = rand() % Range + i;
            gffop += (double(m)/1000.*double(m)/1000.*double(n)/1000.0);

            FFLAS::FFLAS_SIDE side = rand()%2?FFLAS::FflasLeft:FFLAS::FflasRight;
            FFLAS::FFLAS_UPLO uplo = rand()%2?FFLAS::FflasLower:FFLAS::FflasUpper;
            FFLAS::FFLAS_TRANSPOSE trans = rand()%2?FFLAS::FflasNoTrans:FFLAS::FflasTrans;
            FFLAS::FFLAS_DIAG diag = rand()%2?FFLAS::FflasNonUnit:FFLAS::FflasUnit;
            k = (side==FFLAS::FflasLeft?m:n);

            for( size_t i = 0; i < m*n; ++i ) Rand.random( *(B+i) );
            for (size_t i=0;i<k;++i) {
                for (size_t j=0;j<i;++j)
                    A[i*k+j]= (uplo == FFLAS::FflasLower)? Rand.random(tmp) : F.zero;
                A[i*k+i]= (diag == FFLAS::FflasNonUnit)? NZRand.random(tmp) : F.one;
                for (size_t j=i+1;j<k;++j)
                    A[i*k+j]= (uplo == FFLAS::FflasUpper)? Rand.random(tmp) : F.zero;
            }

            chrono.clear(); chrono.start();
            FFLAS::ForceCheck_ftrsm<Field> checker2(Rand, m, n, alpha, B, n);
            chrono.stop(); time1 += chrono.usertime();

            chrono.clear(); chrono.start();
            FFLAS::ftrsm(F, side, uplo, trans, diag, m, n, alpha, A, k, B, n);
            chrono.stop(); time2 += chrono.usertime();

            chrono.clear(); chrono.start();
            pass += checker2.check(side, uplo, trans, diag, m, n, A, k, B, n);
            chrono.stop(); time1 += chrono.usertime();
        }
        time1 /= NR_TESTS;
        time2 /= NR_TESTS;
        stats_f << "     " << i << "-" << i+Range << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2
        << "\t\t" << time1 << endl;
    }
    stats_f << endl;



    // #####   INVERT   #####
    stats_f << "INVERT:\n";
    int nullity;
    for (size_t i=0; i<MAX_SIZE_MATRICES; i+=Range) {
        pass = 0; time1 = 0.0; time2 = 0.0;
        for (size_t j=0; j<NR_TESTS; ++j) {
            m = rand() % Range + i;
            gffop += 2*(double(m)/1000.*double(m)/1000.*double(m)/1000.0);

            FFPACK::RandomMatrixWithRankandRandomRPM(F,m,m,m,A,m);

            try {
                chrono.clear(); chrono.start();
                FFPACK::ForceCheck_invert<Field> checker3(Rand,m,A,m);
                chrono.stop(); time1 += chrono.usertime();

                chrono.clear(); chrono.start();
                FFPACK::Invert(F,m,A,m,nullity);
                chrono.stop(); time2 += chrono.usertime();

                chrono.clear(); chrono.start();
                pass += checker3.check(A,nullity);
                chrono.stop(); time1 += chrono.usertime();
            } catch(FailureInvertCheck &e) {
                stats_f << " invert verification failed! " << nullity << std::endl;
            } catch(FailurePLUQCheck &e) {
                stats_f << " internal PLUQ verification failed! " << std::endl;
            }
        }
        time1 /= NR_TESTS;
        time2 /= NR_TESTS;
        stats_f << "     " << i << "-" << i+Range << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2
        << "\t\t" << time1 << endl;
    }
    stats_f << endl;




    // #####   PLUQ   #####
    stats_f << "PLUQ:\n";
    for (size_t i=0; i<MAX_SIZE_MATRICES; i+=Range) {
        pass = 0; time1 = 0.0; time2 = 0.0;
        for (size_t j=0; j<NR_TESTS; ++j) {
            m = rand() % Range + i;
            n = rand() % Range + i;

            PAR_BLOCK { FFLAS::pfrand(F,Rand, m,n,A,m/MAX_THREADS); }

            size_t *P = FFLAS::fflas_new<size_t>(m);
            size_t *Q = FFLAS::fflas_new<size_t>(n);

            chrono.clear(); chrono.start();
            FFPACK::ForceCheck_PLUQ<Field> checker4 (Rand,m,n,A,n);
            chrono.stop(); time1 += chrono.usertime();

            chrono.clear(); chrono.start();
            k = FFPACK::PLUQ(F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);
            chrono.stop(); time2 += chrono.usertime();

#define CUBE(x) ((x)*(x)*(x))
            gffop += 2.0/3.0*CUBE(double(k)/1000.0) +2*m/1000.0*n/1000.0*double(k)/1000.0  - double(k)/1000.0*double(k)/1000.0*(m+n)/1000;

            chrono.clear(); chrono.start();
            pass += checker4.check(A,n,FFLAS::FflasNonUnit, k,P,Q);
            chrono.stop(); time1 += chrono.usertime();

            FFLAS::fflas_delete(P,Q);
        }
        time1 /= NR_TESTS;
        time2 /= NR_TESTS;
        stats_f << "     " << i << "-" << i+Range << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2
        << "\t\t" << time1 << endl;
    }
    stats_f << endl;
    global.stop();



    // #####   CharPoly   #####
    stats_f << "CharPoly:\n";
    PolRing R(F);

    for (size_t i=0; i<MAX_SIZE_MATRICES; i+=Range) {
        pass = 0; time1 = 0.0; time2 = 0.0;
        for (size_t j=0; j<NR_TESTS; ++j) {
            n = rand() % Range + i;

            PAR_BLOCK { FFLAS::pfrand(F,Rand, n,n,A,n/MAX_THREADS); }

            try {
                Polynomial g(n);

                chrono.clear(); chrono.start();
                FFPACK::ForceCheck_charpoly<Field,Polynomial> checker5(Rand,n,A,n);
                chrono.stop(); time1 += chrono.usertime();

                chrono.clear(); chrono.start();
                FFPACK::CharPoly(R,g,n,A,n,FFPACK::FfpackLUK);
                chrono.stop(); time2 += chrono.usertime();

                chrono.clear(); chrono.start();
                pass += checker5.check(g);
                chrono.stop(); time1 += chrono.usertime();
            } catch(FailureCharpolyCheck &e) {
                stats_f << " charpoly verification failed! " << std::endl;
            } catch(FailurePLUQCheck &e) {
                stats_f << " internal PLUQ verification failed! " << std::endl;
            }
        }
        time1 /= NR_TESTS;
        time2 /= NR_TESTS;
        stats_f << "     " << i << "-" << i+Range << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2
        << "\t\t" << time1 << endl;
    }


    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C);

    std::cout << "Time: " << global.realtime()
    << " Gfops: " << gffop/global.realtime() << std::endl;

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
