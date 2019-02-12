/*
 * Copyright (C) 2016 FFLAS-FFPACK group.
 *
 * Written by Clément Pernet <clement.pernet@imag.fr>
 *
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



#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include <iostream>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/ffpack/ffpack.h"


#ifdef __GIVARO_USE_OPENMP
typedef Givaro::OMPTimer TTimer;
#else
typedef Givaro::Timer TTimer;
#endif

#include <ctime>
#define CUBE(x) ((x)*(x)*(x))
#define GFOPS(m,n,r,t) (2.0/3.0*CUBE(double(n)/1000.0) +2*m/1000.0*n/1000.0*double(r)/1000.0  - double(r)/1000.0*double(r)/1000.0*(m+n)/1000)/t

int main () {
    using namespace std;

    typedef Givaro::ModularBalanced<double> Field;
    Field F(131071);
    typedef Field::Element Element ;
    size_t n=128, nmax=1000, prec=64, nbest=0, count=0;
    TTimer chrono,tim;
    bool bound=false;

    Element * A = FFLAS::fflas_new (F, nmax, nmax);
    Element * B = FFLAS::fflas_new (F, nmax, nmax);
    size_t lda = nmax;
    size_t r;
    size_t * P = new size_t[nmax];
    size_t * Q = new size_t[nmax];
    FFPACK::RandomMatrix (F, nmax, nmax, A,nmax);
    FFLAS::fassign (F, nmax,nmax, A, nmax, B, nmax);
    time_t result = std::time(NULL);
    cerr << std::endl
    << "---------------------------------------------------------------------"
    << std::endl << std::asctime(std::localtime(&result))
    << std::endl
    << "Threshold for PLUQ base case" ;
    F.write(cerr << " (using ") << ')' << endl << endl;

    cerr << "PLUQ:  n                   Base case                        Recursive 1 level" << std::endl;
    cerr << "                    seconds            Gfops          seconds            Gfops" << std::endl;
    do {
        double BCTime, RecTime;
        int iter=10;

        //warm up computation
        r=FFPACK::PLUQ_basecaseCrout (F, FFLAS::FflasNonUnit, n, n, A, lda, P, Q);
        FFLAS::fassign (F, n, n, B, lda, A, lda);
        chrono.clear();tim.clear();
        for (int i=0;i<iter;i++){
            chrono.start();
            r=FFPACK::PLUQ_basecaseCrout (F, FFLAS::FflasNonUnit, n, n, A, lda, P, Q);
            chrono.stop();
            tim+=chrono;
            FFLAS::fassign (F, n, n, B, lda, A, lda);
        }
        BCTime = tim.realtime()/iter;

        tim.clear();chrono.clear();
        for (int i=0;i<iter;i++){
            chrono.start();
            r=FFPACK::_PLUQ (F, FFLAS::FflasNonUnit, n, n, A, lda, P, Q, n);
            chrono.stop();
            tim+=chrono;
            FFLAS::fassign (F, n, n, B, lda, A, lda);
        }
        //FFLAS::fflas_delete(A);
        //FFLAS::fflas_delete(B);
        RecTime = tim.realtime()/iter;

        cerr << "      ";
        cerr.width(4);
        cerr << n;
        cerr << "  ";
        cerr.width(15);
        cerr << BCTime;
        cerr << "  ";
        cerr.width(15);
        cerr << GFOPS(n,n,r, BCTime) << "  ";
        cerr.width(15);
        cerr << RecTime;
        cerr << "  ";
        cerr.width(15);
        cerr << GFOPS(n,n,r, RecTime) << endl;

        if (BCTime > RecTime){
            count++;
            if (count > 1){
                nbest = n;
                bound = true;
                prec = prec >> 1;
                n -= prec;
            }
        }
        else{
            count=0;
            if (bound)
                prec=prec>>1;
            n+=prec;
        }
    } while ((prec > 4 ) && (n < nmax));

    cerr<<endl;
    if (nbest != 0 ) {
        cout << "#ifndef __FFLASFFPACK_PLUQ_THRESHOLD"  << endl;
        cout << "#define __FFLASFFPACK_PLUQ_THRESHOLD" << ' ' <<  nbest << endl;
        cerr << "defined __FFLASFFPACK_PLUQ_THRESHOLD to " << nbest << "" << std::endl;
        std::cout << "#endif" << endl  << endl;
    }
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
