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

// Expect that the following macos are defined:
// VARIANT1
// VARIANT2
// NSTART
// NFIRSTSTEP
// NMAX
// NPREC
// ITER

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
#define GFOPS(m,n,r,t) (2.7*CUBE(double(n)/1000.0))/t // approximative flop count

int main () {
    using namespace std;
    using namespace FFPACK;
    typedef Givaro::ModularBalanced<double> Field;
    Field F(4194301);
    typedef Field::Element Element ;
    size_t n=NSTART, nmax=NMAX, prec=NFIRSTSTEP, nbest=0, count=0;
    TTimer chrono,tim;
    bool bound=false;

    Element * A = FFLAS::fflas_new (F, nmax, nmax);
    Element * B = FFLAS::fflas_new (F, nmax, nmax);
    size_t lda = nmax;
    RandomMatrix (F, nmax, nmax, A, nmax);
    FFLAS::fassign (F, nmax,nmax, A, nmax, B, nmax);
    typedef typename Givaro::Poly1Dom<Field>::Element Polynomial;
    time_t result = std::time(NULL);
    string var1 = (VARIANT1 == FfpackLUK)?"LUKrylov":"Danilevskii";
    string var2 = (VARIANT2 == FfpackLUK)?"LUKrylov":"ArithProg";
    cerr << std::endl
    << "---------------------------------------------------------------------"
    << std::endl << std::asctime(std::localtime(&result))
    << std::endl
    << "Threshold between CharPoly algorithms "<<var1<<" <-> "<<var2 ;
    F.write(cerr << " (using ") << ')' << endl << endl;

    cerr << "Charpoly:  n                   "<<var1<<"                        "<<var2 << std::endl;
    cerr << "                    seconds            Gfops          seconds            Gfops" << std::endl;
    do {
        double Var1Time, Var2Time;
        int iter=ITER;
        //warm up computation
        Polynomial charp(n+1);
        Givaro::Poly1Dom<Field> PolDom(F);
        CharPoly (PolDom, charp, n, A, lda, VARIANT1);
        FFLAS::fassign (F, n, n, B, lda, A, lda);
        chrono.clear();tim.clear();
        for (int i=0;i<iter;i++){
            chrono.start();
            CharPoly (PolDom, charp, n, A, lda, VARIANT1);
            chrono.stop();
            tim+=chrono;
            FFLAS::fassign (F, n, n, B, lda, A, lda);
        }
        Var1Time = tim.realtime()/iter;

        tim.clear();chrono.clear();
        for (int i=0;i<iter;i++){
            chrono.start();
            CharPoly (PolDom, charp, n, A, lda, VARIANT2);
            chrono.stop();
            tim+=chrono;
            FFLAS::fassign (F, n, n, B, lda, A, lda);
        }
        //FFLAS::fflas_delete(A);
        //FFLAS::fflas_delete(B);
        Var2Time = tim.realtime()/iter;

        cerr << "      ";
        cerr.width(4);
        cerr << n;
        cerr << "  ";
        cerr.width(15);
        cerr << Var1Time;
        cerr << "  ";
        cerr.width(15);
        cerr << GFOPS(n,n,r, Var1Time) << "  ";
        cerr.width(15);
        cerr << Var2Time;
        cerr << "  ";
        cerr.width(15);
        cerr << GFOPS(n,n,r, Var2Time) << endl;

        if (Var1Time > Var2Time){
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
    } while ((prec > NPREC ) && (n < nmax));

    cerr<<endl;
    if (nbest != 0 ) {
        cout << "#ifndef __FFLASFFPACK_CHARPOLY_"<<var1<<"_"<<var2<<"_THRESHOLD"  << endl;
        cout << "#define __FFLASFFPACK_CHARPOLY_"<<var1<<"_"<<var2<<"_THRESHOLD" << ' ' <<  nbest << endl;
        cerr << "defined __FFLASFFPACK_CHARPOLY_"<<var1<<"_"<<var2<<"_THRESHOLD to " << nbest << "" << std::endl;
        std::cout << "#endif" << endl  << endl;
    }
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
