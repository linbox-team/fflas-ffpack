/*
 * Copyright (C) 2012 FFLAS-FFPACK group.
 *
 * Extirpé form a m4 macro by Brice Boyer (briceboyer) <boyer.brice@gmail.com>.
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


//#define LinBoxSrcOnly
#define DOUBLE_TO_FLOAT_CROSSOVER 0

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include <iostream>
#include <fstream>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"

template<class Field>
bool balanced(const Field & )
{
    return false;
}

template <class T>
bool balanced(const Givaro::ModularBalanced<T>&)
{
    return true;
}

#ifdef __GIVARO_USE_OPENMP
typedef Givaro::OMPTimer TTimer;
#else
typedef Givaro::Timer TTimer;
#endif

#define GFOPS(n,t) (2.0/t*(double)n/1000.0*(double)n/1000.0*(double)n/1000.0)

#include <ctime>

int main () {
    using namespace std;

    typedef FIELD Field;
    Field F(17);
    typedef Field::Element Element ;
    size_t n=512, nmax=4000, prec=512, nbest=0, count=0;
    TTimer chrono;
    bool bound=false;

    Element * A = FFLAS::fflas_new (F,nmax,nmax);
    Element * B = FFLAS::fflas_new (F,nmax,nmax);
    Element * C = FFLAS::fflas_new (F,nmax,nmax);
    FFPACK::RandomMatrix (F, nmax, nmax, A, nmax);
    FFPACK::RandomMatrix (F, nmax, nmax, B, nmax);
    FFPACK::RandomMatrix (F, nmax, nmax, C, nmax);

    time_t result = std::time(NULL);
    cerr << std::endl
    << "---------------------------------------------------------------------"
    << std::endl << std::asctime(std::localtime(&result))
    << std::endl
    << "Threshold for finite field Strassen-Winograd matrix multiplication" ;
    F.write(cerr << " (using ") << ')' << endl << endl;

    cerr << "fgemm:  n                   Classic                        Winograd 1 level" << std::endl;
    cerr << "                    seconds            Gfops          seconds            Gfops" << std::endl;
    FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> ClassicH(F,0, FFLAS::ParSeqHelper::Sequential());
    FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WinogradH(F,1, FFLAS::ParSeqHelper::Sequential());
    //warm up computation
    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, F.mOne, A, n, B, n, F.one, C, n, ClassicH);
    do {
        double classicTime, winogradTime;

        int iter=3;
        chrono.start();
        for (int i=0;i<iter;i++)
            FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, F.mOne, A, n, B, n, F.one, C, n, ClassicH);
        chrono.stop();

        classicTime = chrono.realtime()/iter;

        chrono.clear(); chrono.start();
        for (int i=0; i<iter; i++)
            FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n, F.mOne, A, n, B,n, F.one, C, n, WinogradH);
        chrono.stop();

        winogradTime = chrono.realtime()/iter;

        cerr << "      ";
        cerr.width(4);
        cerr << n;
        cerr << "  ";
        cerr.width(15);
        cerr << classicTime;
        cerr << "  ";
        cerr.width(15);
        cerr << GFOPS(n, classicTime) << "  ";
        cerr.width(15);
        cerr << winogradTime;
        cerr << "  ";
        cerr.width(15);
        cerr << GFOPS(n, winogradTime) << endl;

        if (classicTime > winogradTime ){
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
    } while ((prec > 32 ) && (n < nmax));

    cerr<<endl;
    if (nbest != 0 ) {
        if (typeid(Element).name() == typeid(double).name()) {
            if ( balanced(F) ) {
                cout << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL"  << endl;
                cout << "#define __FFLASFFPACK_WINOTHRESHOLD_BAL" << ' ' <<  nbest << endl;
                cerr << "defined __FFLASFFPACK_WINOTHRESHOLD_BAL to " << nbest << "" << std::endl;
            }
            else {
                cout << "#ifndef __FFLASFFPACK_WINOTHRESHOLD"  << endl;
                cout << "#define __FFLASFFPACK_WINOTHRESHOLD" << ' ' <<  nbest << endl;
                cerr << "defined __FFLASFFPACK_WINOTHRESHOLD to " << nbest << "" << std::endl;

            }
            std::cout << "#endif" << endl  << endl;
        }

        if (typeid(Element).name() == typeid(float).name()) {
            if ( balanced(F) ) {
                cout << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT"  << endl;
                cout << "#define __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT" << ' ' << nbest << endl;
                cerr << "defined __FFLASFFPACK_WINOTHRESHOLD_BAL_FLT to " << nbest << "" << std::endl;

            }
            else {
                cout << "#ifndef __FFLASFFPACK_WINOTHRESHOLD_FLT"  << endl;
                cout << "#define __FFLASFFPACK_WINOTHRESHOLD_FLT" << ' ' << nbest << endl;
                cerr << "defined __FFLASFFPACK_WINOTHRESHOLD_FLT to " << nbest << "" << std::endl;
            }
            cout << "#endif" << endl << endl;
        }
    }

    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C);

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
