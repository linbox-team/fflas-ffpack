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

// Expect that the follwing macos are defined:
// NSTART
// NMAX
// NSTEP
// DIM
// ITER
// CONFIRM

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

int main (int argc, char** argv) {
    using namespace std;
    using namespace FFPACK;
    typedef Givaro::ModularBalanced<double> Field;
    Field F(2097169);
    Givaro::Poly1Dom<Field> PolDom(F);
    typedef Field::Element Element ;
    size_t n=atoi(argv[1]); // starting value for the block size
    size_t nmax=atoi(argv[2]); // max value for the block size
    size_t step=atoi(argv[3]); // step for the search
    size_t confirm = atoi(argv[4]); // number of increasing values in a row to trigger termination
    size_t dim = atoi(argv[5]); // matrix dimension
    size_t iter = atoi(argv[6]); // number of repetitions
    size_t nbest=0;
    TTimer chrono,tim;

    Element * A = FFLAS::fflas_new (F, dim, dim);
    Element * B = FFLAS::fflas_new (F, dim, dim);
    size_t lda = dim;
    RandomMatrix (F, dim, dim, A, lda);
    FFLAS::fassign (F, dim, dim, A, dim, B, dim);
    typedef typename Givaro::Poly1Dom<Field>::Element Polynomial;
    time_t result = std::time(NULL);
    cerr << std::endl
    << "---------------------------------------------------------------------"
    << std::endl << std::asctime(std::localtime(&result))
    << std::endl
    << "Threshold for ArithProg CharPoly algorithm";
    F.write(cerr << " (using ") << ')' << endl << endl;

    cerr << "Charpoly:  n = "<<dim<<"  degree    seconds            Gfops"<< std::endl;
    double CurrTime;
    double PrevTime = 1000000;
    double BestTime = 1000000;
    size_t increasing = 0;
    Field::RandIter g(F);
    do {
        //warm up computation
        std::list<Polynomial> charp_list;
        try{
            FFPACK::CharPoly (PolDom, charp_list, dim, A, lda, g, FFPACK::FfpackArithProgKrylovPrecond, n);
        }
        catch (CharpolyFailed){}
        FFLAS::fassign (F, dim, dim, B, lda, A, lda);
        chrono.clear();tim.clear();
        for (size_t i=0;i<iter;i++){
            chrono.start();
            try{
                FFPACK::CharPoly (PolDom, charp_list, dim, A, lda, g, FFPACK::FfpackArithProgKrylovPrecond, n);
            }
            catch (CharpolyFailed){	i--;}
            chrono.stop();
            tim+=chrono;
            FFLAS::fassign (F, dim, dim, B, lda, A, lda);
        }
        CurrTime = tim.realtime()/iter;

        cerr << "                           ";
        cerr.width(4);
        cerr << n;
        cerr << "  ";
        cerr.width(15);
        cerr << CurrTime;
        cerr << "  ";
        cerr.width(15);
        cerr << GFOPS(dim,dim,r, CurrTime) << "  "<<std::endl;

        if (BestTime > CurrTime){ // time decreasing -> keep increasing the block size
            nbest = n;
            BestTime = CurrTime;
            increasing=0;
        } else {
            if (CurrTime > PrevTime) increasing ++;
            else increasing = 0;
        }
        PrevTime = CurrTime;
        n += step;
    } while ( (increasing < confirm ) && (n < nmax));

    cerr<<endl;
    if (nbest != 0 ) {
        cout << "#ifndef __FFLASFFPACK_ARITHPROG_THRESHOLD"<< endl;
        cout << "#define __FFLASFFPACK_ARITHPROG_THRESHOLD" << ' ' << nbest << endl;
        cerr << "defined __FFLASFFPACK_ARITHPROG_THRESHOLD to " << nbest << std::endl;
        std::cout << "#endif" << endl  << endl;
    }
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
