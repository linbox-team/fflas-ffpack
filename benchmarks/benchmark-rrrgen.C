/* Copyright (c) FFLAS-FFPACK
 * Written by Hippolyte Signargout <hippolyte.signargout@ens-lyon.fr>
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

// Template from benchmark-quasisep.C

#define __FFLASFFPACK_OPENBLAS_NT_ALREADY_SET 1

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>
#include <fstream>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/fflas_randommatrix.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/ffpack/ffpack_rrrgen.inl"



using namespace std;
using namespace FFLAS;
using namespace FFPACK;




template<class Field>
void run_with_field(int q, size_t n, size_t m, size_t t, size_t r, size_t iter, uint64_t seed, fstream& my_file){

    Field F(q);
    typedef typename Field::Element_ptr Element_ptr;
    typename Field::RandIter G (F, seed);
    FFLAS::Timer chrono;
    
    size_t lda = n;
    size_t ldts = m;


    Element_ptr A, B, TS;
    A = FFLAS::fflas_new (F, n, n);
    B = FFLAS::fflas_new (F, n, n);
    TS = FFLAS::fflas_new (F, n, m);
    Element_ptr Res = fflas_new(F, n, m); // Inadequate name
    RRRgen<Field>* RRRA;
    RRRgen<Field>* RRRB;
    RRRgen<Field>* RRRres;


    Element_ptr A2 = fflas_new (F, n, n);
    Element_ptr B2 = fflas_new (F, n, n);


    double time_RRRxTS = 0, time_RRRxRRR = 0,time_gen_qs = 0, time_gen_rrr = 0;
    size_t * p = FFLAS::fflas_new<size_t> (ceil(n/2.));
    for (size_t i = 0; i < ceil(n/2.); i++)
    {
        p[i] = n - i - 1;
    }
    cout << "Temps de calcul pour les opérations RRRGen, RRRxTS, RRRxRRR et avec n : " << n << endl;
    for (size_t i=0; i<iter;++i){

        cout << "Generation start"<< endl;
        // generates random matrices
        typename Field::RandIter G (F, seed);
        RandomMatrix(F, n, m, TS, ldts, G);
        
        // generate A a t qsmatrix
        
        chrono.clear();
        chrono.start(); 

        RandomLTQSMatrixWithRankandQSorder (F,n,r,t,A, n,G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t, A2, n,G);
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A2, n, p);
        faddin (F, n, n, A2, n, A, n);
        chrono.stop();
        time_gen_qs = chrono.usertime();
        
        cout << "Generation of A ended in "<< time_gen_qs << endl;

        // generate A a t qsmatrix
        chrono.clear();
        chrono.start(); 

        RandomLTQSMatrixWithRankandQSorder (F,n,r,t, B, n,G);
        RandomLTQSMatrixWithRankandQSorder (F,n,r,t, B2, n,G);
        applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B, n, p);
        applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), B2, n, p);
        faddin (F, n, n, B2, n, B, n);
        chrono.stop();
        time_gen_qs += chrono.usertime();
        
        cout << "Generation of B ended in "<< time_gen_qs << endl;

        // create RRR
        chrono.clear();
        chrono.start(); 
        RRRA = new RRRgen<Field>(F, n, t, A, lda,true,true);
        chrono.stop();
        time_gen_rrr += chrono.usertime();

        // create RRR
        chrono.clear();
        chrono.start(); 
        RRRB = new RRRgen<Field>(F, n, t, B, lda,true,true);
        chrono.stop();
        time_gen_rrr += chrono.usertime();

        // RRRxTS product
        chrono.clear();
        chrono.start(); 
        RRRxTS(F,n,m,RRRA,TS,m, Res,m);
        chrono.stop();
        time_RRRxTS += chrono.usertime();
        
        // RRRxRRR
        chrono.clear();
        chrono.start(); 
        RRRres = RRRxRRR(F,RRRA,RRRB);
        chrono.stop();
        time_RRRxRRR += chrono.usertime();
        delete RRRres;
        
        // // Could add RRR invert when Matrix with QSorder and full rank can be created. Here some matrix can't be inverted then raise error.
        // chrono.clear();
        // chrono.start(); 
        // RRRres = RRRinvert(F,RRRA);
        // // Can add an other inverse operation of B to be more precise on time measure. 
        // chrono.stop();
        // time_invert += chrono.usertime();
        // delete RRRres;

        delete RRRA;
        delete RRRB;
    }
    
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(TS);
    FFLAS::fflas_delete(Res); 
    FFLAS::fflas_delete(A2,B2, p);
    
    double mean_time_RRRxTS = time_RRRxTS / double(iter);
    double mean_time_gen_RRR = time_gen_rrr / double(2*iter);
    double mean_time_RRRxRRR = time_RRRxRRR / double(iter);

    cout << "n : " << n << ", t : " << t << ", RRR gen  : " << mean_time_gen_RRR << ", RRRxRRR : " << mean_time_RRRxRRR << ", RRRxTS : " << mean_time_RRRxTS <<endl;
    my_file << "n : " << n << ", t : " << t << ", RRR gen  : " << mean_time_gen_RRR << ", RRRxRRR : " << mean_time_RRRxRRR << ", RRRxTS : " << mean_time_RRRxTS <<endl;
}

int main(int argc, char** argv) {

#ifdef __FFLASFFPACK_OPENBLAS_NUM_THREADS
    openblas_set_num_threads(__FFLASFFPACK_OPENBLAS_NUM_THREADS);
#endif

    size_t iter = 10;
    int    q    = 131071;
    size_t    n    = 2167;
    size_t    m    = 455;
    size_t    t    = 236;
    size_t    r    = 1100;
    uint64_t seed = FFLAS::getSeed();

    Argument as[] = {
                     { 'q', "-q Q", "Set the field characteristic (-1 for the ring ZZ).",     TYPE_INT , &q },
                     { 'n', "-n N", "Set the order of the square matrix A.",                  TYPE_INT , &n },
                     { 'm', "-m M", "Set the column dimension of n x m RHS matrix B.",        TYPE_INT , &m },
                     { 't', "-t T", "Set the quasiseparability order of A.",                  TYPE_INT , &t },
                     { 'r', "-r R", "Set the rank of each upper/lower triangular part of A.", TYPE_INT , &r },
                     { 'i', "-i R", "Set number of repetitions.",                             TYPE_INT , &iter },
                     { 's', "-s S", "Sets seed.",                                             TYPE_INT , &seed },
                     END_OF_ARGUMENTS
    };

    fstream my_file;
	my_file.open("fin RRR measures ", ios::out);
	if (!my_file) {
		cout << "File not created!"<< std::endl;
	}
	else {
		cout << "File created successfully!"<< std::endl;
	}

    // FFLAS::parseArguments(argc,argv,as);
    // for (n = 500; n<6000;n = n+100){
        
    //     run_with_field<Givaro::ModularBalanced<double> >(q, n, m, n/9, n/2, iter, seed,my_file);
    // }
    for (n = 9000; n<=10001;n = n+1000){
        
        run_with_field<Givaro::ModularBalanced<double> >(q, n, m, n/9, n/2, iter, seed,my_file);
    }
    my_file << endl;

    my_file << "t : " << t << ", m : " << m << endl;

    my_file << endl;
    my_file << endl;
    my_file << endl;

    // n = 2167;
    // for (t = 50; t < r; t+=10){
    //     run_with_field<Givaro::ModularBalanced<double> >(q, n, m, t, r, iter, seed,my_file);
    // }
    // my_file << endl;

    // my_file << "n : " << n << ", m : " << m << ", r : " << r << endl;

    my_file.close();
    // std::cout << "( ";
    // FFLAS::writeCommandString(std::cout, as) << ")" << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
