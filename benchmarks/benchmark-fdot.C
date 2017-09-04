/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* Copyright (c) FFLAS-FFPACK
 * Written by Philippe LEDENT <philippe.ledent@etu.univ-grenoble-alpes.fr>
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

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>
#include <givaro/givrational.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/paladin/parallel.h"
#include "fflas-ffpack/paladin/fflas_plevel1.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;
using namespace FFLAS;
using namespace FFPACK;

template<class Field>
typename Field::Element run_with_field(int q, size_t iter, size_t N, const size_t BS, const size_t p, const size_t threads){
    Field F(q);
    typename Field::RandIter G(F, BS);

	typename Field::Element_ptr A, B;
	typename Field::Element d; F.init(d);

#ifdef__GIVARO_USE_OPENMP
	Givaro::OMPTimer chrono, time;
#else
	Givaro::Timer chrono, time;
#endif
	time.clear();
  
	for (size_t i=0;i<iter;++i){
		A = fflas_new(F, N);
		B = fflas_new(F, N);

        PAR_BLOCK { pfrand(F, G, N, 1, A); pfrand(F, G, N, 1, B); }

//         FFLAS::WriteMatrix(std::cerr, F, 1, N, A, 1);
//         FFLAS::WriteMatrix(std::cerr, F, 1, N, B, 1);
        
        F.assign(d, F.zero);


        FFLAS::ParSeqHelper::Parallel<
            FFLAS::CuttingStrategy::Block,
            FFLAS::StrategyParameter::Threads> ParHelper(threads);

		chrono.clear();
        if (p){
            chrono.start();
            F.assign(d, fdot(F, N, A, 1U, B, 1U, ParHelper));
            chrono.stop();
        } else {
            chrono.start();
            F.assign(d, fdot(F, N, A, 1U, B, 1U, FFLAS::ParSeqHelper::Sequential()));
            chrono.stop();
        }


        std::cerr << chrono 
                  << " Gfops: " << ((double(2*N)/1000.)/1000.)/(1000.*chrono.realtime())
                  << std::endl;

		time+=chrono;
		FFLAS::fflas_delete(A);
		FFLAS::fflas_delete(B);
	}
		// -----------
		// Standard output for benchmark
	std::cout << "Time * " << iter << ": " << time 
			  << " Gfops: " << ((double(2*N)/1000.)/1000.)/(1000.*time.realtime())* double(iter);

// 	F.write(std::cerr, d) << std::endl;
    return d;
}

int main(int argc, char** argv) {

	size_t iter = 3;
	size_t N    = 5000;
    size_t BS   = 5000;
    int q		= 131071101;
	size_t p	=0;
    size_t maxallowed_threads; PAR_BLOCK { maxallowed_threads=NUM_THREADS; } 
    size_t threads=maxallowed_threads;

	Argument as[] = {
		{ 'n', "-n N", "Set the dimension of the matrix C.",TYPE_INT , &N },
        { 'q', "-q Q", "Set the field characteristic (0 for the integers).",         TYPE_INT , &q },
		{ 'i', "-i R", "Set number of repetitions.",		TYPE_INT , &iter },
		{ 'b', "-b B", "Set the bitsize of the random elements.",         TYPE_INT , &BS},
		{ 'p', "-p P", "0 for sequential, 1 for parallel.",	TYPE_INT , &p },
		{ 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &threads },
		END_OF_ARGUMENTS
	};

    FFLAS::parseArguments(argc,argv,as);

    if (q > 0){
        BS = Givaro::Integer(q).bitsize();
        double d = run_with_field<Givaro::ModularBalanced<double> >(q, iter, N, BS, p, threads);
        std::cout << ", d: " << d;
    } else {
        auto d = run_with_field<Givaro::ZRing<Givaro::Integer> > (q, iter, N, BS, p, threads);
        std::cout << ", size: " << logtwo(d>0?d:-d);
    }

    FFLAS::writeCommandString(std::cerr, as) << std::endl;
    return 0;
}
