/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


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

// #define ENABLE_ALL_CHECKINGS 1 
#include "fflas-ffpack/fflas-ffpack-config.h"
#include <iostream>
#include <givaro/modular.h>

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"


using namespace std;

int main(int argc, char** argv) {
  
	size_t iter = 3;
    bool p=false;
	int    q    = 1009;
	size_t    n    = 2000;
	std::string file = "";
	int t=MAX_THREADS;
	int NBK = -1;
	bool forcecheck = false;
  
	Argument as[] = {
		{ 'q', "-q Q", "Set the field characteristic (-1 for random).",  TYPE_INT , &q },
		{ 'n', "-n N", "Set the dimension of the matrix.",               TYPE_INT , &n },
		{ 'i', "-i R", "Set number of repetitions.",                     TYPE_INT , &iter },
		{ 'f', "-f FILE", "Set the input file (empty for random).",  TYPE_STR , &file },
		{ 't', "-t T", "number of virtual threads to drive the partition.", TYPE_INT , &t },
		{ 'b', "-b B", "number of numa blocks per dimension for the numa placement", TYPE_INT , &NBK },
		{ 'p', "-p P", "multi-threaded.", TYPE_BOOL , &p },
		{ 'c', "-c C", "force checkers.", TYPE_BOOL , &forcecheck },
		END_OF_ARGUMENTS
	};

	if (NBK==-1) NBK = t;
	FFLAS::parseArguments(argc,argv,as);

    typedef Givaro::ModularBalanced<double> Field;
    typedef Field::Element Element;
    
    Field F(q);
    Field::RandIter G(F);
    Field::Element * A;
    
    FFLAS::Timer chrono;
    double time=0.0;
    
    for (size_t i=0;i<iter;++i){
        if (!file.empty()){
            A = read_field(F, file.c_str(),  &n, &n);
        }
        else {
            A = FFLAS::fflas_new<Element>(n*n);
            PAR_BLOCK { FFLAS::pfrand(F,G,n,n,A,n/size_t(NBK)); }	
        }

		FFPACK::ForceCheck_invert<Field> checker(G,n,A,n);
        
        int nullity=0;
        if (p) {
            chrono.clear(); chrono.start();
            FFLAS::ParSeqHelper::Parallel<
                FFLAS::CuttingStrategy::Block,
                FFLAS::StrategyParameter::Threads> PSH(t);

//             FFLAS::ParSeqHelper::Parallel<
//                 FFLAS::CuttingStrategy::Recursive,
//                 FFLAS::StrategyParameter::TwoDAdaptive> PSH(t);

            PAR_BLOCK { FFPACK::Invert (F, n, A, n, nullity, PSH); }
            chrono.stop();
        } else {
            chrono.clear(); chrono.start();
            FFPACK::Invert (F, n, A, n, nullity);
            chrono.stop();
        }

        if (forcecheck) checker.check(A,nullity);

        time+=chrono.realtime();
        FFLAS::fflas_delete( A);
    }
    
        // -----------
        // Standard output for benchmark - Alexis Breust 2014/11/14
#define CUBE(x) ((x)*(x)*(x))
	std::cout << "Time: " << time / double(iter)
			  << " Gflops: " << 2. * CUBE(double(n)/1000.) / time * double(iter);
	FFLAS::writeCommandString(std::cout, as) << std::endl;
    

  return 0;
}




