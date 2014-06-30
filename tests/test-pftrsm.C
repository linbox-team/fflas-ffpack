/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Ziad Sultan
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
 *.
 */

//--------------------------------------------------------------------------
//                        Test for pftrsm 
//
//--------------------------------------------------------------------------
// Ziad Sultan
//-------------------------------------------------------------------------

#define DEBUG 1
#define TIME 1

#define __FFLASFFPACK_USE_OPENMP
//#define __FFLASFFPACK_USE_KAAPI

#define __FFLAS__TRSM_READONLY

#include <iomanip>
#include <iostream>
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "time.h"


using namespace std;
using namespace FFPACK;

typedef Modular<double> Field;
//typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;




//int main(int argc, char** argv)
//{
BEGIN_PARALLEL_MAIN(int argc, char** argv)
{

	int k,n,m;
	cerr<<setprecision(10);
	Field::Element zero, one;

	if (argc != 10)	{
		cerr<<"Usage : test-ftrsm <p> <A> <B> <iter> <alpha> <left/right> <Up/Low> <NoTrans/Trans> <NonUnit/Unit>"
		    <<endl;
		exit(-1);
	}
	int nbit=atoi(argv[4]); // number of times the product is performed
	Field F(atof(argv[1]));
	F.init(zero,0.0);
	F.init(one,1.0);
	Field::Element * A, *B, *B2;
	A = read_field(F,argv[2],&k,&k);
	B = read_field(F,argv[3],&m,&n);
	B2 = new Field::Element[m*n];


	for (int i=0; i<m;++i){
		for(int j=0; j<n; ++j)
			F.assign(*(B2+i*n+j),*(B+i*n+j));
	}

	Field::Element alpha;
	F.init (alpha, atof(argv[5]));

	FFLAS::FFLAS_SIDE side = (atoi(argv[6])) ? FFLAS::FflasRight :  FFLAS::FflasLeft;
	FFLAS::FFLAS_UPLO uplo = (atoi(argv[7])) ? FFLAS::FflasLower :  FFLAS::FflasUpper;
	FFLAS::FFLAS_TRANSPOSE trans = (atoi(argv[8])) ? FFLAS::FflasTrans :  FFLAS::FflasNoTrans;
	FFLAS::FFLAS_DIAG diag = (atoi(argv[9])) ? FFLAS::FflasUnit :  FFLAS::FflasNonUnit;

	if (   ((side == FFLAS::FflasRight) &&(k != n))
	    || ((side == FFLAS::FflasLeft)&&(k != m))) {
		cerr<<"Error in the dimensions of the input matrices"<<endl;
		exit(-1);
	}

	struct timespec t0,t1;
        double delay, avrg;
        double t_total=0;

	//	Timer t; t.clear();
	//	double time=0.0;
// write_field(F, cerr<<"A="<<endl, A, k,k,k);
// write_field(F, cerr<<"B="<<endl, B, m,n,n);
        const FFLAS::CuttingStrategy Strategy = FFLAS::BLOCK_THREADS;
	//	size_t numThreads;
	for(int i = 0;i<nbit;++i){
		//		t.clear();
		//		t.start();		
		clock_gettime(CLOCK_REALTIME, &t0);
		HPAC_PAR_REGION{
			//		numThreads=HPAC_NUM_THREADS;
			FFLAS::TRSMHelper<FFLAS::StructureHelper::Iterative,
					  FFLAS::ParSeqHelper::Parallel> PH (FFLAS::ParSeqHelper::Parallel(omp_get_max_threads(),Strategy));	
			FFLAS::pftrsm (F, side, uplo, trans, diag, m, n, alpha, A, k, B, n, PH);
		}
		BARRIER;
		clock_gettime(CLOCK_REALTIME, &t1);
		delay = (double)(t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)/1000000000;

		if (i)
                        t_total+=delay;
		//		t.stop();
		//		time+=t.usertime();
		if (i+1<nbit)
                        for (int i=0; i<m*n;++i)
                                F.assign(*(B+i),*(B2+i));

	}
        avrg = t_total/(nbit-1);

#if DEBUG
// write_field (F,cerr<<"S="<<endl,B,m,n,n);
	Field::Element invalpha;
	F.inv(invalpha, alpha);

	FFLAS::ftrmm (F, side, uplo, trans, diag, m, n, invalpha, A, k, B, n);
	bool wrong = false;

	for (int i=0;i<m;++i)
		for (int j=0;j<n;++j)
			if ( !F.areEqual(*(B2+i*n+j), *(B+i*n+j))){
				cerr<<"B2 ["<<i<<", "<<j<<"] = "<<(*(B2+i*n+j))
				    <<" ; B ["<<i<<", "<<j<<"] = "<<(*(B+i*n+j))
				    <<endl;
				wrong = true;
			}

	if ( wrong ){
		cerr<<"FAIL"<<endl;
// write_field (F,cerr<<"B2="<<endl,B2,m,n,n);
// write_field (F,cerr<<"B="<<endl,B,m,n,n);
	} else
		cerr<<"PASS"<<endl;
#endif

	delete[] A;
	delete[] B;
	delete[] B2;

#if TIME
	//	double mflops = m*n/1000000.0*nbit*n/time;
	double mflops = m*n/1000000.0*n/avrg;
	cerr<<"m,n = "<<m<<" "<<n<<". ftrsm "
	    <<((side == FFLAS::FflasLeft)?" Left ":" Right ")
	    <<((uplo == FFLAS::FflasLower)?" Lower ":" Upper ")
	    <<((diag == FFLAS::FflasUnit)?" Unit ":" NonUnit ")
	    <<((trans == FFLAS::FflasTrans)?" Trans ":" NoTrans ")
	    <<"over Z/"<<atoi(argv[1])<<"Z :"
	    <<endl
	    <<"t= "
	    << avrg
	    << " s, Mffops = "<<mflops
	    << endl;

	//	cout<<m<<" "<<n<<" "<<mflops<<" "<<time/nbit<<endl;
	cout<<m<<" "<<n<<" "<<mflops<<" "<<avrg<<endl;
#endif
	//	 }

	/*
	size_t BLOCKSIZE=std::max(m/numThreads,(size_t)1); // There is always 2 TRSM taking place in parallel                                                        
	size_t NBlocks = m/BLOCKSIZE;
	size_t LastBlockSize = m % BLOCKSIZE;
	if (LastBlockSize)
		NBlocks++;
	else
		LastBlockSize=BLOCKSIZE;
	//#pragma omp parallel for default (none) shared(A,B,F,NBlocks, LastBlockSize, BLOCKSIZE)                                                                      
	for (size_t t = 0; t < NBlocks; ++t){
		size_t i = t % NBlocks;
		size_t BlockDim = BLOCKSIZE;
		if (i == NBlocks-1)
			BlockDim = LastBlockSize;


		cout<<"blockDim "<<BlockDim<<"  BLOCKSIZE*i*ldb "<< (BLOCKSIZE*i*n)<<endl;
	}
	FFLAS::ForStrategy1D iter(m, FFLAS::BLOCK_THREADS, numThreads);
	for (iter.begin(); ! iter.end(); ++iter)
		cout<<"iter.iend-iter.ibeg "<<iter.iend-iter.ibeg<<"   iter.ibeg*ldb "<<  (iter.ibeg*n)<<endl;
	*/
}
END_PARALLEL_MAIN()
