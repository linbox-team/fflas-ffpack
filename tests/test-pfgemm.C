/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet
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
//                        Test for fgemm : 1 computation
//
//--------------------------------------------------------------------------
// Clement Pernet, Ziad sultan
//-------------------------------------------------------------------------

#define NEWWINO
#ifndef TIME
#define TIME 1
#endif


//#define DEBUG
#include <iomanip>
#include <iostream>
using namespace std;

#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/utils/args-parser.h"

using namespace FFPACK;

template<class T>
T& myrand (T& r, long size) {
    if (size < 0)
        return r = T( (lrand48() % (-size-size)) + size );
    else
        return r = T(  lrand48() % size ) ;
}

typedef Modular<double> Field;
//typedef Modular<float> Field;
//typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;
//typedef Modular<int> Field;

//#ifdef  __FFLASFFPACK_USE_OPENMP
Field::Element* makemat(const Field::RandIter& RF,int m, int n){
	Field::Element * res = new Field::Element[m*n];
#pragma omp parallel for
	for (long i = 0; i < m; ++i)
		for (long j = 0; j < n; ++j) {
			RF.random(res[j+i*n]);
		}
	return res;
}
//#endif


int main(int argc, char** argv){

	//#ifdef  __FFLASFFPACK_USE_OPENMP
	/*
	if (argc < 8 || argc > 9)	{
		cerr<<"Usage : test-fgemm <p> <nA> <nB> <w> <i> <alpha> <beta> [cutting]"
		    <<endl;
		exit(-1);
	}


	
	int m=atoi(argv[2]),n=atoi(argv[3]),k=m;
	int nbw=atoi(argv[4]); // number of winograd levels
	int pnbw = nbw;
	int nbit=atoi(argv[5]); // number of times the product is performed
	*/

	srand((int)time(NULL));
        srand48(time(NULL));

	size_t m,n,k,p, nbw, nbit, iters;
	iters=5;
	p=65521;
	n=100;
	nbw=2;
	n=20+(size_t)random()%n;
	static Argument as[] = {
                { 'p', "-p P", "Set the field characteristic.",         TYPE_INT , &p },
                { 'n', "-n N", "Set the dimension of the matrix.",      TYPE_INT , &n },
                { 'w', "-w N", "Set the number of winograd levels.",    TYPE_INT , &nbw },
                { 'i', "-i R", "Set number of repetitions.",            TYPE_INT , &iters },
                END_OF_ARGUMENTS
        };

	FFLAS::parseArguments(argc,argv,as);

	m = n ;
        k = m ;
        size_t lda = k;
        size_t ldb = n;
	//        size_t ldc = n;
	nbit = iters;
	//        bool ok = true ;

        srand48(BaseTimer::seed());

	cerr<<setprecision(10);
	Field::Element alpha,beta;


	Field F(p);

	F.init( alpha, Field::Element(1.0));
	F.init( beta, Field::Element(0.0));

	const FFLAS::CuttingStrategy Strategy = FFLAS::BLOCK_THREADS;

	Field::RandIter RF(F);

	Field::Element * A = makemat(RF,m,k);
	Field::Element * B = makemat(RF,k,n);


	enum FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans;
	enum FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans;


	Field::Element * C=NULL;

// 	write_field (F, cerr<<"A = "<<endl, A, m, k, lda);
// 	write_field (F, cerr<<"B = "<<endl, B, k, n, ldb);

    size_t r,c; FFLAS::BlockCuts(r,c,m,n,Strategy, omp_get_max_threads() );
    //    std::cerr << "pfgemm: " << m << 'x' << n << ' ' << r << ':' << c << "  <--  " << omp_get_max_threads() << ':' << (m/r) << 'x' << (n/c) << std::endl;
//     if (nbw <0) {
//         FFLAS::FFLAS_BASE base; size_t winolev, kmax;
//         FFLAS::Protected::MatMulParameters (F, m, n, k, beta, kmax, base, winolev);
//         nbw=winolev;
//         pnbw=0;
//         std::cerr << "Winolevel: " << nbw << '(' << pnbw << ')' << std::endl;
//     }
    FFLAS::MMHelper<Field,
	    FFLAS::MMHelperAlgo::Winograd, 
	    FFLAS::FieldTraits<Field>::value,
	    FFLAS::ParSeqHelper::Parallel> pWH (F,nbw,FFLAS::ParSeqHelper::Parallel(omp_get_max_threads(),Strategy));	

	OMPTimer tim,t; t.clear();tim.clear();
	for(size_t i = 0;i<nbit+1;++i){
        C = new Field::Element[m*n];
		t.clear();
		t.start();
        PAR_REGION{
            
                TASK(READ(A, B),NOWRITE(), READWRITE(C), fgemm, F, ta, tb,m,n,k,alpha, A,lda, B,ldb,
                               beta,C,n, pWH );
            
        }
        t.stop();
		if (i) tim+=t;
//         if (i<nbit) delete[] C;
	}

	//#ifdef DEBUG
	//    cerr<<"Debugging ... ";
	bool wrong = false;
	Field::Element zero;
	F.init(zero, 0.0);
	Field::Element * Cd;
		Cd  = new Field::Element[m*n];
		for (size_t i=0; i<m*n; ++i)
			F.assign (*(Cd+i), zero);

	Field::Element aij, bij,  tmp;
	// Field::Element beta_alpha;
	//F.div (beta_alpha, beta, alpha);
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j){
			F.mulin(*(Cd+i*n+j),beta);
			F.assign (tmp, zero);
			for (size_t l = 0; l < k ; ++l ){
				if ( ta == FFLAS::FflasNoTrans )
					aij = *(A+i*lda+l);
				else
					aij = *(A+l*lda+i);
				if ( tb == FFLAS::FflasNoTrans )
					bij = *(B+l*ldb+j);
				else
					bij = *(B+j*ldb+l);
				//F.mul (tmp, aij, bij);
				//F.axpyin( *(Cd+i*n+j), alpha, tmp );
				F.axpyin (tmp, aij, bij);
			}
			F.axpyin (*(Cd+i*n+j), alpha, tmp);
			//F.mulin( *(Cd+i*n+j),alpha );
			if ( !F.areEqual( *(Cd+i*n+j), *(C+i*n+j) ) ) {
				wrong = true;
			}
		}
	if ( wrong ){
		cerr<<"FAIL"<<endl;

		for (size_t i=0; i<m; ++i)
			for (size_t j =0; j<n; ++j)
                cerr<<"C["<<i<<","<<j<<"]="<<(*(C+i*n+j))<<std::endl;

		for (size_t i=0; i<m; ++i)
			for (size_t j =0; j<n; ++j)
                cerr<<"Cd["<<i<<","<<j<<"]="<<(*(Cd+i*n+j))<<std::endl;

		for (size_t i=0; i<m; ++i){
			for (size_t j =0; j<n; ++j)
				if (!F.areEqual( *(C+i*n+j), *(Cd+i*n+j) ) )
					 cerr<<"Erreur C["<<i<<","<<j<<"]="
					     <<(*(C+i*n+j))<<" Cd["<<i<<","<<j<<"]="
					     <<(*(Cd+i*n+j))<<endl;
		}
	}
	// else{
	// 	cerr<<"PASS"<<endl;
	// }
	delete[] Cd;
	//#endif
	/*
#if TIME
	double mflops = (2.0*(m*k-((!F.isZero(beta))?m:0))/1000000.0)*nbit*n/tim.realtime();
	cerr << pnbw << " Winograd's level over Z/"<<atoi(argv[1])<<"Z : t= "
	     << tim.realtime()/nbit
	     << " s, Mffops = "<<mflops
	     << endl;

	cerr<<"m,n,k,pnbw = "<<m<<", "<<n<<", "<<k<<", "<<alpha
	    <<", "<<beta<<", "<<pnbw<<endl
	    <<alpha
	    <<((ta==FFLAS::FflasNoTrans)?".Ax":".A^Tx")
	    <<((tb==FFLAS::FflasNoTrans)?"B + ":"B^T + ")
	    <<beta<<".C"<<endl;
	cout<<m<<" "<<n<<" "<<k<<" "<<pnbw<<" "<<alpha<<" "<<beta<<" "
	    <<mflops<<" "<<tim.realtime()/nbit<<endl;

	FFLAS::MMHelper<Field,FFLAS::MMHelperAlgo::Winograd> WH(F,nbw, FFLAS::ParSeqHelper::Sequential());

	OMPTimer tims,ts; ts.clear();tims.clear();
	for(int i = 0;i<nbit;++i){
        C = new Field::Element[m*n];
		ts.clear();
		ts.start();        
		FFLAS::fgemm (F, ta, tb,m,n,k,alpha, A,lda, B,ldb,
			      beta,C,n,WH);
		ts.stop();
		//if (i)
			tims+=ts;
//         if (i<nbit) delete[] C;
	}


	mflops = (2.0*(m*k-((!F.isZero(beta))?m:0))/1000000.0)*nbit*n/tims.realtime();
	cerr << nbw << " Winograd's level over Z/"<<atoi(argv[1])<<"Z : t= "
	     << tims.realtime()/nbit
	     << " s, Mffops = "<<mflops
	     << endl;

	cerr<<"m,n,k,nbw = "<<m<<", "<<n<<", "<<k<<", "<<alpha
	    <<", "<<beta<<", "<<nbw<<endl
	    <<alpha
	    <<((ta==FFLAS::FflasNoTrans)?".Ax":".A^Tx")
	    <<((tb==FFLAS::FflasNoTrans)?"B + ":"B^T + ")
	    <<beta<<".C"<<endl;
	cout<<m<<" "<<n<<" "<<k<<" "<<nbw<<" "<<alpha<<" "<<beta<<" "
	    <<mflops<<" "<<tims.realtime()/nbit<<endl;

#endif

	
    std::cerr << "Speed-up: " << tims.realtime()/tim.realtime() << std::endl;
	*/
	delete[] C;
	delete[] A;
	delete[] B;
	//#else
	//	std::cerr << "no openmp available" << std::endl;

	//#endif // __FFLASFFPACK_USE_OPENMP
	return wrong;
}
