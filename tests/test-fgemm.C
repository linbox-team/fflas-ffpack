/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Clément Pernet
 *            Brice Boyer <bbboyer@ncsu.edu>
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
// Clement Pernet
//-------------------------------------------------------------------------

//#define DEBUG 1
#define NEWWINO
#define TIME 1

#include <iomanip>
#include <iostream>
using namespace std;

//#include "fflas-ffpack/field/modular-positive.h"
//#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-int32.h"
#include "fflas-ffpack/utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"

using namespace FFPACK;


//typedef Modular<double> Field;
//typedef Modular<float> Field;
typedef ModularBalanced<double> Field;
//typedef ModularBalanced<float> Field;
//typedef Modular<int> Field;

int main(int argc, char** argv){

	int m,n,k;
	int nbw=atoi(argv[4]); // number of winograd levels
	int nbit=atoi(argv[5]); // number of times the product is performed
	cerr<<setprecision(10);
	Field::Element alpha,beta;


	if (argc != 11)	{
		cerr<<"Usage : test-fgemm <p> <A> <B> <w> <i>"
		    <<" <alpha> <beta> <C> <ta> <tb>"<<endl
		    <<"         to do i computations of C <- alpha AB + beta C"
		    <<" using w recursive levels of Winograd's algorithm"
		    <<endl
		    <<"         if ta=1 (resp tb=1), A (resp B) is transposed."
		    <<endl;
		exit(-1);
	}
	Field F(atoi(argv[1]));

	F.init( alpha, Field::Element(atoi(argv[6])));
	F.init( beta, Field::Element(atoi(argv[7])));

	Field::Element * A;
	Field::Element * B;
	size_t lda;
	size_t ldb;

	enum FFLAS::FFLAS_TRANSPOSE ta = FFLAS::FflasNoTrans;
	enum FFLAS::FFLAS_TRANSPOSE tb = FFLAS::FflasNoTrans;
	if (atoi(argv[9])){
		ta = FFLAS::FflasTrans;
		A = read_field(F,argv[2],&k,&m);
			lda = m;
	}
	else{
		A = read_field(F,argv[2],&m,&k);
		lda = k;
	}
	if (atoi(argv[10])){
		tb = FFLAS::FflasTrans;
		B = read_field(F,argv[3],&n,&k);
		ldb = k;
	}
	else{
		B = read_field(F,argv[3],&k,&n);
		ldb = n;
	}

	Field::Element * C=NULL;

// 	write_field (F, cerr<<"A = "<<endl, A, m, k, lda);
// 	write_field (F, cerr<<"B = "<<endl, B, k, n, ldb);
	Timer tim,t; t.clear();tim.clear();
	for(int i = 0;i<nbit;++i){
		if (!F.isZero(beta)){
			C = read_field(F,argv[8],&m,&n);
		}else
			C = new Field::Element[m*n];
		t.clear();
		t.start();
		FFLAS::fgemm (F, ta, tb,m,n,k,alpha, A,lda, B,ldb,
			      beta,C,n,nbw);
		t.stop();
		tim+=t;
		// if (i<nbit-1) delete[] C; /*  deleted after... */
	}

#if DEBUG
	bool wrong = false;
	Field::Element zero;
	F.init(zero, 0.0);
	Field::Element * Cd;
	if (!F.isZero(beta))
 		Cd = read_field(F,argv[8],&m,&n);
	else{
		Cd  = new Field::Element[m*n];
		for (int i=0; i<m*n; ++i)
			F.assign (*(Cd+i), zero);
	}
	Field::Element aij, bij,  tmp;
	// Field::Element beta_alpha;
	//F.div (beta_alpha, beta, alpha);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j){
			F.mulin(*(Cd+i*n+j),beta);
			F.assign (tmp, zero);
			for ( int l = 0; l < k ; ++l ){
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
		for (int i=0; i<m; ++i){
			for (int j =0; j<n; ++j)
				if (!F.areEqual( *(C+i*n+j), *(Cd+i*n+j) ) )
					 cerr<<"Erreur C["<<i<<","<<j<<"]="
					     <<(*(C+i*n+j))<<" C[d"<<i<<","<<j<<"]="
					     <<(*(Cd+i*n+j))<<endl;
		}
	}
	else{
		cerr<<"PASS"<<endl;
	}
	delete[] Cd;
#endif

	delete[] C;
	delete[] A;
	delete[] B;
#if TIME
	double mflops = (2.0*(m*k-((!F.isZero(beta))?m:0))/1000000.0)*nbit*n/tim.usertime();
	cerr << nbw << " Winograd's level over Z/"<<atoi(argv[1])<<"Z : t= "
	     << tim.usertime()/nbit
	     << " s, Mffops = "<<mflops
	     << endl;

	cerr<<"m,n,k,nbw = "<<m<<", "<<n<<", "<<k<<", "<<alpha
	    <<", "<<beta<<", "<<nbw<<endl
	    <<alpha
	    <<((ta==FFLAS::FflasNoTrans)?".Ax":".A^Tx")
	    <<((tb==FFLAS::FflasNoTrans)?"B + ":"B^T + ")
	    <<beta<<".C"<<endl;
	cout<<m<<" "<<n<<" "<<k<<" "<<nbw<<" "<<alpha<<" "<<beta<<" "
	    <<mflops<<" "<<tim.usertime()/nbit<<endl;
#endif
}


