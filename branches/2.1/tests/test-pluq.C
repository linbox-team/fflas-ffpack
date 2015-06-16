/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * Copyright (C) 2015 the FFLAS-FFPACK group
 * Written by
 *
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
 *
 */

//--------------------------------------------------------------------------
//          Test for the lqup factorisation
//--------------------------------------------------------------------------
// usage: test-lqup p A n, for n lqup factorization
// of A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//#define DEBUG 0
#define __FFLAS__TRSM_READONLY
// Debug option  0: no debug
//               1: check A = LQUP
//-------------------------------------------------------------------------


#define __FFPACK_LUDIVINE_CUTOFF 60
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/timer.h"
#include "givaro/modular-integer.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "test-utils.h"

using namespace std;
using namespace FFPACK;


typedef Givaro::Modular<Givaro::Integer> Field;


int main(int argc, char** argv){
	    //cerr<<setprecision(20);
	int m,n;
	size_t R;

	if (argc!=4){
		cerr<<"usage : test-plup <p> <A>  <i>"<<endl
		    <<"        to do i PLUQ factorisation of A"
		    <<endl;
		exit(-1);
	}
	Field F(atof(argv[1]));
	Field::Element * A;

	A = read_field(F,argv[2],&m,&n);

	size_t maxP, maxQ;

	//	size_t cutoff = atoi(argv[3]);
	size_t nbf = atoi(argv[3]);

 FFLAS::Timer tim,timc, timlud,timludc;
	timc.clear();
	timludc.clear();

	enum FFLAS::FFLAS_DIAG diag = FFLAS::FflasNonUnit;
	enum FFLAS::FFLAS_TRANSPOSE trans = FFLAS::FflasNoTrans;
	if (trans == FFLAS::FflasNoTrans){
		maxP = m;
		maxQ = n;
	} else{
		maxP = n;
		maxQ = m;
	}
	size_t *P = FFLAS::fflas_new<size_t>(maxP);
	size_t *Q = FFLAS::fflas_new<size_t>(maxQ);

	//write_field (F,cerr<<"A = "<<endl, A, m,n,n);
	size_t * RRP, *CRP;
	for ( size_t i=0;i<nbf;i++){
		if (i) {
			FFLAS::fflas_delete( A);
			FFLAS::fflas_delete( RRP);
			FFLAS::fflas_delete( CRP);
			A = read_field(F,argv[2],&m,&n);
		}

		for (size_t j=0;j<maxP;j++)
			P[j]=0;
		for (size_t j=0;j<maxQ;j++)
			Q[j]=0;
		tim.clear();
		tim.start();

		R = FFPACK::PLUQ_basecaseCrout (F, diag, m, n, A, n, P, Q);
		tim.stop();
		timc+=tim;
		FFLAS::fflas_delete( A);
		A = read_field(F,argv[2],&m,&n);
		timlud.clear();
		timlud.start();
		R = FFPACK::LUdivine (F, diag, FFLAS::FflasNoTrans, m, n, A, n, P, Q);
		timlud.stop();
		timludc+=timlud;
//		std::cerr<<"Fini LUdivine"<<std::endl;
		RRP = FFLAS::fflas_new<size_t>(R);
		CRP = FFLAS::fflas_new<size_t>(R);
//		RankProfilesFromPLUQ(RRP, CRP, P, Q, m, n, R);
	}
	    // cerr<<"Row Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<RRP[i]<<" ";
	// cerr<<endl;
	// cerr<<"Column Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<CRP[i]<<" ";
	// cerr<<endl;
	// std::sort(CRP,(CRP+R));
	// std::sort(RRP,(RRP+R));
	// cerr<<"Sorted Row Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<RRP[i]<<" ";
	// cerr<<endl;
	// cerr<<"Sorted Column Rank Profile = ";
	// for (size_t i=0;i<R;++i)
	// 	cerr<<CRP[i]<<" ";
	// cerr<<endl;

	if (nbf){
		FFLAS::fflas_delete( RRP);
		FFLAS::fflas_delete( CRP);
	}
//	write_field (F,cerr<<"Result = "<<endl, A, m,n,n);

// 	cerr<<"P = [";
// 	for (size_t i=0; i<maxP; ++i)
// 		cerr<<P[i]<<" ";
// 	cerr<<"]"<<endl;
	// cerr<<"Q = [";
	// for (size_t i=0; i<maxQ; ++i)
	// 	cerr<<Q[i]<<" ";
	// cerr<<"]"<<endl;
#if DEBUG
	Field::Element * X = FFLAS::fflas_new<Field::Element>(m*n);
	Field::Element * L, *U;
	L = FFLAS::fflas_new<Field::Element>(m*R);
	U = FFLAS::fflas_new<Field::Element>(R*n);

	Field::Element zero,one;
	F.init(zero,0.0);
	F.init(one,1.0);
	for (size_t  i=0; i<R; ++i){
		for (size_t j=0; j<i; ++j)
			F.assign ( *(U + i*n + j), zero);
		for (int j=i; j<n; ++j)
			F.assign (*(U + i*n + j), *(A+ i*n+j));
	}
	for ( size_t j=0; j<R; ++j ){
		for (size_t i=0; i<=j; ++i )
			F.assign( *(L+i*R+j), zero);
		F.assign(*(L+j*R+j), one);
		for (size_t i=j+1; i<(size_t)m; i++)
			F.assign( *(L + i*R+j), *(A+i*n+j));
	}

	    //write_field(F,cerr<<"L = "<<endl,L,m,R,R);
	    //write_field(F,cerr<<"U = "<<endl,U,R,n,n);
	// cerr<<endl;
	FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,m, L, R, P);

	    //write_field(F,cerr<<"L = "<<endl,L,m,m,m);
	    //write_field(F,cerr<<"U = "<<endl,U,m,n,n);
// 		write_field(F,cerr<<"L = "<<endl,L,m,m,m);
// 		write_field(F,cerr<<"U = "<<endl,U,m,n,n);

	FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,R,
		      1.0, L,R, U,n, 0.0, X,n);
	    //FFLAS::fflas_delete( A);

//////
	    //write_field(F,cerr<<"L = "<<endl,L,m ,n,n);
	    //write_field(F,cerr<<"U = "<<endl,U,n,n,n);

	 // cerr<<"P = ";
	 // for (int i=0; i<m; ++i)
	 // 	cerr<<P[i]<<" ";
	 // cerr<<endl;
	 // cerr<<"Q = ";
	 // for (int i=0; i<n; ++i)
	 // 	cerr<<Q[i]<<" ";
	 // cerr<<endl;

	Field::Element * B =  read_field(F,argv[2],&m,&n);

	bool fail = false;
	for (size_t i=0; i<(size_t)m; ++i)
		for (size_t j=0; j<(size_t)n; ++j)
			if (!F.areEqual (*(B+i*n+j), *(X+i*n+j))){
				std::cerr << " B["<<i<<","<<j<<"] = " << (*(B+i*n+j))
					  << " X["<<i<<","<<j<<"] = " << (*(X+i*n+j))
					  << endl;
				fail=true;
			}
	// write_field(F,cerr<<"X = "<<endl,X,m,n,n);
	// write_field(F,cerr<<"B = "<<endl,B,m,n,n);
	FFLAS::fflas_delete( B);
	if (fail)
		cerr<<"FAIL"<<endl;


	else
		cerr<<"PASS"<<endl;
	FFLAS::fflas_delete( U);
	FFLAS::fflas_delete( L);
	FFLAS::fflas_delete( X);
#endif
	FFLAS::fflas_delete( A);
	FFLAS::fflas_delete( P);
	FFLAS::fflas_delete( Q);

	double t = timc.usertime();
	double tlud = timludc.usertime();
	const int sm = std::min(m,n);
	const int sn = std::max(m,n);

	double numops = sm*sm/1000.0*(sn-sm/3.0);

	cout<<m<<"x"<< n
	    << " Trans = "<<trans
	    << " Diag = "<<diag
	    << " : rank = " << R << "  ["
	    << ((double)nbf/1000.0*(double)numops / t)
	    << " MFops "
	    << " in "
	    << t/nbf<<"s"
	    <<" LUdivine : "<<((double)nbf/1000.0*(double)numops / tlud)
	    << " MFops "
	    << " in "
	    <<tlud/nbf
	    <<"s]"<< endl;
	cerr<<m
	    <<" "<<((double)nbf/1000.0*(double)numops / t)
	    <<" "<<t/nbf
	    <<" "<<((double)nbf/1000.0*(double)numops / tlud)
	    <<" "<<tlud/nbf
	    <<" "<<R
	    <<endl;

	return 0;
}
