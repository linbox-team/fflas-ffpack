/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
//--------------------------------------------------------------------------
//                        Test for fgesv : 1 computation
//
//--------------------------------------------------------------------------
// Clement Pernet
//-------------------------------------------------------------------------

#define DEBUG 1
#define TIME 1

#include <iomanip>
#include <iostream>
using namespace std;

#include "fflas-ffpack/field/modular-balanced.h"
#include "timer.h"
#include "Matio.h"
#include "fflas-ffpack/ffpack/ffpack.h"



typedef Modular<double> Field;

int main(int argc, char** argv){

	int n,m,mb,nb;
	cerr<<setprecision(10);
	Field::Element zero, one;

	if (argc != 6)	{
		cerr<<"Usage : test-fgesv <p> <A> <b> <iter> <left/right>"
		    <<endl;
		exit(-1);
	}
	int nbit=atoi(argv[4]); // number of times the product is performed
	Field F(atoi(argv[1]));
	F.init(zero,0.0);
	F.init(one,1.0);
	Field::Element * A, *B, *B2, *X=NULL;
	A = read_field(F,argv[2],&m,&n);
	B = read_field(F,argv[3],&mb,&nb);

	FFLAS::FFLAS_SIDE side = (atoi(argv[5])) ? FFLAS::FflasRight :  FFLAS::FflasLeft;

	size_t ldx=0;
	size_t rhs = (side == FFLAS::FflasLeft) ? nb : mb;
	if (m != n) {
		if (side == FFLAS::FflasLeft){
			X = new Field::Element[n*nb];
			ldx = nb;
		}
		else {
			X = new Field::Element[mb*m];
			ldx = m;
		}
	}

	if ( ((side == FFLAS::FflasRight) && (n != nb))
	     || ((side == FFLAS::FflasLeft)&&(m != mb)) ) {
		cerr<<"Error in the dimensions of the input matrices"<<endl;
		exit(-1);
	}
	int info=0;
	Timer t; t.clear();
	double time=0.0;
	//write_field(F, cerr<<"A="<<endl, A, k,k,k);
	size_t R=0;
	for (int i = 0;i<nbit;++i){
		t.clear();
		t.start();
		if (m == n)
			R = FFPACK::fgesv (F, side, mb, nb, A, n, B, nb, &info);
		else
			R = FFPACK::fgesv (F, side, m, n, rhs, A, n, X, ldx, B, nb, &info);
		if (info > 0){
			std::cerr<<"System is inconsistent"<<std::endl;
			exit(-1);
		}

		t.stop();
		time+=t.usertime();
		if (i+1<nbit){
			delete[]A;
			A = read_field(F,argv[2],&m,&n);
			delete[] B;
			B = read_field(F,argv[3],&mb,&nb);
		}
	}

#if DEBUG
	delete[] A;

	if (info > 0){
		std::cerr<<"System inconsistent"<<std::endl;
		exit (-1);
	}

	A = read_field(F,argv[2],&m,&n);

	B2 = new Field::Element[mb*nb];


	if (m==n)
		if (side == FFLAS::FflasLeft)
			FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, nb, n,
				      one, A, n, B, nb, zero, B2, nb);
		else
			FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mb, n, m,
				      one, B, nb, A, n, zero, B2, nb);
	else
		if (side == FFLAS::FflasLeft)
			FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, nb, n,
				      one, A, n, X, ldx, zero, B2, nb);
		else
			FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mb, n, m,
				      one, X, ldx, A, n, zero, B2, nb);
	delete[] B;
	delete[] X;

	B = read_field(F,argv[3],&mb,&nb);

	bool wrong = false;
	for (int i=0;i<mb;++i)
		for (int j=0;j<nb;++j)
			if ( !F.areEqual(*(B2+i*nb+j), *(B+i*nb+j))){
				cerr<<"B2 ["<<i<<", "<<j<<"] = "<<(*(B2+i*nb+j))
				    <<" ; B ["<<i<<", "<<j<<"] = "<<(*(B+i*nb+j))
				    <<endl;
				wrong = true;
			}

	if (wrong) {
		cerr<<"FAIL"<<endl;
		//write_field (F,cerr<<"B2="<<endl,B2,m,n,n);
		//write_field (F,cerr<<"B="<<endl,B,m,n,n);
	}else{

		cerr<<"PASS"<<endl;
	}


#endif

	delete[] A;
	delete[] B;
	delete[] B2;
#if TIME
	double mflops;
	if (side == FFLAS::FflasLeft)
		mflops = (n*m*m-m*m*m/3+2*R*R*n)/1000000.0*nbit/time;
	else
		mflops = (n*m*m-m*m*m/3+2*R*R*m)/1000000.0*nbit/time;
	cerr<<"m,n,mb,nb = "<<m<<" "<<n<<" "<<mb<<" "<<nb<<". fgesv "
	    <<((side == FFLAS::FflasLeft)?" Left ":" Right ")
	    <<"over Z/"<<atoi(argv[1])<<"Z :"
	    <<endl
	    <<"t= "
	    << time/nbit
	    << " s, Mffops = "<<mflops
	    << endl;

	cout<<m<<" "<<n<<" "<<mflops<<" "<<time/nbit<<endl;
#endif
}
