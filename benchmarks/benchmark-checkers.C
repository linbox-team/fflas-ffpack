#define ENABLE_ALL_CHECKINGS 1 // DO NOT CHANGE
#define NR_TESTS 10
#define MAX_SIZE_MATRICES 3000
#define NR_THREADS 4

#include "fflas-ffpack/config-blas.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/checkers/checkers.h"
#include <fstream>

using namespace std;

int main(int argc, char** argv) {

	if (argc != 2) {
		cout << "./benchmark-checkers <path to stat file>" << endl;
		return -1;
	}

        std::ofstream stats_f(argv[1]);

	srand (time(NULL));
	typedef Givaro::Modular<double> Field;
	typedef std::vector<Field::Element> Polynomial;

	Givaro::Integer q = 131071;
	Field F(q);
	Field::RandIter Rand(F);
	Field::NonZeroRandIter NZRand(Rand);

	size_t pass;
	FFLAS::Timer chrono;
	double time1, time2;

	Field::Element_ptr A = FFLAS::fflas_new(F,MAX_SIZE_MATRICES,MAX_SIZE_MATRICES);
	Field::Element_ptr B = FFLAS::fflas_new(F,MAX_SIZE_MATRICES,MAX_SIZE_MATRICES);
	Field::Element_ptr C = FFLAS::fflas_new(F,MAX_SIZE_MATRICES,MAX_SIZE_MATRICES);
	typename Field::Element alpha,beta,tmp;
	F.init(alpha, rand()%1000+1);
	F.init(beta,  rand()%1000+1);
	size_t m,n,k,lda,ldb,ldc;
	FFLAS::FFLAS_TRANSPOSE ta,tb;

	stats_f << "     Matrix size\tSuccess rate\t\tTime comput.\t\tTime checker\n\n";

	// #####   FGEMM   #####
	stats_f << "FGEMM:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;
			n = rand() % 500 + i;
			k = rand() % 500 + i;
			ta = FFLAS::FflasNoTrans;//rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans,
			tb = FFLAS::FflasNoTrans;//rand()%2 ? FFLAS::FflasNoTrans : FFLAS::FflasTrans;
			lda = ta == FFLAS::FflasNoTrans ? k : m,
			ldb = tb == FFLAS::FflasNoTrans ? n : k,
			ldc = n;

                        PAR_BLOCK { FFLAS::pfrand(F,Rand, m,k,A,m/MAX_THREADS); }
                        PAR_BLOCK { FFLAS::pfrand(F,Rand, k,n,B,k/MAX_THREADS); }
                        PAR_BLOCK { FFLAS::pfrand(F,Rand, m,n,C,n/MAX_THREADS); }
//                         for( size_t i = 0; i < m*k; ++i ) Rand.random( *(A+i) );
// 			for( size_t i = 0; i < k*n; ++i ) Rand.random( *(B+i) );
// 			for( size_t i = 0; i < m*n; ++i ) Rand.random( *(C+i) );

			chrono.clear(); chrono.start();
                        FFLAS::Checker_fgemm<Field> checker1(Rand,m,n,k,beta,C,ldc);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			FFLAS::fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker1.check(ta,tb,alpha,A,lda,B,ldb,C) ? 1 : 0;
			chrono.stop(); time1 += chrono.usertime();
		}
		time1 /= (NR_THREADS*NR_TESTS);
		time2 /= (NR_THREADS*NR_TESTS);
		stats_f << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;



	// #####   FTRSM   #####
	stats_f << "FTRSM:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;
			n = rand() % 500 + i;
			FFLAS::FFLAS_SIDE side = rand()%2?FFLAS::FflasLeft:FFLAS::FflasRight;
			FFLAS::FFLAS_UPLO uplo = rand()%2?FFLAS::FflasLower:FFLAS::FflasUpper;
			FFLAS::FFLAS_TRANSPOSE trans = rand()%2?FFLAS::FflasNoTrans:FFLAS::FflasTrans;
			FFLAS::FFLAS_DIAG diag = rand()%2?FFLAS::FflasNonUnit:FFLAS::FflasUnit;
			k = (side==FFLAS::FflasLeft?m:n);

			for( size_t i = 0; i < m*n; ++i ) Rand.random( *(B+i) );
			for (size_t i=0;i<k;++i) {
				for (size_t j=0;j<i;++j)
					A[i*k+j]= (uplo == FFLAS::FflasLower)? Rand.random(tmp) : F.zero;
				A[i*k+i]= (diag == FFLAS::FflasNonUnit)? NZRand.random(tmp) : F.one;
				for (size_t j=i+1;j<k;++j)
					A[i*k+j]= (uplo == FFLAS::FflasUpper)? Rand.random(tmp) : F.zero;
			}

			chrono.clear(); chrono.start();
                        FFLAS::Checker_ftrsm<Field> checker2(Rand, m, n, alpha, B, n);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			FFLAS::ftrsm(F, side, uplo, trans, diag, m, n, alpha, A, k, B, n);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker2.check(side, uplo, trans, diag, m, n, A, k, B, n);
			chrono.stop(); time1 += chrono.usertime();
		}
		time1 /= (NR_THREADS*NR_TESTS);
		time2 /= (NR_THREADS*NR_TESTS);
		stats_f << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;



	// #####   INVERT   #####
	stats_f << "INVERT:\n";
	int nullity;
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;

                        PAR_BLOCK { FFLAS::pfrand(F,Rand, m,m,A,m/MAX_THREADS); }

			chrono.clear(); chrono.start();
                        FFPACK::Checker_invert<Field> checker3(Rand,m,A,m);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			FFPACK::Invert(F,m,A,m,nullity);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker3.check(A,nullity);
			chrono.stop(); time1 += chrono.usertime();
		}
		time1 /= (NR_THREADS*NR_TESTS);
		time2 /= (NR_THREADS*NR_TESTS);
		stats_f << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;
	



	// #####   PLUQ   #####
	stats_f << "PLUQ:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			m = rand() % 500 + i;
			n = rand() % 500 + i;

                        PAR_BLOCK { FFLAS::pfrand(F,Rand, m,n,A,m/MAX_THREADS); }

			size_t *P = FFLAS::fflas_new<size_t>(m);
			size_t *Q = FFLAS::fflas_new<size_t>(n);

			chrono.clear(); chrono.start();
                        FFPACK::Checker_PLUQ<Field> checker4 (Rand,m,n,A,n);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			k = FFPACK::PLUQ(F, FFLAS::FflasNonUnit, m, n, A, n, P, Q);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker4.check(A,n,k,P,Q);
			chrono.stop(); time1 += chrono.usertime();

			FFLAS::fflas_delete(P,Q);
		}
		time1 /= (NR_THREADS*NR_TESTS);
		time2 /= (NR_THREADS*NR_TESTS);
		stats_f << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}
	cout << endl;




	// #####   CharPoly   #####
	stats_f << "CharPoly:\n";
	for (size_t i=0; i<MAX_SIZE_MATRICES; i+=500) {
		pass = 0; time1 = 0.0; time2 = 0.0;
		for (size_t j=0; j<NR_TESTS; ++j) {
			n = rand() % 500 + i;

                        PAR_BLOCK { FFLAS::pfrand(F,Rand, n,n,A,n/MAX_THREADS); }

			Polynomial g(n);

			chrono.clear(); chrono.start();
                        FFPACK::Checker_charpoly<Field,Polynomial> checker5(Rand,n,A);
			chrono.stop(); time1 += chrono.usertime();

			chrono.clear(); chrono.start();
			FFPACK::CharPoly(F,g,n,A,n,FFPACK::FfpackLUK);
			chrono.stop(); time2 += chrono.usertime();

			chrono.clear(); chrono.start();
			pass += checker5.check(g);
			chrono.stop(); time1 += chrono.usertime();
		}
		time1 /= (NR_THREADS*NR_TESTS);
		time2 /= (NR_THREADS*NR_TESTS);
		stats_f << "     " << i << "-" << i+500 << "\t\t" << pass << "/" << NR_TESTS << "\t\t\t" << time2 
				<< "\t\t" << time1 << endl;
	}


        stats_f.close();

	return 0;
}
