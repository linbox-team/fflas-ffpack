#include "fflas-ffpack/utils/timer.h"
#include "Matio.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "test-utils.h"
#include "assert.h"
#include "fflas-ffpack/utils/args-parser.h"

// using namespace FFPACK;
#define NEWWINO
// #define NOTRANDOM

const int algos = 6 ;
using FFPACK::Modular;
using FFPACK::ModularBalanced;

const size_t selec[] = {
	0
	,1
	,2
	,3
	,4
	,5
};

const char * descr[] = {
	"322 low mem"
	, "322 first 1"
	, "322 4 tmp  "
	, "223 low mem"
	, "232 first 1"
	, "232 all tmp"
};

namespace FFLAS {


	void finit_fuzzy(Modular<double> F, size_t m, size_t n, double * C, size_t ldc)
	{

		double p, invp;
		p=(double)F.cardinality();
		invp=1./p;

		if (n == ldc)
			FFLAS::vectorised::modp<true,true>(C,C,m*n,p,invp,0,p-1);
		else
			for (size_t i = 0 ; i < m ; ++i)
				FFLAS::vectorised::modp<true,true>(C+i*ldc,C+i*ldc,n,p,invp,0,p-1);
	}

	void finit_fuzzy(ModularBalanced<double> F, size_t m, size_t n, double * C, size_t ldc)
	{

		double p, invp;
		p=(double)F.cardinality();
		invp=1./p;
		double pmax = (p-1)/2 ;
		double pmin = pmax-p+1;


		if (n == ldc)
			FFLAS::vectorised::modp<false,true>(C,C,m*n,p,invp,pmin, pmax);
		else
			for (size_t i = 0 ; i < m ; ++i)
				FFLAS::vectorised::modp<false,true>(C+i*ldc,C+i*ldc,n,p,invp,pmin,pmax);
	}


	// C = a*A + B
	void add(const size_t m, const size_t n,
		 double a,
		 const double *A, const size_t lda,
		 const double *B, const size_t ldb,
		 double *C, const size_t ldc)
	{
		const double *Ai = A,*Bi = B;
		double *Ci       = C;
		for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
			for (size_t j = 0 ; j < n ; ++j)
				Ci[j] = a * Ai[j] + Bi[j];
	}

	// C = C-(A+B)
	void subadd(const size_t m, const size_t n,
		    const double *A, const size_t lda,
		    const double *B, const size_t ldb,
		    double *C, const size_t ldc)
	{
		const double *Ai = A,*Bi = B;
		double *Ci       = C;
		for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
			for (size_t j = 0 ; j < n ; ++j) {
				Ci[j] = Ci[j] - Ai[j] - Bi[j] ;
			}

	}

	// C = -(A+B)
	void negadd(const size_t m, const size_t n,
		    const double *A, const size_t lda,
		    const double *B, const size_t ldb,
		    double *C, const size_t ldc)
	{
		const double *Ai = A,*Bi = B;
		double *Ci       = C;
		for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
			for (size_t j = 0 ; j < n ; ++j) {
				Ci[j] =  - Ai[j] - Bi[j] ;
			}

	}


	// C = C+A-B
	void addsub(const size_t m, const size_t n,
		    const double *A, const size_t lda,
		    const double *B, const size_t ldb,
		    double *C, const size_t ldc)
	{
		const double *Ai = A,*Bi = B;
		double *Ci       = C;
		for (;Ai < A+m*lda ; Ai+=lda,Bi+=ldb,Ci+=ldc)
			for (size_t j = 0 ; j < n ; ++j) {
				Ci[j] = Ci[j] + Ai[j] - Bi[j] ;
			}

	}


	// C = (C+B)/e
	template<class Field>
	void addscalinf(const Field & F, const size_t m, const size_t n,
			const double *B, const size_t ldb,
			double e,
			double *C, const size_t ldc)
	{
		const double * Bi = B;
		double * Ci = C;
		for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb)
			for (size_t j = 0 ; j < n ; ++j)
				Ci[j]= (Ci[j]+Bi[j])*e ;
		// F.init( Ci[j], (Ci[j]+Bi[j])/e );

	}

	// C = (C-B)/e
	template<class Field>
	void subscalinf(const Field & F, const size_t m, const size_t n,
			const double *B, const size_t ldb,
			double e,
			double *C, const size_t ldc)
	{
		const double * Bi = B;
		double * Ci = C;
		for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb)
			for (size_t j = 0 ; j < n ; ++j)
				Ci[j]= (Ci[j]-Bi[j])*e ;
		// F.init( Ci[j], (Ci[j]-Bi[j])/e );

	}

	// C = (D-B)/e
	template<class Field>
	void subscal(const Field & F, const size_t m, const size_t n,
		     const double *D, const size_t ldd,
		     const double *B, const size_t ldb,
		     double e,
		     double *C, const size_t ldc)
	{
		const double * Bi = B;
		const double * Di = D;
		double * Ci = C;
		for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb, Di += ldd)
			for (size_t j = 0 ; j < n ; ++j)
				Ci[j] = (Di[j]-Bi[j])*e ;

	}

	// C = (D+B)/e
	template<class Field>
	void addscal(const Field & F, const size_t m, const size_t n,
		     const double *D, const size_t ldd,
		     const double *B, const size_t ldb,
		     double e,
		     double *C, const size_t ldc)
	{
		const double * Bi = B;
		const double * Di = D;
		double * Ci = C;
		for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb, Di += ldd)
			for (size_t j = 0 ; j < n ; ++j)
				Ci[j] = (Di[j]+Bi[j])*e ;

	}

	// C = C + (D-B)/e
	template<class Field>
	void subscalacc(const Field & F, const size_t m, const size_t n,
			const double *D, const size_t ldd,
			const double *B, const size_t ldb,
			double e,
			double *C, const size_t ldc)
	{
		const double * Bi = B;
		const double * Di = D;
		double * Ci = C;
		for (;Bi < B+m*ldb ; Ci+=ldc, Bi += ldb, Di += ldd)
			for (size_t j = 0 ; j < n ; ++j)
				Ci[j] += (Di[j]-Bi[j])*e ;

	}

#ifndef TRE
	// #ifndef NDEBUG
	#define TRE 1
	// #else
// #define TRE (size_t)(__FFLASFFPACK_WINOTHRESHOLD)
	// #define TRE (size_t)(__FFLASFFPACK_WINOTHRESHOLD*0.9)
	// #endif
#endif
	template<class Field>
	double * gemm_fflas(const Field & F,
			    const size_t m, const size_t n, const size_t k,
			    const double *A, size_t lda,
			    const double *B, size_t ldb,
			    double * C, size_t ldc,
			    int rec =  0)
	{
		FFLAS::DoubleDomain R;
		FFLAS::fgemm(R,
			     FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
			     m,n,k,
			     1,
			     A,lda, B,ldb,
			     0,
			     C, ldc);


		return C;
	}


	namespace Protected { namespace Rec {

		// Field must be Modular<double>
		template<class Field>
		double *
		gemm_bini_322_0(const Field & F
				, const size_t m
				, const size_t n
				, const size_t k
				, const double *A , const size_t lda
				, const double *B , const size_t ldb
				, double *C , const size_t ldc
				, int rec
				, const double & epsilon
			       )
		{
			FFPACK::UnparametricField<double>   NoField;
			// const double p = (double)F.characteristic();
			size_t M = (n>m)?std::min(k,m):std::min(k,n);
			// std::cout << rec << ',' <<  M  << std::endl;
			// Field G(p*p);

			if ( M < TRE  || rec <= 0) {
				return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
			}

			assert(k/2*2==k); // k divisible par 2
			assert(n/2*2==n); // n divisible par 2
			assert(m/3*3==m); // m divisible par 3


			size_t n2 = n/2;
			size_t k2 = k/2;
			size_t m3 = m/3;

			// std::cout << "€ = " << epsilon << std::endl;

			// sub matrices in A
			const double * A11 = A;
			const double * A12 = A   +k2;
			const double * A21 = A   +lda*m3;
			const double * A22 = A21 +k2;
			const double * A31 = A21 +lda*m3;
			const double * A32 = A31 +k2;

			// sub matrices in C
			double * C11 = C;
			double * C12 = C   +n2;
			double * C21 = C   +ldc*m3;
			double * C22 = C21 +n2;
			double * C31 = C21 +ldc*m3;
			double * C32 = C31 +n2;

			// sub matrices in B
			const double * B11 = B;
			const double * B12 = B   +n2;
			const double * B21 = B   +ldb*k2;
			const double * B22 = B21 +n2;

			FFLAS::fzero(NoField,m,n,C,ldc);

			/*
			 * Algo :
			 * S1  := A11  +A22;
			 * S4  := e*A12+A22;
			 * S5  := A11  +e*A12;
			 * S6  := A21  +A32;
			 * S9  := A21  +e*A31;
			 * S10 := e*A31+A32;
			 *
			 * T1  := e*B11 +B22;
			 * T2  := B21   +B22;
			 * T4  := -e*B11+B21;
			 * T5  := e*B12 +B22;
			 * T6  := B11   +e*B22;
			 * T7  := B11   +B12;
			 * T9  := B12   -e*B22;
			 * T10 := B11   +e*B21;
			 *
			 * P1 := S1 *T1;
			 * P2 := A22*T2;
			 * P3 := A11*B22;
			 * P4 := S4 *T4;
			 * P5 := S5 *T5;
			 * P6 := S6 *T6;
			 * P7 := A21*T7;
			 * P8 := A32*B11;
			 * P9 := S9 *T9;
			 * P10:= S10*T10;
			 *
			 * C11 := (P1-P2-P3+P4)/e;
			 * C12 := (P3-P5)/(-e) ;
			 * C21 := P4+P6-P10 ;
			 * C22 := P1-P5+P9;
			 * C31 := (-P8+P10)/e;
			 * C32 := (P6-P7-P8+P9)/e;
			 *
			 */

			double * S1 = new double[m3*k2] ;
			// double * C11t = new double[n2*m3] ;
			// S1  := A11  +A22;
			FFLAS::fadd(NoField,m3,k2,A11,lda,A22,lda,S1,k2);
			// T1  := e*B11 +B22;
			double * T1 = new double[n2*k2] ; // ou aire
			add(k2,n2,epsilon,B11,ldb,B22,ldb,T1,n2);
			// P1 := S1 *T1; (dans C22)
			gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,C22,ldc,rec-1,epsilon);
			// S4  := e*A12+A22;
			double * eA12 = new double [m3*k2];
			FFLAS::fscal(NoField,m3,k2,epsilon,A12,lda,eA12,k2) ;
			FFLAS::fadd(NoField,m3,k2,eA12,k2,A22,lda,S1,k2);
			// T4  := -e*B11+B21;
			add(k2,n2,-epsilon,B11,ldb,B21,ldb,T1,n2);
			// P4 := S4 *T4; (dans C21)
			gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,C21,ldc,rec-1,epsilon);
			// C11 = P1+P4
			FFLAS::fadd(NoField,m3,n2,C21,ldc,C22,ldc,C11,ldc);
			// T2  := B21  +B22;
			FFLAS::fadd(NoField,k2,n2,B21,ldb,B22,ldb,T1,n2);
			// P2 := A22*T2;
			double * P1 = new double[n2*m3] ; // ou aire
			gemm_bini_322_0(F,m3,n2,k2,A22,lda,T1,n2,P1,n2,rec-1,epsilon);
			// P3 := A11*B22; (dans C12)
			gemm_bini_322_0(F,m3,n2,k2,A11,lda,B22,ldb,C12,ldc,rec-1,epsilon);
			// C11 -= (P2+P3)
			subadd(m3,n2,P1,n2,C12,ldc,C11,ldc);
			// S5  := A11  +e*A12;
			FFLAS::fadd(NoField,m3,k2,eA12,k2,A11,lda,S1,k2);
			// T5  := e*B12 +B22;
			add(k2,n2,epsilon,B12,ldb,B22,ldb,T1,n2);
			// P5 := S5 *T5;
			double * P2 = new double[n2*m3] ; // ou aire
			gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,P2,n2,rec-1,epsilon);
			// C12 -= P5
			subscalinf(NoField,m3,n2,P2,n2,-(double)1/epsilon,C12,ldc);
			// S6  := A21  +A32;
			FFLAS::fadd(NoField,m3,k2,A21,lda,A32,lda,S1,k2);
			// T6  := B11   +e*B22;
			add(k2,n2,epsilon,B22,ldb,B11,ldb,T1,n2);
			// P6 := S6 *T6; dans C32
			gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,C32,ldc,rec-1,epsilon);
			// C21+= P6
			FFLAS::faddin(NoField,m3,n2,C32,ldc,C21,ldc);
			// T7  := B11  +B12;
			FFLAS::fadd(NoField,k2,n2,B11,ldb,B12,ldb,T1,n2);
			// P7 := A21*T7; !signe
			gemm_bini_322_0(F,m3,n2,k2,A21,lda,T1,n2,P1,n2,rec-1,epsilon);
			// P8 := A32*B11; dans C31 !signe
			gemm_bini_322_0(F,m3,n2,k2,A32,lda,B11,ldb,C31,ldc,rec-1,epsilon);
			// C32 -= P8+P7
			subadd(m3,n2,P1,n2,C31,ldc,C32,ldc);
			// S9  := A21  +e*A31;
			double * eA31 = eA12 ;
			FFLAS::fscal(NoField,m3,k2,epsilon,A31,lda,eA31,k2);
			FFLAS::fadd(NoField,m3,k2,eA31,k2,A21,lda,S1,k2);
			// T9  := B12   -e*B22;
			add(k2,n2,-epsilon,B22,ldb,B12,ldb,T1,n2);
			// P9 := S9 *T9;
			gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,P1,n2,rec-1,epsilon);
			// C32= (C32+P9)/p
			addscalinf(NoField,m3,n2,P1,n2,(double)1/epsilon,C32,ldc);
			// C22+= P9-P5
			addsub(m3,n2,P1,n2,P2,n2,C22,ldc);
			delete[] P2;
			// S10 := e*A31+A32;
			FFLAS::fadd(NoField,m3,k2,eA31,k2,A32,lda,S1,k2);
			delete[] eA12 ;
			// T10 := B11   +e*B21;
			add(k2,n2,epsilon,B21,ldb,B11,ldb,T1,n2);
			// P10:= S10*T10;
			gemm_bini_322_0(F,m3,n2,k2,S1,k2,T1,n2,P1,n2,rec-1,epsilon);
			delete[] S1;
			delete[] T1;
			// C21-= P10
			FFLAS::fsubin(NoField,m3,n2,P1,n2,C21,ldc);
			// C31= (C31-P10)/(-epsilon)
			subscalinf(NoField,m3,n2,P1,n2,-(double)1/epsilon,C31,ldc);
			delete[] P1;
			// C11 := (P1+P-P3+P4)/e;
			FFLAS::fscalin(NoField,m3,n2,(double)1/epsilon,C11,ldc);

			return C;

		}

		// Field must be Modular<double>
		template<class Field>
		double *
		gemm_bini_322_mem(const Field & F
				  , const size_t m
				  , const size_t n
				  , const size_t k
				  , const double *A , const size_t lda
				  , const double *B , const size_t ldb
				  , double *C , const size_t ldc
				  , int rec
				  , const double & epsilon
				 )
		{
			FFPACK::UnparametricField<double>   NoField;
			// const double p = (double)F.characteristic();
			size_t M = (n>m)?std::min(k,m):std::min(k,n);
			// std::cout << rec << ',' <<  M  << std::endl;
			// Field G(p*p);

			if ( M < TRE  || rec <= 0) {
				// std::cout << "ffw" << std::endl;
				return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
				// return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
			}

			assert(k/2*2==k); // k divisible par 2
			assert(n/2*2==n); // n divisible par 2
			assert(m/3*3==m); // m divisible par 3

			// std::cout << "tested" << std::endl;

			size_t n2 = n/2;
			size_t k2 = k/2;
			size_t m3 = m/3;

			// std::cout << "€ = " << epsilon << std::endl;

			// sub matrices in A
			const double * A11 = A;
			const double * A12 = A   +k2;
			const double * A21 = A   +lda*m3;
			const double * A22 = A21 +k2;
			const double * A31 = A21 +lda*m3;
			const double * A32 = A31 +k2;

			// sub matrices in C
			double * C11 = C;
			double * C12 = C   +n2;
			double * C21 = C   +ldc*m3;
			double * C22 = C21 +n2;
			double * C31 = C21 +ldc*m3;
			double * C32 = C31 +n2;

			// sub matrices in B
			const double * B11 = B;
			const double * B12 = B   +n2;
			const double * B21 = B   +ldb*k2;
			const double * B22 = B21 +n2;

			FFLAS::fzero(F,m,n,C,ldc);

			/*
			 * Algo :
			 * S1  := A11  +A22;
			 * S4  := e*A12+A22;
			 * S5  := A11  +e*A12;
			 * S6  := A21  +A32;
			 * S9  := A21  +e*A31;
			 * S3 := e*A31+A32;
			 *
			 * T1  := e*B11 +B22;
			 * T2  := B21   +B22;
			 * T4  := -e*B11+B21;
			 * T5  := e*B12 +B22;
			 * T6  := B11   +e*B22;
			 * T7  := B11   +B12;
			 * T9  := B12   -e*B22;
			 * T3 := B11   +e*B21;
			 *
			 * P1 := S1 *T1;
			 * P2 := A22*T2;
			 * P10 := A11*B22;
			 * P4 := S4 *T4;
			 * P5 := S5 *T5;
			 * P6 := S6 *T6;
			 * P7 := A21*T7;
			 * P8 := A32*B11;
			 * P9 := S9 *T9;
			 * P3:= S3*T3;
			 *
			 * C11 := (P1-P2-P10+P4)/e;
			 * C12 := (P10-P5)/(-e) ;
			 * C21 := P4+P6-P3 ;
			 * C22 := P1-P5+P9;
			 * C31 := (-P8+P3)/e;
			 * C32 := (P6-P7-P8+P9)/e;
			 *
			 */


			// P10
			gemm_bini_322_mem(F,m3,n2,k2,A11,lda,B22,ldb,C11,ldc,rec-1,epsilon);
			// S5
			double * X = new double[m3*k2];
			add(m3,k2,epsilon,A12,lda,A11,lda,X,k2);
			// T5
			// double * Y = new double[std::max(k2,m3)*n2];
			double * Y = new double[k2*n2];
			add(k2,n2,epsilon,B12,ldb,B22,ldb,Y,n2);
			// P5
			gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C22,ldc,rec-1,epsilon);
			// C12
			subscal(NoField,m3,n2,C22,ldc,C11,ldc,(double)1/epsilon,C12,ldc);
			// T2
			FFLAS::fadd(NoField,k2,n2,B21,ldb,B22,ldb,Y,n2);
			// P2
			gemm_bini_322_mem(F,m3,n2,k2,A22,lda,Y,n2,C31,ldc,rec-1,epsilon);
			// C11
			FFLAS::faddin(NoField,m3,n2,C31,ldc,C11,ldc);
			// S1
			FFLAS::fadd(NoField,m3,k2,A11,lda,A22,lda,X,k2);
			// T1
			add(k2,n2,epsilon,B11,ldb,B22,ldb,Y,n2);
			// P1
			gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C21,ldc,rec-1,epsilon);
			// C22
			FFLAS::fsub(NoField,m3,n2,C21,ldc,C22,ldc,C22,ldc);
			// C11
			FFLAS::fsub(NoField,m3,n2,C21,ldc,C11,ldc,C11,ldc);
			// S4
			add(m3,k2,epsilon,A12,lda,A22,lda,X,k2);
			// T4
			add(k2,n2,-epsilon,B11,ldb,B21,ldb,Y,n2);
			// P4
			gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C21,ldc,rec-1,epsilon);
			// C11
			addscalinf(NoField,m3,n2,C21,ldc,(double)1/epsilon,C11,ldc);
			// S9
			add(m3,k2,epsilon,A31,lda,A21,lda,X,k2);
			// T9
			add(k2,n2,-epsilon,B22,ldb,B12,ldb,Y,n2);
			// P9
			gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C32,ldc,rec-1,epsilon);
			//  C22
			FFLAS::faddin(NoField,m3,n2,C32,ldc,C22,ldc);
			// S6
			FFLAS::fadd(NoField,m3,k2,A21,lda,A32,lda,X,k2);
			// T6
			add(k2,n2,epsilon,B22,ldb,B11,ldb,Y,n2);
			// P6
			gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C31,ldc,rec-1,epsilon);
			// C21
			FFLAS::faddin(NoField,m3,n2,C31,ldc,C21,ldc);
			// C32
			FFLAS::faddin(NoField,m3,n2,C31,ldc,C32,ldc);
			// T7
			FFLAS::fadd(NoField,k2,n2,B11,ldb,B12,ldb,Y,n2);
			// P7
			gemm_bini_322_mem(F,m3,n2,k2,A21,lda,Y,n2,C31,ldc,rec-1,epsilon);
			// if (epsilon > 1 && rec == 2) { FFLAS::finit(G,m3,n2,C31,ldc);}
			// C32
			FFLAS::fsubin(NoField,m3,n2,C31,ldc,C32,ldc);
			// S3
			add(m3,k2,epsilon,A31,lda,A32,lda,X,k2);
			// T3
			add(k2,n2,epsilon,B21,ldb,B11,ldb,Y,n2);
			// P3
			gemm_bini_322_mem(F,m3,n2,k2,X,k2,Y,n2,C31,ldc,rec-1,epsilon);
			delete[] X;
			delete[] Y ;
			// C21
			FFLAS::fsubin(NoField,m3,n2,C31,ldc,C21,ldc);
			// P8
			Y = new double[m3*n2];
			gemm_bini_322_mem(F,m3,n2,k2,A32,lda,B11,ldb,Y,n2,rec-1,epsilon);
			// C31
			subscalinf(NoField,m3,n2,Y,n2,(double)1/epsilon,C31,ldc);
			// FFLAS::fsubin(NoField,m3,n2,Y,n2,C31,ldc);
			// C32
			subscalinf(NoField,m3,n2,Y,n2,(double)1/epsilon,C32,ldc);
			// FFLAS::fsubin(NoField,m3,n2,Y,n2,C32,ldc);
			// FFLAS::fscalin(NoField,m3,n,(double)1/epsilon,C31,ldc);
			delete[] Y ;


			return C;

		}

		// Field must be Modular<double>
		template<class Field>
		double *
		gemm_bini_223_mem(const Field & F
				  , const size_t m
				  , const size_t n
				  , const size_t k
				  , const double *A , const size_t lda
				  , const double *B , const size_t ldb
				  , double *C , const size_t ldc
				  , int rec
				  , const double & epsilon
				 )
		{
			FFPACK::UnparametricField<double>   NoField;
			// const double p = (double)F.characteristic();
			size_t M = (n>m)?std::min(k,m):std::min(k,n);
			// std::cout << rec << ',' <<  M  << std::endl;
			// Field G(p*p);

			if ( M < TRE  || rec <= 0) {
				// std::cout << "ffw" << std::endl;
				return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
				// return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
			}

			assert(k/2*2==k); // k divisible par 2
			assert(n/3*3==n); // n divisible par 2
			assert(m/2*2==m); // m divisible par 3

			// std::cout << "tested" << std::endl;

			size_t m2 = m/2;
			size_t k2 = k/2;
			size_t n3 = n/3;

			// std::cout << "€ = " << epsilon << std::endl;

			// sub matrices in A
			const double * A11 = A;
			const double * A12 = A   +k2;
			const double * A21 = A   +lda*m2;
			const double * A22 = A21 +k2;

			// sub matrices in C
			double * C11 = C;
			double * C12 = C   +n3;
			double * C13 = C   +2*n3;
			double * C21 = C   +ldc*m2;
			double * C22 = C21 +n3;
			double * C23 = C21 +2*n3;



			// sub matrices in B
			const double * B11 = B;
			const double * B12 = B   +n3;
			const double * B13 = B   +2*n3;
			const double * B21 = B   +ldb*k2;
			const double * B22 = B21 +n3;
			const double * B23 = B21 +2*n3;


			FFLAS::fzero(F,m,n,C,ldc);

			/*
			 * Algo :
			 * S1  := B11  +B22;
			 * S4  := e*B21+B22;
			 * S5  := B11  +e*B21;
			 * S6  := B12  +B23;
			 * S9  := B12  +e*B13;
			 * S3 := e*B13+B23;
			 *
			 * T1  := e*A11 +A22;
			 * T2  := A12   +A22;
			 * T4  := -e*A11+A12;
			 * T5  := e*A21 +A22;
			 * T6  := A11   +e*A22;
			 * T7  := A11   +A21;
			 * T9  := A21   -e*A22;
			 * T3  := A11   +e*A12;
			 *
			 * P1 := S1 *T1;
			 * P2 := T2 * B22;
			 * P10 := A22 * B11;
			 * P4 := S4 *T4;
			 * P5 := S5 *T5;
			 * P6 := S6 *T6;
			 * P7 := T7*B12;
			 * P8 := A11*B23;
			 * P9 := S9 *T9;
			 * P3 := S3*T3;
			 *
			 * C11 := (P1-P2-P10+P4)/e;
			 * C21 := (P10-P5)/(-e) ;
			 * C12 := P4+P6-P3 ;
			 * C22 := P1-P5+P9;
			 * C13 := (-P8+P3)/e;
			 * C23 := (P6-P7-P8+P9)/e;
			 *
			 */


			// P10
			gemm_bini_223_mem(F,m2,n3,k2,A22,lda,B11,ldb,C11,ldc,rec-1,epsilon);
			// S5
			double * Y = new double[k2*n3];
			add(k2,n3,epsilon,B21,ldb,B11,ldb,Y,n3);
			// T5
			double * X = new double[m2*k2];
			add(m2,k2,epsilon,A21,lda,A22,lda,X,k2);
			// P5
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C22,ldc,rec-1,epsilon);
			// C12
			subscal(NoField,m2,n3,C22,ldc,C11,ldc,(double)1/epsilon,C21,ldc);
			// T2
			FFLAS::fadd(NoField,m2,k2,A12,lda,A22,lda,X,k2);
			// P2
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,B22,ldb,C13,ldc,rec-1,epsilon);
			// C11
			FFLAS::faddin(NoField,m2,n3,C13,ldc,C11,ldc);
			// S1
			FFLAS::fadd(NoField,k2,n3,B11,ldb,B22,ldb,Y,n3);
			// T1
			add(m2,k2,epsilon,A11,lda,A22,lda,X,k2);
			// P1
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C12,ldc,rec-1,epsilon);
			// C22
			FFLAS::fsub(NoField,m2,n3,C12,ldc,C22,ldc,C22,ldc);
			// C11
			FFLAS::fsub(NoField,m2,n3,C12,ldc,C11,ldc,C11,ldc);
			// S4
			add(k2,n3,epsilon,B21,ldb,B22,ldb,Y,n3);
			// T4
			add(m2,k2,-epsilon,A11,lda,A12,lda,X,k2);
			// P4
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C12,ldc,rec-1,epsilon);
			// C11
			addscalinf(NoField,m2,n3,C12,ldc,(double)1/epsilon,C11,ldc);
			// S9
			add(k2,n3,epsilon,B13,ldb,B12,ldb,Y,n3);
			// T9
			add(m2,k2,-epsilon,A22,lda,A21,lda,X,k2);
			// P9
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C23,ldc,rec-1,epsilon);
			//  C22
			FFLAS::faddin(NoField,m2,n3,C23,ldc,C22,ldc);
			// S6
			FFLAS::fadd(NoField,k2,n3,B12,ldb,B23,ldb,Y,n3);
			// T6
			add(m2,k2,epsilon,A22,lda,A11,lda,X,k2);
			// P6
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C13,ldc,rec-1,epsilon);
			// C21
			FFLAS::faddin(NoField,m2,n3,C13,ldc,C12,ldc);
			// C32
			FFLAS::faddin(NoField,m2,n3,C13,ldc,C23,ldc);
			// T7
			FFLAS::fadd(NoField,m2,k2,A11,lda,A21,lda,X,k2);
			// P7
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,B12,ldb,C13,ldc,rec-1,epsilon);
			// if (epsilon > 1 && rec == 2) { FFLAS::finit(G,m2,n3,C31,ldc);}
			// C32
			FFLAS::fsubin(NoField,m2,n3,C13,ldc,C23,ldc);
			// S3
			add(k2,n3,epsilon,B13,ldb,B23,ldb,Y,n3);
			// T3
			add(m2,k2,epsilon,A12,lda,A11,lda,X,k2);
			// P3
			gemm_bini_223_mem(F,m2,n3,k2,X,k2,Y,n3,C13,ldc,rec-1,epsilon);
			delete[] Y ;
			delete[] X ;
			// C21
			FFLAS::fsubin(NoField,m2,n3,C13,ldc,C12,ldc);
			// P8
			Y = new double[m2*n3];
			gemm_bini_223_mem(F,m2,n3,k2,A11,lda,B23,ldb,Y,n3,rec-1,epsilon);
			// C31
			subscalinf(NoField,m2,n3,Y,n3,(double)1/epsilon,C13,ldc);
			// C32
			subscalinf(NoField,m2,n3,Y,n3,(double)1/epsilon,C23,ldc);
			delete[] Y ;


			return C;

		}

		// Field must be Modular<double>
		template<class Field>
		double *
		gemm_bini_322_2(const Field & F
				, const size_t m
				, const size_t n
				, const size_t k
				, const double *A , const size_t lda
				, const double *B , const size_t ldb
				, double *C , const size_t ldc
				, int rec
				, const double & epsilon
			       )
		{
			FFPACK::UnparametricField<double>   NoField;
			// const double p = (double)F.characteristic();
			size_t M = (n>m)?std::min(k,m):std::min(k,n);
			// std::cout << rec << ',' <<  M  << std::endl;
			// Field G(p*p);

			if ( M < TRE  || rec <= 0) {
				// std::cout << "ffw" << std::endl;
				return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
				// return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
			}

			assert(k/2*2==k); // k divisible par 2
			assert(n/2*2==n); // n divisible par 2
			assert(m/3*3==m); // m divisible par 3

			// std::cout << "tested" << std::endl;

			size_t n2 = n/2;
			size_t k2 = k/2;
			size_t m3 = m/3;

			// std::cout << "€ = " << epsilon << std::endl;

			// sub matrices in A
			const double * A11 = A;
			const double * A12 = A   +k2;
			const double * A21 = A   +lda*m3;
			const double * A22 = A21 +k2;
			const double * A31 = A21 +lda*m3;
			const double * A32 = A31 +k2;

			// sub matrices in C
			double * C11 = C;
			double * C12 = C   +n2;
			double * C21 = C   +ldc*m3;
			double * C22 = C21 +n2;
			double * C31 = C21 +ldc*m3;
			double * C32 = C31 +n2;

			// sub matrices in B
			const double * B11 = B;
			const double * B12 = B   +n2;
			const double * B21 = B   +ldb*k2;
			const double * B22 = B21 +n2;

			FFLAS::fzero(F,m,n,C,ldc);

			/*
			 * Algo :
			 * S1 := A11  +A22;
			 * S4 := e*A12+A22;
			 * S5 := A11  +e*A12;
			 * S3 := e*A31+A32;
			 * S6 := A21  +A32;
			 * S9 := A21  +e*A31;
			 *
			 * T1 := e*B11 +B22;
			 * T2 := B21   +B22;
			 * T3 := B11   +e*B21;
			 * T4 := -e*B11+B21;
			 * T5 := e*B12 +B22;
			 * T6 := B11   +e*B22;
			 * T7 := B11   +B12;
			 * T9 := B12   -e*B22;
			 *
			 * P1 := S1 *T1;
			 * P2 := A22*T2;
			 * P10 := A11*B22;
			 * P4 := S4 *T4;
			 * P5 := S5 *T5;
			 * P6 := S6 *T6;
			 * P7 := A21*T7;
			 * P8 := A32*B11;
			 * P9 := S9 *T9;
			 * P3:= S3*T3;
			 *
			 * C11 := (P1-P2-P10+P4)/e;
			 * C12 := (P10-P5)/(-e) ;
			 * C21 := P4+P6-P3 ;
			 * C22 := P1-P5+P9;
			 * C31 := (-P8+P3)/e;
			 * C32 := (P6-P7-P8+P9)/e;
			 *
			 */

			double * U = new double[m3*n2];
			double * V = new double[m3*n2];
			double * X = new double[m3*std::max(k2,n2)];
			double * Y = new double[std::max(k2,m3)*n2];

			// S4
			add(m3,k2,epsilon,A12,lda,A22,lda,X,k2);
			// T4
			add(k2,n2,-epsilon,B11,ldb,B21,ldb,Y,n2);
			// P4
			gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,U,n2,rec-1,epsilon);
			// S9
			add(m3,k2,epsilon,A31,lda,A21,lda,X,k2);
			// T9
			add(k2,n2,-epsilon,B22,ldb,B12,ldb,Y,n2);
			// P9
			gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,V,n2,rec-1,epsilon);
			// S5
			add(m3,k2,epsilon,A12,lda,A11,lda,X,k2);
			// T5
			add(k2,n2,epsilon,B12,ldb,B22,ldb,Y,n2);
			// P5
			gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,C12,ldc,rec-1,epsilon);
			// S3
			add(m3,k2,epsilon,A31,lda,A32,lda,X,k2);
			// T3
			add(k2,n2,epsilon,B21,ldb,B11,ldb,Y,n2);
			// P3
			gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,C31,ldc,rec-1,epsilon);
			// C22 = P9-P5
			FFLAS::fsub(NoField,m3,n2,V,n2,C12,ldc,C22,ldc);
			// C21 = P4-P3
			FFLAS::fsub(NoField,m3,n2,U,n2,C31,ldc,C21,ldc);
			// T2
			FFLAS::fadd(NoField,k2,n2,B21,ldb,B22,ldb,Y,n2);
			// P2
			gemm_bini_322_2(F,m3,n2,k2,A22,lda,Y,n2,X,n2,rec-1,epsilon);
			// XXX approximate
			// C11 = (P4 - P2) / e
			subscal(NoField,m3,n2,U,n2,X,n2,1./epsilon,C11,ldc);
			// T7
			FFLAS::fadd(NoField,k2,n2,B11,ldb,B12,ldb,Y,n2);
			// P7
			gemm_bini_322_2(F,m3,n2,k2,A21,lda,Y,n2,X,n2,rec-1,epsilon);
			// XXX approximate
			// C32 = (P9-P7) / e
			subscal(NoField,m3,n2,V,n2,X,n2,1./epsilon,C32,ldc);
			// S1
			FFLAS::fadd(NoField,m3,k2,A11,lda,A22,lda,X,k2);
			// T1
			add(k2,n2,epsilon,B11,ldb,B22,ldb,Y,n2);
			// P1
			gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,U,n2,rec-1,epsilon);
			// C22 += P1
			FFLAS::faddin(NoField,m3,n2,U,n2,C22,ldc);
			// P10
			gemm_bini_322_2(F,m3,n2,k2,A11,lda,B22,ldb,V,n2,rec-1,epsilon);
			// C12 = (P5-P10)/e
			subscalinf(NoField,m3,n2,V,n2,1./epsilon,C12,ldc);
			// XXX approximate
			// C11 = C11 + (P1-P10)/e
			subscalacc(NoField,m3,n2,U,n2,V,n2,1./epsilon,C11,ldc);
			// S6
			FFLAS::fadd(NoField,m3,k2,A21,lda,A32,lda,X,k2);
			// T6
			add(k2,n2,epsilon,B22,ldb,B11,ldb,Y,n2);
			// P6
			gemm_bini_322_2(F,m3,n2,k2,X,k2,Y,n2,U,n2,rec-1,epsilon);
			// C21 += P6
			FFLAS::faddin(NoField,m3,n2,U,n2,C21,ldc);
			// P8
			gemm_bini_322_2(F,m3,n2,k2,A32,lda,B11,ldb,V,n2,rec-1,epsilon);
			// C31 = (P3-P8)/2
			subscalinf(NoField,m3,n2,V,n2,1./epsilon,C31,ldc);
			// XXX approximate
			// C32 = C32 + (P6-P8)/e
			subscalacc(NoField,m3,n2,U,n2,V,n2,1./epsilon,C32,ldc);


			delete[] X;
			delete[] Y ;
			delete[] U;
			delete[] V;


			return C;

		}


		// Field must be Modular<double>
		template<class Field>
		double *
		gemm_bini_232_2(const Field & F
				, const size_t m
				, const size_t n
				, const size_t k
				, const double *A , const size_t lda
				, const double *B , const size_t ldb
				, double *C , const size_t ldc
				, int rec
				, const double & epsilon
			       )
		{
			FFPACK::UnparametricField<double>   NoField;
			// const double p = (double)F.characteristic();
			size_t M = (n>m)?std::min(k,m):std::min(k,n);
			// Field G(p*p);

			if ( M < TRE  || rec <= 0) {
				// std::cout << "ffw" << std::endl;
				return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
				// return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
			}

			assert(k/3*3==k); // k divisible par 3
			assert(n/2*2==n); // n divisible par 2
			assert(m/2*2==m); // m divisible par 2

			// std::cout << "tested" << std::endl;

			size_t n2 = n/2;
			size_t k3 = k/3;
			size_t m2 = m/2;

			// std::cout << "€ = " << epsilon << std::endl;

			// sub matrices in B
			const double * B11 = B;
			const double * B12 = B   +n2;
			const double * B21 = B   +ldb*k3;
			const double * B22 = B21 +n2;
			const double * B31 = B21 +ldb*k3;
			const double * B32 = B31 +n2;

			// sub matrices in C
			double * C11 = C;
			double * C12 = C   +n2;
			double * C21 = C   +ldc*m2;
			double * C22 = C21 +n2;

			// sub matrices in A

			const double * A11 = A;
			const double * A12 = A   +k3;
			const double * A13 = A   +2*k3;
			const double * A21 = A   +lda*m2;
			const double * A22 = A21 +k3;
			const double * A23 = A21 +2*k3;


			FFLAS::fzero(F,m,n,C,ldc);

			/*
			 * Algo :
			 *
			 * S1  := A11  +A22*e;
			 * S3  := -(A11+A21);
			 * S4  := A11+A12*e;
			 * S5  := A21 - A22*e;
			 * S6  := A12*e  + A23;
			 * S8 := -(A13+A23):
			 * S9  := A22*e  + A23;
			 * S10 := -A12*e+A13;
			 *
			 * T1  := B11 +B22;
			 * T4  := e*B12+B22;
			 * T5  := B11 +e*B12;
			 * T6  := B21   +B32;
			 * T9  := B21   + e*B31;
			 * T10 := e*B31   +B32;
			 *
			 * P1 := Bini232(S1,T1 ,e);
			 * P2 := Bini232(A11,B22 ,e);
			 * P3 := Bini232(S3,B11,e);
			 * P4 := Bini232(S4,T4 ,e);
			 * P5 := Bini232(S5,T5 ,e);
			 * P6 := Bini232(S6,T6 ,e);
			 * P7 := Bini232(A23,B21 ,e);
			 * P8 := Bini232(S8,B32,e);
			 * P9 := Bini232(S9,T9 ,e);
			 * P10:= Bini232(S10,T10,e);
			 *
			 *
			 * C11 := evalm(P1-P4+(P6-P7+P8+P10)/e);
			 * C12 := evalm((-P2+P4)/e+P10) ;
			 * C21 := evalm(P5+(-P7+P9)/e) ;
			 * C22 := evalm((P1-P2+P3+P5)/e+P6-P9);
			 *
			 */

			double * U = new double[m2*n2];
			double * V = new double[m2*n2];
			double * X = new double[m2*k3];
			double * Y = new double[k3*n2];

			// S1
			add(m2,k3,epsilon,A22,lda,A11,lda,X,k3);
			// T1
			FFLAS::fadd(NoField,k3,n2,B11,ldb,B22,ldb,Y,n2);
			// P1 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// S3
			negadd(m2,k3,A11,lda,A21,lda,X,k3);
			// P3 (in V)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,B11,ldb,V,n2,rec-1,epsilon);
			// C22 = (P1+P3)/e
			// FFLAS::fadd(NoField,m2,n2,U,n2,V,n2,C22,ldc); // XXX acc
			addscal(NoField,m2,n2,U,n2,V,n2,(double)1/epsilon,C22,ldc);
			// S6
			add(m2,k3,epsilon,A12,lda,A23,lda,X,k3);
			// T6
			FFLAS::fadd(NoField,k3,n2,B21,ldb,B32,ldb,Y,n2);
			// P6 (in V)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
			// C22 += P6
			FFLAS::faddin(NoField,m2,n2,V,n2,C22,ldc);
			// S8
			negadd(m2,k3,A13,lda,A23,lda,X,k3);
			// P8 (in C11)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,B32,ldb,C11,ldc,rec-1,epsilon);
			// C11 = (P8+P6)/e
			addscalinf(NoField,m2,n2,V,n2,(double)1/epsilon,C11,ldc);
			// C11 += P1
			FFLAS::faddin(NoField,m2,n2,U,n2,C11,ldc);
			// S4
			add(m2,k3,epsilon,A12,lda,A11,lda,X,k3);
			// T4
			add(k3,n2,epsilon,B12,ldb,B22,ldb,Y,n2);
			// P4 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// C11 -= P4
			FFLAS::fsubin(NoField,m2,n2,U,n2,C11,ldc);
			// P2 (in C12)
			gemm_bini_232_2(F,m2,n2,k3,A11,lda,B22,ldb,C12,ldc,rec-1,epsilon);
			// S5
			add(m2,k3,-epsilon,A22,lda,A21,lda,X,k3);
			// T5
			add(k3,n2,epsilon,B12,ldb,B11,ldb,Y,n2);
			// P5 (in V)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
			// C22 += (P5-P2)/e
			subscalacc(NoField,m2,n2,V,n2,C12,ldc,(double)1/epsilon,C22,ldc);
			// C12 = (P4-P2)/e
			subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C12,ldc);
			// S9
			add(m2,k3,epsilon,A22,lda,A23,lda,X,k3);
			// T9
			add(k3,n2,epsilon,B31,ldb,B21,ldb,Y,n2);
			// P9 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// C22 -= P9
			FFLAS::fsubin(NoField,m2,n2,U,n2,C22,ldc);
			// P7 (in C21)
			gemm_bini_232_2(F,m2,n2,k3,A23,lda,B21,ldb,C21,ldc,rec-1,epsilon);
			// C11 = C11  - P7/e
			add(m2,n2,-(double)1/epsilon,C21,ldc,C11,ldc,C11,ldc);
			// C21 =  (P9-P7)/e
			subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C21,ldc);
			// C21 += P5
			FFLAS::faddin(NoField,m2,n2,V,n2,C21,ldc);
			// S10
			add(m2,k3,-epsilon,A12,lda,A13,lda,X,k3);
			// T10
			add(k3,n2,epsilon,B31,ldb,B32,ldb,Y,n2);
			// P10 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// C12 += P10
			FFLAS::faddin(NoField,m2,n2,U,n2,C12,ldc);
			// C11 += P10/e
			add(m2,n2,(double)1/epsilon,U,n2,C11,ldc,C11,ldc);


			delete[] X ;
			delete[] Y ;
			delete[] U ;
			delete[] V ;


			return C;

		}

		template<class Field>
		double *
		gemm_bini_232_3_acc(const Field & F
				, const size_t m
				, const size_t n
				, const size_t k
				, const double *A , const size_t lda
				, const double *B , const size_t ldb
				, double *C , const size_t ldc
				, int rec
				, const double & epsilon
			       )
		{
			if (rec != 0)
				exit(-1);
			FFLAS::DoubleDomain R;
			FFLAS::fgemm(R,
			     FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
			     m,n,k,
			     1,
			     A,lda, B,ldb,
			     1,
			     C, ldc);


		}

		template<class Field>
		double *
		gemm_bini_232_3(const Field & F
				, const size_t m
				, const size_t n
				, const size_t k
				, const double *A , const size_t lda
				, const double *B , const size_t ldb
				, double *C , const size_t ldc
				, int rec
				, const double & epsilon
			       )
		{
			FFPACK::UnparametricField<double>   NoField;
			// const double p = (double)F.characteristic();
			size_t M = (n>m)?std::min(k,m):std::min(k,n);
			// Field G(p*p);

			if ( M < TRE  || rec <= 0) {
				// std::cout << "ffw" << std::endl;
				return gemm_fflas(F, m,n,k, A,lda, B,ldb, C, ldc);
				// return gemm_fflas(NoField, m,n,k, A,lda, B,ldb, C, ldc);
			}

			assert(k/3*3==k); // k divisible par 3
			assert(n/2*2==n); // n divisible par 2
			assert(m/2*2==m); // m divisible par 2

			// std::cout << "tested" << std::endl;

			size_t n2 = n/2;
			size_t k3 = k/3;
			size_t m2 = m/2;

			// std::cout << "€ = " << epsilon << std::endl;

			// sub matrices in B
			const double * B11 = B;
			const double * B12 = B   +n2;
			const double * B21 = B   +ldb*k3;
			const double * B22 = B21 +n2;
			const double * B31 = B21 +ldb*k3;
			const double * B32 = B31 +n2;

			// sub matrices in C
			double * C11 = C;
			double * C12 = C   +n2;
			double * C21 = C   +ldc*m2;
			double * C22 = C21 +n2;

			// sub matrices in A

			const double * A11 = A;
			const double * A12 = A   +k3;
			const double * A13 = A   +2*k3;
			const double * A21 = A   +lda*m2;
			const double * A22 = A21 +k3;
			const double * A23 = A21 +2*k3;


			FFLAS::fzero(F,m,n,C,ldc);

			/*
			 * Algo :
			 *
			 * S1  := A11  +A22*e;
			 * S3  := -(A11+A21);
			 * S4  := A11+A12*e;
			 * S5  := A21 - A22*e;
			 * S6  := A12*e  + A23;
			 * S8 := -(A13+A23):
			 * S9  := A22*e  + A23;
			 * S10 := -A12*e+A13;
			 *
			 * T1  := B11 +B22;
			 * T4  := e*B12+B22;
			 * T5  := B11 +e*B12;
			 * T6  := B21   +B32;
			 * T9  := B21   + e*B31;
			 * T10 := e*B31   +B32;
			 *
			 * P1 := Bini232(S1,T1 ,e);
			 * P2 := Bini232(A11,B22 ,e);
			 * P3 := Bini232(S3,B11,e);
			 * P4 := Bini232(S4,T4 ,e);
			 * P5 := Bini232(S5,T5 ,e);
			 * P6 := Bini232(S6,T6 ,e);
			 * P7 := Bini232(A23,B21 ,e);
			 * P8 := Bini232(S8,B32,e);
			 * P9 := Bini232(S9,T9 ,e);
			 * P10:= Bini232(S10,T10,e);
			 *
			 *
			 * C11 := evalm(P1-P4+(P6-P7+P8+P10)/e);
			 * C12 := evalm((-P2+P4)/e+P10) ;
			 * C21 := evalm(P5+(-P7+P9)/e) ;
			 * C22 := evalm((P1-P2+P3+P5)/e+P6-P9);
			 *
			 */

			// could be just one band for the scalings



			double * U = new double[m2*n2];
			double * V = new double[std::max(k3,m2)*n2];
			double * X = new double[m2*k3];
			double * Y = new double[k3*n2];

			// S1
			double * eA22 = new double[std::max(m2,n2)*k3];
			FFLAS::fscal(NoField,m2,k3,epsilon,A22,lda,eA22,k3);
			FFLAS::fadd(NoField,m2,k3,eA22,k3,A11,lda,X,k3);
			// T1
			FFLAS::fadd(NoField,k3,n2,B11,ldb,B22,ldb,Y,n2);
			// P1 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// S3
			negadd(m2,k3,A11,lda,A21,lda,X,k3);
			// P3 (in V)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,B11,ldb,V,n2,rec-1,epsilon);
			// C22 = (P1+P3)/e
			addscal(NoField,m2,n2,U,n2,V,n2,(double)1/epsilon,C22,ldc);
			// S6
			double * eA12 = new double[m2*k3];
			FFLAS::fscal(NoField,m2,k3,epsilon,A12,lda,eA12,k3);
			FFLAS::fadd(NoField,m2,k3,eA12,k3,A23,lda,X,k3);
			// T6
			FFLAS::fadd(NoField,k3,n2,B21,ldb,B32,ldb,Y,n2);
			// P6 (in V)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
			// C22 += P6
			FFLAS::faddin(NoField,m2,n2,V,n2,C22,ldc);
			// S8
			negadd(m2,k3,A13,lda,A23,lda,X,k3);
			// P8 (in C11)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,B32,ldb,C11,ldc,rec-1,epsilon);
			// C11 = (P8+P6)/e
			addscalinf(NoField,m2,n2,V,n2,(double)1/epsilon,C11,ldc);
			// C11 += P1
			FFLAS::faddin(NoField,m2,n2,U,n2,C11,ldc);
			// S4
			FFLAS::fadd(NoField,m2,k3,eA12,k3,A11,lda,X,k3);
			// T4
			double * eB12 = V ; // new double[n2*k3];
			FFLAS::fscal(NoField,k3,n2,epsilon,B12,ldb,eB12,n2);
			FFLAS::fadd(NoField,k3,n2,eB12,n2,B22,ldb,Y,n2);
			// P4 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// C11 -= P4
			FFLAS::fsubin(NoField,m2,n2,U,n2,C11,ldc);
			// P2 (in C12)
			gemm_bini_232_2(F,m2,n2,k3,A11,lda,B22,ldb,C12,ldc,rec-1,epsilon);
			// S5
			FFLAS::fsub(NoField,m2,k3,A21,lda,eA22,k3,X,k3);
			// T5
			FFLAS::fadd(NoField,k3,n2,eB12,n2,B11,ldb,Y,n2);
			// delete[] eB12;
			// P5 (in V)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,V,n2,rec-1,epsilon);
			// C22 += (P5-P2)/e
			subscalacc(NoField,m2,n2,V,n2,C12,ldc,(double)1/epsilon,C22,ldc);
			// C12 = (P4-P2)/e
			subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C12,ldc);
			// S9
			FFLAS::fadd(NoField,m2,k3,eA22,k3,A23,lda,X,k3);
			double * eB31 = eA22 ;
			FFLAS::fscal(NoField,k3,n2,epsilon,B31,ldb,eB31,n2);
			// T9
			FFLAS::fadd(NoField,k3,n2,eB31,n2,B21,ldb,Y,n2);
			// P9 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// C22 -= P9
			FFLAS::fsubin(NoField,m2,n2,U,n2,C22,ldc);
			// P7 (in C21)
			gemm_bini_232_2(F,m2,n2,k3,A23,lda,B21,ldb,C21,ldc,rec-1,epsilon);
			// C11 = C11  - P7/e
			add(m2,n2,-(double)1/epsilon,C21,ldc,C11,ldc,C11,ldc);
			// C21 =  (P9-P7)/e
			subscalinf(NoField,m2,n2,U,n2,-(double)1/epsilon,C21,ldc);
			// C21 += P5
			FFLAS::faddin(NoField,m2,n2,V,n2,C21,ldc);
			// S10
			FFLAS::fsub(NoField,m2,k3,A13,lda,eA12,k3,X,k3);
			delete[] eA12;
			// T10
			FFLAS::fadd(NoField,k3,n2,eB31,n2,B32,ldb,Y,n2);
			delete[] eA22;
			// P10 (in U)
			gemm_bini_232_2(F,m2,n2,k3,X,k3,Y,n2,U,n2,rec-1,epsilon);
			// C12 += P10
			FFLAS::faddin(NoField,m2,n2,U,n2,C12,ldc);
			// C11 += P10/e
			add(m2,n2,(double)1/epsilon,U,n2,C11,ldc,C11,ldc);


			delete[] X ;
			delete[] Y ;
			delete[] U ;
			delete[] V ;


			return C;

		}


	} // Rec


	template<class Field>
	typename Field::Element *
	gemm_bini_322_p(const Field &F
			, const size_t m
			, const size_t n
			, const size_t k
			, const typename Field::Element *A
			, const size_t lda
			, const typename Field::Element *B
			, const size_t ldb
			, typename Field::Element *C
			, const size_t ldc
			, int rec
			, size_t algo
		       )
	{

		assert(k/6*6==k); // k divisible par 6
		assert(n/6*6==n); // n divisible par 6
		assert(m/6*6==m); // m divisible par 6

		// e-formule
		double epsilon = (double) F.characteristic() ;
		switch(algo) {
		case 0 :
				Rec::gemm_bini_322_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
				FFLAS::finit_fuzzy(F,m,n,C,ldc);
				// FFLAS::finit(F,m,n,C,ldc);
				break;
		case 1 :
				Rec::gemm_bini_322_0(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
				FFLAS::finit_fuzzy(F,m,n,C,ldc);
				// FFLAS::finit(F,m,n,C,ldc);
				break;
		case 2 :
				Rec::gemm_bini_322_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
				FFLAS::finit_fuzzy(F,m,n,C,ldc);
				break;
		case 3 :
				Rec::gemm_bini_223_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
				FFLAS::finit_fuzzy(F,m,n,C,ldc);
				// FFLAS::finit(F,m,n,C,ldc);
				break;
		case 4 :
				Rec::gemm_bini_232_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
				FFLAS::finit_fuzzy(F,m,n,C,ldc);
				break;
		case 5 :
				Rec::gemm_bini_232_3(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
				FFLAS::finit_fuzzy(F,m,n,C,ldc);
				break;
		default :
			std::cout << " not an algo :" << algo << std::endl;;
			exit(-1);
		}



		return C;

	}

	template<class Field>
	typename Field::Element *
	gemm_bini_322_e(const Field &F
			, const size_t m
			, const size_t n
			, const size_t k
			, const typename Field::Element *A
			, const size_t lda
			, const typename Field::Element *B
			, const size_t ldb
			, typename Field::Element *C
			, const size_t ldc
			, int rec
			, size_t algo
		       )
	{

		assert(k/2*2==k); // k divisible par 2
		assert(n/2*2==n); // n divisible par 2
		assert(m/3*3==m); // m divisible par 3

		// e-formule
		double epsilon = 1./(1<<27);
		switch(algo) {
		case 0 :
			Rec::gemm_bini_322_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
			break;
		case 1 :
			Rec::gemm_bini_322_0(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
			break;
		case 2 :
			Rec::gemm_bini_322_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
			break;
		case 3 :
			Rec::gemm_bini_223_mem(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
			break;
		case 4 :
			Rec::gemm_bini_232_2(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
			break;
		case 5 :
			Rec::gemm_bini_232_3(F,m,n,k,A,lda,B,ldb,C,ldc,rec,epsilon);
			break;
		default :
			std::cout << " not an algo :" << algo << std::endl;;
			exit(-1);
		}


		// vire les e.
		// FFLAS::finit_fuzzy(F,m,n,C,ldc);
		FFLAS::finit_fuzzy(F,m,n,C,ldc);

		return C;

	}
	} // Protected
} // FFLAS



template<class Field>
void test_algo(const Field &F, size_t m, size_t n, size_t k
	       , const typename Field::Element * A, size_t lda
	       , const typename Field::Element * B, size_t ldb
	       , typename Field::Element * C, size_t ldc
	       , typename Field::Element * D, size_t ldd
	       , typename Field::Element * E, size_t lde
	       , int r, size_t algo, Timer & tim_e, Timer & tim_p, Timer & tom, int & ok_e, int & ok_p,
	       bool with_e)
{
	Timer tmp ;

	if (with_e) {
	tmp.clear();tmp.start();
	FFLAS::Protected::gemm_bini_322_e(F,m,n,k,A,k,B,n,C,n,r,algo);
	tmp.stop(); tim_e += tmp ;
	}


	tmp.clear();tmp.start();
	fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
	      m,n,k, 1, A,k, B,n, 0, D, n);
	tmp.stop(); tom += tmp ;

	tmp.clear();tmp.start();
	FFLAS::Protected::gemm_bini_322_p(F,m,n,k,A,k,B,n,E,n,r,algo);
	tmp.stop(); tim_p += tmp ;

	size_t faux = 0 ;

	if (with_e) {
		for (size_t i = 0 ; i < m ; ++i) {
			for (size_t j = 0 ; j < n ; ++j) {
				if (!F.areEqual(C[i*n+j],D[i*n+j])) {
					++faux ;
				}
			}
		}
		if (faux) {
			std::cout << "bini_e " << descr[algo] << " : bad/all = " << faux << '/' << m*n << " ~~  " << (double)faux/(double)(m*n) << std::endl;
		}
		else ok_e ++ ;


#if 1
		if (faux && (n<20)) {
			std::cout << "Diff" << std::endl;
			for (size_t i = 0 ; i < m ; ++i) {
				for (size_t j = 0 ; j < n ; ++j)
					std::cout << C[i*n+j]-D[i*n+j] << ' ';
				std::cout << std::endl;
			}
		}

#endif
	}


	faux = 0 ;
	for (size_t i = 0 ; i < m ; ++i) {
		for (size_t j = 0 ; j < n ; ++j) {
			if (!F.areEqual(D[i*n+j],E[i*n+j])) {
				++faux ;
			}
		}
	}
	if (faux) {
		std::cout << "bini_p " << descr[algo] << " : bad/all = " << faux << '/' << m*n << " ~~  " << (double)faux/(double)(m*n) << std::endl;
	}
	else ok_p ++ ;


#if 1
	if (faux && (n<20)) {
		std::cout << "Diff" << std::endl;
		for (size_t i = 0 ; i < m ; ++i) {
			for (size_t j = 0 ; j < n ; ++j)
				std::cout << D[i*n+j]-E[i*n+j] << ' ';
			std::cout << std::endl;
		}
	}
#endif
}

template<class T>
struct changeField {
	typedef T other ;
};

template<>
struct changeField<Modular<double> > {
typedef Modular<float> other;
};

template<>
struct changeField<ModularBalanced<double> > {
typedef ModularBalanced<float> other;
};

template<class Field>
void test(size_t m, size_t k, size_t n, size_t p, int r, bool with_e)
{

	typedef typename Field::Element Element;

	Element * A = new Element[m*k];
	Element * B = new Element[n*k];
	Element * C = new Element[m*n];
	Element * D = new Element[m*n];
	Element * E = new Element[m*n];


	Field F(p);
	F.write(std::cout<< " * Field " ) << std::endl;

	typedef typename changeField<Field>::other Field_f  ;
	typedef typename Field_f::Element Element_f ;
	Field_f F_f(p);
	Element_f * A_f = new Element_f[m*k];
	Element_f * B_f = new Element_f[n*k];
	Element_f * C_f = new Element_f[m*n];

#if defined(NOTRANDOM)
	size_t i0 ;
	size_t j0 ;
	Element p2 ; F.init(p2,(size_t)F.mOne/2);
	std::cout << p2 << std::endl;
#warning "not random"
	for (size_t i = 0 ; i < m ; ++i)
		for (size_t j = 0 ; j < k ; ++j) {
			i0 = i/(m/3);
			j0 = j/(k/2);
			if      (i0 == 0 and j0 == 0) A[i*k+j] = F.mOne ;
			else if (i0 == 0 and j0 == 1) A[i*k+j] = F.zero ;
			else if (i0 == 1 and j0 == 0) A[i*k+j] = F.mOne ;
			else if (i0 == 1 and j0 == 1) A[i*k+j] = F.mOne ;
			else if (i0 == 2 and j0 == 0) A[i*k+j] = F.mOne ;
			else if (i0 == 2 and j0 == 1) A[i*k+j] = F.mOne ;
			else A[i*k+j] = F.mOne ;
		}
	for (size_t i = 0 ; i < k ; ++i)
		for (size_t j = 0 ; j < n ; ++j) {
			i0 = i/(k/2);
			j0 = j/(n/2);
			if      (i0 == 0 and j0 == 0) B[i*n+j] = F.mOne ;
			else if (i0 == 0 and j0 == 1) B[i*n+j] = F.mOne ;
			else if (i0 == 1 and j0 == 0) B[i*n+j] = F.mOne ;
			else if (i0 == 1 and j0 == 1) B[i*n+j] = F.zero ;
			else B[i*n+j] = F.mOne ;

		}
#endif

	std::vector<Timer> tim_e(algos);
	std::vector<Timer> tim_p(algos);
	Timer tom; tom.clear();
	Timer tam; tam.clear();
	Timer tmp;
	for (size_t i = 0 ; i < algos ; ++i) {
		tim_e[i].clear();
		tim_p[i].clear();
	}

	size_t iters = 4 ;

	std::vector<int> ok_p(algos,0);
	std::vector<int> ok_e(algos,0);

	for (size_t b = 0 ; b < iters ; ++b) {
		std::cout << b+1 << " of " << iters << std::endl;
#if not defined(NOTRANDOM)
		RandomMatrix(F,A,m,k,k);
		RandomMatrix(F,B,k,n,n);
#endif
		FFLAS::finit(F_f,m,k,A,k,A_f,k);
		FFLAS::finit(F_f,k,n,B,n,B_f,n);

		tmp.clear();tmp.start();
		fgemm(F_f,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
		      m,n,k, 1, A_f,k, B_f,n, 0, C_f, n);
		tmp.stop(); tam += tmp ;

		for (size_t i = 0 ; i < algos ; ++i) {
			test_algo(F,m,n,k,A,k,B,n,C,n,D,n,E,n,r,selec[i],tim_e[i],tim_p[i],tom,ok_e[i],ok_p[i],with_e);
		}


	}
	std::cout  << std::endl << "results" << std::endl;

	int min_e = -1 ;
	double bini_e = -1 ;
	for (size_t i = 0 ; i < algos ; ++i){
		if (ok_e[i] == iters) {
			double bini1 = tim_e[i].usertime()/(double)ok_e[i] ;
			if (bini_e <  0) {
				bini_e = bini1;
				min_e = (int) i ;
			}
			else if (bini1 < bini_e) {
				min_e  = (int)i ;
				bini_e = bini1 ;
			}
		}
	}
	for (size_t i = 0 ; i < algos ; ++i){
		if (ok_e[i] == iters) {
		double bini1 = tim_e[i].usertime()/(double)ok_e[i] ;
		std::cout << "Bini_e ( " << descr[i] << " ) : " ;
		if ((int)i == min_e) std::cout << " * " ;
		else std::cout << "   ";
		std::cout << bini1  << 's'<<  std::endl;
	}
	}
		int min_p = -1 ;
		double bini_p = -1 ;
		for (size_t i = 0 ; i < algos ; ++i){
			if (ok_p[i] == iters) {
				double bini1 = tim_p[i].usertime()/(double)ok_p[i];
				if (bini_p <  0) {
					bini_p = bini1;
					min_p = (int)i ;
				}
				else if (bini1 < bini_p) {
					min_p  = (int)i ;
					bini_p = bini1 ;
				}

			}
		}
		for (size_t i = 0 ; i < algos ; ++i){
			if (ok_p[i] == iters) {
				double bini1 = tim_p[i].usertime()/(double)ok_p[i];
				std::cout << "Bini_p ( " << descr[i] << " ) : " ;
				if ((int)i == min_p) std::cout << " * " ;
				else std::cout << "   ";
				std::cout << bini1  << 's'<<  std::endl;
			}
		}

		std::cout << "Wino d : " << tom.usertime()/(double)(iters*algos) << 's'<<  std::endl;
		std::cout << "Wino f : " << tam.usertime()/(double)iters << 's'<<  std::endl;
		double wino =  std::min(tom.usertime()/(double)algos,tam.usertime())/(double)iters ;
		if (bini_e>=0)
			std::cout << "Gain e: " << ((bini_e-wino)/wino)*100 << '%' << std::endl;
		if (bini_p>=0)
			std::cout << "Gain p: " << ((bini_p-wino)/wino)*100 << '%' << std::endl;



		delete[] A ;
		delete[] B;
		delete[] C;
		delete[] D;
		delete[] E;

		delete[] A_f ;
		delete[] B_f;
		delete[] C_f;
	}



int main(int ac, char **av) {
	static size_t m = 36 ;
	static size_t n = 12 ;
	static size_t k = 18 ;
	static size_t p = 101;
	bool  eps = false ;
	int r = 1 ;
	int seed = (int) time(NULL);

	static Argument as[] = {
		{ 'p', "-p P", "Set the field characteristic.",  TYPE_INT , &p },
		{ 'n', "-n N", "Set the number of cols in C.",   TYPE_INT , &n },
		{ 'm', "-m N", "Set the number of rows in C.",   TYPE_INT , &m },
		{ 'k', "-k N", "Set the number of rows in B.",   TYPE_INT , &k },
		{ 'r', "-k N", "Set the recursive number Bini.", TYPE_INT , &r },
		{ 's', "-s N", "Set the seed                 .", TYPE_INT , &seed },
		{ 'e', "-e " , "epsilon                 .", TYPE_NONE , &eps },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(ac,av,as);

	srand(seed);
	srand48(seed);

	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	std::cout << "size: " << m << ',' << k << ',' << n << std::endl;
	std::cout << "mod : " << p << std::endl;
	std::cout << "rec : " << r << std::endl;
	std::cout << "seed: " << seed << std::endl;
	std::cout << "thre: " << TRE << std::endl;
	std::cout << "=====================================================" << std::endl;
	test<Modular<double> > (m,k,n,p,r,eps);
	std::cout << "=====================================================" << std::endl;
	test<ModularBalanced<double> > (m,k,n,p,r,eps);
	std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

	return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

