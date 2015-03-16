
#define SIMD_INT

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/fflas/fflas_igemm/igemm.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/timer.h"

// COL_MAJOR true not supported in test. To be updated.
#define COL_MAJOR false
#define LEAD_GEN true
#define DISPLAY false
#define TRUST_FGEMM false


int test_igemm(size_t m, size_t n, size_t k, enum CBLAS_TRANSPOSE tA, enum CBLAS_TRANSPOSE tB, int a_scal, int b_scal, bool timing)
{


	 FFLAS::Timer tim;

	srand((unsigned int)time(NULL));
	typedef Givaro::ModularBalanced<Givaro::Integer> IField ;
	IField Z(1UL<<63);

	size_t ra = (tA==CblasNoTrans) ? m : k ;
	size_t ca = (tA==CblasNoTrans) ? k : m ;
	size_t rb = (tB==CblasNoTrans) ? k : n;
	size_t cb = (tB==CblasNoTrans) ? n : k;

	size_t lda = ca ;
	size_t ldb = cb ; // n
	size_t ldc = n  ; // n

#if COL_MAJOR
	size_t ldA = m;//+rand() % 3 ; // m
	size_t ldB = k;//+rand() % 3 ; // k
	size_t ldC = m;//+rand() % 3 ; // m
#else
	size_t ldA = ca ; // k
	size_t ldB = cb ; // n
	size_t ldC = n  ; // n
#endif

#if LEAD_GEN
	lda += rand() % 5;
	ldb += rand() % 5;
	ldc += rand() % 5;
	ldA += rand() % 5;
	ldB += rand() % 5;
	ldC += rand() % 5;
#endif


	int seed=0;
	typename IField::RandIter Rand(Z,seed);
	// typename IField::RandIter Rand(Z,seed);

	IField::Element_ptr A,B,C,D;
	C= FFLAS::fflas_new(Z,m,ldc);
	D= FFLAS::fflas_new(Z,m,n);
	A= FFLAS::fflas_new(Z,ra,lda);
	B= FFLAS::fflas_new(Z,rb,ldb) ;

	for (size_t i=0;i<ra;++i)
		for (size_t j=0;j<ca;++j)
			// Rand.random(A[i*lda+j]);
			A[i*lda+j] = rand() % 10;
	for (size_t i=0;i<rb;++i)
		for (size_t j=0;j<cb;++j)
			// Rand.random(B[i*ldb+j]);
			B[i*ldb+j] = rand() % 10;
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			// Rand.random(C[i*ldc+j]);
			D[i*n+j]=C[i*ldc+j] = 0 ; //rand() % 10;

#if DISPLAY
	write_field(Z,std::cout << "A:=", A, (int)ra, (int)ca, (int)lda,true,false) <<';' <<std::endl;
	// write_field(Z,std::cout << "A:=", A, (int)ra, (int)ca, (int)lda,true,(tA==CblasTrans))  <<';'<<std::endl;
	write_field(Z,std::cout << "B:=", B, (int)rb, (int)cb, (int)ldb,true,false) <<';' <<std::endl;
	// write_field(Z,std::cout << "B:=", B, (int)rb, (int)cb, (int)ldb,true,(tB==CblasTrans)) <<';' <<std::endl;
#endif



	IField::Element alpha,beta ;
	alpha = (IField::Element) a_scal ;
	beta = (IField::Element) b_scal ;

	tim.clear(); tim.start();
#if TRUST_FGEMM
	FFLAS::fgemm(Z,(FFLAS::FFLAS_TRANSPOSE)tA,(FFLAS::FFLAS_TRANSPOSE)tB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
#else
	if (!timing) {
	IField::Element_ptr ail,blj;
	IField::Element tmp;
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j){
			Z.mulin(*(C+i*ldc+j),beta);
			Z.assign (tmp, Z.zero);
			for ( size_t l = 0; l < k ; ++l ){
				if ( tA == CblasNoTrans )
					ail = A+i*lda+l;
				else
					ail = A+l*lda+i;
				if ( tB == CblasNoTrans )
					blj = B+l*ldb+j;
				else
					blj = B+j*ldb+l;
				Z.axpyin (tmp, *ail, *blj);
			}
			Z.axpyin (*(C+i*ldc+j), alpha, tmp);
		}
	}
#endif
	tim.stop();
	if (timing) std::cout << "zgemm time: " << tim << std::endl;


#if DISPLAY
	write_field(Z,std::cout << "C:=", C, (int)m, (int)n, (int)ldc,true,false) <<';' <<std::endl;
	std::cout << ((tA == CblasTrans) ? "LinearAlgebra:-Transpose":"") << "(A).";
	std::cout << ((tB == CblasTrans) ? "LinearAlgebra:-Transpose":"") << "(B);" << std::endl;;
#endif

	std::cout << "---------------------------------------------" << std::endl;


	typedef Givaro::UnparametricRing<int64_t> FField ;
	FField F ;

	FField::Element_ptr Ci,Ai,Bi;
#if COL_MAJOR
	Ci= FFLAS::fflas_new(F,ldC,n);
	Ai= FFLAS::fflas_new(F,ldA,k);
	Bi= FFLAS::fflas_new(F,ldB,n);

	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<k;++j)
			F.init(Ai[j*ldA+i],A[i*lda+j]);
	for (size_t i=0;i<k;++i)
		for (size_t j=0;j<n;++j)
			F.init(Bi[j*ldB+i],B[i*ldb+j]);
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			F.init(Ci[j*ldC+i],D[i*n+j]);
#else
	Ci= FFLAS::fflas_new(F,m,ldC);
	Ai= FFLAS::fflas_new(F,ra,ldA);
	Bi= FFLAS::fflas_new(F,rb,ldB);

	for (size_t i=0;i<ra;++i)
		for (size_t j=0;j<ca;++j)
			F.init(Ai[i*ldA+j],A[i*lda+j]);
	for (size_t i=0;i<rb;++i)
		for (size_t j=0;j<cb;++j)
			F.init(Bi[i*ldB+j],B[i*ldb+j]);
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			F.init(Ci[i*ldC+j],D[i*n+j]);

#endif

#if DISPLAY
	write_field(F,std::cout << "A:=", Ai, (int)ra, (int)ca, (int)ldA,true,COL_MAJOR)  <<';'<<std::endl;
	// write_field(F,std::cout << "A:=", Ai, (int)ra, (int)ca, (int)ldA,true,(tA==CblasTrans))  <<';'<<std::endl;
	write_field(F,std::cout << "B:=", Bi, (int)rb, (int)cb, (int)ldB,true,COL_MAJOR) <<';' <<std::endl;
	// write_field(F,std::cout << "B:=", Bi, (int)rb, (int)cb, (int)ldB,true,(tB==CblasTrans)) <<';' <<std::endl;
#endif

	FField::Element a,b ;
	a= (FField::Element) a_scal;
	b= (FField::Element) b_scal;
#if 0
	FFLAS::igemm_(CblasColMajor, tA, tB, m, n, k, a, Ai, ldA, Bi, ldB, b, Ci, ldC);
#else
	tim.clear(); tim.start();
	FFLAS::igemm_(CblasRowMajor,tA,tB, (int)m, (int)n, (int)k, a, Ai, (int)ldA, Bi, (int)ldB, b, Ci, (int)ldC);
	tim.stop();
	if (timing) std::cout << "igemm time: " << tim << std::endl;
#endif


#if DISPLAY
	write_field(F,std::cout << "C:=", Ci, (int)m, (int)n, (int)ldC,true,COL_MAJOR) <<';' <<std::endl;
	std::cout << ((tA == CblasTrans) ? "LinearAlgebra:-Transpose":"") << "(A).";
	std::cout << ((tB == CblasTrans) ? "LinearAlgebra:-Transpose":"") << "(B);" << std::endl;;
#endif

	bool pass = true ;
	if (!timing) {
#if DISPLAY
	for (size_t i = 0 ; i < m ; ++i) {
		for (size_t j = 0 ; j < n  ; ++j) {
			if (Ci[i*ldC+j] != (typename IField::Element) C[i*ldc+j]) {
				pass = false;
				std::cout << 'x' ;
			}
			else
				std::cout << 'o' ;
		}
		std::cout << std::endl;
	}
#else
	for (size_t i = 0 ; i < m && pass; ++i) {
		for (size_t j = 0 ; j < n && pass ; ++j) {
			if (Ci[i*ldC+j] != (typename IField::Element) C[i*ldc+j]) {
				pass = false;
			}
		}
	}
#endif
	}

	if (!pass) {
		std::cout << "***       *** " << std::endl;
		std::cout << "*** error *** " << std::endl;
		std::cout << "***       *** " << std::endl;
	}
	else {
		std::cout << "+++ pass  +++" << std::endl;
	}

	if  (timing) {
		FFLAS::Timer tom;
		// Givaro::Modular<double> G(65537);
		Givaro::UnparametricRing<double> G;
		double af, bf ;
		G.init(af,alpha);
		G.init(bf,beta);

		double *Cf,*Af,*Bf;
		Cf= FFLAS::fflas_new(G,m,ldC);
		Af= FFLAS::fflas_new(G,ra,ldA);
		Bf= FFLAS::fflas_new(G,rb,ldB);

		for (size_t i=0;i<ra;++i)
			for (size_t j=0;j<ca;++j)
				G.init(Af[i*ldA+j],A[i*lda+j]);
		for (size_t i=0;i<rb;++i)
			for (size_t j=0;j<cb;++j)
				G.init(Bf[i*ldB+j],B[i*ldb+j]);
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				G.init(Cf[i*ldC+j],D[i*n+j]);


		tom.clear(); tom.start();
		FFLAS::fgemm(G,(FFLAS::FFLAS_TRANSPOSE)tA,(FFLAS::FFLAS_TRANSPOSE)tB,m,n,k,af,Af,ldA,Bf,ldB,bf,Cf,ldC);
		tom.stop();
		std::cout << "fgemm time: " << tom << std::endl;

		FFLAS::fflas_delete(Af);
		FFLAS::fflas_delete(Bf);
		FFLAS::fflas_delete(Cf);

	}

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(B);
	FFLAS::fflas_delete(C);

	FFLAS::fflas_delete(Ai);
	FFLAS::fflas_delete(Bi);
	FFLAS::fflas_delete(Ci);

	FFLAS::fflas_delete(D);

	return pass;
}

int main(int argc, char** argv)
{

	static size_t m = 9 ;
	static size_t n = 10 ;
	static size_t k = 11 ;
	static bool timing = false ;

	static Argument as[] = {
		{ 'm', "-m m", "m.",      TYPE_INT , &m },
		{ 'n', "-n n", "n.",      TYPE_INT , &n },
		{ 'k', "-k k", "k.",      TYPE_INT , &k },
		{ 't', "-timing", "Output timing"            , TYPE_NONE, &timing},
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);


	for (int i = -1 ; i < 2 ; ++i) {
		for (int j = -1 ; j < 2 ; ++j) {
			std::cout << "===================================================" << std::endl;
			std::cout << " C = " << i << " A B + " << j << "C" << std::endl;
			test_igemm(m,n,k,CblasNoTrans, CblasNoTrans,i,j,timing);
			std::cout << " C = " << i << " A tB + " << j << "C" << std::endl;
			test_igemm(m,n,k,CblasTrans, CblasNoTrans,i,j,timing);
			std::cout << " C = " << i << " tA B + " << j << "C" << std::endl;
			test_igemm(m,n,k,CblasNoTrans, CblasTrans,i,j,timing);
			std::cout << " C = " << i << " tA tB + " << j << "C" << std::endl;
			test_igemm(m,n,k,CblasTrans, CblasTrans,i,j,timing);
		}
	}
	for (size_t a = 0 ; a < 4 ; ++a) {
		int i = rand() % 25 ;
		int j = rand() % 25 ;
		if (rand()%2) i = -i ;
		if (rand()%2) j = -j ;
		std::cout << " C = " << i << " A B + " << j << "C" << std::endl;
		test_igemm(m,n,k,CblasNoTrans, CblasNoTrans,i,j,timing);
		std::cout << " C = " << i << " A tB + " << j << "C" << std::endl;
		test_igemm(m,n,k,CblasTrans, CblasNoTrans,i,j,timing);
		std::cout << " C = " << i << " tA B + " << j << "C" << std::endl;
		test_igemm(m,n,k,CblasNoTrans, CblasTrans,i,j,timing);
		std::cout << " C = " << i << " tA tB + " << j << "C" << std::endl;
		test_igemm(m,n,k,CblasTrans, CblasTrans,i,j,timing);

	}


	return 0;
}
