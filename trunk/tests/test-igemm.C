
#define SIMD_INT

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/fflas/fflas_igemm/igemm.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"

#define COL_MAJOR false


int test_igemm(size_t m, size_t n, size_t k, enum CBLAS_TRANSPOSE tA, enum CBLAS_TRANSPOSE tB)
{
	srand(time(NULL));
	typedef FFPACK::Modular<FFPACK::Integer> IField ;
	IField Z(std::pow(2,63));


	size_t lda = k;//+rand() % 3 ; // k
	size_t ldb = n;//+rand() % 3 ; // n
	size_t ldc = n;//+rand() % 3 ; // n
#if COL_MAJOR
	size_t ldA = m;//+rand() % 3 ; // m
	size_t ldB = k;//+rand() % 3 ; // k
	size_t ldC = m;//+rand() % 3 ; // m
#else
	size_t ldA = k;//+rand() % 3 ; // k
	size_t ldB = n;//+rand() % 3 ; // n
	size_t ldC = n;//+rand() % 3 ; // n
#endif

	int seed=0;
	typename IField::RandIter Rand(Z,seed);

	IField::Element_ptr A,B,C,D;
	C= FFLAS::fflas_new(Z,m,ldc);
	D= FFLAS::fflas_new(Z,m,n);
	A= FFLAS::fflas_new(Z,m,lda);
	B= FFLAS::fflas_new(Z,k,ldb);

	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<k;++j)
			// Rand.random(A[i*lda+j]);
			A[i*lda+j] = rand() % 10 ;
	for (size_t i=0;i<k;++i)
		for (size_t j=0;j<n;++j)
			// Rand.random(B[i*ldb+j]);
			B[i*ldb+j] = rand() %10;
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			// Rand.random(C[i*ldc+j]);
			D[i*n+j]=C[i*ldc+j] = 0 ; //rand() % 10;

	write_field(Z,std::cout << "A:=", A, m, k, lda,true,false) <<';' <<std::endl;
	write_field(Z,std::cout << "B:=", B, k, n, ldb,true,false) <<';' <<std::endl;



	IField::Element alpha,beta ;
	alpha = Z.one ;
	beta = Z.one ;

	FFLAS::fgemm(Z,(FFLAS::FFLAS_TRANSPOSE)tA,(FFLAS::FFLAS_TRANSPOSE)tB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);

	write_field(Z,std::cout << "C:=", C, m, n, ldc,true,false) <<';' <<std::endl;
	std::cout << "==========================" << std::endl;


	typedef FFPACK::UnparametricField<int64_t> FField ;
	FField F ;

	FField::Element_ptr Cc,Aa,Bb;
#if COL_MAJOR
	Cc= FFLAS::fflas_new(F,ldC,n);
	Aa= FFLAS::fflas_new(F,ldA,k);
	Bb= FFLAS::fflas_new(F,ldB,n);

	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<k;++j)
			F.init(Aa[j*ldA+i],A[i*lda+j]);
	for (size_t i=0;i<k;++i)
		for (size_t j=0;j<n;++j)
			F.init(Bb[j*ldB+i],B[i*ldb+j]);
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			F.init(Cc[j*ldC+i],D[i*n+j]);
#else
	Cc= FFLAS::fflas_new(F,m,ldC);
	Aa= FFLAS::fflas_new(F,m,ldA);
	Bb= FFLAS::fflas_new(F,k,ldB);

	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<k;++j)
			F.init(Aa[i*ldA+j],A[i*lda+j]);
	for (size_t i=0;i<k;++i)
		for (size_t j=0;j<n;++j)
			F.init(Bb[i*ldB+j],B[i*ldb+j]);
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			F.init(Cc[i*ldC+j],D[i*n+j]);

#endif

	write_field(F,std::cout << "A:=", Aa, m, k, ldA,true,COL_MAJOR)  <<';'<<std::endl;
	write_field(F,std::cout << "B:=", Bb, k, n, ldB,true,COL_MAJOR) <<';' <<std::endl;

	FField::Element a,b ;
	a=F.one;
	b=F.one;
#if 0
	FFLAS::igemm_(CblasColMajor, tA, tB, m, n, k, a, Aa, ldA, Bb, ldB, b, Cc, ldC);
#else
	FFLAS::igemm_(CblasRowMajor,tA,tB, m, n, k, a, Aa, ldA, Bb, ldB, b, Cc, ldC);
#endif


	write_field(F,std::cout << "C:=", Cc, m, n, ldC,true,COL_MAJOR) <<';' <<std::endl;

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(B);
	FFLAS::fflas_delete(C);

	FFLAS::fflas_delete(Aa);
	FFLAS::fflas_delete(Bb);
	FFLAS::fflas_delete(Cc);

	return true;
}

int main(int argc, char** argv)
{

	static size_t m = 9 ;
	static size_t n = 10 ;
	static size_t k = 11 ;

	static Argument as[] = {
		{ 'm', "-m m", "m.",      TYPE_INT , &m },
		{ 'n', "-n n", "n.",      TYPE_INT , &n },
		{ 'k', "-k k", "k.",      TYPE_INT , &k },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(argc,argv,as);


	test_igemm(m,n,k,CblasNoTrans, CblasNoTrans);
	std::cout << "xxxxxxxxxxxxxxxxxxxxxx" << std::endl;
	test_igemm(m,n,k,CblasTrans, CblasNoTrans);
	std::cout << "xxxxxxxxxxxxxxxxxxxxxx" << std::endl;
	test_igemm(m,n,k,CblasNoTrans, CblasTrans);
	std::cout << "xxxxxxxxxxxxxxxxxxxxxx" << std::endl;
	test_igemm(m,n,k,CblasTrans, CblasTrans);


	return 0;
}

