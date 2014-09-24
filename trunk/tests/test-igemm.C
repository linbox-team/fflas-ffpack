
#define SIMD_INT

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/fflas/fflas_igemm/igemm.h"
#include "fflas-ffpack/utils/Matio.h"
#include "fflas-ffpack/utils/args-parser.h"


int main(int argc, char** argv)
{
	typedef FFPACK::Modular<FFPACK::Integer> IField ;
	IField Z(std::pow(2,63));

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


	size_t lda = k+rand() % 3 ; // k
	size_t ldb = n+rand() % 3 ; // n
	size_t ldc = n+rand() % 3 ; // n
	size_t ldA = m+rand() % 3 ; // m
	size_t ldB = k+rand() % 3 ; // k
	size_t ldC = m+rand() % 3 ; // m

	int seed=0;
	typename IField::RandIter Rand(Z,seed);

	IField::Element_ptr A,B,C;
	C= FFLAS::fflas_new(Z,m,ldc);
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
			C[i*ldc+j] = rand() % 10;

	write_field(Z,std::cout << "A", A, m, k, lda,true,false) <<std::endl;
	write_field(Z,std::cout << "B", B, k, n, ldb,true,false) <<std::endl;



	IField::Element alpha,beta ;
	alpha = Z.one ;
	beta = Z.zero ;

	FFLAS::fgemm(Z,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);

	write_field(Z,std::cout << "C", C, m, n, ldc,true,false) <<std::endl;
	std::cout << "==========================" << std::endl;


	typedef FFPACK::UnparametricField<int64_t> FField ;
	FField F ;

	FField::Element_ptr Cc,Aa,Bb;
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
			F.init(Cc[j*ldC+i],0UL);

	write_field(F,std::cout << "A", Aa, m, k, ldA,true,true) <<std::endl;
	write_field(F,std::cout << "B", Bb, k, n, ldB,true,true) <<std::endl;

	FFLAS::igemm_(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, Aa, ldA, Bb, ldB, 0, Cc, ldC);


	write_field(F,std::cout << "C", Cc, m, n, ldC,true,true) <<std::endl;

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(B);
	FFLAS::fflas_delete(C);

	FFLAS::fflas_delete(Aa);
	FFLAS::fflas_delete(Bb);
	FFLAS::fflas_delete(Cc);



	return 0;
}

