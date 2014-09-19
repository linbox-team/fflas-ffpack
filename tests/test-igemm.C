
#define SIMD_INT

#include "fflas-ffpack/fflas-ffpack.h"
#include "fflas-ffpack/fflas/fflas_igemm/igemm.h"
#include "fflas-ffpack/utils/Matio.h"


int main()
{
	typedef FFPACK::Modular<FFPACK::Integer> IField ;
	IField Z(std::pow(2,63));

	size_t m = 5 ;
	size_t n = 5 ;
	size_t k = 5 ;
	size_t lda = 5 ;
	size_t ldb = 5 ;
	size_t ldc = 5 ;

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
	// for (size_t i=0;i<m;++i)
		// for (size_t j=0;j<n;++j)
			// Rand.random(C[i*ldc+j]);
			// C[i*ldc+j] = rand() % 10;

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
	Cc= FFLAS::fflas_new(F,ldc,n);
	Aa= FFLAS::fflas_new(F,lda,k);
	Bb= FFLAS::fflas_new(F,ldb,n);

	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<k;++j)
			F.init(Aa[j*lda+i],A[i*lda+j]);
	for (size_t i=0;i<k;++i)
		for (size_t j=0;j<n;++j)
			F.init(Bb[j*ldb+i],B[i*ldb+j]);
	for (size_t i=0;i<m;++i)
		for (size_t j=0;j<n;++j)
			F.init(Cc[j*ldc+i],0UL);

	write_field(F,std::cout << "A", Aa, m, k, lda,true,true) <<std::endl;
	write_field(F,std::cout << "B", Bb, k, n, ldb,true,true) <<std::endl;

	FFLAS::igemm_(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, Aa, lda, Bb, ldb, 0, Cc, ldc);


	write_field(F,std::cout << "C", Cc, m, n, ldc,true,true) <<std::endl;

	FFLAS::fflas_delete(A);
	FFLAS::fflas_delete(B);
	// FFLAS::fflas_delete(C);

	FFLAS::fflas_delete(Aa);
	FFLAS::fflas_delete(Bb);
	// FFLAS::fflas_delete(Cc);



	return 0;
}

