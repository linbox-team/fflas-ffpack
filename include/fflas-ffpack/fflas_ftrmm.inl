/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* fflas/fflas_ftrmm.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

//---------------------------------------------------------------------
// ftrmm: TRiangular Matrix Multiply
// Computes  B <- alpha.op(A).B,  B <- alpha.B.op(A)
// B is M*N, A is M*M if Side==FflasLeft, N*N if Side==FflasRight
// Warning :Assumes alpha==1
//---------------------------------------------------------------------
template<class Field>
inline void
FFLAS::ftrmm(const Field& F, const enum FFLAS_SIDE Side,
	     const enum FFLAS_UPLO Uplo, 
	     const enum FFLAS_TRANSPOSE TransA,
	     const enum FFLAS_DIAG Diag, 
	     const size_t M, const size_t N,
	     const typename Field::Element alpha,
	     typename Field::Element * A, const size_t lda,
	     typename Field::Element * B, const size_t ldb){
	
	if (!M || !N ) return; 
	
	typename Field::Element zero;
	F.init(zero, 0.0);
	size_t nmax = DotProdBound (F, 0, zero);
	if ( Side==FflasLeft ){
		if ( Uplo==FflasUpper){
			if (TransA == FflasNoTrans){
				ftrmmLeftUpNoTrans(F,Diag,M,N,A,lda,B,ldb, nmax);
			}
			else{
				ftrmmLeftUpTrans(F,Diag,M,N,A,lda,B,ldb, nmax);
			}
		}
		else{
			if (TransA == FflasNoTrans){
				ftrmmLeftLowNoTrans(F,Diag,M,N,A,lda,B,ldb, nmax);
			}
			else{
				ftrmmLeftLowTrans(F,Diag,M,N,A,lda,B,ldb, nmax);
			}
		}
	}
	else{
	if ( Uplo==FflasUpper){
			if (TransA == FflasNoTrans){
				ftrmmRightUpNoTrans(F,Diag,M,N,A,lda,B,ldb,nmax);
			}
			else{
				ftrmmRightUpTrans(F,Diag,M,N,A,lda,B,ldb, nmax);
			}
		}
		else{
			if (TransA == FflasNoTrans){
				ftrmmRightLowNoTrans(F,Diag,M,N,A,lda,B,ldb, nmax);
			}
			else{
				ftrmmRightLowTrans(F,Diag,M,N,A,lda,B,ldb, nmax);
			}
		}
	}
	if (!F.isOne(alpha))
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				F.mulin (*(B+i*ldb+j), alpha);

}
template<class Field>
inline void 
FFLAS::ftrmmLeftUpNoTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			   const size_t M, const size_t N,
			   const typename Field::Element * A, const size_t lda,
			   typename Field::Element * B, const size_t ldb, 
			   const size_t nmax) {
	callFtrmmLeftUpNoTrans<typename Field::Element> () (F,Diag,M,N,A,lda,B,ldb,nmax);
}
template <class Element>
class FFLAS::callFtrmmLeftUpNoTrans{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, 
			 const size_t nmax) {
		typename Field::Element one;
		F.init(one, 1.0);
		if (M < nmax) {
			DoubleDomain::Element * Ad = new DoubleDomain::Element[M*M];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD (F, Ad, M, A, lda, M, M);
			MatF2MatD (F, Bd, N, B, ldb, M, N);
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, M, Bd, N);
			delete[] Ad;
			MatD2MatF( F, B, ldb, Bd, N, M, N );
			delete[] Bd;
		} else{
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator()(F, Diag, Mup, N, A, lda, B, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, Mup, N, Mdown, one,
			       A+Mup, lda, B+Mup*ldb, ldb, one, B, ldb);
			this->operator()(F, Diag, Mdown, N, A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
		}
	}
};
	
template <>
class FFLAS::callFtrmmLeftUpNoTrans<double>{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, 
			 const size_t nmax) {
		if (M < nmax) {
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else{
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator()(F, Diag, Mup, N, A, lda, B, ldb, nmax);
			fgemm( F, FflasNoTrans, FflasNoTrans, Mup, N, Mdown, 1.0,
			       A+Mup, lda, B+Mup*ldb, ldb, 1.0, B, ldb);
			this->operator()(F, Diag, Mdown, N, A+Mup*(lda+1), lda, 
					    B+Mup*ldb, ldb, nmax);
		}
	}
};
	
template <>
class FFLAS::callFtrmmLeftUpNoTrans<float>{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, 
			 const size_t nmax) {
		if (M < nmax) {
			cblas_strmm (CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else{
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator()(F, Diag, Mup, N, A, lda, B, ldb, nmax);
			fgemm( F, FflasNoTrans, FflasNoTrans, Mup, N, Mdown, 1.0,
			       A+Mup, lda, B+Mup*ldb, ldb, 1.0, B, ldb);
			this->operator()(F, Diag, Mdown, N, A+Mup*(lda+1), lda, 
					    B+Mup*ldb, ldb, nmax);
		}
	}
};
	
template<class Field>
inline void
FFLAS::ftrmmLeftUpTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 const size_t nmax) 
{
	callFtrmmLeftUpTrans<typename Field::Element> () (F,Diag,M,N,A,lda,B,ldb,nmax);
}

template <class Element>
class FFLAS::callFtrmmLeftUpTrans{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 const size_t nmax) {
		typename Field::Element one;
		F.init(one, 1.0);
		if (M < nmax) {
			DoubleDomain::Element * Ad = new DoubleDomain::Element[M*M];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD (F, Ad, M, A, lda, M, M);
			MatF2MatD (F, Bd, N, B, ldb, M, N);
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasUpper, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, M, Bd, N);
			delete[] Ad;
			MatD2MatF( F, B, ldb, Bd, N, M, N );
			delete[] Bd;
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator ()(F, Diag, Mup, N, A, lda, B, ldb, nmax);
			fgemm( F, FflasTrans, FflasNoTrans, Mup, N, Mdown, one, 
			       A+Mup*lda, lda, B+Mup*ldb, ldb, one, B, ldb);
			this->operator () (F, Diag, Mdown, N, A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
		}
	}
};
	
template <>
class FFLAS::callFtrmmLeftUpTrans<double>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax) {
		if (M < nmax) {
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasUpper, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator () (F, Diag, Mup, N, A, lda, B, ldb, nmax);
			fgemm( F, FflasTrans, FflasNoTrans, Mup, N, Mdown, 1.0, 
			       A+Mup*lda, lda, B+Mup*ldb, ldb, 1.0, B, ldb);
			this->operator() (F, Diag, Mdown, N, A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
		}
	}
};
template <>
class FFLAS::callFtrmmLeftUpTrans<float>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax) {
		if (M < nmax) {
			cblas_strmm (CblasRowMajor, CblasLeft, CblasUpper, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator () (F, Diag, Mup, N, A, lda, B, ldb, nmax);
			fgemm( F, FflasTrans, FflasNoTrans, Mup, N, Mdown, 1.0, 
			       A+Mup*lda, lda, B+Mup*ldb, ldb, 1.0, B, ldb);
			this->operator() (F, Diag, Mdown, N, A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
		}
	}
};
	
template<class Field>
inline void
FFLAS::ftrmmLeftLowNoTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			    const size_t M, const size_t N,
			    const typename Field::Element * A, const size_t lda,
			    typename Field::Element * B, const size_t ldb, 
			    const size_t nmax){
	callFtrmmLeftLowNoTrans<typename Field::Element>()(F,Diag,M,N,A,lda,B,ldb,nmax);
}

template <class Element>
class FFLAS::callFtrmmLeftLowNoTrans{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, 
			 const size_t nmax){
		typename Field::Element one;
		F.init (one, 1.0);
		if ( M < nmax ){
			DoubleDomain::Element * Ad = new DoubleDomain::Element[M*M];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD( F, Ad, M, A, lda, M, M );
			MatF2MatD( F, Bd, N, B, ldb, M, N );
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, M, Bd, N);
			delete[] Ad;
			MatD2MatF (F, B, ldb, Bd, N, M, N);
			delete[] Bd;
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			ftrmmLeftLowNoTrans (F, Diag, Mdown, N, A+Mup*(lda+1), lda,
					     B+Mup*ldb, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
			       one, A+Mup*lda, lda, B, ldb, one, B+Mup*ldb, ldb);
			ftrmmLeftLowNoTrans (F, Diag, Mup, N, A, lda, B, ldb, nmax);
		}
	}
};
template<>
class FFLAS::callFtrmmLeftLowNoTrans<double>{
public:
	template<class Field>
	void operator ()(const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, const size_t nmax){
		if ( M < nmax ){
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator()(F, Diag, Mdown, N, A+Mup*(lda+1), lda,
					 B+Mup*ldb, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
			       1.0, A+Mup*lda, lda, B, ldb, 1.0, B+Mup*ldb, ldb);
			this->operator()(F, Diag, Mup, N, A, lda, B, ldb, nmax);
		}
	}
};

template<>
class FFLAS::callFtrmmLeftLowNoTrans<float>{
public:
	template<class Field>
	void operator ()(const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, const size_t nmax){
		if ( M < nmax ){
			cblas_strmm (CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator()(F, Diag, Mdown, N, A+Mup*(lda+1), lda,
					 B+Mup*ldb, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
			       1.0, A+Mup*lda, lda, B, ldb, 1.0, B+Mup*ldb, ldb);
			this->operator()(F, Diag, Mup, N, A, lda, B, ldb, nmax);
		}
	}
};
	
template<class Field>
inline void 
FFLAS::ftrmmLeftLowTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax){
	callFtrmmLeftLowTrans<typename Field::Element> () (F,Diag,M,N,A,lda,B,ldb,nmax);
}

template <class Element>
class FFLAS::callFtrmmLeftLowTrans{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax){
		typename Field::Element one;
		F.init(one, 1.0);
		if ( M < nmax ) {
			DoubleDomain::Element * Ad = new DoubleDomain::Element[M*M];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD( F, Ad, M, A, lda, M, M );
			MatF2MatD( F, Bd, N, B, ldb, M, N );
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasLower, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, M, Bd, N);
			delete[] Ad;
			MatD2MatF (F, B, ldb, Bd, N, M, N);
			delete[] Bd;
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			this->operator()(F, Diag, Mdown, N, A+Mup*(lda+1), lda,
					   B+Mup*ldb, ldb, nmax);
			fgemm (F, FflasTrans, FflasNoTrans, Mdown, N, Mup,
			       one, A+Mup, lda, B, ldb, one, B+Mup*ldb, ldb);
			this->operator()( F, Diag, Mup, N, A, lda, B, ldb, nmax);
		}
	}
};
	
template <>
class FFLAS::callFtrmmLeftLowTrans<double>{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 const size_t nmax){
		if ( M < nmax ) {
			cblas_dtrmm (CblasRowMajor, CblasLeft, CblasLower, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			ftrmmLeftLowTrans (F, Diag, Mdown, N, A+Mup*(lda+1), lda,
					   B+Mup*ldb, ldb, nmax);
			fgemm (F, FflasTrans, FflasNoTrans, Mdown, N, Mup,
			       1.0, A+Mup, lda, B, ldb, 1.0, B+Mup*ldb, ldb);
			ftrmmLeftLowTrans( F, Diag, Mup, N, A, lda, B, ldb, nmax);
		}
	}
};
template <>
class FFLAS::callFtrmmLeftLowTrans<float>{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 const size_t nmax){
		if ( M < nmax ) {
			cblas_strmm (CblasRowMajor, CblasLeft, CblasLower, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Mup=M>>1;
			size_t Mdown = M-Mup;
			ftrmmLeftLowTrans (F, Diag, Mdown, N, A+Mup*(lda+1), lda,
					   B+Mup*ldb, ldb, nmax);
			fgemm (F, FflasTrans, FflasNoTrans, Mdown, N, Mup,
			       1.0, A+Mup, lda, B, ldb, 1.0, B+Mup*ldb, ldb);
			ftrmmLeftLowTrans( F, Diag, Mup, N, A, lda, B, ldb, nmax);
		}
	}
};

template<class Field>
inline void 
FFLAS::ftrmmRightUpNoTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			    const size_t M, const size_t N,
			    const typename Field::Element * A, const size_t lda,
			    typename Field::Element * B, const size_t ldb, 
			    const size_t nmax) {
	callFtrmmRightUpNoTrans<typename Field::Element>()(F,Diag,M,N,A,lda,B,ldb,nmax);
}

template <class Element>
class FFLAS::callFtrmmRightUpNoTrans{
public:
	
	template<class Field>
	void operator ()(const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, 
			 const size_t nmax) {
		typename Field::Element one;
		F.init(one, 1.0);
		if (N < nmax) {
			DoubleDomain::Element * Ad = new DoubleDomain::Element[N*N];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD( F, Ad, N, A, lda, N, N );
			MatF2MatD( F, Bd, N, B, ldb, M, N );
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, N, Bd, N);
			delete[] Ad;
			MatD2MatF (F, B, ldb, Bd, N, M, N);
			delete[] Bd;
		} else {
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator ()(F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					     B+Nup, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
			       one, B, ldb, A+Nup, lda, one, B+Nup, ldb);
			this->operator ()( F, Diag, M, Nup, A, lda, B, ldb, nmax);
		}
	}
};

template <>
class FFLAS::callFtrmmRightUpNoTrans<double>{
public:
	
	template<class Field>
	void operator()(const Field& F, const enum FFLAS_DIAG Diag, 
			const size_t M, const size_t N,
			const typename Field::Element * A, const size_t lda,
			typename Field::Element * B, const size_t ldb, const size_t nmax){
		if ( N < nmax ){
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
				     (CBLAS_DIAG)Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init(*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator ()( F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					     B+Nup, ldb, nmax);
			fgemm( F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
			       1.0, B, ldb, A+Nup, lda, 1.0, B+Nup, ldb);
			this->operator ()( F, Diag, M, Nup, A, lda, B, ldb, nmax);
		}
	}
};
template <>
class FFLAS::callFtrmmRightUpNoTrans<float>{
public:
	
	template<class Field>
	void operator()(const Field& F, const enum FFLAS_DIAG Diag, 
			const size_t M, const size_t N,
			const typename Field::Element * A, const size_t lda,
			typename Field::Element * B, const size_t ldb, const size_t nmax){
		if ( N < nmax ){
			cblas_strmm (CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
				     (CBLAS_DIAG)Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init(*(B+i*ldb+j),*(B+i*ldb+j));
		} else {
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator ()( F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					     B+Nup, ldb, nmax);
			fgemm( F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
			       1.0, B, ldb, A+Nup, lda, 1.0, B+Nup, ldb);
			this->operator ()( F, Diag, M, Nup, A, lda, B, ldb, nmax);
		}
	}
};

template<class Field>
inline void
FFLAS::ftrmmRightUpTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax){
	callFtrmmRightUpTrans<typename Field::Element>() (F,Diag,M,N,A,lda,B,ldb,nmax);
}

template <class Element>
class FFLAS::callFtrmmRightUpTrans{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 const size_t nmax){
		typename Field::Element one;
		F.init(one, 1.0);
		if (N < nmax) {
			DoubleDomain::Element * Ad = new DoubleDomain::Element[N*N];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD( F, Ad, N, A, lda, N, N );
			MatF2MatD( F, Bd, N, B, ldb, M, N );
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasUpper, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, N, Bd, N);
			delete[] Ad;
			MatD2MatF (F, B, ldb, Bd, N, M, N);
			delete[] Bd;
		}
		else{	
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator()( F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					   B+Nup, ldb, nmax);
			fgemm( F, FflasNoTrans, FflasTrans, M, Ndown, Nup, 
			       one, B, ldb, A+Nup*lda, lda, one, B+Nup, ldb);
			this->operator()( F, Diag, M, Nup, A, lda, B, ldb, nmax);
		}
	}
};

template <>
class FFLAS::callFtrmmRightUpTrans<double>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax){
		if (N < nmax) {
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasUpper, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		}
		else{
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator()( F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					   B+Nup, ldb, nmax);
			fgemm( F, FflasNoTrans, FflasTrans, M, Ndown, Nup, 
			       1.0, B, ldb, A+Nup*lda, lda, 1.0, B+Nup, ldb);
			this->operator()( F, Diag, M, Nup, A, lda, B, ldb, nmax);
		}
	}
};
template <>
class FFLAS::callFtrmmRightUpTrans<float>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax){
		if (N < nmax) {
			cblas_strmm (CblasRowMajor, CblasRight, CblasUpper, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init (*(B+i*ldb+j),*(B+i*ldb+j));
		}
		else{
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator()( F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					   B+Nup, ldb, nmax);
			fgemm( F, FflasNoTrans, FflasTrans, M, Ndown, Nup, 
			       1.0, B, ldb, A+Nup*lda, lda, 1.0, B+Nup, ldb);
			this->operator()( F, Diag, M, Nup, A, lda, B, ldb, nmax);
		}
	}
};
	
template<class Field>
inline void
FFLAS::ftrmmRightLowNoTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			     const size_t M, const size_t N,
			     const typename Field::Element * A, const size_t lda,
			     typename Field::Element * B, const size_t ldb, 
			     const size_t nmax){
	callFtrmmRightLowNoTrans<typename Field::Element>() (F,Diag,M,N,A,lda,B,ldb,nmax);
}
template <class Element>
class FFLAS::callFtrmmRightLowNoTrans{
public:
	template<class Field>
	void operator() (const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, 
			 const size_t nmax){
		typename Field::Element one;
		F.init(one, 1.0);
		if (N < nmax) {
			DoubleDomain::Element * Ad = new DoubleDomain::Element[N*N];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD( F, Ad, N, A, lda, N, N );
			MatF2MatD( F, Bd, N, B, ldb, M, N );
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasLower, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, N, Bd, N);
			delete[] Ad;
			MatD2MatF (F, B, ldb, Bd, N, M, N);
			delete[] Bd;
		} else {
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator() (F, Diag, M, Nup, A, lda, B, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, M, Nup, Ndown,
			       one, B+Nup, ldb, A+Nup*lda, lda, one, B, ldb);
			this->operator() (F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					      B+Nup, ldb, nmax);
		}
	}
};

template <>
class FFLAS::callFtrmmRightLowNoTrans<double>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax){
		typename Field::Element one;
		F.init(one, 1.0);
		if (N < nmax) {
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasLower, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init(*(B+i*ldb+j),*(B+i*ldb+j));			
			
		}
		else{
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator() (F, Diag, M, Nup, A, lda, B, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, M, Nup, Ndown,
			       one, B+Nup, ldb, A+Nup*lda, lda, one, B, ldb);
			this->operator() (F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					      B+Nup, ldb, nmax);
		}
	}
};
template <>
class FFLAS::callFtrmmRightLowNoTrans<float>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb,
			  const size_t nmax){
		typename Field::Element one;
		F.init(one, 1.0);
		if (N < nmax) {
			cblas_strmm (CblasRowMajor, CblasRight, CblasLower, CblasNoTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init(*(B+i*ldb+j),*(B+i*ldb+j));			
			
		}
		else{
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator() (F, Diag, M, Nup, A, lda, B, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasNoTrans, M, Nup, Ndown,
			       one, B+Nup, ldb, A+Nup*lda, lda, one, B, ldb);
			this->operator() (F, Diag, M, Ndown, A+Nup*(lda+1), lda, 
					      B+Nup, ldb, nmax);
		}
	}
};

template<class Field>
inline void
FFLAS::ftrmmRightLowTrans (const Field& F, const enum FFLAS_DIAG Diag, 
			   const size_t M, const size_t N,
			   const typename Field::Element * A, const size_t lda,
			   typename Field::Element * B, const size_t ldb, 
			   const size_t nmax){
	callFtrmmRightLowTrans<typename Field::Element>() (F,Diag,M,N,A,lda,B,ldb,nmax);
}

template <class Element>
class FFLAS::callFtrmmRightLowTrans{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb, 
			  const size_t nmax){
	
		typename Field::Element one;
		F.init(one, 1.0);
		if (N < nmax) {
			DoubleDomain::Element * Ad = new DoubleDomain::Element[N*N];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
			MatF2MatD( F, Ad, N, A, lda, N, N );
			MatF2MatD( F, Bd, N, B, ldb, M, N );
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasLower, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, Ad, N, Bd, N);
			delete[] Ad;
			MatD2MatF (F, B, ldb, Bd, N, M, N);
			delete[] Bd;
		} else{
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator() (F, Diag, M, Nup, A, lda, B, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasTrans, M, Nup, Ndown, one, 
			       B+Nup, ldb, A+Nup, lda, one, B, ldb);
			this->operator() (F, Diag, M, Ndown, A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
		}
	}
};
template<>
class FFLAS::callFtrmmRightLowTrans<double>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb, 
			  const size_t nmax){
	
		if (N < nmax) {
			cblas_dtrmm (CblasRowMajor, CblasRight, CblasLower, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init(*(B+i*ldb+j),*(B+i*ldb+j));			
			
		} else{
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator() (F, Diag, M, Nup, A, lda, B, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasTrans, M, Nup, Ndown, 1.0, 
			       B+Nup, ldb, A+Nup, lda, 1.0, B, ldb);
			this->operator() (F, Diag, M, Ndown, A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
		}
	}
};
template<>
class FFLAS::callFtrmmRightLowTrans<float>{
public:
	template<class Field>
	void operator () (const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb, 
			  const size_t nmax){
	
		if (N < nmax) {
			cblas_strmm (CblasRowMajor, CblasRight, CblasLower, CblasTrans,
				     (CBLAS_DIAG) Diag, M, N, 1.0, A, lda, B, ldb);
			for (size_t i=0; i< M; ++i)
				for (size_t j=0; j<N; ++j)
					F.init(*(B+i*ldb+j),*(B+i*ldb+j));			
			
		} else{
			size_t Nup=N>>1;
			size_t Ndown = N-Nup;
			this->operator() (F, Diag, M, Nup, A, lda, B, ldb, nmax);
			fgemm (F, FflasNoTrans, FflasTrans, M, Nup, Ndown, 1.0, 
			       B+Nup, ldb, A+Nup, lda, 1.0, B, ldb);
			this->operator() (F, Diag, M, Ndown, A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
		}
	}
};


