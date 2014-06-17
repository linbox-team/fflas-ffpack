


#ifdef __FFLASFFPACK_USE_KAAPI

	template<class Field>			
	struct Taskfgemm : public ka::Task<15>::Signature<
		Field,
		FFLAS_TRANSPOSE,
		FFLAS_TRANSPOSE,
		size_t ,
		size_t ,
		size_t ,
		typename Field::Element,
		ka::R<typename Field::Element>,
		size_t ,
		ka::R<typename Field::Element>,
		size_t ,
		typename Field::Element,
		ka::RW<typename Field::Element>,
		size_t ,
		size_t
		>{};
	template<class Field>			
	struct Taskfgemmw : public ka::Task<14>::Signature<
		Field,
		FFLAS_TRANSPOSE,
		FFLAS_TRANSPOSE,
		size_t ,
		size_t ,
		size_t ,
		typename Field::Element,
		ka::R<typename Field::Element>,
		size_t ,
		ka::R<typename Field::Element>,
		size_t ,
		typename Field::Element,
		ka::RW<typename Field::Element>,
		size_t 
		>{};
}

template<class Field>
struct TaskBodyCPU<FFLAS::Taskfgemm<Field> >{
	void operator()(const Field& F,
			const FFLAS::FFLAS_TRANSPOSE ta,
			const FFLAS::FFLAS_TRANSPOSE tb,
			size_t BlockRowDim,
			size_t BlockColDim,
			size_t k,
			const typename Field::Element alpha,
			ka::pointer_r<typename Field::Element> A,
			const size_t lda,
			ka::pointer_r<typename Field::Element> B,
			const size_t ldb,
			const typename Field::Element beta,
			ka::pointer_rw<typename Field::Element> C, const size_t ldc,
			const size_t w)
	{
		FFLAS::fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A.ptr(), lda, B.ptr() , ldb,
			      beta, C.ptr(), ldc, w);
	}
};

template<class Field>
struct TaskBodyCPU<FFLAS::Taskfgemmw<Field> >{
	void operator()(const Field& F,
			const FFLAS::FFLAS_TRANSPOSE ta,
			const FFLAS::FFLAS_TRANSPOSE tb,
			size_t BlockRowDim,
			size_t BlockColDim,
			size_t k,
			const typename Field::Element alpha,
			ka::pointer_r<typename Field::Element> A,
			const size_t lda,
			ka::pointer_r<typename Field::Element> B,
			const size_t ldb,
			const typename Field::Element beta,
			ka::pointer_rw<typename Field::Element> C, const size_t ldc)
	{
		FFLAS::fgemm( F, ta, tb, BlockRowDim, BlockColDim, k, alpha, A.ptr(), lda, B.ptr() , ldb,
			      beta, C.ptr(), ldc);
	}
};

#endif
