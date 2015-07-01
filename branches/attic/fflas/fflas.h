#ifdef LB_TRTR
	// BB
	/** @brief ftrtr: Triangular-Triangular matrix multiplication.
	 * \f$B \gets \alpha \mathrm{op}(A) \times B\f$ (for FFLAS_SIDE::FflasLeft)
	 * A and B are triangular, with B UpLo
	 * and op(A) = A, A^T according to TransA
	 * A and B can be (non)unit
	 *
	 */
	template<class Field>
	typename Field::Element* ftrtr (const Field& F, const FFLAS_SIDE Side,
					       const FFLAS_UPLO Uplo,
					       const FFLAS_TRANSPOSE TransA,
					       const FFLAS_DIAG ADiag,
					       const FFLAS_DIAG BDiag,
					       const size_t M,
					       const typename Field::Element alpha,
					       typename Field::Element * A, const size_t lda,
					       typename Field::Element * B, const size_t ldb);
#endif

	/** faddm.
	 * A <- A+op(B)
	 * with op(B) = B or B^T
	 * @bug not tested
	 */
	template<class Field>
	void faddm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb);

	/** faddm.
	 * C <- op(A)+op(B)
	 * with op(B) = B or B^T
	 * @bug not tested
	 */
	template<class Field>
	void faddm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const FFLAS_TRANSPOSE transB,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  const typename Field::Element * B, const size_t ldb,
			  typename Field::Element * C, const size_t ldc );

	/** fsubm.
	 * A <- A-op(B)
	 * with op(B) = B or B^T
	 * @bug not tested
	 */
	template<class Field>
	void fsubm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb) ;

	/** fsubm.
	 * C <- op(A)-op(B)
	 * with op(B) = B or B^T
	 * @bug not tested
	 */
	template<class Field>
	void fsubm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const FFLAS_TRANSPOSE transB,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  const typename Field::Element * B, const size_t ldb,
			  typename Field::Element * C, const size_t ldc );

#if 0
		template<class Element>
		class faddmTrans;
		template<class Element>
		class faddmNoTrans;
		template<class Element>
		class fsubmTrans;
		template<class Element>
		class fsubmNoTrans;
		template<class Element>
		class faddmTransTrans;
		template<class Element>
		class faddmNoTransTrans;
		template<class Element>
		class faddmTransNoTrans;
		template<class Element>
		class faddmNoTransNoTrans;
		template<class Element>
		class fsubmTransTrans;
		template<class Element>
		class fsubmNoTransTrans;
		template<class Element>
		class fsubmTransNoTrans;
		template<class Element>
		class fsubmNoTransNoTrans;
#endif

#if 0
		// BB : ça peut servir...
#ifdef LB_TRTR
		template <class Element>
		class ftrtrLeftUpperNoTransNonUnitNonUnit;
		template <class Element>
		class ftrtrLeftUpperNoTransUnitNonUnit;
		template <class Element>
		class ftrtrLeftUpperTransNonUnitNonUnit;
		template <class Element>
		class ftrtrLeftUpperTransUnitNonUnit;
		template <class Element>
		class ftrtrLeftLowerNoTransNonUnitNonUnit;
		template <class Element>
		class ftrtrLeftLowerNoTransUnitNonUnit;
		template <class Element>
		class ftrtrLeftLowerTransNonUnitNonUnit;
		template <class Element>
		class ftrtrLeftLowerTransUnitNonUnit;
		template <class Element>
		class ftrtrLeftUpperNoTransNonUnitUnit;
		template <class Element>
		class ftrtrLeftUpperNoTransUnitUnit;
		template <class Element>
		class ftrtrLeftUpperTransNonUnitUnit;
		template <class Element>
		class ftrtrLeftUpperTransUnitUnit;
		template <class Element>
		class ftrtrLeftLowerNoTransNonUnitUnit;
		template <class Element>
		class ftrtrLeftLowerNoTransUnitUnit;
		template <class Element>
		class ftrtrLeftLowerTransNonUnitUnit;
		template <class Element>
		class ftrtrLeftLowerTransUnitUnit;
		template <class Element>
		class ftrtrRightUpperNoTransNonUnitNonUnit;
		template <class Element>
		class ftrtrRightUpperNoTransUnitNonUnit;
		template <class Element>
		class ftrtrRightUpperTransNonUnitNonUnit;
		template <class Element>
		class ftrtrRightUpperTransUnitNonUnit;
		template <class Element>
		class ftrtrRightLowerNoTransNonUnitNonUnit;
		template <class Element>
		class ftrtrRightLowerNoTransUnitNonUnit;
		template <class Element>
		class ftrtrRightLowerTransNonUnitNonUnit;
		template <class Element>
		class ftrtrRightLowerTransUnitNonUnit;
		template <class Element>
		class ftrtrRightUpperNoTransNonUnitUnit;
		template <class Element>
		class ftrtrRightUpperNoTransUnitUnit;
		template <class Element>
		class ftrtrRightUpperTransNonUnitUnit;
		template <class Element>
		class ftrtrRightUpperTransUnitUnit;
		template <class Element>
		class ftrtrRightLowerNoTransNonUnitUnit;
		template <class Element>
		class ftrtrRightLowerNoTransUnitUnit;
		template <class Element>
		class ftrtrRightLowerTransNonUnitUnit;
		template <class Element>
		class ftrtrRightLowerTransUnitUnit;
#endif
#endif

