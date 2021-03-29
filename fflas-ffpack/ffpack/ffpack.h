/* ffpack.h
 * Copyright (C) 2005 Clement Pernet
 *               2014 FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file ffpack.h
 * \brief Set of elimination based routines for dense linear algebra.
 * Matrices are supposed over finite prime field of characteristic less than 2^26.

*/

#ifndef __FFLASFFPACK_ffpack_H
#define __FFLASFFPACK_ffpack_H

#include "givaro/givpoly1.h"
#include <fflas-ffpack/fflas-ffpack-config.h>

#ifdef __FFLASFFPACK_USE_OPENMP
#include <omp.h>
#endif

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_helpers.inl"

//#include "parallel.h"
#include <list>
#include <vector>
#include <iostream> // std::cout
#include <algorithm>

#define  __FFLASFFPACK_FTRSTR_THRESHOLD 64
#define  __FFLASFFPACK_FTRSSYR2K_THRESHOLD 64

/** @brief <b>F</b>inite <b>F</b>ield <b>PACK</b>
 * Set of elimination based routines for dense linear algebra.
 *
 * This namespace enlarges the set of BLAS routines of the class FFLAS, with higher
 * level routines based on elimination.
 \ingroup ffpack
 */
namespace FFPACK  { /* tags */

    /* \cond */
    enum FFPACK_LU_TAG
    {
        FfpackSlabRecursive = 1,
        FfpackTileRecursive = 2,
        FfpackSingular = 3,
        FfpackGaussJordanSlab = 4,
        FfpackGaussJordanTile = 5
    };

    enum FFPACK_CHARPOLY_TAG
    {
        FfpackAuto = 0,
        FfpackDanilevski = 1,
        FfpackLUK = 2,
        FfpackArithProgKrylovPrecond = 3,
        FfpackArithProg = 4,
        FfpackKG = 5,
        FfpackKGFast = 6,
        FfpackHybrid = 7,
        FfpackKGFastG = 8
    };
    /* \endcond */
    class CharpolyFailed{};

    /* \cond */
    enum FFPACK_MINPOLY_TAG
    {
        FfpackDense=1,
        FfpackKGF=2
    };
    /* \endcond */

}
namespace FFPACK { /* Permutations */

    /*****************/
    /* PERMUTATIONS  */
    /*****************/


    void LAPACKPerm2MathPerm (size_t * MathP, const size_t * LapackP,
                              const size_t N);

    void MathPerm2LAPACKPerm (size_t * LapackP, const size_t * MathP,
                              const size_t N);

    /* \cond */
    template <class Field>
    void MatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda, const size_t width,
                       const size_t M2,
                       const size_t R1, const size_t R2,
                       const size_t R3, const size_t R4);

    template <class Field>
    void MatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda,
                       const size_t width, const size_t M2,
                       const size_t R1, const size_t R2,
                       const size_t R3, const size_t R4,
                       const FFLAS::ParSeqHelper::Sequential seq);

    template <class Field, class Cut, class Param>
    void MatrixApplyS (const Field& F, typename Field::Element_ptr A, const size_t lda,
                       const size_t width, const size_t M2,
                       const size_t R1, const size_t R2,
                       const size_t R3, const size_t R4,
                       const FFLAS::ParSeqHelper::Parallel<Cut, Param> par);

    template <class Element>
    void PermApplyS (Element* A, const size_t lda, const size_t width,
                     const size_t M2,
                     const size_t R1, const size_t R2,
                     const size_t R3, const size_t R4);

    template <class Field>
    void MatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda, const size_t width,
                       const size_t N2,
                       const size_t R1, const size_t R2,
                       const size_t R3, const size_t R4);

    template <class Field>
    void MatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda,
                       const size_t width, const size_t N2,
                       const size_t R1, const size_t R2,
                       const size_t R3, const size_t R4,
                       const FFLAS::ParSeqHelper::Sequential seq);

    template <class Field, class Cut, class Param>
    void MatrixApplyT (const Field& F, typename Field::Element_ptr A, const size_t lda,
                       const size_t width, const size_t N2,
                       const size_t R1, const size_t R2,
                       const size_t R3, const size_t R4,
                       const FFLAS::ParSeqHelper::Parallel<Cut, Param> par);

    template <class Element>
    void PermApplyT (Element* A, const size_t lda, const size_t width,
                     const size_t N2,
                     const size_t R1, const size_t R2,
                     const size_t R3, const size_t R4);
    /* \endcond */

    /**
     * @brief Computes P1 x Diag (I_R, P2) where P1 is a LAPACK and P2 a LAPACK permutation
     * and store the result in P1 as a LAPACK permutation
     * @param [inout] P1 a LAPACK permutation of size N
     * @param P2 a LAPACK permutation of size N-R
     */

    /* \cond */
    inline void composePermutationsLLL (size_t * P1,
                                        const size_t * P2,
                                        const size_t R, const size_t N);

    /**
     * @brief Computes P1 x Diag (I_R, P2) where P1 is a LAPACK and P2 a LAPACK permutation
     * and store the result in MathP as a MathPermutation format.
     * @param [out] MathP a MathPermutation of size N
     * @param P1 a LAPACK permutation of size N
     * @param P2 a LAPACK permutation of size N-R
     */
    inline void composePermutationsLLM (size_t * MathP,
                                        const size_t * P1,
                                        const size_t * P2,
                                        const size_t R, const size_t N);

    /**
     * @brief Computes MathP1 x Diag (I_R, P2) where MathP1 is a MathPermutation and P2 a LAPACK permutation
     * and store the result in MathP1 as a MathPermutation format.
     * @param [inout] MathP1 a MathPermutation of size N
     * @param P2 a LAPACK permutation of size N-R
     */
    inline void composePermutationsMLM (size_t * MathP1,
                                        const size_t * P2,
                                        const size_t R, const size_t N);

    void cyclic_shift_mathPerm (size_t * P,  const size_t s);
    template<typename Base_t>
    void cyclic_shift_row_col(Base_t * A, size_t m, size_t n, size_t lda);
    template<class Field>
    void cyclic_shift_row(const Field& F, typename Field::Element_ptr A, size_t m, size_t n, size_t lda);
    template<class Field>
    void cyclic_shift_col(const Field& F, typename Field::Element_ptr A, size_t m, size_t n, size_t lda);
    /* \endcond */

    /**
     * @brief Applies a permutation P to the matrix A.
     * Apply a permutation P, stored in the LAPACK format (a sequence of transpositions)
     * between indices ibeg and iend of P to (iend-ibeg) vectors of size M stored in A (as column for NoTrans and rows for Trans).
     * Side==FFLAS::FflasLeft for row permutation Side==FFLAS::FflasRight for a column permutation
     * Trans==FFLAS::FflasTrans for the inverse permutation of P
     * @param F base field
     * @param Side  decides if rows (FflasLeft) or columns (FflasRight) are permuted
     * @param Trans decides if the matrix is seen as columns (FflasTrans) or rows (FflasNoTrans)
     * @param M size of the elements to permute
     * @param ibeg first index to consider in P
     * @param iend last index to consider in P
     * @param A input matrix
     * @param lda leading dimension of A
     * @param P permutation in LAPACK format
     * @param psh (optional): a sequential or parallel helper, to choose between sequential or parallel execution
     * @warning not sure the submatrix is still a permutation and the one we expect in all cases... examples for iend=2, ibeg=1 and P=[2,2,2]
     */
    template<class Field>
    void applyP( const Field& F,
                 const FFLAS::FFLAS_SIDE Side,
                 const FFLAS::FFLAS_TRANSPOSE Trans,
                 const size_t M, const size_t ibeg, const size_t iend,
                 typename Field::Element_ptr A, const size_t lda, const size_t * P );

    template<class Field>
    void applyP( const Field& F,
                 const FFLAS::FFLAS_SIDE Side,
                 const FFLAS::FFLAS_TRANSPOSE Trans,
                 const size_t m, const size_t ibeg, const size_t iend,
                 typename Field::Element_ptr A, const size_t lda, const size_t * P,
                 const FFLAS::ParSeqHelper::Sequential seq);

    template<class Field, class Cut, class Param>
    void applyP( const Field& F,
                 const FFLAS::FFLAS_SIDE Side,
                 const FFLAS::FFLAS_TRANSPOSE Trans,
                 const size_t m, const size_t ibeg, const size_t iend,
                 typename Field::Element_ptr A, const size_t lda, const size_t * P,
                 const FFLAS::ParSeqHelper::Parallel<Cut, Param> par);

    /** Apply a R-monotonically increasing permutation P, to the matrix A.
     * The permutation represented by P is defined as follows:
     *  - the first R values of P is a LAPACK reprensentation (a sequence of transpositions)
     *  - the remaining iend-ibeg-R values of the permutation are in a monotonically increasing progression
     * Side==FFLAS::FflasLeft for row permutation Side==FFLAS::FflasRight for a column permutation
     * Trans==FFLAS::FflasTrans for the inverse permutation of P
     * @param F	base field
     * @param Side selects if it is a row (FflasLeft) or column (FflasRight) permutation
     * @param Trans inverse permutation (FflasTrans/NoTrans)
     * @param M
     * @param ibeg
     * @param iend
     * @param A input matrix
     * @param lda leading dimension of A
     * @param P LAPACK permuation
     * @param R first values of P
     */
    template<class Field>
    void
    MonotonicApplyP (const Field& F,
                     const FFLAS::FFLAS_SIDE Side,
                     const FFLAS::FFLAS_TRANSPOSE Trans,
                     const size_t M, const size_t ibeg, const size_t iend,
                     typename Field::Element_ptr A, const size_t lda, const size_t * P, const size_t R);
    /* \cond */
    template<class Field>
    void
    MonotonicCompress (const Field& F,
                       const FFLAS::FFLAS_SIDE Side,
                       const size_t M,
                       typename Field::Element_ptr A, const size_t lda, const size_t incA, const size_t * P,
                       const size_t R, const size_t maxpiv, const size_t rowstomove,
                       const std::vector<bool> &ispiv);
    template<class Field>
    void
    MonotonicCompressMorePivots (const Field& F, const FFLAS::FFLAS_SIDE Side, const size_t M,
                                 typename Field::Element_ptr A, const size_t lda, const size_t incA,
                                 const size_t * MathP, const size_t R, const size_t rowstomove, const size_t lenP);
    template<class Field>
    void
    MonotonicCompressCycles (const Field& F, const FFLAS::FFLAS_SIDE Side, const size_t M,
                             typename Field::Element_ptr A, const size_t lda, const size_t incA,
                             const size_t * MathP, const size_t lenP);

    template<class Field>
    void
    MonotonicExpand (const Field& F, const FFLAS::FFLAS_SIDE Side, const size_t M,
                     typename Field::Element_ptr A, const size_t lda, const size_t incA,
                     const size_t * MathP, const size_t R, const size_t maxpiv,
                     const size_t rowstomove, const std::vector<bool> &ispiv);
    /* \endcond */

} // FFPACK permutations
// #include "ffpack_permutation.inl"

namespace FFPACK { /* fgetrs, fgesv */

    /** Solve the system \f$A X = B\f$ or \f$X A = B\f$.
     * Solving using the \c PLUQ decomposition of \p A
     * already computed inplace with \c PLUQ (FFLAS::FflasNonUnit).
     * Version for A square.
     * If A is rank deficient, a solution is returned if the system is consistent,
     * Otherwise an info is 1
     *
     * @param F base field
     * @param Side Determine wheter the resolution is left (FflasLeft) or right (FflasRight) looking.
     * @param M row dimension of \p B
     * @param N col dimension of \p B
     * @param R rank of \p A
     * @param A input matrix
     * @param lda leading dimension of \p A
     * @param P row permutation of the \c PLUQ decomposition of \p A
     * @param Q column permutation of the \c PLUQ decomposition of \p A
     * @param B Right/Left hand side matrix. Initially stores \p B, finally stores the solution \p X.
     * @param ldb leading dimension of \p B
     * @param info Success of the computation: 0 if successfull, >0 if system is inconsistent
     */
    template <class Field>
    void
    fgetrs (const Field& F,
            const FFLAS::FFLAS_SIDE Side,
            const size_t M, const size_t N, const size_t R,
            typename Field::Element_ptr A, const size_t lda,
            const size_t *P, const size_t *Q,
            typename Field::Element_ptr B, const size_t ldb,
            int * info);

    /** Solve the system A X = B or X A = B.
     * Solving using the PLUQ decomposition of A
     * already computed inplace with PLUQ(FFLAS::FflasNonUnit).
     * Version for A rectangular.
     * If A is rank deficient, a solution is returned if the system is consistent,
     * Otherwise an info is 1
     *
     * @param F base field
     * @param Side Determine wheter the resolution is left (FflasLeft) or right (FflasRight) looking.
     * @param M row dimension of A
     * @param N col dimension of A
     * @param NRHS number of columns (if Side = FFLAS::FflasLeft) or row (if Side = FFLAS::FflasRight) of the matrices X and B
     * @param R rank of A
     * @param A input matrix
     * @param lda leading dimension of A
     * @param P row permutation of the PLUQ decomposition of A
     * @param Q column permutation of the PLUQ decomposition of A
     * @param X solution matrix
     * @param ldx leading dimension of X
     * @param B Right/Left hand side matrix.
     * @param ldb leading dimension of B
     * @param info Succes of the computation: 0 if successfull, >0 if system is inconsistent
     */
    template <class Field>
    typename Field::Element_ptr
    fgetrs (const Field& F,
            const FFLAS::FFLAS_SIDE Side,
            const size_t M, const size_t N, const size_t NRHS, const size_t R,
            typename Field::Element_ptr A, const size_t lda,
            const size_t *P, const size_t *Q,
            typename Field::Element_ptr X, const size_t ldx,
            typename Field::ConstElement_ptr B, const size_t ldb,
            int * info);

    /** @brief Square system solver
     * @param F The computation domain
     * @param Side Determine wheter the resolution is left (FflasLeft) or right (FflasRight) looking
     * @param M row dimension of B
     * @param N col dimension of B
     * @param A input matrix
     * @param lda leading dimension of A
     * @param B Right/Left hand side matrix. Initially contains B, finally contains the solution X.
     * @param ldb leading dimension of B
     * @param info Success of the computation: 0 if successfull, >0 if system is inconsistent
     * @return the rank of the system
     *
     * Solve the system A X = B or X A = B.
     * Version for A square.
     * If A is rank deficient, a solution is returned if the system is consistent,
     * Otherwise an info is 1
     */
    template <class Field>
    size_t
    fgesv (const Field& F,
           const FFLAS::FFLAS_SIDE Side,
           const size_t M, const size_t N,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::Element_ptr B, const size_t ldb,
           int * info);

    /**  @brief Rectangular system solver
     * @param F The computation domain
     * @param Side Determine wheter the resolution is left (FflasLeft) or right (FflasRight) looking
     * @param M row dimension of A
     * @param N col dimension of A
     * @param NRHS number of columns (if Side = FFLAS::FflasLeft) or row (if Side = FFLAS::FflasRight) of the matrices X and B
     * @param A input matrix
     * @param lda leading dimension of A
     * @param B Right/Left hand side matrix. Initially contains B, finally contains the solution X.
     * @param ldb leading dimension of B
     * @param X
     * @param ldx
     * @param info Success of the computation: 0 if successfull, >0 if system is inconsistent
     * @return the rank of the system
     *
     * Solve the system A X = B or X A = B.
     * Version for A square.
     * If A is rank deficient, a solution is returned if the system is consistent,
     * Otherwise an info is 1
     */
    template <class Field>
    size_t
    fgesv (const Field& F,
           const FFLAS::FFLAS_SIDE Side,
           const size_t M, const size_t N, const size_t NRHS,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::Element_ptr X, const size_t ldx,
           typename Field::ConstElement_ptr B, const size_t ldb,
           int * info);
} // FFPACK fgesv, fgetrs
// #include "ffpack_fgesv.inl"
// #include "ffpack_fgetrs.inl"

namespace FFPACK { /* ftrtr */

    /** Compute the inverse of a triangular matrix.
     * @param F base field
     * @param Uplo whether the matrix is upper or lower triangular
     * @param Diag whether the matrix is unit diagonal (FflasUnit/NoUnit)
     * @param N input matrix order
     * @param A the input matrix
     * @param lda leading dimension of A
     *
     */
    template<class Field>
    void
    ftrtri (const Field& F, const FFLAS::FFLAS_UPLO Uplo, const FFLAS::FFLAS_DIAG Diag,
            const size_t N, typename Field::Element_ptr A, const size_t lda,
            const size_t threshold = __FFLASFFPACK_FTRTRI_THRESHOLD);


    template<class Field>
    void trinv_left( const Field& F, const size_t N, typename Field::ConstElement_ptr L, const size_t ldl,
                     typename Field::Element_ptr X, const size_t ldx );

    /**  @brief Compute the product of two triangular matrices of opposite shape.
     * Product UL or LU of the upper, resp lower triangular matrices U and L
     * stored one above the other in the square matrix A.
     * @param F base field
     * @param Side set to FflasLeft to compute the product UL, FflasRight to compute LU
     * @param diag whether the matrix U is unit diagonal (FflasUnit/NoUnit)
     * @param N input matrix order
     * @param A the input matrix
     * @param lda leading dimension of A
     *
     */
    template<class Field>
    void
    ftrtrm (const Field& F, const FFLAS::FFLAS_SIDE side, const FFLAS::FFLAS_DIAG diag,
            const size_t N,	typename Field::Element_ptr A, const size_t lda);

    /** @brief Solve a triangular system with a triangular right hand side of the same shape.
     * @param F base field
     * @param Side set to FflasLeft to compute U1^-1*U2 or L1^-1*L2, FflasRight to compute U1*U2^-1 or L1*L2^-1
     * @param Uplo whether the matrix A is upper or lower triangular
     * @param diag1 whether the matrix U1 or L2 is unit diagonal (FflasUnit/NoUnit)
     * @param diag2 whether the matrix U2 or L2 is unit diagonal (FflasUnit/NoUnit)
     * @param N order of the input matrices
     * @param A the input matrix to be inverted (U1 or L1)
     * @param lda leading dimension of A
     * @param B the input right hand side (U2 or L2)
     * @param ldb leading dimension of B
     */
    template<class Field>
    void
    ftrstr (const Field& F, const FFLAS::FFLAS_SIDE side, const FFLAS::FFLAS_UPLO Uplo,
            const FFLAS::FFLAS_DIAG diagA, const FFLAS::FFLAS_DIAG diagB, const size_t N,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::Element_ptr B, const size_t ldb, const size_t threshold=__FFLASFFPACK_FTRSTR_THRESHOLD);

    /** @brief Solve a triangular system in a symmetric sum: find B upper/lower triangular such that A^T B + B^T A = C
     * where C is symmetric. C is overwritten by B.
     * @param F base field
     * @param Side set to FflasLeft to compute U1^-1*U2 or L1^-1*L2, FflasRight to compute U1*U2^-1 or L1*L2^-1
     * @param Uplo whether the matrix A is upper or lower triangular
     * @param diagA whether the matrix A is unit diagonal (FflasUnit/NoUnit)
     * @param N order of the input matrices
     * @param A the input matrix
     * @param lda leading dimension of A
     * @param [inout] B the input right hand side where the output is written
     * @param ldb leading dimension of B
     */
    template<class Field>
    void
    ftrssyr2k (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
               const FFLAS::FFLAS_DIAG diagA, const size_t N,
               typename Field::ConstElement_ptr A, const size_t lda,
               typename Field::Element_ptr B, const size_t ldb, const size_t threshold=__FFLASFFPACK_FTRSSYR2K_THRESHOLD);

} // FFPACK ftrtr
// #include "ffpack_ftrtr.inl"

namespace FFPACK {

    /* LDLT or UTDU factorizations */

    /** @brief Triangular factorization of symmetric matrices
     * @param F The computation domain
     * @param UpLo Determine wheter to store the upper (FflasUpper) or lower (FflasLower) triangular factor
     * @param N order of the matrix A
     * @param [inout] A input matrix
     * @param lda leading dimension of A
     * @return false if the \p A does not have generic rank profile, making the computation fail.
     *
     * Compute the a triangular factorization of the matrix A: \f$ A = L \times D \times  L^T\f$ if UpLo = FflasLower or
     * \f$ A = U^T \times D \times  U\f$ otherwise. \p D is a diagonal matrix. The matrices \p L and \p U are unit
     * diagonal lower (resp. upper) triangular and overwrite the input matrix \p A.
     * The matrix \p D is stored on the diagonal of \p A, as the diagonal of \p L or \p U is known to be all ones.
     * If A does not have generic rank profile, the LDLT or UTDU factorizations is not defined, and the algorithm
     * returns false.
     */
    template <class Field>
    bool fsytrf (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 const size_t threshold = __FFLASFFPACK_FSYTRF_THRESHOLD);

    template <class Field>
    bool fsytrf (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 const FFLAS::ParSeqHelper::Sequential seq,
                 const size_t threshold = __FFLASFFPACK_FSYTRF_THRESHOLD);

    template <class Field, class Cut, class Param>
    bool fsytrf (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 const FFLAS::ParSeqHelper::Parallel<Cut,Param> par,
                 const size_t threshold = __FFLASFFPACK_FSYTRF_THRESHOLD);

    /* LDLT or UTDU factorizations */

    /** @brief Triangular factorization of symmetric matrices
     * @param F The computation domain
     * @param UpLo Determine wheter to store the upper (FflasUpper) or lower (FflasLower) triangular factor
     * @param N order of the matrix A
     * @param [inout] A input matrix
     * @param [inout] D
     * @param lda leading dimension of A
     * @return false if the \p A does not have generic rank profile, making the computation fail.
     *
     * Compute the a triangular factorization of the matrix A: \f$ A = L \times Dinv \times  L^T\f$ if UpLo = FflasLower
     * or \f$ A = U^T \times D \times  U\f$ otherwise. \p D is a diagonal matrix. The matrices \p L and \p U are
     * lower (resp. upper) triangular and overwrite the input matrix \p A.
     * The matrix \p D need to be stored separately, as the diagonal of \p L or \p U are not unit.
     * If A does not have generic rank profile, the LDLT or UTDU factorizations is not defined, and the algorithm
     * returns false.
     */
    template <class Field>
    bool fsytrf_nonunit (const Field& F, const FFLAS::FFLAS_UPLO UpLo, const size_t N,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr D, const size_t incD,
                         const size_t threshold = __FFLASFFPACK_FSYTRF_THRESHOLD);
    /* PLUQ */

    /** @brief Compute a PLUQ factorization of the given matrix.
     * Return its rank.
     * The permutations P and Q are represented
     * using LAPACK's convention.
     * @param F base field
     * @param Diag   whether U should have a unit diagonal (FflasUnit) or not (FflasNoUnit)
     * @param M matrix row dimension
     * @param N matrix column dimension
     * @param A input matrix
     * @param lda leading dimension of \p A
     * @param P the row permutation
     * @param Q the column permutation

     * @return the rank of \p A
     * @bib
     * - Dumas J-G.,  Pernet C., and Sultan Z. <i>\c Simultaneous computation of the row and column rank profiles </i>, ISSAC'13, 2013
     * .
     */
    template<class Field>
    size_t PLUQ (const Field& F, const FFLAS::FFLAS_DIAG Diag,
                 const size_t M, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 size_t*P, size_t *Q);

    template<class Field>
    size_t pPLUQ (const Field& F, const FFLAS::FFLAS_DIAG Diag,
                 const size_t M, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 size_t*P, size_t *Q);

    template<class Field>
    size_t PLUQ (const Field& F, const FFLAS::FFLAS_DIAG Diag,
                 const size_t M, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 size_t*P, size_t *Q, const FFLAS::ParSeqHelper::Sequential& PSHelper,
                 size_t BCThreshold = __FFLASFFPACK_PLUQ_THRESHOLD);

    template<class Field, class Cut, class Param>
    size_t PLUQ (const Field& F, const FFLAS::FFLAS_DIAG Diag,
                 const size_t M, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 size_t*P, size_t *Q, const FFLAS::ParSeqHelper::Parallel<Cut,Param>& PSHelper);

} // FFPACK PLUQ
// #include "ffpack_pluq.inl"

namespace FFPACK { /* ludivine */

    /** @brief Compute the CUP or PLE factorization of the given matrix.
     * Using
     * a block algorithm and return its rank.
     * The permutations P and Q are represented
     * using LAPACK's convention.
     * @param F base field
     * @param Diag  whether the transformation matrix (U of the CUP, L of the PLE) should have a unit diagonal (FflasUnit)
     * or not (FflasNoUnit)
     * @param trans whether to compute the CUP decomposition (FflasNoTrans) or the PLE decomposition (FflasTrans)
     * @param M matrix row dimension
     * @param N matrix column dimension
     * @param A input matrix
     * @param lda leading dimension of \p A
     * @param P the factor of CUP or PLE
     * @param Q a permutation indicating the pivot position in the echelon form C or E in its first r positions
     * @param LuTag flag for setting the earling termination if the matrix
     * is singular
     * @param cutoff threshold to basecase
     *
     * @return the rank of \p A
     * @bib
     * - Jeannerod C-P, Pernet, C. and Storjohann, A. <i>\c Rank-profile revealing Gaussian elimination and the CUP matrix decomposition  </i>, J. of Symbolic Comp., 2013
     * - Pernet C, Brassel M <i>\c LUdivine, une divine factorisation \c LU</i>, 2002
     * .
     */
    template <class Field>
    size_t
    LUdivine (const Field& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
              const size_t M, const size_t N,
              typename Field::Element_ptr A, const size_t lda,
              size_t* P, size_t* Qt,
              const FFPACK_LU_TAG LuTag = FfpackSlabRecursive,
              const size_t cutoff=__FFLASFFPACK_LUDIVINE_THRESHOLD);

    /* \cond */
    template<class Element>
    class callLUdivine_small;

    //! LUdivine small case
    template <class Field>
    size_t
    LUdivine_small (const Field& F, const FFLAS::FFLAS_DIAG Diag,  const FFLAS::FFLAS_TRANSPOSE trans,
                    const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda,
                    size_t* P, size_t* Q,
                    const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    //! LUdivine gauss
    template <class Field>
    size_t
    LUdivine_gauss (const Field& F, const FFLAS::FFLAS_DIAG Diag,
                    const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda,
                    size_t* P, size_t* Q,
                    const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    /* \endcond */
    namespace Protected {



        //---------------------------------------------------------------------
        // LUdivine_construct: (Specialisation of LUdivine)
        // LUP factorisation of X, the Krylov base matrix of A^t and v, in A.
        // X contains the nRowX first vectors v, vA, .., vA^{nRowX-1}
        // A contains the LUP factorisation of the nUsedRowX first row of X.
        // When all rows of X have been factorized in A, and rank is full,
        // then X is updated by the following scheme: X <= ( X; X.B ), where
        // B = A^2^i.
        // This enables to make use of Matrix multiplication, and stop computing
        // Krylov vector, when the rank is not longer full.
        // P is the permutation matrix stored in an array of indexes
        //---------------------------------------------------------------------

        template <class Field>
        size_t
        LUdivine_construct( const Field& F, const FFLAS::FFLAS_DIAG Diag,
                            const size_t M, const size_t N,
                            typename Field::ConstElement_ptr A, const size_t lda,
                            typename Field::Element_ptr X, const size_t ldx,
                            typename Field::Element_ptr u, const size_t incu, size_t* P,
                            bool computeX, const FFPACK_MINPOLY_TAG MinTag= FfpackDense
                            , const size_t kg_mc =0
                            , const size_t kg_mb =0
                            , const size_t kg_j  =0);

    } // Protected

} //FFPACK ludivine, turbo
// #include "ffpack_ludivine.inl"

namespace FFPACK { /* echelon */
    /*****************/
    /* ECHELON FORMS */
    /*****************/

    /** Compute the Column Echelon form of the input matrix in-place.
     *
     * If LuTag == FfpackTileRecursive, then after the computation A = [ M \ V ]
     * such that AU = C is a column echelon decomposition of A,
     * with U = P^T [   V    ] and C = M + Q [ Ir ]
     *              [ 0 In-r ]               [ 0  ]
     * If LuTag == FfpackTileRecursive then A = [ N \ V ] such that the same holds with M = Q N
     *
     * Qt = Q^T
     * If transform=false, the matrix V is not computed.
     * See also test-colechelon for an example of use
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param[in] A input matrix
     * @param lda leading dimension of A
     * @param P the column permutation
     * @param Qt the row position of the pivots in the echelon form
     * @param transform decides whether V is computed
     * @param LuTag chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
     */
    template <class Field>
    size_t
    ColumnEchelonForm (const Field& F, const size_t M, const size_t N,
                       typename Field::Element_ptr A, const size_t lda,
                       size_t* P, size_t* Qt, bool transform=false,
                       const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    template <class Field>
    size_t
    pColumnEchelonForm (const Field& F, const size_t M, const size_t N,
                        typename Field::Element_ptr A, const size_t lda,
                        size_t* P, size_t* Qt, bool transform=false,
                        size_t numthreads = 0, const FFPACK_LU_TAG LuTag=FfpackTileRecursive);

    template <class Field, class PSHelper>
    size_t
    ColumnEchelonForm (const Field& F, const size_t M, const size_t N,
                              typename Field::Element_ptr A, const size_t lda,
                              size_t* P, size_t* Qt, const bool transform,
                              const FFPACK_LU_TAG LuTag, const PSHelper& psH);



    /**  Compute the Row Echelon form of the input matrix in-place.
     *
     * If LuTag == FfpackTileRecursive, then after the computation A = [ L \ M ]
     * such that X A = R is a row echelon decomposition of A,
     * with X =  [ L  0   ] P  and R = M + [Ir 0] Q^T
     *           [    In-r]
     * If LuTag == FfpackTileRecursive then A = [ L \ N ] such that the same holds with M =  N Q^T
     * Qt = Q^T
     * If transform=false, the matrix L is not computed.
     * See also test-rowechelon for an example of use
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param[in] A the input matrix
     * @param lda leading dimension of A
     * @param P the row permutation
     * @param Qt the column position of the pivots in the echelon form
     * @param transform decides whether L is computed
     * @param LuTag chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
     */
    template <class Field>
    size_t
    RowEchelonForm (const Field& F, const size_t M, const size_t N,
                    typename Field::Element_ptr A, const size_t lda,
                    size_t* P, size_t* Qt, const bool transform=false,
                    const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    template <class Field>
    size_t
    pRowEchelonForm (const Field& F, const size_t M, const size_t N,
                     typename Field::Element_ptr A, const size_t lda,
                     size_t* P, size_t* Qt, const bool transform=false,
                     size_t numthreads = 0, const FFPACK_LU_TAG LuTag=FfpackTileRecursive);

    template <class Field, class PSHelper>
    size_t
    RowEchelonForm (const Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           size_t* P, size_t* Qt, const bool transform,
                           const FFPACK_LU_TAG LuTag, const PSHelper& psH);


    /** Compute the Reduced Column Echelon form of the input matrix in-place.
     *
     * After the computation A = [ V   ] such that AX = R is a reduced col echelon
     *                           [ M 0 ]
     * decomposition of A, where X = P^T [ V      ] and R = Q [ Ir   ]
     *                                   [ 0 In-r ]           [ M  0 ]
     * Qt = Q^T
     * If transform=false, the matrix X is not computed and the matrix A = R
     *
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param[in] A input matrix
     * @param lda leading dimension of A
     * @param P the column permutation
     * @param Qt the row position of the pivots in the echelon form
     * @param transform decides whether X is computed
     * @param LuTag chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
     */
    template <class Field>
    size_t
    ReducedColumnEchelonForm (const Field& F, const size_t M, const size_t N,
                              typename Field::Element_ptr A, const size_t lda,
                              size_t* P, size_t* Qt, const bool transform = false,
                              const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    template <class Field>
    size_t
    pReducedColumnEchelonForm (const Field& F, const size_t M, const size_t N,
                               typename Field::Element_ptr A, const size_t lda,
                               size_t* P, size_t* Qt, const bool transform = false,
                               size_t numthreads = 0, const FFPACK_LU_TAG LuTag=FfpackTileRecursive);

    template <class Field, class PSHelper>
    size_t
    ReducedColumnEchelonForm (const Field& F, const size_t M, const size_t N,
                              typename Field::Element_ptr A, const size_t lda,
                              size_t* P, size_t* Qt, const bool transform,
                              const FFPACK_LU_TAG LuTag, const PSHelper& psH);



    /** Compute the Reduced Row Echelon form of the input matrix in-place.
     *
     * After the computation A = [ V1 M ] such that X A = R is a reduced row echelon
     *                           [ V2 0 ]
     * decomposition of A, where X =  [ V1  0   ] P and R =  [ Ir M  ] Q^T
     *                                [ V2 In-r ]            [ 0     ]
     * Qt = Q^T
     * If transform=false, the matrix X is not computed and the matrix A = R
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param[in] A input matrix
     * @param lda leading dimension of A
     * @param P the row permutation
     * @param Qt the column position of the pivots in the echelon form
     * @param transform decides whether X is computed
     * @param LuTag chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
     */
    template <class Field>
    size_t
    ReducedRowEchelonForm (const Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           size_t* P, size_t* Qt, const bool transform = false,
                           const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    template <class Field>
    size_t
    pReducedRowEchelonForm (const Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           size_t* P, size_t* Qt, const bool transform = false,
                           size_t numthreads = 0, const FFPACK_LU_TAG LuTag=FfpackTileRecursive);

    template <class Field, class PSHelper>
    size_t
    ReducedRowEchelonForm (const Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           size_t* P, size_t* Qt, const bool transform,
                           const FFPACK_LU_TAG LuTag, const PSHelper& psH);




    namespace Protected {
        /**  @brief Gauss-Jordan algorithm computing the Reduced Row echelon form and its transform matrix.
         * @bib
         *  - Algorithm 2.8 of A. Storjohann Thesis 2000,
         *  - Algorithm 11 of Jeannerod C-P., Pernet, C. and Storjohann, A. <i>\c Rank-profile revealing Gaussian elimination and the CUP matrix decomposition  </i>, J. of Symbolic Comp., 2013
         * @param M row dimension of A
         * @param N column dimension of A
         * @param [inout] A an m x n matrix
         * @param lda leading dimension of A
         * @param P row permutation
         * @param Q column permutation
         * @param LuTag set the base case to a Tile (FfpackGaussJordanTile)  or Slab (FfpackGaussJordanSlab) recursive RedEchelon
         */
        template <class Field>
        size_t
        GaussJordan (const Field& F, const size_t M, const size_t N,
                     typename Field::Element_ptr A, const size_t lda,
                     const size_t colbeg, const size_t rowbeg, const size_t colsize,
                     size_t* P, size_t* Q, const FFPACK::FFPACK_LU_TAG LuTag);
    } // Protected
} // FFPACK
// #include "ffpack_echelonforms.inl"

namespace FFPACK { /* invert */
    /*****************/
    /*   INVERSION   */
    /*****************/
    /**  @brief Invert the given matrix in place
     * or computes its nullity if it is singular.
     *
     * An inplace \f$2n^3\f$ algorithm is used.
     * @param F The computation domain
     * @param M order of the matrix
     * @param [in,out] A input matrix (\f$M \times M\f$)
     * @param lda leading dimension of A
     * @param nullity dimension of the kernel of A
     * @return pointer to \f$A\f$ and \f$A \gets A^{-1}\f$
     */
    template <class Field>
    typename Field::Element_ptr
    Invert (const Field& F, const size_t M,
            typename Field::Element_ptr A, const size_t lda,
            int& nullity);

    /** @brief Invert the given matrix
     * or computes its nullity if it is singular.
     *
     * @pre \p X is preallocated and should be large enough to store the
     * \f$ m \times m\f$ matrix \p A.
     *
     * @param F The computation domain
     * @param M order of the matrix
     * @param [in] A input matrix (\f$M \times M\f$)
     * @param lda leading dimension of \p A
     * @param [out] X this is the inverse of \p A if \p A is invertible
     * (non \c NULL and \f$ \mathtt{nullity} = 0\f$). It is untouched
     * otherwise.
     * @param ldx leading dimension of \p X
     * @param nullity dimension of the kernel of \p A
     * @return pointer to \f$X = A^{-1}\f$
     */
    template <class Field>
    typename Field::Element_ptr
    Invert (const Field& F, const size_t M,
            typename Field::ConstElement_ptr A, const size_t lda,
            typename Field::Element_ptr X, const size_t ldx,
            int& nullity);

    /** @brief Invert the given matrix or computes its nullity if it is singular.
     *
     * An \f$2n^3f\f$ algorithm is used.
     * This routine can be \% faster than FFPACK::Invert but is not totally inplace.
     *
     * @pre \p X is preallocated and should be large enough to store the
     * \f$ m \times m\f$ matrix \p A.
     *
     * @warning A is overwritten here !
     * @bug not tested.
     * @param F the computation domain
     * @param M order of the matrix
     * @param [in,out] A input matrix (\f$M \times M\f$). On output, \p A
     * is modified and represents a "psycological" factorisation \c LU.
     * @param lda leading dimension of A
     * @param [out] X this is the inverse of \p A if \p A is invertible
     * (non \c NULL and \f$ \mathtt{nullity} = 0\f$). It is untouched
     * otherwise.
     * @param ldx leading dimension of \p X
     * @param nullity dimension of the kernel of \p A
     * @return pointer to \f$X = A^{-1}\f$
     */
    template <class Field>
    typename Field::Element_ptr
    Invert2( const Field& F, const size_t M,
             typename Field::Element_ptr A, const size_t lda,
             typename Field::Element_ptr X, const size_t ldx,
             int& nullity);

} // FFPACK invert
// #include "ffpack_invert.inl"

namespace FFPACK { /* charpoly */
    /*****************************/
    /* CHARACTERISTIC POLYNOMIAL */
    /*****************************/


    /**
     * @brief Compute the characteristic polynomial of the matrix A.
     * @param R the polynomial ring of charp (contains the base field)
     * @param [out] charp the characteristic polynomial of \p as a list of factors
     * @param N order of the matrix \p A
     * @param [in] A the input matrix (\f$ N \times N\f$) (could be overwritten in some algorithmic variants)
     * @param lda leading dimension of \p A
     * @param CharpTag the algorithmic variant
     * @param G a random iterator (required for the randomized variants LUKrylov and ArithProg)
     */
    template <class PolRing>
    inline std::list<typename PolRing::Element>&
    CharPoly (const PolRing& R, std::list<typename PolRing::Element>& charp, const size_t N,
              typename PolRing::Domain_t::Element_ptr A, const size_t lda,
              typename PolRing::Domain_t::RandIter& G,
              const FFPACK_CHARPOLY_TAG CharpTag= FfpackAuto,
              const size_t degree = __FFLASFFPACK_ARITHPROG_THRESHOLD);

    /**
     * @brief Compute the characteristic polynomial of the matrix A.
     * @param R the polynomial ring of charp (contains the base field)
     * @param [out] charp the characteristic polynomial of \p as a single polynomial
     * @param N order of the matrix \p A
     * @param [in] A the input matrix (\f$ N \times N\f$) (could be overwritten in some algorithmic variants)
     * @param lda leading dimension of \p A
     * @param CharpTag the algorithmic variant
     * @param G a random iterator (required for the randomized variants LUKrylov and ArithProg)
     */
    template <class PolRing>
    inline typename PolRing::Element&
    CharPoly (const PolRing& R, typename PolRing::Element& charp, const size_t N,
              typename PolRing::Domain_t::Element_ptr A, const size_t lda,
              typename PolRing::Domain_t::RandIter& G,
              const FFPACK_CHARPOLY_TAG CharpTag= FfpackAuto,
              const size_t degree = __FFLASFFPACK_ARITHPROG_THRESHOLD);

    /**
     * @brief Compute the characteristic polynomial of the matrix A.
     * @param R the polynomial ring of charp (contains the base field)
     * @param [out] charp the characteristic polynomial of \p as a single polynomial
     * @param N order of the matrix \p A
     * @param [in] A the input matrix (\f$ N \times N\f$) (could be overwritten in some algorithmic variants)
     * @param lda leading dimension of \p A
     * @param CharpTag the algorithmic variant
     */
    template <class PolRing>
    inline typename PolRing::Element&
    CharPoly (const PolRing& R, typename PolRing::Element& charp, const size_t N,
              typename PolRing::Domain_t::Element_ptr A, const size_t lda,
              const FFPACK_CHARPOLY_TAG CharpTag= FfpackAuto,
              const size_t degree = __FFLASFFPACK_ARITHPROG_THRESHOLD){
        typename PolRing::Domain_t::RandIter G(R.getdomain());
        return CharPoly (R, charp, N, A, lda, G, CharpTag, degree);
    }


    namespace Protected {
        template <class Field, class Polynomial>
        std::list<Polynomial>&
        KellerGehrig( const Field& F, std::list<Polynomial>& charp, const size_t N,
                      typename Field::ConstElement_ptr A, const size_t lda );

        template <class Field, class Polynomial>
        int
        KGFast ( const Field& F, std::list<Polynomial>& charp, const size_t N,
                 typename Field::Element_ptr A, const size_t lda,
                 size_t * kg_mc, size_t* kg_mb, size_t* kg_j );

        template <class Field, class Polynomial>
        std::list<Polynomial>&
        KGFast_generalized (const Field& F, std::list<Polynomial>& charp,
                            const size_t N,
                            typename Field::Element_ptr A, const size_t lda);


        template<class Field>
        void
        fgemv_kgf( const Field& F,  const size_t N,
                   typename Field::ConstElement_ptr A, const size_t lda,
                   typename Field::ConstElement_ptr X, const size_t incX,
                   typename Field::Element_ptr Y, const size_t incY,
                   const size_t kg_mc, const size_t kg_mb, const size_t kg_j );

        template <class Field, class Polynomial, class RandIter>
        std::list<Polynomial>&
        LUKrylov( const Field& F, std::list<Polynomial>& charp, const size_t N,
                  typename Field::Element_ptr A, const size_t lda,
                  typename Field::Element_ptr U, const size_t ldu, RandIter& G);

        template <class Field, class Polynomial>
        std::list<Polynomial>&
        Danilevski (const Field& F, std::list<Polynomial>& charp,
                    const size_t N, typename Field::Element_ptr A, const size_t lda);


        template <class PolRing>
        inline void
        RandomKrylovPrecond (const PolRing& PR, std::list<typename PolRing::Element>& completedFactors, const size_t N,
                             typename PolRing::Domain_t::Element_ptr A, const size_t lda,
                             size_t& Nb, typename PolRing::Domain_t::Element_ptr& B, size_t& ldb,
                             typename PolRing::Domain_t::RandIter& g, const size_t degree=__FFLASFFPACK_ARITHPROG_THRESHOLD);
        
        template <class PolRing>
        inline std::list<typename PolRing::Element>&
        ArithProg (const PolRing& PR, std::list<typename PolRing::Element>& frobeniusForm,
                   const size_t N, typename PolRing::Domain_t::Element_ptr A, const size_t lda,
                   const size_t degree);

        template <class Field, class Polynomial>
        std::list<Polynomial>&
        LUKrylov_KGFast( const Field& F, std::list<Polynomial>& charp, const size_t N,
                         typename Field::Element_ptr A, const size_t lda,
                         typename Field::Element_ptr X, const size_t ldx);
    } // Protected
} // FFPACK charpoly
// #include "ffpack_charpoly_kglu.inl"
// #include "ffpack_charpoly_kgfast.inl"
// #include "ffpack_charpoly_kgfastgeneralized.inl"
// #include "ffpack_charpoly_danilevski.inl"
// #include "ffpack_charpoly.inl"

namespace FFPACK { /* frobenius, charpoly */


} // FFPACK frobenius
// #include "ffpack_frobenius.inl"

namespace FFPACK { /* minpoly */


    /**********************/
    /* MINIMAL POLYNOMIAL */
    /**********************/

    /**
     * @brief Compute the minimal polynomial of the matrix A.
     * The algorithm is randomized probabilistic, and computes the minimal polynomial of
     * the Krylov iterates of a random vector: (v, Av, .., A^kv)
     * @param F the base field
     * @param [out] minP the minimal polynomial of \p A
     * @param N order of the matrix \p A
     * @param [in] A the input matrix (\f$ N \times N\f$)
     * @param lda leading dimension of \p A
     */
    template <class Field, class Polynomial>
    Polynomial&
    MinPoly (const Field& F, Polynomial& minP, const size_t N,
             typename Field::ConstElement_ptr A, const size_t lda);

    /**
     * @brief Compute the minimal polynomial of the matrix A.
     * The algorithm is randomized probabilistic, and computes the minimal polynomial of
     * the Krylov iterates of a random vector: (v, Av, .., A^kv)
     * @param F the base field
     * @param [out] minP the minimal polynomial of \p A
     * @param N order of the matrix \p A
     * @param [in] A the input matrix (\f$ N \times N\f$)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     */
    template <class Field, class Polynomial, class RandIter>
    Polynomial&
    MinPoly (const Field& F, Polynomial& minP, const size_t N,
             typename Field::ConstElement_ptr A, const size_t lda, RandIter& G);

    /**
     * @brief Compute the minimal polynomial of the matrix A and a vector v, namely the first linear dependency relation in the Krylov basis \f$(v,Av, ..., A^Nv)\f$.
     * @param F the base field
     * @param [out] minP the minimal polynomial of \p A and v
     * @param N order of the matrix \p A
     * @param [in] A the input matrix (\f$ N \times N\f$)
     * @param lda leading dimension of \p A
     * @param K an \f$ N \times (N+1)\f$ matrix containing the vector v on its first row
     * @param ldk leading dimension of \p K
     * @param P [out] (optional) the permutation used in the elimination of the Krylov matrix \p K
     */
    template <class Field, class Polynomial>
    Polynomial&
    MatVecMinPoly (const Field& F, Polynomial& minP, const size_t N,
                   typename Field::ConstElement_ptr A, const size_t lda,
                   typename Field::ConstElement_ptr v, const size_t incv);

    namespace Protected{
        template <class Field, class Polynomial>
        Polynomial&
        MatVecMinPoly (const Field& F, Polynomial& minP, const size_t N,
                       typename Field::ConstElement_ptr A, const size_t lda,
                       typename Field::Element_ptr v, const size_t incv,
                       typename Field::Element_ptr K, const size_t ldk,
                       size_t * P);

        template <class Field, class Polynomial>
        Polynomial&
        Hybrid_KGF_LUK_MinPoly (const Field& F, Polynomial& minP, const size_t N,
                                typename Field::ConstElement_ptr A, const size_t lda,
                                typename Field::Element_ptr X, const size_t ldx, size_t* P,
                                const FFPACK_MINPOLY_TAG MinTag= FFPACK::FfpackDense,
                                const size_t kg_mc=0, const size_t kg_mb=0, const size_t kg_j=0);
    } // Protected
} // FFPACK minpoly
// #include "ffpack_minpoly.inl"

namespace FFPACK { /* Krylov Elim */

    /* \cond */
    template <class Field>
    size_t KrylovElim( const Field& F, const size_t M, const size_t N,
                       typename Field::Element_ptr A, const size_t lda, size_t*P,
                       size_t *Q, const size_t deg, size_t *iterates, size_t * inviterates, const size_t maxit,size_t virt);

    template <class Field>
    size_t  SpecRankProfile (const Field& F, const size_t M, const size_t N,
                             typename Field::Element_ptr A, const size_t lda, const size_t deg, size_t *rankProfile);
    /* \endcond */

} // FFPACK KrylovElim
// #include "ffpack_krylovelim.inl"

namespace FFPACK { /* Solutions */
    /********/
    /* RANK */
    /********/



    /** Computes the rank of the given matrix using a PLUQ factorization.
     * The input matrix is modified.
     * @param F base field
     * @param M row dimension of the matrix
     * @param N column dimension of the matrix
     * @param [in] A input matrix
     * @param lda leading dimension of A
     * @param psH (optional) a ParSeqHelper to choose between sequential and parallel execution
     */

    template <class Field>
    size_t
    Rank( const Field& F, const size_t M, const size_t N,
          typename Field::Element_ptr A, const size_t lda);

    template <class Field>
    size_t
    pRank (const Field& F, const size_t M, const size_t N,
           typename Field::Element_ptr A, const size_t lda, size_t numthreads = 0);

    template <class Field, class PSHelper>
    size_t
    Rank( const Field& F, const size_t M, const size_t N,
          typename Field::Element_ptr A, const size_t lda, const PSHelper& psH) ;


    /********/
    /* DET  */
    /********/


    /**  Returns true if the given matrix is singular.
     * The method is a block elimination with early termination
     *
     * using LQUP factorization  with early termination.
     * If <code>M != N</code>,
     * then the matrix is virtually padded with zeros to make it square and
     * it's determinant is zero.
     * @warning The input matrix is modified.
     * @param F base field
     * @param M row dimension of the matrix
     * @param N column dimension of the matrix.
     * @param [in,out] A input matrix
     * @param lda leading dimension of A
     */
    template <class Field>
    bool
    IsSingular( const Field& F, const size_t M, const size_t N,
                typename Field::Element_ptr A, const size_t lda);

    /** @brief Returns the determinant of the given square matrix.
     * @details The method is a block elimination
     * using PLUQ factorization. The input matrix A is overwritten.
     * @warning The input matrix is modified.
     * @param F base field
     * @param [out] det the determinant of A
     * @param N the order of the square matrix A.
     * @param [in,out] A input matrix
     * @param lda leading dimension of A
     * @param psH (optional) a ParSeqHelper to choose between sequential and parallel execution
     * @param P,Q (optional) row and column permutations to be used by the PLUQ factorization. randomized checkers (see cherckes/checker_det.inl) need them for certification
     */

    template <class Field>
    typename Field::Element&
    Det (const Field& F, typename Field::Element& det, const size_t N,
         typename Field::Element_ptr A, const size_t lda,
         size_t * P = NULL, size_t * Q = NULL);

    template <class Field>
    typename Field::Element&
    pDet (const Field& F, typename Field::Element& det, const size_t N,
          typename Field::Element_ptr A, const size_t lda,
          size_t numthreads = 0, size_t * P = NULL, size_t * Q = NULL);

    template <class Field, class PSHelper>
    typename Field::Element&
    Det(const Field& F, typename Field::Element& det, const size_t N,
        typename Field::Element_ptr A, const size_t lda, const PSHelper& psH,
        size_t * P = NULL, size_t * Q = NULL);

    /*********/
    /* SOLVE */
    /*********/


    /**
     * @brief Solves a linear system AX = b using PLUQ factorization.
     * @oaram F base field
     * @oaram M matrix order
     * @param [in] A input matrix
     * @param lda leading dimension of A
     * @param [out] x output solution vector
     * @param incx increment of x
     * @param b input right hand side of the system
     * @param incb increment of b
     */
    template <class Field>
    typename Field::Element_ptr
    Solve( const Field& F, const size_t M,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::Element_ptr x, const int incx,
           typename Field::ConstElement_ptr b, const int incb );

    template <class Field, class PSHelper>
    typename Field::Element_ptr
    Solve( const Field& F, const size_t M,
           typename Field::Element_ptr A, const size_t lda,
           typename Field::Element_ptr x, const int incx,
           typename Field::ConstElement_ptr b, const int incb, PSHelper& psH);

    template <class Field>
    typename Field::Element_ptr
    pSolve (const Field& F, const size_t M,
            typename Field::Element_ptr A, const size_t lda,
            typename Field::Element_ptr x, const int incx,
            typename Field::ConstElement_ptr b, const int incb, size_t numthreads = 0);

    //! Solve L X = B or X L = B in place.
    //! L is M*M if Side == FFLAS::FflasLeft and N*N if Side == FFLAS::FflasRight, B is M*N.
    //! Only the R non trivial column of L are stored in the M*R matrix L
    //! Requirement :  so that L could  be expanded in-place
    /* \cond */
    template<class Field>
    void
    solveLB( const Field& F, const FFLAS::FFLAS_SIDE Side,
             const size_t M, const size_t N, const size_t R,
             typename Field::Element_ptr L, const size_t ldl,
             const size_t * Q,
             typename Field::Element_ptr B, const size_t ldb );

    //! Solve L X = B in place.
    //! L is M*M or N*N, B is M*N.
    //! Only the R non trivial column of L are stored in the M*R matrix L
    template<class Field>
    void
    solveLB2( const Field& F, const FFLAS::FFLAS_SIDE Side,
              const size_t M, const size_t N, const size_t R,
              typename Field::Element_ptr L, const size_t ldl,
              const size_t * Q,
              typename Field::Element_ptr B, const size_t ldb );
    /* \endcond */

    /*************/
    /* NULLSPACE */
    /*************/

    /**  Computes a vector of the Left/Right nullspace of the matrix A.
     *
     * @param F The computation domain
     * @param Side decides whether it computes the left (FflasLeft) or right (FflasRight) nullspace
     * @param M number of rows
     * @param N number of columns
     * @param[in,out] A input matrix of dimension M x N, A is modified to its LU version
     * @param lda leading dimension of A
     * @param[out] X output vector
     * @param incX increment of X
     *
     */
    template <class Field>
    void RandomNullSpaceVector (const Field& F, const FFLAS::FFLAS_SIDE Side,
                                const size_t M, const size_t N,
                                typename Field::Element_ptr A, const size_t lda,
                                typename Field::Element_ptr X, const size_t incX);

    /**  Computes a basis of the Left/Right nullspace of the matrix A.
     * return the dimension of the nullspace.
     *
     * @param F The computation domain
     * @param Side decides whether it computes the left (FflasLeft) or right (FflasRight) nullspace
     * @param M number of rows
     * @param N number of columns
     * @param[in,out] A input matrix of dimension M x N, A is modified
     * @param lda leading dimension of A
     * @param[out] NS output matrix of dimension N x NSdim (allocated here)
     * @param[out] ldn leading dimension of NS
     * @param[out] NSdim the dimension of the Nullspace (N-rank(A))
     *
     */
    template <class Field>
    size_t NullSpaceBasis (const Field& F, const FFLAS::FFLAS_SIDE Side,
                           const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           typename Field::Element_ptr& NS, size_t& ldn,
                           size_t& NSdim);

    /*****************/
    /* RANK PROFILES */
    /*****************/

    /** @brief Computes the row rank profile of A.
     *
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param [in] A input matrix of dimension M x N
     * @param lda leading dimension of A
     * @param [out] rkprofile return the rank profile as an array of row indexes, of dimension r=rank(A)
     * @param LuTag: chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
     *
     * A is modified
     * rkprofile is allocated during the computation.
     * @returns R
     */
    template <class Field>
    size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           size_t* &rkprofile, const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    template <class Field>
    size_t pRowRankProfile (const Field& F, const size_t M, const size_t N,
                           typename Field::Element_ptr A, const size_t lda,
                           size_t* &rkprofile, size_t numthreads = 0, const FFPACK_LU_TAG LuTag=FfpackTileRecursive);

    template <class Field, class PSHelper>
    size_t RowRankProfile (const Field& F, const size_t M, const size_t N,
                                  typename Field::Element_ptr A, const size_t lda,
                                  size_t* &rkprofile, const FFPACK_LU_TAG LuTag, PSHelper& psH);

    /**  @brief Computes the column rank profile of A.
     *
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param [in] A input matrix of dimension
     * @param lda leading dimension of A
     * @param [out] rkprofile return the rank profile as an array of row indexes, of dimension r=rank(A)
     * @param LuTag: chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
     *
     * A is modified
     * rkprofile is allocated during the computation.
     * @returns R
     */
    template <class Field>
    size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
                              typename Field::Element_ptr A, const size_t lda,
                              size_t* &rkprofile, const FFPACK_LU_TAG LuTag=FfpackSlabRecursive);

    template <class Field>
    size_t pColumnRankProfile (const Field& F, const size_t M, const size_t N,
                               typename Field::Element_ptr A, const size_t lda,
                               size_t* &rkprofile, size_t numthreads = 0,
                               const FFPACK_LU_TAG LuTag=FfpackTileRecursive);

    template <class Field, class PSHelper>
    size_t ColumnRankProfile (const Field& F, const size_t M, const size_t N,
                              typename Field::Element_ptr A, const size_t lda,
                              size_t* &rkprofile, const FFPACK_LU_TAG LuTag, PSHelper& psH);


    /**  @brief Recovers the column/row rank profile from the permutation of an LU decomposition.
     *
     * Works with both the CUP/PLE decompositions (obtained by LUdivine) or the PLUQ decomposition.
     * Assumes that the output vector containing the rank profile is already allocated.
     * @param P the permutation carrying the rank profile information
     * @param N the row/col dimension for a row/column rank profile
     * @param R the rank of the matrix
     * @param [out] rkprofile return the rank profile as an array of indices
     * @param LuTag: chooses the elimination algorithm. SlabRecursive for LUdivine, TileRecursive for PLUQ
     *
     *
     */
    void RankProfileFromLU (const size_t* P, const size_t N, const size_t R,
                            size_t* rkprofile, const FFPACK_LU_TAG LuTag);

    /**  @brief Recovers the row and column rank profiles of any leading submatrix from the PLUQ decomposition.
     *
     * Only works with the PLUQ decomposition
     * Assumes that the output vectors containing the rank profiles are already allocated.
     *
     * @param P the permutation carrying the rank profile information
     * @param M the row dimension of the initial matrix
     * @param N the column dimension of the initial matrix
     * @param R the rank of the initial matrix
     * @param LSm the row dimension of the leading submatrix considered
     * @param LSn the column dimension of the leading submatrix considered
     * @param P the row permutation of the PLUQ decomposition
     * @param Q the column permutation of the PLUQ decomposition
     * @param RRP return the row rank profile of the leading submatrix
     * @return the rank of the LSm x LSn leading submatrix
     *
     * A is modified
     * @bib
     * - Dumas J-G., Pernet C., and Sultan Z. <i>\c Simultaneous computation of the row and column rank profiles </i>, ISSAC'13.
     */
    size_t LeadingSubmatrixRankProfiles (const size_t M, const size_t N, const size_t R,
                                         const size_t LSm, const size_t LSn,
                                         const size_t* P, const size_t* Q,
                                         size_t* RRP, size_t* CRP);
    /** RowRankProfileSubmatrixIndices.
     * Computes the indices of the submatrix r*r X of A whose rows correspond to
     * the row rank profile of A.
     *
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param [in] A input matrix of dimension
     * @param rowindices array of the row indices of X in A
     * @param colindices array of the col indices of X in A
     * @param lda leading dimension of A
     * @param[out] R list of indices
     *
     * rowindices and colindices are allocated during the computation.
     * A is modified
     * @returns R
     */
    template <class Field>
    size_t RowRankProfileSubmatrixIndices (const Field& F,
                                           const size_t M, const size_t N,
                                           typename Field::Element_ptr A,
                                           const size_t lda,
                                           size_t*& rowindices,
                                           size_t*& colindices,
                                           size_t& R);

    /** Computes the indices of the submatrix r*r X of A whose columns correspond to
     * the column rank profile of A.
     *
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param [in] A input matrix of dimension
     * @param rowindices array of the row indices of X in A
     * @param colindices array of the col indices of X in A
     * @param lda leading dimension of A
     * @param[out] R list of indices
     *
     * rowindices and colindices are allocated during the computation.
     * @warning A is modified
     * \return R
     */
    template <class Field>
    size_t ColRankProfileSubmatrixIndices (const Field& F,
                                           const size_t M, const size_t N,
                                           typename Field::Element_ptr A,
                                           const size_t lda,
                                           size_t*& rowindices,
                                           size_t*& colindices,
                                           size_t& R);

    /** Computes the r*r submatrix X of A, by picking the row rank profile rows of A.
     *
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param [in] A input matrix of dimension M x N
     * @param lda leading dimension of A
     * @param [out] X the output matrix
     * @param[out] R list of indices
     *
     * A is not modified
     * X is allocated during the computation.
     * @return R
     */
    template <class Field>
    size_t RowRankProfileSubmatrix (const Field& F,
                                    const size_t M, const size_t N,
                                    typename Field::Element_ptr A,
                                    const size_t lda,
                                    typename Field::Element_ptr& X, size_t& R);

    /** Compute the \f$ r\times r\f$ submatrix X of A, by picking the row rank profile rows of A.
     *
     *
     * @param F base field
     * @param M number of rows
     * @param N number of columns
     * @param[in] A input matrix of dimension M x N
     * @param lda leading dimension of A
     * @param[out] X the output matrix
     * @param[out] R list of indices
     *
     * A is not modified
     * X is allocated during the computation.
     * \returns R
     */
    template <class Field>
    size_t ColRankProfileSubmatrix (const Field& F, const size_t M, const size_t N,
                                    typename Field::Element_ptr A, const size_t lda,
                                    typename Field::Element_ptr& X, size_t& R);

    /*********************************************/
    /* Accessors to Triangular and Echelon forms */
    /*********************************************/

    /** Extracts a triangular matrix from a compact storage A=L\U of rank R.
     * if OnlyNonZeroVectors is false, then T and A have the same dimensions
     * Otherwise, T is R x N if UpLo = FflasUpper, else T is  M x R
     * @param F: base field
     * @param UpLo: selects if the upper (FflasUpper) or lower (FflasLower) triangular matrix is returned
     * @param diag: selects if the triangular matrix unit-diagonal (FflasUnit/NoUnit)
     * @param M: row dimension of T
     * @param N: column dimension of T
     * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
     * @param[in] A: input matrix
     * @param lda: leading dimension of A
     * @param[out] T: output matrix
     * @param ldt: leading dimension of T
     * @param OnlyNonZeroVectors: decides whether the last zero rows/columns should be ignored
     */
    template <class Field>
    void
    getTriangular (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                   const FFLAS::FFLAS_DIAG diag,
                   const size_t M, const size_t N, const size_t R,
                   typename Field::ConstElement_ptr A, const size_t lda,
                   typename Field::Element_ptr T, const size_t ldt,
                   const bool OnlyNonZeroVectors = false);

    /** Cleans up a compact storage A=L\U to reveal a triangular matrix of rank R.
     * @param F: base field
     * @param UpLo: selects if the upper (FflasUpper) or lower (FflasLower) triangular matrix is revealed
     * @param diag: selects if the triangular matrix unit-diagonal (FflasUnit/NoUnit)
     * @param M: row dimension of A
     * @param N: column dimension of A
     * @param R: rank of the triangular matrix
     * @param[inout] A: input/output matrix
     * @param lda: leading dimension of A
     */
    template <class Field>
    void
    getTriangular (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                   const FFLAS::FFLAS_DIAG diag,
                   const size_t M, const size_t N, const size_t R,
                   typename Field::Element_ptr A, const size_t lda);

    /** Extracts a matrix in echelon form from a compact storage A=L\U of rank R obtained by
     * RowEchelonForm or ColumnEchelonForm.
     * Either L or U is in Echelon form (depending on Uplo)
     * The echelon structure is defined by the first R values of the array P.
     * row and column dimension of T are greater or equal to that of A
     * @param F: base field
     * @param UpLo: selects if the upper (FflasUpper) or lower (FflasLower) triangular matrix is returned
     * @param diag: selects if the echelon matrix has unit pivots (FflasUnit/NoUnit)
     * @param M: row dimension of T
     * @param N: column dimension of T
     * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
     * @param P: positions of the R pivots
     * @param[in] A: input matrix
     * @param lda: leading dimension of A
     * @param[out] T: output matrix
     * @param ldt: leading dimension of T
     * @param OnlyNonZeroVectors: decides whether the last zero rows/columns should be ignored
     * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
     */
    template <class Field>
    void
    getEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                    const FFLAS::FFLAS_DIAG diag,
                    const size_t M, const size_t N, const size_t R, const size_t* P,
                    typename Field::ConstElement_ptr A, const size_t lda,
                    typename Field::Element_ptr T, const size_t ldt,
                    const bool OnlyNonZeroVectors = false,
                    const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

    /** Cleans up a compact storage A=L\U obtained by RowEchelonForm or ColumnEchelonForm
     * to reveal an echelon form of rank R.
     * Either L or U is in Echelon form (depending on Uplo)
     * The echelon structure is defined by the first R values of the array P.
     * @param F: base field
     * @param UpLo: selects if the upper (FflasUpper) or lower (FflasLower) triangular matrix is returned
     * @param diag: selects if the echelon matrix has unit pivots (FflasUnit/NoUnit)
     * @param M: row dimension of A
     * @param N: column dimension of A
     * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
     * @param P: positions of the R pivots
     * @param[inout] A: input/output matrix
     * @param lda: leading dimension of A
     * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
     */
    template <class Field>
    void
    getEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                    const FFLAS::FFLAS_DIAG diag,
                    const size_t M, const size_t N, const size_t R, const size_t* P,
                    typename Field::Element_ptr A, const size_t lda,
                    const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

    /** Extracts a transformation matrix to echelon form from a compact storage A=L\U
     * of rank R obtained by RowEchelonForm or ColumnEchelonForm.
     * If Uplo == FflasLower:
     *   T is N x N (already allocated) such that A T = C is a transformation of A in
     *   Column echelon form
     * Else
     *   T is M x M (already allocated) such that T A = E is a transformation of A in
     *   Row Echelon form
     * @param F: base field
     * @param UpLo: Lower (FflasLower) means Transformation to Column Echelon Form, Upper (FflasUpper), to Row Echelon Form
     * @param diag: selects if the echelon matrix has unit pivots (FflasUnit/NoUnit)
     * @param M: row dimension of A
     * @param N: column dimension of A
     * @param R: rank of the triangular matrix
     * @param P: permutation matrix
     * @param[in] A: input matrix
     * @param lda: leading dimension of A
     * @param[out] T: output matrix
     * @param ldt: leading dimension of T
     * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
     */
    template <class Field>
    void
    getEchelonTransform (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                         const FFLAS::FFLAS_DIAG diag,
                         const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
                         typename Field::ConstElement_ptr A, const size_t lda,
                         typename Field::Element_ptr T, const size_t ldt,
                         const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);
    /** Extracts a matrix in echelon form from a compact storage A=L\U of rank R obtained by
     * ReducedRowEchelonForm or ReducedColumnEchelonForm with transform = true.
     * Either L or U is in Echelon form (depending on Uplo)
     * The echelon structure is defined by the first R values of the array P.
     * row and column dimension of T are greater or equal to that of A
     * @param F: base field
     * @param UpLo: selects if the upper (FflasUpper) or lower (FflasLower) triangular matrix is returned
     * @param diag: selects if the echelon matrix has unit pivots (FflasUnit/NoUnit)
     * @param M: row dimension of T
     * @param N: column dimension of T
     * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
     * @param P: positions of the R pivots
     * @param[in] A: input matrix
     * @param lda: leading dimension of A
     * @param ldt: leading dimension of T
     * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
     * @param OnlyNonZeroVectors: decides whether the last zero rows/columns should be ignored
     */
    template <class Field>
    void
    getReducedEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                           const size_t M, const size_t N, const size_t R, const size_t* P,
                           typename Field::ConstElement_ptr A, const size_t lda,
                           typename Field::Element_ptr T, const size_t ldt,
                           const bool OnlyNonZeroVectors = false,
                           const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

    /** Cleans up a compact storage A=L\U of rank R obtained by ReducedRowEchelonForm or
     * ReducedColumnEchelonForm with transform = true.
     * Either L or U is in Echelon form (depending on Uplo)
     * The echelon structure is defined by the first R values of the array P.
     * @param F: base field
     * @param UpLo: selects if the upper (FflasUpper) or lower (FflasLower) triangular matrix is returned
     * @param diag: selects if the echelon matrix has unit pivots (FflasUnit/NoUnit)
     * @param M: row dimension of A
     * @param N: column dimension of A
     * @param R: rank of the triangular matrix (how many rows/columns need to be copied)
     * @param P: positions of the R pivots
     * @param[inout] A: input/output matrix
     * @param lda: leading dimension of A
     * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
     */
    template <class Field>
    void
    getReducedEchelonForm (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                           const size_t M, const size_t N, const size_t R, const size_t* P,
                           typename Field::Element_ptr A, const size_t lda,
                           const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);

    /** Extracts a transformation matrix to echelon form from a compact storage A=L\U
     * of rank R obtained by RowEchelonForm or ColumnEchelonForm.
     * If Uplo == FflasLower:
     *   T is N x N (already allocated) such that A T = C is a transformation of A in
     *   Column echelon form
     * Else
     *   T is M x M (already allocated) such that T A = E is a transformation of A in
     *   Row Echelon form
     * @param F: base field
     * @param UpLo: selects Col (FflasLower) or Row (FflasUpper) Echelon Form
     * @param diag: selects if the echelon matrix has unit pivots (FflasUnit/NoUnit)
     * @param M: row dimension of A
     * @param N: column dimension of A
     * @param R: rank of the triangular matrix
     * @param P: permutation matrix
     * @param[in] A: input matrix
     * @param lda: leading dimension of A
     * @param[out] T: output matrix
     * @param ldt: leading dimension of T
     * @param LuTag: which factorized form (CUP/PLE if FfpackSlabRecursive, PLUQ if FfpackTileRecursive)
     */
    template <class Field>
    void
    getReducedEchelonTransform (const Field& F, const FFLAS::FFLAS_UPLO Uplo,
                                const size_t M, const size_t N, const size_t R, const size_t* P, const size_t* Q,
                                typename Field::ConstElement_ptr A, const size_t lda,
                                typename Field::Element_ptr T, const size_t ldt,
                                const FFPACK_LU_TAG LuTag = FfpackSlabRecursive);
    /** Auxiliary routine: determines the permutation that changes a PLUQ decomposition
     * into a echelon form revealing PLUQ decomposition
     */
    void
    PLUQtoEchelonPermutation (const size_t N, const size_t R, const size_t * P, size_t * outPerm);

} // FFPACK
// #include "ffpack.inl"


namespace FFPACK { /* Quasi-separable matrices*/
        // TODO add signatures of ffpack_bruhat.inl
    /** LTBruhatGen
     * Suppose A is Left Triangular Matrix
     * This procedure computes the Bruhat Representation of A and return the rank of A
     * @param Fi base Field
     * @param diag
     * @param N size of A
     * @param A the matrix we search the Bruhat representation
     * @param lda the leading dimension of A
     * @param P a permutation matrix
     * @param Q a permutation matrix
     */
  template<class Field>
  size_t LTBruhatGen (const Field& Fi, const FFLAS::FFLAS_DIAG diag,
                             const size_t N,
                             typename Field::Element_ptr A, const size_t lda,
                             size_t * P, size_t * Q);
    /** GetLTBruhatGen
     * This procedure Computes the Rank Revealing Matrix based on the Bruhta representation of a Matrix
     * @param Fi base Field
     * @param N size of the matrix
     * @param r the rank of the matrix
     * @param P a permutation matrix
     * @param Q a permutation matrix
     * @param R the matrix that will contain the rank revealing matrix
     * @param ldr the leading fimension of R
     */
  template<class Field>
  void getLTBruhatGen(const Field& Fi, const size_t N, const size_t r,const size_t * P, const size_t * Q, typename Field::Element_ptr R, const size_t ldr);
    /** GetLTBruhatGen
     * This procedure computes the matrix L or U f the Bruhat Representation
     * Suppose that A is the bruhat representation of a matrix
     * @param Fi base Field
     * @param Uplo choose if the procedure return L or U
     * @param diag 
     * @param N size of A
     * @param r rank of A
     * @param P permutaion matrix
     * @param Q permutation matrix
     * @param A a bruhat representation
     * @param lda leading dimension of A
     * @param T matrix that will contains L or U
     * @param ldt leading dimension of T
     */
  template<class Field>
   void getLTBruhatGen(const Field& Fi, const FFLAS::FFLAS_UPLO Uplo,const FFLAS::FFLAS_DIAG diag ,const size_t N, const size_t r, const size_t *P, const size_t * Q, typename Field::ConstElement_ptr A, const size_t lda, typename Field::Element_ptr T, const size_t ldt);

    /** LTQSorder
     * This procedure computes the order of quasiseparability of a matrix
     * @param N size of the matrix
     * @param r rank of the matrix
     * @param P permutation matrix
     * @param Q permutation matrix
     */
  size_t LTQSorder(const size_t N, const size_t r,const size_t * P, const size_t * Q);

    /**CompressToBlockBiDiagonal
     * This procedure compress a compact representation of a row echelon form or column echelon form
     * @param Fi base Field
     * @param Uplo chosse if the procedure is based on row or column
     * @param N size of the matrix
     * @param s order of qausiseparability
     * @param r rank
     * @param P permutation matrix
     * @param Q permutation matrix
     * @param A the matrix to compact
     * @param lda leading dimension of A
     * @param X matrix that will stock the representation
     * @param ldx leading dimension of X
     * @param K stock the position of the blocks in A
     * @param M permutation matrix
     * @param T stock the operation done in the procedure
     */
  template<class Field>
  size_t CompressToBlockBiDiagonal(const Field&Fi, const FFLAS::FFLAS_UPLO Uplo, size_t N, size_t s, size_t r, const size_t *P, const size_t *Q,  typename Field::Element_ptr A, size_t lda, typename Field::Element_ptr X, size_t ldx, size_t *K, size_t *M, size_t *T);

  /**ExpandBlockBiDiagonal
     * This procedure expand a compact representation of a row echelon form or column echelon form
     * @param Fi base Field
     * @param Uplo chosse if the procedure is based on row or column
     * @param N size of the matrix
     * @param s order of qausiseparability
     * @param r rank
     * @param A the matrix that will sotck the expanded representation
     * @param lda leading dimension of A
     * @param X matrix to expand
     * @param ldx leading dimension of X
     * @param K stock the position of the blocks in A
     * @param M permutation matrix
     * @param T stock the operation done in the procedure
     */
  template<class Field>
  void ExpandBlockBiDiagonalToBruhat(const Field&Fi, const FFLAS::FFLAS_UPLO Uplo, size_t N, size_t s, size_t r, typename Field::Element_ptr A, size_t lda, typename Field::Element_ptr X, size_t ldx,size_t NbBlocks,size_t *K, size_t *M, size_t *T);

   /**Bruhat2EchelonPermutation (N,R,P,Q)
    * Compute M such that LM or MU is in echelon form where L or U are factors of the Bruhat Rpresentation
    * @param[in] N size of the matrix
    * @param[in] R rank
    * @param[in] P permutation Matrix
    * @param[in] Q permutation Matrix
    * @param[out] M output permutation matrix
    */
    void Bruhat2EchelonPermutation (size_t N,size_t R, const size_t* P,const size_t *Q, size_t * M);


  size_t * TInverter (size_t * T, size_t r);

  template<class Field>
  void ComputeRPermutation (const Field&Fi, size_t N, size_t r, const size_t * P, const size_t * Q, size_t * R,size_t * MU, size_t * ML);
  /**productBruhatxTS
   *Comput the product between the CRE compact representation of a matrix A and B a tall matrix
   */

  template<class Field>
  void productBruhatxTS (const Field&Fi, size_t N, size_t s, size_t r, const size_t *P, const size_t *Q,  const typename Field::Element_ptr Xu,size_t ldu, size_t NbBlocksU, size_t * Ku, size_t *Tu, size_t * MU,const typename Field::Element_ptr Xl, size_t ldl, size_t NbBlocksL,size_t *Kl, size_t * Tl, size_t * ML,typename Field::Element_ptr B,size_t t, size_t ldb, typename Field::Element_ptr C, size_t ldc);
}

namespace FFPACK { /* not used */

    /** LQUPtoInverseOfFullRankMinor.
     * Suppose A has been factorized as L.Q.U.P, with rank r.
     * Then Qt.A.Pt has an invertible leading principal r x r submatrix
     * This procedure efficiently computes the inverse of this minor and puts it into X.
     * @note It changes the lower entries of A_factors in the process (NB: unless A was nonsingular and square)
     *
     * @param F base field
     * @param rank       rank of the matrix.
     * @param A_factors  matrix containing the L and U entries of the factorization
     * @param lda leading dimension of A
     * @param QtPointer  theLQUP->getQ()->getPointer() (note: getQ returns Qt!)
     * @param X          desired location for output
     * @param ldx leading dimension of X
     */
    template <class Field>
    typename Field::Element_ptr
    LQUPtoInverseOfFullRankMinor( const Field& F, const size_t rank,
                                  typename Field::Element_ptr A_factors, const size_t lda,
                                  const size_t* QtPointer,
                                  typename Field::Element_ptr X, const size_t ldx);

} // FFPACK
// include precompiled instantiation headers (avoiding to recompile them)
#ifdef FFPACK_COMPILED
#include "fflas-ffpack/interfaces/libs/ffpack_inst.h"
#endif

//---------------------------------------------------------------------
// Checkers
#include "fflas-ffpack/checkers/checkers_ffpack.h"
//---------------------------------------------------------------------

#include "ffpack_fgesv.inl"
#include "ffpack_fgetrs.inl"
//---------------------------------------------------------------------
// Checkers
#include "fflas-ffpack/checkers/checkers_ffpack.inl"
//---------------------------------------------------------------------
#include "ffpack_pluq.inl"
#include "ffpack_pluq_mp.inl"
#include "ffpack_ppluq.inl"
#include "ffpack_ludivine.inl"
#include "ffpack_ludivine_mp.inl"
#include "ffpack_echelonforms.inl"
#include "ffpack_fsytrf.inl"
#include "ffpack_invert.inl"
#include "ffpack_ftrtr.inl"
#include "ffpack_ftrstr.inl"
#include "ffpack_ftrssyr2k.inl"
#include "ffpack_charpoly_kglu.inl"
#include "ffpack_charpoly_kgfast.inl"
#include "ffpack_charpoly_kgfastgeneralized.inl"
#include "ffpack_charpoly_danilevski.inl"
#include "ffpack_charpoly.inl"
#include "ffpack_frobenius.inl"
#include "ffpack_minpoly.inl"
#include "ffpack_krylovelim.inl"
#include "ffpack_permutation.inl"
#include "ffpack_rankprofiles.inl"
#include "ffpack_det_mp.inl"
#include "ffpack_bruhatgen.inl"
#include "ffpack.inl"

#endif // __FFLASFFPACK_ffpack_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

