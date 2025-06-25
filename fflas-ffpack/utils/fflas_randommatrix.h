/*
 * Copyright (C) FFLAS-FFPACK
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * This file is Free Software and part of FFLAS-FFPACK.
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

/* @file utils/fflas_randommatrix.h
 * @ingroup tests
 * @brief Utilities to create matrices with prescribed shapes, properties,...
 * To be used in benchmarks/tests
 */

#ifndef __FFLASFFPACK_randommatrix_H
#define __FFLASFFPACK_randommatrix_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/debug.h"
#include "fflas-ffpack/fflas/fflas.h"
#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include <givaro/givranditer.h>

namespace FFPACK {

    /** @brief  Random non-zero Matrix.
     * Creates a \c m x \c n matrix with random entries, and at least one of them is non zero.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    NonZeroRandomMatrix(const Field & F, size_t m, size_t n, typename Field::Element_ptr A, size_t lda, RandIter& G) {
        bool ok=false;
        while (!ok)
            for (size_t i=0 ; i<m ; ++i)
                for (size_t j= 0; j<n ;++j)
                    if (!F.isZero(G.random (A[i*lda+j])))
                        ok = true;
        return A;
    }

    /** @brief  Random non-zero Matrix.
     * Creates a \c m x \c n matrix with random entries, and at least one of them is non zero.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    NonZeroRandomMatrix(const Field & F, size_t m, size_t n,
                        typename Field::Element_ptr A, size_t lda) {
        typename Field::RandIter G(F);
        return NonZeroRandomMatrix(F, m, n, A, lda, G);
    }

    /** @brief  Random Matrix.
     * Creates a \c m x \c n matrix with random entries.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomMatrix(const Field & F, size_t m, size_t n, typename Field::Element_ptr A, size_t lda, RandIter& G) {
        for (size_t i=0 ; i<m ; ++i)
            for (size_t j= 0; j<n ;++j)
                G.random (A[i*lda+j]);
        return A;
    }

    /** @brief  Random Matrix.
     * Creates a \c m x \c n matrix with random entries.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomMatrix(const Field & F, size_t m, size_t n, typename Field::Element_ptr A, size_t lda) {
        typename Field::RandIter G(F);
        return RandomMatrix (F, m, n, A, lda, G);
    }

    /** @brief  Random Triangular Matrix.
     * Creates a \c m x \c n triangular matrix with random entries. The \c UpLo parameter defines wether it is upper or lower triangular.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param UpLo whether \c A is upper or lower triangular
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomTriangularMatrix (const Field & F, size_t m, size_t n,
                            const FFLAS::FFLAS_UPLO UpLo, const FFLAS::FFLAS_DIAG Diag, bool nonsingular,
                            typename Field::Element_ptr A, size_t lda, RandIter& G) {

        if (UpLo == FFLAS::FflasUpper){
            for (size_t i=0 ; i<m ; ++i){
                FFLAS::fzero(F, std::min(i,n), A + i*lda, 1);
                for (size_t j= i; j<n ;++j)
                    G.random (A[i*lda+j]);
            }
        } else { // FflasLower
            for (size_t i=0 ; i<m ; ++i){
                for (size_t j=0; j<=i ;++j)
                    G.random (A[i*lda+j]);
                FFLAS::fzero (F, n-1-std::min(i,n-1), A + i*lda + i+1, 1);
            }
        }
        if (Diag == FFLAS::FflasUnit){
            for (size_t i=0; i< std::min(m,n); ++i)
                F.assign(A[i*(lda+1)], F.one);
        } else { // NonUnit
            if (nonsingular){
                Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzG (G);
                for (size_t i=0; i< std::min(m,n); ++i)
                    nzG.random (A[i*(lda+1)]);
            }
        }
        return A;
    }
    /** @brief  Random Triangular Matrix.
     * Creates a \c m x \c n triangular matrix with random entries. The \c UpLo parameter defines wether it is upper or lower triangular.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param UpLo whether \c A is upper or lower triangular
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomTriangularMatrix (const Field & F, size_t m, size_t n,
                            const FFLAS::FFLAS_UPLO UpLo, const FFLAS::FFLAS_DIAG Diag, bool nonsingular,
                            typename Field::Element_ptr A, size_t lda) {
        typename Field::RandIter G(F);
        return RandomTriangularMatrix (F, m, n, UpLo, Diag, nonsingular, A, lda, G);
    }

    /* Random integer in range.
     * @param a min bound
     * @param b max bound
     * @return a random integer in [a,b[  */
    inline size_t RandInt(size_t a, size_t b)
    {
        size_t x = a ;
        x += (size_t)rand()%(b-a);
        FFLASFFPACK_check(x<b && x>=a);
        return x ;
    }
    /** @brief  Random Symmetric Matrix.
     * Creates a \c m x \c n triangular matrix with random entries. The \c UpLo parameter defines wether it is upper or lower triangular.
     * @param F field
     * @param n order of \p A
     * @param [out] A the matrix (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomSymmetricMatrix (const Field & F,size_t n, bool nonsingular,
                           typename Field::Element_ptr A, size_t lda, RandIter& G) {
        RandomTriangularMatrix (F, n, n, FFLAS::FflasUpper, FFLAS::FflasNonUnit, nonsingular, A, lda, G);
        for (size_t i=0; i<n; i++){
            typename Field::Element piv = A[i*(lda+1)];
            if (!F.isZero(piv)){
                typename Field::Element inv;
                F.init(inv);
                F.inv(inv, A[i*(lda+1)]);
                FFLAS::fscal(F, n-i-1, inv, A+i*(lda+1)+1, 1, A+i*(lda+1)+lda, lda);
            }
        }
        ftrtrm (F, FFLAS::FflasRight, FFLAS::FflasNonUnit, n, A, lda);
        return A;
    }
} // FFPACK

#include "fflas-ffpack/ffpack/ffpack.h"

namespace FFPACK{
    /** @brief  Random Matrix with prescribed rank.
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomMatrixWithRank (const Field & F, size_t m, size_t n, size_t r,
                          typename Field::Element_ptr A, size_t lda, RandIter& G){
        FFLASFFPACK_check(r <= std::min(m,n));
        FFLASFFPACK_check(n <= lda);
        typedef typename Field::Element_ptr  Element_ptr;

        Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzG (G);

        size_t * P = FFLAS::fflas_new<size_t>(n);
        size_t * Q = FFLAS::fflas_new<size_t>(m);
        for (size_t i = 0 ; i < m ; ++i ) Q[i] = 0;
        for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

        Element_ptr U = FFLAS::fflas_new(F,m,n);
        Element_ptr L = FFLAS::fflas_new(F,m,m);
        /*  Create L, lower invertible */
        RandomTriangularMatrix (F, m, m, FFLAS::FflasLower, FFLAS::FflasNonUnit, true, L, m, G);
        /*  Create U, upper or rank r */
        RandomTriangularMatrix (F, m, n, FFLAS::FflasUpper, FFLAS::FflasNonUnit, true, U, n, G);

        /*  Create a random P,Q */
        for (size_t i = 0 ; i < n ; ++i)
            P[i] = i + RandInt(0U,n-i);
        for (size_t i = 0 ; i < m ; ++i)
            Q[i] = i + RandInt(0U,m-i);

        /*  compute product */

        FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
                        m,0,(int)n, U, n, P);
        FFPACK::applyP (F, FFLAS::FflasLeft,  FFLAS::FflasNoTrans,
                        m,0,(int)m, L, m, Q);
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                      m, n, m, F.one, L, m, U, n, F.zero, A, lda);
        // @todo compute LU with ftrtr

        FFLAS::fflas_delete(P);
        FFLAS::fflas_delete(L);
        FFLAS::fflas_delete(U);
        FFLAS::fflas_delete(Q);

        return A;
    }

    /** @brief  Random Matrix with prescribed rank.
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param [out] A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomMatrixWithRank (const Field & F, size_t m, size_t n, size_t r,
                          typename Field::Element_ptr A, size_t lda){
        typename Field::RandIter G(F);
        return RandomMatrixWithRank(F, m, n, r, A, lda, G);
    }

    /** @brief Pick uniformly at random a sequence of \c R distinct elements from the set \f$ \{0,\dots, N-1\}\f$  using Knuth's shuffle.
     * @param N the cardinality of the sampling set
     * @param R the number of elements to sample
     * @param [out] P the output sequence (pre-allocated to at least R indices)
     */
    inline size_t * RandomIndexSubset (size_t N, size_t R, size_t* P){
        size_t * Q = FFLAS::fflas_new<size_t>(N);
        for (size_t i=0; i<N; ++i)
            Q[i] = i;
        for (size_t i=0; i<R; ++i){
            size_t j = RandInt(i,N);
            P[i] = Q[j];
            Q[j] = Q[i];
        }
        FFLAS::fflas_delete(Q);
        return P;
    }

    /** @brief Pick uniformly at random a permutation of size \c N stored in LAPACK format using Knuth's shuffle.
     * @param N the length  of the permutation
     * @param [out] P the output permutation (pre-allocated to at least N indices)
     */
    inline size_t * RandomPermutation (size_t N, size_t* P){
        for (size_t i = 0 ; i < N ; ++i)
            P[i] = i + RandInt(0U,N-i);
        return P;
    }

    /** @brief Pick uniformly at random an R-subpermutation of dimension \c M x \c N : a matrix with only R non-zeros equal to one, in a random rook placement.
     * @param M row dimension
     * @param N column dimension
     * @param [out] rows the row position of each non zero element (pre-allocated)
     * @param [out] cols the column position of each non zero element (pre-allocated)
     */
    inline void RandomRankProfileMatrix (size_t M, size_t N, size_t R, size_t* rows, size_t* cols){
        RandomIndexSubset (M, R, rows);
        RandomIndexSubset (N, R, cols);
    }

    inline void swapval(size_t k, size_t N, size_t * P, size_t val){
        size_t i = k;
        int found =-1;
        do {
            if (P[i] == val)
                found = i;
            i++;
        } while(found<0);
        P[found] = P[k];
    }
    /** @brief Pick uniformly at random a symmetric R-subpermutation of dimension \c N x \c N : a symmetric matrix with only R non-zeros, all equal to one, in a random rook placement.
     * @param N matrix order
     * @param [out] rows the row position of each non zero element (pre-allocated)
     * @param [out] cols the column position of each non zero element (pre-allocated)
     */
    inline void RandomSymmetricRankProfileMatrix (size_t N, size_t R, size_t* rows, size_t* cols){

        size_t * rr = FFLAS::fflas_new<size_t>(N);
        size_t * cc = FFLAS::fflas_new<size_t>(N);
        for (size_t i=0; i<N; ++i)
            rr[i] = cc[i] = i;
        for (size_t k=0; k<R; k+=2){
            size_t i = RandInt(k,N);
            size_t j = RandInt(k,N);
            cols[k] = cc[j];
            cc[j] = cc[k];
            rows[k] = rr[i];
            rr[i] = rr[k];
            if (rows[k] != cols[k] && k < R-1){
                // adding the symmetric element
                rows[k+1] = cols[k];
                cols[k+1] = rows[k];
                swapval(k+1,N,rr,cols[k]);
                swapval(k+1,N,cc,rows[k]);
            } else {
                // we need to add a diagonal pivot since
                // - either k==R-1 and there is only one pivot left to be added
                // - or we just added a diagonal pivot. We need to pick another one so
                //   that they appear with the same probability 2/N^2 as off-diagonal pivots
                if (k<R-1) k++; //
                size_t l, co;
                int found =-1;
                do{
                    l = RandInt(k,N);
                    co = cc[l];
                    for (size_t m=k; m<N; m++)
                        if (rr[m] == co) // l is valid as row co still available
                            found = m;
                    // TODO: Write a variant for when k < N/2
                } while(found<0);
                cols[k] = co;
                cc[l] = cc[k];
                rows[k] = co;
                rr[found] = rr[k];
                if (k<R) k--; //
            }
        }
        FFLAS::fflas_delete(rr,cc);
    }

    inline void RandomLTQSRankProfileMatrix (size_t n, size_t r, size_t t, size_t * rows, size_t *cols){

        size_t * leadingRk = FFLAS::fflas_new<size_t>(n-1); //leadingRk[i] = rank of the leading (i+1) x (n-i-1) submatrix
        size_t fullRkBlocks, p, max_t;
        size_t* randrows = FFLAS::fflas_new<size_t>(n-1);
        size_t* randcols = FFLAS::fflas_new<size_t>(n-1);
        do {
                //std::cerr<<"Iteration n,r,t="<<n<<" "<<r<<" "<<t<<std::endl;
            std::vector<bool> pivot_in_row (n-1,false);
            std::vector<bool> pivot_in_col (n-1,false);
            RandomIndexSubset (n-1, n-1, randrows);
            RandomIndexSubset (n-1, n-1, randcols);
            p=0; max_t=0;
            for (size_t i=0; i<n-1; ++i) leadingRk[i]=0;
            fullRkBlocks = 0;
            for (size_t k=0; k<n-1 && p<r; ++k){
                size_t i = randrows[k];
                size_t j = randcols[k];
                if (j < n-i-1){
                    bool admissible = true;
                    for (size_t l=i; l<n-j-1; l++)
                        admissible = admissible && (leadingRk[l] < t);
                    if (admissible){
                        rows[p] = i; pivot_in_row[i]=true;
                        cols[p] = j; pivot_in_col[j]=true;
                            //std::cerr<<"Pick ("<<i<<", "<<j<<")"<<std::endl;
                        p++;
                        for (size_t l=i; l<n-j-1; l++){
                            leadingRk[l]++;
                            if (max_t < leadingRk[l]) max_t = leadingRk[l];
                            if (leadingRk[l] == std::min(std::min(t,l+1),n-1-l)) fullRkBlocks++;
                        }
                    }
                }
            }
            
                //std::cerr<<"done with Step 1: p = "<<p<<" max_t = "<<max_t<<" fullRkBlocks = "<<fullRkBlocks<<std::endl;
            std::vector<size_t> available_rows;
            std::vector<size_t> available_cols;
            for (size_t i=0; i<n-1; i++) {
                if (!pivot_in_row[i]) available_rows.push_back(i);
                if (!pivot_in_col[i]) available_cols.push_back(i);
            }
            size_t tries=0;
            size_t maxtries=n*t;
            while (p<r && fullRkBlocks < n-1 && tries<maxtries){ // finish the last pivots by random samples
                size_t i = available_rows[RandInt (0,n-1-r)];
                size_t j = available_cols[RandInt (0,n-1-r)];
                if (j<n-i-1){
                        //std::cerr<<"Pick ("<<i<<", "<<j<<")"<<std::endl;
                    tries++;
                    if (!pivot_in_row[i] && !pivot_in_col[j]){
                        bool admissible = true;
                        for (size_t k=i; k<n-j-1; k++)
                            admissible = admissible && (leadingRk[k] < t);
                        if (admissible){
                            rows[p] = i; pivot_in_row[i]=true;
                            cols[p] = j; pivot_in_col[j]=true;
                            p++;
                            pivot_in_col[j]=true;
                            for (size_t l=i; l<n-j-1; l++){
                                leadingRk[l]++;
                                if (max_t< leadingRk[l]) max_t=leadingRk[l];
                                if (leadingRk[l] == std::min(std::min(t,l+1),n-1-l)) fullRkBlocks++;
                            }
                        }
                    }
                }
            }
                //std::cerr<<"Before looping: p < r : "<<p<<" < "<<r<<" max_t < t: "<<max_t<<" < "<<t<<std::endl;
        } while (p<r || max_t <t);
        FFLAS::fflas_delete (leadingRk);
        FFLAS::fflas_delete (randcols,randrows);
    }

    // O(l*len(L)+M*len(M))
    // compute all the valid moves in the array valid_moves depending on the free areas described in L and M 
    inline void compute_valid_moves(std::vector<size_t>* valid_moves, size_t* start_up, size_t l, size_t* L, size_t lengthL, int llim, size_t m, size_t* M, size_t lengthM, int mlim, std::vector<bool> availrows, std::vector<bool> availcols) {
        valid_moves->clear();

        // Base case: Identity
        valid_moves->push_back(l);

        // Compute all the empty columns/lines
        for (size_t j = (size_t) llim+1; j < l; ++j) {
            if (availcols[j]) {
                valid_moves->push_back(j);
            }
        }

        // start_up point where the valid moves change direction
        *start_up = valid_moves->size();

        for (size_t i = (size_t) mlim+1; i < m; ++i) {
            if (availrows[i]) {
                valid_moves->push_back(i);
            }
        }
    }

    // O(1)
    // estimate what rank a submatrix is able to reach given the current configuration of the algorithm
    size_t valuate(int ti, int r, int n, int i, int Is_pivot_not_moved,int nb_pivot_after_i, int nb_pivot_bf_i, int nb_pivot_after_i_not_moved, int nb_pivot_before_i_not_moved) {
        // Compute the valuation of a given leading submatrix ti
            return ti + std::min(i + 1 - nb_pivot_bf_i, nb_pivot_after_i_not_moved - Is_pivot_not_moved) +
                        std::min(n - (i + 1) - nb_pivot_after_i, nb_pivot_before_i_not_moved - Is_pivot_not_moved);
    }

    /** @brief genenration of an Rank Profile Matrix randomly
     * @param n        Dimension
     * @param r        rank of the n*n matrix
     * @param t        quasi-separability order
     * @param rows
     * @param cols
    */
    void RandomLTQSRankProfileMatrix_Tom(size_t n, size_t r, size_t t, size_t* rows, size_t* cols) {

        if (r <= 0 || n <= 0 || t <= 0) {
            std::cout << "care, all the arguments must be positive" << std::endl;
            return;
        }

        if (r >= n) {
            std::cout << "care, you called the function with invalid arguments : r >= n IMPOSSIBLE" << std::endl;
            return;
        }

        if (t > r) {
            std::cout << "care, you called the function with invalid arguments : t >= r IMPOSSIBLE" << std::endl;
            return;
        }

        if (r + t > n) {
            std::cout << "care, you called the function with invalid arguments : r+t > n IMPOSSIBLE" << std::endl;
            return;
        }

        /// Initialisation

        size_t* H = FFLAS::fflas_new<size_t>(n-1);
        size_t* T = FFLAS::fflas_new<size_t>(n-1);
        size_t* nb_pivot_before_index = FFLAS::fflas_new<size_t>(n-1); // pivots that are on lines before line i
        size_t* nb_pivot_after_index = FFLAS::fflas_new<size_t>(n-1); // pivots that are on columns before the last column of the line i
        size_t* nb_pivot_before_index_not_moved = FFLAS::fflas_new<size_t>(n-1);  // same idea with pivots that didn't move yet
        size_t* nb_pivot_after_index_not_moved = FFLAS::fflas_new<size_t>(n-1);

        std::vector<bool> availablerows (n-1, true); // indicate wether the row i is available or not
        std::vector<bool> availablecols (n-1, true); // same idea with col j

        for (size_t i = 0; i<n-1; ++i) {
            H[i]=0;
            T[i]=0;
        }

        RandomIndexSubset(n-1,r,rows);

        for (size_t i = 0; i < r; ++i) {
            cols[i] = n - rows[i] - 2;
            T[rows[i]] = 1;
            H[rows[i]] = 1;

            availablerows[rows[i]] = false;
            availablecols[cols[i]] = false;
        }


        nb_pivot_before_index[0] = H[0];
        nb_pivot_after_index[0] = r;
        nb_pivot_after_index_not_moved[0]=r;
        nb_pivot_before_index_not_moved[0]=H[0];
        for (size_t i = 1; i < n - 1; ++i) {
            nb_pivot_before_index[i] = nb_pivot_before_index[i - 1] + H[i];
            nb_pivot_after_index[i] = r-nb_pivot_before_index[i-1];
            nb_pivot_after_index_not_moved[i] = r-nb_pivot_before_index[i-1];
            nb_pivot_before_index_not_moved[i] = nb_pivot_before_index[i];
        }

        /// loop initialisation
        size_t loop = 0;
        std::vector<size_t> valid_moves;
        size_t start_up;
        int ilim, jlim;
        int h = 0;
        int prev_h = -1;

        while ((size_t) h < r) {
            size_t prev_i = rows[h];
            size_t prev_j = cols[h];

            size_t new_i = prev_i;
            size_t new_j = prev_j;


            /// computating valid_moves
            // complex : 2n
            if (prev_h != h) {
                // compute impossible ones (those which make the configuration having a greater t than asked)
                ilim = -1;
                jlim = -1;
                int k = rows[h]+1;
                while ((size_t) k <= n-2) {
                    if (T[k] == t) {
                        jlim = n-k-2;
                        break;
                    }
                    ++k;
                }
                k = rows[h]-1;
                while (k >= 0) {
                    if (T[k] == t) {
                        ilim = k;
                        break;
                    }
                    --k;
                }
                compute_valid_moves(& valid_moves, & start_up, prev_j, cols, r, jlim, prev_i, rows, r, ilim, availablerows, availablecols);
            }

            /// choose and make the chosen move
            size_t rand_index;
            if (valid_moves.size()>1){
                rand_index = RandInt(0, valid_moves.size());
            }
            else {
                rand_index = 0;
            }

            if (rand_index < start_up) {
                new_j = valid_moves[rand_index];
                cols[h] = new_j;
                // updates of the arrays for the valuation function
                for (size_t k = prev_i + 1; k <= prev_i + (prev_j - new_j); ++k) {
                    T[k] += 1;
                    nb_pivot_after_index[k] += 1;
                }
            } else {
                new_i = valid_moves[rand_index];
                rows[h] = new_i;
                // updates of the arrays for the valuation function
                for (size_t k = new_i; k < prev_i; ++k) {
                    T[k] += 1;
                    nb_pivot_before_index[k] += 1;
                }
            }

            // updates of the arrays for the valuation function
            // O(n)
            for (size_t k = 0; k<n-1; k++){
                if (k <= prev_i) {
                    nb_pivot_after_index_not_moved[k]--;
                }
                if (k >= prev_i) {
                    nb_pivot_before_index_not_moved[k]--;
                }
            }

            H[prev_i] = 0;

            /// valutation
            // O(n)
            bool undo = true;
            for (size_t k = 0; k < n - 1; ++k) {
                if (valuate(T[k], r, n, k, H[k], nb_pivot_after_index[k], nb_pivot_before_index[k], nb_pivot_after_index_not_moved[k], nb_pivot_before_index_not_moved[k]) >= t) {
                    undo = false;
                    break;
                }
            }

            /// undo if needed
            if (undo) {
                if (rand_index >= start_up) {
                    // undo the move and erase it from the valid moves
                    rows[h] = prev_i;
                    valid_moves.erase(valid_moves.begin()+rand_index);
                    // undo of the arrays for valuation func
                    for (size_t k = new_i; k < prev_i; ++k) {
                        T[k] -= 1;
                        nb_pivot_before_index[k] -= 1;
                    }
                } else {
                    // undo the move and erase it from the valid moves
                    cols[h] = prev_j;
                    valid_moves.erase(valid_moves.begin()+rand_index);
                    // undo of the arrays for valuation func
                    start_up--;
                    for (size_t k = prev_i + 1; k <= prev_i + (prev_j - new_j); ++k) {
                        T[k] -= 1;
                        nb_pivot_after_index[k] -= 1;
                    }
                }

                // undo of the arrays for valuation func
                // O(n)
                for (size_t k = 0; k<n-1; k++){
                    if (k <= prev_i) {
                        nb_pivot_after_index_not_moved[k]--;
                    }
                    if (k >= prev_i) {
                        nb_pivot_before_index_not_moved[k]--;
                    }
                }

                H[prev_i] = 0;

                ++loop;
                if (prev_h != h) {
                    prev_h = h;
                }
                // let's try another move
                continue;
            }

            /// set up for next loop turn
            if (prev_j != new_j){
                availablecols[prev_j] = true;
                availablecols[new_j] = false;
            }
            else if (prev_i != new_i){
                availablerows[prev_i] = true;
                availablerows[new_i] = false;
            }
            ++loop;
            ++h;
        }
        FFLAS::fflas_delete (H);
        FFLAS::fflas_delete (T);
        FFLAS::fflas_delete(nb_pivot_before_index);
        FFLAS::fflas_delete(nb_pivot_after_index);
        FFLAS::fflas_delete(nb_pivot_before_index_not_moved);
        FFLAS::fflas_delete(nb_pivot_after_index_not_moved);
    }


    /** @brief  Random Matrix with prescribed rank and rank profile matrix
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param RRP the R dimensional array with row positions of the rank profile matrix' pivots
     * @param CRP the R dimensional array with column positions of the rank profile matrix' pivots
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field,class RandIter>
    inline typename Field::Element_ptr
    RandomMatrixWithRankandRPM (const Field& F, size_t M, size_t N, size_t R,
                                typename Field::Element_ptr A, size_t lda,
                                const size_t * RRP, const size_t * CRP, RandIter& G){

        Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzG(G);

        typename Field::Element_ptr L= FFLAS::fflas_new(F,M,N);

        FFLAS::fzero(F, M, N, L, N);
        // Disabling the  parallel loop, as there is no way to declare G as SHARED in paladin
        //FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads> H;
        //SYNCH_GROUP (FOR1D(k, R, H,
        for (size_t k=0; k<R; ++k){
            size_t i = RRP[k];
            size_t j = CRP[k];
            nzG.random (L [i*N+j]);
            for (size_t l=i+1; l < M; ++l)
                G.random (L [l*N+j]);
        }
        //));

        typename Field::Element_ptr U= FFLAS::fflas_new(F,N,N);
        RandomTriangularMatrix (F, N, N, FFLAS::FflasUpper, FFLAS::FflasNonUnit, true, U, N, G);

        // auto sp=SPLITTER(); //CP: broken with Modular<Integer>. Need to reorganize  the helper behaviour with ParSeq and ModeTraits
        auto sp=NOSPLIT();
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M, N, N,
                      F.one, L, N, U, N, F.zero, A, lda, sp);

        FFLAS::fflas_delete(L);
        FFLAS::fflas_delete(U);
        return A;
    }

    /** @brief  Random Matrix with prescribed rank and rank profile matrix
     * Creates an \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param RRP the R dimensional array with row positions of the rank profile matrix' pivots
     * @param CRP the R dimensional array with column positions of the rank profile matrix' pivots
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomMatrixWithRankandRPM (const Field& F, size_t M, size_t N, size_t R,
                                typename Field::Element_ptr A, size_t lda,
                                const size_t * RRP, const size_t * CRP){
        typename Field::RandIter G(F);
        return RandomMatrixWithRankandRPM (F, M, N, R, A, lda, RRP, CRP, G);
    }

    /** @brief  Random Symmetric Matrix with prescribed rank and rank profile matrix
     * Creates an \c n x \c n symmetric matrix with random entries and rank \c r.
     * @param F field
     * @param n order of \p A
     * @param r rank of \p A
     * @param A the matrix (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param RRP the R dimensional array with row positions of the rank profile matrix' pivots
     * @param CRP the R dimensional array with column positions of the rank profile matrix' pivots
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field,class RandIter>
    inline typename Field::Element_ptr
    RandomSymmetricMatrixWithRankandRPM (const Field& F,  size_t N, size_t R,
                                         typename Field::Element_ptr A, size_t lda,
                                         const size_t * RRP, const size_t * CRP, RandIter& G){

        typename Field::Element_ptr U= FFLAS::fflas_new(F,N,N);
        typename Field::Element_ptr L= FFLAS::fflas_new(F,N,N);
        // U <- $
        RandomTriangularMatrix (F, N, N, FFLAS::FflasUpper, FFLAS::FflasNonUnit, true, U, N, G);
        // L <-  U^T x R
        FFLAS::fzero(F, N, N, L, N);
        for (size_t k=0; k<R; ++k){
            size_t i = RRP[k];
            size_t j = CRP[k];
            FFLAS::fassign (F, N-i, U+i*(N+1), 1, L+j+i*N, N);
        }

        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N,N,N, F.one, L, N, U, N, F.zero, A, lda);

        FFLAS::fflas_delete(L);
        FFLAS::fflas_delete(U);
        return A;
    }

    /** @brief  Random Symmetric Matrix with prescribed rank and rank profile matrix
     * Creates an \c n x \c n symmetric matrix with random entries and rank \c r.
     * @param F field
     * @param n order of \p A
     * @param r rank of \p A
     * @param A the matrix (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param RRP the R dimensional array with row positions of the rank profile matrix' pivots
     * @param CRP the R dimensional array with column positions of the rank profile matrix' pivots
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomSymmetricMatrixWithRankandRPM (const Field& F, size_t M, size_t N, size_t R,
                                         typename Field::Element_ptr A, size_t lda,
                                         const size_t * RRP, const size_t * CRP){
        typename Field::RandIter G(F);
        return RandomSymmetricMatrixWithRankandRPM (F, N, R, A, lda, RRP, CRP, G);
    }

    /** @brief  Random Matrix with prescribed rank, with random rank profile matrix
     * Creates an \c m x \c n matrix with random entries, rank \c r and with a
     * rank profile matrix chosen uniformly at random.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomMatrixWithRankandRandomRPM (const Field& F, size_t M, size_t N, size_t R,
                                      typename Field::Element_ptr A, size_t lda, RandIter& G){
        // generate the r pivots in the rank profile matrix E
        size_t * pivot_r = FFLAS::fflas_new<size_t> (R);
        size_t * pivot_c = FFLAS::fflas_new<size_t> (R);
        RandomRankProfileMatrix (M, N, R, pivot_r, pivot_c);
        RandomMatrixWithRankandRPM (F, M, N, R, A, lda, pivot_r, pivot_c, G);
        FFLAS::fflas_delete(pivot_r);
        FFLAS::fflas_delete(pivot_c);
        return A;
    }

    /** @brief  Random Matrix with prescribed rank, with random  rank profile matrix
     * Creates an \c m x \c n matrix with random entries, rank \c r and with a
     * rank profile matrix chosen uniformly at random.
     * @param F field
     * @param m number of rows in \p A
     * @param n number of cols in \p A
     * @param r rank of the matrix to build
     * @param A the matrix (preallocated to at least \c m x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomMatrixWithRankandRandomRPM (const Field& F, size_t M, size_t N, size_t R,
                                      typename Field::Element_ptr A, size_t lda){
        typename Field::RandIter G(F);
        return RandomMatrixWithRankandRandomRPM (F, M, N, R, A, lda, G);
    }

    /** @brief Random Symmetric Matrix with prescribed rank, with random rank profile matrix
     * Creates an \c n x \c n matrix with random entries, rank \c r and with a
     * rank profile matrix chosen uniformly at random.
     * @param F field
     * @param n order of \p A
     * @param r rank of \p A
     * @param A the matrix (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @param G a random iterator
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomSymmetricMatrixWithRankandRandomRPM (const Field& F, size_t N, size_t R,
                                               typename Field::Element_ptr A, size_t lda, RandIter& G){
        // generate the r pivots in the rank profile matrix E
        size_t * pivot_r = FFLAS::fflas_new<size_t> (R);
        size_t * pivot_c = FFLAS::fflas_new<size_t> (R);
        RandomSymmetricRankProfileMatrix (N, R, pivot_r, pivot_c);
        RandomSymmetricMatrixWithRankandRPM (F, N, R, A, lda, pivot_r, pivot_c, G);
        FFLAS::fflas_delete(pivot_r);
        FFLAS::fflas_delete(pivot_c);
        return A;
    }

    /** @brief Random Symmetric Matrix with prescribed rank, with random rank profile matrix
     * Creates an \c n x \c n matrix with random entries, rank \c r and with a
     * rank profile matrix chosen uniformly at random.
     * @param F field
     * @param n order of \p A
     * @param r rank of \p A
     * @param A the matrix (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomSymmetricMatrixWithRankandRandomRPM (const Field& F, size_t N, size_t R,
                                               typename Field::Element_ptr A, size_t lda){
        typename Field::RandIter G(F);
        return RandomSymmetricMatrixWithRankandRandomRPM (F, N, R, A, lda, G);
    }

    /** @brief  Random Matrix with prescribed det.
     * Creates a \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param d the prescribed value for the determinant of A
     * @param n number of cols in \p A
     * @param A the matrix to be generated (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field>
    inline typename Field::Element_ptr
    RandomMatrixWithDet(const Field & F, size_t n, const typename Field::Element d,
                        typename Field::Element_ptr A, size_t lda) {
        typename Field::RandIter G(F);
        return RandomMatrixWithDet (F, n, d, A, lda, G);
    }
    /** @brief  Random Matrix with prescribed det.
     * Creates a \c m x \c n matrix with random entries and rank \c r.
     * @param F field
     * @param d the prescribed value for the determinant of A
     * @param n number of cols in \p A
     * @param A the matrix to be generated (preallocated to at least \c n x \c lda field elements)
     * @param lda leading dimension of \p A
     * @return \c A.
     */
    template<class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomMatrixWithDet(const Field & F, size_t n, const typename Field::Element d,
                        typename Field::Element_ptr A, size_t lda, RandIter& G){
        FFLASFFPACK_check(n <= lda);
        typedef typename Field::Element  Element ;
        Givaro::GeneralRingNonZeroRandIter<Field,RandIter> nzG (G);

        size_t * P = FFLAS::fflas_new<size_t>(n);
        size_t * Q = FFLAS::fflas_new<size_t>(n);
        for (size_t i = 0 ; i < n ; ++i ) Q[i] = 0;
        for (size_t i = 0 ; i < n ; ++i ) P[i] = 0;

        Element * U = FFLAS::fflas_new<Element>(n*n);
        Element * L = FFLAS::fflas_new<Element>(n*n);

        /*  Create a random P,Q */
        RandomPermutation (n, P);
        RandomPermutation (n, Q);

        /*  det of P,Q */
        int d1 =1 ;
        for (size_t i = 0 ; i < n ; ++i)
            if (P[i] != i)
                d1 = -d1;
        for (size_t i = 0 ; i < n ; ++i)
            if (Q[i] != i)
                d1 = -d1;

        /*  Create L, lower det d */
        RandomTriangularMatrix (F, n, n, FFLAS::FflasLower, FFLAS::FflasNonUnit, true, L, n, G);

        Element dd = F.one;
        for (size_t i=0 ; i<n-1 ; ++i)
            F.mulin(dd,L[i*n+i]);

        F.div(dd,d,dd);
        if (d1<0) F.negin(dd);
        F.assign (L[n*n-1],dd);

        /*  Create U, upper unit*/
        RandomTriangularMatrix (F, n, n, FFLAS::FflasUpper, FFLAS::FflasUnit, true, U, n, G);

        /*  compute product */
        FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
                        n,0,(int)n, U, n, P);
        FFPACK::applyP (F, FFLAS::FflasLeft,  FFLAS::FflasNoTrans,
                        n,0,(int)n, L, n, Q);
        FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                      n,n,n, 1.0, L, n, U, n, 0.0, A, lda);
        // @todo compute LU with ftrtr

        FFLAS::fflas_delete( P);
        FFLAS::fflas_delete( L);
        FFLAS::fflas_delete( U);
        FFLAS::fflas_delete( Q);

        return A;
    }

    template <class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomLTQSMatrixWithRankandQSorder (Field& F, size_t n, size_t r, size_t t,
					typename Field::Element_ptr A, size_t lda, RandIter& G){

        size_t * pivot_r = FFLAS::fflas_new<size_t> (r);
        size_t * pivot_c = FFLAS::fflas_new<size_t> (r);

        RandomLTQSRankProfileMatrix_Tom (n, r, t, pivot_r, pivot_c);
        // typename Field::Element_ptr R =FFLAS::fflas_new(F,n,n);
        // getLTBruhatGen(F, n, r, pivot_r, pivot_c, R, n);
        // FFLAS:: WriteMatrix (std::cerr<<"R = "<<std::endl,F,n,n,R,n);
        // FFLAS::fflas_delete(R);
        // FFLAS::WritePermutation(std::cerr<<"generating a matrix with P=",pivot_r,r);
        // FFLAS::WritePermutation(std::cerr<<"generating a matrix with Q=",pivot_c,r);
    
        RandomMatrixWithRankandRPM (F, n, n, r, A, lda, pivot_r, pivot_c, G);

        FFLAS::fzero (F, FFLAS::FflasRightTri, FFLAS::FflasNonUnit, n, A, lda);

        FFLAS::fflas_delete(pivot_r,pivot_c);
        return A;
    }
} // FFPACK
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
