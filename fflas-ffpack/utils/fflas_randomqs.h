/*
 * Copyright (C) FFLAS-FFPACK
 * Written by 
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * Hippolyte Signargout <hippolyte.signargout@ens-lyon.fr>
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

#ifndef __FFLASFFPACK_randomqs_H
#define __FFLASFFPACK_randomqs_H

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/debug.h"
#include "fflas-ffpack/fflas/fflas.h"
#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include <givaro/givranditer.h>

namespace FFPACK {

    /** @brief Random RPM for left triangular quasi-separable matrix
     * @param n dimension of the matrix
     * @param r rank of the matrix
     * @param t order of quasiseparability
     * @param[inout] rows row indices of the pivots
     * @param[inout] cols column indices of the pivots
     */
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

    /** @brief Random RPM for left triangular quasi-separable matrix using a
     * recursive strategy
     * @param n dimension of the matrix
     * @param r number of pivots to put in the triangle
     * @param t order of quasiseparability
     * @param[inout] ranks array of ranks of the current (n-i * i) submatrices
     * @param row_index row index of the most topleft entry of the subtriangle
     * to consider
     * @param col_index column index of the most topleft entry of the 
     * subtriangle to consider
     * @param[inout] rows row indices of the pivots
     * @param[inout] cols column indices of the pivots
     * @param pivots number of pivots already set
     */
    inline void RandomRecursiveLTQSRankProfileMatrix (size_t n, size_t r,
                                                      size_t t, size_t *ranks,
                                                      size_t row_index,
                                                      size_t col_index,
                                                      size_t * rows,
                                                      size_t *cols,
                                                      size_t pivots)
    {
        size_t triangle_size = n - row_index - col_index;
        if (triangle_size != 0) // Otherwise do nothing: we consider r = 0
            {
                if (triangle_size == 1) // No division
                    {
                        if (r == 1) // Otherwise r = 0, do nothing
                            {
                                // Add entry to pivots
                                rows[pivots] = row_index;
                                cols[pivots] = col_index;
                                //pivots++ //This does not seem useful
                            }
                    }
                else // n \geq 2, division
                    {
                        // Choose division (could be in a subfunction?)
                        if (triangle_size == r)
                            size_t rc = 1;
                        else
                            {
                                size_t ceil_ts = ceil (triangle_size / 2.)
                                size_t rc_bound = min(r, t, triangle_size
                                                      + 1 - r, ceil_ts);
                                /* Does this work?
                                 * r <= n
                                 */
                                size_t rc = RandInt (0, rc_bound);
                            }
                        size_t rd_lbound = max (0, r - ceil_ts, r - ceil_ts - rc + 1);
                        size_t rd_ubound = min (r - rc, triangle_size/2,
                                                triangle_size/2 - rc + 1);
                        size_t rd = RandInt (rd_lbound, rd_ubound);
                        size_t rb = r - rc - rd;
                        // Recursive calls
                        /* Place pivots in bottom triangular part */
                        RandomRecursiveLTQSRankProfileMatrix(n, rb, t, ranks, row_index +
                                                             triangle_size + 1 - ceil_ts,
                                                             col_index,rows,cols,pivots);
                        pivots += rb;
                        /* Place pivots in  triangular part */
                        RandomRecursiveLTQSRankProfileMatrix(n, rd, t, ranks, row_index,
                                                             col_index + ceil_ts, rows,
                                                             cols, pivots);
                        pivots += rd;
                        // Preliminaries for rectangle
                        ranks[ceil_ts - 1] = t;
                        //findimin
                        size_t imin = 0;
                        for (size_t l =col_index;(imin == 0) && (l < n - row_index); l++)
                                if (ranks[l] == 0)
                                    imin = n - l + 2;
                        //findjmin
                        size_t jmin = 0;
                        for (size_t l =n - row_index;(imin == 0) && (l >= col_index);l--)
                                if (ranks[l] == 0)
                                    imin = l + 1;
                        // Fill rectangle -> Placer
                        for (size_t k = 0; k < rc; k++)
                            {
                                bool correcti = false;
                                while (correcti == false)
                                    {
                                        size_t i = RandInt (imin, triangle_size/2 +1);
                                        correcti = true;
                                        for (size_t l = 0; correcti && (l < pivots); l++)
                                            if (rows[l] == i)
                                                correcti = false;
                                    }
                                bool correctj = false;
                                while (correctj == false)
                                    {
                                        size_t j = RandInt (jmin, ceil_ts);
                                        correctj = true;
                                        for (size_t l = 0; correctj && (l < pivots); l++)
                                            if (rows[l] == j)
                                                correctj = false;
                                    }
                                for(size_t l = col_index + j; l < n -row_index-i+2; l++)
                                    {
                                        ranks[l]--;
                                        //Rest of Placer (might miss a "change_imin" var)
                                    }
                            }
                    }
            }
    }
    
    template <class Field, class RandIter>
    inline typename Field::Element_ptr
    RandomLTQSMatrixWithRankandQSorder (Field& F, size_t n, size_t r, size_t t,
					typename Field::Element_ptr A, size_t lda, RandIter& G){

        size_t * pivot_r = FFLAS::fflas_new<size_t> (r);
        size_t * pivot_c = FFLAS::fflas_new<size_t> (r);

        RandomLTQSRankProfileMatrix (n, r, t, pivot_r, pivot_c);
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
