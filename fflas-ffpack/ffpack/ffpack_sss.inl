/* ffpack/ffpack_sss.inl
 * Copyright (C) 2021 Hippolyte Signargout
 *
 * Written by Hippolyte Signargout <hippolyte.signargout@ens-lyon.fr>
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

#ifndef __FFLASFFPACK_ffpack_sss_inl
#define __FFLASFFPACK_ffpack_sss_inl

namespace FFPACK{
    template<class Field>
    inline void productSSSxTS (const Field& Fi, size_t N, size_t t, size_t s,
                               const typename Field::Element alpha,
                               typename Field::ConstElement_ptr P, size_t ldp,
                               typename Field::ConstElement_ptr Q, size_t ldq,
                               typename Field::ConstElement_ptr R, size_t ldr,
                               typename Field::ConstElement_ptr U, size_t ldu,
                               typename Field::ConstElement_ptr V, size_t ldv,
                               typename Field::ConstElement_ptr W, size_t ldw,
                               typename Field::ConstElement_ptr D, size_t ldd,
                               typename Field::ConstElement_ptr B, size_t ldb,
                               const typename Field::Element beta,
                               typename Field::Element_ptr C, size_t ldc)
    {
        /*        +--------+------+------+--------+---  +--+         +--+
         *        |   D1   | U1V2 |U1W2V3|U1W2W3V4|     |B1|         |C1|
         *        +--------+------+------+--------+---  +--+         +--+
         *        |  P2Q1  |  D2  | U2V3 | U2W3V4 |     |B2|         |C2|
         * alpha  +--------+------+------+--------+---  +--+  + beta +--+
         *        | P3R2Q1 | P3Q2 |  D3  |  U3V4  |     |B3|         |C3|
         *        +--------+------+------+--------+---  +--+         +--+
         *        |P4R3R2Q1|P4R3Q2| P4Q3 |   D4   |     |B4|         |C4|
         *        +--------+------+------+--------+---  +--+         +--+
         *        |        |      |      |        |     |  |         |  | */
         
        /* Block division */
        size_t kf = N/s;  // Nb of full slices of dimension s
        size_t rs = N%s; //size of the partial block
        size_t k = rs ? kf+1 : kf; // Total number of blocks
        size_t ls = (rs)? rs: s;   // Size of the last block:

        /* Correspondence between matrix and table indices
         * First block   Last block
         * D -> D_1,     D + ((n - ls) * ldd)      -> D_k
         * P -> P_2,     P + ((n - s - ls) * ldp)  -> P_k
         * Q -> Q_1,     Q + ((n - s - ls) * ldq)  -> Q_{k-1}
         * R -> R_2,     R + ((n - 2s - ls) * ldr) -> R_{k-1}
         * U -> U_1,     U + ((n - s - ls) * ldu)  -> U_{k-1}
         * V -> V_2,     V + ((n - s - ls) * ldv)  -> V_k
         * W -> W_2,     W + ((n - 2s - ls) * ldw) -> W_{k-1}
         */

        /* Unused space when s does not divide n:
         *
         * |       |  |       |
         * +---+---+  +---+---+            Code readability is preferred to efficiency:
         * | D | * |  |   |   |        One could store D_last and V_last in new variables
         * +---+---+  | V | * |
         *            |   |   |
         *            +---+---+ */
  
        /************** Diagonal Blocks ************
         * C <- beta * C + alpha D B */
        for (size_t block = 0; block < kf; block++) // Full blocks
            fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
                   s, t, s, alpha, D + block * s * ldd, //s x s by s x t, alpha DB + beta C
                   ldd, B + block * s * ldb, ldb, beta, // D_block is block * s rows under D
                   C + block * s * ldc, ldc);
        if (rs) // Last block
            fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
                   rs, t, rs, alpha, D + kf * s * ldd, //rs x rs by rs x t, alpha DB + beta C
                   ldd, B + kf * s * ldb, ldb, beta, // D_kf is kf * s rows under D
                   C + kf * s * ldc, ldc);

        /************ Lower Triangular Part ***********/
        typename Field::Element_ptr Temp1 = FFLAS::fflas_new(Fi, s, t);
        typename Field::Element_ptr Temp2 = FFLAS::fflas_new(Fi, s, t);
        if (k > 1)
            /* Temp1 <- alpha * Q_1 * B_1 */
            fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                   s, t, s, alpha, Q, ldq, B, ldb, Fi.zero, Temp1, t);
        /* unrolling by step of 2 to avoid swapping temporaries */
        size_t bs = s;
        for (size_t block = 0; (block) < kf; block+=2)
            {
                if ((block + 2) == k)
                    bs = ls;
                if (((block + 2) < kf) || ((kf % 2) == 0) || rs)
                    {
                        // Step a
                        /* C_{block + 2} += P_{block + 2} * Temp1 */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               bs, t, s, Fi.one, P + (block) * s * ldp,
                               ldp, Temp1, t, Fi.one, C + (block + 1) * s * ldc, ldc);
                    }
                if (((block + 2) < kf) || ((rs) && (kf%2 == 0)))
                    {
                        /* Temp2 <- alpha * Q_{block + 2} * B_{block + 2} */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               s, t, s, alpha, Q + (block + 1) * s * ldq, ldq,
                               B + (block + 1) * s * ldb, ldb, Fi.zero, Temp2, t);
                        /* Temp2 += R_{block + 2} * Temp1 */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               s, t, s, Fi.one, R + (block) * s * ldr,
                               ldr, Temp1, t, Fi.one, Temp2, t);
                        if ((block + 3) == k)
                            bs = ls;
                        // Step b
                        /* C_{block + 3} += P_{block + 3} * Temp2 */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               bs, t, s, Fi.one, P + (block + 1) * s * ldp,
                               ldp, Temp2, t, Fi.one,           
                               C + (block + 2) * s * ldc, ldc);
                    }
                /* Only needs to be done if the results is useful :
                 * - Not last loop with instructions
                 */
                if (((block + 2) < kf) && (((block + 4) < kf) || (kf%2 == 0) || rs))
                    {
                        /* Temp1 <- alpha * Q_{block + 3} * B_{block + 3} */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               s, t, s, alpha, Q + (block + 2) * s * ldq,
                               ldq, B + (block + 2) * s * ldb, ldb, Fi.zero,
                               Temp1, t);
                        /* Temp1 += R_{block + 3} * Temp2 */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               s, t, s, Fi.one, R + (block + 1) * s * ldr,
                               ldr, Temp2, t, Fi.one, Temp1, t);
                    }
            }
                /*********** Upper ****************/
                /* Copy, but the other way: partial blocks are taken care of first*/
        
                /* Temp1 <- alpha * V_lastv * B_lastb */
                if (kf > 1 || rs){
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, ls, alpha,
                           V + (k - 2) * s * ldv, ldv, B + (k - 1) * s * ldb, ldb, Fi.zero, Temp1, t);
                }
                /* Starting block: n - ls - s */
                for (size_t block = 0; (block + 2) < k; block+=2) // Block increases but index decreases
                    {
                        /* C_{last - 1 - Block} += U_{last' - block} * Temp1 
                         * last' = last - 1: the last block of U */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s,
                               Fi.one, U + (k - 2 - block) * s * ldu, ldu, Temp1, t,
                               Fi.one, C + (k - 2 - block) * s * ldc, ldc);
                        /* Temp2 <- alpha * V_{last - 1 - block} * B_{last - 1 - block} */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               s, t, s, alpha, V + (k - 3 - block) * s * ldv,
                               ldv, B + (k - 2 - block) * s * ldb, ldb, Fi.zero, Temp2, t);
                        /* Temp2 += W_{last - 1 - block} * Temp1 */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                               s, t, s, Fi.one, W + (k - 3 - block) * s * ldw,
                               ldw, Temp1, t, Fi.one, Temp2, t);

                        /* C_{last - 2 - Block} += U_{last' - 1 - block} * Temp2
                         * last' = last - 1 */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s,
                               Fi.one, U + (k - 3 - block) * s * ldu, ldu, Temp2, t,
                               Fi.one, C + (k - 3 - block) * s * ldc, ldc);
                        /* Only needs to be done if the results is useful :
                         * - Not last loop
                         * - Odd amount of U-blocks: C was increased an even number of times at this point
                         * Partial blocks are taken care of at the beginning */
                        /* Check first condition depending on 'for' loop */
                        if (((block + 4) < k)||((k - 1)%2)) // next loop condition passes / (k - 1) U-blocks
                            {
                                /* Temp1 <- alpha * V_{last - 2 - block} * B_{last - 2 - block} */
                                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                                       s, t, s, alpha, V + (k - 4 - block) * s * ldv,
                                       ldv, B + (k - 3 - block) * s * ldb, ldb, Fi.zero, Temp1, t);
                                /* Temp1 += W_{last - 2 - block} * Temp2 */
                                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                                       s, t, s, Fi.one, W + (k - 4 - block) * s * ldw,
                                       ldw, Temp2, t, Fi.one, Temp1, t);
                            }
                    }
                if ((k - 1)%2) // Odd amount of U-blocks, could be included in the for loop but not necessary
                    /* C_{1} += U_{1} * Temp1 */
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s,
                           Fi.one, U, ldu, Temp1, t, Fi.one, C, ldc);

                FFLAS::fflas_delete(Temp1);
                FFLAS::fflas_delete(Temp2);
            }

        template<class Field>
            inline void SSSToDense (const Field& Fi, size_t N, size_t s,
                                    typename Field::ConstElement_ptr P, size_t ldp,
                                    typename Field::ConstElement_ptr Q, size_t ldq,
                                    typename Field::ConstElement_ptr R, size_t ldr,
                                    typename Field::ConstElement_ptr U, size_t ldu,
                                    typename Field::ConstElement_ptr V, size_t ldv,
                                    typename Field::ConstElement_ptr W, size_t ldw,
                                    typename Field::ConstElement_ptr D, size_t ldd,
                                    typename Field::Element_ptr A, size_t lda)
        {
            /*      +--------+------+------+--------+---
             *      |   D1   | U1V2 |U1W2V3|U1W2W3V4|
             *      +--------+------+------+--------+---
             *      |  P2Q1  |  D2  | U2V3 | U2W3V4 |
             * A =  +--------+------+------+--------+---
             *      | P3R2Q1 | P3Q2 |  D3  |  U3V4  |
             *      +--------+------+------+--------+---
             *      |P4R3R2Q1|P4R3Q2| P4Q3 |   D4   |
             *      +--------+------+------+--------+---
             *      |        |      |      |        |
             */
  
            /* Block division */
        size_t kf = N/s;           // Nb of full slices of dimension s
        size_t rs = N%s;           // Size of the partial block
        size_t k = rs ? kf+1 : kf; // Total number of blocks
        size_t ls = (rs)? rs: s;   // Size of the last block:
        /*   First block    Last block
         * D -> D_1,     D + ((n - ls) * ldd)      -> D_k
         * P -> P_2,     P + ((n - s - ls) * ldp)  -> P_k
         * Q -> Q_1,     Q + ((n - s - ls) * ldq)  -> Q_{k-1}
         * R -> R_2,     R + ((n - 2s - ls) * ldr) -> R_{k-1}
         * U -> U_1,     U + ((n - s - ls) * ldu)  -> U_{k-1}
         * V -> V_2,     V + ((n - s - ls) * ldv)  -> V_k
         * W -> W_2,     W + ((n - 2s - ls) * ldw) -> W_{k-1}
         */

            /* Unused space:
             *
             * |       |  |       |
             * +---+---+  +---+---+            Code readability is preferred to efficiency:
             * | D | * |  |   |   |        One could store D_last and V_last in new variables
             * +---+---+  | V | * |
             *            |   |   |
             *            +---+---+ */

            /******************* Diagonal Blocks **************
             * A <- diag (D) */
            for (size_t block = 0; block < kf; block++) // Full blocks
                /* Diagonal block 'block' starts on row s*block, column s*block */
                FFLAS::fassign (Fi, s, s, D + block * s * ldd, ldd, A + block * s * (lda + 1), lda);
            if (rs) // Last block
                FFLAS::fassign (Fi, rs, rs, D + kf * s * ldd, ldd, A + kf * s * (lda + 1), lda);

        /************** Lower triangular part **********************/
        /* Blocks are computed row by row, by successively applying the R and Q to the P(RRR...) */
        typename Field::Element_ptr Temp1 = FFLAS::fflas_new(Fi, s, s);
        typename Field::Element_ptr Temp2 = FFLAS::fflas_new(Fi, s, s);
        size_t bsize;
        // Loop on rows which have LT part
        for (size_t row = 0; row < k - 1; row++) 
        {
            // In the last iteration, the block may not be full
            bsize = (row == k - 2)? ls: s;
            /* A_{row + 2, row + 1} <- P_{row + 2} * Q_{row + 1} */
            fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                   bsize, s, s, Fi.one, P + row * s * ldp, // In this loop, all blocks are s x s
                   ldp, Q + row * s * ldq, ldq, Fi.zero, A + (row + 1) * s * lda + row * s, lda);
            if (row > 0) /* After row 2, R is also applied 
                          * Next inner loop needs Temp1, which is first set here
                          * Temp1 <- P_{row + 2} * R_{row + 1} */
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bsize, s, s, Fi.one, P + row * s * ldp, 
                       ldp, R + (row - 1) * s * ldr, ldr, Fi.zero, Temp1, s);
            /* unrolling by step of 2 to avoid swapping temporaries */
            for (size_t block = 1; block < row; block+=2)
                {
                    /* A_{row + 2, row + 1 - block} <- Temp1 * Q_{row + 1 - block} */
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bsize, s, s, Fi.one, Temp1, s,
                           Q + (row - block) * s * ldq, ldq, Fi.zero,
                           A + (row + 1) * s * lda + (row - block) * s, lda);
                    /* Temp2 <- Temp1 * R_{row + 1 - block} */
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bsize, s, s, Fi.one, Temp1,
                           s, R + (row - 1 - block) * s * ldr, ldr, Fi.zero, Temp2, s);
                    /* A_{row + 2, row - block} <- Temp2 * Q_{row - block} */
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bsize, s, s, Fi.one, Temp2, s,
                           Q + (row - block - 1) * s * ldq, ldq, Fi.zero,
                           A + (row + 1) * s * lda + (row - block - 1) * s, lda);
                    /* If necessary:
                     * - Not last loop
                     * - One more block
                     * -> At least one more block */
                    if (block + 1 < row)
                        /* Temp1 <- Temp2 * R_{row - block} */
                        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bsize, s, s, Fi.one, Temp2,
                               s, R + (row - block - 2) * s * ldr, ldr, Fi.zero, Temp1, s);
                }
            /* First column if not done already */
            if (row%2)
                /* A_{row + 2, 1} <- Temp1 * Q_{1} */
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, bsize, s, s, Fi.one, Temp1, s,
                       Q, ldq, Fi.zero, A + (row + 1) * s * lda, lda);
          }
       /******************* Upper triangular part *****************/
        /* Symmetrically identical to the lower part (Could be merged with transposes) */
        for (size_t column = 0; column < k - 1; column++) // Loop on columns which have a W
        {
            // In the last iteration, the block may not be full
            bsize = (column == k - 2)? ls: s;
            /* A_{column + 1, column + 2} <- U_{column + 1} * V_{column + 2} */
            fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                   s, bsize, s, Fi.one, U + column * s * ldu, // In this loop, all blocks are s x s
                   ldu, V + column * s * ldv, ldv, Fi.zero, A + (column) * s * lda + (column + 1) * s, lda);
            /* After column 2, R is also applied */
            if (column > 0)
                /* Block loop needs Temp1, which is first updated here */
                /* Temp1 <- W_{column + 1} * V_{column + 2} */
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, bsize, s, Fi.one,
                       W + (column - 1) * s * ldw, ldw, V + column * s * ldv, ldv, Fi.zero, Temp1, s);
            /* Instructions are doubled in the loop in order to avoid using more than two temporary blocks */
            for (size_t block = 1; block < column; block+=2)
            {
                /* A_{column + 1 - block, column + 2} <- U_{column + 1 - block} * Temp1 */
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, bsize, s, Fi.one,
                       U + (column - block) * s * ldu, ldu, Temp1, s, Fi.zero,
                       A + (column - block) * s * lda + (column + 1) * s, lda);
                /* Temp2 <- W_{column + 1 - block} * Temp1 */
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, bsize, s, Fi.one,
                       W + (column - 1 - block) * s * ldw, ldw, Temp1, s, Fi.zero, Temp2, s);
                /* A_{column - block, column + 2} <- U_{column - block} * Temp2 */
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, bsize, s, Fi.one,
                       U + (column - block - 1) * s * ldu,
                       ldu, Temp2, s, Fi.zero, A + (column - block - 1) * s * lda + (column + 1) * s, lda);
                    /* If necessary:
                     * - Not last loop
                     * - One more block
                     * -> At least one more block */
                if (block + 1 < column)
                    /* Temp1 <- W_{column - block} * Temp2 */
                    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, bsize, s, Fi.one,
                           W + (column - block - 2) * s * ldw, ldw, Temp2, s, Fi.zero, Temp1, s);
            }
            /* First row if not done already */
            if (column%2)
                /* A_{1, column + 2} <- U_{1} * Temp1 */
                fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, bsize, s, Fi.one, U,
                       ldu, Temp1, s, Fi.zero, A + (column + 1) * s, lda);
        }    
        FFLAS::fflas_delete(Temp1);
        FFLAS::fflas_delete(Temp2);
    }
        template<class Field>
    inline void DenseToSSS (const Field& Fi, size_t N, size_t s,
                            typename Field::Element_ptr P, size_t ldp,
                            typename Field::Element_ptr Q, size_t ldq,
                            typename Field::Element_ptr R, size_t ldr,
                            typename Field::Element_ptr U, size_t ldu,
                            typename Field::Element_ptr V, size_t ldv,
                            typename Field::Element_ptr W, size_t ldw,
                            typename Field::Element_ptr D, size_t ldd,
                            typename Field::ConstElement_ptr A, size_t lda)
    {
            /*       +--------+------+------+--------+---
             *       |   D1   | U1V2 |U1W2V3|U1W2W3V4|
             *       +--------+------+------+--------+---
             *       |  P2Q1  |  D2  | U2V3 | U2W3V4 |
             * A <-  +--------+------+------+--------+---
             *       | P3R2Q1 | P3Q2 |  D3  |  U3V4  |
             *       +--------+------+------+--------+---
             *       |P4R3R2Q1|P4R3Q2| P4Q3 |   D4   |
             *       +--------+------+------+--------+---
             *       |        |      |      |        |
             */
  
            /* Block division */
        size_t kf = N/s;           // Nb of full slices of dimension s
        size_t rs = N%s;           // Size of the partial block
        size_t k = rs ? kf+1 : kf; // Total number of blocks
        size_t ls = (rs)? rs: s;   // Size of the last block:
        /*   First block    Last block
         * D -> D_1,     D + ((n - ls) * ldd)      -> D_k
         * P -> P_2,     P + ((n - s - ls) * ldp)  -> P_k
         * Q -> Q_1,     Q + ((n - s - ls) * ldq)  -> Q_{k-1}
         * R -> R_2,     R + ((n - 2s - ls) * ldr) -> R_{k-1}
         * U -> U_1,     U + ((n - s - ls) * ldu)  -> U_{k-1}
         * V -> V_2,     V + ((n - s - ls) * ldv)  -> V_k
         * W -> W_2,     W + ((n - 2s - ls) * ldw) -> W_{k-1}
         */

            /* Unused space:
             *
             * |       |  |       |
             * +---+---+  +---+---+            Code readability is preferred to efficiency:
             * | D | * |  |   |   |        One could store D_last and V_last in new variables
             * +---+---+  | V | * |
             *            |   |   |
             *            +---+---+ */

        /******************* Diagonal Blocks **************
             * A <- diag (D) */
        for (size_t block = 0; block < kf; block++) // Full blocks
                /* Diagonal block 'block' starts on row s*block, column s*block */
            FFLAS::fassign (Fi, s, s, A + block * s * (lda + 1), lda, D + block * s * ldd, ldd);
        if (rs) // Last block
            FFLAS::fassign (Fi, rs, rs, A + kf * s * (lda + 1), lda, D + kf * s * ldd, ldd);

        if (N > s) // Otherwise we're done
            {
                size_t fs = s;
                if (N - s < s)
                    {
                        fs = N - s;
                        FFLAS::fzero(Fi, s, s, U, ldu); // U, Q^T will only be updated on their fs first columns
                        FFLAS::fzero(Fi, s, s, Q, ldq);
                    }

                /******************* Upper triangular part *****************/
                // Temporary submatrix, copied to be pluqed
                typename  Field::Element_ptr H = FFLAS::fflas_new (Fi, N, N);
                FFLAS::fassign (Fi, N, N, A, lda, H, N);
                size_t * p = FFLAS::fflas_new<size_t> (N - ls); /* - ls instead of -s so that the size remains 
                                                                   greater than s */
                size_t * q = FFLAS::fflas_new<size_t> (N - ls);
                size_t r = FFPACK::PLUQ (Fi, FFLAS::FflasNonUnit, s, N - s, H + s, N, p, q);

                // pL -> U_1
                FFPACK::getTriangular(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, s, fs, r, H + s,
                                      N, U, ldu);       
                FFPACK::applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans,
                                 fs, 0, s, U, ldu, p);
                // Uq -> V_2
                FFPACK::getTriangular(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit, s, N - s, r,
                                      H + s, N); // Remove L
                FFPACK::applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, s,
                                0, N - s, H + s, N, q); // Apply permutation q to H
                FFLAS::fassign (Fi, s, fs, H + s, N, V, ldv);

                // Temporary matrix for storing lower triangular of pluq
                typename  Field::Element_ptr Temp = FFLAS::fflas_new (Fi, 2 * s, s);
        
                for (size_t brow = 0; brow < k - 2; brow++)
                    {
                        r = FFPACK::PLUQ (Fi, FFLAS::FflasNonUnit, s*(2),
                                          N - s*(brow + 2), H + N * s * brow + s * (brow + 2), N,
                                          p, q);
                        // pL -> [W_{brow + 2} \\ U_{brow + 2}]
                        // 1) L -> Temp
                        FFLAS::fzero (Fi, 2 * s, s, Temp, s);           
                        FFPACK::getTriangular(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, 2*s,
                                              N - s * (brow + 2), r,
                                              H + s * N * brow + s * (brow + 2),
                                              N, Temp, s, true);
                        // 2) p * Temp -> Temp
                        FFPACK::applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, s,
                                        0, 2*s, Temp, s, p);
                        // 3) Temp1 -> W_{brow + 2}
                        FFLAS::fassign (Fi, s, s, Temp, s, W + ldw * s * brow, ldw);
                        // 4) Temp2 -> U_{brow + 2}
                        FFLAS::fassign (Fi, s, s, Temp + s * s, s, U + ldu * s * (brow + 1), ldu);
                        // Uq -> [V_{brow + 3} & H]
                        FFPACK::getTriangular(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit, s, N - s * (brow + 2),
                                              r,
                                              H + s * N * brow + s * (brow + 2), N,
                                              H + s * N * (brow + 1) + s * (brow + 2), N); // Remove L
                        FFPACK::applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, s,
                                        0, N - s* (brow + 2), H + s * N * (brow + 1) + s * (brow + 2), N,
                                        q); // Apply permutation q to H
                        FFLAS::fassign (Fi, s, s, H + s * N * (brow + 1) + s * (brow + 2), N,
                                        V + ldv * s * (brow + 1), ldv); /* Sould not cause any trouble even if
                                                                           last block*/
                    }
                FFLAS::fflas_delete (Temp);
        
                /******************* Lower triangular part *****************/
                // Does it need to be copied? It seems easier than to play with transposes
                r = FFPACK::PLUQ (Fi, FFLAS::FflasNonUnit, N - s, s, H + s * N, N, p, q);
                // Uq -> Q_1 []
                FFPACK::getTriangular(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit, fs, s, r, H + s * N,
                                      N, Q, ldq);       
                FFPACK::applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, fs,
                                0, s, Q, ldq, q);
                // pL -> [P_2 \\ H]
                // Remove L
                FFPACK::getTriangular(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, N - s, s, r,
                                      H + s * N, N);
                // Apply permutation p to H
                FFPACK::applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, s,
                                0, N - s, H + s * N, N, p);
                FFLAS::fassign (Fi, fs, s, H + s * N, N, P, ldp);

                // Temporary matrix for storing upper triangular of pluq
                typename  Field::Element_ptr TempFlat = FFLAS::fflas_new (Fi, s, 2*s);

                for (size_t brow = 0; brow < k - 2; brow++)
                    {
                        r = FFPACK::PLUQ (Fi, FFLAS::FflasNonUnit,
                                          N - s*(brow + 2), s*(2), H + N * s * (brow + 2) + s * (brow), N,
                                          p, q);
                        // Uq -> [R_{brow + 2} & Q_{brow + 2}]
                        // 1) U -> TempFlat
                        FFLAS::fzero (Fi, s, 2 * s, TempFlat, 2 * s);
                        FFPACK::getTriangular(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit,
                                              N - s * (brow + 2), 2*s, r,
                                              H + s * N * (brow + 2) + s * (brow),
                                              N, TempFlat, 2 * s, true);
                        // 2) TempFlat * q -> TempFlat
                        FFPACK::applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, s,
                                        0, 2*s, TempFlat, 2 * s, q);
                        // 3) TempFlat_1 -> R_{brow + 2}
                        FFLAS::fassign (Fi, s, s, TempFlat, 2 * s, R + s * ldr * brow, ldr);
                        // 4) TempFlat_2 -> Q_{brow + 2}
                        FFLAS::fassign (Fi, s, s, TempFlat + s, 2 * s, Q + ldq * s * (brow + 1), ldq);

                        // pL -> [P_{brow + 3} \\ H]
                        // Remove L
                        FFPACK::getTriangular(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, N - s * (brow + 2), s, r,
                                              H + s * N * (brow + 2) + s * (brow), N,
                                              H + s * N * (brow + 2) + s * (brow + 1), N);
                        // Apply permutation p to H
                        FFPACK::applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, s,
                                        0, N - s* (brow + 2), H + s * N * (brow + 2) + s * (brow + 1), N,
                                        p);
                        // H -> P_{brow + 3}
                        FFLAS::fassign (Fi, ((brow == k - 3)? ls: s), s, H + s * N * (brow + 2) + s * (brow + 1), N,
                                        P + ldp * s * (brow + 1), ldp);
                    }
                FFLAS::fflas_delete (H, p, q, TempFlat);
            }
    }
}
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
