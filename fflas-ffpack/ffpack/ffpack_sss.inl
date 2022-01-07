/* ffpack/ffpack_bruhatgen.inl
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

  /**
   * @brief Compute the product of a quasi-separable matrix A, 
   *represented by a sequentially semi-separable generator, 
   * with a dense rectangular matrix B:  \f$ C \gets \alpha * A \times B + beta C \f$
   *
   * @param Fi the base field
   * @param N the order of \p A
   * @param t the width of \p B and \p C
   * @param s the order of quasiseparability of \p A
   * @param alpha scaling product constant
   * @param D an \f$ N \times s\f$ dense matrix
   * @param ldd leading dimension of \p D
   * @param P an \f$ (N - s) \times s\f$ dense matrix
   * @param ldp leading dimension of \p P
   * @param Q an \f$ (N - ls) \times s\f$ dense matrix
   * @param ldq leading dimension of \p Q
   * @param R an \f$ (N - s - ls) \times s\f$ dense matrix
   * @param ldr leading dimension of \p R
   * @param U an \f$ (N - ls) \times s\f$ dense matrix
   * @param ldu leading dimension of \p U
   * @param V an \f$ (N - ls) \times s\f$ dense matrix
   * @param ldv leading dimension of \p V
   * @param W an \f$ (N - s - ls) \times s\f$ dense matrix
   * @param ldw leading dimension of \p W
   * @param B an \f$ N \times t\f$ dense matrix
   * @param ldb leading dimension of \p B
   * @param beta scaling constant
   * @param [inout] C output matrix
   * @param ldc leading dimension of \p C
   *
   * @bib S. Chandrasekaran et al. “Fast Stable Solver for Sequentially Semi-separable Linear
   * Systems of Equations”. In : High Performance Computing — HiPC 2002. Sous la dir.
   * de Sartaj Sahni, Viktor K. Prasanna et Uday Shukla. Berlin, Heidelberg : Springer
   * Berlin Heidelberg, 2002, p. 545-554. isbn : 978-3-540-36265-4.
   */
  template<class Field>
  inline  void productSSSxTS (const Field& Fi, size_t N, size_t t, size_t s,
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
    /*     +--------+------+------+--------+---  +--+         +--+
	   |   D1   | U1V2 |U1W2V3|U1W2W3V4|     |B1|         |C1|
	   +--------+------+------+--------+---  +--+         +--+
	   |  P2Q1  |  D2  | U2V3 | U2W3V4 |     |B2|         |C2|
    alpha  +--------+------+------+--------+---  +--+  + beta +--+
	   | P3R2Q1 | P3Q2 |  D3  |  U3V4  |     |B3|         |C3|
	   +--------+------+------+--------+---  +--+         +--+
	   |P4R3R4Q1|P4R3Q2| P4Q3 |   D4   |     |B4|         |C4|
	   +--------+------+------+--------+---  +--+         +--+
	   |        |      |      |        |     |  |         |  | */
  
    /* Correspondence between matrix and table indices
     * 
     *   First            Last
     * D.0 -> D_1,   D.((n - ls) * ldd)  
     * P.0 -> P_2,   P.((n - s - ls) * ldp)  
     * Q.0 -> Q_1,   Q.((n - s - ls) * ldq)  
     * R.0 -> R_1,   R.((n - 2s - ls) * ldr)  
     * U.0 -> U_1,   U.((n - s - ls) * ldu)  
     * V.0 -> V_2,   V.((n - s - ls) * ldv)  
     * W.0 -> W_1,   W.((n - 2s - ls) * ldw)   */

    /* Unused space when s does not divide n:
     *
     * |       |  |       |
     * +---+---+  +---+---+            Code readability is preferred to efficiency:
     * | D | * |  |   |   |        One could store D_last and V_last in new variables
     * +---+---+  | V | * |
     *            |   |   |
     *            +---+---+ */
  
    /* Block division */
    size_t kf = N/s;  // Nb of full slices of dimension s
    size_t rs = N%s; //size of the partial block
    size_t k = kf;
    if (rs) k++; // Total number of blocks
    // In the preceeding, ls is the size of the last block:
    //  size_t ls = (rs)? rs: s; Not used in the code
  
    /* Diagonal Blocks :
     * C <- beta * C + alpha D B */
    for (size_t block = 0; block < kf; block++) // Full blocks
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
	     s, t, s, alpha, D + block * s * ldd, //s x s by s x t, alpha DB + beta C
	     ldd, B + block * s * ldb, ldb, beta, // D_block is block * s lines under D
	     C + block * s * ldc, ldc);

    if (rs) // Last block
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
	     rs, t, rs, alpha, D + kf * s * ldd, //rs x rs by s x t, alpha DB + beta C
	     ldd, B + kf * s * ldb, ldb, beta, // D_block is block * s lines under D
	     C + kf * s * ldc, ldc);

    /* Lower Triangular Part */
    typename Field::Element_ptr Temp1 = FFLAS::fflas_new(Fi, s, t);
    typename Field::Element_ptr Temp2 = FFLAS::fflas_new(Fi, s, t);

    if (k > 1)
      /* Temp1 <- alpha * Q_1 * B_1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
	     s, t, s, alpha, Q,          
	     ldq, B, ldb, Fi.zero,       // No addition
	     Temp1, t);                                    // Leading dimension s
  
    /* Instructions are doubled in the loop 
     * in order to avoid using more than two temporary blocks */
    for (size_t block = 0; (block + 2) < kf; block+=2)
      {
	/* C_{block + 1} += P_{block + 1} * Temp1 */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
	       s, t, s, Fi.one, P + (block) * s * ldp,   //s x s by s x t, alpha VB 
	       ldp, Temp1, t, Fi.one,                        // No mult
	       C + (block + 1) * s * ldc, ldc);
	/* Temp2 <- alpha * Q_{block + 1} * B_{block + 1} */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	       s, t, s, alpha, Q + (block + 1) * s * ldq,          
	       ldq, B + (block + 1) * s * ldb, ldb, Fi.zero,       // No addition
	       Temp2, t);                                    // Leading dimension s
	/* Temp2 += R_{block + 1} * Temp1 */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	       s, t, s, Fi.one, R + (block) * s * ldr,          
	       ldr, Temp1, t, Fi.one,    
	       Temp2, t);                                    // Leading dimension s

	/* C_{block + 2} += P_{block + 2} * Temp2 */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	       s, t, s, Fi.one, P + (block + 1) * s * ldp,   
	       ldp, Temp2, t, Fi.one,                        
	       C + (block + 2) * s * ldc, ldc);
	/* Only needs to be done if the results is useful :
	 * - Not last loop
	 * - Even amount of blocks (why?)
	 * - Partial blocks */
	if (((block + 4) < kf)||((kf%2 == 0)||rs))
	  {
	    /* Temp1 <- alpha * Q_{block + 2} * B_{block + 2} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		   s, t, s, alpha, Q + (block + 2) * s * ldq,          
		   ldq, B + (block + 2) * s * ldb, ldb, Fi.zero,       // No addition
		   Temp1, t);                                    // Leading dimension s
	    /* Temp1 += R_{block + 2} * Temp2 /!\ */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		   s, t, s, Fi.one, R + (block + 1) * s * ldr,          
		   ldr, Temp2, t, Fi.one,    
		   Temp1, t);                                    // Leading dimension s
	  }
      }

    /* Remaining blocks */
    if (kf%2 == 0) // Even amount of full blocks (why?)
      {
	/* C_{kf - 1} += P_{kf - 1} * Temp1 */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
	       s, t, s, Fi.one, P + (kf - 2) * s * ldp,   //s x s by s x t, alpha VB 
	       ldp, Temp1, t, Fi.one,                        // No mult
	       C + (kf - 1) * s * ldc, ldc);
	if (rs) // Last small block
	  {
	    /* Temp2 <- alpha * Q_{kf - 1} * B_{kf - 1} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		   s, t, s, alpha, Q + (kf - 1) * s * ldq,          
		   ldq, B + (kf - 1) * s * ldb, ldb, Fi.zero,       // No addition
		   Temp2, t);                                    // Leading dimension s
	    /* Temp2 += R_{kf - 1} * Temp1 */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		   s, t, s, Fi.one, R + (kf - 2) * s * ldr,          
		   ldr, Temp1, t, Fi.one,    
		   Temp2, t);                                    // Leading dimension s

	    /* C_{kf} += P_{kf} * Temp2 */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		   rs, t, s, Fi.one, P + (kf - 1) * s * ldp,   
		   ldp, Temp2, t, Fi.one,                        
		   C + (kf) * s * ldc, ldc);
	  }
      }
  
    else
      if (rs)
	{
	  /* C_{kf} += P_{kf} * Temp1 */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 rs, t, s, Fi.one, P + (kf - 1) * s * ldp,   //s x s by s x t, alpha VB 
		 ldp, Temp1, t, Fi.one,                        // No mult
		 C + (kf) * s * ldc, ldc);
	}

    /* Upper */
    /* Copy, but the other way: partial blocks are taken care of first*/

    /* Temp1 <- alpha * V_lastv * B_lastb */
    if (rs)
      {
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, rs, // s x rs by rs x t
	        alpha,V + (k - 2) * s * ldv, ldv, B + (k - 1) * s * ldb, ldb, Fi.zero, 
	       Temp1, t);    // Last block is full in height 
      }
    else
      {
	if (kf > 1) //No need if no full anti-diagonal block
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, // s x s by s x t
		 alpha, V + (k - 2) * s * ldv, ldv, B + (k - 1) * s * ldb, ldb, Fi.zero,    
		 Temp1, t); // Last block is same as above 
      }
  
    /* Sarting block: n - ls - s */
    for (size_t block = 0; (block + 2) < k; block+=2) // Block increases but index decreases
      {
	/* C_{last - 1 - Block} += U_{last' - block} * Temp1 
	 * last' = last - 1 */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, // s x s by s x t
	       Fi.one, U + (k - 2 - block) * s * ldu, ldu, Temp1, t,
	       Fi.one, C + (k - 2 - block) * s * ldc, ldc);  
	/* Temp2 <- alpha * V_{last - 1 - block} * B_{last - 1 - block} */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	       s, t, s, alpha, V + (k - 3 - block) * s * ldv,          
	       ldv, B + (k - 2 - block) * s * ldb, ldb, Fi.zero,    
	       Temp2, t);                                    // Leading dimension s
	/* Temp2 += W_{last - 1 - block} * Temp1 */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	       s, t, s, Fi.one, W + (k - 3 - block) * s * ldw,          
	       ldw, Temp1, t, Fi.one, Temp2, t);                     // Leading dimension s

	/* C_{last - 2 - Block} += U_{last' - 1 - block} * Temp2
	 * last' = last - 1 */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, // s x s by s x t
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
		   ldv, B + (k - 3 - block) * s * ldb, ldb, Fi.zero,       // No addition
		   Temp1, t);                                    // Leading dimension s
	    /* Temp1 += W_{last - 2 - block} * Temp1 */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		   s, t, s, Fi.one, W + (k - 4 - block) * s * ldw,          
		   ldw, Temp2, t, Fi.one, Temp1, t);                // Leading dimension s
	  }
      }
    if ((k - 1)%2) // Odd amount of U-blocks
      /* C_{1} += U_{1} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, // s x s by s x t
	     Fi.one, U, ldu, Temp1, t, Fi.one, C, ldc);
  
    /* Free memory */
    FFLAS::fflas_delete(Temp1);
    FFLAS::fflas_delete(Temp2);
  }
  

  /**
   * @brief Computes a quasi-separable matrix A from its SSS generators
   *
   * @param Fi the base field
   * @param N the order of \p A
   * @param s the order of quasiseparability of \p A
   * @param P an \f$ (N - s) \times s\f$ dense matrix
   * @param ldp leading dimension of \p P
   * @param Q an \f$ (N - ls) \times s\f$ dense matrix
   * @param ldq leading dimension of \p Q
   * @param R an \f$ (N - s - ls) \times s\f$ dense matrix
   * @param ldr leading dimension of \p R
   * @param U an \f$ (N - ls) \times s\f$ dense matrix
   * @param ldu leading dimension of \p U
   * @param V an \f$ (N - ls) \times s\f$ dense matrix
   * @param ldv leading dimension of \p V
   * @param W an \f$ (N - s - ls) \times s\f$ dense matrix
   * @param ldw leading dimension of \p W
   * @param D an \f$ N \times s\f$ dense matrix
   * @param ldd leading dimension of \p D
   * @param [inout] A the \f$ N \times N \f$ output matrix
   * @param lda leading dimension of \p A
   *
   * @bib S. Chandrasekaran et al. “Fast Stable Solver for Sequentially Semi-separable Linear
   * Systems of Equations”. In : High Performance Computing — HiPC 2002. Sous la dir.
   * de Sartaj Sahni, Viktor K. Prasanna et Uday Shukla. Berlin, Heidelberg : Springer
   * Berlin Heidelberg, 2002, p. 545-554. isbn : 978-3-540-36265-4.
   */
  template<class Field>
  inline  void SSSToDense (const Field& Fi, size_t N, size_t s,
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
     *      |P4R3R4Q1|P4R3Q2| P4Q3 |   D4   |   
     *      +--------+------+------+--------+---
     *      |        |      |      |        |    */
  
    /*   First            Last
     * D.0 -> D_1,   D.((n - ls) * ldd)  
     * P.0 -> P_2,   P.((n - s - ls) * ldp)  
     * Q.0 -> Q_1,   Q.((n - s - ls) * ldq)  
     * R.0 -> R_1,   R.((n - 2s - ls) * ldr)  
     * U.0 -> U_1,   U.((n - s - ls) * ldu)  
     * V.0 -> V_2,   V.((n - s - ls) * ldv)  
     * W.0 -> W_1,   W.((n - 2s - ls) * ldw)   */

    /* Unused space:
     *
     * |       |  |       |
     * +---+---+  +---+---+            Code readability is preferred to efficiency:
     * | D | * |  |   |   |        One could store D_last and V_last in new variables
     * +---+---+  | V | * |
     *            |   |   |
     *            +---+---+ */
  
    /* Block division */
    size_t kf = N/s;  // Nb of full slices of dimension s
    size_t rs = N%s; //size of the partial block
    size_t k = kf;
    if (rs) k++; // Total number of blocks
    // In the preceeding, ls is the size of the last block:
    //  size_t ls = (rs)? rs: s;
  
    /* Diagonal Blocks :
     * A <- diag (D) */
    for (size_t block = 0; block < kf; block++) // Full blocks
      /* Diagonal block 'block' starts on row s*block, column s*block */
      FFLAS::fassign(Fi, s, s, D + block * s * ldd, ldd, A + block * s * (lda + 1), lda); 
    if (rs) // Last block
      FFLAS::fassign(Fi, rs, rs, D + kf * s * ldd, ldd, A + kf * s * (lda + 1), lda); // Last rs x rs block

    /* Lower triangular part */
    /* Blocks are computed line by line, by successively applying the R and Q to the P(RRR...) */
    typename Field::Element_ptr Temp1 = FFLAS::fflas_new(Fi, s, s);
    typename Field::Element_ptr Temp2 = FFLAS::fflas_new(Fi, s, s);

    /* Line 2 */
    if (kf > 1)
      /* A_{2, 1} <- P_{2} * Q_{1} */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, P, 
	     ldp, Q, ldq, Fi.zero, A +  s * lda, lda);
    
    /* After line 2, R is also applied */
    for (size_t line = 1; line < kf - 1; line++) // Loop on lines which have LT part
      {
	/* Loop needs Temp1, which is first set here */
	/* A_{line + 2, line + 1} <- P_{line + 2} * Q_{line + 1} */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
	       s, s, s, Fi.one, P + line * s * ldp, // In this loop, all blocks are s x s
	       ldp, Q + line * s * ldq, ldq, Fi.zero, A + (line + 1) * s * lda + line * s, lda);
	/* Temp1 <- P_{line + 2} * R_{line + 1} */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, P + line * s * ldp, 
	       ldp, R + (line - 1) * s * ldr, ldr, Fi.zero, Temp1, s);
	
	/* Instructions are doubled in the loop in order to avoid using more than two temporary blocks */
	for (size_t block = 1; block < line; block+=2)
	  {
	    /* A_{line + 2, line + 1 - block} <- Temp1 * Q_{line + 1 - block} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, Temp1, s, 
		   Q + (line - block) * s * ldq, ldq, Fi.zero, A + (line + 1) * s * lda + (line - block) * s, lda);
	    /* Temp2 <- Temp1 * R_{line + 1 - block} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, Temp1, 
		   s, R + (line - 1 - block) * s * ldr, ldr, Fi.zero, Temp2, s);

	    /* A_{line + 2, line - block} <- Temp2 * Q_{line - block} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, Temp2, s,
		   Q + (line - block - 1) * s * ldq, ldq, Fi.zero, A + (line + 1) * s * lda + (line - block - 1) * s, lda);
	    /* If necessary:
	     * - Not last loop
	     * - One more block
	     * -> At least one more block */
	    if (block + 1 < line)
	    /* Temp1 <- Temp2 * R_{line - block} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, Temp2, 
		   s, R + (line - block - 2) * s * ldr, ldr, Fi.zero, Temp1, s);
	  }
	/* First column if not done already */
	if (line%2)
	  /* A_{line + 2, 1} <- Temp1 * Q_{1} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, Temp1, s,
		   Q, ldq, Fi.zero, A + (line + 1) * s * lda, lda);
      }
    
    /* Last line if partial */
    /* rs lines in the block */
    if (rs && (k > 1))
            {
	/* Loop needs Temp1, which is first set here */
	/* A_{kf + 1, kf} <- P_{kf + 1} * Q_{kf} */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
	       rs, s, s, Fi.one, P + (kf - 1) * s * ldp, 
	       ldp, Q + (kf - 1) * s * ldq, ldq, Fi.zero, A + kf * s * lda + (kf - 1) * s, lda);
	if (kf > 1) // Otherwise no R
	  /* Temp1 <- P_{(kf + 1)} * R_{kf} */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rs, s, s, Fi.one, P + (kf - 1) * s * ldp, 
		 ldp, R + (kf - 2) * s * ldr, ldr, Fi.zero, Temp1, s); // rs x s by s x s
	
	/* Instructions are doubled in the loop in order to avoid using more than two temporary blocks */
	for (size_t block = 1; block < (kf - 1); block+=2)
	  {
	    /* A_{kf + 1, kf - block} <- Temp1 * Q_{kf - block} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rs, s, s, Fi.one, Temp1, s, //Temp1 is only up to date on rs lines
		   Q + ((kf - 1) - block) * s * ldq, ldq, Fi.zero, A + (kf) * s * lda + ((kf - 1) - block) * s, lda);
	    /* Temp2 <- Temp1 * R_{kf - block} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rs, s, s, Fi.one, Temp1, 
		   s, R + ((kf - 2) - block) * s * ldr, ldr, Fi.zero, Temp2, s);

	    /* A_{kf + 1, (kf - 1) - block} <- Temp2 * Q_{(kf - 1) - block} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rs, s, s, Fi.one, Temp2, s,
		   Q + ((kf - 2) - block) * s * ldq, ldq, Fi.zero, A + (kf) * s * lda + ((kf - 2) - block) * s, lda);
	    /* If necessary:
	     * - Not last loop
	     * - One more block
	     * -> At least one more block */
	    if (block + 1 < kf - 1)
	    /* Temp1 <- Temp2 * R_{(kf - 1) - block} */
	      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rs, s, s, Fi.one, Temp2, 
		     s, R + ((kf - 3) - block) * s * ldr, ldr, Fi.zero, Temp1, s);
	  }
	/* First column if not done already */
	if ((kf - 1)%2)
	  /* A_{(kf + 1), 1} <- Temp1 * Q_{1} */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, rs, s, s, Fi.one, Temp1, s,
		   Q, ldq, Fi.zero, A + kf * s * lda, lda);
      }

    /* Upper tirangular part */
    /* Symmetrically identical to the lower part */

    /* Column 2 */
    if (kf > 1)
      /* A_{1, 2} <- U_{2} * V_{1} */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, U, 
	     ldu, V, ldv, Fi.zero, A + s, lda);
    
    /* After column 2, R is also applied */
    for (size_t column = 1; column < kf - 1; column++) // Loop on columns which have a W
      {
	/* Block loop needs Temp1, which is first updated here */
	/* A_{column + 1, column + 2} <- U_{column + 1} * V_{column + 2} */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
	       s, s, s, Fi.one, U + column * s * ldu, // In this loop, all blocks are s x s
	       ldu, V + column * s * ldv, ldv, Fi.zero, A + (column) * s * lda + (column + 1) * s, lda);
	/* Temp1 <- W_{column + 1} * V_{column + 2} */
	fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, W + (column - 1) * s * ldw, 
	       ldw, V + column * s * ldv, ldv, Fi.zero, Temp1, s);
	
	/* Instructions are doubled in the loop in order to avoid using more than two temporary blocks */
	for (size_t block = 1; block < column; block+=2)
	  {
	    /* A_{column + 1 - block, column + 2} <- U_{column + 1 - block} * Temp1 */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, U + (column - block) * s * ldu,
		   ldu, Temp1, s, Fi.zero, A + (column - block) * s * lda + (column + 1) * s, lda);
	    /* Temp2 <- W_{column + 1 - block} * Temp1 */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one,
		   W + (column - 1 - block) * s * ldw, ldw, Temp1, s, Fi.zero, Temp2, s);

	    /* A_{column - block, column + 2} <- U_{column - block} * Temp2 */
	    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, U + (column - block - 1) * s * ldu,
		   ldu, Temp2, s, Fi.zero, A + (column - block - 1) * s * lda + (column + 1) * s, lda);
	    /* If necessary:
	     * - Not last loop
	     * - One more block
	     * -> At least one more block */
	    if (block + 1 < column)
	      /* Temp1 <- W_{column - block} * Temp2 */
	      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one,
		     W + (column - block - 2) * s * ldw, ldw, Temp2, s, Fi.zero, Temp1, s);
	  }
	/* First line if not done already */
	if (column%2)
	  /* A_{1, column + 2} <- U_{1} * Temp1 */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, s, s, Fi.one, U,
		 ldu, Temp1, s, Fi.zero, A + (column + 1) * s, lda);
      }
    
    /* Last column if partial 
     * and partial column is not just a diagonal */
    /* rs columns in the block */
    if (rs && (k > 1))
            {
	      /* Block loop needs Temp1, which is first updated here */
	      /* A_{kf, kf + 1} <- U_{kf} * V_{kf + 1} */
	      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
		     s, rs, s, Fi.one, U + (kf - 1) * s * ldu, // In this loop, all blocks are s x rs
		     ldu, V + (kf - 1) * s * ldv, ldv, Fi.zero, A + (kf - 1) * s * lda + (kf) * s, lda);
	      if (kf > 1) // Or W does not exist
		/* Temp1 <- W_{kf} * V_{kf + 1} */
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, rs, s, Fi.one, W + (kf - 2) * s * ldw, 
		       ldw, V + (kf - 1) * s * ldv, ldv, Fi.zero, Temp1, s);
	
	      /* Instructions are doubled in the loop in order to avoid using more than two temporary blocks */
	      for (size_t block = 1; block < kf - 1; block+=2)
		{
		  /* A_{kf - block, kf + 1} <- U_{kf - block} * Temp1 */
		  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, rs, s, Fi.one, U + (kf - 1 - block) * s * ldu,
			 ldu, Temp1, s, Fi.zero, A + (kf - 1 - block) * s * lda + (kf) * s, lda);
		  /* Temp2 <- W_{kf - block} * Temp1 */
		  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, rs, s, Fi.one,
			 W + (kf - 2 - block) * s * ldw, ldw, Temp1, s, Fi.zero, Temp2, s);

		  /* A_{kf - 1 - block, kf + 1} <- U_{kf - 1 - block} * Temp2 */
		  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, rs, s, Fi.one, U + (kf - block - 2) * s * ldu,
			 ldu, Temp2, s, Fi.zero, A + (kf - block - 2) * s * lda + (kf) * s, lda);
		  /* If necessary:
		   * - Not last loop
		   * - One more block
		   * -> At least one more block */
		  if (block + 1 < kf - 1)
		    /* Temp1 <- W_{kf - 1 - block} * Temp2 */
		    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, rs, s, Fi.one,
			   W + (kf - block - 3) * s * ldw, ldw, Temp2, s, Fi.zero, Temp1, s);
		}
	      /* First column if not done already */
	      if ((kf%2) == 0)
		/* A_{1, kf + 1} <- U_{1} * Temp1 */
		fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, rs, s, Fi.one, U,
		       ldu, Temp1, s, Fi.zero, A + (kf) * s, lda);
      }
    
    /* Free memory */
    FFLAS::fflas_delete(Temp1);
    FFLAS::fflas_delete(Temp2);
  }
}

#endif
