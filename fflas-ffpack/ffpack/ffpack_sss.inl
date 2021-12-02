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
 * @brief Compute the product of a quasi-separable matrix A, represented by a sequentially semi-separable generator, 
 * with a dense rectangular matrix B:  \f$ C \gets A \times B + beta C \f$
 *
 * @param F the base field
 * @param N the order of \p A
 * @param s the order of quasiseparability of \p A
 * @param D an \f$ N \times s\f$ dense matrix
 * @param ldp leading dimension of \p D
 * @param P an \f$ (N - s) \times s\f$ dense matrix
 * @param ldp leading dimension of \p P
 * @param R an \f$ (N - s - ls) \times s\f$ dense matrix
 * @param ldr leading dimension of \p R
 * @param Q an \f$ (N - ls) \times s\f$ dense matrix
 * @param ldq leading dimension of \p Q
 * @param U an \f$ (N - ls) \times s\f$ dense matrix
 * @param ldu leading dimension of \p U
 * @param V an \f$ (N - ls) \times s\f$ dense matrix
 * @param ldv leading dimension of \p V
 * @param W an \f$ (N - s - ls) \times s\f$ dense matrix
 * @param ldw leading dimension of \p W
 * @param t the number of columns of \p B
 * @param B an \f$ N \times t\f$ dense matrix
 * @param ldb leading dimension of \p B
 * @param beta scaling constant
 * @param [inout] C output matrix
 * @param ldc leading dimension of \p C
 *
 * @bib Missing
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
			    typename  Field::Element_ptr B, size_t ldb,
			    const typename Field::Element beta,
			    typename Field::Element_ptr C, size_t ldc){
  /* Prototype copied from ffpack_bruhatgen.inl */
  
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

  /* Lower */
  typename Field::Element_ptr Temp1 = FFLAS::fflas_new(Fi, s, s);
  typename Field::Element_ptr Temp2 = FFLAS::fflas_new(Fi, s, s);

  /* Temp1 <- alpha * Q_1 * B_1 */
  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
	 s, t, s, alpha, Q,          
	 ldq, B, ldb, Fi.zero,       // No addition
	 Temp1, s);                                    // Leading dimension s
  
  /* Instructions are doubled in the loop in order to avoid using more than two temporary blocks */
  for (size_t block = 0; (block + 2) < kf; block+=2)
    {
      /* C_{block + 1} += P_{block + 1} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
	     s, t, s, Fi.one, P + (block) * s * ldp,   //s x s by s x t, alpha VB 
	     ldp, Temp1, s, Fi.one,                        // No mult
	     C + (block + 1) * s * ldc, ldc);
      /* Temp2 <- alpha * Q_{block + 1} * B_{block + 1} */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     s, t, s, alpha, Q + (block + 1) * s * ldq,          
	     ldq, B + (block + 1) * s * ldb, ldb, Fi.zero,       // No addition
	     Temp2, s);                                    // Leading dimension s
      /* Temp2 += R_{block + 1} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     s, t, s, Fi.one, R + (block) * s * ldr,          
	     ldr, Temp1, s, Fi.one,    
	     Temp2, s);                                    // Leading dimension s

      /* C_{block + 2} += P_{block + 2} * Temp2 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     s, t, s, Fi.one, P + (block + 1) * s * ldp,   
	     ldp, Temp2, s, Fi.one,                        
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
		 Temp1, s);                                    // Leading dimension s
	  /* Temp1 += R_{block + 2} * Temp2 /!\ */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		 s, t, s, Fi.one, R + (block + 1) * s * ldr,          
		 ldr, Temp2, s, Fi.one,    
		 Temp1, s);                                    // Leading dimension s
	}
    }

  /* Remaining blocks */
  if (kf%2 == 0) // Even amount of full blocks (why?)
    {
       /* C_{kf - 1} += P_{kf - 1} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
	     s, t, s, Fi.one, P + (kf - 2) * s * ldp,   //s x s by s x t, alpha VB 
	     ldp, Temp1, s, Fi.one,                        // No mult
	     C + (kf - 1) * s * ldc, ldc);
      if (rs) // Last small block
	{
      /* Temp2 <- alpha * Q_{kf - 1} * B_{kf - 1} */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     s, t, s, alpha, Q + (kf - 1) * s * ldq,          
	     ldq, B + (kf - 1) * s * ldb, ldb, Fi.zero,       // No addition
	     Temp2, s);                                    // Leading dimension s
      /* Temp2 += R_{kf - 1} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     s, t, s, Fi.one, R + (kf - 2) * s * ldr,          
	     ldr, Temp1, s, Fi.one,    
	     Temp2, s);                                    // Leading dimension s

	  /* C_{kf} += P_{kf} * Temp2 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     rs, t, s, Fi.one, P + (kf - 1) * s * ldp,   
	     ldp, Temp2, s, Fi.one,                        
	     C + (kf) * s * ldc, ldc);
	}
    }
  
  else
    if (rs)
      {
	/* C_{kf} += P_{kf} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
	     rs, t, s, Fi.one, P + (kf - 1) * s * ldp,   //s x s by s x t, alpha VB 
	     ldp, Temp1, s, Fi.one,                        // No mult
	     C + (kf) * s * ldc, ldc);
      }

  /* Upper */
  /* Copy, but the other way: partial blocks are taken care of first*/

  /* Temp1 <- alpha * V_lastv * B_lastb */
  if (rs)
    {
    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, rs, alpha,// s x rs by rs x t
	   V + (k - 2) * s * ldv, ldv, B + (k - 1) * s * ldb, ldb, Fi.zero, 
	 Temp1, s);    // Last block is full in height 
    }
  else
    {
    fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, alpha,// s x s by s x t
	   V + (k - 2) * s * ldv, ldv, B + (k - 1) * s * ldb, ldb, Fi.zero,    
	 Temp1, s); // Last block is same as above 
    }
  
  /* Sarting block: n - ls - s */
  for (size_t block = 0; (block + 2) < k; block+=2) // Block increases but index decreases
    {
      /* C_{last - 1 - Block} += U_{last' - block} * Temp1 
       * last' = last - 1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, // s x s by s x t
	     Fi.one, U + (k - 2 - block) * s * ldu, ldu, Temp1, s,
	     Fi.one, C + (k - 2 - block) * s * ldc, ldc);  
      /* Temp2 <- alpha * V_{last - 1 - block} * B_{last - 1 - block} */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     s, t, s, alpha, V + (k - 3 - block) * s * ldv,          
	     ldv, B + (k - 2 - block) * s * ldb, ldb, Fi.zero,    
	     Temp2, s);                                    // Leading dimension s
      /* Temp2 += W_{last - 1 - block} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
	     s, t, s, Fi.one, W + (k - 3 - block) * s * ldw,          
	     ldw, Temp1, s, Fi.one, Temp2, s);                     // Leading dimension s

      /* C_{last - 2 - Block} += U_{last' - 1 - block} * Temp2
       * last' = last - 1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, // s x s by s x t
	     Fi.one, U + (k - 3 - block) * s * ldu, ldu, Temp2, s,
	     Fi.one, C + (k - 3 - block) * s * ldc, ldc);  
      /* Only needs to be done if the results is useful :
       * - Not last loop
       * - Odd amount of U-blocks: C was increased an even number of times at this point
       * Partial blocks are taken care of at the beginning */
      /* Check first condition depending on 'for' loop */
      if (((block + 4) < k)||((k - 1)%2)) // next loop condition passes  -  (k - 1) U-blocks
	{
	  /* Temp1 <- alpha * V_{last - 2 - block} * B_{last - 2 - block} */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		 s, t, s, alpha, V + (k - 4 - block) * s * ldv,          
		 ldv, B + (k - 3 - block) * s * ldb, ldb, Fi.zero,       // No addition
		 Temp1, s);                                    // Leading dimension s
	  /* Temp1 += W_{last - 2 - block} * Temp1 */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, 
		 s, t, s, Fi.one, W + (k - 4 - block) * s * ldw,          
		 ldw, Temp2, s, Fi.one, Temp1, s);                    // Leading dimension s
	}
    }
  if ((k - 1)%2) // Odd amount of U-blocks
    /* C_{1} += U_{1} * Temp1 */
      fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, s, t, s, // s x s by s x t
	     Fi.one, U, ldu, Temp1, s, Fi.one, C, ldc);
  
  /* Free memory */
  FFLAS::fflas_delete(Temp1);
  FFLAS::fflas_delete(Temp2);
}
}

#endif
