/* ffpack/ffpack_bruhatgen.inl
 * Copyright (C) 2020 Quentin Houssier
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
 * @brief Compute the product of a lower-triangular quasi-separable matrix A, represented by a sequentially semi-separable generator, 
 * with a dense rectangular matrix B:  \f$ C \gets A \times B + beta C \f$
 *
 * @param F the base field
 * @param N the order of \p A
 * @param s the order of quasiseparability of \p A
 * @param D an \f$ N \times s\f$ dense matrix
 * @param ldp leading dimension of \p D
 * @param P an \f$ N \times s\f$ dense matrix
 * @param ldp leading dimension of \p P
 * @param R an \f$ N \times s\f$ dense matrix
 * @param ldr leading dimension of \p R
 * @param Q an \f$ N \times s\f$ dense matrix
 * @param ldq leading dimension of \p Q
 * @param U an \f$ N \times s\f$ dense matrix
 * @param ldu leading dimension of \p U
 * @param V an \f$ N \times s\f$ dense matrix
 * @param ldv leading dimension of \p V
 * @param W an \f$ N \times s\f$ dense matrix
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
inline  void productSSSxTS (const Field& Fi, size_t N, size_t s,
			    typename Field::ConstElement_ptr P, size_t ldp,
			    typename Field::ConstElement_ptr Q, size_t ldq,
			    typename Field::ConstElement_ptr R, size_t ldr,
			    typename Field::ConstElement_ptr U, size_t ldu,
			    typename Field::ConstElement_ptr V, size_t ldv,
			    typename Field::ConstElement_ptr W, size_t ldw,
			    typename Field::ConstElement_ptr D, size_t ldd,
			    size_t t, const typename Field::Element alpha,
			    typename  Field::Element_ptr B, size_t ldb,
			    const typename Field::Element beta,
			    typename Field::Element_ptr C, size_t ldc){

  /*     +--------+------+------+--------+---  +--+         +--+
         |   D1   | U1V2 |U1W2V3|U1W2W3V4|     |B1|         |C1|
	 +--------+------+------+--------+---  +--+         +--+
	 |  P2Q1  |  D2  | U2V3 | U2W3V4 |     |B2|         |C2|
  alpha  +--------+------+------+--------+---  +--+  + beta +--+
	 | P3R2Q1 | P3Q2 |  D3  |  U3V4  |     |B3|         |C3|
	 +--------+------+------+--------+---  +--+         +--+
	 |P4R3R4Q1|P4R3Q2| P4Q3 |   D4   |     |B4|         |C4|
	 +--------+------+------+--------+---  +--+         +--+
	 |        |      |      |       |      |  |         |  |
  */

  /* Block division */
      size_t kf = N/s;  // Nb of full slices of dimension s
      size_t rs = N%s; //size of the last block
      size_t k = kf;
      if (rs) k++; // Total number of blocks

      TempU = fflas_new (Fi, N, s); // Initialise a temporary N x s matrix
      TempL = fflas_new (Fi, N, s);

      for (size_t block = 0; block < kf; block++) // Full blocks
	{
	  /* Diagonal Blocks */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 s, s, t, alpha, D + block * s * ldd, //s x s by s x t, alpha DB + beta C
		 ldd, B + block * s * ldb, ldb, beta, // D_block is block * s lines under D
		 C + block * s * ldc, ldc);
	  /* Right */
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 s, s, t, alpha, V + block * s * ldv, //s x s by s x t, alpha VB 
		 ldv, B + block * s * ldb, ldb, Fi.zero, // No addition
		 TempU + block * s * s, s); // Leading dimension s
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 s, s, t, alpha, Q + block * s * ldq, //s x s by s x t, alpha VB 
		 ldq, B + block * s * ldb, ldb, Fi.zero, // No addition
		 TempL + block * s * s, s); // Leading dimension s
	}
      if (rs)
	{
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 rs, rs, t, alpha, D + kf * s * ldd, //rs x rs by s x t, alpha DB + beta C
		 ldd, B + kf * s * ldb, ldb, beta, // D_block is block * s lines under D
		 C + kf * s * ldc, ldc);
	  
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 rs, rs, t, alpha, V + kf * s * ldv, //rs x rs by s x t, alpha DB + beta C
		 ldv, B + kf * s * ldb, ldb, Fi.zero, // D_block is block * s lines under D
		 TempU + kf * s * s, s);
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 rs, rs, t, alpha, Q + kf * s * ldq, //rs x rs by s x t, alpha DB + beta C
		 ldq, B + kf * s * ldb, ldb, Fi.zero, // D_block is block * s lines under D
		 Templ + kf * s * s, s);
	}
      /* At this point:
       * C <- beta * C + alpha D B
       * TempU <- alpha * diag(V) * B
       * TempL <- alpha * diag(Q) * B

       /* R and W */
      for (size_t block = 0; block < kf; block++) // Full blocks
	{
	  // C_{block + 1} += P_{block + 1} * TempL_{block}
	  fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, //Element field, no transpose
		 s, s, t, Fi.one, P + (block + 1) * s * ldp, //s x s by s x t, alpha VB 
		 ldp, TempL + block * s * s, s, Fi.one, // No mult
		 C + (block + 1) * s * ldc, ldc);
	  // TempL_{block} *= R_{block + 1}
	  /* Needs to be written */
}
