/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack/ffpack_pluq_check.inl
 * Copyright (C) 2016 Ashley Lesdalons
 *
 * Written by Ashley Lesdalons <Ashley.Lesdalons@e.ujf-grenoble.fr>
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

#ifndef __FFLASFFPACK_ffpack_pluq_check_INL
#define __FFLASFFPACK_ffpack_pluq_check_INL


#if  defined(DEBUG) || defined(PLUQ_check)
class FailurePLUQcheck {};
#define PLUQ_check_init(Fi,A,M,N) 					\
 	typename Field::Element_ptr v,w; 				\
	v = FFLAS::fflas_new(Fi,N,1); 					\
	w = FFLAS::fflas_new(Fi,N,1); 					\
	FFPACK::init_check_pluq(Fi,A,M,N,v,w);
#define PLUQ_check(Fi,M,N,R,P,A,Q,v,w) 				\
	bool p = FFPACK::check_pluq(Fi,M,N,R,P,A,Q,v,w); \
	FFLAS::fflas_delete(v); 						\
	FFLAS::fflas_delete(w); 						\
	if (!p) throw FailurePLUQcheck();
#endif


namespace FFPACK {


	/** check if a matrix is a zero matrix
	 * @param F
     * @param A
     * @param m
     * @param n
     */
    template <class Field>
    inline bool check_zero_matrix(const Field &F, typename Field::Element_ptr A, size_t m, size_t n) {
    	bool p = true;
    	for (size_t i=0; i < m*n; ++i)
    		p &= F.isZero(A[i]);

    	return p;
    }

	/** check if the PLUQ factorization is correct.
	 *  Returns true if w - P(L(U(Q.v))) == 0
	 * @param F
	 * @param m
	 * @param n
	 * @param r
     * @param A
     * @param P
     * @param Q
     * @param v
     * @param w
     */
    template <class Field>
	inline bool check_pluq(const Field &F,
					size_t m, size_t n, size_t r,
					size_t *P, typename Field::Element_ptr A, size_t *Q,
					typename Field::Element_ptr v, typename Field::Element_ptr w) {

		// v1 <-- Q.v
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, n, v, 1, Q);

		// v1 <-- U.v
		FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, n, 1, F.one, A, n, v, 1);

		// v1 <-- L.v
		FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, n, 1, F.one, A, n, v, 1);

		// v1 <-- P.v
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, m, v, 1, P);

		// is v == w ?
		FFLAS::fsub(F, n, 1, w, 1, v, 1, v, 1);
		bool pass = check_zero_matrix(F, v, n, 1);
		
		return pass;
	} // check_pluq

	/** initialize v with random values and compute w <- A.v
	 * @param F
	 * @param A
	 * @param m
     * @param n
     * @param v unitialized vector of size n
     * @param w unitialized vector of size n
     */
	template <class Field>
	inline void init_check_pluq(const Field &F,
						 typename Field::Element_ptr A,
						 size_t m, size_t n,
						 typename Field::Element_ptr v, typename Field::Element_ptr w) {
		// v is a random vector
		typename Field::RandIter G(F);
        FFLAS::frand(F,G,n,v,1);

        // w <-- A.v
        FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, A, n, v, 1, F.zero, w, 1);	
    } // init_check_pluq


} // FFPACK


 #endif // __FFLASFFPACK_ffpack_pluq_check_INL