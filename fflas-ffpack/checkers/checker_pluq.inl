/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* checkers/checker_pluq.inl
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

#ifdef ENABLE_CHECKER_PLUQ

class FailurePLUQcheck {};

template <class Field> 
class Checker_PLUQ {
public:
	const Field F;	// add & (BUG)
	typename Field::Element_ptr v,w;
	const size_t m,n;

//public:
	Checker_PLUQ(Field F, typename Field::Element_ptr A, size_t m, size_t n) 
				: F(F), v(FFLAS::fflas_new(F,n,1)), w(FFLAS::fflas_new(F,m,1)), m(m), n(n)
	{
		// v is a random vector
		typename Field::RandIter G(F);
		init(G,A);
	}

	Checker_PLUQ(typename Field::RandIter &G, typename Field::Element_ptr A, size_t m, size_t n)
				: F(G.field()), v(FFLAS::fflas_new(F,n,1)), w(FFLAS::fflas_new(F,m,1)), m(m), n(n)
	{
		init(G,A);
	}

	~Checker_PLUQ() {
		FFLAS::fflas_delete(v,w);
	}

	/** check if the PLUQ factorization is correct.
	 *  Returns true if w - P(L(U(Q.v))) == 0
	 * @param r
	 * @param P
	 * @param Q
	 */
	inline bool check(typename Field::Element_ptr A, size_t r, size_t *P, size_t *Q) {
		typename Field::Element_ptr R = FFLAS::fflas_new(F,r,1),
									_w = FFLAS::fflas_new(F,m,1);
		for (size_t i=0; i<m; ++i)
			F.assign(*(_w+i),F.zero);

		// v <-- Q.v
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, n, v, 1, Q);

		// R <- V1
 		FFLAS::fassign(F, r, 1, v, 1, R, 1);

 		// R <- U1.R
 		FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, r, 1, F.one, A, n, R, 1);

 		// R <- U2.V2 + R
 		if (r < n)
 			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, r, 1, n-r, F.one, A+r, n, v+r, 1, F.one, R, 1);

 		// w2 <- L2.R
 		if (r < m)
 			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m-r, 1, r, F.one, A+r*n, n, R, 1, F.zero, _w+r, 1);

 		// R <- L1.R
 		FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, r, 1, F.one, A, n, R, 1);

 		// w1 <- R
 		FFLAS::fassign(F, r, 1, R, 1, _w, 1);

 		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, m, _w, 1, P);

 		// is _w == w ?
		FFLAS::fsub(F, m, 1, w, 1, _w, 1, _w, 1);
		bool pass = FFLAS::fiszero(F,m,1,_w,1);

		FFLAS::fflas_delete(R,_w);

		if (!pass) throw FailurePLUQcheck();

		return pass;
	}

private:	
	inline void init(typename Field::RandIter &G, typename Field::Element_ptr A) {
		FFLAS::frand(F,G,n,v,1);
    	// w <-- A.v
    	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, A, n, v, 1, F.zero, w, 1);
	}
};

#endif // ENABLE_CHECKER_PLUQ

#endif // __FFLASFFPACK_ffpack_pluq_check_INL
