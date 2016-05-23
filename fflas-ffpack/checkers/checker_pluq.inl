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


#ifdef ENABLE_CHECKING

class FailurePLUQcheck {};

template <class Field> 
class Checker_PLUQ {

	Field F;
	typename Field::Element_ptr A,v,w;
	size_t m,n;

public:
	Checker_PLUQ(Field F, typename Field::Element_ptr A, size_t m, size_t n): F(F), A(A), m(m), n(n) {
		v = FFLAS::fflas_new(F,n,1);
		w = FFLAS::fflas_new(F,n,1);

		// v is a random vector
		typename Field::RandIter G(F);
    	FFLAS::frand(F,G,n,v,1);

    	// w <-- A.v
    	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, A, n, v, 1, F.zero, w, 1);
	}

	~Checker_PLUQ() {
		FFLAS::fflas_delete(v);
		FFLAS::fflas_delete(w);
	}

	/** check if the PLUQ factorization is correct.
	 *  Returns true if w - P(L(U(Q.v))) == 0
	 * @param r
	 * @param P
	 * @param Q
	 */
	inline bool check(size_t r, size_t *P, size_t *Q) {

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
		bool pass = FFLAS::fiszero(F,n,1,v,1);
	
		if (!pass) throw FailurePLUQcheck();

		return pass;
	}
};

#else

template <class Field> 
class Checker_PLUQ {
public:
	Checker_PLUQ(Field F, typename Field::Element_ptr A, size_t m, size_t n) {}
	~Checker_PLUQ() {}
	inline bool check(size_t r, size_t *P, size_t *Q) { return true; }
};

#endif // ENABLE_CHECKING


 #endif // __FFLASFFPACK_ffpack_pluq_check_INL
