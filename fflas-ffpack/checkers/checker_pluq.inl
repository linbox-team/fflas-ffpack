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
	typename Field::Element_ptr v,w,v1,v2;
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
		typename Field::Element_ptr  L, U, v1, v2;
		L = FFLAS::fflas_new(F,m,r);
		U = FFLAS::fflas_new(F,r,n);

		// decompose A into L and U
		for (size_t  i=0; i<r; ++i) {
			for (size_t j=0; j<i; ++j)
				F.assign (*(U + i*n + j), F.zero);
			for (size_t j=i; j<n; ++j)
				F.assign (*(U + i*n + j), *(A+ i*n + j));
		}
		for ( size_t j=0; j<r; ++j) {
			for (size_t i=0; i<=j; ++i)
				F.assign(*(L + i*r + j), F.zero);
			F.assign(*(L + j*r + j), F.one);
			for (size_t i=j+1; i<m; ++i)
				F.assign( *(L + i*r + j), *(A + i*n + j));
		}

		//FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, r,0,m, L, r, P);
		//FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, r,0,n, U, n, Q);

		// v <-- Q.v
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, n, v, 1, Q);

		// v1 <-- U.v
		v1 = FFLAS::fflas_new(F,r,1);
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, r, n, F.one, U, n, v, 1, F.zero, v1, 1);
		//FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, 1, n, F.one, A, r, v, n);

		// v2 <-- L.v1
		v2 = FFLAS::fflas_new(F,m,1);
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, r, F.one, L, r, v1, 1, F.zero, v2, 1);
		//FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, r, 1, F.one, A, m, v, r);

		// v2 <-- P.v2
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, m, v2, 1, P);

		// is v2 == w ?
		FFLAS::fsub(F, m, 1, w, 1, v2, 1, v2, 1);
		bool pass = FFLAS::fiszero(F,n,1,v2,1);

		FFLAS::fflas_delete(L,U,v1,v2);

		//if (!pass) throw FailurePLUQcheck();

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
