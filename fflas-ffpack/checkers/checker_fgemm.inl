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

#ifndef __FFLASFFPACK_checker_fgemm_INL
#define __FFLASFFPACK_checker_fgemm_INL

#ifdef ENABLE_CHECKER_fgemm

class FailureFgemmCheck {};

template <class Field> 
class Checker_fgemm {

	const Field& F;	// add & (BUG)
	const size_t m,n,k,ldc;
	typename Field::Element_ptr v,w1;

public:
	Checker_fgemm(const Field &F_,
	       		  const size_t m_, const size_t n_, const size_t k_,
	       		  const typename Field::Element beta,
	       		  typename Field::Element_ptr C, const size_t ldc_)
		: F(F_), m(m_), n(n_), k(k_), ldc(ldc_), v(FFLAS::fflas_new(F_,n,1)),w1(FFLAS::fflas_new(F_,m,1))
	{			
			//std::cout << "Verifing...";
			typename Field::RandIter G(F);
			init(G,beta,C);
	}

	Checker_fgemm(typename Field::RandIter &G,
	       		  const size_t m_, const size_t n_, const size_t k_,
	       		  const typename Field::Element beta,
	      		  typename Field::Element_ptr C, const size_t ldc_)
		: F(G.field()), m(m_), n(n_), k(k_), ldc(ldc_), v(FFLAS::fflas_new(F,n,1)),w1(FFLAS::fflas_new(F,m,1))
	{
		init(G,beta,C);
	}

	~Checker_fgemm() {
		FFLAS::fflas_delete(v,w1);
	}

	inline bool check(const FFLAS::FFLAS_TRANSPOSE ta,
	    			  const FFLAS::FFLAS_TRANSPOSE tb,
	    			  const typename Field::Element alpha,
	    			  typename Field::ConstElement_ptr A, const size_t lda,
	    			  typename Field::ConstElement_ptr B, const size_t ldb,
	    			  typename Field::ConstElement_ptr C)
	{	
		// w1 <- C.v - w1
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, C, ldc, v, 1, F.mOne, w1, 1);

		// w2 <- B.v
		typename Field::Element_ptr w2 = FFLAS::fflas_new(F,k,1);
		FFLAS::fgemv(F, tb, k, n, F.one, B, ldb, v, 1, F.zero, w2, 1);

		// w1 <- alpha.A.w2 - w1
		FFLAS::fgemv(F, ta, m, k, alpha, A, lda, w2, 1, F.mOne, w1, 1);

		//FFLAS::fflas_delete(w2);

		// is w1 == O ?
		bool pass = FFLAS::fiszero(F, m, w1, 1);
		if (!pass) throw FailureFgemmCheck();
		return pass;
	}

private:
	inline void init(typename Field::RandIter &G, const typename Field::Element beta, typename Field::Element_ptr C) {
		FFLAS::frand(F,G,n,v,1);

		// w1 <- beta.C.v
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, beta, C, ldc, v, 1, F.zero, w1, 1);
	}

};

#endif // ENABLE_CHECKER_fgemm
#endif // __FFLASFFPACK_checker_fgemm_INL
