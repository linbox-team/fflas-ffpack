/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* checkers/Checker_ftrsm.inl
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

#ifndef __FFLASFFPACK_checker_ftrsm_INL
#define __FFLASFFPACK_checker_ftrsm_INL

#ifdef ENABLE_CHECKER_ftrsm

class FailureTrsmCheck {};

template <class Field> 
class Checker_ftrsm {

	const Field& F;	// add & (BUG)
	typename Field::Element_ptr v,w;
	size_t m,n,k,ldb;

public:
	Checker_ftrsm(const Field& F_, const size_t m_, const size_t n_, typename Field::ConstElement_ptr B, const size_t ldb_) 
				: F(F_), v(FFLAS::fflas_new(F_,n_,1)), w(FFLAS::fflas_new(F_,m_,1)), m(m_), n(n_), ldb(ldb_)
	{
		std::cout << "Verifing...\n";
		typename Field::RandIter G(F);
		init(G,B);
	}

	Checker_ftrsm(typename Field::RandIter &G, const size_t m_, const size_t n_, typename Field::ConstElement_ptr B, const size_t ldb_)
				: F(G.field()), v(FFLAS::fflas_new(F,n_,1)), w(FFLAS::fflas_new(F,m_,1)), m(m_), n(n_), ldb(ldb_)
	{
		init(G,B);
	}

	~Checker_ftrsm() {
		FFLAS::fflas_delete(v,w);
	}

	inline bool check(const FFLAS::FFLAS_SIDE side,
					  const FFLAS::FFLAS_UPLO uplo,
					  const FFLAS::FFLAS_TRANSPOSE trans,
					  const FFLAS::FFLAS_DIAG diag,
					  const typename Field::Element alpha,
					  typename Field::ConstElement_ptr A, size_t lda,
					  typename Field::ConstElement_ptr X) {
		k = (side==FFLAS::FflasLeft?m:n);
		typename Field::Element_ptr v1 = FFLAS::fflas_new(F,k,1);
		write_field(F,std::cerr<<"X:=",X,m,n,ldb,true) <<std::endl;

		// (Left) v1 <- X.v OR (Right) v1 <- alpha.A.v
		if (side==FFLAS::FflasLeft)
			FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, X, ldb, v, 1, F.zero, v1, 1);
		else
			FFLAS::fgemv(F, trans, k, k, alpha, A, lda, v, 1, F.zero, v1, 1);
		write_field(F,std::cerr<<"v1:=",v1,k,1,1,true) <<std::endl;

		// (Left) w <- alpha.A.v1 - w OR (Right) w <- X.v1 - w
		if (side==FFLAS::FflasLeft)
			FFLAS::fgemv(F, trans, k, k, alpha, A, lda, v1, 1, F.mOne, w, 1);
		else 
			FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, X, ldb, v1, 1, F.mOne, w, 1);
		write_field(F,std::cerr<<"w:=",w,m,1,1,true) <<std::endl;

		//FFLAS::fflas_delete(v1);

		bool pass = FFLAS::fiszero(F,m,1,w,1);
		if (!pass) throw FailureTrsmCheck();
		return pass;
	}

private:	
	inline void init(typename Field::RandIter &G, typename Field::ConstElement_ptr B) {
		FFLAS::frand(F,G,n,v,1);
		write_field(F,std::cerr<<"v:=",v,n,1,1,true) <<std::endl;

		// w <- B.v
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, B, ldb, v, 1, F.zero, w, 1);
		write_field(F,std::cerr<<"w:=",w,m,1,1,true) <<std::endl;
	}
};

#endif // ENABLE_CHECKER_ftrsm

#endif // __FFLASFFPACK_checker_ftrsm_INL
