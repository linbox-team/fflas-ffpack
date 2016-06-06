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

#ifndef __FFLASFFPACK_checker_pluq_INL
#define __FFLASFFPACK_checker_pluq_INL

#ifdef ENABLE_CHECKER_PLUQ

class FailurePLUQcheck {};

template <class Field> 
class Checker_PLUQ {

	const Field& F;
	typename Field::Element_ptr v,w;
	const size_t m,n;

public:
	Checker_PLUQ(const Field& F_, typename Field::Element_ptr A, size_t m_, size_t n_) 
				: F(F_), v(FFLAS::fflas_new(F_,n_,1)), w(FFLAS::fflas_new(F_,m_,1)), m(m_), n(n_)
	{
		typename Field::RandIter G(F);
		init(G,A);
	}

	Checker_PLUQ(typename Field::RandIter &G, typename Field::Element_ptr A, size_t m_, size_t n_)
				: F(G.field()), v(FFLAS::fflas_new(F,n_,1)), w(FFLAS::fflas_new(F,m_,1)), m(m_), n(n_)
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
		//std::cerr << "check, r: " << r << std::endl;
		//write_perm(std::cerr<<"P:=",P,m)<<std::endl;
		//write_field(F,std::cerr<<"L,U:=",A,m,n,n,true) <<std::endl;
		//write_perm(std::cerr<<"Q:=",Q,n)<<std::endl;

		typename Field::Element_ptr _w = FFLAS::fflas_new(F,m,1); // _w = [w1|w2]

		//write_field(F,std::cerr<<"v0:=",v,n,1,1,true) <<std::endl;

		// v <-- Q.v
		FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, 1, 0, n, v, 1, Q);
		//write_field(F,std::cerr<<"v1:=",v,n,1,1,true) <<std::endl;

		// w1 <- V1 && w2 <- 0
 		FFLAS::fassign(F, r, 1, v, 1, _w, 1);
 		for (size_t i=r; i<m; ++i)
			F.assign(*(_w+i),F.zero);

		//write_field(F,std::cerr<<"w0:=",_w,m,1,1,true) <<std::endl;

 		// w1 <- U1.w1
 		FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, r, 1, F.one, A, n, _w, 1);
		
		//write_field(F,std::cerr<<"w1:=",_w,m,1,1,true) <<std::endl;

 		// w1 <- U2.V2 + w1
 		if (r < n)
 			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, r, 1, n-r, F.one, A+r, n, v+r, 1, F.one, _w, 1);
		//write_field(F,std::cerr<<"w2:=",_w,m,1,1,true) <<std::endl;

 		// w2 <- L2.w1
 		if (r < m)
 			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m-r, 1, r, F.one, A+r*n, n, _w, 1, F.zero, _w+r, 1);
		//write_field(F,std::cerr<<"w3:=",_w,m,1,1,true) <<std::endl;
 		

 		// w1 <- L1.w1
 		FFLAS::ftrmm(F, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, FFLAS::FflasUnit, r, 1, F.one, A, n, _w, 1);
		//write_field(F,std::cerr<<"w4:=",_w,m,1,1,true) <<std::endl;

 		// _w <- P._w
 		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, 1, 0, m, _w, 1, P);
		//write_field(F,std::cerr<<"w5:=",_w,m,1,1,true) <<std::endl;

 		// is _w == w ?
		FFLAS::fsub(F, m, 1, w, 1, _w, 1, _w, 1);

		//write_field(F,std::cerr<<"_w:=",_w,m,1,1,true) <<std::endl;

		bool pass = FFLAS::fiszero(F,m,1,_w,1);

		FFLAS::fflas_delete(_w);

		if (!pass) throw FailurePLUQcheck();

		return pass;
	}

private:	
	inline void init(typename Field::RandIter &G, typename Field::Element_ptr A) {
		FFLAS::frand(F,G,n,v,1);
		//write_field(F,std::cerr<<"v:=",v,n,1,1,true) <<std::endl;
    	
    	// w <-- A.v
    	FFLAS::fgemv(F, FFLAS::FflasNoTrans, m, n, F.one, A, n, v, 1, F.zero, w, 1);
		//write_field(F,std::cerr<<"w:=",w,m,1,1,true) <<std::endl;

	}
};

#endif // ENABLE_CHECKER_PLUQ

#endif // __FFLASFFPACK_checker_pluq_INL
