/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* checkers/Checker_charpoly.inl
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

#ifndef __FFLASFFPACK_checker_charpoly_INL
#define __FFLASFFPACK_checker_charpoly_INL

#ifdef ENABLE_CHECKER_charpoly

class FailureCharpolyCheck {};

template <class Field, class Polynomial> 
class Checker_charpoly {

	const Field& F;
	typename Field::Element lambda;
	typename Field::Element_ptr Ac;
	size_t n;

public:
	Checker_charpoly(const Field& F_, const size_t n_, typename Field::ConstElement_ptr A) 
			: F(F_), Ac(FFLAS::fflas_new(F_,n_,n_)), n(n_)
	{
		//std::cout << "Verifing...";
		typename Field::RandIter G(F);
		//Checker_charpoly(G,n_,A);

		// random lambda
		G.random(lambda);
		//std::cout << "lambda= " << lambda << std::endl;

		// Ac <- A - lambda.I
		for (size_t i=0; i<n; ++i)
			for (size_t j=0; j<n; ++j) {
				if (i==j) F.sub(*(Ac+i*n+j),*(A+i*n+j),lambda);
				else F.assign(*(Ac+i*n+j),*(A+i*n+j));
			}
		//write_field(F,std::cerr<<"Ac=",Ac,n,n,n,true) <<std::endl;

		typename Field::Element_ptr w,	Av = FFLAS::fflas_new(F,n,1),
										v = FFLAS::fflas_new(F,n,1);
		FFLAS::frand(F,G,n,v,1);
		//write_field(F,std::cerr<<"v:=",v,n,1,1,true) <<std::endl;

		// Av <- -A.v
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.mOne, A, n, v, 1, F.zero, Av, 1);

		// Av <- lambda.v + w
		FFLAS::faxpy(F, n, lambda, v, 1, Av, 1);

		// w <- PLUQ.v
		Checker_PLUQ<Field> checker (F,Ac,n,n);
		size_t *P = FFLAS::fflas_new<size_t>(n);
		size_t *Q = FFLAS::fflas_new<size_t>(n);
		size_t r = FFPACK::PLUQ(F, FFLAS::FflasNonUnit, n, n, Ac, n, P, Q);
		//write_field(F,std::cerr<<"L,U:=",Ac,n,n,n,true) <<std::endl;
		w = checker.computePLUQ(Ac,r,P,Q,v);
		//write_field(F,std::cerr<<"w:=",w,n,1,1,true) <<std::endl;

		// w <- Av + w
		FFLAS::faxpy(F, n, F.one, Av, 1, w, 1);
		//write_field(F,std::cerr<<"w:=",w,n,1,1,true) <<std::endl;

		// is w == 0 ?
		bool pass = FFLAS::fiszero(F,n,1,w,1);
		if (!pass) throw FailureCharpolyCheck();
	}

	Checker_charpoly(typename Field::RandIter &G, const size_t n_, typename Field::ConstElement_ptr A)
			: F(G.ring()), Ac(FFLAS::fflas_new(F,n,n)), n(n_)
	{

	}

	~Checker_charpoly() {
	}

	inline bool check(typename Field::ConstElement_ptr A) {

	}
};

#endif // ENABLE_CHECKER_charpoly

#endif // __FFLASFFPACK_checker_charpoly_INL
