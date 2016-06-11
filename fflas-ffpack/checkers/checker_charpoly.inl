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
	const size_t n;
	typename Field::Element lambda, det;

public:
	Checker_charpoly(const Field& F_, const size_t n_, typename Field::Element_ptr A) 
			: F(F_), n(n_)
	{
		typename Field::RandIter G(F);
		init(G,A);
	}

	Checker_charpoly(typename Field::RandIter &G, const size_t n_, typename Field::Element_ptr A)
			: F(G.ring()), n(n_)
	{
		init(G,A);
	}

	~Checker_charpoly() {
	}

	inline bool check(Polynomial &g) {
		//std::cout << "det= "; F.write(std::cout,det); std::cout << std::endl;

		typename Field::Element h = F.zero,
								t = F.one,
								u;
		for (size_t i=0; i < g.size(); ++i) {
			F.mul(u,g[i],t);
			F.add(h,h,u);
			F.mul(t,t,lambda);
		}
		//std::cout << "h= "; F.write(std::cout,h); std::cout << std::endl;

		// is h == det ?
		bool pass = F.areEqual(h,det);
		if (!pass) throw FailureCharpolyCheck();

		return pass;
	}

	//inline typename Field:Element eval_poly(Polynomial &g, typename Field::Element lambda) {

	//}

	// TODO: optimize this function
	inline void pluq(typename Field::Element_ptr A, size_t *P, size_t *Q, size_t R) {
		typename Field::Element_ptr L,U;
		L = FFLAS::fflas_new(F,n,R);
		U = FFLAS::fflas_new(F,R,n);

		for (size_t  i=0; i<R; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign ( *(U + i*n + j), F.zero);
			for (size_t j=i; j<n; ++j)
				F.assign (*(U + i*n + j), *(A+ i*n+j));
		}
		for ( size_t j=0; j<R; ++j ){
			for (size_t i=0; i<=j; ++i )
				F.assign( *(L+i*R+j), F.zero);
			F.assign(*(L+j*R+j), F.one);
			for (size_t i=j+1; i<n; i++)
				F.assign( *(L + i*R+j), *(A+i*n+j));
		}

		FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasTrans, R,0,n, L, R, P);
		FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, R,0,n, U, n, Q);
		FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n,n,R,1.0, L,R, U,n, 0.0, A,n);

		FFLAS::fflas_delete(L,U);
	}

private:
	inline void init(typename Field::RandIter &G, typename Field::Element_ptr A) {
		// random lambda
		G.random(lambda);
		//std::cout << "lambda= " << lambda << std::endl;

		typename Field::Element_ptr v = FFLAS::fflas_new(F,n,1),
									w = FFLAS::fflas_new(F,n,1);
		FFLAS::frand(F,G,n,v,1);
		//write_field(F,std::cerr<<"v:=",v,n,1,1,true) <<std::endl;

		// w <- -A.v
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.mOne, A, n, v, 1, F.zero, w, 1);

		if (!F.isZero(lambda)) {
			// w <- lambda.v + w
			FFLAS::faxpy(F, n, lambda, v, 1, w, 1);

			// A <- A - lambda.I
			for (size_t i=0; i<n; ++i)
				F.sub(*(A+i*n+i),*(A+i*n+i),lambda);
		}

		// P,A,Q <- PLUQ(A)
		size_t *P = FFLAS::fflas_new<size_t>(n);
		size_t *Q = FFLAS::fflas_new<size_t>(n);
		size_t R = FFPACK::PLUQ(F, FFLAS::FflasNonUnit, n, n, A, n, P, Q);
		//std::cout << "rang= " << R << std::endl;
		//write_perm(std::cout<<"P= ",P,n);
		//write_perm(std::cout<<"Q= ",Q,n);

		// compute the determinant of A
		F.init(det,*A);
		for (size_t i=1; i<n; ++i)
			F.mul(det,det,*(A+i*n+i));
		if (n%2 == 1) F.neg(det,det);

		// count the number of permutations
		int t = 0;
		for (size_t i=0; i<n; ++i) {
			if (P[i] != i) t++;
			if (Q[i] != i) t++;
		}
		if (t%2 == 1) F.neg(det,det);
		//std::cout << "det= "; F.write(std::cout,det); std::cout << std::endl;

		// A <- P.L.U.Q (get the initial A-lambda.I)
		pluq(A,P,Q,R);

		// w <- A.v + w
		FFLAS::fgemv(F, FFLAS::FflasNoTrans, n, n, F.one, A, n, v, 1, F.one, w, 1);

		// A <- A + lambda.I (restore A to initial matrix)
		if (!F.isZero(lambda))
			for (size_t i=0; i<n; ++i)
					F.add(*(A+i*n+i),*(A+i*n+i),lambda);

		// is w == 0 ?
		bool pass = FFLAS::fiszero(F,n,1,w,1);

		FFLAS::fflas_delete(v,w);

		if (!pass) throw FailureCharpolyCheck();
	}
};

#endif // ENABLE_CHECKER_charpoly

#endif // __FFLASFFPACK_checker_charpoly_INL
