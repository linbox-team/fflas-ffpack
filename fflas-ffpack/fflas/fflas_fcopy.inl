/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_fcopy.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __FFLAFLAS_fcopy_INL
#define __FFLAFLAS_fcopy_INL

namespace FFLAS {

	template<class Field>
	inline void
	fcopy (const Field& F, const size_t N,
	       typename Field::Element * X, const size_t incX,
	       const typename Field::Element * Y, const size_t incY )
	{

		if (incY == 1 && incY == 1) {
			memcpy(X,Y,N*sizeof(typename Field::Element)); // much faster (hopefully)
			return;
		}
		typename Field::Element * Xi = X;
		const typename Field::Element * Yi=Y;
		for (; Xi < X+N*incX; Xi+=incX, Yi+=incY )
			F.assign(*Xi,*Yi);
		return;
	}

	template<class Field>
	void fcopy (const Field& F, const size_t m, const size_t n,
		    typename Field::Element * A, const size_t lda,
		    const typename Field::Element * B, const size_t ldb )
	{
		// if possible, copy one big block
		if (lda == n && ldb == n) {
			fcopy(F,m*n,A,1,B,1);
			return ;
		}
		// else, copy row after row
		for (size_t i = 0 ; i < m ; ++i) {
			fcopy(F,n,A+i*lda,1,B+i*ldb,1);
		}
		return;

	}


}


#endif // __FFLAFLAS_fcopy_INL
