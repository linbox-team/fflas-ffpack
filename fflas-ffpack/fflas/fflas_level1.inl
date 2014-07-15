/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by BB <bbboyer@ncsu.edu>
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

/** @file fflas/fflas_level1.h
 * @brief  Vector operations
 * or anything of \f$n\f$ complexity
 */

#ifndef __FFLASFFPACK_fflas_fflas_level1_INL
#define __FFLASFFPACK_fflas_fflas_level1_INL

namespace FFLAS {



	//---------------------------------------------------------------------
	// Level 1 routines
	//---------------------------------------------------------------------

	/** finit
	 * \f$x \gets  x mod F\f$.
	 * @param F field
	 * @param n size of the vectors
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 * @bug use cblas_(d)scal when possible
	 */
	template<class Field>
	void
	finit (const Field& F, const size_t n,
	       typename Field::Element * X, const size_t incX);

	/** finit
	 * \f$x \gets  y mod F\f$.
	 * @param F field
	 * @param n size of the vectors
	 * \param Y vector of \p OtherElement
	 * \param incY stride of \p Y
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 * @bug use cblas_(d)scal when possible
	 */
	template<class Field, class OtherElement>
	void
	finit (const Field& F, const size_t n,
	       const OtherElement * Y, const size_t incY,
	       typename Field::Element * X, const size_t incX);

	/** fconvert
	 * \f$x \gets  y mod F\f$.
	 * @param F field
	 * @param n size of the vectors
	 * \param Y vector of \p F
	 * \param incY stride of \p Y
	 * \param X vector in \p OtherElement
	 * \param incX stride of \p X
	 * @bug use cblas_(d)scal when possible
	 */
	template<class Field, class OtherElement>
	void
	fconvert (const Field& F, const size_t n,
		  OtherElement * X, const size_t incX,
		  const typename  Field::Element* Y, const size_t incY)
	{
		OtherElement * Xi = X ;
		const typename Field::Element * Yi = Y ;
		for (; Xi < X+n*incX; Xi+=incX, Yi += incY )
			F.convert( *Xi , *Yi);
	}

	/** fnegin
	 * \f$x \gets - x\f$.
	 * @param F field
	 * @param n size of the vectors
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 * @bug use cblas_(d)scal when possible
	 */
	template<class Field>
	void
	fnegin (const Field& F, const size_t n,
		typename Field::Element * X, const size_t incX)
	{
		typename Field::Element * Xi = X ;
		for (; Xi < X+n*incX; Xi+=incX )
			F.negin( *Xi );
	}

	/** fneg
	 * \f$x \gets - y\f$.
	 * @param F field
	 * @param n size of the vectors
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 * \param Y vector in \p F
	 * \param incY stride of \p Y
	 * @bug use cblas_(d)scal when possible
	 */
	template<class Field>
	void
	fneg (const Field& F, const size_t n,
	      const typename Field::Element * Y, const size_t incY,
	      typename Field::Element * X, const size_t incX)
	{
		typename Field::Element * Xi = X ;
		const typename Field::Element * Yi = Y ;
		for (; Xi < X+n*incX; Xi+=incX,Yi+=incY  )
			F.neg( *Xi, *Yi );
	}

	/** \brief fzero : \f$A \gets 0 \f$.
	 * @param F field
	 * @param n number of elements to zero
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 */
	template<class Field>
	void
	fzero (const Field& F, const size_t n,
	       typename Field::Element *X, const size_t incX)
	{
		if (incX == 1) { // contigous data
			// memset(X,(int)F.zero,n); // might be bogus ?
			for (size_t i = 0 ; i < n ; ++i)
				F.assign(*(X+i), F.zero);

		}
		else { // not contiguous (strided)
			for (size_t i = 0 ; i < n ; ++i)
				F.assign(*(X+i*incX), F.zero);
		}
	}

	/** \brief fcopy : \f$x \gets y \f$.
	 * X is preallocated
	 * @todo variant for triagular matrix
	 * @param F field
	 * @param N size of the vectors
	 * \param [out] X vector in \p F
	 * \param incX stride of \p X
	 * \param [in] Y vector in \p F
	 * \param incY stride of \p Y
	 */
	template<class Field>
	void
	fcopy (const Field& F, const size_t N,
	       const typename Field::Element * Y, const size_t incY ,
	       typename Field::Element * X, const size_t incX);


	/** fscalin
	 * \f$x \gets a \cdot x\f$.
	 * @param F field
	 * @param n size of the vectors
	 * @param alpha homotÃ©ti scalar
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 * @bug use cblas_(d)scal when possible
	 * @internal
	 * @todo check if comparison with +/-1,0 is necessary.
	 */
	template<class Field>
	void
	fscalin (const Field& F, const size_t n, const typename Field::Element alpha,
		 typename Field::Element * X, const size_t incX);


	/** fscal
	 * \f$y \gets a \cdot x\f$.
	 * @param F field
	 * @param n size of the vectors
	 * @param alpha homotÃ©ti scalar
	 * \param[in] X vector in \p F
	 * \param incX stride of \p X
	 * \param[out] Y vector in \p F
	 * \param incY stride of \p Y
	 * @bug use cblas_(d)scal when possible
	 * @internal
	 * @todo check if comparison with +/-1,0 is necessary.
	 */
	template<class Field>
	void
	fscal (const Field& F, const size_t n
	       , const typename Field::Element alpha
	       , const typename Field::Element * X, const size_t incX
	       , typename Field::Element * Y, const size_t incY);



	/** \brief faxpy : \f$y \gets \alpha \cdot x + y\f$.
	 * @param F field
	 * @param N size of the vectors
	 * @param alpha scalar
	 * \param[in] X vector in \p F
	 * \param incX stride of \p X
	 * \param[in,out] Y vector in \p F
	 * \param incY stride of \p Y
	 */
	template<class Field>
	void
	faxpy (const Field& F, const size_t N,
	       const typename Field::Element alpha,
	       const typename Field::Element * X, const size_t incX,
	       typename Field::Element * Y, const size_t incY );

	/** \brief faxpby : \f$y \gets \alpha \cdot x + \beta \cdot y\f$.
	 * @param F field
	 * @param N size of the vectors
	 * @param alpha scalar
	 * \param[in] X vector in \p F
	 * \param incX stride of \p X
	 * \param beta scalar
	 * \param[in,out] Y vector in \p F
	 * \param incY stride of \p Y
	 * \note this is a catlas function
	 */
	template<class Field>
	void
	faxpby (const Field& F, const size_t N,
		const typename Field::Element alpha,
		const typename Field::Element * X, const size_t incX,
		const typename Field::Element beta,
		typename Field::Element * Y, const size_t incY );


	/** \brief fdot: dot product \f$x^T  y\f$.
	 * @param F field
	 * @param N size of the vectors
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 * \param Y vector in \p F
	 * \param incY stride of \p Y
	 */
	template<class Field>
	typename Field::Element
	fdot (const Field& F, const size_t N,
	      const typename Field::Element * X, const size_t incX,
	      const typename Field::Element * Y, const size_t incY );

	/** \brief fswap: \f$ X \leftrightarrow Y\f$.
	 * @bug use cblas_dswap when double
	 * @param F field
	 * @param N size of the vectors
	 * \param X vector in \p F
	 * \param incX stride of \p X
	 * \param Y vector in \p F
	 * \param incY stride of \p Y
	 */
	template<class Field>
	void
	fswap (const Field& F, const size_t N, typename Field::Element * X, const size_t incX,
	       typename Field::Element * Y, const size_t incY )
	{

		typename Field::Element tmp;
		typename Field::Element * Xi = X;
		typename Field::Element * Yi=Y;
		for (; Xi < X+N*incX; Xi+=incX, Yi+=incY ){
			F.assign( tmp, *Xi );
			F.assign( *Xi, *Yi );
			F.assign( *Yi, tmp );
		}
	}

	template <class Field>
	void
	fadd (const Field& F,  const size_t N,
	      const typename Field::Element* A, const size_t inca,
	      const typename Field::Element* B, const size_t incb,
	      typename Field::Element* C, const size_t incc);

	template <class Field>
	void
	fsub (const Field& F,  const size_t N,
	      const typename Field::Element* A, const size_t inca,
	      const typename Field::Element* B, const size_t incb,
	      typename Field::Element* C, const size_t incc);

	template <class Field>
	void
	faddin (const Field& F,  const size_t N,
		const typename Field::Element* B, const size_t incb,
		typename Field::Element* C, const size_t incc);

	template <class Field>
	void
	fsubin (const Field& F,  const size_t N,
		typename Field::Element* C, const size_t incc);


	template <class Field>
	void
	fadd (const Field& F,  const size_t N,
	      const typename Field::Element* A, const size_t inca,
	      const typename Field::Element alpha,
	      const typename Field::Element* B, const size_t incb,
	      typename Field::Element* C, const size_t incc);

} // FFLAS


#endif // __FFLASFFPACK_fflas_fflas_level1_INL
