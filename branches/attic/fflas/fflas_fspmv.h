/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
 *              Bastien Vialla <bastien.vialla@lirmm.fr>
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

/** @file fflas/fflas_fspmv.h
*/

#ifndef __FFLASFFPACK_fflas_fflas_fspmv_H
#define __FFLASFFPACK_fflas_fflas_fspmv_H


#ifndef index_t
#define index_t size_t
#endif

#ifndef DIVIDE_INTO
#define DIVIDE_INTO(x,y) (((x) + (y) - 1)/(y))
#endif

#include "fflas-ffpack/config.h"
#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/field/field-traits.h"
#include <type_traits>
 // #include "fflas-ffpack/fflas/fflas_helpers.inl"

#ifdef __FFLASFFPACK_HAVE_MKL
#ifndef _MKL_H_ // temporary
#error "MKL (mkl.h) not present, while you have MKL enabled"
#endif
#undef index_t
#define index_t MKL_INT
#endif

namespace FFLAS { /*  DNS */

	template<class Field>
	struct VECT ;

	template<class Fiedl>
	struct DNS ;

}

namespace FFLAS { /*  COO */

	template<class Field>
	struct COO ;

	template<class Field>
	struct COO_sub ;

	template<class Field>
	struct COO_ZO ;


	template<class Field>
	void fspmv(
		      const Field& F,
		      const COO<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmv(
		      const Field& F,
		      const COO_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmv(
		      const Field & F,
		      const COO_ZO<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const COO<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const COO_sub<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const COO_ZO<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

} // FFLAS

namespace FFLAS { /*  CSR */

	template<class Field>
	struct CSR ;

	template<class Field>
	struct CSR_sub ;

	template<class Field>
	struct CSR_ZO ;


	template<class Field>
	void fspmv(
		      const Field& F,
		      const CSR<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );


	// y = A.x + b y
	template<class Field>
	void fspmv(
		      const Field& F,
		      const CSR_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmv(
		      const Field & F,
		      const CSR_ZO<Field> & A,
		      const VECT<Field > & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const CSR<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const CSR_sub<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const CSR_ZO<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );


}

namespace FFLAS { /*  ELL */

	 template<class Field>
	 struct ELL;

	template<class Field>
	struct ELL_sub;

	template<class Field>
	struct ELL_ZO;

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	// y = A.x + b y
	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELL_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmv(
		      const Field & F,
		      const ELL_ZO<Field> & A,
		      const VECT<Field > & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const ELL<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const ELL_sub<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

	template<class Field>
	void fspmm(
		      const Field& F,
		      const ELL_ZO<Field> & A,
		      const int blockSize,
		      const typename Field::Element_ptr & x,
		      const int ldx,
		      const typename Field::Element & b,
		      typename Field::Element_ptr & y,
		      const int ldy
		     );

} // FFLAS

namespace FFLAS{ /* ELLR */

	template<class Field>
	struct ELLR;

	template<class Field>
	struct ELLR_sub;

	template<class Field>
	struct ELLR_ZO;

	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELLR<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	// y = A.x + b y
	template<class Field>
	void fspmv(
		      const Field& F,
		      const ELLR_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmv(
		      const Field & F,
		      const ELLR_ZO<Field> & A,
		      const VECT<Field > & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

} // FFLAS

namespace FFLAS { /* SELL */

	template<class Field>
	struct SELL;

	template<class Field>
	struct SELL_sub;

	template<class Field>
	struct SELL_ZO;

	template<class Field>
	void fspmv(
		      const Field& F,
		      const SELL<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	// y = A.x + b y
	template<class Field>
	void fspmv(
		      const Field& F,
		      const SELL_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void fspmv(
		      const Field & F,
		      const SELL_ZO<Field> & A,
		      const VECT<Field > & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );
}

namespace FFLAS{
	namespace details{
		template<class Field>
		inline void init_y(const Field & F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y, FieldCategories::ModularTag);

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y, FieldCategories::UnparametricTag);

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y, FieldCategories::GenericTag);		

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const size_t n,
						  const typename Field::Element b, typename Field::Element_ptr y,
						  const int lda, FieldCategories::UnparametricTag);

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const size_t n,
						  const typename Field::Element b, typename Field::Element_ptr y,
						  const int lda, FieldCategories::GenericTag);

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const size_t n,
						  const typename Field::Element b, typename Field::Element_ptr y,
						  const int lda, FieldCategories::ModularTag);
	}	
}


#include "fflas-ffpack/fflas/fflas_fspmv.inl"
// #include "fflas-ffpack/fflas/fflas_fspmv/ellr.inl"
// #include "fflas-ffpack/fflas/fflas_fspmv/sell.inl"

#endif // __FFLASFFPACK_fflas_fflas_fspmv_H
