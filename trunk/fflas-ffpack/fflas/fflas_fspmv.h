/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
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
	void sp_fgemv(
		      const Field& F,
		      const COO<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );


	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      const COO_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      const COO_ZO<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
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
	void sp_fgemv(
		      const Field& F,
		      const CSR<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );


	// y = A.x + b y
	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      const CSR_sub<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );

	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      const CSR_ZO<Field> & A,
		      const VECT<Field > & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     );


}

namespace FFLAS { /*  CSC */

} // FFLAS

namespace FFLAS { /*  ELL */

	 template<class Field, bool simd>
	 struct ELL;

	 template<class Field, bool simd>
	 struct ELL_sub;

	 template<class Field, bool simd>
	 struct ELL_ZO;

} // FFLAS

namespace FFLAS{ /* ELLR */

	template<class Field>
	struct ELLR;

	template<class Field>
	struct ELLR_sub;

	template<class Field>
	struct ELLR_ZO;

} // FFLAS

namespace FFLAS { /* SELL */

	template<class Element>
	struct SELL;

	template<class Element>
	struct SELL_sub;

	template<class Element>
	struct SELL_ZO;
}

namespace FFLAS { /*  BCSR */

} // FFLAS

namespace FFLAS { /*  DIA */

} // FFLAS

namespace FFLAS { /*  SKY */

} // FFLAS

namespace FFLAS { /*  JAG */

} // FFLAS


#include "fflas-ffpack/fflas/fflas_fspmv.inl"
#include "fflas-ffpack/fflas/fflas_fspmv/coo.inl"
#include "fflas-ffpack/fflas/fflas_fspmv/csr.inl"
#include "fflas-ffpack/fflas/fflas_fspmv/ell.inl"
#include "fflas-ffpack/fflas/fflas_fspmv/ellr.inl"
// #include "fflas-ffpack/fflas/fflas_fspmv/sell.inl"

#endif // __FFLASFFPACK_fflas_fflas_fspmv_H
