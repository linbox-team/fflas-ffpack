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

/** @file fflas/fflas_sparse_fgemv.h
*/

#ifndef __FFLASFFPACK_fflas_fflas_sparse_fgemv_H
#define __FFLASFFPACK_fflas_fflas_sparse_fgemv_H


#define index_t size_t

#ifdef __FFLASFFPACK_HAVE_MKL
// #include <mkl.h>
// #include <mkl_spblas.h>
#undef index_t
#define index_t MKL_INT
#endif

namespace FFLAS { /*  DNS */

	template<class Element>
	struct VECT ;

	template<class Element>
	struct DNS ;

}

namespace FFLAS { /*  CSR */

	template<class Element>
	struct CSR ;

	template<class Element>
	struct CSR_sub ;

	template<class Element>
	struct CSR_ZO ;


	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const CSR<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     );


	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      // const FFLAS_TRANSPOSE tA,
		      const CSR_sub<typename Field::Element> & A,
		      const VECT<typename Field::Element> & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element> & y
		     );

	template<class Field>
	void sp_fgemv(
		      const Field & F,
		      // const FFLAS_TRANSPOSE tA,
		      const CSR_ZO<typename Field::Element> & A,
		      const VECT<typename Field::Element > & x,
		      const typename Field::Element & b,
		      VECT<typename Field::Element > & y
		     );


}

namespace FFLAS { /*  ELL */

	template<class Element>
	struct ELL ;

	template<class Element>
	struct ELLR ;
} // FFLAS

namespace FFLAS { /*  COO */

} // FFLAS

namespace FFLAS { /*  BCSR */

} // FFLAS

namespace FFLAS { /*  DIA */

} // FFLAS

namespace FFLAS { /*  SKY */

} // FFLAS

namespace FFLAS { /*  JAG */

} // FFLAS


#include "fflas-ffpack/fflas/fflas_sparse_fgemv.inl"

#endif // __FFLASFFPACK_fflas_fflas_sparse_fgemv_H
