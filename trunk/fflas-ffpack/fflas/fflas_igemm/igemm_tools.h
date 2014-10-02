/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2013,2014  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 * the code is inspired and adapted from the Eigen library
 * modified by BB <brice.boyer@lip6.fr>
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

#ifndef __FFLASFFPACK_fflas_igemm_igemm_tools_H
#define __FFLASFFPACK_fflas_igemm_igemm_tools_H


/* ***** */
/* TOOLS */
/* ***** */

namespace FFLAS { namespace details { /*  tools */

	// duplicate each entry into vector register
	template<size_t N>
	inline void duplicate_vect (int64_t* XX, const int64_t* X, size_t n){}

	template<size_t k>
	void pack_lhs(int64_t* XX, const int64_t* X, size_t ldx, size_t rows, size_t cols);

	template<size_t k>
	void pack_rhs(int64_t* XX, const int64_t* X, size_t ldx, size_t rows, size_t cols);

	void gebp(size_t rows, size_t cols, size_t depth,int64_t* C, size_t ldc, const int64_t* blockA, size_t lda,
		  const int64_t* BlockB, size_t ldb, int64_t* BlockW);

	void BlockingFactor(size_t& m, size_t& n, size_t& k);


} // details
} // FFLAS

#include "igemm_tools.inl"

#endif // __FFLASFFPACK_fflas_igemm_igemm_tools_H

