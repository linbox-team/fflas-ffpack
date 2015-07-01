/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   BB <bbboyer@ncsu.edu>
 *		Bastien Vialla <bastien.vialla@lirmm.fr>
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

/** @file fflas/fflas_fspmv.inl
*/

#ifndef __FFLASFFPACK_fflas_fflas_fspmv_INL
#define __FFLASFFPACK_fflas_fflas_fspmv_INL

#include "fflas-ffpack/fflas/fflas_bounds.inl"
#include "fflas-ffpack/fflas/fflas_helpers.inl"

namespace FFLAS { /*  DNS */

	template<class Field>
	struct VECT {
		size_t m;
		size_t inc;
		typename Field::Element_ptr dat;

		inline typename Field::Element_ptr data() {return dat;}
	};

	template<class Field>
	struct DNS {
		size_t n;
		size_t ld;
		typename Field::Element_ptr dat;

		inline typename Field::Element_ptr data() {return dat;}
	};

} // FFLAS

#if defined(__FFLASFFPACK_HAVE_CXX11) && defined(__FFLASFFPACK_USE_SIMD)

namespace FFLAS{ /* ELL */
	
	template<class Field>
	struct ELL {
		size_t m = 0;
		size_t n = 0;
		size_t ld = 0;
		index_t  * col = nullptr;
		typename Field::Element_ptr dat = nullptr;
	};

	template<class Field>
	struct ELL_sub : public ELL<Field> {
	};

	template<class Field>
	struct ELL_ZO : public ELL<Field> {
		typename Field::Element cst = 1;
	};
}

// #include "fflas-ffpack/fflas/fflas_fspmv/ell.inl"

#ifdef __FFLASFFPACK_USE_SIMD

namespace FFLAS{

	template<class Field>
	struct ELL_simd
	{
		size_t m = 0;
		size_t n = 0;
		size_t ld = 0;
		size_t chunk = 0;
		index_t  * col = nullptr;
		typename Field::Element_ptr dat = nullptr;
	};

	template<class Field>
	struct ELL_simd_sub : public ELL_simd<Field> {
	};

	template<class Field>
	struct ELL_simd_ZO : public ELL<Field> {
		typename Field::Element cst = 1;
	};
}

#include "fflas-ffpack/fflas/fflas_fspmv/ell_simd.inl"

namespace FFLAS{ /* SELL */

	template<class Field>
	struct SELL
	{
		size_t m = 0;
		size_t n = 0;
		size_t chunk = 0;
		size_t nChunks = 0;
		uint32_t * col = nullptr;
		uint64_t * perm = nullptr;
		uint64_t * ptr = nullptr;
		uint32_t * chs = nullptr;
		typename Field::Element_ptr dat;
	};

	template<class Field>
	struct SELL_sub : public SELL<Field>
	{};

	template<class Field>
	struct SELL_ZO : public SELL<Field>
	{
		typename Field::Element cst = 1;
	};

}

#endif // __FFLASFFPACK_USE_SIMD

namespace FFLAS{ /* CSR */
	template<class Field>
	struct CSR {
		index_t m = 0;
		index_t n = 0;
		index_t maxrow = 0;
		index_t  * st = nullptr;
		index_t  * col = nullptr;
		typename Field::Element_ptr dat ;
	};

	template<class Field>
	struct CSR_sub : public CSR<Field> {
	};

	template<class Field>
	struct CSR_ZO : public CSR<Field> {
		typename Field::Element cst = 1;
	};
}

#include "fflas-ffpack/fflas/fflas_fspmv/csr.inl"

namespace FFLAS{ /* COO */

	template<class Field>
	struct COO {
		index_t m  = 0;
		index_t n  = 0;
		uint64_t z = 0;
		index_t maxrow = 0;
		index_t  * row  = nullptr;
		index_t  * col = nullptr;
		typename Field::Element_ptr dat;
	};

	template<class Field>
	struct COO_sub : public COO<Field> {
	};

	template<class Field>
	struct COO_ZO : public COO<Field >{
		typename Field::Element cst = 1;
	};
}

#include "fflas-ffpack/fflas/fflas_fspmv/coo.inl"



namespace FFLAS { /* HYB */
#if 0
	template<class Element>
	struct SPADD {
		size_t ncsr;
		CSR<Element> * csr;
		size_t ncoo;
		COO<Element> * coo;
		size_t ndns;
		DNS<Element> * dns;
		size_t nell;
		ELL<Element> * ell;
		size_t nellr ;
		ELLR<Element> * ellr ;
		size_t ndia ;
		DIA<Element> * dia;

		SPADD() :
			ncsr(0)  ,csr(NULL)
			,ncoo(0) ,coo(NULL)
			,ndns(0) ,dns(NULL)
			,ndia(0) ,dia(NULL)
			,nell(0) ,ell(NULL)
			,nellr(0),ellr(NULL)
		{}
	};
#endif

} // FFLAS


namespace FFLAS{
	namespace details {

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y, FieldCategories::ModularTag)
		{
			if(b != 1)
			{
				if(b == 0)
				{
					for(size_t i = 0 ; i < m ; ++i)
						y[i] = 0;
				}
				else if(b == -1)
				{
					for(size_t i = 0 ; i < m ; ++i)
						y[i] *= -1;
				}
				else
				{
					fscalin(F, m, b, y, 1);
				}
			}
		}

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const size_t n,
						  const typename Field::Element b, typename Field::Element_ptr y,
						  const int lda, FieldCategories::ModularTag)
		{
			if(b != 1){
				if(b == 0){
					for(size_t i = 0 ; i < m ; ++i){
						for(size_t j = 0 ; j < n ; ++j){
							y[i*lda+j] = 0;
						}
					}
				}else if(b == -1){
					for(size_t i = 0 ; i < m ; ++i){
						for(size_t j = 0 ; j < n ; ++j){
							y[i*lda+j] *= -1;
						}
					}
				}else{
					fscalin(F, m, n, y, lda);
				}
			}
		}

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y, FieldCategories::UnparametricTag)
		{
			if(b != 1)
			{
				if(b == 0)
				{
					for(size_t i = 0 ; i < m ; ++i)
						y[i] = 0;
				}
				else if(b == -1)
				{
					for(size_t i = 0 ; i < m ; ++i)
						y[i] *= -1;
				}
				else
				{
					fscalin(F, m, b, y, 1);
				}
			}
		}

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const size_t n,
						  const typename Field::Element b, typename Field::Element_ptr y,
						  const int lda, FieldCategories::UnparametricTag)
		{
			if(b != 1){
				if(b == 0){
					for(size_t i = 0 ; i < m ; ++i){
						for(size_t j = 0 ; j < n ; ++j){
							y[i*lda+j] = 0;
						}
					}
				}else if(b == -1){
					for(size_t i = 0 ; i < m ; ++i){
						for(size_t j = 0 ; j < n ; ++j){
							y[i*lda+j] *= -1;
						}
					}
				}else{
					fscalin(F, m, n, y, lda);
				}
			}
		}

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const typename Field::Element b, typename Field::Element_ptr y, FieldCategories::GenericTag)
		{
			if(!F.isOne(b))
			{
				if(F.isZero(b))
				{
					for(size_t i = 0 ; i < m ; ++i)
						F.assign(y[i], F.zero);
				}
				else if(F.isMOne(b))
				{
					for(size_t i = 0 ; i < m ; ++i)
						F.negin(y[i]);
				}
				else
				{
					fscalin(F, m, b, y, 1);
				}
			}
		}

		template<class Field>
		inline void init_y(const Field & F, const size_t m, const size_t n,
						  const typename Field::Element b, typename Field::Element_ptr y,
						  const int lda, FieldCategories::GenericTag)
		{
			if(b != 1){
				if(b == 0){
					for(size_t i = 0 ; i < m ; ++i){
						for(size_t j = 0 ; j < n ; ++j){
							F.assign(y[i*lda+j], F.zero);
						}
					}
				}else if(b == -1){
					for(size_t i = 0 ; i < m ; ++i){
						for(size_t j = 0 ; j < n ; ++j){
							F.negin(y[i*lda+j]);
						}
					}
				}else{
					fscalin(F, m, n, y, lda);
				}
			}
		}

	}/* details */
}/* FFLAS */
#else // __FFLASFFPACK_HAVE_CXX11
namespace FFLAS{

	template<class Field>
	struct COO {
		index_t m;
		index_t n;
		uint64_t z;
		index_t maxrow;
		index_t  * row  ;
		index_t  * col ;
		typename Field::Element_ptr dat;
	};

	template<class Field>
	void sp_delete(COO<Field> & M){
		fflas_delete(M.row);
		fflas_delete(M.col);
		fflas_delete(M.dat);
	}

	template<class Field>
	void sp_fgemv(const Field& F, const COO<Field> & A, const VECT<Field> & x, const typename Field::Element & b, VECT<Field> & y){
		if(!F.isOne(b)){
			if(F.isMOne(b))
			{
				fnegin(F, A.m, y.dat, 1);
			}else if(F.isZero(b)){
				for(size_t i = 0 ; i < A.m ; ++i)
					F.assign(y.dat[i], F.zero);
			}
			else{
				fscalin(F, A.m, b,  y.dat, 1);
			}
		}
		for (index_t j = 0 ; j < A.z ; ++j)
			F.axpyin(y.dat[A.row[j]],A.dat[j],x.dat[A.col[j]]);
	}
}
#endif // __FFLASFFPACK_HAVE_CXX11

#endif // __FFLASFFPACK_fflas_fflas_fspmv_INL
