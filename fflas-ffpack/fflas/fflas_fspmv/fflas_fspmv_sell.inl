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

/** @file fflas/fflas_fspmv_sell.inl
 * NO DOC
*/

#ifndef __FFLASFFPACK_fflas_fflas_spmv_sell_INL
#define __FFLASFFPACK_fflas_fflas_spmv_sell_INL

namespace FFLAS {
 
	template<class Element, bool bSimd>
	struct SELL
	{
		Element * val = nullptr;
		uint32_t * col = nullptr;
		uint64_t * ptr = nullptr;
		uint32_t * chs = nullptr;
		uint32_t n = 0;
		uint32_t m = 0;
		uint32_t chunk = Simd<Element>::vect_size;
	};

	template<class Element, bool simd>
	struct SELL_sub : public SELL<Element, bSimd>
	{};

	template<class Element, bool bSimd>
	struct SELL_ZO : public SELL<Element, bSimd>
	{};

	namespace details{

		// y = Ax + by (generic)
		template<class Field>
		inline std::enable_if<!std::is_same<FieldTraits<Field>::value, FieldCategories::FloatingPointTag>::value, void>::type
		sp_fgemv(const Field & F, const uint32_t m, const uint32_t n, const uint32_t chunk, const uint64_t * ptr, const uint32_t * chs, const uint32_t * col, const typename Field::Element * dat,
				 const typename Field::Element *x, const typename Field::Element b, typename Field::Element * y)
		{
			for(size_t i = 0 ; i < nbChunk ; ++i)
			{
				for(size_t j = 0 ; j < chs[i] ; ++j)
				{
					for(size_t k = 0 ; k < chunk ; ++k)
					{
						F.axpyin(y[i*Chunk+k], dat[ptr[i]+j*chunk+k], x[col[ptr[i]+j*chunk+k]]);
					}
				}
			}
		}

		// y = Ax + by (numeric)
		template<class Field, bool bSimd>
		inline std::enable_if<std::is_same<FieldTraits<Field>::value, FieldCategories::FloatingPointTag>::value, void>::type
		sp_fgemv(const Field & F, const size_t m, const size_t n, const size_t ld, const index_t * col, const typename Field::Element * dat,
				 const typename Field::Element *x, const typename Field::Element b, typename Field::Element * y)
		{
			if (!bsimd) {
				for(size_t i = 0 ; i < nbChunk ; ++i)
				{
					for(size_t j = 0 ; j < chs[i] ; ++j)
					{
						for(size_t k = 0 ; k < chunk ; ++k)
						{
							y[i*Chunk+k]+=dat[ptr[i]+j*chunk+k]*x[col[ptr[i]+j*chunk+k]];
						}
					}
				}
			}
			else {
				size_t i = 0 ;
				using simd = Simd<typename Field::Element>;
				using vect_t = typename simd::vect_t;
				vect_t X,Y,D ;
				for ( ; i < m ; i += simd::vect_size) {
					for (index_t j = 0 ; j < ld ; ++j) {
						D = simd::load(dat+i*simd::vect_size*ld+j*simd::vect_size);
						X = simd::gather(x,col+i*simd::vect_size*ld+j*simd::vect_size);
						Y = simd::madd(Y,D,X);
					}
					simd::store(y+i,Y);
				}
				size_t deplacement = m -i*simd::vect_size ;
				for (index_t j = 0 ; j < ld ; ++j) {
					for (size_t ii = 0 ; ii < deplacement ; ++ii) {
						y[i*simd::vect_size+ii] += dat[i*simd::vect_size*ld+j*deplacement + ii]*x[col[i*simd::vect_size*ld+j*deplacement + ii]];
					}
				}
			}
		}

		template<class Field, bool bSimd>
		inline std::enable_if<!std::is_same<FieldTraits<Field>::value, FieldCategories::FloatingPointTag>::value, void>::type
		sp_fgemv_zo(const Field & F, const size_t m, const size_t n, const size_t ld, const index_t * col,
				 	const typename Field::Element *x, const typename Field::Element b, typename Field::Element * y)
		{
			if(!bSimd)
			{
				for (size_t i = 0 ; i < m ; ++i) {
					for (index_t j = 0 ; j < ld ; ++j) {
						F.addin(y[i],x[col[i*ld+j]]);
					}
				}
			}
			else
			{
				size_t i = 0 ;
				using simd = Simd<typename Field::Element>;
				using vect_t = typename simd::vect_t;
				for ( ; i < m ; i += simd::vect_size) {
					for (index_t j = 0 ; j < ld ; ++j) {
						for(size_t k = 0 ; k < simd::vect_size ; ++k)
								F.addin(y[i*simd::vect_size+k], x[col[i*simd::vect_size*ld+j*simd::vect_size+k]]);
					}
				}
				size_t deplacement = m -i*simd::vect_size ;
				for (index_t j = 0 ; j < ld ; ++j) {
					for (size_t ii = 0 ; ii < deplacement ; ++ii) {
						F.addin(y[i*simd::vect_size+ii], x[col[i*simd::vect_size*ld+j*deplacement + ii]]);
					}
				}
			}
		}

		// y = Ax + by (numeric)
		template<class Field, bool bSimd>
		inline std::enable_if<std::is_same<FieldTraits<Field>::value, FieldCategories::FloatingPointTag>::value, void>::type
		sp_fgemv_zo(const Field & F, const size_t m, const size_t n, const size_t ld, const index_t * col,
				 const typename Field::Element *x, const typename Field::Element b, typename Field::Element * y)
		{
			if (!bsimd) {
				for ( size_t i = 0 ;  i < m   ; ++i ) {
					for (index_t j = 0 ; j < ld ; ++j) {
						y[i] += x[col[i*ld+j]];
					}
				}
			}
			else {
				size_t i = 0 ;
				using simd = Simd<double>;
				using vect_t = typename simd::vect_t;
				vect_t X,Y,D ;
				for ( ; i < m ; i += simd::vect_size) {
					for (index_t j = 0 ; j < ld ; ++j) {
						X = simd::gather(x,col+i*simd::vect_size*ld+j*simd::vect_size);
						Y = simd::add(Y,X);
					}
					simd::store(y+i,Y);
				}
				size_t deplacement = m -i*simd::vect_size ;
				for (index_t j = 0 ; j < ld ; ++j) {
					for (size_t ii = 0 ; ii < deplacement ; ++ii) {
						y[i*simd::vect_size+ii] += x[col[i*simd::vect_size*ld+j*deplacement + ii]];
					}
				}
			}
		}

	}/* details */

} /* FFLAS */

#endif // __FFLASFFPACK_fflas_fflas_spmv_sell_INL

