/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla <bastien.vialla@lirmm.fr>
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

#ifndef __FFLASFFPACK_fflas_sparse_ELL_simd_pspmv_INL
#define __FFLASFFPACK_fflas_sparse_ELL_simd_pspmv_INL

#ifdef __FFLASFFPACK_USE_TBB
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#endif

namespace FFLAS{
	namespace sparse_details_impl{
	template<class Field>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::GenericTag){
#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nbChunks, 2),
				[&F, &A, &x, &y](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						index_t j = 0;
						for(; i < A.ld ; ++j){
							for(index_t k = 0 ; k < A.chunk ; ++k){
								F.axpyin(y[i*A.chunk+k], A.dat[i*A.ld*A.chunk+j*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
							}
						}
					}
				}
		);
#else
#pragma omp parallel for
		for(index_t i = 0 ; i < A.nChunks ; ++i){
			index_t j = 0;
			for(; i < A.ld ; ++j){
				for(index_t k = 0 ; k < A.chunk ; ++k){
					F.axpyin(y[i*A.chunk+k], A.dat[i*A.ld*A.chunk+j*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
				}
			}
		}
#endif
	}

#ifdef __FFLASFFPACK_USE_SIMD
	template<class Field>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::UnparametricTag){
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;
		
#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nChunks, 2),
				[&F, &A, &x, &y](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						index_t j = 0;
						vect_t y1, y2, x1, x2, dat1, dat2, yy;
						y1 = simd::zero();
						y2 = simd::zero();
						for(; j < ROUND_DOWN(A.ld, 2) ; j+=2){
							dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
							dat2 = simd::load(A.dat+i*A.ld*A.chunk+(j+1)*A.chunk);
							x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
							x2 = simd::gather(x, A.col+i*A.ld*A.chunk+(j+1)*A.chunk);
							y1 = simd::fmadd(y1, dat1, x1);
							y2 = simd::fmadd(y2, dat2, x2);
						}
						for(; j < A.ld ; ++j){
							dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
							x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
							y1 = simd::fmadd(y1, dat1, x1);
					}
					yy = simd::load(y+i*A.chunk);
					simd::store(y+i*A.chunk, simd::add(yy, simd::add(y1, y2)));
				}
			}
		);
#else
#pragma omp parallel for
		for(index_t i = 0 ; i < A.nChunks ; ++i){
			index_t j = 0;
			vect_t y1, y2, x1, x2, dat1, dat2, yy;
			y1 = simd::zero();
			y2 = simd::zero();
			for(; j < ROUND_DOWN(A.ld, 2) ; j+=2){
				dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
				dat2 = simd::load(A.dat+i*A.ld*A.chunk+(j+1)*A.chunk);
				x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
				x2 = simd::gather(x, A.col+i*A.ld*A.chunk+(j+1)*A.chunk);
				y1 = simd::fmadd(y1, dat1, x1);
				y2 = simd::fmadd(y2, dat2, x2);
			}
			for(; j < A.ld ; ++j){
				dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
				x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
				y1 = simd::fmadd(y1, dat1, x1);
			}
			yy = simd::load(y+i*A.chunk);
			simd::store(y+i*A.chunk, simd::add(yy, simd::add(y1, y2)));
		}
#endif
	}
#else
	template<class Field>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::UnparametricTag){
#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nChunks, 2),
				[&F, &A, &x, &y](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						for(index_t j = 0 ; j < A.ld ; ++j){
							int k = 0;
							for( ; k < ROUND_DOWN(A.chunk, 4) ; k+=4){
								y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k]];
								y[i*A.chunk+k+1] += A.dat[i*A.ld*A.chunk+j*A.chunk+k+1]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k+1]];
								y[i*A.chunk+k+2] += A.dat[i*A.ld*A.chunk+j*A.chunk+k+2]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k+2]];
								y[i*A.chunk+k+3] += A.dat[i*A.ld*A.chunk+j*A.chunk+k+3]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k+3]];
							}
							for(; k < A.chunk ; ++k)
								y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k]];
						}
					}
				}
		);
#else
#pragma omp parallel for 
		for(index_t i = 0 ; i < A.nChunks ; ++i){
			for(index_t j = 0 ; j < A.ld ; ++j){
				int k = 0;
				for( ; k < ROUND_DOWN(A.chunk, 4) ; k+=4){
					y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k]];
					y[i*A.chunk+k+1] += A.dat[i*A.ld*A.chunk+j*A.chunk+k+1]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k+1]];
					y[i*A.chunk+k+2] += A.dat[i*A.ld*A.chunk+j*A.chunk+k+2]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k+2]];
					y[i*A.chunk+k+3] += A.dat[i*A.ld*A.chunk+j*A.chunk+k+3]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k+3]];
				}
				for(; k < A.chunk ; ++k)
					y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+J*A.chunk+k]];
			}
		}
#endif // TBB
	}
#endif // SIMD

#ifdef __FFLASFFPACK_USE_SIMD
	template<class Field>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, const uint64_t kmax){
		index_t block = (A.ld)/kmax ; // use DIVIDE_INTO from pfspmvgpu
		index_t chunk = A.chunk;
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;

		vect_t X, Y, D, C, Q, TMP, NEGP, INVP, MIN, MAX, P;
		double p = (typename Field::Element)F.characteristic();

		P = simd::set1(p);
		NEGP = simd::set1(-p);
		INVP = simd::set1(1/p);
		MIN = simd::set1(F.minElement());
		MAX = simd::set1(F.maxElement());

#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nChunks, 2),
				[&F, &A, &x, &y, P, NEGP, INVP, MAX, MIN](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						index_t j = 0 ;
						index_t j_loc = 0 ;
						Y = simd::load(y+i*chunk);
						for (size_t l = 0 ; l < block ; ++l) {
							j_loc += kmax ;
							for ( ; j < j_loc ; ++j) {
								D = simd::load(A.dat+i*A.chunk*A.ld+j*A.chunk);
								X = simd::gather(x, A.col+i*A.chunk*A.ld+j*A.chunk);
								Y = simd::fmadd(Y,D,X);
							}
							simd::mod(Y,P, INVP, NEGP, MIN, MAX, Q, TMP);
						}
						for ( ; j < A.ld ; ++j) {
							D = simd::load(A.dat+i*A.chunk*A.ld+j*A.chunk);
							X = simd::gather(x,A.col+i*A.chunk*A.ld+j*A.chunk);
							Y = simd::fmadd(Y,D,X);
						}
						simd::mod(Y,P, INVP, NEGP, MIN, MAX, Q, TMP);
						simd::store(y+i*A.chunk,Y);
					}
				}
		);
#else

#pragma omp parallel for
		for ( size_t i = 0; i < A.nChunks ; ++i ) {
			index_t j = 0 ;
			index_t j_loc = 0 ;
			Y = simd::load(y+i*chunk);
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;

				for ( ; j < j_loc ; ++j) {
					D = simd::load(A.dat+i*A.chunk*A.ld+j*A.chunk);
					X = simd::gather(x, A.col+i*A.chunk*A.ld+j*A.chunk);
					Y = simd::fmadd(Y,D,X);
				}
				simd::mod(Y,P, INVP, NEGP, MIN, MAX, Q, TMP);
			}
			for ( ; j < A.ld ; ++j) {
				D = simd::load(A.dat+i*A.chunk*A.ld+j*A.chunk);
				X = simd::gather(x,A.col+i*A.chunk*A.ld+j*A.chunk);
				Y = simd::fmadd(Y,D,X);
			}
			simd::mod(Y,P, INVP, NEGP, MIN, MAX, Q, TMP);
			simd::store(y+i*A.chunk,Y);
		}
#endif // TBB
	}

#else
	template<class Field>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, const uint64_t kmax){
		index_t block = (A.ld)/kmax ; // use DIVIDE_INTO from pfspmvgpu
		index_t chunk = A.chunk;
		size_t end = (A.m%chunk == 0)? A.m : A.m+A.m%chunk;
#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nChunks, 2),
				[&F, &A, &x, &y, P, NEGP, INVP, MAX, MIN](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						index_t j = 0 ;
						index_t j_loc = 0 ;
						for (size_t l = 0 ; l < block ; ++l) {
							j_loc += kmax ;
							for ( ; j < j_loc ; ++j) {
								for(int k = 0 ; k < A.chunk ; ++k){
									y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+j*A.chunk+k]];
								}
							}
							for(int k = 0 ; k < A.chunk ; ++k)
								F.reduce(y[i*A.chunk+k], y[i*A.chunk+k]);
						}
						for ( ; j < A.ld ; ++j) {
							for(int k = 0 ; k < A.chunk ; ++k){
								y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+j*A.chunk+k]];
							}
						}
						for(int k = 0 ; k < A.chunk ; ++k)
							F.reduce(y[i*A.chunk+k], y[i*A.chunk+k]);
					}
				}
		);
#else
#pragma omp parallel for
		for ( size_t i = 0; i < A.nChunks ; ++i ) {
			index_t j = 0 ;
			index_t j_loc = 0 ;
			for (size_t l = 0 ; l < block ; ++l) {
				j_loc += kmax ;

				for ( ; j < j_loc ; ++j) {
					for(int k = 0 ; k < A.chunk ; ++k){
						y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+j*A.chunk+k]];
					}
				}
				for(int k = 0 ; k < A.chunk ; ++k)
					F.reduce(y[i*A.chunk+k], y[i*A.chunk+k]);
			}
			for ( ; j < A.ld ; ++j) {
				for(int k = 0 ; k < A.chunk ; ++k){
					y[i*A.chunk+k] += A.dat[i*A.ld*A.chunk+j*A.chunk+k]*x[A.col[i*A.ld*A.chunk+j*A.chunk+k]];
				}
			}
			for(int k = 0 ; k < A.chunk ; ++k)
			F.reduce(y[i*A.chunk+k], y[i*A.chunk+k]);
		}
#endif // TBB
	}
#endif // SIMD

	template<class Field, class Func>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd_ZO> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, Func && func, FieldCategories::GenericTag){
#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nChunks, 2),
				[&F, &A, &x, &y, P, NEGP, INVP, MAX, MIN](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						index_t j = 0;
						for(; i < A.ld ; ++j){
							index_t k = 0;
							for( ; k < ROUND_DOWN(A.chunk, 4) ; k+=4){
								func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
								func(y[i*A.chunk+k+1], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+1]]);
								func(y[i*A.chunk+k+2], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+2]]);
								func(y[i*A.chunk+k+3], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+3]]);
							}
							for(; k < A.chunk ; ++k)
								func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
						}
					}
				}
		);
#else
#pragma omp parallel for
		for(index_t i = 0 ; i < A.nChunks ; ++i){
			index_t j = 0;
			for(; i < A.ld ; ++j){
				index_t k = 0;
				for( ; k < ROUND_DOWN(A.chunk, 4) ; k+=4){
					func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
					func(y[i*A.chunk+k+1], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+1]]);
					func(y[i*A.chunk+k+2], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+2]]);
					func(y[i*A.chunk+k+3], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+3]]);
				}
				for(; k < A.chunk ; ++k)
					func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
			}
		}
#endif
	}

#ifdef __FFLASFFPACK_USE_SIMD
	template<class Field, class Func>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd_ZO> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, Func && func, FieldCategories::UnparametricTag){
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;
#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nChunks, 2),
				[&F, &A, &x, &y, P, NEGP, INVP, MAX, MIN](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						index_t j = 0;
						vect_t y1, y2, x1, x2, dat1, dat2, yy;
						y1 = simd::zero();
						y2 = simd::zero();
						for(; j < ROUND_DOWN(A.ld, 2) ; j+=2){
							dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
							dat2 = simd::load(A.dat+i*A.ld*A.chunk+(j+1)*A.chunk);
							x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
							x2 = simd::gather(x, A.col+i*A.ld*A.chunk+(j+1)*A.chunk);
							y1 = func(y1, x1);
							y1 = func(y2, x2);
						}
						for(; j < A.ld ; ++j){
							dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
							x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
							y1 = func(y1, dat1, x1);
						}
						yy = simd::load(y+i*A.chunk);
						simd::store(y+i*A.chunk, func(yy, func(y1, y2)));
					}
				}
		);
#else
#pragma omp parallel for
		for(index_t i = 0 ; i < A.nChunks ; ++i){
			index_t j = 0;
			vect_t y1, y2, x1, x2, dat1, dat2, yy;
			y1 = simd::zero();
			y2 = simd::zero();
			for(; j < ROUND_DOWN(A.ld, 2) ; j+=2){
				dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
				dat2 = simd::load(A.dat+i*A.ld*A.chunk+(j+1)*A.chunk);
				x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
				x2 = simd::gather(x, A.col+i*A.ld*A.chunk+(j+1)*A.chunk);
				y1 = func(y1, x1);
				y1 = func(y2, x2);
			}
			for(; j < A.ld ; ++j){
				dat1 = simd::load(A.dat+i*A.ld*A.chunk+j*A.chunk);
				x1 = simd::gather(x, A.col+i*A.ld*A.chunk+j*A.chunk);
				y1 = func(y1, dat1, x1);
			}
			yy = simd::load(y+i*A.chunk);
			simd::store(y+i*A.chunk, func(yy, func(y1, y2)));
		}
#endif
	}

#else // SIMD

	template<class Field, class Func>
	inline void pfspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_simd_ZO> & A, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, Func && func, FieldCategories::UnparametricTag){
#ifdef __FFLASFFPACK_USE_TBB
		tbb::parallel_for(
			tbb::blocked_range<index_t>(0, A.nChunks, 2),
				[&F, &A, &x, &y, P, NEGP, INVP, MAX, MIN](const tbb::blocked_range<index_t> & r){
					for(index_t i = r.begin(), end = r.end() ; i < end ; ++i){
						index_t j = 0;
						for(; i < A.ld ; ++j){
							index_t k = 0;
							for( ; k < ROUND_DOWN(A.chunk, 4) ; k+=4){
								func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
								func(y[i*A.chunk+k+1], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+1]]);
								func(y[i*A.chunk+k+2], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+2]]);
								func(y[i*A.chunk+k+3], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+3]]);
							}
							for(; k < A.chunk ; ++k)
								func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
						}
					}
				}
		);
#else
#pragma omp parallel for
		for(index_t i = 0 ; i < A.nChunks ; ++i){
			index_t j = 0;
			for(; i < A.ld ; ++j){
				index_t k = 0;
				for( ; k < ROUND_DOWN(A.chunk, 4) ; k+=4){
					func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
					func(y[i*A.chunk+k+1], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+1]]);
					func(y[i*A.chunk+k+2], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+2]]);
					func(y[i*A.chunk+k+3], x[A.col[i*A.ld*A.chunk+j*A.chunk+k+3]]);
				}
				for(; k < A.chunk ; ++k)
					func(y[i*A.chunk+k], x[A.col[i*A.ld*A.chunk+j*A.chunk+k]]);
			}
		}
#endif // TBB
	}
#endif // SIMD
}// ELL_simd_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_ELL_simd_pspmv_INL