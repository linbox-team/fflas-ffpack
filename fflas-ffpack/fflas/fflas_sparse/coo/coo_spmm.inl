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

#ifndef __FFLASFFPACK_fflas_sparse_coo_spmm_INL
#define __FFLASFFPACK_fflas_sparse_coo_spmm_INL

namespace FFLAS{
	namespace sparse_details_impl{
	
	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::GenericTag){
		for(index_t i = 0 ; i < A.nnz ; ++i){
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, 4) ; k+=4){
				F.axpyin(y[A.row[i]*blockSize+k], A.dat[i], x[A.col[i]*blockSize+k]);
				F.axpyin(y[A.row[i]*blockSize+k+1], A.dat[i], x[A.col[i]*blockSize+k+1]);
				F.axpyin(y[A.row[i]*blockSize+k+2], A.dat[i], x[A.col[i]*blockSize+k+2]);
				F.axpyin(y[A.row[i]*blockSize+k+3], A.dat[i], x[A.col[i]*blockSize+k+3]);
			}
			for(; k < blockSize ; ++k)
				F.axpyin(y[A.row[i]*blockSize+k], A.dat[i], x[A.col[i]*blockSize+k]);
		}
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, FieldCategories::GenericTag){
		for(index_t i = 0 ; i < A.nnz ; ++i){
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, 4) ; k+=4){
				F.axpyin(y[A.row[i]*ldy+k], A.dat[i], x[A.col[i]*ldx+k]);
				F.axpyin(y[A.row[i]*ldy+k+1], A.dat[i], x[A.col[i]*ldx+k+1]);
				F.axpyin(y[A.row[i]*ldy+k+2], A.dat[i], x[A.col[i]*ldx+k+2]);
				F.axpyin(y[A.row[i]*ldy+k+3], A.dat[i], x[A.col[i]*ldx+k+3]);
			}
			for(; k < blockSize ; ++k)
				F.axpyin(y[A.row[i]*ldy+k], A.dat[i], x[A.col[i]*ldx+k]);
		}
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FieldCategories::UnparametricTag){
		for(index_t i = 0 ; i < A.nnz ; ++i){
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, 4) ; k+=4){
				y[A.row[i]*blockSize+k] += A.dat[i]*x[A.col[i]*blockSize+k];
				y[A.row[i]*blockSize+k+1] += A.dat[i]*x[A.col[i]*blockSize+k+1];
				y[A.row[i]*blockSize+k+2] += A.dat[i]*x[A.col[i]*blockSize+k+2];
				y[A.row[i]*blockSize+k+3] += A.dat[i]*x[A.col[i]*blockSize+k+3];
			}
			for(; k < blockSize ; ++k)
				y[A.row[i]*blockSize+k] += A.dat[i]*x[A.col[i]*blockSize+k];
		}
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, FieldCategories::UnparametricTag){
		for(index_t i = 0 ; i < A.nnz ; ++i){
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, 4) ; k+=4){
				y[A.row[i]*ldy+k] += A.dat[i]*x[A.col[i]*ldx+k];
				y[A.row[i]*ldy+k+1] += A.dat[i]*x[A.col[i]*ldx+k+1];
				y[A.row[i]*ldy+k+2] += A.dat[i]*x[A.col[i]*ldx+k+2];
				y[A.row[i]*ldy+k+3] += A.dat[i]*x[A.col[i]*ldx+k+3];
			}
			for(; k < blockSize ; ++k)
				y[A.row[i]*ldy+k] += A.dat[i]*x[A.col[i]*ldx+k];
		}
	}

#ifdef __FFLASFFPACK_USE_SIMD
	
	template<class Field, class LFunc, class SFunc>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, LFunc && lfunc, SFunc && sfunc, FieldCategories::UnparametricTag){
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;
		for(index_t i = 0 ; i < A.nnz ; ++i){	
			vect_t vy, vx, vdat;
			vdat = simd::set1(A.dat[i]);
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, simd::vect_size) ; k+=simd::vect_size){
				vy = lfunc(y+A.row[i]*blockSize+k);
				vx = lfunc(x+A.col[i]*blockSize+k);
				sfunc(y+A.row[i]*blockSize+k, simd::fmadd(vy, vdat, vx));
			}
			for(; k < blockSize ; ++k)
				y[A.row[i]*blockSize+k] += A.dat[i]*x[A.col[i]*blockSize+k];
		}
	}

	template<class Field, class LFunc, class SFunc>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, LFunc && lfunc, SFunc && sfunc, FieldCategories::UnparametricTag){
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_t;
		for(index_t i = 0 ; i < A.nnz ; ++i){	
			vect_t vy, vx, vdat;
			vdat = simd::set1(A.dat[i]);
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, simd::vect_size) ; k+=simd::vect_size){
				vy = lfunc(y+A.row[i]*ldy+k);
				vx = lfunc(x+A.col[i]*ldx+k);
				sfunc(y+A.row[i]*ldy+k, simd::fmadd(vy, vdat, vx));
			}
			for(; k < blockSize ; ++k)
				y[A.row[i]*ldy+k] += A.dat[i]*x[A.col[i]*ldx+k];
		}
	}

#endif // SIMD

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, const int64_t kmax){
		// TODO
	}

	template<class Field>
	inline void fspmm(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, const int64_t kmax){
		// TODO
	}

	template<class Field, class Func>
	inline void	fspmm_zo(const Field & F, const Sparse<Field, SparseMatrix_t::COO_ZO> & A, int blockSize, typename Field::ConstElement_ptr x,
			typename Field::Element_ptr y, Func && func){
		for(index_t i = 0 ; i < A.nnz ; ++i){
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, 4) ; k+=4){
				func(y[A.row[i]*blockSize+k], x[A.col[i]*blockSize+k]);
				func(y[A.row[i]*blockSize+k+1],x[A.col[i]*blockSize+k+1]);
				func(y[A.row[i]*blockSize+k+2],x[A.col[i]*blockSize+k+2]);
				func(y[A.row[i]*blockSize+k+3],x[A.col[i]*blockSize+k+3]);
			}
			for(; k < blockSize ; ++k)
				func(y[A.row[i]*blockSize+k], x[A.col[i]*blockSize+k]);
		}
	}

	template<class Field, class Func>
	inline void	fspmm_zo(const Field & F, const Sparse<Field, SparseMatrix_t::COO_ZO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx, 
			typename Field::Element_ptr y, int ldy, Func && func){
		for(index_t i = 0 ; i < A.nnz ; ++i){
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, 4) ; k+=4){
				func(y[A.row[i]*ldy+k], x[A.col[i]*ldx+k]);
				func(y[A.row[i]*ldy+k+1],x[A.col[i]*ldx+k+1]);
				func(y[A.row[i]*ldy+k+2],x[A.col[i]*ldx+k+2]);
				func(y[A.row[i]*ldy+k+3],x[A.col[i]*ldx+k+3]);
			}
			for(; k < blockSize ; ++k)
				func(y[A.row[i]*ldy+k], x[A.col[i]*ldx+k]);
		}
	}

#ifdef __FFLASFFPACK_USE_SIMD

	template<class Field, class LFunc, class SFunc, class FuncVect, class FuncScal>
	inline void fspmm_zo(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x,
		      typename Field::Element_ptr y, FuncVect && vecfunc, FuncScal && scalfunc, LFunc && lfunc, SFunc && sfunc){
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_size;
		for(index_t i = 0 ; i < A.nnz ; ++i){	
			vect_t vy, vx;
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, simd::vect_size) ; k+=simd::vect_size){
				vy = lfunc(y+A.row[i]*blockSize+k);
				vx = lfunc(x+A.col[i]*blockSize+k);
				sfunc(y+A.row[i]*blockSize+k, vecfunc(vy,  vx));
			}
			for(; k < blockSize ; ++k)
				scalfunc(y[A.row[i]*blockSize+k], x[A.col[i]*blockSize+k]);
		}
	}

	template<class Field, class LFunc, class SFunc, class FuncVect, class FuncScal>
	inline void fspmm_zo(const Field & F, const Sparse<Field, SparseMatrix_t::COO> & A, int blockSize, typename Field::ConstElement_ptr x, int ldx,
		      typename Field::Element_ptr y, int ldy, FuncVect && vecfunc, FuncScal && scalfunc, LFunc && lfunc, SFunc && sfunc){
		using simd = Simd<typename Field::Element>;
		using vect_t = typename simd::vect_size;
		for(index_t i = 0 ; i < A.nnz ; ++i){	
			vect_t vy, vx;
			int k = 0;
			for( ; k < ROUND_DOWN(blockSize, simd::vect_size) ; k+=simd::vect_size){
				vy = lfunc(y+A.row[i]*ldy+k);
				vx = lfunc(x+A.col[i]*ldx+k);
				sfunc(y+A.row[i]*ldy+k, vecfunc(vy,  vx));
			}
			for(; k < blockSize ; ++k)
				scalfunc(y[A.row[i]*ldy+k], x[A.col[i]*ldx+k]);
		}
	}

#endif

}// coo_details

} // FFLAS

#endif //  __FFLASFFPACK_fflas_coo_spmm_INL