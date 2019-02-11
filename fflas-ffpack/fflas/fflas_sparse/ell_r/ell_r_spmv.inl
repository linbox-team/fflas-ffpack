/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#ifndef __FFLASFFPACK_fflas_sparse_ELL_R_spmv_INL
#define __FFLASFFPACK_fflas_sparse_ELL_R_spmv_INL

namespace FFLAS{
	namespace ell_r_details{
		template<class Field>
		inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_R> & A, typename Field::ConstElement_ptr x,
				  typename Field::Element_ptr y, FieldCategories::GenericTag){
			index_t start = 0;
			for(index_t i = 0 ; i < A.mRow ; ++i, start+=A.ld){
				index_t j = 0;
				typename Field::Element y1, y2, y3, y4;
				F.assign(y1, F.zero);
				F.assign(y2, F.zero);
				F.assign(y3, F.zero);
				F.assign(y4, F.zero);
				for(; j < ROUND_DOWN(A.ld, 4) ; j+=4){
					F.axpyin(y1,A.dat[start+j],x[A.col[start+j]]);
					F.axpyin(y2,A.dat[start+j+1],x[A.col[start+j+1]]);
					F.axpyin(y3,A.dat[start+j+2],x[A.col[start+j+2]]);
					F.axpyin(y4,A.dat[start+j+3],x[A.col[start+j+3]]);
				}
				for(; j < A.ld ; ++j){
					F.axpyin(y1,A.dat[start+j],x[A.col[start+j]]);
				}
				F.addin(y[A.row[i]], y1);
				F.addin(y[A.row[i]], y2);
				F.addin(y[A.row[i]], y3);
				F.addin(y[A.row[i]], y4);
			}
		}

		template<class Field>
		inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_R> & A, typename Field::ConstElement_ptr x,
				  typename Field::Element_ptr y, FieldCategories::UnparametricTag){
			index_t start = 0;
			for(index_t i = 0 ; i < A.mRow ; ++i, start+=A.ld){
				index_t j = 0;
				typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
				for(; j < ROUND_DOWN(A.ld, 4) ; j+=4){
					y1 += A.dat[start+j] * x[A.col[start+j]];
					y2 += A.dat[start+j+1] * x[A.col[start+j+1]];
					y3 += A.dat[start+j+2] * x[A.col[start+j+2]];
					y4 += A.dat[start+j+3] * x[A.col[start+j+3]];
				}
				for(; j < A.ld ; ++j){
					y1 += A.dat[start+j] * x[A.col[start+j]];
				}
				y[A.row[i]] += y1+y2+y3+y4;
			}
		}

		template<class Field>
		inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_R> & A, typename Field::ConstElement_ptr x,
				  typename Field::Element_ptr y, const int64_t kmax){
			index_t start = 0;
			index_t block = (A.ld)/kmax ;
			for (index_t i = 0 ; i < A.mRow ; ++i, start+=A.ld) {
				index_t j_loc = 0, j = 0;
				for (index_t l = 0 ; l < (index_t) block ; ++l) {
					j_loc += kmax ;
					for ( ; j < j_loc ; ++j) {
						y[A.row[i]] += A.dat[start+j] * x[A.col[start+j]];
					}
					F.reduce(y[A.row[i]]);
				}
				for ( ; j < A.ld ; ++j) {
					y[A.row[i]] += A.dat[start+j] * x[A.col[start+j]];
				}
				F.reduce(y[A.row[i]];
					 }
					 }

					 template<class Field, class Func>
					 inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A, typename Field::ConstElement_ptr x,
							   typename Field::Element_ptr y, Func && func, FieldCategories::GenericTag){
					 index_t start = 0;
					 for(index_t i = 0 ; i < A.mRow ; ++i, start+=A.ld){
					 index_t j = 0;
					 typename Field::Element y1, y2, y3, y4;
					 F.assign(y1, F.zero);
					 F.assign(y2, F.zero);
					 F.assign(y3, F.zero);
					 F.assign(y4, F.zero);
					 for(; j < ROUND_DOWN(A.ld, 4) ; j+=4){
					 func(y1,x[A.col[start+j]]);
					 func(y2,x[A.col[start+j+1]]);
					 func(y3,x[A.col[start+j+2]]);
					 func(y4,x[A.col[start+j+3]]);
					 }
					 for(; j < A.ld ; ++j){
						 func(y1,x[A.col[start+j]]);
					 }
					 F.addin(y[A.row[i]], y1);
					 F.addin(y[A.row[i]], y2);
					 F.addin(y[A.row[i]], y3);
					 F.addin(y[A.row[i]], y4);
					 }
					 }

				template<class Field, class Func>
				inline void fspmv(const Field & F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A, typename Field::ConstElement_ptr x,
						  typename Field::Element_ptr y, Func && func, FieldCategories::UnparametricTag){
					index_t start = 0;
					for(index_t i = 0 ; i < A.mRow ; ++i, start+=A.ld){
						index_t j = 0;
						typename Field::Element y1 = 0, y2 = 0, y3 = 0, y4 = 0;
						for(; j < ROUND_DOWN(A.ld, 4) ; j+=4){
							func(y1,x[A.col[start+j]]);
							func(y2,x[A.col[start+j+1]]);
							func(y3,x[A.col[start+j+2]]);
							func(y4,x[A.col[start+j+3]]);
						}
						for(; j < A.ld ; ++j){
							func(y1,x[A.col[start+j]]);
						}
						y[A.row[i]] += y1+y2+y3+y4;
					}
				}
	}// ELL_R_details

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R> & A, typename Field::ConstElement_ptr x,
			  const typename Field::Element & beta, typename Field::Element_ptr y){
		sparse_details::init_y(F, A.m, beta, y, typename FieldTraits<Field>::category());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R> & A, typename Field::ConstElement_ptr x,
			  typename Field::Element_ptr y, FieldCategories::GenericTag){
		ell_r_details::fspmv(F, A, x, y, FieldCategories::GenericTag());
	}

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R> & A, typename Field::ConstElement_ptr x,
			  typename Field::Element_ptr y, FieldCategories::UnparametricTag){
		ell_r_details::fspmv(F, A, x, y, FieldCategories::UnparametricTag());
	}

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R> & A, typename Field::ConstElement_ptr x,
			  typename Field::Element_ptr y, FieldCategories::ModularTag){
		if(A.delayed){
			ell_r_details::fspmv(F, A, x, y, FieldCategories::UnparametricTag());
			freduce(F, A.m, y, 1);
		}else{
			ell_r_details::fspmv(F, A, x, y, A.kmax);
		}
	}

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A, typename Field::ConstElement_ptr x,
			  const typename Field::Element & beta, typename Field::Element_ptr y){
		sparse_details::init_y(F, A.m, beta, y, typename FieldTraits<Field>::category());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A, typename Field::ConstElement_ptr x,
			  typename Field::Element_ptr y, FieldCategories::GenericTag){
		using Element = typename Field::Element;
		if(A.cst == 1){
			ell_r_details::fspmv(F, A, x, y, [&F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::GenericTag());
		}else if(A.cst == -1){
			ell_r_details::fspmv(F, A, x, y, [&F](Element & a, const Element & b){F.subin(a, b);}, FieldCategories::GenericTag());
		}else{
			auto x1 = fflas_new(F, A.n, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x, 1, x1, 1);
			ell_r_details::fspmv(F, A, x, y, [&F](Element & a, const Element & b){F.addin(a, b);}, FieldCategories::GenericTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A, typename Field::ConstElement_ptr x,
			  typename Field::Element_ptr y, FieldCategories::UnparametricTag){
		using Element = typename Field::Element;
		if(A.cst == 1){
			ell_r_details::fspmv(F, A, x, y, [](Element & a, const Element & b){a += b;}, FieldCategories::UnparametricTag());
		}else if(A.cst == -1){
			ell_r_details::fspmv(F, A, x, y, [](Element & a, const Element & b){a -= b;}, FieldCategories::UnparametricTag());
		}else{
			auto x1 = fflas_new(F, A.n, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x, 1, x1, 1);
			ell_r_details::fspmv(F, A, x, y, [](Element & a, const Element & b){a += b;}, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(const Field& F, const Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A, typename Field::ConstElement_ptr x,
			  typename Field::Element_ptr y, FieldCategories::ModularTag){
		fspmv(F, A, x, y, FieldCategories::UnparametricTag());
		freduce(F, A.m, y, 1);
	}

	template<class Field>
	inline void sparse_delete(const Sparse<Field, SparseMatrix_t::ELL_R> & A){
		fflas_delete(A.dat);
		fflas_delete(A.col);
	}

	template<class Field>
	inline void sparse_delete(const Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A){
		fflas_delete(A.col);
	}

	template<class Field, class IndexT>
	inline void sparse_init(const Field & F, Sparse<Field, SparseMatrix_t::ELL_R> & A,
				const IndexT * row, const IndexT * col, typename Field::ConstElement_ptr dat,
				uint64_t rowdim, uint64_t coldim, uint64_t nnz){
		// TODO

		// A.kmax = Protected::DotProdBoundClassic(F,F.one);
		// A.m = rowdim;
		//       A.n = coldim;
		// A.nnz = nnz;
		// std::vector<uint64_t> rows(A.m, 0);
		// for(uint64_t i = 0 ; i < A.nnz ; ++i)
		// 	rows[row[i]]++;
		// A.maxrow = *(std::max_element(rows.begin(), rows.end()));
		// A.ld = A.maxrow;
		// for(auto & x : rows)
		// 	if(x != 0)
		// 		A.mRow++;

		// if(A.kmax > A.maxrow)
		// 	A.delayed = true;

		// A.col = fflas_new<index_t>(A.mRow*A.ld, Alignment::CACHE_LINE);
		//       A.dat = fflas_new(F, rowdim*A.ld, 1, Alignment::CACHE_LINE);

		//       for(size_t i = 0 ; i < rowdim*A.ld ; ++i){
		//       	A.col[i] = 0;
		//       	F.assign(A.dat[i], F.zero);
		//       }

		//       size_t currow = row[0], it = 0;

		//       for(size_t i = 0 ; i < nnz ; ++i){
		//       	if(row[i] != currow){
		//       		it = 0;
		//       		currow = row[i];
		//       	}
		//       	A.col[row[i]*A.ld + it] = col[i];
		//       	A.dat[row[i]*A.ld + it] = dat[i];
		//       	++it;
		//       }
	}

	template<class Field, class IndexT>
	inline void sparse_init(const Field & F, Sparse<Field, SparseMatrix_t::ELL_R_ZO> & A,
				const IndexT * row, const IndexT * col, typename Field::ConstElement_ptr dat,
				uint64_t rowdim, uint64_t coldim, uint64_t nnz){
		// TODO

		// A.kmax = Protected::DotProdBoundClassic(F,F.one);
		// A.m = rowdim;
		//       A.n = coldim;
		// A.nnz = nnz;
		// std::vector<uint64_t> rows(A.m, 0);
		// for(uint64_t i = 0 ; i < A.nnz ; ++i)
		// 	rows[row[i]]++;
		// A.maxrow = *(std::max_element(rows.begin(), rows.end()));
		// A.ld = A.maxrow;
		// if(A.kmax > A.maxrow)
		// 	A.delayed = true;

		// A.col = fflas_new<index_t>(rowdim*A.ld, Alignment::CACHE_LINE);

		//       for(size_t i = 0 ; i < rowdim*A.ld ; ++i){
		//       	A.col[i] = 0;
		//       }

		//       size_t currow = row[0], it = 0;

		//       for(size_t i = 0 ; i < nnz ; ++i){
		//       	if(row[i] != currow){
		//       		it = 0;
		//       		currow = row[i];
		//       	}
		//       	A.col[row[i]*A.ld + it] = col[i];
		//       	++it;
		//       }
	}

} // FFLAS

#endif //  __FFLASFFPACK_fflas_ELL_R_spmv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
