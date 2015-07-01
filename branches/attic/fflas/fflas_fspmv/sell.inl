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

namespace FFLAS {	namespace sell_details{

		// y = A x + b y ; (generic)
		template<class Field>
		inline void fspmv(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t chunk,
			      const uint32_t nChunks,
			      const uint32_t * col,
			      const uint64_t * ptr,
			      const uint32_t * chs,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      typename Field::Element * y,
			      FieldCategories::GenericTag
			     )
		{
			for(size_t i = 0 ; i < nChunks ; ++i){
				for(size_t j = 0 ; j < chs[i] ; ++j){
					size_t start = ptr[i];
					for(size_t k = 0 ; k < chunk ; ++k){
						F.axpyin(y[i*chunk], dat[start+j*chunk+k], x[col[start+j*chunk+k]]);
					}
				}
			}
		}

		template<class Field>
		inline void fspmv(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t chunk,
			      const uint32_t nChunks,
			      const uint32_t * col,
			      const uint64_t * ptr,
			      const uint32_t * chs,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      typename Field::Element * y,
			      FieldCategories::UnparametricTag
			     )
		{
			using simd = Simd<typename Field::Element>;
			using vect_t = typename simd::vect_t;

			vect_t vx, vy, vdat;

			for(size_t i = 0 ; i < nChunks ; ++i){
				vy = simd::load(y+i*chunk);
				for(size_t j = 0 ; j < chs[i] ; ++j){
					size_t start = ptr[i];
					vdat = simd::load(dat+start+j*chunk);
					vx = simd::gather(x, col+start+j*chunk);
					vy = simd::fmadd(vy, vx, vdat);
				}
				simd::stream(y+i*chunk, vy);
			}
		}

		// y = A x + b y ; (generic)
		template<class Field, bool add>
		inline void fspmv_zo(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t chunk,
			      const uint32_t nChunks,
			      const uint32_t * col,
			      const uint64_t * ptr,
			      const uint32_t * chs,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      typename Field::Element * y,
			      FieldCategories::GenericTag
			     )
		{
			if(add){
				for(size_t i = 0 ; i < nChunks ; ++i){
					for(size_t j = 0 ; j < chs[i] ; ++j){
						size_t start = ptr[i];
						for(size_t k = 0 ; k < chunk ; ++k){
							F.addin(y[i*chunk], x[col[start+j*chunk+k]]);
						}
					}
				}
			}else{
				for(size_t i = 0 ; i < nChunks ; ++i){
					for(size_t j = 0 ; j < chs[i] ; ++j){
						size_t start = ptr[i];
						for(size_t k = 0 ; k < chunk ; ++k){
							F.subin(y[i*chunk], x[col[start+j*chunk+k]]);
						}
					}
				}
			}
		}

		template<class Field, bool add>
		inline void fspmv_zo(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t chunk,
			      const uint32_t nChunks,
			      const uint32_t * col,
			      const uint64_t * ptr,
			      const uint32_t * chs,
			      const typename Field::Element *  dat,
			      const typename Field::Element * x ,
			      typename Field::Element * y,
			      FieldCategories::UnparametricTag
			     )
		{
			using simd = Simd<typename Field::Element>;
			using vect_t = typename simd::vect_t;

			vect_t vx, vy, vdat;
			if(add)
			{
				for(size_t i = 0 ; i < nChunks ; ++i){
					vy = simd::load(y+i*chunk);
					for(size_t j = 0 ; j < chs[i] ; ++j){
						size_t start = ptr[i];
						vdat = simd::load(dat+start+j*chunk);
						vx = simd::gather(x, col+start+j*chunk);
						vy = simd::add(vy, vx);
					}
					simd::store(y+i*chunk, vy);
				}
			}else{
				for(size_t i = 0 ; i < nChunks ; ++i){
					vy = simd::load(y+i*chunk);
					for(size_t j = 0 ; j < chs[i] ; ++j){
						size_t start = ptr[i];
						vdat = simd::load(dat+start+j*chunk);
						vx = simd::gather(x, col+start+j*chunk);
						vy = simd::sub(vy, vx);
					}
					simd::store(y+i*chunk, vy);
				}
			}
		}


		// delayed by kmax
		template<class Field >
		inline void fspmv(
			      const Field& F,
			      const size_t m,
			      const size_t n,
			      const size_t chunk,
			      const size_t nChunks,
			      const uint32_t * col,
			      const uint64_t * ptr,
			      const uint32_t * chs,
			      const typename Field::Element * dat,
			      const typename Field::Element * x ,
			      typename Field::Element * y,
			      const index_t & kmax,
			      FieldCategories::ModularFloatingPointTag
			      )
		{

		}

	} // details

	/* ******* */
	/* SELL_sub */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field& F,
			     const SELL_sub<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, y.m, b, y.dat, typename FieldTraits<Field>::category());
		fspmv( F, A, x, y, typename FieldTraits<Field>::category());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const SELL_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag)
	{
		sell_details::fspmv(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
		freduce (F,A.m,y.dat,1);
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const SELL_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag)
	{
		sell_details::fspmv(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	inline void fspmv(const Field & F,
			     const SELL_sub<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag)
	{
		sell_details::fspmv(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}


	/* ******* */
	/* SELL_ZO  */
	/* ******* */

	// y = A x + b y ; (generic)
	// it is supposed that no reduction is needed.
	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const SELL_ZO<Field> & A,
			     const VECT<Field> & x,
			     const typename Field::Element & b,
			     VECT<Field> & y
			    )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F, A, x, y, typename FieldTraits<Field>::category() );
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const SELL_ZO<Field > & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::GenericTag
			    )
	{
		if (A.cst == F.one) {
			sell_details::fspmv<Field, true>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else if (A.cst == F.mOne) {
			sell_details::fspmv<Field, false>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::GenericTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			sell_details::fspmv<Field, true>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::GenericTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const SELL_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::UnparametricTag
			    )
	{
		if (A.cst == F.one) {
			sell_details::fspmv<Field, true>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else if (A.cst == F.mOne) {
			sell_details::fspmv<Field, false>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::UnparametricTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			sell_details::fspmv<Field, true>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::UnparametricTag());
			fflas_delete(x1);
		}
	}

	template<class Field>
	inline void fspmv(
			     const Field & F,
			     const ELL_ZO<Field> & A,
			     const VECT<Field> & x,
			     VECT<Field> & y,
			     FieldCategories::ModularTag
			    )
	{
		if (A.cst == F.one) {
			sell_details::fspmv<Field, true>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::ModularTag());
		}
		else if (A.cst == F.mOne) {
			sell_details::fspmv<Field, false>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::ModularTag());
		}
		else {
			auto x1 = fflas_new(F, A.n, 1, Alignment::CACHE_LINE);
			fscal(F, A.n, A.cst, x.dat, 1, x1, 1);
			sell_details::fspmv<Field, true>(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,x.dat,y.dat, FieldCategories::ModularTag());
			fflas_delete(x1);
		}
		freduce (F,A.m,y.dat,1);
	}

	/* *** */
	/* SELL */
	/* *** */

	// y = A x + b y ; (generic)
	// reductions are delayed.
	template<class Field>
	void fspmv(
		      const Field& F,
		      const SELL<Field> & A,
		      const VECT<Field> & x,
		      const typename Field::Element & b,
		      VECT<Field> & y
		     )
	{
		details::init_y(F, y.m, b, y.dat,  typename FieldTraits<Field>::category());
		fspmv(F,A,x,y,typename FieldTraits<Field>::category());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const SELL<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::GenericTag
		     )
	{
		sell_details::fspmv(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,A.dat,x.dat,y.dat, FieldCategories::GenericTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const SELL<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::UnparametricTag
		     )
	{
		sell_details::fspmv(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,A.dat,x.dat,y.dat, FieldCategories::UnparametricTag());
	}

	template<class Field>
	void fspmv(
		      const Field& F,
		      const SELL<Field> & A,
		      const VECT<Field> & x,
		      VECT<Field> & y,
		      FieldCategories::ModularTag
		     )
	{
		size_t kmax = Protected::DotProdBoundClassic(F,F.one) ;
		sell_details::fspmv(F,A.m,A.n,A.chunk,A.nChunks,A.col,A.ptr,A.chs,A.dat,x.dat,y.dat, kmax);
	}

namespace details{
	struct Infos{
		size_t begin = 0;
		size_t size = 0;
		size_t rowIndex = 0;
	};

    template<class Element, class Idx>
    struct Coo
    {
        Element val = 0;
        Idx col = 0;
        Idx row = 0;
    };
}
	template<class Field, class ValT, class IdxT>
	inline void sp_sell_from_coo(const Field & F, size_t rowdim, size_t coldim, size_t nnz, ValT * COO_val, IdxT * COO_col, IdxT * COO_row, size_t chunk, size_t nbChunk)
	{
		using namespace details;
		using Coo = details::COO<ValT, IdxT>;

		std::vector<Coo> data(nnz);

		for(size_t i = 0 ; i < nnz ; ++i)
		{
			data[i].val = COO_val[i];
			data[i].col = COO_col[i];
			data[i].row = COO_row[i];
		}

		// sort nnz by row major
		std::sort(data.begin(), data.end(), [](const Coo & a, const Coo & b)
			{
				return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col));
			});

		size_t currow = 0; // current row
		Ctn::vector<Infos> rowsInfos(Rowdim);
        rowsInfos[0].begin = 0;

        for(size_t i = 0, end = data.size() ; i < end ; ++i)
		{
			if(data[i].row > currow)
			{
				rowsInfos[data[i].row].begin = i;
				currow = data[i].row;
			}
			rowsInfos[data[i].row].size++;
		}

        for(size_t i = 0 ; i < Rowdim; ++i){
            rowsInfos[i].rowIndex = i;
        }

		std::sort(rowsInfos.begin(), rowsInfos.end(), [](const Infos & a, const Infos & b)
		{
			return a.size > b.size;
		});

        // compute row permutation
        perm.resize(Rowdim);
        perm = fflas_new<uint64_t>(Rowdim);
        for(i = 0 ; i < Rowdim ; ++i)
            perm[i] = rowsInfos[i].rowIndex;

		// max row size by chunk

   		chs.resize(Rowdim/Chunk);
   		chs = fflas_new<uint32_t>(Rowdim/Chunk+1);
		i = 0;
		for(; i < Rowdim-Chunk ; i+=Chunk)
		{
			chs[i/Chunk] = (*(std::max_element(rowsInfos.begin()+i, rowsInfos.begin()+i+1, [](const Infos & a, const Infos & b)
			{
				return a.size < b.size;
			}))).size;
		}
        std::cout <<  Rowdim << " " << Chunk << " " << i << std::endl;
		if(i < Rowdim - 1)
		{
			chs[Rowdim/Chunk] = (*(std::max_element(rowsInfos.begin()+i, rowsInfos.end(), [](const Infos & a, const Infos & b)
			{
				return a.size < b.size;
			}))).size);
		}
		else
		{
			chs[Rowdim/Chunk] = 0;
		}

        cout << "kmax : " << (*(std::max_element(rowsInfos.begin(), rowsInfos.end(), [](Infos & a, Infos & b){
                return a.size < b.size;
            }))).size << endl;

        // compute the size of the vetor of indexes
        size_t vectorIndexesSize = 0;
        for(size_t i = 0 ; i < nbChunk ; ++i)
        {
        	vectorIndexesSize += chs[i]*Chunk;
        }

        // store indexes
//        idx_.resize(vectorIndexesSize);
        i = 0;
        for(; i < Rowdim ; i+=Chunk)
        {
            for(size_t j = 0 ;  j < chs_[i/Chunk] ; ++j)
            {
                for(size_t k = 0 ; k < Chunk ; ++k)
                {
                    if(j > rowsInfos[i+k].size)
                    {
                        idx_.push_back(0);
                    }
                    else
                    {
                        idx_.push_back(data[rowsInfos[i+k].begin+j].col);
                    }
                }
            }
        }
        if(i < Rowdim - 1)
        {
            for(size_t j = 0 ; j < chs_[Rowdim+1] ; ++j)
            {
                for(size_t k = 0 ; k < Rowdim - i ; ++k)
                {
                    if(j > rowsInfos[i+k].size)
                    {
                        idx_.push_back(0);
                    }
                    else
                    {
                        idx_.push_back(data[rowsInfos[i+k].begin+j].col);
                    }
                }
            }
        }
        idx_.shrink_to_fit();
        std::cout << "idx_ size " << idx_.size() << std::endl;

        // store pointers to chunks
        ptr_.resize(chs_.size()+1);
        ptr_[0] = 0;
        for(i = 1 ; i < ptr_.size() ; ++i)
        {
            ptr_[i] = ptr_[i-1]+chs_[i]*Chunk;
        }

#ifdef DEBUG
        for(i = 0 ; i < Chunk ; ++i)
        {
            if(i % Sigma == 0)
            {
                std::cout << std::endl;
                std::cout << "----------- Sigma ----------- ";
            }
            if(i % Chunk == 0)
            {
                std::cout << std::endl;
                std::cout << "----------- Chunk ----------- " << std::endl;
            }
            std::cout << rowsInfos[i].rowIndex << " (" << rowsInfos[i].size << ") : ";
            for(size_t j = rowsInfos[i].begin ; j < rowsInfos[i].begin+rowsInfos[i].size ; ++j)
                std:cout << data[j].col << " ";
            std::cout << std::endl;
        }
        for(i = 0 ; i < ptr_.size() ; ++i)
        {
            auto start = ptr_[i];
            for(size_t k = 0 ; k < Chunk ; ++k)
            {
                for(size_t j = 0 ; j < chs_[i] ; ++j)
                {
                    std::cout << idx_[start+k+j*Chunk] << " ";
                }
                std::cout << std::endl;
            }
        }
#endif
	}


} /* FFLAS */

#endif // __FFLASFFPACK_fflas_fflas_spmv_sell_INL


