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

/** @file fflas/fflas_sparse_fgemv.inl
*/

#ifndef __FFLASFFPACK_fflas_fflas_sparse_fgemv_INL
#define __FFLASFFPACK_fflas_fflas_sparse_fgemv_INL

namespace FFLAS {
	template<class Element>
	struct VECT {
		size_t m ;
		size_t inc ;
		Element * dat ;
	}

	template<class Element>
	struct DNS {

		size_t n ;
		size_t ld ;
		Element * dat ;
	};

	template<class Element>
	struct CSR {
		size_t m ;
		size_t n ;
		size_t  * st  ;
		size_t  * col ;
		Element * dat ;
		// int mc ;
		// int ml ;
	};

	template<class Element>
	struct ELL {
		size_t m ;
		size_t n ;
		size  ld ;
		size_t  * col ;
		Element * dat ;
	};

	template<class Element>
	struct ELLR {
		size_t m ;
		size_t n ;
		size  ld ;
		size_t  * row ;
		size_t  * col ;
		Element * dat ;
	};

	template<class Element>
	struct SPADD {
		size_t ncsr;
		CSR * csr;
		size_t ncoo;
		COO * coo;
		size_t ndns;
		DNS * dns;
		size_t nell;
		ELL * ell;
		size_t nellr ;
		ELLR * ellr ;
		size_t ndia ;
		DIA * dia;

		SPADD() :
			ncsr(0)  ,csr(NULL)
			,ncoo(0) ,coo(NULL)
			,ndns(0) ,dns(NULL)
			,ndia(0) ,dia(NULL)
			,nell(0) ,ell(NULL)
			,nellr(0),ellr(NULL)
		{}
	};

	template<class Field>
	void sp_fgemv(
		      const Field& F,
		      const FFLAS_TRANSPOSE tA,
		      CSR & A,
		      const VECT<typename Field::Element> & x,
		      VECT<typename Field::Element> & y
		     )
	{
		if (tA == FflasNoTrans) {
			for (size_t i = 0 ; i < M.m ; ++i) {
				F.assign(y.dat[i],F.zero);
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					F.axpyin(y.dat[i],M.dat[j],x.dat[M.col[j]]);
			}
		}
		else {
			for (size_t i = 0 ; i < M.m ; ++i) {
				F.assign(y.dat[i],F.zero);
			}
			for (size_t i = 0 ; i < M.m ; ++i) {
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					F.axpyin(y.dat[M.col[j]],M.dat[j],x.dat[j]);
			}
		}

	}

	template<>
	void sp_fgemv(
		      const DoubleField & F,
		      const FFLAS_TRANSPOSE tA,
		      CSR & A,
		      const VECT<double> & x,
		      VECT<double> & y
		     )
	{
#ifdef __FFLASFFPACK_HAVE_MKL
		char * transa = (ta==FflasNoTrans)?'n':'t';
		mkl_cspblas_dcsrgemv (ransa, &M.m, M.dat, M.st , M.col, x.dat, y.dat);
#else
		if (tA == FflasNoTrans) {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i] = 0;
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[i] += M.dat[j] * x.dat[M.col[j]];
			}
		}
		else {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i]=0;
			}
			for (size_t i = 0 ; i < M.m ; ++i) {
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[M.col[j]] += M.dat[j] * x.dat[j];
			}
		}

#endif
	}

	template<>
	void sp_fgemv(
		      const FloatField & F,
		      const FFLAS_TRANSPOSE tA,
		      CSR & A,
		      const VECT<float> & x,
		      VECT<float> & y
		     )
	{
#ifdef __FFLASFFPACK_HAVE_MKL
		char * transa = (ta==FflasNoTrans)?'n':'t';
		mkl_cspblas_scsrgemv (ransa, &M.m, M.dat, M.st , M.col, x.dat, y.dat);
#else
		if (tA == FflasNoTrans) {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i] = 0;
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[i] += M.dat[j] * x.dat[M.col[j]];
			}
		}
		else {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i]=0;
			}
			for (size_t i = 0 ; i < M.m ; ++i) {
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[M.col[j]] += M.dat[j] * x.dat[j];
			}
		}

#endif
	}

	template<>
	void sp_fgemv(
		      const Modular<double>& F,
		      const FFLAS_TRANSPOSE tA,
		      CSR & A,
		      const VECT<double> & x,
		      VECT<double> & y
		     )
	{
		int kmax = DotProdBoundClassic(F,F.one,FflasDouble) ;

		if (tA == FflasNoTrans) {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i] = 0;
				size_t j = M.st[i];
				size_t j_loc = j;
				size_t j_end = M.st[i+1];
				size_t block = (j_end - j_loc)/kmax ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
					for ( ; j < j_loc ; ++j) {
						y.dat[i] += M.dat[j] * x.dat[M.col[j]];
					}
					F.init(y.dat[i]);
				}
				for ( ; j < j_end ; ++j) {
					y.dat[i] += M.dat[j] * x.dat[M.col[j]];
				}
				F.init(y.dat[i]);
			}
		}
		else {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i]=0;
			}
			for (size_t i = 0 ; i < M.m ; ++i) {
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[M.col[j]] += M.dat[j] * x.dat[j];
			}
		}
	}

	template<>
	void sp_fgemv(
		      const Modular<float>& F,
		      const FFLAS_TRANSPOSE tA,
		      CSR & A,
		      const VECT<float> & x,
		      VECT<float> & y
		     )
	{
		int kmax = DotProdBoundClassic(F,F.one,FflasFloat) ;

		if (tA == FflasNoTrans) {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i] = 0;
				size_t j = M.st[i];
				size_t j_loc = j;
				size_t j_end = M.st[i+1];
				size_t block = (j_end - j_loc)/kmax ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
					for ( ; j < j_loc ; ++j) {
						y.dat[i] += M.dat[j] * x.dat[M.col[j]];
					}
					F.init(y.dat[i]);
				}
				for ( ; j < j_end ; ++j) {
					y.dat[i] += M.dat[j] * x.dat[M.col[j]];
				}
				F.init(y.dat[i]);
			}
		}
		else {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i]=0;
			}
			for (size_t i = 0 ; i < M.m ; ++i) {
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[M.col[j]] += M.dat[j] * x.dat[j];
			}
		}
	}

	template<>
	void sp_fgemv(
		      const ModularBalanced<double>& F,
		      const FFLAS_TRANSPOSE tA,
		      CSR & A,
		      const VECT<double> & x,
		      VECT<double> & y
		     )
	{
		int kmax = DotProdBoundClassic(F,F.one,FflasDouble) ;

		if (tA == FflasNoTrans) {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i] = 0;
				size_t j = M.st[i];
				size_t j_loc = j;
				size_t j_end = M.st[i+1];
				size_t block = (j_end - j_loc)/kmax ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
					for ( ; j < j_loc ; ++j) {
						y.dat[i] += M.dat[j] * x.dat[M.col[j]];
					}
					F.init(y.dat[i]);
				}
				for ( ; j < j_end ; ++j) {
					y.dat[i] += M.dat[j] * x.dat[M.col[j]];
				}
				F.init(y.dat[i]);
			}
		}
		else {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i]=0;
			}
			for (size_t i = 0 ; i < M.m ; ++i) {
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[M.col[j]] += M.dat[j] * x.dat[j];
			}
		}
	}

	template<>
	void sp_fgemv(
		      const ModularBalanced<float>& F,
		      const FFLAS_TRANSPOSE tA,
		      CSR & A,
		      const VECT<float> & x,
		      VECT<float> & y
		     )
	{
		int kmax = DotProdBoundClassic(F,F.one,FflasFloat) ;

		if (tA == FflasNoTrans) {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i] = 0;
				size_t j = M.st[i];
				size_t j_loc = j;
				size_t j_end = M.st[i+1];
				size_t block = (j_end - j_loc)/kmax ;
				for (size_t l = 0 ; l < block ; ++l) {
					j_loc += block ;
					for ( ; j < j_loc ; ++j) {
						y.dat[i] += M.dat[j] * x.dat[M.col[j]];
					}
					F.init(y.dat[i]);
				}
				for ( ; j < j_end ; ++j) {
					y.dat[i] += M.dat[j] * x.dat[M.col[j]];
				}
				F.init(y.dat[i]);
			}
		}
		else {
			for (size_t i = 0 ; i < M.m ; ++i) {
				y.dat[i]=0;
			}
			for (size_t i = 0 ; i < M.m ; ++i) {
				for (size_t j = M.st[i] ; j < M.st[i+1] ; ++j)
					y.dat[M.col[j]] += M.dat[j] * x.dat[j];
			}
		}
	}


} // FFLAS

#endif // __FFLASFFPACK_fflas_fflas_sparse_fgemv_INL
