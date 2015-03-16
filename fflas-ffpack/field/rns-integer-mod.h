/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
/*! @file field/rns-integer-mod.h
 * @ingroup field
 * @brief  representation of <code>Z/pZ</code> using RNS representation (note: fixed precision)
 */


#ifndef __FFPACK_rns_integer_mod_H
#define __FFPACK_rns_integer_mod_H

#include <vector>
#include <cmath> 
// extern "C" {
// 	#include <quadmath.h>
// }
//using namespace std; // NO WAY! A.B. - 2014-12-18

#include <givaro/modular-integer.h>
#include <givaro/givinteger.h>

#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-integer.h"

#include "fflas-ffpack/fflas/fflas_level1.inl"
#include "fflas-ffpack/fflas/fflas_level2.inl"
#include "fflas-ffpack/fflas/fflas_level3.inl"
#include "fflas-ffpack/fflas/fflas_fscal_mp.inl"
 
#if defined(BENCH_PERF_FGEMM_MP) || defined(BENCH_PERF_TRSM_MP) || defined(BENCH_PERF_LQUP_MP)
#define BENCH_PERF_SCAL_MP
#define BENCH_MODP
#endif

namespace FFPACK {

	template<typename RNS>
	class RNSIntegerMod {
	public:
		typedef typename RNS::Element                   Element;
		typedef typename RNS::Element_ptr           Element_ptr;
		typedef typename RNS::ConstElement_ptr ConstElement_ptr;

	protected:
		typedef typename RNS::BasisElement BasisElement;
		typedef Givaro::Modular<BasisElement> ModField;
		typedef Givaro::Integer integer;
		
		integer                              _p;
		std::vector<BasisElement>       _Mi_modp_rns;
		std::vector<BasisElement>       _iM_modp_rns;
		const RNS                         *_rns;
		Givaro::Modular<Givaro::Integer>     _F;
		RNSInteger<RNS>                  _RNSdelayed;
	public:
		Element                one, mOne,zero;

#ifdef BENCH_MODP
		mutable double t_modp, t_igemm, t_scal,t_trsm; 
		mutable size_t n_modp;
#endif


		RNSIntegerMod(const integer& p, const RNS& myrns) : _p(p),
								  _Mi_modp_rns(myrns._size*myrns._size),
								  _iM_modp_rns(myrns._size*myrns._size),
								  _rns(&myrns),
								    _F(p),
								    _RNSdelayed(myrns){
			init(one,1);
			init(zero,0);
			init(mOne,-1);
			integer iM=0;
			size_t mysize=myrns._size;
			integer sum=0;
			for (size_t i=0;i<mysize;i++){
				integer Mi = myrns._Mi[i] % _p;
				for (size_t j=0;j<mysize;j++){
					_iM_modp_rns[i+j*mysize]=  iM % myrns._basis[j];
					_Mi_modp_rns[i+j*mysize]=  Mi % myrns._basis[j];
				}
				iM+=myrns._M;iM%=_p;
				sum+=myrns._basis[i];
			}
#ifdef BENCH_MODP
		t_modp=t_igemm=t_scal=t_trsm=0.;
		n_modp=0;
#endif
		}

		const rns_double& rns() const {return *_rns;}
		const RNSInteger<RNS>& delayed() const {return _RNSdelayed;}
		
		size_t size() const {return _rns->_size;}

		bool isOne(const Element& x) const {
			bool isone=true;
			for (size_t i=0;i<_rns->_size;i++)
				isone&= (one._ptr[i]== x._ptr[i]);
			return isone;
		}

		bool isMOne(const Element& x) const {
			bool ismone=true;
			for (size_t i=0;i<_rns->_size;i++)
				ismone&= (mOne._ptr[i]== x._ptr[i]);
			return ismone;
		}

		bool isZero(const Element& x) const {
			//write(std::cout,x)<<" == ";
			//write(std::cout,zero)<<std::endl;
			integer t1; 
			t1=convert(t1,x)%_p;
			//std::cout<<"t1="<<t1<<std::endl;
			bool iszero=true;
			for (size_t i=0;i<_rns->_size;i++)
				iszero&= (zero._ptr[i]==x._ptr[i]);
			//std::cout<<(iszero || (t1==integer(0))?"zero":"nonzero")<<std::endl;
			return iszero || (t1==integer(0));
		}

		integer characteristic(integer &p) const { return p=_p;}

		integer cardinality(integer &p) const { return p=_p;}

		Element& init(Element& x) const{
			if (x._ptr == NULL){
				x._ptr = FFLAS::fflas_new<BasisElement>(_rns->_size);
				x._stride=1;
				x._alloc=true;
			}
			return x;
		}
		Element& init(Element& x, const Givaro::Integer& y) const{
			init(x);
			size_t k =(_p.bitsize())/16+((_p.bitsize())%16?1:0);
			_rns->init(1,1,x._ptr,x._stride, &y,1,k);
			return x;
		}

		// assume this is the mod p operation
		Element& reduce (Element& x, const Element& y) const{
			Givaro::Integer tmp;
			convert(tmp,y);
			tmp %= _p;
			init (x,tmp);
			return x;
		}

		Element& reduce (Element& x) const{
			Givaro::Integer tmp;
			convert (tmp, x);
			tmp %= _p;
			return init (x, tmp);
		}

		Element& init(Element& x, const Element& y) const{
			return reduce (x, y);
		}

       
		Givaro::Integer convert(Givaro::Integer& x, const Element& y)const {
			_rns->convert(1,1,integer(0),&x,1,y._ptr,y._stride);
			return x;
		}

		Element& assign(Element& x, const Element& y) const {
			for(size_t i=0;i<_rns->_size;i++)
				x._ptr[i*x._stride] = y._ptr[i*y._stride];
			return x;
		}

		Element& add(Element& x, const Element& y, const Element& z) const {
			for(size_t i=0;i<_rns->_size;i++)
				_rns->_field_rns[i].add((x._ptr)[i*x._stride],
							(y._ptr)[i*y._stride],
							(z._ptr)[i*z._stride]);
			return x;
		}

		Element& sub(Element& x, const Element& y, const Element& z) const {
			for(size_t i=0;i<_rns->_size;i++)
				_rns->_field_rns[i].sub((x._ptr)[i*x._stride],
							(y._ptr)[i*y._stride],
							(z._ptr)[i*z._stride]);
			return x;
		}

		Element& neg(Element& x, const Element& y) const {
			for(size_t i=0;i<_rns->_size;i++)
				_rns->_field_rns[i].neg((x._ptr)[i*x._stride],
							(y._ptr)[i*y._stride]);						
			return x;
		}

		Element& mul(Element& x, const Element& y, const Element& z) const {
			for(size_t i=0;i<_rns->_size;i++)
				_rns->_field_rns[i].mul((x._ptr)[i*x._stride],
							(y._ptr)[i*y._stride],
							(z._ptr)[i*z._stride]);
			return x;
		}


		Element& axpyin(Element& x, const Element& y, const Element& z) const {
			for(size_t i=0;i<_rns->_size;i++)
				_rns->_field_rns[i].axpyin((x._ptr)[i*x._stride],
							   (y._ptr)[i*y._stride],
							   (z._ptr)[i*z._stride]);
			return x;
		}

		Element& inv(Element& x, const Element& y) const {
			Givaro::Integer tmp;
			convert(tmp,y);
			_F.invin(tmp);
			init(x,tmp);
			return x;
		}

		bool areEqual(const Element& x, const Element& y) const {
			for(size_t i=0;i<_rns->_size;i++)
				if (!_rns->_field_rns[i].areEqual((x._ptr)[i*x._stride],(y._ptr)[i*y._stride]))
					return false;
			return true;
		}
		std::ostream& write(std::ostream& os, const Element& y) const {
			os<<"[ "<< (long) (y._ptr)[0];
			for(size_t i=1;i<_rns->_size;i++)
				os<<" , "<< (long) ((y._ptr)[i*y._stride]);
			return os<<" ]";
		}


		std::ostream& write(std::ostream& os) const {
			os<<"M:=[ "<< (long) _rns->_basis[0];
			for(size_t i=1;i<_rns->_size;i++)
				os<<" , "<< (long) _rns->_basis[i];
			return os<<" ]"<<std::endl;
		}


		void reduce_modp(size_t n, Element_ptr B) const{
#ifdef BENCH_MODP
			FFLAS::Timer chrono; chrono.start();
#endif			
			size_t _size= _rns->_size;
			BasisElement *Gamma, *alpha, *A;
			A=B._ptr;
			size_t rda = B._stride;
			Givaro::UnparametricRing<BasisElement> D;
			Gamma = FFLAS::fflas_new(D,_size,n);
			alpha = FFLAS::fflas_new(D,n,1);

			// compute Gamma
			//for(size_t i=0;i<_size;i++)
			//	FFLAS::fscal(_rns->_field_rns[i], n, _rns->_MMi[i], A+i*rda, 1, Gamma+i*n,1);
			typename RNS::Element mmi(const_cast<typename RNS::BasisElement*>(_rns->_MMi.data()),1);
			FFLAS::fscal(_RNSdelayed, n, mmi, B, 1, typename RNS::Element_ptr(Gamma,n), 1);
			
			// compute A = _Mi_modp_rns.Gamma (note must be reduced mod m_i, but this is postpone to the end)
			FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans, _size, n, _size, D.one, _Mi_modp_rns.data(), _size, Gamma, n, D.zero, A, rda);
			
			//std::cout<<"fgemv (Y)...";
			//std::cout<<"fgemv (Y)..."<<n<<" -> "<<_size<<endl;;
			// compute alpha = _invbase.Gamma
			FFLAS::fgemv(D,FFLAS::FflasTrans, _size, n, D.one, Gamma, n, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);
			//std::cout<<"done"<<std::endl;
                        			
			// compute ((z-(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
			for(size_t i=0;i<_size;i++){
				for(size_t j=0;j<n;j++){
					//long aa=floor(alpha[j]+0.5);
					long aa= (long)rint(alpha[j]);
					A[j+i*rda]-=_iM_modp_rns[aa+i*_size];
				}
			}

			// reduce each row of A modulo m_i
			for (size_t i=0;i<_size;i++)
				FFLAS::freduce (_rns->_field_rns[i], n, A+i*rda, 1);

			FFLAS::fflas_delete(Gamma);
			FFLAS::fflas_delete(alpha);

#ifdef BENCH_MODP
			chrono.stop();
			t_modp+=chrono.usertime();
#endif
	 	}

		std::ostream& write_matrix(std::ostream& c,
					double* E,
					int n, int m, int lda) const
		{
			c<<std::endl<<"***********************"<<std::endl;
			for (int i = 0; i<n;++i){
				for (int j=0; j<m;++j)
					c << (long)*(E+j+lda*i) << " ";
				c << std::endl;
			}
			c<<"***********************"<<std::endl;
			return c << std::endl;
		}
		
		void reduce_modp(size_t m, size_t n, Element_ptr B, size_t lda) const{
#ifdef BENCH_MODP
			FFLAS::Timer chrono; chrono.start();
#endif
			//cout<<"REDUCE MOD WITH LDA!=N"<<endl;
			size_t _size= _rns->_size;
			size_t mn=m*n;
			BasisElement *Gamma, *alpha, *z, *A;
			A=B._ptr;
			size_t rda=B._stride;
			Gamma = FFLAS::fflas_new<BasisElement>(mn*_size);
			alpha = FFLAS::fflas_new<BasisElement>(mn);
			z     = FFLAS::fflas_new<BasisElement>(mn*_size);

			// compute Gamma
			//for(size_t i=0;i<_size;i++)
			//	FFLAS::fscal(_rns->_field_rns[i], m, n, _rns->_MMi[i], A+i*rda, lda, Gamma+i*mn,n);
			typename RNS::Element mmi(const_cast<typename RNS::BasisElement*>(_rns->_MMi.data()),1);
			FFLAS::fscal(_RNSdelayed, m, n, mmi, B, lda, typename RNS::Element_ptr(Gamma,mn), n);
			
			// compute Gamma = _Mi_modp_rns.Gamma (note must be reduced mod m_i, but this is postpone to the end)
			Givaro::UnparametricRing<BasisElement> D;
		
			FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,_size, mn, _size, D.one, _Mi_modp_rns.data(), _size, Gamma, mn, D.zero, z, mn);

			//write_matrix(std::cout,Gamma, mn, _size, mn);
					
			// compute alpha = _invbase.Gamma
			//std::cout<<"fgemv (X)..."<<m<<"x"<<n<<" -> "<<_size<<"  "<<lda<<endl;;
			FFLAS::fgemv(D, FFLAS::FflasTrans, _size, mn, D.one, Gamma, mn, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);
			//std::cout<<"done"<<std::endl; 

			// compute A=((Gamma--(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
 			for(size_t k=0;k<_size;k++){
				for(size_t i=0;i<m;i++)
					for(size_t j=0;j<n;j++){
						long aa=(long)floor(alpha[j+i*n]+0.5);
						A[j+i*lda+k*rda]= z[j+i*n+k*mn]-_iM_modp_rns[aa+k*_size];
					}
			}
			
			// reduce each row of A modulo m_i
			for (size_t i=0;i<_size;i++)
				FFLAS::freduce (_rns->_field_rns[i], m, n, A+i*rda, lda);
			FFLAS::fflas_delete(Gamma);
			FFLAS::fflas_delete(alpha);
			FFLAS::fflas_delete(z);
#ifdef BENCH_MODP
			chrono.stop();
			t_modp+=chrono.usertime();
#endif
		}

#ifdef __DLP_CHALLENGE		
		void reduce_modp_rnsmajor_scal(size_t n, Element_ptr B) const {
			std::cout << "modp scalar" << std::endl;
			size_t _size= _rns->_size;
            BasisElement *Gamma, *A;
            const BasisElement *Mi_modp = _Mi_modp_rns.data();
            const BasisElement *inv_basis = _rns->_invbasis.data();
			A=B._ptr;
			Givaro::UnparametricRing<BasisElement> D;
			FFLAS::Timer T;
			
            Gamma = FFLAS::fflas_new(D,n,1);
            for(size_t i = 0 ; i < n; ++i){
            	// Compute Gamma
            	for(size_t k = 0; k < _size ; ++k){
                	Gamma[k] = A[i*_size+k]*_rns->_MMi[k];
                    Gamma[k] -= std::floor(Gamma[k]*_rns->_invbasis[k])*_rns->_basis[k];
					if(Gamma[k] >= _rns->_basis[k]){
						Gamma[k] -= _rns->_basis[k]; 
					}else if(Gamma[k] < 0){
						Gamma[k] += _rns->_basis[k]; 
					}
                }
                
                // FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasTrans, n, _size, _size, D.one, Gamma, _size, _Mi_modp_rns.data(), _size, D.zero, A, _size);
                // Mul by Mi_modp
                for(size_t k = 0 ; k < _size ; ++k){
                	A[i*_size+k] = 0;
                	for(size_t j = 0 ; j < _size ; ++j){
                		A[i*_size+k] += Gamma[k] * Mi_modp[j*_size+k];
                	}
                }

                // compute alpha
                // FFLAS::fgemv(D,FFLAS::FflasNoTrans, n, _size, D.one, Gamma, _size, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);
                BasisElement alpha = 0;
                for(size_t k = 0 ; k < _size ; ++k){
                	alpha += Gamma[k]*inv_basis[k];
                }

                // -= alpha
                long aa= (long)rint(alpha);
                for(size_t k = 0; k < _size ; k++){
                	A[i*_size+k] -= _iM_modp_rns[aa+k*_size];
                }
            }
            _rns->reduce(n,A,1,true);
            FFLAS::fflas_delete(Gamma);
		}

		void reduce_modp_rnsmajor_scal_quad(size_t n, Element_ptr B) const {
			std::cout << "modp scalar quad" << std::endl;
			FFLAS::Timer T;
			size_t _size= _rns->_size;

			// A=B._ptr;
			Givaro::UnparametricRing<BasisElement> D;		
#pragma omp parallel for schedule(static, 32)
            for(size_t i = 0 ; i < n; ++i){
            	__int128_t *A;
            	BasisElement *Gamma;
            	Gamma = FFLAS::fflas_new<BasisElement>(_size);
            	A = FFLAS::fflas_new<__int128_t>(_size);
            	// Compute Gamma
            	for(size_t k = 0; k < _size ; ++k){
                	Gamma[k] = B._ptr[i*_size+k]*_rns->_MMi[k];
                    Gamma[k] -= std::floor(Gamma[k]*_rns->_invbasis[k])*_rns->_basis[k];
					if(Gamma[k] >= _rns->_basis[k]){
						Gamma[k] -= _rns->_basis[k]; 
					}else if(Gamma[k] < 0){
						Gamma[k] += _rns->_basis[k]; 
					}
                }
                // FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasTrans, n, _size, _size, D.one, Gamma, _size, _Mi_modp_rns.data(), _size, D.zero, A, _size);
                // Mul by Mi_modp
                for(size_t k = 0 ; k < _size ; ++k){
                	A[k] = 0;
                	for(size_t j = 0 ; j < _size ; ++j){
                		A[k] += (__int128_t)(Gamma[k])*(int64_t)_Mi_modp_rns[j*_size+k];
                		// A[k], mul_q(QFloat(Gamma[k], 0), _Mi_modp_rns[j*_size+k]));
                		// A[i*_size+k] += Gamma[k] * Mi_modp[j*_size+k];
                	}
                }

                // compute alpha
                // FFLAS::fgemv(D,FFLAS::FflasNoTrans, n, _size, D.one, Gamma, _size, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);
                __int128_t alpha = 0;
                for(size_t k = 0 ; k < _size ; ++k){
                	alpha += (__int128_t)Gamma[k]*(int64_t)_rns->_invbasis[k];
                }

                // -= alpha
                long aa= (long)alpha;
                for(size_t k = 0; k < _size ; k++){
                	A[k] -= (int64_t)_iM_modp_rns[aa+k*_size];
                	A[k] %= (int64_t)_rns->_basis[k];
                	B._ptr[k] = (double)((int64_t)A[k]);
                }
                FFLAS::fflas_delete(Gamma);
            	FFLAS::fflas_delete(A);
            }
            // _rns->reduce(n,A,1,true);
		}
#endif // __DLP_CHALLENGE		

		void reduce_modp_rnsmajor(size_t n, Element_ptr B) const{
			// std::cout << "modp BLAS" << std::endl;
#ifdef BENCH_MODP
                        FFLAS::Timer chrono; chrono.start();
#endif
                        size_t _size= _rns->_size;
                        BasisElement *Gamma, *alpha, *A;
			A=B._ptr;
			Givaro::UnparametricRing<BasisElement> D;
			FFLAS::Timer T;
			// T.start();
                        Gamma = FFLAS::fflas_new(D,n,_size);
                        alpha = FFLAS::fflas_new(D,n,1);
                        // T.stop();
                        // std::cout << "Alloc: " << T << std::endl;
                        // compute Gamma (NOT EFFICIENT)
                        //for(size_t i=0;i<_size;i++)
                        // 
			       // FFLAS::fscal(_rns->_field_rns[i], n, _rns->_MMi[i], A+i, _size, Gamma+i,_size);
                        T.start();
#ifdef __FFLASFFPACK_USE_SIMD
                        using simd = Simd<BasisElement>;
                        using vect_t = typename simd::vect_t;

                       if(((int64_t)A%simd::alignment == 0) && (_size%simd::vect_size==0)){
                        	auto MMi = _rns->_MMi.data();
                        	for(size_t i = 0 ; i < n ; ++i){
                        		vect_t vA1, vA2, vMi1, vMi2, tmp1, tmp2, tmp3, v, max, basis, inv, neg;
                        		size_t k = 0;
                        		for( ; k < ROUND_DOWN(_size, simd::vect_size) ; k+=simd::vect_size){
                        			basis = simd::load(_rns->_basis.data()+k);
						inv   = simd::load(_rns->_invbasis.data()+k);
						max   = simd::load(_rns->_basisMax.data()+k);
						neg   = simd::load(_rns->_negbasis.data()+k);
						vA1  = simd::load(A+i*_size+k);
                        			vMi1 = simd::load(MMi+k);
                        			v = simd::mul(vA1, vMi1);
                        			tmp1  = simd::floor(simd::mul(v, inv));				
						tmp2  = simd::fnmadd(v, tmp1, basis);
						tmp1  = simd::greater(tmp2, max);
						tmp3  = simd::lesser(tmp2, simd::zero());
						tmp1  = simd::vand(tmp1, neg);
						tmp3  = simd::vand(tmp3, basis);
						tmp1  = simd::vor(tmp1, tmp3);
						tmp2  = simd::add(tmp2, tmp1);
                        			simd::store(Gamma+i*_size+k, tmp2);
                        		}
                        	}
                        }else{
                        	vect_t vA1, vA2, vMi1, vMi2, tmp1, tmp2, tmp3, v, max, basis, inv, neg;
                        	auto MMi = _rns->_MMi.data();
                        	for(size_t i = 0 ; i < n ; ++i){
                        		size_t k = 0;
                        		for( ; k < ROUND_DOWN(_size, simd::vect_size) ; k+=simd::vect_size){
                        			basis = simd::load(_rns->_basis.data()+k);
									inv   = simd::load(_rns->_invbasis.data()+k);
									max   = simd::load(_rns->_basisMax.data()+k);
									neg   = simd::load(_rns->_negbasis.data()+k);                        			
                        			vA1 = simd::loadu(A+i*_size+k);
                        			vMi1 = simd::loadu(MMi+k);
                        			v = simd::mul(vA1, vMi1);
                        			tmp1  = simd::floor(simd::mul(v, inv));				
						tmp2  = simd::fnmadd(v, tmp1, basis);
						tmp1  = simd::greater(tmp2, max);
						tmp3  = simd::lesser(tmp2, simd::zero());
						tmp1  = simd::vand(tmp1, neg);
						tmp3  = simd::vand(tmp3, basis);
						tmp1  = simd::vor(tmp1, tmp3);
						tmp2  = simd::add(tmp2, tmp1);
                        			simd::storeu(Gamma+i*_size+k, tmp2);
                        		}
                        		for(; k < _size ; ++k){
                        			Gamma[i*_size+k] = A[i*_size+k]*MMi[k];
                        			Gamma[i*_size+k] -= std::floor(Gamma[i*_size+k]*_rns->_invbasis[k])*_rns->_basis[k];
						if(Gamma[i*_size+k] >= _rns->_basis[k]){
							Gamma[i*_size+k] -= _rns->_basis[k]; 
						}else if(Gamma[i*_size+k] < 0){
							Gamma[i*_size+k] += _rns->_basis[k]; 
						}
                        		}
                        	}
                        }
                        // _rns->reduce(n,Gamma,1,true);
#else
			typename RNS::Element mmi(const_cast<typename RNS::BasisElement*>(_rns->_MMi.data()),1);
			FFLAS::fscal(_RNSdelayed, n, mmi, B, 1, typename RNS::Element_ptr(Gamma,1), 1);
#endif	
					T.stop();
					// std::cout << "Gamma: " << T << std::endl;
			
			
                        // compute A = Gamma._Mi_modp_rns^T (note must be reduced mod m_i, but this is postpone to the end)
					T.start();
                        FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasTrans, n, _size, _size, D.one, Gamma, _size, _Mi_modp_rns.data(), _size, D.zero, A, _size);
					T.stop();
					// std::cout << "fgemm: " << T << std::endl;
                        // std::cout<<"fgemv (Y)...";
                        //std::cout<<"fgemv (Y)..."<<n<<" -> "<<_size<<endl;;
                        // compute alpha = Gamma._invbasis
					T.start();
                        FFLAS::fgemv(D,FFLAS::FflasNoTrans, n, _size, D.one, Gamma, _size, _rns->_invbasis.data(), 1 , D.zero, alpha, 1);
                        T.stop();
					// std::cout << "fgemv: " << T << std::endl;
                        //std::cout<<"done"<<std::endl;
					T.start();
                        // compute ((z-(alpha.M mod p)) mod m_i (perform the subtraction over Z and reduce at the end)
                        for(size_t j=0;j<n;j++){
                        	long aa= (long)rint(alpha[j]);
                        	for(size_t i=0;i<_size;i++){
                                        //long aa=floor(alpha[j]+0.5);
                                        A[j*_size+i]-=_iM_modp_rns[aa+i*_size];
                                }
                        }
                        // vect_t viM;
                        // auto iM = _iM_modp_rns.data();
                        // for(size_t j = 0 ; j < n ; ++j){
                        // 	long aa= (long)rint(alpha[j]);
                        // 	for(int i = 0 ; i < ROUND_DOWN(_size, simd::vect_size) ; i+=simd::vect_size){
                        // 		vA = simd::load(A+j*_size+i);
                        // 		viM = simd::load(iM+aa)
                        // 	}
                        // }
                        T.stop();
					// std::cout << "last: " << T << std::endl;
					T.start();
                        // reduce each column of A modulo m_i 
					_rns->reduce(n,A,1,true);
					T.stop();
					// std::cout << "reduce: "<< T << std::endl;

						// T.start();
                        FFLAS::fflas_delete(Gamma);
                        FFLAS::fflas_delete(alpha);
                        // T.stop();
                        // std::cout << "delete: " << T << std::endl;
#ifdef BENCH_MODP
                        chrono.stop();
                        t_modp+=chrono.usertime();
#endif

                }


	}; // end of class RNSIntegerMod

} // end of namespace FFPACK


namespace FFLAS {

	// specialization for the fflas alloc function
	template<>
	inline FFPACK::rns_double_elt_ptr
	fflas_new(const FFPACK::RNSIntegerMod<FFPACK::rns_double> &F,  const size_t m,  const size_t n, const Alignment align){
		return FFPACK::rns_double_elt_ptr(FFLAS::fflas_new<double>(m*n*F.size(),align),m*n);
	}

	// function to convert from integer to RNS (note: this is not the finit function from FFLAS, extra k)
	template<typename RNS>
	void finit_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n, size_t k,
		   const Givaro::Integer *B, const size_t ldb, typename RNS::Element_ptr A)
	{
			F.rns().init(m,n,A._ptr,A._stride, B,ldb,k);
	}
	template<typename RNS>
	void finit_trans_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n, size_t k,
			     const Givaro::Integer *B, const size_t ldb, typename RNS::Element_ptr A)
	{
		F.rns().init_transpose(m,n,A._ptr,A._stride, B,ldb,k);
	}

	// function to convert from RNS to integer (note: this is not the fconvert function from FFLAS, extra alpha)
	template<typename RNS>
	void fconvert_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n,
		      Givaro::Integer alpha, Givaro::Integer *B, const size_t ldb, typename RNS::ConstElement_ptr A)
	{
		F.rns().convert(m,n,alpha,B,ldb,A._ptr,A._stride);
	}
	template<typename RNS>
	void fconvert_trans_rns(const FFPACK::RNSIntegerMod<RNS> &F, const size_t m, const size_t n,
				Givaro::Integer alpha, Givaro::Integer *B, const size_t ldb, typename RNS::ConstElement_ptr A)
	{
		F.rns().convert_transpose(m,n,alpha,B,ldb,A._ptr,A._stride);
	}

} // end of namespace FFLAS

#endif
