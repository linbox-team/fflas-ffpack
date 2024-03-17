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


#ifndef __FFLASFFPACK_field_rns_double_INL
#define __FFLASFFPACK_field_rns_double_INL

#include "fflas-ffpack/fflas/fflas_freduce.h"
#include "fflas-ffpack/field/rns-common.h"


  // Arns must be an array of m*n*_size
  // abs(||A||) < 2^(16k)
  inline void rns_double::init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR, const FFLAS::FFLAS_TRANSPOSE trans) const
  {
    if (k>_ldm){
      FFPACK::failure()(__func__,__FILE__,__LINE__,"rns_double [init] -> rns basis is too small to handle integers with 2^(16*k) values ");
      std::cerr<<"with k="<<k<<" _ldm="<<_ldm<<std::endl;
    }
    const size_t mn=m*n;
    Givaro::ZRing<double> ZD;

    if (mn) {
      // split A into A_beta according to a Kronecker transform in base 2^16
      double *A_beta = FFLAS::fflas_new<double >(mn*k);
      RNS_COMMON::Kronecker_base16(trans,m,n,A_beta,A,lda,k);

      // Using Helper for potential parallelism -> need to be activated by hand
      //using ParallelStrategy= FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>;
      using ParallelStrategy= FFLAS::ParSeqHelper::Sequential;
      FFLAS::MMHelper<Givaro::ZRing<double>, FFLAS::MMHelperAlgo::Winograd>  MatMulHelper (ZD, -1, ParallelStrategy());

      // Pseudo reduction of the matrix entries
      if (RNS_MAJOR==false) 
        // Arns = _crt_in x A_beta^T
        FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in.data(),_ldm,A_beta,k,0.,Arns,mn,MatMulHelper);
      else 
        // Arns =  A_beta x _crt_in^T
        FFLAS::fgemm (ZD, FFLAS::FflasNoTrans,FFLAS::FflasTrans,mn,_size,k,1.0,A_beta, k, _crt_in.data(),_ldm,0.,Arns,_size,MatMulHelper);

      // Reduce the pseudo reduce matrix modulo RNS basis
      reduce(mn,Arns,rda,RNS_MAJOR);

#ifdef CHECK_RNS // check the result
      RNS_COMMON::check_integers_to_rns(trans, m, n, Arns, rda, A, lda, _basis, RNS_MAJOR);
#endif
      FFLAS::fflas_delete(A_beta);
    }
  }


  inline void rns_double::convert(size_t m, size_t n, integer gamma, integer* A, size_t lda,
                                  const double* Arns, size_t rda, bool RNS_MAJOR, const FFLAS::FFLAS_TRANSPOSE trans) const
  {
    const size_t  mn= m*n;
    Givaro::ZRing<double> ZD;
    Givaro::ZRing<Givaro::Integer> ZZ;
    if (mn) {
#ifdef CHECK_RNS
      integer* Acopy=new integer[m*n];
      FFLAS::fassign(ZZ,m,n,Acopy,n,A,lda);
#endif 

      double *A_beta= FFLAS::fflas_new<double>(mn*_ldm);

      // Using Helper for potential parallelism -> need to be activated by hand
      //using ParallelStrategy= FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>;
      using ParallelStrategy= FFLAS::ParSeqHelper::Sequential;
      FFLAS::MMHelper<Givaro::ZRing<double>, FFLAS::MMHelperAlgo::Winograd>  MatMulHelper (ZD, -1, ParallelStrategy());

      if (RNS_MAJOR==false) // compute A_beta = Ap^T x M_beta
        FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans, mn, _ldm, _size, 1.0 , Arns, rda, _crt_out.data(), _ldm, 0., A_beta,_ldm, MatMulHelper);
      else  // compute A_beta = Ap x M_Beta
        FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mn, _ldm, _size, 1.0 , Arns, _size, _crt_out.data(), _ldm, 0., A_beta, _ldm, MatMulHelper);


      RNS_COMMON::inverse_Kronecker_base16(trans,m,n,A_beta,A,lda,_ldm, _M, gamma);

      //FFLAS::WriteMatrix(std::cout,Givaro::ZRing<Givaro::Integer>(),m,n, A, lda, FFLAS::FflasSageMath)<<std::endl;
      //FFLAS::WriteMatrix(std::cout,Givaro::ZRing<double>(),1,_size, _basis.data(), _size, FFLAS::FflasSageMath)<<std::endl;      
      // if (RNS_MAJOR)
      //   FFLAS::WriteMatrix(std::cout,Givaro::ZRing<double>(),mn,_size, Arns, rda,FFLAS::FflasSageMath)<<std::endl;
      // else
      //   FFLAS::WriteMatrix(std::cout,Givaro::ZRing<double>(),_size, mn, Arns, rda,FFLAS::FflasSageMath)<<std::endl;

      FFLAS::fflas_delete(A_beta);
#ifdef CHECK_RNS
      RNS_COMMON::check_rns_to_integers(trans, m, n, Arns, rda, A, lda, gamma, Acopy, n, _basis, RNS_MAJOR);
      FFLAS::fflas_delete(Acopy);
#endif

}
  }



  // reduce entries of Arns to be less than the rns basis elements
  inline void rns_double::reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR) const{

    if (RNS_MAJOR) {
#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
      using simd = Simd<double>;
      using vect_t = typename simd::vect_t;

      if(_size % simd::vect_size == 0){
        for(size_t i = 0 ; i < n ; i++){
          vect_t tmp1, tmp2, tmp3, v, max, basis, inv, neg;
          for(size_t j = 0 ; j < _size ; j+=simd::vect_size){
            basis = simd::load(_basis.data()+j);
            inv   = simd::load(_invbasis.data()+j);
            max   = simd::load(_basisMax.data()+j);
            neg   = simd::load(_negbasis.data()+j);
            v     = simd::load(Arns+i*_size+j);
            tmp1  = simd::floor(simd::mul(v, inv));
            tmp2  = simd::fnmadd(v, tmp1, basis);
            tmp1  = simd::greater(tmp2, max); 
            tmp3  = simd::lesser(tmp2, simd::zero());
            tmp1  = simd::vand(tmp1, neg);
            tmp3  = simd::vand(tmp3, basis);
            tmp1  = simd::vor(tmp1, tmp3);
            tmp2  = simd::add(tmp2, tmp1);
            simd::store(Arns+i*_size+j, tmp2);
          }
        }
      }else{
        for(size_t i = 0 ; i < n ; i++){
          vect_t tmp1, tmp2, tmp3, v, max, basis, inv, neg;
          size_t j = 0;
          for( ; j < ROUND_DOWN(_size, simd::vect_size) ; j+=simd::vect_size){
            basis = simd::load(_basis.data()+j);
            inv   = simd::load(_invbasis.data()+j);
            max   = simd::load(_basisMax.data()+j);
            neg   = simd::load(_negbasis.data()+j);
            v     = simd::loadu(Arns+i*_size+j);
            tmp1  = simd::floor(simd::mul(v, inv));
            tmp2  = simd::fnmadd(v, tmp1, basis);
            tmp1  = simd::greater(tmp2, max);
            tmp3  = simd::lesser(tmp2, simd::zero());
            tmp1  = simd::vand(tmp1, neg);
            tmp3  = simd::vand(tmp3, basis);
            tmp1  = simd::vor(tmp1, tmp3);
            tmp2  = simd::add(tmp2, tmp1);
            simd::storeu(Arns+i*_size+j, tmp2);
          }
          for( ; j < _size ; ++j){
            // std::cout << j << std::endl;
            // auto x = std::floor(Arns[i*_size+j] * _invbasis[j]);
            Arns[i*_size+j] -= std::floor(Arns[i*_size+j]*_invbasis[j])*_basis[j];
            // Arns[i*_size+j] = std::fma(Arns[i*_size+j], -x, _basis[j]);
            if(Arns[i*_size+j] >= _basis[j]){
              Arns[i*_size+j] -= _basis[j];
            }else if(Arns[i*_size+j] < 0){
              Arns[i*_size+j] += _basis[j];
            }
          }
        }
      }
#else
      for(size_t i = 0 ; i < n ; i++){
        for(size_t j = 0 ; j < _size ; ++j){
          //_field_rns.reduce(Arns+i*_size+j);
          _field_rns[j].reduce(Arns[i*_size+j]);
        }
      }
#endif
    }
    else { // NOT IN RNS MAJOR
      // #ifndef __FFLASFFPACK_SEQUENTIAL
      // 			auto sp=SPLITTER(MAX_THREADS);
      // #else
      // 			auto sp=SPLITTER(1);
      // #endif
      PARFOR1D(i,_size,SPLITTER(NUM_THREADS),
               //for(size_t i=0;i<_size;i++)
               FFLAS::freduce (_field_rns[i],n,Arns+i*rda,1);
               );
    }

  }


} // FFPACK

#endif // __FFLASFFPACK_field_rns_double_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
 
