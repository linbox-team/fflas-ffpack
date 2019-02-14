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

namespace FFPACK {

    // Arns must be an array of m*n*_size
    // abs(||A||) < 2^(16k)
    inline void rns_double::init(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR) const
    {
        if (k>_ldm){
            FFPACK::failure()(__func__,__FILE__,__LINE__,"rns_double [init] -> rns basis is too small to handle integers with 2^(16*k) values ");
            std::cerr<<"with k="<<k<<" _ldm="<<_ldm<<std::endl;
        }
        const size_t mn=m*n;
        if (mn) {
        double *A_beta = FFLAS::fflas_new<double >(mn*k);
        const integer* Aiter=A;
        // split A into A_beta according to a Kronecker transform in base 2^16
        //		auto sp=SPLITTER(MAX_THREADS,FFLAS::CuttingStrategy::Column,FFLAS::StrategyParameter::Threads);

        Givaro::Timer tkr; tkr.start();
        // #ifndef __FFLASFFPACK_SEQUENTIAL
        // 			auto sp=SPLITTER(MAX_THREADS);
        // #else
        // 			auto sp=SPLITTER(1);
        // #endif
        // FOR2D(i,j,m,n,sp,
        //       TASK(MODE(READ(Aiter[0]) READWRITE(A_beta[0])),
        //for(size_t i=0;i<m;i++)
        //PAR_BLOCK{
        //			FOR1D(i,m,sp,
        PARFOR1D(i,m,SPLITTER(NUM_THREADS),

                 for(size_t j=0;j<n;j++){
                 size_t idx=j+i*n;
                 const mpz_t*    m0     = reinterpret_cast<const mpz_t*>(Aiter+j+i*lda);
                 const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(m0[0]->_mp_d);
                 size_t l=0;
                 //size_t maxs=std::min(k,(Aiter[j+i*lda].size())<<2);
                 size_t maxs=std::min(k,(Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2);// to ensure 32 bits portability

#ifdef __FFLASFFPACK_HAVE_LITTLE_ENDIAN
                 if (m0[0]->_mp_size >= 0)
                 for (;l<maxs;l++)
                 A_beta[l+idx*k]=  m0_ptr[l];
                 else
                 for (;l<maxs;l++)
                 A_beta[l+idx*k]= - double(m0_ptr[l]);
#else
                 if (m0[0]->_mp_size >= 0)
                 for (;l<maxs;l++) {
                 size_t l2 = l ^ ((sizeof(mp_limb_t)/2U) - 1U);
                 A_beta[l+idx*k]=  m0_ptr[l2];
                 }
                 else
                     for (;l<maxs;l++) {
                         size_t l2 = l ^ ((sizeof(mp_limb_t)/2U) - 1U);
                         A_beta[l+idx*k]= - double(m0_ptr[l2]);
                     }
#endif
                 for (;l<k;l++)
                     A_beta[l+idx*k]=  0.;

                 // 	   );
                 }
        );

        tkr.stop();
        //if(m>1 && n>1) std::cerr<<"Kronecker : "<<tkr.realtime()<<std::endl;
        if (RNS_MAJOR==false) {
            // Arns = _crt_in x A_beta^T
            Givaro::Timer tfgemm; tfgemm.start();
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in.data(),_ldm,A_beta,k,0.,Arns,rda,
                          FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)_size,(int)mn,(int)k,1.0,_crt_in.data(),(int)_ldm,A_beta,(int)k,0.,Arns,(int)rda);
#endif
            tfgemm.stop();
            //if(m>1 && n>1) 	std::cerr<<"fgemm : "<<tfgemm.realtime()<<std::endl;
        }
        else {
            // Arns =  A_beta x _crt_in^T
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,mn,_size,k,1.0,A_beta, k, _crt_in.data(),_ldm,0.,Arns,_size,
                          FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)mn,(int)_size,(int)k,1.0,A_beta,(int)k,_crt_in.data(),(int)_ldm,0.,Arns,(int)_size);
#endif
        }
        Givaro::Timer tred; tred.start();

        reduce(mn,Arns,rda,RNS_MAJOR);
        tred.stop();
        //if(m>1 && n>1) 			std::cerr<<"Reduce : "<<tred.realtime()<<std::endl;

        FFLAS::fflas_delete( A_beta);

#ifdef CHECK_RNS
        bool ok=true;
        for (size_t i=0;i<m;i++)
            for(size_t j=0;j<n;j++)
                for(size_t k=0;k<_size;k++){
                    ok&= (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0)) == (int64_t) Arns[i*n+j+k*rda]);
                    if (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
                        != (int64_t) Arns[i*n+j+k*rda])
                    {
                        std::cout<<((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
                        <<" != "
                        <<(int64_t) Arns[i*n+j+k*rda]
                        <<std::endl;
                    }
                }
        std::cout<<"RNS freduce ... "<<(ok?"OK":"ERROR")<<std::endl;
#endif
        }
    }

    // Arns must be an array of m*n*_size
    // abs(||A||) < 2^(16k)
    inline void rns_double::init_transpose(size_t m, size_t n, double* Arns, size_t rda, const integer* A, size_t lda, size_t k, bool RNS_MAJOR) const
    {
        if (k>_ldm)
            FFPACK::failure()(__func__,__FILE__,__LINE__,"rns_struct: init (too large entry)");

        const size_t mn=m*n;
        if (mn) {
        double *A_beta = FFLAS::fflas_new<double >(mn*k);
        const integer* Aiter=A;
        // split A into A_beta according to a Kronecker transform in base 2^16
        for(size_t j=0;j<n;j++){
            for(size_t i=0;i<m;i++){
                size_t idx=i+j*m;
                const mpz_t*    m0     = reinterpret_cast<const mpz_t*>(Aiter+j+i*lda);
                const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(m0[0]->_mp_d);
                size_t l=0;
                //size_t maxs=std::min(k,(Aiter[j+i*lda].size())<<2);
                size_t maxs=std::min(k,(Aiter[j+i*lda].size())*sizeof(mp_limb_t)/2); // to ensure 32 bits portability
#ifdef __FFLASFFPACK_HAVE_LITTLE_ENDIAN
                if (m0[0]->_mp_size >= 0)
                    for (;l<maxs;l++)
                        A_beta[l+idx*k]=  m0_ptr[l];
                else
                    for (;l<maxs;l++)
                        A_beta[l+idx*k]= - double(m0_ptr[l]);
#else
                if (m0[0]->_mp_size >= 0)
                    for (;l<maxs;l++) {
                        size_t l2 = l ^ ((sizeof(mp_limb_t)/2U) - 1U);
                        A_beta[l+idx*k]=  m0_ptr[l2];
                    }
                else
                    for (;l<maxs;l++) {
                        size_t l2 = l ^ ((sizeof(mp_limb_t)/2U) - 1U);
                        A_beta[l+idx*k]= - double(m0_ptr[l2]);
                    }
#endif
                for (;l<k;l++)
                    A_beta[l+idx*k]=  0.;
            }
        }
        if (RNS_MAJOR==false) {
            // Arns = _crt_in x A_beta^T
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in.data(),_ldm, A_beta, k, 0.,Arns,rda,
                          FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)_size,(int)mn,(int)k,1.0,_crt_in.data(),(int)_ldm,A_beta,(int)k,0.,Arns,(int)rda);
#endif
        }
        else {
            // Arns =  A_beta x _crt_in^T
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,mn,_size,k,1.0,A_beta, k, _crt_in.data(),_ldm,0.,Arns,_size,
                          FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)mn,(int)_size,(int)k,1.0,A_beta,(int)k,_crt_in.data(),(int)_ldm,0.,Arns,(int)_size);
#endif
        }
        reduce(mn,Arns,rda,RNS_MAJOR);

        FFLAS::fflas_delete( A_beta);

#ifdef CHECK_RNS
        bool ok=true;
        for (size_t i=0;i<m;i++)
            for(size_t j=0;j<n;j++)
                for(size_t k=0;k<_size;k++)
                    ok&= (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]<0?(int64_t)_basis[k]:0))
                          == (int64_t) Arns[j*m+i+k*rda]);
        std::cout<<"RNS freduce ... "<<(ok?"OK":"ERROR")<<std::endl;
#endif
        }
    }

    inline void rns_double::convert(size_t m, size_t n, integer gamma, integer* A, size_t lda,
                                    const double* Arns, size_t rda, bool RNS_MAJOR) const
    {
        const size_t  mn= m*n;
        if (mn) {
#ifdef CHECK_RNS
        integer* Acopy=new integer[m*n];
        for(size_t i=0;i<m;i++)
            for(size_t j=0;j<n;j++)
                Acopy[i*n+j]=A[i*lda+j];

#endif

        integer hM= (_M-1)>>1;
        double *A_beta= FFLAS::fflas_new<double>(mn*_ldm);
        Givaro::Timer tfgemmc;tfgemmc.start();
        if (RNS_MAJOR==false) {// compute A_beta = Ap^T x M_beta
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans, mn, _ldm, _size, 1.0 , Arns, rda, _crt_out.data(), _ldm, 0., A_beta,_ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
#else
            cblas_dgemm(CblasRowMajor,CblasTrans, CblasNoTrans, (int)mn, (int)_ldm, (int)_size, 1.0 , Arns, (int)rda, _crt_out.data(), (int)_ldm, 0., A_beta,(int)_ldm);
#endif
        }
        else  {// compute A_beta = Ap x M_Beta
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mn, _ldm, _size, 1.0 , Arns, _size, _crt_out.data(), _ldm, 0., A_beta, _ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, (int)mn, (int)_ldm, (int)_size, 1.0 , Arns, (int)_size, _crt_out.data(), (int)_ldm, 0., A_beta,(int)_ldm);
#endif
        }
        tfgemmc.stop();
        //if(m>1 && n>1) std::cerr<<"fgemm Convert : "<<tfgemmc.realtime()<<std::endl;
        // compute A using inverse Kronecker transform of A_beta expressed in base 2^log_beta
        integer* Aiter= A;
        size_t k=_ldm;
        size_t k4=((k+3)>>2)+ (((k+3)%4==0)?0:1);
        std::vector<uint16_t> A0(k4<<2,0),A1(k4<<2,0),A2(k4<<2,0),A3(k4<<2,0);
        integer a0,a1,a2,a3,res;
        mpz_t *m0,*m1,*m2,*m3;
        m0= reinterpret_cast<mpz_t*>(&a0);
        m1= reinterpret_cast<mpz_t*>(&a1);
        m2= reinterpret_cast<mpz_t*>(&a2);
        m3= reinterpret_cast<mpz_t*>(&a3);
        mp_limb_t *m0_d,*m1_d,*m2_d,*m3_d;
        m0_d = m0[0]->_mp_d;
        m1_d = m1[0]->_mp_d;
        m2_d = m2[0]->_mp_d;
        m3_d = m3[0]->_mp_d;
        m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc = m3[0]->_mp_alloc = (int) (k4*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
        m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size  = m3[0]->_mp_size  = (int) (k4*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
        Givaro::Timer tkroc;
        tkroc.start();
        //		auto sp=SPLITTER();
        //		PARFOR1D(i,m,sp,
        for(size_t i=0;i<m;i++)
            for (size_t j=0;j<n;j++){
                size_t idx=i*n+j;
                for (size_t l=0;l<k;l++){
                    uint64_t tmp=(uint64_t)A_beta[l+idx*k];
                    uint16_t* tptr= reinterpret_cast<uint16_t*>(&tmp);
#ifdef __FFLASFFPACK_HAVE_LITTLE_ENDIAN
                    A0[l  ]= tptr[0];
                    A1[l+1]= tptr[1];
                    A2[l+2]= tptr[2];
                    A3[l+3]= tptr[3];
#else
                    A0[l     ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[3];
                    A1[(l+1) ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[2];
                    A2[(l+2) ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[1];
                    A3[(l+3) ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[0];
#endif
                }
                // see A0,A1,A2,A3 as a the gmp integers a0,a1,a2,a3
                m0[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A0[0]);
                m1[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A1[0]);
                m2[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A2[0]);
                m3[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A3[0]);
                res = a0;res+= a1;res+= a2;res+= a3;
                res%=_M;

                // get the correct result according to the expected sign of A
                if (res>hM)
                    res-=_M;
                if (gamma==0)
                    Aiter[j+i*lda]=res;
                else
                    if (gamma==integer(1))
                        Aiter[j+i*lda]+=res;
                    else
                        if (gamma==integer(-1))
                            Aiter[j+i*lda]=res-Aiter[j+i*lda];
                        else{
                            Aiter[j+i*lda]*=gamma;
                            Aiter[j+i*lda]+=res;
                        }

            }
        tkroc.stop();
        //if(m>1 && n>1) std::cerr<<"Kronecker Convert : "<<tkroc.realtime()<<std::endl;

        m0[0]->_mp_d = m0_d;
        m1[0]->_mp_d = m1_d;
        m2[0]->_mp_d = m2_d;
        m3[0]->_mp_d = m3_d;
        m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc= m3[0]->_mp_alloc = 1;
        m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size = m3[0]->_mp_size  = 0;
        FFLAS::fflas_delete( A_beta);

#ifdef CHECK_RNS
        bool ok=true;
        for (size_t i=0;i<m;i++)
            for(size_t j=0;j<n;j++)
                for(size_t k=0;k<_size;k++){
                    int64_t _p =(int64_t) _basis[k];
                    integer curr=A[i*lda+j] - gamma*Acopy[i*n+j];
                    ok&= ( curr% _p +(curr%_p<0?_p:0) == (int64_t) Arns[i*n+j+k*rda]);
                    //std::cout<<A[i*lda+j]<<" mod "<<(int64_t) _basis[k]<<"="<<(int64_t) Arns[i*n+j+k*rda]<<";"<<std::endl;
                }
        std::cout<<"RNS convert ... "<<(ok?"OK":"ERROR")<<std::endl;
#endif
        }
    }

    inline void rns_double::convert_transpose(size_t m, size_t n, integer gamma, integer* A, size_t lda,
                                              const double* Arns, size_t rda, bool RNS_MAJOR) const
    {
        const size_t  mn= m*n;
        if (mn) {
        integer hM= (_M-1)>>1;
        double *A_beta= FFLAS::fflas_new<double>(mn*_ldm);

        if (RNS_MAJOR==false){
            // compute A_beta = Ap^T x M_beta
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans, mn, _ldm, _size, 1.0 , Arns, rda, _crt_out.data(), _ldm, 0., A_beta, _ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
#else
            cblas_dgemm(CblasRowMajor,CblasTrans, CblasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , Arns,(int) rda, _crt_out.data(),(int) _ldm, 0., A_beta,(int)_ldm);
#endif
        }
        else { // compute A_beta = Ap x M_Beta
#ifndef ENABLE_CHECKER_fgemm
            FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mn, _ldm, _size, 1.0 , Arns, rda, _crt_out.data(), _ldm, 0., A_beta, _ldm, FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
#else
            cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, (int)mn, (int)_ldm, (int)_size, 1.0 , Arns, (int)rda, _crt_out.data(), (int)_ldm, 0., A_beta,(int)_ldm);
#endif
        }
        // compute A using inverse Kronecker transform of A_beta expressed in base 2^log_beta
        integer* Aiter= A;
        size_t k=_ldm;
        size_t k4=((k+3)>>2)+ (((k+3)%4==0)?0:1);
        std::vector<uint16_t> A0(k4<<2,0),A1(k4<<2,0),A2(k4<<2,0),A3(k4<<2,0);
        integer a0,a1,a2,a3,res;
        mpz_t *m0,*m1,*m2,*m3;
        m0= reinterpret_cast<mpz_t*>(&a0);
        m1= reinterpret_cast<mpz_t*>(&a1);
        m2= reinterpret_cast<mpz_t*>(&a2);
        m3= reinterpret_cast<mpz_t*>(&a3);
        mp_limb_t *m0_d,*m1_d,*m2_d,*m3_d;
        m0_d = m0[0]->_mp_d;
        m1_d = m1[0]->_mp_d;
        m2_d = m2[0]->_mp_d;
        m3_d = m3[0]->_mp_d;
        m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc = m3[0]->_mp_alloc = (int32_t)(k4*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
        m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size  = m3[0]->_mp_size  = (int32_t)(k4*8/sizeof(mp_limb_t)); // to ensure 32 bits portability
        for (size_t j=0;j<n;j++)
            for(size_t i=0;i<m;i++){


                size_t idx=i+j*m;
                for (size_t l=0;l<k;l++){
                    uint64_t tmp=(uint64_t)A_beta[l+idx*k];
                    uint16_t* tptr= reinterpret_cast<uint16_t*>(&tmp);
#ifdef __FFLASFFPACK_HAVE_LITTLE_ENDIAN
                    A0[l  ]= tptr[0];
                    A1[l+1]= tptr[1];
                    A2[l+2]= tptr[2];
                    A3[l+3]= tptr[3];
#else
                    A0[l     ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[3];
                    A1[(l+1) ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[2];
                    A2[(l+2) ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[1];
                    A3[(l+3) ^ ((sizeof(mp_limb_t)/2U) - 1U)] = tptr[0];
#endif
                }
                // see A0,A1,A2,A3 as a the gmp integers a0,a1,a2,a3
                m0[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A0[0]);
                m1[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A1[0]);
                m2[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A2[0]);
                m3[0]->_mp_d= reinterpret_cast<mp_limb_t*>(&A3[0]);
                res = a0;res+= a1;res+= a2;res+= a3;
                res%=_M;

                // get the correct result according to the expected sign of A
                if (res>hM)
                    res-=_M;
                if (gamma==0)
                    Aiter[j+i*lda]=res;
                else
                    if (gamma==integer(1))
                        Aiter[j+i*lda]+=res;
                    else
                        if (gamma==integer(-1))
                            Aiter[j+i*lda]=res-Aiter[j+i*lda];
                        else{
                            Aiter[j+i*lda]*=gamma;
                            Aiter[j+i*lda]+=res;
                        }

            }
        m0[0]->_mp_d = m0_d;
        m1[0]->_mp_d = m1_d;
        m2[0]->_mp_d = m2_d;
        m3[0]->_mp_d = m3_d;
        m0[0]->_mp_alloc = m1[0]->_mp_alloc = m2[0]->_mp_alloc= m3[0]->_mp_alloc = 1;
        m0[0]->_mp_size  = m1[0]->_mp_size  = m2[0]->_mp_size = m3[0]->_mp_size  = 0;
        FFLAS::fflas_delete( A_beta);
#ifdef CHECK_RNS
        bool ok=true;
        for (size_t i=0;i<m;i++)
            for(size_t j=0;j<n;j++)
                for(size_t k=0;k<_size;k++){
                    ok&= (((A[i*lda+j] % (int64_t) _basis[k])+(A[i*lda+j]% (int64_t) _basis[k]<0?(int64_t)_basis[k]:0)) == (int64_t) Arns[i+j*m+k*rda]);
                    //std::cout<<A[i*lda+j]<<" mod "<<(int64_t) _basis[k]<<"="<<(int64_t) Arns[i*n+j+k*rda]<<";"<<std::endl;
                }
        std::cout<<"RNS convert ... "<<(ok?"OK":"ERROR")<<std::endl;
#endif // CHECK_RNS
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
            for(size_t i = 0 ; i < n ; i+= _size){
                for(size_t j = 0 ; j < _size ; ++j){
                    //_field_rns.reduce(Arns+i*_size+j);
                    _field_rns[i].reduce(Arns[i*_size+j]);
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


    // TODO: less naive implementation
    inline void rns_double_extended::init(size_t m, double* Arns, const integer* A, size_t lda) const{
        for(size_t i = 0 ; i < m ; ++i){
            for(size_t j = 0 ; j < _size ; ++j){
                Arns[i*_size+j] = (double)((A[i*lda]%integer(_basis[j]))[0]);
            }
        }
    }

    // TODO: less naive implementation
    inline void rns_double_extended::convert(size_t m, integer *A, const double *Arns) const{
        integer hM= (_M-1)/2;
        for(size_t i = 0 ; i < m ; ++i){
            A[i] = 0;
            integer tmp;
            for(size_t j = 0 ; j < _size ; ++j){
                A[i] += ((integer(Arns[i*_size+j])*integer(_MMi[j]))%integer(_basis[j]))*integer(_Mi[j]);
            }
            A[i] %= _M;
            if(A[i] > hM)
                A[i] -= _M;
        }
    }

    // reduce entries of Arns to be less than the rns basis elements
    inline void rns_double_extended::reduce(size_t n, double* Arns, size_t rda, bool RNS_MAJOR) const{

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
        using simd = Simd<double>;
        using vect_t = typename simd::vect_t;

        if(_size % simd::vect_size == 0){
            //#pragma omp parallel for schedule(static, 256)
            for(size_t i = 0 ; i < n ; i++){
                vect_t tmp1, tmp2, tmp3, v, max, basis, inv, neg;
                for(size_t j = 0 ; j < _size ; j+=simd::vect_size){
                    basis = simd::load(_basis.data()+j);
                    inv   = simd::load(_invbasis.data()+j);
                    max   = simd::load(_basisMax.data()+j);
                    neg   = simd::load(_negbasis.data()+j);
                    v     = simd::load(Arns+i*_size+j);
                    tmp2 = modSimd(v, basis, inv, neg);
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
            //#pragma omp parallel for schedule(static, 256)
            for(size_t i = 0 ; i < n ; i++){
                vect_t tmp1, tmp2, tmp3, v, max, basis, inv, neg;
                size_t j = 0;
                for( ; j < ROUND_DOWN(_size, simd::vect_size) ; j+=simd::vect_size){
                    basis = simd::load(_basis.data()+j);
                    inv   = simd::load(_invbasis.data()+j);
                    max   = simd::load(_basisMax.data()+j);
                    neg   = simd::load(_negbasis.data()+j);
                    v     = simd::loadu(Arns+i*_size+j);
                    tmp2 = modSimd(v, basis, inv, neg);
                    tmp1  = simd::greater(tmp2, max);
                    tmp3  = simd::lesser(tmp2, simd::zero());
                    tmp1  = simd::vand(tmp1, neg);
                    tmp3  = simd::vand(tmp3, basis);
                    tmp1  = simd::vor(tmp1, tmp3);
                    tmp2  = simd::add(tmp2, tmp1);
                    simd::storeu(Arns+i*_size+j, tmp2);
                }
                for( ; j < _size ; ++j){
                    _field_rns[j].reduce(Arns[i*_size+j]);
                }
            }
        }
#else

        // TODO : SIMD version
        for(size_t i = 0 ; i < n ; i+= _size){
            for(size_t j = 0 ; j < _size ; ++j){
                //_field_rns.reduce(Arns+i*_size+j);
                _field_rns[i].reduce(Arns[i*_size+j]);
            }
        }

#endif

    }

    } // FFPACK

#endif // __FFLASFFPACK_field_rns_double_INL
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
