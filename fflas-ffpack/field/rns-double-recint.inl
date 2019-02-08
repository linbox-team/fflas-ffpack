/*
 * Copyright (C) 2016 the FFLAS-FFPACK group
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


#ifndef __FFLASFFPACK_field_rns_double_recint_INL
#define __FFLASFFPACK_field_rns_double_recint_INL

#include "fflas-ffpack/fflas/fflas_freduce.h"

namespace FFPACK {

    // Arns must be an array of m*n*_size
    // abs(||A||) < 2^(16k)

    template<size_t K>
    inline void rns_double::init(size_t m, size_t n, double* Arns, size_t rda, const RecInt::ruint<K>* A, size_t lda, size_t k, bool RNS_MAJOR) const
    {
        if (k>_ldm){
            FFPACK::failure()(__func__,__FILE__,__LINE__,"rns_double [init] -> rns basis is too small to handle ruint<K> with 2^(16*k) values ");
            std::cerr<<"K="<<K<<" k="<<k<<" _ldm="<<_ldm<<std::endl;
        }
        size_t mn=m*n;
        double *A_beta = FFLAS::fflas_new<double >(mn*k);
        const RecInt::ruint<K>* Aiter=A;
        // split A into A_beta according to a Kronecker transform in base 2^16
        //		auto sp=SPLITTER(MAX_THREADS,FFLAS::CuttingStrategy::Column,FFLAS::StrategyParameter::Threads);

        Givaro::Timer tkr; tkr.start();
#ifndef __FFLASFFPACK_SEQUENTIAL
        //auto sp=SPLITTER(MAX_THREADS);
#endif
        // FOR2D(i,j,m,n,sp,
        //       TASK(MODE(READ(Aiter[0]) READWRITE(A_beta[0])),
        for(size_t i=0;i<m;i++)
            //PAR_BLOCK{
            //			FOR1D(i,m,sp,
            //PARFOR1D(i,m,sp,
            for(size_t j=0;j<n;j++){
                size_t idx=j+i*n;
                const uint16_t* m0_ptr = reinterpret_cast<const uint16_t*>(Aiter+j+i*lda);
                size_t l=0;
                size_t maxs=std::min(k,size_t(1UL<<(K-4)));

                for (;l<maxs;l++){
#ifdef __FFLASFFPACK_HAVE_LITTLE_ENDIAN
                    A_beta[l+idx*k]= m0_ptr[l];
#else
                    A_beta[l+idx*k]= m0_ptr[l^((__RECINT_LIMB_BITS/16U)-1U)]; /* big endian: from __RECINT_LIMB_BITS-bit word to 16-bit words */
#endif
                }
                for (;l<k;l++)
                    A_beta[l+idx*k]=  0.;

                // 	   );
            }


        tkr.stop();
        //if(m>1 && n>1) std::cerr<<"Kronecker : "<<tkr.realtime()<<std::endl;
        if (RNS_MAJOR==false) {
            // Arns = _crt_in x A_beta^T
            Givaro::Timer tfgemm; tfgemm.start();
            PAR_BLOCK{
                FFLAS::fgemm (Givaro::ZRing<double>(), FFLAS::FflasNoTrans,FFLAS::FflasTrans,_size,mn,k,1.0,_crt_in.data(),_ldm,A_beta,k,0.,Arns,rda,
                              //			      FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads>());
                FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive>());

            }
            tfgemm.stop();
            //if(m>1 && n>1) 	std::cerr<<"fgemm : "<<tfgemm.realtime()<<std::endl;
            //			cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)_size,(int)mn,(int)k,1.0,_crt_in.data(),(int)_ldm,A_beta,(int)k,0.,Arns,(int)rda);
            // reduce each row i of Arns modulo moduli[i]
            //for(size_t i=0;i<_size;i++)
            //	FFLAS::freduce (_field_rns[i],mn,Arns+i*rda,1);
        }
        else {
            // Arns =  A_beta x _crt_in^T
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)mn,(int)_size,(int)k,1.0,A_beta,(int)k,_crt_in.data(),(int)_ldm,0.,Arns,(int)_size);
            // reduce each column j of Arns modulo moduli[i]
            //for(size_t i=0;i<_size;i++)
            //	FFLAS::freduce (_field_rns[i],mn,Arns+i,_size);
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


    template<size_t K>
    inline void rns_double::convert(size_t m, size_t n, integer gamma, RecInt::ruint<K>* A, size_t lda,
                                    const double* Arns, size_t rda, integer p,bool RNS_MAJOR) const
    {
        if (p==0 && _M > integer(1)<<(1<<K)){
            std::cerr<<"RNS convert [error] : ruint<"<<K<<"> too small for the rns basis log[2](M)="<<_M.bitsize()<<std::endl;
            std::terminate();
        }

#ifdef CHECK_RNS
        integer* Acopy=new integer[m*n];
        for(size_t i=0;i<m;i++)
            for(size_t j=0;j<n;j++)
                Acopy[i*n+j]=A[i*lda+j];

#endif

        integer hM= (_M-1)>>1;
        size_t  mn= m*n;
        double *A_beta= FFLAS::fflas_new<double>(mn*_ldm);
        Givaro::Timer tfgemmc;tfgemmc.start();
        if (RNS_MAJOR==false)
            // compute A_beta = Ap^T x M_beta
            PAR_BLOCK{
                FFLAS::fgemm(Givaro::ZRing<double>(),FFLAS::FflasTrans, FFLAS::FflasNoTrans,(int) mn,(int) _ldm,(int) _size, 1.0 , Arns,(int) rda, _crt_out.data(),(int) _ldm, 0., A_beta,(int)_ldm,
                             FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Recursive,FFLAS::StrategyParameter::TwoDAdaptive >());
                //				FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Block,FFLAS::StrategyParameter::Threads >());
            }
        else // compute A_beta = Ap x M_Beta
            cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, (int)mn, (int)_ldm, (int)_size, 1.0 , Arns, (int)_size, _crt_out.data(), (int)_ldm, 0., A_beta,(int)_ldm);

        tfgemmc.stop();
        //if(m>1 && n>1) std::cerr<<"fgemm Convert : "<<tfgemmc.realtime()<<std::endl;
        // compute A using inverse Kronecker transform of A_beta expressed in base 2^log_beta

#ifdef CHECK_RNS
        //std::cout<<"CHECKING RNS CONVERT : ruint<"<<K<<"> with log[2](M)="<<_M.bitsize()<<std::endl;
        //std::cout<<"RNS : _ldm*16="<<(_ldm+2)*16<<std::endl;
        bool ok=true;
#endif
        Givaro::Modular<RecInt::ruint<K>,RecInt::ruint<K+1> > Fp(p);
        RecInt::ruint<K>* Aiter= A;
        size_t k=_ldm;
        if ((_ldm+3)*16 > (1<<K) || p!=0){
            //std::cerr<<"ERROR: RNS with recint<"<<K<<"> -> convert needs "<<(_ldm+3)*16<<"bits ...aborting"<<std::endl;
            //std::terminate();
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

#ifdef CHECK_RNS
                    for(size_t k=0;k<_size;k++){
                        int64_t _p =(int64_t) _basis[k];
                        integer curr=res;
                        if ( curr% _p +(curr%_p<0?_p:0) != (int64_t) Arns[i*n+j+k*rda])
                            std::cout<<A[i*lda+j]<<" mod "<<(int64_t) _basis[k]<<"="<<(int64_t) Arns[i*n+j+k*rda]<<";"<<std::endl;
                        ok&= ( curr% _p +(curr%_p<0?_p:0) == (int64_t) Arns[i*n+j+k*rda]);
                    }
#endif

                    if (p!=0){
                        res%=p;
                        if (gamma==0)
                            Aiter[j+i*lda]=RecInt::ruint<K>(res);
                        else
                            if (gamma==integer(1))
                                Fp.addin(Aiter[j+i*lda],RecInt::ruint<K>(res));
                            else
                                if (gamma==p-1)
                                    Fp.sub(Aiter[j+i*lda],RecInt::ruint<K>(res),Aiter[j+i*lda]);
                                else{
                                    Fp.mulin(Aiter[j+i*lda],RecInt::ruint<K>(gamma));
                                    Fp.addin(Aiter[j+i*lda],RecInt::ruint<K>(res));
                                }
                    }else {

                        // get the correct result according to the expected sign of A
                        if (res>hM)
                            res-=_M;

                        if (gamma==0)
                            Aiter[j+i*lda]=RecInt::ruint<K>(res);
                        else
                            if (gamma==integer(1))
                                Aiter[j+i*lda]+=RecInt::ruint<K>(res);
                            else
                                if (gamma==integer(-1))
                                    Aiter[j+i*lda]=RecInt::ruint<K>(res)-Aiter[j+i*lda];
                                else{
                                    Aiter[j+i*lda]*=RecInt::ruint<K>(gamma);
                                    Aiter[j+i*lda]+=RecInt::ruint<K>(res);
                                }
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


        }
        else {
            //size_t k4=((k+3)>>2)+ (((k+3)%4==0)?0:1);

            std::vector<uint16_t> A0(1<<(K-4),0),A1(1<<(K-4),0),A2(1<<(K-4),0),A3(1<<(K-4),0);
            RecInt::ruint<K> *a0,*a1,*a2,*a3,res;
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
                    a0= reinterpret_cast<RecInt::ruint<K>*>(&A0[0]);
                    a1= reinterpret_cast<RecInt::ruint<K>*>(&A1[0]);
                    a2= reinterpret_cast<RecInt::ruint<K>*>(&A2[0]);
                    a3= reinterpret_cast<RecInt::ruint<K>*>(&A3[0]);

                    res = *a0;res+= *a1;res+= *a2;res+= *a3;
                    res%= RecInt::ruint<K>(_M);

#ifdef CHECK_RNS
                    for(size_t k=0;k<_size;k++){
                        int64_t _p =(int64_t) _basis[k];
                        integer curr=res;
                        if ( curr% _p +(curr%_p<0?_p:0) != (int64_t) Arns[i*n+j+k*rda])
                            std::cout<<A[i*lda+j]<<" mod "<<(int64_t) _basis[k]<<"="<<(int64_t) Arns[i*n+j+k*rda]<<";"<<std::endl;
                        ok&= ( curr% _p +(curr%_p<0?_p:0) == (int64_t) Arns[i*n+j+k*rda]);
                    }
#endif


                    // get the correct result according to the expected sign of A
                    //if (res>hM)
                    //	res-=_M;
                    if (gamma==0)
                        Aiter[j+i*lda]=res;
                    else
                        if (gamma==1)
                            Aiter[j+i*lda]+=res;
                        else
                            if (gamma==-1)
                                Aiter[j+i*lda]=res-Aiter[j+i*lda];
                            else{
                                Aiter[j+i*lda]*=RecInt::ruint<K>(gamma);
                                Aiter[j+i*lda]+=res;
                            }

                }
            tkroc.stop();
        }
        //if(m>1 && n>1) std::cerr<<"Kronecker Convert : "<<tkroc.realtime()<<std::endl;

        FFLAS::fflas_delete( A_beta);

#ifdef CHECK_RNS
        std::cout<<"CHECKING RNS CONVERT : ruint<"<<K<<"> with log[2](M)="<<_M.bitsize()<<std::endl;
        // std::cout<<"RNS : _ldm*16="<<(_ldm+2)*16<<std::endl;
        // bool ok=true;
        // for (size_t i=0;i<m;i++)
        // 	for(size_t j=0;j<n;j++)
        // 		for(size_t k=0;k<_size;k++){
        // 			int64_t _p =(int64_t) _basis[k];
        // 			integer curr=integer(A[i*lda+j]) - gamma*Acopy[i*n+j]; if (p!=0) curr%=p;
        // 			if ( curr% _p +(curr%_p<0?_p:0) != (int64_t) Arns[i*n+j+k*rda])
        // 				std::cout<<A[i*lda+j]<<" mod "<<(int64_t) _basis[k]<<"="<<(int64_t) Arns[i*n+j+k*rda]<<";"<<std::endl;
        // 			ok&= ( curr% _p +(curr%_p<0?_p:0) == (int64_t) Arns[i*n+j+k*rda]);

        // 		}

        std::cout<<"RNS convert ... "<<(ok?"OK":"ERROR")<<std::endl;
#endif
    }

    } // FFPACK

#endif // __FFLASFFPACK_field_rns_double_recint_INL
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
