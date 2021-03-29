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
/** @file fflas_fgemm/fgemm_classical_mp.inl
 * @brief matrix multiplication with multiprecision input (either over Z or over Z/pZ)
 */


#ifndef __FFPACK_fgemm_classical_INL
#define __FFPACK_fgemm_classical_INL

#include <givaro/modular-integer.h>
#include <givaro/zring.h>
#ifdef PROFILE_FGEMM_MP
#include "fflas-ffpack/utils/timer.h"
#endif
#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/field/rns-integer-mod.h"
#include "fflas-ffpack/field/field-traits.h"
#include "fflas-ffpack/fflas/fflas_helpers.inl"
#include "fflas-ffpack/fflas/fflas_bounds.inl"
namespace FFLAS {

    template<typename Field,
    typename AlgoTrait,
    typename ParSeqTrait>
    struct MMHelper<Field, AlgoTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeqTrait> {
        typedef MMHelper<Field, AlgoTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeqTrait>  Self_t;
        Givaro::Integer normA,normB;
        int recLevel;
        ParSeqTrait parseq;
        MMHelper() : normA(0), normB(0), recLevel(-1) {}
        template <class F2, class A2, class M2, class PS2>
        MMHelper(MMHelper<F2, A2, M2, PS2> H2) :
            normA(H2.normA), normB(H2.normB), recLevel(H2.recLevel), parseq(H2.parseq) {}
        MMHelper(Givaro::Integer Amax, Givaro::Integer Bmax) : normA(Amax), normB(Bmax), recLevel(-1) {}
        MMHelper(const Field& F, size_t m, size_t n, size_t k, ParSeqTrait PS=ParSeqTrait())
        : recLevel(-1), parseq(PS)
        {F.characteristic(normA);F.characteristic(normB); }
        MMHelper(const Field& F, int wino, ParSeqTrait PS=ParSeqTrait()) : recLevel(wino), parseq(PS)
        {F.characteristic(normA);F.characteristic(normB);}
        void setNorm(Givaro::Integer p){normA=normB=p;}

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
        {
            return out <<"Helper: "
            <<typeid(AlgoTrait).name()<<' '
            <<typeid(ModeCategories::ConvertTo<ElementCategories::RNSElementTag>).name()<< ' '
            << M.parseq <<std::endl
            <<"  recLevel = "<<M.recLevel<<std::endl;
        }
    };
    template<typename E,
    typename AlgoTrait,
    typename ParSeqTrait>
    struct MMHelper<FFPACK::RNSInteger<E>, AlgoTrait,ModeCategories::DefaultTag, ParSeqTrait> {
        typedef  MMHelper<FFPACK::RNSInteger<E>, AlgoTrait,ModeCategories::DefaultTag, ParSeqTrait> Self_t;
        Givaro::Integer normA,normB;
        int recLevel;
        ParSeqTrait parseq;
        MMHelper() : normA(0), normB(0), recLevel(-1) {}
        MMHelper(Givaro::Integer Amax, Givaro::Integer Bmax) : normA(Amax), normB(Bmax), recLevel(-1) {}
        MMHelper(const FFPACK::RNSInteger<E>& F, size_t m, size_t n, size_t k, ParSeqTrait PS=ParSeqTrait())
        : recLevel(-1), parseq(PS)
        {F.characteristic(normA);F.characteristic(normB);}
        MMHelper(const FFPACK::RNSInteger<E>& F, int wino, ParSeqTrait PS=ParSeqTrait()) : recLevel(wino), parseq(PS)
        {F.characteristic(normA);F.characteristic(normB);}
        template <class F2, class A2, class M2, class PS2>
        MMHelper(MMHelper<F2, A2, M2, PS2> H2) :
            normA(H2.normA), normB(H2.normB), recLevel(H2.recLevel), parseq(H2.parseq) {}
        void setNorm(Givaro::Integer p){normA=normB=p;}

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
        {
            return out <<"Helper: "
            <<typeid(AlgoTrait).name()<<' '
            <<typeid(ModeCategories::DefaultTag).name()<< ' '
            << M.parseq <<std::endl
            <<"  recLevel = "<<M.recLevel<<std::endl;
        }
    };

    template<typename E,
    typename AlgoTrait,
    typename ParSeqTrait>
    struct MMHelper<FFPACK::RNSIntegerMod<E>, AlgoTrait,ModeCategories::DefaultTag, ParSeqTrait> {
        typedef MMHelper<FFPACK::RNSIntegerMod<E>, AlgoTrait,ModeCategories::DefaultTag, ParSeqTrait> Self_t;
        Givaro::Integer normA,normB;
        int recLevel;
        ParSeqTrait parseq;
        MMHelper() : normA(0), normB(0), recLevel(-1) {}
        MMHelper(Givaro::Integer Amax, Givaro::Integer Bmax) : normA(Amax), normB(Bmax), recLevel(-1) {}
        MMHelper(const FFPACK::RNSIntegerMod<E>& F, size_t m, size_t n, size_t k, ParSeqTrait PS=ParSeqTrait())
        : recLevel(-1), parseq(PS)
        {F.characteristic(normA);F.characteristic(normB);}
        MMHelper(const FFPACK::RNSIntegerMod<E>& F, int wino, ParSeqTrait PS=ParSeqTrait()) : recLevel(wino), parseq(PS)
        {F.characteristic(normA);F.characteristic(normB);}
        // copy constructor from other Field and Algo Traits
        template<class F2, typename AlgoT2, typename FT2, typename PS2>
        MMHelper(MMHelper<F2, AlgoT2, FT2, PS2>& WH) : recLevel(WH.recLevel), parseq(WH.parseq) {}

        void setNorm(Givaro::Integer p){normA=normB=p;}

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
        {
            return out <<"Helper: "
            <<typeid(AlgoTrait).name()<<' '
            <<typeid(ModeCategories::DefaultTag).name()<< ' '
            << M.parseq <<std::endl
            <<"  recLevel = "<<M.recLevel<<std::endl;
        }
    };

    /***********************************
     *** MULTIPRECISION FGEMM OVER Z ***
     ***********************************/

    // fgemm for RnsInteger: simple loop over the moduli
    template<typename RNS, typename ParSeqTrait>
    inline  typename FFPACK::RNSInteger<RNS>::Element_ptr
    fgemm (const FFPACK::RNSInteger<RNS> &F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m, const size_t n,const size_t k,
           const typename FFPACK::RNSInteger<RNS>::Element alpha,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Ad, const size_t lda,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Bd, const size_t ldb,
           const typename FFPACK::RNSInteger<RNS>::Element beta,
           typename FFPACK::RNSInteger<RNS>::Element_ptr Cd, const size_t ldc,
           MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Classic, ModeCategories::DefaultTag, ParSeqHelper::Compose<ParSeqHelper::Sequential, ParSeqTrait> > & H)
    {
#ifdef PROFILE_FGEMM_MP
        Givaro::Timer t;t.start();
#endif
        for(size_t i=0;i<F.size();i++){
            MMHelper<typename RNS::ModField, MMHelperAlgo::Winograd, typename ModeTraits<typename RNS::ModField>::value, ParSeqTrait> Hsub (F.rns()._field_rns[i], H.recLevel, H.parseq.second_component());
            FFLAS::fgemm(F.rns()._field_rns[i],ta,tb,
                         m, n, k, alpha._ptr[i*alpha._stride],
                         Ad._ptr+i*Ad._stride, lda,
                         Bd._ptr+i*Bd._stride, ldb,
                         beta._ptr[i*beta._stride],
                         Cd._ptr+i*Cd._stride, ldc, Hsub);
        }
#ifdef PROFILE_FGEMM_MP
        t.stop();
        std::cerr<<"=========================================="<<std::endl
        <<"Pointwise fgemm : "<<t.realtime()<<" ("<<F.size()<<") moduli "<<std::endl
        <<"=========================================="<<std::endl;
#endif
        return Cd;
    }

    // fgemm for RnsInteger: Sequential => Sequential, Sequential
    template<typename RNS>
    inline  typename FFPACK::RNSInteger<RNS>::Element_ptr
    fgemm (const FFPACK::RNSInteger<RNS> &F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m, const size_t n,const size_t k,
           const typename FFPACK::RNSInteger<RNS>::Element alpha,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Ad, const size_t lda,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Bd, const size_t ldb,
           const typename FFPACK::RNSInteger<RNS>::Element beta,
           typename FFPACK::RNSInteger<RNS>::Element_ptr Cd, const size_t ldc,
           MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Classic, ModeCategories::DefaultTag, ParSeqHelper::Sequential> & H)
    {
        MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Classic, ModeCategories::DefaultTag, ParSeqHelper::Compose<ParSeqHelper::Sequential, ParSeqHelper::Sequential> > H2(H);
        return fgemm (F, ta, tb, m, n, k, alpha, Ad, lda, Bd, ldb, beta, Cd, ldc, H2);
    }

    // fgemm for RnsInteger: handle the moduli in parallel
    template<typename RNS, typename ParSeqTrait>
    inline  typename FFPACK::RNSInteger<RNS>::Element_ptr
    fgemm (const FFPACK::RNSInteger<RNS> &F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m, const size_t n,const size_t k,
           const typename FFPACK::RNSInteger<RNS>::Element alpha,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Ad, const size_t lda,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Bd, const size_t ldb,
           const typename FFPACK::RNSInteger<RNS>::Element beta,
           typename FFPACK::RNSInteger<RNS>::Element_ptr Cd, const size_t ldc,
           MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Classic, ModeCategories::DefaultTag, ParSeqHelper::Compose<ParSeqHelper::Parallel<CuttingStrategy::RNSModulus, StrategyParameter::Threads>, ParSeqTrait> > & H)
    {
#ifdef PROFILE_FGEMM_MP
        Givaro::Timer t;t.start();
#endif
        size_t rns_size = F.size();
        typedef MMHelper<typename RNS::ModField, MMHelperAlgo::Winograd, typename ModeTraits<typename RNS::ModField>::value, ParSeqTrait> SubHelper;
        FORBLOCK1D(iter, rns_size, H.parseq.first_component(),
                   TASK(MODE(CONSTREFERENCE(F,H)),
                        {
                        for(auto i=iter.begin(); i!=iter.end(); ++i)
                        {
                        SubHelper Hsub (F.rns()._field_rns[i], H.recLevel, H.parseq.second_component());
                        FFLAS::fgemm(F.rns()._field_rns[i],ta,tb,
                                     m, n, k, alpha._ptr[i*alpha._stride],
                                     Ad._ptr+i*Ad._stride, lda,
                                     Bd._ptr+i*Bd._stride, ldb,
                                     beta._ptr[i*beta._stride],
                                     Cd._ptr+i*Cd._stride, ldc, Hsub);
                        }
                        })
                  );
#ifdef PROFILE_FGEMM_MP
        t.stop();
        std::cerr<<"=========================================="<<std::endl
        <<"Pointwise fgemm : "<<t.realtime()<<" ("<<F.size()<<") moduli "<<std::endl
        <<"=========================================="<<std::endl;
#endif
        return Cd;
    }

    // fgemm for RnsInteger default parallel version
    template<typename RNS, typename Cut, typename Param>
    inline  typename FFPACK::RNSInteger<RNS>::Element_ptr
    fgemm (const FFPACK::RNSInteger<RNS> &F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m, const size_t n,const size_t k,
           const typename FFPACK::RNSInteger<RNS>::Element alpha,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Ad, const size_t lda,
           typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Bd, const size_t ldb,
           const typename FFPACK::RNSInteger<RNS>::Element beta,
           typename FFPACK::RNSInteger<RNS>::Element_ptr Cd, const size_t ldc,
           MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Classic, ModeCategories::DefaultTag, ParSeqHelper::Parallel<Cut,Param> > & H)
    {
        // compute each fgemm componentwise
        size_t rns_size = F.size();
        size_t nt = H.parseq.numthreads();
        size_t loop_nt = std::min (rns_size, nt);
        size_t iter_nt = nt / loop_nt;
        size_t leftover_nt = nt % loop_nt;
        ParSeqHelper::Parallel<CuttingStrategy::RNSModulus, StrategyParameter::Threads> Hloop (loop_nt);
#ifdef PROFILE_FGEMM_MP
        Givaro::Timer t;t.start();
#endif
        typedef MMHelper<typename RNS::ModField,
                MMHelperAlgo::Winograd,
                typename ModeTraits<typename RNS::ModField>::value,
                ParSeqHelper::Parallel<Cut,Param> > SubPar;

        typedef MMHelper<typename RNS::ModField,
                MMHelperAlgo::Winograd,
                typename ModeTraits<typename RNS::ModField>::value,
                ParSeqHelper::Sequential> SubSeq;

        FORBLOCK1D(iter, rns_size, Hloop,
                   TASK(MODE(CONSTREFERENCE(F,H)),
                        {
                        for(auto i=iter.begin(); i!=iter.end(); ++i)
                        {
                        size_t fgemm_nt = iter_nt;
                        if (i < leftover_nt)
                        fgemm_nt++;
                        if (fgemm_nt>1) // Running a parallel fgemm
                        {
                        SubPar H2(F.rns()._field_rns[i], H.recLevel, ParSeqHelper::Parallel<Cut,Param>(fgemm_nt));
                        fgemm(F.rns()._field_rns[i], ta, tb, m, n, k,
                              alpha._ptr[i*alpha._stride], Ad._ptr+i*Ad._stride,
                              lda, Bd._ptr+i*Bd._stride, ldb,
                              beta._ptr[i*beta._stride], Cd._ptr+i*Cd._stride,
                              ldc, H2);
                        }
                        else // Running a sequential fgemm
                        {
                        SubSeq H2(F.rns()._field_rns[i], H.recLevel, ParSeqHelper::Sequential());
                        fgemm(F.rns()._field_rns[i], ta, tb, m, n, k,
                              alpha._ptr[i*alpha._stride], Ad._ptr+i*Ad._stride,
                              lda, Bd._ptr+i*Bd._stride, ldb,
                              beta._ptr[i*beta._stride], Cd._ptr+i*Cd._stride,
                              ldc, H2);
                        }
                        }
                        }); // TASK
        ); // FLORBLOCK1D

#ifdef PROFILE_FGEMM_MP
        t.stop();

        std::cerr<<"=========================================="<<std::endl
        <<"Pointwise fgemm : "<<t.realtime()<<" ("<<rns_size<<") moduli "<<std::endl
        <<"=========================================="<<std::endl;
#endif
        return Cd;
    }


    template<class ParSeq>
    inline Givaro::Integer*
    fgemm (const Givaro::ZRing<Givaro::Integer>& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m, const size_t n,const size_t k,
           const Givaro::Integer alpha,
           const Givaro::Integer* A, const size_t lda,
           const Givaro::Integer* B, const size_t ldb,
           Givaro::Integer beta,
           Givaro::Integer* C, const size_t ldc,
           MMHelper<Givaro::ZRing<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq >  & H)
    {
        //std::cerr<<"Entering fgemm<ZRing<Integer>> ParSeq"<<std::endl;
#ifdef PROFILE_FGEMM_MP
        Timer chrono;
        chrono.start();
#endif
        if (alpha == 0){
            fscalin(F,m,n,beta,C,ldc);
            return C;
        }

        if (k==0) return C;
        // compute bit size of feasible prime for FFLAS
        size_t _k=k,lk=0;
        while ( _k ) {_k>>=1; ++lk;}
        size_t prime_bitsize= (53-lk)>>1;

        // compute bound on the output
        Givaro::Integer  mA,mB,mC;
        size_t logA,logB;
        mA=H.normA;
        mB=H.normB;
        if (H.normA==0)
            H.normA = InfNorm ((ta==FflasNoTrans)?m:k,(ta==FflasNoTrans)?k:m,A,lda);
        logA = H.normA.bitsize();
        if (H.normB==0)
            H.normB = InfNorm ((tb==FflasNoTrans)?k:n,(tb==FflasNoTrans)?n:k,B,ldb);
        logB = H.normB.bitsize();

        mC = 2*uint64_t(k)*H.normA*H.normB*abs(alpha); // need to use 2x bound to reach both positive and negative

        // A or B is zero, no need to modify C
        if (mC == 0) return C;

        // construct an RNS structure and its associated Domain
        FFPACK::rns_double RNS(mC, prime_bitsize);

        typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
        RnsDomain Zrns(RNS);

        size_t Acold,Arowd,Bcold,Browd;
        if (ta == FFLAS::FflasNoTrans){Arowd=m; Acold = k; }
        else { Arowd=k; Acold = m;}
        if (tb == FFLAS::FflasNoTrans){Browd=k; Bcold = n; }
        else { Browd=n; Bcold = k;}

        // allocate data for RNS representation
        typename RnsDomain::Element_ptr Ap,Bp,Cp;
        Ap = FFLAS::fflas_new(Zrns,Arowd,Acold);
        Bp = FFLAS::fflas_new(Zrns,Browd,Bcold);
        Cp = FFLAS::fflas_new(Zrns,m,n);

#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"-------------------------------"<<std::endl;
        std::cout<<"FGEMM_MP: nb prime: "<<RNS._size<<std::endl;
        std::cout<<"FGEMM_MP:     init: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        chrono.start();
#endif

        // convert the input matrices to RNS representation
        finit_rns(Zrns,Arowd,Acold,(logA/16)+((logA%16)?1:0),A,lda,Ap);
        finit_rns(Zrns,Browd,Bcold,(logB/16)+((logB%16)?1:0),B,ldb,Bp);

#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"FGEMM_MP:   to RNS: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        chrono.start();
#endif

        // perform the fgemm in RNS
        // Classic as no Winograd over ZZ available for the moment
        MMHelper<RnsDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag, ParSeq> H2(Zrns,H.recLevel,H.parseq);

        // compute alpha and beta in RNS
        typename RnsDomain::Element alphap, betap;
        Zrns.init(alphap, alpha);
        Zrns.init(betap, F.zero);

        // call  fgemm
        fgemm(Zrns,ta,tb,m,n,k,alphap,Ap,Acold,Bp,Bcold,betap,Cp,n,H2);

#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"FGEMM_MP:  RNS Mul: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        chrono.start();
#endif


        // convert the RNS output to integer representation (C=beta.C+ RNS^(-1)(Cp) )
        fconvert_rns(Zrns,m,n,beta,C,ldc,Cp);

        FFLAS::fflas_delete(Ap);
        FFLAS::fflas_delete(Bp);
        FFLAS::fflas_delete(Cp);
#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"FGEMM_MP: from RNS: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        std::cout<<"-------------------------------"<<std::endl;
#endif

        return C;
    }



    // Simple switch Winograd -> Classic (waiting for Winograd's algorithm to be generic wrt ModeTrait)
    template<typename RNS, class ModeT>
    inline typename RNS::Element_ptr fgemm (const FFPACK::RNSInteger<RNS> &F,
                                            const FFLAS_TRANSPOSE ta,
                                            const FFLAS_TRANSPOSE tb,
                                            const size_t m, const size_t n,const size_t k,
                                            const typename RNS::Element alpha,
                                            typename RNS::ConstElement_ptr Ad, const size_t lda,
                                            typename RNS::ConstElement_ptr Bd, const size_t ldb,
                                            const typename RNS::Element beta,
                                            typename RNS::Element_ptr Cd, const size_t ldc,
                                            MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Winograd, ModeT, ParSeqHelper::Sequential> & H)
    {
        MMHelper<FFPACK::RNSInteger<RNS>, MMHelperAlgo::Classic, ModeT, ParSeqHelper::Sequential> H2(F, H.recLevel,H.parseq);
        return fgemm(F,ta,tb,m,n,k,alpha,Ad,lda,Bd,ldb,beta,Cd,ldc,H2);
    }

    // template<class ParSeq>
    // inline Givaro::Integer*
    // fgemm (const Givaro::ZRing<Givaro::Integer>& F,
    //        const FFLAS_TRANSPOSE ta,
    //        const FFLAS_TRANSPOSE tb,
    //        const size_t m, const size_t n,const size_t k,
    //        const Givaro::Integer alpha,
    //        const Givaro::Integer* A, const size_t lda,
    //        const Givaro::Integer* B, const size_t ldb,
    //        Givaro::Integer beta,
    //        Givaro::Integer* C, const size_t ldc,
    //        MMHelper<Givaro::ZRing<Givaro::Integer>, MMHelperAlgo::Winograd, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq >  & H)
    // {
    // 	MMHelper<Givaro::ZRing<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq> H2(F, H.recLevel,H.parseq);
    // 	return fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H2);

    // }
    /************************************
     *** MULTIPRECISION FGEMM OVER Fp ***
     ************************************/

    // fgemm for RNSIntegerMod  with Winograd Helper
    template<typename RNS>
    inline typename RNS::Element_ptr fgemm (const FFPACK::RNSIntegerMod<RNS> &F,
                                            const FFLAS_TRANSPOSE ta,
                                            const FFLAS_TRANSPOSE tb,
                                            const size_t m, const size_t n,const size_t k,
                                            const typename RNS::Element alpha,
                                            typename RNS::ConstElement_ptr Ad, const size_t lda,
                                            typename RNS::ConstElement_ptr Bd, const size_t ldb,
                                            const typename RNS::Element beta,
                                            typename RNS::Element_ptr Cd, const size_t ldc,
                                            MMHelper<FFPACK::RNSIntegerMod<RNS>, MMHelperAlgo::Winograd> & H)
    {
        // compute the product over Z
        typedef FFPACK::RNSInteger<RNS> RnsDomain;
        RnsDomain Zrns(F.rns());
        MMHelper<RnsDomain, MMHelperAlgo::Classic> H2(Zrns, H.recLevel,H.parseq);
#ifdef BENCH_PERF_FGEMM_MP
        FFLAS::Timer chrono;chrono.start();
#endif
        fgemm(Zrns,ta,tb,m,n,k,alpha,Ad,lda,Bd,ldb,beta,Cd,ldc,H2);
        // reduce the product mod p (note that entries are larger than p, due to RNS modulo reduction)
        freduce (F, m, n, Cd, ldc);
#ifdef BENCH_PERF_FGEMM_MP
        chrono.stop();
        F.t_igemm+=chrono.realtime();
#endif
        return Cd;
    }


    // fgemm for IntegerDomain with Winograd Helper

    inline Givaro::Integer* fgemm (const Givaro::Modular<Givaro::Integer>& F,
                                   const FFLAS_TRANSPOSE ta,
                                   const FFLAS_TRANSPOSE tb,
                                   const size_t m, const size_t n,const size_t k,
                                   const Givaro::Integer alpha,
                                   const Givaro::Integer *A, const size_t lda,
                                   const Givaro::Integer *B, const size_t ldb,
                                   const Givaro::Integer beta,
                                   Givaro::Integer* C, const size_t ldc,
                                   MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> > & H)
    {
        // compute the product over Z
        //std::cerr<<"Entering fgemm<Modular<Integer>>"<<std::endl;
        typedef Givaro::ZRing<Givaro::Integer> IntegerDomain;
        Givaro::Integer p;
        F.cardinality(p);
        IntegerDomain Z;
        MMHelper<IntegerDomain,MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag> > H2(Z,H.recLevel,H.parseq);
        H2.setNorm(p);

        fgemm(Z,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H2);

        // reduce the product mod p
        freduce (F, m, n, C, ldc);

        return C;
    }
    template<class ParSeq>
    inline Givaro::Integer* fgemm (const Givaro::Modular<Givaro::Integer>& F,
                                   const FFLAS_TRANSPOSE ta,
                                   const FFLAS_TRANSPOSE tb,
                                   const size_t m, const size_t n,const size_t k,
                                   const Givaro::Integer alpha,
                                   const Givaro::Integer *A, const size_t lda,
                                   const Givaro::Integer *B, const size_t ldb,
                                   const Givaro::Integer beta,
                                   Givaro::Integer* C, const size_t ldc,
                                   MMHelper<Givaro::Modular<Givaro::Integer>, MMHelperAlgo::Auto, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq > & H)
    {
        // compute the product over Z
        //std::cerr<<"Entering fgemm<Modular<Integer>> PArSeq"<<std::endl;
        typedef Givaro::ZRing<Givaro::Integer> IntegerDomain;
        Givaro::Integer p;
        F.cardinality(p);
        IntegerDomain Z;
        MMHelper<IntegerDomain,MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq > H2(Z,H.recLevel,H.parseq);
        H2.setNorm(p);

        fgemm(Z,ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,H2);

        // reduce the product mod p
        freduce (F, m, n, C, ldc);

        return C;
    }


    // 	// PARALLEL VERSION (NOT PARALLEL YET)
    // template<class Cut, class Param>
    // inline Givaro::Integer* fgemm (const Givaro::ZRing<Givaro::Integer>& F,
    // 							   const FFLAS_TRANSPOSE ta,
    // 							   const FFLAS_TRANSPOSE tb,
    // 							   const size_t m, const size_t n,const size_t k,
    // 							   const Givaro::Integer alpha,
    // 							   const Givaro::Integer* A, const size_t lda,
    // 							   const Givaro::Integer* B, const size_t ldb,
    // 							   Givaro::Integer beta,
    // 							   Givaro::Integer* C, const size_t ldc,
    // 							   MMHelper<Givaro::ZRing<Givaro::Integer>,MMHelperAlgo::Winograd,FieldCategories::UnparametricTag,ParSeqHelper::Parallel<Cut,Param> > & H){
    // 	MMHelper<Givaro::ZRing<Givaro::Integer>,MMHelperAlgo::Winograd> H2(F,H.recLevel);
    // 	return fgemm(F,ta,tb,m,n,k,alpha,A,lda,B,lda,beta,C,ldc,H2);
    // }





    // fgemm with Modular<ruint<K1,K2> >
    template<size_t K1, size_t K2, class ParSeq>
    inline  RecInt::ruint<K1>*
    fgemm (const Givaro::Modular<RecInt::ruint<K1>,RecInt::ruint<K2> >& F,
           const FFLAS_TRANSPOSE ta,
           const FFLAS_TRANSPOSE tb,
           const size_t m, const size_t n,const size_t k,
           const RecInt::ruint<K1> alpha,
           const RecInt::ruint<K1>* A, const size_t lda,
           const RecInt::ruint<K1>* B, const size_t ldb,
           RecInt::ruint<K1> beta,
           RecInt::ruint<K1>* C, const size_t ldc,
           MMHelper<Givaro::Modular<RecInt::ruint<K1>,RecInt::ruint<K2> >, MMHelperAlgo::Classic, ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq >  & H)
    {
        //std::cerr<<"Entering fgemm<Modular<ruint<K1,K2> > ParSeq with alpha="<<alpha<<" and beta="<<beta<<std::endl;
        //std::cerr << alpha << "*"<< *A <<" * "<< *B <<" + "<<beta<<" * "<<*C<<std::endl;
#ifdef PROFILE_FGEMM_MP
        Timer chrono;
        chrono.start();
#endif
        if (alpha == 0){
            fscalin(F,m,n,beta,C,ldc);
            return C;
        }

        if (k==0) return C;
        // compute bit size of feasible prime for FFLAS
        size_t _k=k,lk=0;
        while ( _k ) {_k>>=1; ++lk;}
        size_t prime_bitsize= (53-lk)>>1;

        // compute bound on the output
        Givaro::Integer  mA,mB,mC;
        F.cardinality(mA);
        mB=mA;
        size_t logA,logB;
        //logA=logB=1<<K1;
	logA=logB=mA.bitsize();

        mC = 2*uint64_t(k)*mA*mB*abs(Givaro::Integer(alpha)); // need to use 2x bound to reach both positive and negative
        // std::cout<<"mA= "<<mA<< "  ("<<mA.bitsize()<<")\n";
        // std::cout<<"mB= "<<mB<< "  ("<<mB.bitsize()<<")\n";
        // std::cout<<"mC= "<<mC<< "  ("<<mC.bitsize()<<")\n";
        // construct an RNS structure and its associated Domain
        FFPACK::rns_double RNS(mC, prime_bitsize);

        typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
        RnsDomain Zrns(RNS);

        size_t Acold,Arowd,Bcold,Browd;
        if (ta == FFLAS::FflasNoTrans){Arowd=m; Acold = k; }
        else { Arowd=k; Acold = m;}
        if (tb == FFLAS::FflasNoTrans){Browd=k; Bcold = n; }
        else { Browd=n; Bcold = k;}

        // allocate data for RNS representation
        typename RnsDomain::Element_ptr Ap,Bp,Cp;
        Ap = FFLAS::fflas_new(Zrns,Arowd,Acold);
        Bp = FFLAS::fflas_new(Zrns,Browd,Bcold);
        Cp = FFLAS::fflas_new(Zrns,m,n);

#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"-------------------------------"<<std::endl;
        std::cout<<"FGEMM_MP: nb prime: "<<RNS._size<<std::endl;
        std::cout<<"FGEMM_MP:     init: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        chrono.start();
#endif

        // convert the input matrices to RNS representation
        RNS.init(Arowd,Acold,Ap._ptr,Ap._stride, A,lda,(logA/16)+((logA%16)?1:0));
        RNS.init(Browd,Bcold,Bp._ptr,Bp._stride, B,ldb,(logB/16)+((logB%16)?1:0));

#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"FGEMM_MP:   to RNS: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        chrono.start();
#endif

        // perform the fgemm in RNS
        // Classic as no Winograd over ZZ available for the moment
        MMHelper<RnsDomain, MMHelperAlgo::Classic, ModeCategories::DefaultTag, ParSeq> H2(Zrns,H.recLevel,H.parseq);

        // compute alpha and beta in RNS
        typename RnsDomain::Element alphap, betap;
        Zrns.init(alphap, alpha);
        Zrns.init(betap, F.zero);

        // call  fgemm
        fgemm(Zrns,ta,tb,m,n,k,alphap,Ap,Acold,Bp,Bcold,betap,Cp,n,H2);

#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"FGEMM_MP:  RNS Mul: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        chrono.start();
#endif


        // convert the RNS output to integer representation (C=beta.C+ RNS^(-1)(Cp) ) modulo mA
        //fconvert_rns(Zrns,m,n,beta,C,ldc,Cp);
        RNS.convert(m,n,beta,C,ldc,Cp._ptr,Cp._stride,mA);


        FFLAS::fflas_delete(Ap);
        FFLAS::fflas_delete(Bp);
        FFLAS::fflas_delete(Cp);
#ifdef PROFILE_FGEMM_MP
        chrono.stop();
        std::cout<<"FGEMM_MP: from RNS: "<<uint64_t(chrono.realtime()*1000)<<"ms"<<std::endl;
        std::cout<<"-------------------------------"<<std::endl;
#endif

        return C;
    }




}// END of namespace FFLAS

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
