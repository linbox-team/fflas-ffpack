/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

/** @file fflas/fflas_mmhelper.h
 * @brief  Matrix-Matrix Helper class
 */

#ifndef __FFLASFFPACK_fflas_fflas_mmhelper_INL
#define __FFLASFFPACK_fflas_fflas_mmhelper_INL

#include "fflas-ffpack/field/field-traits.h"
#include "fflas-ffpack/paladin/parallel.h"
#include "fflas-ffpack/utils/flimits.h"

#include <algorithm> // std::max

namespace FFLAS{ namespace Protected{
    /** \brief Computes the number of recursive levels to perform.
     *
     * \param m the common dimension in the product AxB
     */
    template<class Field>
    int WinogradSteps (const Field & F, const size_t & m);

}//Protected
}//FFLAS

namespace FFLAS {

    namespace Protected{
        template <class DFE> inline size_t min_types(const DFE& k) {return static_cast<size_t>(k);}
#if __FFLASFFPACK_SIZEOF_LONG == 4
        template <> inline size_t min_types(const double& k) {return static_cast<size_t>(std::min(k,double(std::numeric_limits<size_t>::max())));}
        template <> inline size_t min_types(const int64_t& k) {return static_cast<size_t>(std::min(k,int64_t(std::numeric_limits<size_t>::max())));}
#endif
        template <> inline size_t min_types(const RecInt::rint<6>& k) {return static_cast<size_t>(uint64_t(std::min(k,RecInt::rint<6>(uint64_t(std::numeric_limits<size_t>::max())))));}
        template <> inline size_t min_types(const RecInt::rint<7>& k) {return static_cast<size_t>(uint64_t(std::min(k,RecInt::rint<7>(uint64_t(std::numeric_limits<size_t>::max())))));}
        template <> inline size_t min_types(const RecInt::rint<8>& k) {return static_cast<size_t>(uint64_t(std::min(k,RecInt::rint<8>(uint64_t(std::numeric_limits<size_t>::max())))));}
        template <> inline size_t min_types(const RecInt::rint<9>& k) {return static_cast<size_t>(uint64_t(std::min(k,RecInt::rint<9>(uint64_t(std::numeric_limits<size_t>::max())))));}
        template <> inline size_t min_types(const RecInt::rint<10>& k) {return static_cast<size_t>(uint64_t(std::min(k,RecInt::rint<10>(uint64_t(std::numeric_limits<size_t>::max())))));}
        template <> inline size_t min_types(const Givaro::Integer& k) {return static_cast<size_t>(uint64_t(std::min(k,Givaro::Integer(uint64_t(std::numeric_limits<size_t>::max())))));}

        template <class T>
        inline bool unfit(T x){return false;}
        template <>
        inline bool unfit(int64_t x){return (x>limits<int32_t>::max());}
        template <size_t K>
        inline bool unfit(RecInt::rint<K> x){return (x > RecInt::rint<K>(limits<RecInt::rint<K-1>>::max()));}
        template <>
        inline bool unfit(RecInt::rint<6> x){return (x > limits<int32_t>::max());}
    }

    namespace MMHelperAlgo{
        struct Auto{};
        struct Classic{};
        struct DivideAndConquer{};
        struct Winograd{};
        struct WinogradPar{};
        struct Bini{};
    }

    template<class ModeT, class ParSeq>
    struct AlgoChooser{typedef MMHelperAlgo::Winograd value;};
    template<class ParSeq>
    struct AlgoChooser<ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq>{typedef MMHelperAlgo::Classic value;};

    template<class Field,
    typename AlgoTrait = MMHelperAlgo::Auto,
    typename ModeTrait = typename ModeTraits<Field>::value,
    typename ParSeqTrait = ParSeqHelper::Sequential >
    struct MMHelper;
    /*! FGEMM  Helper for Default and ConvertTo modes of operation
    */
    template<class Field,
    typename AlgoTrait,
    typename ParSeqTrait >
    struct MMHelper<Field, AlgoTrait, ModeCategories::DefaultTag, ParSeqTrait>
    {
        typedef MMHelper<Field,AlgoTrait, ModeCategories::DefaultTag,ParSeqTrait> Self_t;
        int recLevel ;
        ParSeqTrait parseq;

        MMHelper(){}
        MMHelper(const Field& F, size_t m, size_t k, size_t n, ParSeqTrait _PS) : recLevel(-1), parseq(_PS) {}
        MMHelper(const Field& F, int w, ParSeqTrait _PS=ParSeqTrait()) : recLevel(w), parseq(_PS) {}

        // copy constructor from other Field and Algo Traits
        template<class F2, typename AlgoT2, typename FT2, typename PS2>
        MMHelper(MMHelper<F2, AlgoT2, FT2, PS2>& WH) : recLevel(WH.recLevel), parseq(WH.parseq) {}

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
        {
            return out <<"Helper: "
            <<typeid(AlgoTrait).name()<<' '
            <<typeid(ModeCategories::DefaultTag).name()<< ' '
            << M.parseq <<std::endl
            <<"  recLevel = "<<M.recLevel<<std::endl;
        }
    };
    template<class Field,
    typename AlgoTrait,
    typename Dest,
    typename ParSeqTrait>
    struct MMHelper<Field, AlgoTrait, ModeCategories::ConvertTo<Dest>, ParSeqTrait>
    {
        typedef MMHelper<Field,AlgoTrait, ModeCategories::ConvertTo<Dest>,ParSeqTrait> Self_t;
        int recLevel ;
        ParSeqTrait parseq;

        MMHelper(){}
        MMHelper(const Field& F, size_t m, size_t k, size_t n, ParSeqTrait _PS) : recLevel(-1), parseq(_PS) {}
        MMHelper(const Field& F, int w, ParSeqTrait _PS=ParSeqTrait()) : recLevel(w), parseq(_PS) {}

        // copy constructor from other Field and Algo Traits
        template<class F2, typename AlgoT2, typename FT2, typename PS2>
        MMHelper(MMHelper<F2, AlgoT2, FT2, PS2>& WH) : recLevel(WH.recLevel), parseq(WH.parseq) {}

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
        {
            return out <<"Helper: "
            <<typeid(AlgoTrait).name()<<' '
            <<typeid(ModeCategories::ConvertTo<Dest>).name()<< ' '
            << M.parseq <<std::endl
            <<"  recLevel = "<<M.recLevel<<std::endl;
        }
    };
    // MMHelper for Delayed and Lazy Modes of operation
    template<class Field,
    typename AlgoTrait,
    typename ModeTrait,
    typename ParSeqTrait>
    struct MMHelper {
        typedef MMHelper<Field,AlgoTrait,ModeTrait,ParSeqTrait> Self_t;
        typedef typename associatedDelayedField<const Field>::type DelayedField_t;
        typedef typename associatedDelayedField<const Field>::field DelayedField;
        typedef typename DelayedField::Element DFElt;
        int recLevel ;
        DFElt FieldMin, FieldMax, Amin, Amax, Bmin, Bmax, Cmin, Cmax, Outmin, Outmax;
        DFElt MaxStorableValue;

        const DelayedField_t delayedField;
        ParSeqTrait parseq;
        void initC(){Cmin = FieldMin; Cmax = FieldMax;}
        void initA(){Amin = FieldMin; Amax = FieldMax;}
        void initB(){Bmin = FieldMin; Bmax = FieldMax;}
        void initOut(){Outmin = FieldMin; Outmax = FieldMax;}


        size_t MaxDelayedDim(DFElt beta)
        {
            if (MaxStorableValue < DFElt(0))
                //Infinte precision delayed field
                return std::numeric_limits<size_t>::max();

            DFElt absbeta;
            delayedField.init(absbeta);
            delayedField.assign(absbeta, beta);
            if (beta < 0) absbeta = -beta;
            // This cast is needed when Cmin base type is int8/16_t,
            // getting -Cmin returns a int, not the same base type.
            DFElt diff = MaxStorableValue - absbeta * std::max(static_cast<const DFElt&>(-Cmin), Cmax);
            DFElt AB = std::max(static_cast<const DFElt&>(-Amin), Amax) * std::max(static_cast<const DFElt&>(-Bmin), Bmax);

            if ((diff < DFElt(0u))||(AB<DFElt(0u))) return 0;

            DFElt kmax = diff/AB;
            return FFLAS::Protected::min_types<DFElt>(kmax);
            // if (kmax > std::numeric_limits<size_t>::max())
            // 	return std::numeric_limits<size_t>::max();
            // else
            // 	return kmax;
        }
        bool Aunfit(){ return Protected::unfit(std::max(static_cast<const DFElt&>(-Amin),Amax));}
        bool Bunfit(){ return Protected::unfit(std::max(static_cast<const DFElt&>(-Bmin),Bmax));}
        void setOutBounds(const size_t k, const DFElt alpha, const DFElt beta)
        {
            if (beta<0){
                Outmin = beta*Cmax;
                Outmax = beta*Cmin;
            } else {
                Outmin = beta*Cmin;
                Outmax = beta*Cmax;
            }
            if (alpha >0){
                Outmin += DFElt(k)*alpha*std::min(Amin*Bmax, Amax*Bmin);
                Outmax += DFElt(k)*alpha*std::max(Amin*Bmin, Amax*Bmax);
            }else{
                Outmin += DFElt(k)*alpha*std::max(Amin*Bmin, Amax*Bmax);
                Outmax += DFElt(k)*alpha*std::min(Amin*Bmax, Amax*Bmin);
            }
        }

        bool checkA(const Field& F, const FFLAS::FFLAS_TRANSPOSE ta, const size_t M, const size_t N,
                    typename Field::ConstElement_ptr A, const size_t lda )
        {
#ifdef __FFLASFFPACK_DEBUG
            for (size_t i=0; i<M;++i)
                for (size_t j=0; j<N;++j){
                    const typename Field::Element x = (ta == FFLAS::FflasNoTrans)? A[i*lda+j] : A[i+j*lda];
                    if (x > Amax || x < Amin){
                        std::cerr<<"Error in "<<Amin<<" < = A["<<i<<", "<<j<<"] ="<<x<<" <= "<<Amax<<std::endl;
                        return false;
                    }
                }
#endif
            return true;
        }

        bool checkB(const Field& F, const FFLAS::FFLAS_TRANSPOSE tb, const size_t M, const size_t N,
                    typename Field::ConstElement_ptr B, const size_t ldb)
        {
#ifdef __FFLASFFPACK_DEBUG
            for (size_t i=0; i<M;++i)
                for (size_t j=0; j<N;++j){
                    const typename Field::Element x = (tb == FFLAS::FflasNoTrans)? B[i*ldb+j] : B[i+j*ldb];
                    if (x > Bmax || x < Bmin){
                        std::cerr<<"Error in "<<Bmin<<" < = B["<<i<<", "<<j<<"] ="<<B[i*ldb+j]<<" <= "<<Bmax<<std::endl;
                        return false;
                    }
                }
#endif
            return true;
        }

        bool checkOut(const Field& F, const size_t M, const size_t N,
                      typename Field::ConstElement_ptr A, const size_t lda ){
#ifdef __FFLASFFPACK_DEBUG
            for (size_t i=0; i<M;++i)
                for (size_t j=0; j<N;++j)
                    if ((A[i*lda+j]>Outmax) || (A[i*lda+j]<Outmin)){
                        std::cerr<<"Error in "<<Outmin<<" <= Out["<<i<<", "<<j<<"] = "<<A[i*lda+j]<<" <= "<<Outmax<<std::endl;
                        return false;
                    }
#endif
            return true;
        }
        bool checkOut(const Field& F, FFLAS_UPLO uplo, const size_t M, const size_t N,
                      typename Field::ConstElement_ptr A, const size_t lda ){
#ifdef __FFLASFFPACK_DEBUG
            if (uplo == FflasUpper){
                for (size_t i=0; i<M;++i)
                    for (size_t j=i; j<N;++j)
                        if ((A[i*lda+j]>Outmax) || (A[i*lda+j]<Outmin)){
                            std::cerr<<"Error in "<<Outmin<<" <= Out["<<i<<", "<<j<<"] = "<<A[i*lda+j]<<" <= "<<Outmax<<std::endl;
                            return false;
                        }
            } else {
                for (size_t i=0; i<M;++i)
                    for (size_t j=0; j<=i;++j)
                        if ((A[i*lda+j]>Outmax) || (A[i*lda+j]<Outmin)){
                            std::cerr<<"Error in "<<Outmin<<" <= Out["<<i<<", "<<j<<"] = "<<A[i*lda+j]<<" <= "<<Outmax<<std::endl;
                            return false;
                        }

            }
#endif
            return true;
        }

        MMHelper(){}
        //TODO: delayedField constructor has a >0 characteristic even when it is a Double/FloatDomain
        // correct but semantically not satisfactory
        MMHelper(const Field& F, size_t m, size_t k, size_t n, ParSeqTrait _PS) :
            recLevel(-1),
            FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
            Amin(FieldMin), Amax(FieldMax),
            Bmin(FieldMin), Bmax(FieldMax),
            Cmin(FieldMin), Cmax(FieldMax),
            Outmin(0), Outmax(0),
            MaxStorableValue ((DFElt)(limits<typename DelayedField::Element>::max())),
            delayedField(F),
            // delayedField((typename Field::Element)F.characteristic()),
            parseq(_PS)
        {
        }

        MMHelper(const Field& F, int w, ParSeqTrait _PS=ParSeqTrait()) :
            recLevel(w),
            FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
            Amin(FieldMin), Amax(FieldMax),
            Bmin(FieldMin), Bmax(FieldMax),
            Cmin(FieldMin), Cmax(FieldMax),
            Outmin(0), Outmax(0),
            MaxStorableValue ((DFElt)(limits<typename DelayedField::Element>::max())),
            delayedField(F),
            parseq(_PS)
        {
        }

        // copy constructor from other Field and Algo Traits
        template<class F2, typename AlgoT2, typename FT2, typename PS2>
        MMHelper(MMHelper<F2, AlgoT2, FT2, PS2>& WH) :
            recLevel(WH.recLevel),
            FieldMin(WH.FieldMin), FieldMax(WH.FieldMax),
            Amin(WH.Amin), Amax(WH.Amax),
            Bmin(WH.Bmin), Bmax(WH.Bmax),
            Cmin(WH.Cmin), Cmax(WH.Cmax),
            Outmin(WH.Outmin), Outmax(WH.Outmax),
            MaxStorableValue(WH.MaxStorableValue),
            delayedField(WH.delayedField),
            parseq(WH.parseq)
        {
        }

        MMHelper(const Field& F, int w,
                 DFElt _Amin, DFElt _Amax,
                 DFElt _Bmin, DFElt _Bmax,
                 DFElt _Cmin, DFElt _Cmax,
                 ParSeqTrait _PS=ParSeqTrait()):
            recLevel(w), FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
            Amin(_Amin), Amax(_Amax),
            Bmin(_Bmin), Bmax(_Bmax),
            Cmin(_Cmin), Cmax(_Cmax),
            Outmin(0),Outmax(0),
            MaxStorableValue(limits<typename DelayedField::Element>::max()),
            delayedField(F),
            parseq(_PS)
        {
        }

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
        {
            return out <<"Helper: "
            <<typeid(AlgoTrait).name()<<' '
            <<typeid(ModeTrait).name()<< ' '
            << M.parseq <<std::endl
            <<" DelayedField = "<<typeid(DelayedField).name()<<std::endl
            <<"  recLevel = "<<M.recLevel<<std::endl
            <<"  FieldMin = "<<M.FieldMin<<" FieldMax = "<<M.FieldMax<<std::endl
            <<"  MaxStorableValue = "<< M.MaxStorableValue <<std::endl
            <<"  Amin = "<<M.Amin<<" Amax = "<<M.Amax<<std::endl
            <<"  Bmin = "<<M.Bmin<<" Bmax = "<<M.Bmax<<std::endl
            <<"  Cmin = "<<M.Cmin<<" Cmax = "<<M.Cmax<<std::endl
            <<"  Outmin = "<<M.Outmin<<" Outmax = "<<M.Outmax<<std::endl;
        }
    }; // MMHelper


    // to be used in the future, when Winograd's algorithm will be made generic wrt the ModeTrait
    // template <class Field, class AlgoT, class ParSeqH>
    // void copyOutBounds(const MMHelper<Field,AlgoT,ModeCategories::DelayedTag, ParSeqH> &Source,
    // 		   MMHelper<Field,AlgoT,ModeCategories::DelayedTag, ParSeqH> & Dest){
    // 	Dest.Outmax = Source.Outmax;
    // 	Dest.Outmin = Source.Outmin;
    // }
    // template <class Field, class AlgoT, class ParSeqH>
    // void copyOutBounds(const MMHelper<Field,AlgoT,ModeCategories::LazyTag, ParSeqH> &Source,
    // 		   MMHelper<Field,AlgoT,ModeCategories::LazyTag, ParSeqH> & Dest){
    // 	Dest.Outmax = Source.Outmax;
    // 	Dest.Outmin = Source.Outmin;
    // }
    // template <class MMH1, class MMH2>
    // void copyOutBounds(const MMH1 &Source, MMH2 & Dest){}
    /*! StructureHelper for ftrsm
    */
    namespace StructureHelper {
        struct Recursive{};
        struct Iterative{};
        struct Hybrid{};
    }

    /*! TRSM Helper
    */
    template<typename RecIterTrait = StructureHelper::Recursive, typename ParSeqTrait = ParSeqHelper::Sequential>
    struct TRSMHelper {
        ParSeqTrait parseq;
        template<class Cut,class Param>
        TRSMHelper(ParSeqHelper::Parallel<Cut,Param> _PS):parseq(_PS){}
        TRSMHelper(ParSeqHelper::Sequential _PS):parseq(_PS){}
        template<typename RIT, typename PST>
        TRSMHelper(TRSMHelper<RIT,PST>& _TH):parseq(_TH.parseq){}

        template<class Dom, class Algo=FFLAS::MMHelperAlgo::Winograd, class ModeT=typename FFLAS::ModeTraits<Dom>::value>
        FFLAS::MMHelper<Dom, Algo, ModeT, ParSeqTrait> pMMH (Dom& D, size_t m, size_t k, size_t n, ParSeqTrait p) const {
            return FFLAS::MMHelper<Dom, Algo, ModeT, ParSeqTrait>(D,m,k,n,p);
        }

        template<class Dom, class Algo=FFLAS::MMHelperAlgo::Winograd, class ModeT=typename FFLAS::ModeTraits<Dom>::value>
        FFLAS::MMHelper<Dom, Algo, ModeT, ParSeqTrait> pMMH (Dom& D, size_t m, size_t k, size_t n) const {
            return pMMH(D,m,k,n,this->parseq);
        }

    };




} // FFLAS
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
