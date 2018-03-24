/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
        struct Auto{
            Auto(){}
            template<class Other>
            Auto(Other& o){}
        };
        struct Classic{
            Classic(){}
            template<class Other>
            Classic(Other & o){}
        };
        struct Winograd{
            int recLevel;
            Winograd():recLevel(-1){}
            Winograd(int r):recLevel(r){}
        };
        struct WinogradPar{
            int recLevel;
            WinogradPar():recLevel(-1){}
            WinogradPar(int r):recLevel(r){}
        };
        struct Bini{};
    }

    template<class ModeT, class ParSeq>
    struct AlgoChooser{typedef MMHelperAlgo::Winograd value;};
    template<class ParSeq>
    struct AlgoChooser<ModeCategories::ConvertTo<ElementCategories::RNSElementTag>, ParSeq>{typedef MMHelperAlgo::Classic value;};

    template<class Field, class ModeTrait,class Enable=void>
    struct Operand{};
    template<class Field, class ModeTrait>
    struct Operand<Field, ModeTrait,
                   typename std::enable_if<std::is_same<ModeTrait,ModeCategories::DelayedTag>::value ||
                                           std::is_same<ModeTrait,ModeCategories::LazyTag>::value ||
                                           std::is_same<ModeTrait,ModeCategories::DefaultBoundedTag>::value ||
                                           std::is_same<ModeTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag> >::value >::type> {
        typedef Operand<Field,ModeTrait,void> Self_t;
        typedef ModeTrait ModeTrait_t;
        typedef typename associatedDelayedField<const Field>::field DelayedField;
        typedef typename DelayedField::Element DFElt;
        template<class OtherOp>
        Operand(const OtherOp& Other) {}
        
        Operand(const DFElt& mi, const DFElt& ma) : min(mi), max(ma){}
        DFElt min, max;
        const DFElt absMax() const{return std::max(static_cast<const DFElt&>(-min),max);}
            //void init(DF){min = F.minElement(); max = F.maxElement();}
        bool unfit(){ return Protected::unfit(absMax());}
        bool check (const Field& F, const FFLAS::FFLAS_TRANSPOSE ta, const size_t M, const size_t N,
                    typename Field::ConstElement_ptr A, const size_t lda ) {
#ifdef __FFLASFFPACK_DEBUG
            for (size_t i=0; i<M;++i)
                for (size_t j=0; j<N;++j){
                    const typename Field::Element x = (ta == FFLAS::FflasNoTrans)? A[i*lda+j] : A[i+j*lda];
                    if (x > max || x < min){
                        std::cerr<<"Error in "<<min<<" < = A["<<i<<", "<<j<<"] ="<<x<<" <= "<<max<<std::endl;
                        return false;
                    }
                }
#endif
            return true;
        }

    };

    template<class Field, class ModeTrait,class Enable=void>
    struct ModeManager_t{
        typedef ModeManager_t<Field,ModeTrait,Enable> Self_t;
        typedef typename associatedDelayedField<const Field>::field_ref DelayedFieldRef;
        typedef typename associatedDelayedField<const Field>::field DelayedField;
        typedef typename DelayedField::Element DFElt;
        typedef typename DelayedField::Element_ptr DFEptr;
        const DelayedFieldRef delayedField;

        Operand<Field,ModeTrait> A;
        Operand<Field,ModeTrait> B;
        Operand<Field,ModeTrait> C;
        Operand<Field,ModeTrait> Out;
  
//        ModeManager_t (){}
        ModeManager_t(const Field&F) : delayedField(F){}
        
        template<class OtherMode>
        ModeManager_t (const Field& F, const OtherMode& OM) : delayedField(F){}

        ModeManager_t (const Field&F, const Operand<Field, ModeTrait>& OA, const Operand<Field,ModeTrait>& OB):
                delayedField(F) {}
    
        ModeManager_t (const Field&F, const Operand<Field, ModeTrait>& OA, const Operand<Field,ModeTrait>& OB, const Operand<Field,ModeTrait>& OC):
                delayedField(F) {}

        template<class OtherMM>
        void updateOutBounds (const OtherMM&M){}

        template<class MM1, class MM2, class MM3, class MM4>
        void setOutBoundsMM(const MM1& M1,const MM2& M2,const MM3& M3,const MM4& M4){}
        
        template<class MMSrc, class MMDest>
        void updateDynPeelHelpers (const MMDest& MMModd, const MMDest& MMNodd, const MMSrc& MMacc){
        }

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M){
                return out <<"ModeManager: "
                           <<typeid(ModeTrait).name()<< ' ';
                }
    };
    
    template<class Field, class ModeTrait>
    struct ModeManager_t <Field, ModeTrait,
                          typename std::enable_if<std::is_same<ModeTrait,ModeCategories::DelayedTag>::value ||
                                                  std::is_same<ModeTrait,ModeCategories::LazyTag>::value ||
                                                  std::is_same<ModeTrait,ModeCategories::DefaultBoundedTag>::value ||
                                                  std::is_same<ModeTrait,ModeCategories::ConvertTo<ElementCategories::RNSElementTag> >::value >::type> {

        typedef ModeManager_t<Field,ModeTrait,void> Self_t;
        typedef typename associatedDelayedField<const Field>::field_ref DelayedFieldRef;
        typedef typename associatedDelayedField<const Field>::field DelayedField;
        typedef typename DelayedField::Element DFElt;
        typedef typename DelayedField::Element_ptr DFEptr;

        const DFElt FieldMin, FieldMax;
        Operand<Field,ModeTrait> A, B, C, Out;
        const DFElt MaxStorableValue;
        const DelayedFieldRef delayedField;
        ModeManager_t (){}
        
        ModeManager_t (const Field& F):
            FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
            A (FieldMin, FieldMax), B (FieldMin, FieldMax), C (FieldMin, FieldMax), Out (0, 0),
            MaxStorableValue ((DFElt)(limits<typename DelayedField::Element>::max())),
            delayedField(F) {}
    
        template<class OtherMode>
        ModeManager_t (const Field& F, const OtherMode& OM):
            FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
            A(OM.A), B(OM.B), C(OM.C), Out(0,0),
            MaxStorableValue ((DFElt)(limits<typename DelayedField::Element>::max())),
            delayedField(F) {}
        
         ModeManager_t (const Field&F, const Operand<Field, ModeTrait>& OA, const Operand<Field,ModeTrait>& OB):
                FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
                A(OA), B(OB), C(0,0), Out(0,0),
                MaxStorableValue ((DFElt)(limits<typename DelayedField::Element>::max())),
                delayedField(F) {}
    
        ModeManager_t (const Field&F, const Operand<Field, ModeTrait>& OA, const Operand<Field,ModeTrait>& OB, const Operand<Field,ModeTrait>& OC):
                FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
                A(OA), B(OB), C(OC), Out(0,0),
                MaxStorableValue ((DFElt)(limits<typename DelayedField::Element>::max())),
                delayedField(F) {}

        size_t MaxDelayedDim(DFElt beta) {
            if (MaxStorableValue < DFElt(0))
                    //Infinte precision delayed field
                return std::numeric_limits<size_t>::max();
            
            DFElt absbeta;
            delayedField.init(absbeta,beta);
            if (beta < 0) absbeta = -beta;
                // This cast is needed when Cmin base type is int8/16_t,
                // getting -Cmin returns a int, not the same base type.
            DFElt diff = MaxStorableValue - absbeta * C.absMax();
            DFElt AB = A.absMax() * B.absMax();
            if ((diff < DFElt(0u))||(AB<DFElt(0u))) return 0;

            DFElt kmax = diff/AB;
            return FFLAS::Protected::min_types<DFElt>(kmax);
        }

        void initC(){C.min = FieldMin; C.max = FieldMax;}
        void initA(){A.min = FieldMin; A.max = FieldMax;}
        void initB(){B.min = FieldMin; B.max = FieldMax;}
        void initOut(){Out.min = FieldMin; Out.max = FieldMax;}

        template<class OtherMM>
        void updateOutBounds (const OtherMM&M){
            Out = M.Out;
            return;
        }
	   template<class MM1, class MM2, class MM3, class MM4>
        void setOutBoundsMM(const MM1& M1,const MM2& M2,const MM3& M3,const MM4& M4){
            Out.min = std::min (M1.ModeManager.Out.min, std::min (M2.ModeManager.Out.min, std::min (M3.ModeManager.Out.min, M4.ModeManager.Out.min)));
            Out.max = std::max (M1.ModeManager.Out.max, std::max (M2.ModeManager.Out.max, std::max (M3.ModeManager.Out.max, M4.ModeManager.Out.max)));
        }
        void setOutBoundsMM(const size_t k, const DFElt alpha, const DFElt beta)
        {
            if (beta<0){
                Out.min = beta*C.max;
                Out.max = beta*C.min;
            } else {
                Out.min = beta*C.min;
                Out.max = beta*C.max;
            }
            if (alpha >0){
                Out.min += DFElt(k)*alpha*std::min(A.min*B.max, A.max*B.min);
                Out.max += DFElt(k)*alpha*std::max(A.min*B.min, A.max*B.max);
            }else{
                Out.min += DFElt(k)*alpha*std::max(A.min*B.min, A.max*B.max);
                Out.max += DFElt(k)*alpha*std::min(A.min*B.max, A.max*B.min);
            }
        }

        void setOutBoundsAdd() { Out.min = A.min+B.min; Out.max = A.max+B.max; }
        void setOutBoundsSub() { Out.min = A.min-B.max; Out.max = A.max-B.min; }
        
        template<class MMSrc, class MMDest>
        void updateDynPeelHelpers (const MMDest& MMModd, const MMDest& MMNodd, const  MMSrc& MMacc){
            Out.min = min4(MMModd.Out.min,MMNodd.Out.min, MMacc.Out.min, Out.min);
            Out.max = max4(MMModd.Out.max,MMNodd.Out.max, MMacc.Out.max, Out.max);
        }
        
        friend std::ostream& operator<<(std::ostream& out, const Self_t& M){
            return out <<"ModeManager: "
                       <<typeid(ModeTrait).name()<< ' '
                       <<" DelayedField = "<<typeid(DelayedField).name()<<std::endl
                       <<"  FieldMin = "<<M.FieldMin<<" FieldMax = "<<M.FieldMax<<std::endl
                       <<"  MaxStorableValue = "<< M.MaxStorableValue <<std::endl
                       <<"  Amin = "<<M.A.min<<" Amax = "<<M.A.max<<std::endl
                       <<"  Bmin = "<<M.B.min<<" Bmax = "<<M.B.max<<std::endl
                       <<"  Cmin = "<<M.C.min<<" Cmax = "<<M.C.max<<std::endl
                       <<"  Outmin = "<<M.Out.min<<" Outmax = "<<M.Out.max<<std::endl;
        }
    };
        //TODO : write generic case with empty functions
    namespace Protected{
        template<class MMDest, class MMSrc>
        void setDynPeelHelpers(MMDest& MMacc, MMDest& MMModd, MMDest& MMNodd, const MMSrc& MMH, const MMSrc& MMHC){
            MMacc.C = MMH.Out;
            MMModd.C = MMHC.C;
            MMModd.A = MMH.B;
            MMModd.B = MMH.A;
            MMNodd.C = MMHC.C;
        }

    } //Protected

	namespace CuttingStrategy{
		struct Single{};
		struct Row{};
		struct Column{};
		struct Block{};
		struct Recursive{};		
	}

	namespace StrategyParameter{
		struct Fixed{};
		struct Threads{};
		struct Grain{};
		struct TwoD{};
		struct TwoDAdaptive{};
		struct ThreeD{};
		struct ThreeDInPlace{};
		struct ThreeDAdaptive{};
	}

	/*! ParSeqHelper for both fgemm and ftrsm
	*/
		/*! ParSeqHelper for both fgemm and ftrsm
	*/
	namespace ParSeqHelper {
		template <typename C=CuttingStrategy::Block, typename P=StrategyParameter::Threads>
		struct Parallel{
			typedef C Cut;
			typedef P Param;
			
			Parallel(size_t n=NUM_THREADS):_numthreads(n){}

			friend std::ostream& operator<<(std::ostream& out, const Parallel& p) {
				return out << "Parallel: " << p.numthreads();
			}
			size_t numthreads() const { return _numthreads; }
			size_t& set_numthreads(size_t n) { return _numthreads=n; }
			// CuttingStrategy method() const { return _method; }
			// StrategyParameter strategy() const { return _param; }
        private:
			size_t _numthreads;
			// CuttingStrategy _method;
			// StrategyParameter _param;
            
		};
		struct Sequential{
			Sequential() {}
			template<class Cut,class Param>
			Sequential(const Parallel<Cut,Param>& ) {}
			friend std::ostream& operator<<(std::ostream& out, const Sequential&) {
				return out << "Sequential";
			}
			size_t numthreads() const { return 1; }
		// 	CuttingStrategy method() const { return SINGLE; }
                // // numthreads==1 ==> a single block
		// 	StrategyParameter strategy() const { return THREADS; }
		};
	}


    template<class Field,
             typename ModeTrait = typename ModeTraits<Field>::value,
             typename ParSeqTrait = ParSeqHelper::Sequential >
     struct AddSubHelper {

        typedef AddSubHelper<Field,ModeTrait,ParSeqTrait> Self_t;
        typedef ModeManager_t<Field,ModeTrait> ModeMgr_t;
        ModeMgr_t ModeManager;
        ParSeqTrait ParSeqManager;

        // Operand<Field,ModeTrait>& A;
        // Operand<Field,ModeTrait>& B;
        // Operand<Field,ModeTrait>& Out;
        
        AddSubHelper(const Field& F, const ParSeqTrait& _PS=ParSeqTrait()):
                ModeManager(F), ParSeqManager(_PS)
                // A(ModeManager.A), B(ModeManager.B), Out(ModeManager.Out)
            {}
         
        
        AddSubHelper(const Field& F, ModeMgr_t _MM,
                  const ParSeqTrait& _PS = ParSeqTrait()):
                ModeManager(_MM), ParSeqManager(_PS)
                // A(ModeManager.A), B(ModeManager.B), Out(ModeManager.Out)
            {}
        
        AddSubHelper(const Field& F,
                     const Operand<Field,ModeTrait>& OA,
                     const Operand<Field,ModeTrait>& OB, 
                     const ParSeqTrait& _PS = ParSeqTrait()):
                ModeManager(F, OA, OB), ParSeqManager(_PS)
                // A(ModeManager.A), B(ModeManager.B), Out(ModeManager.Out)
            {}

        template <class OtherHelper>
        AddSubHelper(const Field& F, const OtherHelper& OH):
                ModeManager(F,OH.ModeManager), ParSeqManager(OH.ParSeqManager)
                // A(ModeManager.A), B(ModeManager.B), Out(ModeManager.Out)
            {}

        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
            {
                return out <<"AddSubHelper: "
                           << M.ModeManager <<std::endl
                           << M.ParSeqManager <<std::endl;
            }
    };

    template<class Field,
             typename AlgoTrait = MMHelperAlgo::Auto,
             typename ModeTrait = typename ModeTraits<Field>::value,
             typename ParSeqTrait = ParSeqHelper::Sequential >
    struct MMHelper {

        typedef MMHelper<Field,AlgoTrait,ModeTrait,ParSeqTrait> Self_t;
        typedef ModeManager_t<Field,ModeTrait> ModeMgr_t;
        AlgoTrait AlgoManager;
        ModeManager_t<Field,ModeTrait> ModeManager;
        ParSeqTrait ParSeqManager;

             //MMHelper(){}
            //TODO: delayedField constructor has a >0 characteristic even when it is a Double/FloatDomain
            // correct but semantically not satisfactory

        MMHelper(const Field& F,
                 const ParSeqTrait& _PS=ParSeqTrait(),
                 const AlgoTrait& _AT = AlgoTrait()) :
                AlgoManager(_AT), ModeManager(F), ParSeqManager(_PS) {}
         
        MMHelper(const Field& F,
                 const AlgoTrait& _AT,
                 const ParSeqTrait& _PS=ParSeqTrait()):
                AlgoManager(_AT), ModeManager(F), ParSeqManager(_PS) {}
        
        MMHelper(const Field& F,
                 ModeManager_t<Field,ModeTrait> _MM,
                 const ParSeqTrait& _PS = ParSeqTrait(),
                 const AlgoTrait& _AT = AlgoTrait()) :
                AlgoManager(_AT), ModeManager(_MM), ParSeqManager(_PS){}
        

        template <class OtherHelper>
        MMHelper(const Field& F, const OtherHelper& OH):
                AlgoManager(OH.AlgoManager), ModeManager(F,OH.ModeManager), ParSeqManager(OH.ParSeqManager){}
        
        // MMHelper(const Field& F, int w, ParSeqTrait _PS=ParSeqTrait()) :
        //         recLevel(w), 
        //         FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
        //         Amin(FieldMin), Amax(FieldMax),
        //         Bmin(FieldMin), Bmax(FieldMax),
        //         Cmin(FieldMin), Cmax(FieldMax),
        //         Outmin(0), Outmax(0),
        //         MaxStorableValue ((DFElt)(limits<typename DelayedField::Element>::max())),
        //         delayedField(F),
        //         parseq(_PS)
        //     {
        //     }

        //     // copy constructor from other Field and Algo Traits
        // template<class F2, typename AT2, typename MT2, typename PS2>
        // MMHelper(MMHelper<F2, AT2, MT2, PS2>& H2) :
        //         AlgoManager(H2.AlgoManager),
        //         ModeManager(H2.ModeManager),
        //         ParSeqManager(H2.ParSeqManager) {}

        template<class MT1, class MT2>
        MMHelper(const Field& F, int w,
                 Operand<Field,MT1>& OA,
                 Operand<Field,MT2>& OB,
                 ParSeqTrait _PS=ParSeqTrait()):
                AlgoManager(w), ModeManager(F,OA,OB),ParSeqManager(_PS) {}

        template<class MT1, class MT2, class MT3>
        MMHelper(const Field& F, int w,
                 Operand<Field,MT1>& OA,
                 Operand<Field,MT2>& OB,
                 Operand<Field,MT3>& OC,
                 ParSeqTrait _PS=ParSeqTrait()):
                AlgoManager(w), ModeManager(F,OA,OB,OC),ParSeqManager(_PS) {}
                
        friend std::ostream& operator<<(std::ostream& out, const Self_t& M)
            {
                return out <<"MMHelper: "
                           << M.AlgoManager <<' '
                           << M.ModeManager <<std::endl
                           << M.ParSeqManager <<std::endl;
            }
    }; // MMHelper


        // to be used in the future, when Winograd's algorithm will be made generic wrt the ModeTrait
        // template <class Field, class AlgoT, class ParSeqH>
        // void copyOutBounds(const MMHelper<Field,AlgoT,ModeCategories::DelayedTag, ParSeqH> &Source,
        //                 MMHelper<Field,AlgoT,ModeCategories::DelayedTag, ParSeqH> & Dest){
        //      Dest.Outmax = Source.Outmax;
        //      Dest.Outmin = Source.Outmin;
        // }
        // template <class Field, class AlgoT, class ParSeqH>
        // void copyOutBounds(const MMHelper<Field,AlgoT,ModeCategories::LazyTag, ParSeqH> &Source,
        //                 MMHelper<Field,AlgoT,ModeCategories::LazyTag, ParSeqH> & Dest){
        //      Dest.Outmax = Source.Outmax;
        //      Dest.Outmin = Source.Outmin;
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
