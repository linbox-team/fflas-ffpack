/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
		template <class T>
		inline bool unfit(T x){return false;}
		template <>
		inline bool unfit(int64_t x){return (x>limits<int32_t>::max());}
	};

	enum CuttingStrategy {
        SINGLE			,
		ROW_FIXED		,
		COLUMN_FIXED	,
		BLOCK_FIXED		,
		ROW_THREADS		,
		COLUMN_THREADS	,
		BLOCK_THREADS	,
        GRAIN_SIZE		,
		TWO_D			,
		THREE_D_INPLACE	,
		THREE_D_ADAPT	,
		TWO_D_ADAPT		,
		THREE_D
	};

	/*! ParSeqHelper for both fgemm and ftrsm
	*/
	namespace ParSeqHelper {
		struct Parallel{
			Parallel(size_t n=MAX_THREADS, CuttingStrategy m=BLOCK_THREADS):_numthreads(n),_method(m){}

			friend std::ostream& operator<<(std::ostream& out, const Parallel& p) {
				return out << "Parallel: " << p.numthreads() << ',' << p.method();
			}
			size_t numthreads() const { return _numthreads; }
            size_t& set_numthreads(size_t n) { return _numthreads=n; }
			CuttingStrategy method() const { return _method; }            
        private:
			size_t _numthreads;
			CuttingStrategy _method;
            
		};
		struct Sequential{
			Sequential() {}
			//template<class T>
			Sequential(Parallel& ) {}
			friend std::ostream& operator<<(std::ostream& out, const Sequential&) {
				return out << "Sequential";
			}
			size_t numthreads() const { return 1; }
			CuttingStrategy method() const { return SINGLE; }            
		};
	}

	namespace MMHelperAlgo{
		struct Auto{};
		struct Classic{};
		struct Winograd{};
		struct WinogradPar{};
		struct Bini{};
	}


	/*! FGEMM Helper
	*/
	template<class Field,
		 typename AlgoTrait = MMHelperAlgo::Auto,
		 typename ModeTrait = typename ModeTraits<Field>::value,
		 typename ParSeqTrait = ParSeqHelper::Sequential >
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
			DFElt absbeta;
			delayedField.init(absbeta,beta);
			if (beta < 0) absbeta = -beta;
			    //std::cout<<"In MaxDelayedDim: beta = "<<beta<<" absbeta = "<<absbeta<<std::endl;
			DFElt diff = MaxStorableValue - absbeta * std::max(-Cmin, Cmax);
			DFElt AB = std::max (-Amin, Amax) * std::max (-Bmin, Bmax);
			    //std::cout<<"diff = "<<diff<<" AB = "<<AB<<" Helper = "<<std::endl<<*this<<std::endl;
			return static_cast<size_t>(((diff < DFElt(0u))||(AB<DFElt(0u)))? DFElt(0u) : diff / AB);
		}
		bool Aunfit(){ return Protected::unfit(std::max(-Amin,Amax));}
		bool Bunfit(){ return Protected::unfit(std::max(-Bmin,Bmax));}
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
#ifdef DEBUG
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
#ifdef DEBUG
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
#ifdef DEBUG
			for (size_t i=0; i<M;++i)
				for (size_t j=0; j<N;++j)
					if ((A[i*lda+j]>Outmax) || (A[i*lda+j]<Outmin)){
						std::cerr<<"Error in "<<Outmin<<" <= Out["<<i<<", "<<j<<"] = "<<A[i*lda+j]<<" <= "<<Outmax<<std::endl;
						return false;
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
			 DFElt _Cmin, DFElt _Cmax):
			recLevel(w), FieldMin((DFElt)F.minElement()), FieldMax((DFElt)F.maxElement()),
			Amin(_Amin), Amax(_Amax),
			Bmin(_Bmin), Bmax(_Bmax),
			Cmin(_Cmin), Cmax(_Cmax),
			Outmin(0),Outmax(0),
			MaxStorableValue(limits<typename DelayedField::Element>::max()),
			delayedField(F)
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
		TRSMHelper(ParSeqHelper::Parallel _PS):parseq(_PS){}
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
