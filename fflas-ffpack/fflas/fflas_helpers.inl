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

#include "parallel.h"

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

	template <class Field>
	struct associatedDelayedField{
		typedef Field field;
		typedef Field& value; // reference to avoid copying heavy fields 	
	};
	template <>
	struct associatedDelayedField<const FFPACK::Modular<float> >{
		typedef FloatDomain field;
		typedef FloatDomain value;
	};
	template <>
	struct associatedDelayedField<const FFPACK::ModularBalanced<float> >{
		typedef FloatDomain field;	
		typedef FloatDomain value;
	};
	template <>
	struct associatedDelayedField<const FFPACK::Modular<double> >{
		typedef DoubleDomain field;
		typedef DoubleDomain value;
};
	template <>
	struct associatedDelayedField<const FFPACK::ModularBalanced<double> >{
		typedef DoubleDomain field;
		typedef DoubleDomain value;
	};

	// Traits and categories will need to be placed in a proper file later
	namespace FieldCategories {
		//! generic ring.
		struct GenericTag{};
		//! If it can init/convert elements to/from floating point types: float, double
		struct FloatingPointConvertibleTag : public  GenericTag{};
		//! If it is a Modular or ModularBalanced templated by float or double
		struct ModularFloatingPointTag : public GenericTag{};
		//! If it is a Modular or ModularBalanced templated by float or double, and result is not reduced
		struct DelayedModularFloatingPointTag : public GenericTag{};
		//! If it is a multiprecision field
		struct MultiPrecisionTag : public  GenericTag{};
		//! If it is DoubleDomain or a FloatDomain
		struct FloatingPointTag : public GenericTag{};
	}

	/*! FieldTrait
	*/
	template <class Field>
	struct FieldTraits {typedef typename FieldCategories::GenericTag value;};
	template<>
	struct FieldTraits<FFPACK::Modular<double> > {typedef  FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<FFPACK::Modular<float> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<FFPACK::ModularBalanced<double> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<FFPACK::ModularBalanced<float> > {typedef FieldCategories::ModularFloatingPointTag value;};
	template<>
	struct FieldTraits<DoubleDomain> {typedef FieldCategories::FloatingPointTag value;};
	template<>
	struct FieldTraits<FloatDomain> {typedef FieldCategories::FloatingPointTag value;};
	template<>
	struct FieldTraits<FFPACK::Modular<int32_t> > {typedef FieldCategories::FloatingPointConvertibleTag value;};
	template<>
	struct FieldTraits<FFPACK::ModularBalanced<int32_t> > {typedef FieldCategories::FloatingPointConvertibleTag value;};
	//template<> struct FieldTraits<Modular<integer> > {typedef FieldCategories::MultiPrecisionTag value;};


	enum CuttingStrategy {
		ROW_FIXED	,
		COLUMN_FIXED	,
		BLOCK_FIXED	,
		ROW_THREADS	,
		COLUMN_THREADS	,
		BLOCK_THREADS	//,
	};

	/*! ParSeqHelper for both fgemm and ftrsm
	 */
	namespace ParSeqHelper {
		struct Parallel{
			const int numthreads;
			const CuttingStrategy method;
			Parallel(int n=NUM_THREADS, CuttingStrategy m=BLOCK_THREADS):numthreads(n),method(m){}

			friend std::ostream& operator<<(std::ostream& out, const Parallel& p) {
				return out << "Parallel: " << p.numthreads << ',' << p.method;
			}

		};
		struct Sequential{
			Sequential() {}
			    //template<class T> 
			Sequential(Parallel& ) {}
			friend std::ostream& operator<<(std::ostream& out, const Sequential& p) {
				return out << "Sequential";
			}
		};
	}

	namespace MMHelperAlgo{
		struct Auto{};
		struct Classic{};
		struct Winograd{};
		struct Bini{};
	}

	// TODO: fieldMin and fieldMax could be offered by the Field object directly ?
	template <class Field>
	inline double getFieldMin (const Field& F){return 0;}
	template <class Field>
	inline double getFieldMax (const Field& F){return (double) F.characteristic()-1;}
	template<class Element>
	inline double getFieldMin (const FFPACK::ModularBalanced<Element>& F){
		return -((double) F.characteristic()-1)/2;
	}
	template<class Element>
	inline double getFieldMax (const FFPACK::ModularBalanced<Element>& F){
		return ((double) F.characteristic()-1)/2;
	}

       /*! FGEMM Helper
	*/
	template<class Field,
		 typename AlgoTrait = MMHelperAlgo::Auto,
		 typename FieldTrait = typename FieldTraits<Field>::value,
		 typename ParSeqTrait = ParSeqHelper::Sequential >
	struct MMHelper {
		typedef MMHelper<Field,AlgoTrait,FieldTrait,ParSeqTrait> Self_t;
		int recLevel ;
		double FieldMin, FieldMax, Amin, Amax, Bmin, Bmax, Cmin, Cmax, Outmin, Outmax;
		double MaxStorableValue;
		typedef typename associatedDelayedField<const Field>::value DelayedField_t;
		typedef typename associatedDelayedField<const Field>::field DelayedField_v;
		
		const DelayedField_t delayedField;
		ParSeqTrait parseq;
		void initC(){Cmin = FieldMin; Cmax = FieldMax;}
		void initA(){Amin = FieldMin; Amax = FieldMax;}
		void initB(){Bmin = FieldMin; Bmax = FieldMax;}
		void initOut(){Outmin = FieldMin; Outmax = FieldMax;}


		size_t MaxDelayedDim(double beta){ return (size_t) std::max (0.0, floor ( (MaxStorableValue - fabs (beta)*std::max (-Cmin, Cmax) ) / (std::max (-Amin, Amax) * std::max (-Bmin, Bmax))));}

		void setOutBounds(const size_t k, const double alpha, const double beta){
			if (beta<0){
				Outmin = beta*Cmax;
				Outmax = beta*Cmin;
			} else {
				Outmin = beta*Cmin;
				Outmax = beta*Cmax;
			}
			if (alpha >0){
				Outmin += double(k)*alpha*std::min(Amin*Bmax, Amax*Bmin);
				Outmax += double(k)*alpha*std::max(Amin*Bmin, Amax*Bmax);
			}else{
				Outmin += double(k)*alpha*std::max(Amin*Bmin, Amax*Bmax);
				Outmax += double(k)*alpha*std::min(Amin*Bmax, Amax*Bmin);
			}
		}

		bool checkA(const Field& F, const size_t M, const size_t N,
			    typename Field::ConstElement_ptr A, const size_t lda ){
#ifdef DEBUG
			for (size_t i=0; i<M;++i)
				for (size_t j=0; j<N;++j)
					if ((A[i*lda+j]>Amax) || (A[i*lda+j]<Amin)){
						std::cerr<<"Error in "<<i<<j<<std::endl;
						return false;
					}
#endif
			return true;
		}
		bool checkB(const Field& F, const size_t M, const size_t N,
			    typename Field::ConstElement_ptr A, const size_t lda ){
#ifdef DEBUG
			for (size_t i=0; i<M;++i)
				for (size_t j=0; j<N;++j)
					if ((A[i*lda+j]>Bmax) || (A[i*lda+j]<Bmin)){
						std::cerr<<"Error in "<<i<<j<<std::endl;
						return false;
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
						std::cerr<<"Error in A["<<i<<", "<<j<<"] = "<<A[i*lda+j]<<std::endl;
						return false;
					}
#endif
			return true;
		}

		MMHelper(){}
		    //TODO: delayedField constructor has a >0 characteristic even when it is a Double/FloatDomain
		    // correct but semantically not satisfactory
		MMHelper(const Field& F, size_t m, size_t k, size_t n, ParSeqTrait _PS) :
		                recLevel(Protected::WinogradSteps (F, min3(m,k,n))),
				FieldMin(getFieldMin(F)), FieldMax(getFieldMax(F)),
				Amin(FieldMin), Amax(FieldMax),
				Bmin(FieldMin), Bmax(FieldMax),
				Cmin(FieldMin), Cmax(FieldMax),
				Outmin(0.0), Outmax(0.0),
				MaxStorableValue ((double)((1ULL << Protected::Mantissa<typename DelayedField_v::Element>())-1)),
				delayedField((typename Field::Element)F.characteristic()),
				parseq(_PS) {}

		MMHelper(const Field& F, int w, ParSeqTrait _PS=ParSeqTrait()) :
				recLevel(w), //base(FflasDouble),
				FieldMin(getFieldMin(F)), FieldMax(getFieldMax(F)),
				Amin(FieldMin), Amax(FieldMax),
				Bmin(FieldMin), Bmax(FieldMax),
				Cmin(FieldMin), Cmax(FieldMax),
				Outmin(0.0), Outmax(0.0),
				MaxStorableValue ((double)((1ULL << Protected::Mantissa<typename DelayedField_v::Element>())-1)),
				delayedField((typename Field::Element)F.characteristic()),
				parseq(_PS) {}

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
				parseq(WH.parseq) {}

		MMHelper(const Field& F, int w,
			 double _Amin, double _Amax,
			 double _Bmin, double _Bmax,
			 double _Cmin, double _Cmax):
				recLevel(w), FieldMin(getFieldMin(F)), FieldMax(getFieldMax(F)),
				Amin(_Amin), Amax(_Amax),
				Bmin(_Bmin), Bmax(_Bmax),
				Cmin(_Cmin), Cmax(_Cmax),
				MaxStorableValue((double)((1ULL << Protected::Mantissa<typename DelayedField_v::Element>())-1)),
				delayedField((typename Field::Element)F.characteristic()) {}

		friend std::ostream& operator<<(std::ostream& out, const Self_t& M)  {
            return out <<"Helper: "
                <<typeid(AlgoTrait).name()<<' '
                <<typeid(FieldTrait).name()<< ' '
                       << M.parseq <<std::endl
                <<"  recLevel = "<<M.recLevel<<std::endl
                <<"  FieldMin = "<<M.FieldMin<<" FieldMax = "<<M.FieldMax<<std::endl
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
	}
	/*! TRSM Helper
	 */
	template<typename RecIterTrait = StructureHelper::Recursive, typename ParSeqTrait = ParSeqHelper::Sequential>
	struct TRSMHelper {
		ParSeqTrait parseq;
		TRSMHelper(ParSeqTrait _PS):parseq(_PS){}

        template<class Dom, class Algo=FFLAS::MMHelperAlgo::Winograd, class FieldT=typename FFLAS::FieldTraits<Dom>::value>
        FFLAS::MMHelper<Dom, Algo, FieldT, ParSeqTrait> pMMH (Dom& D, size_t m, size_t k, size_t n, ParSeqTrait p) const {
            return FFLAS::MMHelper<Dom, Algo, FieldT, ParSeqTrait>(D,m,k,n,p);
        }
        template<class Dom, class Algo=FFLAS::MMHelperAlgo::Winograd, class FieldT=typename FFLAS::FieldTraits<Dom>::value>
        FFLAS::MMHelper<Dom, Algo, FieldT, ParSeqTrait> pMMH (Dom& D, size_t m, size_t k, size_t n) const {
            return pMMH(D,m,k,n,this->parseq);
        }

	};




} // FFLAS
#endif
