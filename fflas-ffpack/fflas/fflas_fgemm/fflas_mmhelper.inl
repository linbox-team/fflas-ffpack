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
	struct associatedDelayedField{typedef Field value;};
	template <>
	struct associatedDelayedField<FFPACK::Modular<float> >{typedef FloatDomain value;};
	template <>
	struct associatedDelayedField<FFPACK::ModularBalanced<float> >{typedef FloatDomain value;};
	template <>
	struct associatedDelayedField<FFPACK::Modular<double> >{typedef DoubleDomain value;};
	template <>
	struct associatedDelayedField<FFPACK::ModularBalanced<double> >{typedef DoubleDomain value;};
        
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
	template<typename  Element>
	struct FieldTraits<FFPACK::Modular<Element> > {typedef FieldCategories::FloatingPointConvertibleTag value;};
	//template<> struct FieldTraits<Modular<integer> > {typedef FieldCategories::MultiPrecisionTag value;};

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

	template<class Field, typename AlgoTrait = MMHelperAlgo::Auto, typename FieldTrait = typename FieldTraits<Field>::value >
	struct MMHelper {
		int recLevel ;
		double FieldMin, FieldMax, Amin, Amax, Bmin, Bmax, Cmin, Cmax, Outmin, Outmax;
		double MaxStorableValue;
		typedef  typename associatedDelayedField<Field>::value DelayedField_t;
		DelayedField_t delayedField;

		void initC(){Cmin = FieldMin; Cmax = FieldMax;}
		void initA(){Amin = FieldMin; Amax = FieldMax;}
		void initB(){Bmin = FieldMin; Bmax = FieldMax;}
		void initOut(){Outmin = FieldMin; Outmax = FieldMax;}


		size_t MaxDelayedDim(double beta){return std::max (0.0, floor ( (MaxStorableValue - abs (beta)*std::max (-Cmin, Cmax) ) / (std::max (-Amin, Amax) * std::max (-Bmin, Bmax))));}

		void setOutBounds(const size_t k, const double alpha, const double beta){
			if (beta<0){
				Outmin = beta*Cmax;
				Outmax = beta*Cmin;
			} else {
				Outmin = beta*Cmin;
				Outmax = beta*Cmax;
			}
			if (alpha >0){
				Outmin += k*alpha*std::min(Amin*Bmax, Amax*Bmin);
				Outmax += k*alpha*std::max(Amin*Bmin, Amax*Bmax);
			}else{
				Outmin += k*alpha*std::max(Amin*Bmin, Amax*Bmax);
				Outmax += k*alpha*std::min(Amin*Bmax, Amax*Bmin);
			}
		}

		bool checkA(const Field& F, const size_t M, const size_t N,
			    const typename Field::Element* A, const size_t lda ){
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
			    const typename Field::Element* A, const size_t lda ){
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
			      const typename Field::Element* A, const size_t lda ){
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
		MMHelper(const Field& F, size_t m, size_t k, size_t n) :
		                recLevel(Protected::WinogradSteps (F, min3(m,k,n))),
				FieldMin(getFieldMin(F)), FieldMax(getFieldMax(F)),
				Amin(FieldMin), Amax(FieldMax),
				Bmin(FieldMin), Bmax(FieldMax),
				Cmin(FieldMin), Cmax(FieldMax),
				MaxStorableValue ((1ULL << Protected::Mantissa<typename DelayedField_t::Element>())-1),
				delayedField(F.characteristic()) {}

		MMHelper(const Field& F, int w) :
				recLevel(w), //base(FflasDouble),
				FieldMin(getFieldMin(F)), FieldMax(getFieldMax(F)),
				Amin(FieldMin), Amax(FieldMax),
				Bmin(FieldMin), Bmax(FieldMax),
				Cmin(FieldMin), Cmax(FieldMax),
				MaxStorableValue ((1ULL << Protected::Mantissa<typename DelayedField_t::Element>())-1),
				delayedField(F.characteristic()){}

		// copy constructor from other Field and Algo Traits
		template<class F2, typename AlgoT2, typename FT2>
		MMHelper(MMHelper<F2, AlgoT2, FT2>& WH) :
		                recLevel(WH.recLevel),
				FieldMin(WH.FieldMin), FieldMax(WH.FieldMax),
				Amin(WH.Amin), Amax(WH.Amax),
				Bmin(WH.Bmin), Bmax(WH.Bmax),
				Cmin(WH.Cmin), Cmax(WH.Cmax),
				Outmin(WH.Outmin), Outmax(WH.Outmax),
				MaxStorableValue(WH.MaxStorableValue),
				delayedField(WH.delayedField) {}

		MMHelper(const Field& F, int w,
			 double _Amin, double _Amax,
			 double _Bmin, double _Bmax,
			 double _Cmin, double _Cmax):
				recLevel(w), FieldMin(getFieldMin(F)), FieldMax(getFieldMax(F)),
				Amin(_Amin), Amax(_Amax),
				Bmin(_Bmin), Bmax(_Bmax),
				Cmin(_Cmin), Cmax(_Cmax),
				MaxStorableValue((1ULL << Protected::Mantissa<typename DelayedField_t::Element>())-1),
				delayedField(F.characteristic()) {}

		void print() const {
			std::cerr<<"Helper: "
				 <<typeid(AlgoTrait).name()<<" "
				 <<typeid(FieldTrait).name()<<std::endl
				 <<"  recLevel = "<<recLevel<<std::endl
				 <<"  FieldMin = "<<FieldMin<<" FieldMax = "<<FieldMax<<std::endl
				 <<"  Amin = "<<Amin<<" Amax = "<<Amax<<std::endl
				 <<"  Bmin = "<<Bmin<<" Bmax = "<<Bmax<<std::endl
				 <<"  Cmin = "<<Cmin<<" Cmax = "<<Cmax<<std::endl
				 <<"  Outmin = "<<Outmin<<" Outmax = "<<Outmax<<std::endl;
		}
	}; // MMHelper
} // FFLAS
#endif
