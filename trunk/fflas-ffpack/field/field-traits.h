/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * BB <brice.boyer@lip6.fr>
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

/** @file field/field-traits.h
 * @brief  Field Traits
 */

#ifndef __FFLASFFPACK_field_field_traits_H
#define __FFLASFFPACK_field_field_traits_H

#include <type_traits> // CXX11
#include "fflas-ffpack/field/rns-double-elt.h"
#include "recint/rmint.h"
// ----- Forward declarations

namespace Givaro {

	template<class T>
	class Modular;

	template<class T>
	class ModularBalanced;

	template<class T>
	class UnparametricRing;

	template<class T>
	class ZRing;

}

namespace FFPACK {

	template<class T>
	class RNSInteger;

	template<class T>
	class RNSIntegerMod;

}

namespace FFLAS { /*  Categories */

	//! Traits and categories will need to be placed in a proper file later

	namespace FieldCategories {
	
		    // Classify
		//! generic ring.
		struct GenericTag{};
		//! This is a <code>Modular<T></code> or <code>ModularBalanced<T></code>
		struct ModularTag{};
		//! If it is a multiprecision field (ie \c Givaro::Integer)
		struct MultiPrecisionTag{};
		//! If the field uses a representation with infix operators
		struct UnparametricTag{};
	}

	    //! Specifies the mode of action for an algorithm w.r.t. its field
	    //! 
	namespace ModeCategories {
		    //! No specific mode of action: use standard field operations
		struct DefaultTag{};
		    //! Force conversion to appropriate element type of ElementCategory T.
		    //! e.g. 
		    //!    - ConvertTo<ElementCategories::MachineFloatTag> tries conversion 
		    //!      of Modular<int> to Modular<double>
		    //!    - ConvertTo<ElementCategories::FixedPrecIntTag> tries conversion 
                    //!      of Modular<Integer> to Modular<RecInt<K> >
		    //!    - ConvertTo<ElementCategories::ArbitraryPrecIntTag> tries conversion 
		    //!     of Modular<Integer> to RNSInteger
		    //! .
		template<class T>
		struct ConvertTo{};
		    //! Performs field operations with delayed mod reductions. Ensures result is reduced.
		struct DelayedTag{};
		    //! Performs field operations with delayed mod only when necessary. Result may not be reduced.
		struct LazyTag{};
	}

	namespace ElementCategories {
		    //! default is generic
		struct GenericTag{};
		    //! float or double
		struct MachineFloatTag{};
		    //! short, int, long, long long, and unsigned variants
		struct MachineIntTag{};
		    //! Fixed precision integers above machine precision: Givaro::recInt
		struct FixedPrecIntTag{};
		    //! Arbitrary precision integers: GMP
		struct ArbitraryPrecIntTag{};
		    //! Representation in a Residue Number System 
		struct RNSElementTag{};
		//- If it can support SIMD operations (ie \c double or \c int32_t, etc)
		// struct SIMDTag : public GenericTag{};
	}

} // FFLAS

namespace FFLAS { /*  Traits */

	/*! ElementTraits
	*/
	template <class Element>
	struct ElementTraits {typedef typename ElementCategories::GenericTag value;};
	
	template<> struct ElementTraits<float> {typedef ElementCategories::MachineFloatTag value;};
	template<> struct ElementTraits<double> {typedef ElementCategories::MachineFloatTag value;};
	template<> struct ElementTraits<int8_t> {typedef ElementCategories::MachineIntTag value;};
	template<> struct ElementTraits<int16_t> {typedef ElementCategories::MachineIntTag value;};
	template<> struct ElementTraits<int32_t> {typedef ElementCategories::MachineIntTag value;};
	template<> struct ElementTraits<int64_t> {typedef ElementCategories::MachineIntTag value;};
	template<> struct ElementTraits<uint8_t> {typedef ElementCategories::MachineIntTag value;};
	template<> struct ElementTraits<uint16_t> {typedef ElementCategories::MachineIntTag value;};
	template<> struct ElementTraits<uint32_t> {typedef ElementCategories::MachineIntTag value;};
	template<> struct ElementTraits<uint64_t> {typedef ElementCategories::MachineIntTag value;};
	template<>
	struct ElementTraits<Givaro::Integer> {typedef ElementCategories::ArbitraryPrecIntTag value;};
	template<size_t K> 
	struct ElementTraits<RecInt::ruint<K> > {typedef ElementCategories::FixedPrecIntTag value;};
	template<size_t K, int MG> 
	struct ElementTraits<RecInt::rmint<K, MG> >{typedef ElementCategories::FixedPrecIntTag value;};
	template<>
	struct ElementTraits<FFPACK::rns_double_elt>{typedef ElementCategories::RNSElementTag value;};


	/*! ModeTraits
	*/
	template <class Field> 
	struct ModeTraits {typedef typename ModeCategories::DefaultTag value;};
	template <> struct ModeTraits<Givaro::Modular<float> >{typedef typename ModeCategories::DelayedTag value;};
	template <> struct ModeTraits<Givaro::Modular<double> > {typedef typename ModeCategories::DelayedTag value;};
	template <> struct ModeTraits<Givaro::Modular<int8_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<int16_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<int32_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<int64_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<uint8_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<uint16_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<uint32_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<uint64_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::Modular<Givaro::Integer> > {typedef typename ModeCategories::ConvertTo<ElementCategories::RNSElementTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<float> >{typedef typename ModeCategories::DelayedTag value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<double> > {typedef typename ModeCategories::DelayedTag value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<int8_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<int16_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<int32_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<int64_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<uint8_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<uint16_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<uint32_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<uint64_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
	template <> struct ModeTraits<Givaro::ModularBalanced<Givaro::Integer> > {typedef typename ModeCategories::ConvertTo<ElementCategories::RNSElementTag> value;};
	template <> struct ModeTraits<Givaro::UnparametricRing<Givaro::Integer> > {typedef typename ModeCategories::ConvertTo<ElementCategories::RNSElementTag> value;};

	/*! FieldTrait
	*/
	template <class Field>
	struct FieldTraits {
		typedef typename FieldCategories::GenericTag category;
		// typedef false_type balanced ;
		static  const bool balanced = false ;
	};
	

	// RecInt
	template<size_t K>
	struct FieldTraits<Givaro::UnparametricRing<RecInt::ruint<K> > > {
		    //typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	
	template<size_t K, int MG>
	struct FieldTraits<Givaro::UnparametricRing<RecInt::rmint<K, MG> > > {
		    //typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};

	// Modular <double|float>
	// ModularBalanced <double|float>
	template<>
	struct FieldTraits<Givaro::Modular<double> > {
//		typedef  FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<float> > {
//		typedef FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<double> > {
//		typedef FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<float> > {
//		typedef FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};


	//Givaro::Modular< intX >
	// ModularBalanced < intX >
	template<>
	struct FieldTraits<Givaro::Modular<int16_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<uint16_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<int32_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<uint32_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<int64_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<uint64_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<int32_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<uint32_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<int64_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<uint64_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<int16_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<uint16_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};

	//Givaro::Modular< intX >
	// ModularBalanced < intX >
	template<>
	struct FieldTraits<Givaro::Modular<Givaro::Integer> > {
//		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		static  const bool balanced = false ;
	};

	template<>
	struct FieldTraits<Givaro::ModularBalanced<Givaro::Integer> > {
//		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};


	// ZRing< float|double >
	template<>
	struct FieldTraits<Givaro::ZRing<double>> {
//		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::ZRing<float>> {
//		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};

	// UnparametricRing< intX >
	template<>
	struct FieldTraits<Givaro::UnparametricRing<int16_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<uint16_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<int32_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<uint32_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<int64_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<uint64_t> > {
//		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<double> > {
//		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<float> > {
//		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};

	// UnparametricRing<Integer>
	template<>
	struct FieldTraits<Givaro::UnparametricRing<Givaro::Integer> >
	{
//		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::UnparametricTag category;
		static  const bool balanced = false ;
	};

	// RNSInteger
	template<typename T>
	struct FieldTraits<FFPACK::RNSInteger<T> > {
//		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = false ;
	};
	// RNSIntegerMod
	template<typename T>
	struct FieldTraits<FFPACK::RNSIntegerMod<T> >{
//		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		// typedef true_type balanced ;
		static  const bool balanced = false ;
	};


} // FFLAS

namespace FFLAS { /* associatedDelayedField */

	template <class Field>
	struct associatedDelayedField{
		typedef Field field;
		typedef Field& type; // reference to avoid copying heavy fields
	};
	template <>
	struct associatedDelayedField<const Givaro::Modular<float> >{
		typedef Givaro::UnparametricRing<float> field;
		typedef Givaro::UnparametricRing<float> type;
	};
	template <>
	struct associatedDelayedField<const Givaro::ModularBalanced<float> >{
		typedef Givaro::UnparametricRing<float> field;
		typedef Givaro::UnparametricRing<float> type;
	};
	template <>
	struct associatedDelayedField<const Givaro::Modular<double> >{
		typedef Givaro::UnparametricRing<double> field;
		typedef Givaro::UnparametricRing<double> type;
	};
	template <>
	struct associatedDelayedField<const Givaro::ModularBalanced<double> >{
		typedef Givaro::UnparametricRing<double> field;
		typedef Givaro::UnparametricRing<double> type;
	};
	template <>
	struct associatedDelayedField<const Givaro::Modular<Givaro::Integer> >{
		typedef Givaro::UnparametricRing<Givaro::Integer> field;
		typedef Givaro::UnparametricRing<Givaro::Integer> type;
	};
	template <typename RNS>
	struct associatedDelayedField<const FFPACK::RNSIntegerMod<RNS> >{
		typedef FFPACK::RNSInteger<RNS> field;
		typedef FFPACK::RNSInteger<RNS> type;
	};
	template <size_t K, int MG>
	struct associatedDelayedField<const Givaro::Modular<RecInt::rmint<K, MG> > >{
		typedef Givaro::UnparametricRing<RecInt::rmint<K, MG> > field;
		typedef Givaro::UnparametricRing<RecInt::rmint<K, MG> > type;
	};
	template <size_t K, int MG>
	struct associatedDelayedField<const Givaro::ModularBalanced<RecInt::rmint<K, MG> > >{
		typedef Givaro::UnparametricRing<RecInt::rmint<K, MG> > field;
		typedef Givaro::UnparametricRing<RecInt::rmint<K, MG> > type;
	};


} // FFLAS

#endif // __FFLASFFPACK_field_field_traits_H

