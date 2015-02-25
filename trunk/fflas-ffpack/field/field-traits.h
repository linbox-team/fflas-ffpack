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
		//! generic ring.
		struct GenericTag{};
		//! This is a <code>Modular<T></code> or <code>ModularBalanced<T></code>
		struct ModularTag : public GenericTag {} ;
		//! This is an <code>UnparametricRing<T></code>
		struct UnparametricTag : public GenericTag {} ;

		//! If it is a machine floating point (ie \c float or \c double)
		struct FloatingPointTag : public GenericTag{};
		//! If it is a machine int (ie  \c intX_t or \c uintX_t)
		struct IntegralTag : public GenericTag{};
		//! If it is a multiprecision field (ie \c Givaro::Integer)
		struct MultiPrecisionTag : public  GenericTag{};
		//- If it can support SIMD operations (ie \c double or \c int32_t, etc)
		// struct SIMDTag : public GenericTag{};

		// this is weird :
		//! If it can init/convert elements to/from floating point types: float, double
		struct FloatingPointConvertibleTag : public  GenericTag{};
		//! If it is a Modular or ModularBalanced templated by float or double
		struct ModularFloatingPointTag : public GenericTag{};
		//! If it is a Modular or ModularBalanced templated by float or double, and result is not reduced
		struct DelayedModularFloatingPointTag : public GenericTag{};
	}

} // FFLAS

namespace FFLAS { /*  Traits */

	/*! FieldTrait
	*/
	template <class Field>
	struct FieldTraits {
		typedef typename FieldCategories::GenericTag value;
		typedef typename FieldCategories::GenericTag category;
		typedef typename FieldCategories::GenericTag rep_t ;
		// typedef false_type balanced ;
		static  const bool balanced = false ;
	};
	

	// RecInt
	template<size_t K>
	struct FieldTraits<Givaro::UnparametricRing<RecInt::ruint<K>> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	
	template<size_t K, int MG>
	struct FieldTraits<Givaro::UnparametricRing<RecInt::rmint<K, MG>> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};

	// Modular <double|float>
	// ModularBalanced <double|float>
	template<>
	struct FieldTraits<Givaro::Modular<double> > {
		typedef  FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<float> > {
		typedef FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<double> > {
		typedef FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<float> > {
		typedef FieldCategories::ModularFloatingPointTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};


	//Givaro::Modular< intX >
	// ModularBalanced < intX >
	template<>
	struct FieldTraits<Givaro::Modular<int16_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<uint16_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<int32_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<uint32_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<int64_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::Modular<uint64_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<int32_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<uint32_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<int64_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<uint64_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<int16_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};
	template<>
	struct FieldTraits<Givaro::ModularBalanced<uint16_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};

	//Givaro::Modular< intX >
	// ModularBalanced < intX >
	template<>
	struct FieldTraits<Givaro::Modular<Givaro::Integer> > {
		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::MultiPrecisionTag rep_t ;
		static  const bool balanced = false ;
	};

	template<>
	struct FieldTraits<Givaro::ModularBalanced<Givaro::Integer> > {
		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::MultiPrecisionTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = true ;
	};


	// ZRing< float|double >
	template<>
	struct FieldTraits<Givaro::ZRing<double>> {
		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::ZRing<float>> {
		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		static  const bool balanced = false ;
	};


	// UnparametricRing< intX >
	template<>
	struct FieldTraits<Givaro::UnparametricRing<int16_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<uint16_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<int32_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<uint32_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<int64_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<uint64_t> > {
		typedef FieldCategories::FloatingPointConvertibleTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::IntegralTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<double>> {
		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		static  const bool balanced = false ;
	};
	template<>
	struct FieldTraits<Givaro::UnparametricRing<float>> {
		typedef FieldCategories::FloatingPointTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::FloatingPointTag rep_t ;
		static  const bool balanced = false ;
	};

	// UnparametricRing<Integer>
	template<>
	struct FieldTraits<Givaro::UnparametricRing<Givaro::Integer> >
	{
		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::UnparametricTag category;
		typedef typename FieldCategories::MultiPrecisionTag rep_t ;
		static  const bool balanced = false ;
	};

	// RNSInteger
	template<typename T>
	struct FieldTraits<FFPACK::RNSInteger<T> > {
		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::MultiPrecisionTag rep_t ;
		// typedef true_type balanced ;
		static  const bool balanced = false ;
	};
	// RNSIntegerMod
	template<typename T>
	struct FieldTraits<FFPACK::RNSIntegerMod<T> >{
		typedef FieldCategories::MultiPrecisionTag value;
		typedef FieldCategories::ModularTag category;
		typedef typename FieldCategories::MultiPrecisionTag rep_t ;
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
	struct associatedDelayedField<const Givaro::Modular<Givaro::Integer>>{
		typedef Givaro::UnparametricRing<Givaro::Integer> field;
		typedef Givaro::UnparametricRing<Givaro::Integer> type;
	};
	template <typename RNS>
	struct associatedDelayedField<const FFPACK::RNSIntegerMod<RNS> >{
		typedef FFPACK::RNSInteger<RNS> field;
		typedef FFPACK::RNSInteger<RNS> type;
	};


} // FFLAS

#endif // __FFLASFFPACK_field_field_traits_H

