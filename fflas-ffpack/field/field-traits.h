/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

// ----- Forward declarations

#include "recint/rmint.h"
#include "givaro/modular-general.h"
#include "givaro/zring.h"

namespace RecInt {

    template<size_t K>
    class rint;

    template<size_t K>
    class ruint;

}

namespace Givaro {

    template<class T>
    class ModularBalanced;

    template<class T>
    class Montgomery;

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
        //! This is a modular field like e.g. <code>Modular<T></code> or <code>ModularBalanced<T></code>
        struct ModularTag{};
        //! If the field uses a representation with infix operators
        struct UnparametricTag{};
    }

    //! Specifies the mode of action for an algorithm w.r.t. its field
    //!
    namespace ModeCategories {
        //! No specific mode of action: use standard field operations
        struct DefaultTag{};

        //! Use standard field operations, but keeps track of bounds on input and output
        struct DefaultBoundedTag{};

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
    struct ElementTraits<RecInt::rint<K> > {typedef ElementCategories::FixedPrecIntTag value;};
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

    template <typename Element, typename Compute>
    struct ModeTraits<Givaro::Modular<Element,Compute> >{typedef typename ModeCategories::DelayedTag value;};
    template<> struct ModeTraits<Givaro::Modular<int64_t,uint64_t> > {typedef typename ModeCategories::DefaultTag value;};

    template<typename Compute> struct ModeTraits<Givaro::Modular<int8_t,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
    template<typename Compute> struct ModeTraits<Givaro::Modular<int16_t,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
    template<typename Compute> struct ModeTraits<Givaro::Modular<int32_t,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
    template<typename Compute> struct ModeTraits<Givaro::Modular<uint8_t,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
    template<typename Compute> struct ModeTraits<Givaro::Modular<uint16_t,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
    template<typename Compute> struct ModeTraits<Givaro::Modular<uint32_t,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};

#ifndef INTEGER_NO_RNS
    template<typename Compute> struct ModeTraits<Givaro::Modular<Givaro::Integer,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::RNSElementTag> value;};
    template<typename Compute, size_t K> struct ModeTraits<Givaro::Modular<RecInt::ruint<K>,Compute> > {typedef typename ModeCategories::ConvertTo<ElementCategories::RNSElementTag> value;};
#endif

    template <typename Element>
    struct ModeTraits<Givaro::ModularBalanced<Element> >{typedef typename ModeCategories::DelayedTag value;};

    template <> struct ModeTraits<Givaro::ModularBalanced<int8_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
    template <> struct ModeTraits<Givaro::ModularBalanced<int16_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
    template <> struct ModeTraits<Givaro::ModularBalanced<int32_t> > {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};

#ifndef INTEGER_NO_RNS
    template <> struct ModeTraits<Givaro::ModularBalanced<Givaro::Integer> > {typedef typename ModeCategories::ConvertTo<ElementCategories::RNSElementTag> value;};
    template <> struct ModeTraits<Givaro::ZRing<Givaro::Integer> > {typedef typename ModeCategories::ConvertTo<ElementCategories::RNSElementTag> value;};
#endif

    // These ones are here temporarily, to ensure
    // In the long term ZRing should be in DefaultTag, and forced to be in DefaultBoundedTag be the caller. However this would prevent these rings to use Winograd's algorithm (extensive use of bounded helpers) in the current implementation. Needs work.
    template <> struct ModeTraits<Givaro::ZRing<float> > {typedef typename ModeCategories::DefaultBoundedTag value;};
    template <> struct ModeTraits<Givaro::ZRing<double> > {typedef typename ModeCategories::DefaultBoundedTag value;};
    template <class T> struct ModeTraits<Givaro::Montgomery<T> > {typedef typename ModeCategories::DefaultBoundedTag value;};

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
    struct FieldTraits<Givaro::ZRing<RecInt::ruint<K> > > {
        //typedef FieldCategories::FloatingPointConvertibleTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };

    // Modular <double|float>
    // ModularBalanced <double|float>
    template<class Element>
    struct FieldTraits<Givaro::Modular<Element> > {
        typedef FieldCategories::ModularTag category;
        static  const bool balanced = false ;
    };

    template<class Element>
    struct FieldTraits<Givaro::ModularBalanced<Element> > {
        typedef FieldCategories::ModularTag category;
        static  const bool balanced = true ;
    };

    // ZRing< float|double >
    template<>
    struct FieldTraits<Givaro::ZRing<double> > {
        //		typedef FieldCategories::FloatingPointTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };
    template<>
    struct FieldTraits<Givaro::ZRing<float> > {
        //		typedef FieldCategories::FloatingPointTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };

    // ZRing< intX >
    template<>
    struct FieldTraits<Givaro::ZRing<int16_t> > {
        //		typedef FieldCategories::FloatingPointConvertibleTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };
    template<>
    struct FieldTraits<Givaro::ZRing<uint16_t> > {
        //		typedef FieldCategories::FloatingPointConvertibleTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };
    template<>
    struct FieldTraits<Givaro::ZRing<int32_t> > {
        //		typedef FieldCategories::FloatingPointConvertibleTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };
    template<>
    struct FieldTraits<Givaro::ZRing<uint32_t> > {
        //		typedef FieldCategories::FloatingPointConvertibleTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };
    template<>
    struct FieldTraits<Givaro::ZRing<int64_t> > {
        //		typedef FieldCategories::FloatingPointConvertibleTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };
    template<>
    struct FieldTraits<Givaro::ZRing<uint64_t> > {
        //		typedef FieldCategories::FloatingPointConvertibleTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };

    // ZRing<Integer>
    template<>
    struct FieldTraits<Givaro::ZRing<Givaro::Integer> >
    {
        //		typedef FieldCategories::MultiPrecisionTag value;
        typedef FieldCategories::UnparametricTag category;
        static  const bool balanced = false ;
    };

    // RNSInteger
    template<typename T>
    struct FieldTraits<FFPACK::RNSInteger<T> > {
        //		typedef FieldCategories::MultiPrecisionTag value;
        typedef FieldCategories::UnparametricTag category;
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
    template <typename T,typename X>
    struct associatedDelayedField<const Givaro::Modular<T,X>> {
        typedef Givaro::ZRing<T> field;
        typedef Givaro::ZRing<T> type;
    };
    template <typename T>
    struct associatedDelayedField<const Givaro::ModularBalanced<T>> {
        typedef Givaro::ZRing<T> field;
        typedef Givaro::ZRing<T> type;
    };
    template <typename T>
    struct associatedDelayedField<const Givaro::ZRing<T>> {
        typedef Givaro::ZRing<T> field;
        typedef Givaro::ZRing<T> type;
    };
    template <typename RNS>
    struct associatedDelayedField<const FFPACK::RNSIntegerMod<RNS>> {
        typedef FFPACK::RNSInteger<RNS> field;
        typedef FFPACK::RNSInteger<RNS> type;
    };

} // FFLAS

namespace FFLAS { /* MaxCadinality */
    template <class Field, class enable=void>
    inline typename Field::Residu_t maxCardinality() {return Field::maxCardinality();}

// Need to override Givaro's default, as Compute_t (uint64_t) is larger than Storage_t
    template<>
    inline uint64_t maxCardinality <Givaro::Modular<int64_t> >(){
        // ceil(2^31.5) such that ab+cd fits in uint64_t with int64_t a,b,c,d and  abs(a,b,c,d) <= (p-1)
        return UINT64_C(3037000500);
    }
    template<>
    inline uint32_t maxCardinality<Givaro::Modular<int32_t> >(){
        // ceil(2^15.5) such that ab+cd fits in uint32_t with a,b,c,d int32_t and abs(a,b,c,d) <= (p-1)
        return UINT32_C(46341);
    }
    template<class Field, typename std::enable_if<is_rint<typename Field::Element>::value || is_ruint<typename Field::Element>::value, typename Field::Element>::type>
    typename Field::Residu_t maxCardinality (){
            // such that ab+cd fits in Field::Element with a,b,c,d Field::Element and abs(a,b,c,d) <= (p-1)
        return typename Field::Element::maxFFLAS();
    }

    template <class Field>
    inline typename Field::Residu_t minCardinality() {return Field::minCardinality();}

} // FFLAS

#endif // __FFLASFFPACK_field_field_traits_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
