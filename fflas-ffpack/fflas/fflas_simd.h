/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
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

#ifndef __FFLASFFPACK_utils_simd_H
#define __FFLASFFPACK_utils_simd_H

#define SIMD_INT 1


//#include <x86intrin.h>
//#include <immintrin.h> -> only define for AVX
#include "fflas-ffpack/utils/fflas_intrinsic.h"
#include <iostream>
#include <type_traits>
#include <limits>

#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/utils/debug.h"

#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <type_traits>

#include "givaro/givtypestring.h"

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define INLINE __attribute__((always_inline)) inline
#else
#define INLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define CONST __attribute__((const))
#else
#define CONST
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define PURE __attribute__((pure))
#else
#define PURE
#endif

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
namespace std { // Why? - A.B. 2015-04-30

    inline
    std::ostream &operator<<(std::ostream &o, const __m128 &v) {
        const float *vArray = (const float *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1];
        o << ',';
        o << vArray[2] << ',' << vArray[3];
        o << '>';
        return o;
    }

    inline
    std::ostream &operator<<(std::ostream &o, const __m128i &v) {
        const int64_t *vArray = (const int64_t *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1];
        o << '>';
        return o;
    }

    inline
    std::ostream &operator<<(std::ostream &o, const __m128d &v) {
        const double *vArray = (const double *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1];
        o << '>';
        return o;
    }
} // std

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
namespace std {

    inline
    std::ostream &operator<<(std::ostream &o, const __m256 &v) {
        const float *vArray = (const float *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1] << ',' << vArray[2] << ',' << vArray[3];
        o << ',';
        o << vArray[4] << ',' << vArray[5] << ',' << vArray[6] << ',' << vArray[7];
        o << '>';
        return o;
    }

    inline
    std::ostream &operator<<(std::ostream &o, const __m256i &v) {
        const int64_t *vArray = (const int64_t *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1] << ',' << vArray[2] << ',' << vArray[3];
        o << '>';
        return o;
    }

    inline
    std::ostream &operator<<(std::ostream &o, const __m256d &v) {
        const double *vArray = (const double *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1] << ',' << vArray[2] << ',' << vArray[3];
        o << '>';
        return o;
    }
} // std

#endif // __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

#ifdef __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS
namespace std {

    inline
    std::ostream &operator<<(std::ostream &o, const __m512 &v) {
        const float *vArray = (const float *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1] << ',' << vArray[2] << ',' << vArray[3];
        o << ',';
        o << vArray[4] << ',' << vArray[5] << ',' << vArray[6] << ',' << vArray[7];
        o << ',';
        o << vArray[8] << ',' << vArray[9] << ',' << vArray[10] << ',' << vArray[11];
        o << ',';
        o << vArray[12] << ',' << vArray[13] << ',' << vArray[14] << ',' << vArray[15];
        o << '>';
        return o;
    }

    inline
    std::ostream &operator<<(std::ostream &o, const __m512i &v) {
        const int64_t *vArray = (const int64_t *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1] << ',' << vArray[2] << ',' << vArray[3];
        o << ',';
        o << vArray[4] << ',' << vArray[5] << ',' << vArray[6] << ',' << vArray[7];
        o << '>';
        return o;
    }

    inline
    std::ostream &operator<<(std::ostream &o, const __m512d &v) {
        const double *vArray = (const double *)(&v);
        o << '<';
        o << vArray[0] << ',' << vArray[1] << ',' << vArray[2] << ',' << vArray[3];
        o << ',';
        o << vArray[4] << ',' << vArray[5] << ',' << vArray[6] << ',' << vArray[7];
        o << '>';
        return o;
    }
} // std

#endif // __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS

namespace FFLAS {
    template <class T> struct support_simd : public std::false_type {};

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
    template <> struct support_simd<float> : public std::true_type {};
    template <> struct support_simd<double> : public std::true_type {};
#ifdef SIMD_INT
    template <> struct support_simd<int64_t> : public std::true_type {};
    template <> struct support_simd<int32_t> : public std::true_type {};
    template <> struct support_simd<int16_t> : public std::true_type {};
#endif
#endif

} // FFLAS

#define NORML_MOD(C, P, NEGP, MIN, MAX, Q, T)                                                                      \
{                                                                                                                  \
    Q = greater(C, MAX);                                                                                           \
    T = lesser(C, MIN);                                                                                            \
    Q = vand(Q, NEGP);                                                                                             \
    T = vand(T, P);                                                                                                \
    Q = vor(Q, T);                                                                                                 \
    C = add(C, Q);                                                                                                 \
}

#define FLOAT_MOD(C, P, INVP, Q)                                                                                   \
{                                                                                                                  \
    Q = mul(C, INVP);                                                                                              \
    Q = floor(Q);                                                                                                  \
    C = fnmadd(C, Q, P);                                                                                           \
}

// to activate SIMD with integers
//#define SIMD_INT

template <class T> struct simdToType;

/*
 * is_simd trait
 */

template <class T> struct is_simd {
    static const constexpr bool value = false;
    using type = std::integral_constant<bool, false>;
};

// SSE
#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS) // SSE or better
#include "fflas-ffpack/fflas/fflas_simd/simd128.inl"

template <> struct simdToType<__m128d> { using type = double; };

template <> struct simdToType<__m128> { using type = float; };

template <> struct is_simd<__m128d> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};

template <> struct is_simd<__m128> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};

#ifdef SIMD_INT
template <> struct is_simd<__m128i> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};
#endif

#endif // SSE

// AVX
#if defined(__FFLASFFPACK_HAVE_AVX_INSTRUCTIONS) or defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
#include "fflas-ffpack/fflas/fflas_simd/simd256.inl"

template <> struct simdToType<__m256d> { using type = double; };

template <> struct simdToType<__m256> { using type = float; };

template <> struct is_simd<__m256d> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};

template <> struct is_simd<__m256> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};

#ifdef SIMD_INT
template <> struct is_simd<__m256i> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};
#endif
#endif // AVX

// AVX512F
#if defined(__FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS)
#include "fflas-ffpack/fflas/fflas_simd/simd512.inl"

template <> struct simdToType<__m512d> { using type = double; };

template <> struct simdToType<__m512> { using type = float; };

template <> struct is_simd<__m512d> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};

template <> struct is_simd<__m512> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};

#ifdef SIMD_INT
template <> struct is_simd<__m512i> {
    static const constexpr bool value = true;
    using type = std::integral_constant<bool, true>;
};
#endif
#endif // AVX512F

/*
 * Simd functors
 */

template<typename T>
struct NoSimd {
    /*
     * alias to 128 bit simd register
     */
    using vect_t = T*;

    /*
     * define the scalar type corresponding to the specialization
     */
    using scalar_t = T;

    /*
     *  number of scalar_t in a simd register
     */
    static const constexpr size_t vect_size = 1;

    /*
     *  alignement for scalar_t pointer to be loaded in a vect_t
     * [ available to have the same interface as SimdXXX classes ]
     */
    static const constexpr size_t alignment = static_cast<size_t>(Alignment::Normal);
    using aligned_allocator = AlignedAllocator<scalar_t, Alignment(alignment)>;
    using aligned_vector = std::vector<scalar_t, aligned_allocator>;

    /* To check compatibility with Modular struct */
    template <class Field>
    using is_same_element = std::is_same<typename Field::Element, T>;

    /* Name of the NoSimd struct */
    static inline const std::string type_string () {
        return "NoSimd<" + Givaro::TypeString<scalar_t>::get() + ">";
    }

    // Test if the pointer p is multiple of alignment
    template <class TT> static constexpr bool valid(TT p) { return false; }

    // Test if n is multiple of vect_size
    template <class TT> static constexpr bool compliant(TT n) { return false; }
};

// #if defined(__FFLASFFPACK_HAVE_AVX_INSTRUCTIONS)

template <class T, bool = std::is_arithmetic<T>::value, bool = std::is_integral<T>::value> struct SimdChooser {};

template <class T, bool b> struct SimdChooser<T, false, b> { using value = NoSimd<T>; };

template <class T>
struct SimdChooser<T, true, false> // floating number
{
#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    using value = Simd512<T>;
#elif defined (__FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS)
    using value = Simd256<T>;
#elif defined(__FFLASFFPACK_HAVE_AVX_INSTRUCTIONS)
    using value = Simd256<T>;
#elif defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
    using value = Simd128<T>;
#else
    using value = NoSimd<T>;
#endif
};

template <class T>
struct SimdChooser<T, true, true> // integral number
{
#ifdef __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS
    using value = Simd512<T>;
#elif __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
    using value = Simd256<T>;
#elif __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    using value = Simd128<T>;
#else
    using value = NoSimd<T>;
#endif
};

template <class T> using Simd = typename SimdChooser<T>::value;

// template <class T> struct SimdChooser<T, true> {
// #if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
//     typedef Simd256<T> value;
// #else
//     typedef Simd128<T> value;
// #endif // __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
// };

// #elif defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS) // not AVX

// template <class T> using Simd = Simd128<T>;

// #endif // __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS) // SSE or better

// template <class T> struct floating_simd;

// template <> struct floating_simd<float> { typedef Simd<float> value; };

// template <> struct floating_simd<double> { typedef Simd<double> value; };

// template <> struct floating_simd<int64_t> {
// #if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
// // typedef Simd256<double> value;
// #else
//     typedef Simd128<double> value;
// #endif
// };

#endif

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

namespace FFLAS { /*  print helper */

    // need friend ?
    template <class simdT>
    inline std::ostream &print(std::ostream &os, const typename simdT::vect_t &P) {
        typename simdT::scalar_t p[simdT::vect_size];
        os << '<';
        simdT::storeu(p, P);
        for (size_t i = 0; i < simdT::vect_size; ++i) {
            os << p[i];
            if (i < simdT::vect_size - 1)
                os << '|';
        }
        os << '>';

        return os;
    }

} // FFLAS

namespace std {
    // cannot be instanciated, T is not deductible
    template <class T>
    inline std::ostream &operator<<(std::ostream &o, const typename Simd128<T>::vect_t &v) {
        FFLAS::print<Simd128<T>>(o, v);
        return o;
    }
} // std

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
namespace std {
    // cannot be instanciated, T is not deductible
    template <class T>
    inline std::ostream &operator<<(std::ostream &o, const typename Simd256<T>::vect_t &v) {
        FFLAS::print(o, v);
        return o;
    }
}
#endif // __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

// Provide simd modular support
#include <fflas-ffpack/fflas/fflas_simd/simd_modular.inl>

#undef INLINE
#undef PURE
#undef CONST
#undef SIMD_INT

#endif /* __FFLASFFPACK_utils_simd_H */
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
