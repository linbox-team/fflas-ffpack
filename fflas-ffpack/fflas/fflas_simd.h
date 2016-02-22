/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifdef __FFLASFFPACK_USE_SIMD
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

#ifdef __FFLASFFPACK_USE_AVX
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
#endif // __FFLASFFPACK_USE_AVX

#endif // __FFLASFFPACK_USE_SIMD

namespace FFLAS {
template <class T> struct support_simd : public std::false_type {};

#if defined(__FFLASFFPACK_USE_SIMD)
template <> struct support_simd<float> : public std::true_type {};
template <> struct support_simd<double> : public std::true_type {};
#ifdef SIMD_INT
template <> struct support_simd<int64_t> : public std::true_type {};
template <> struct support_simd<int32_t> : public std::true_type {};
template <> struct support_simd<int16_t> : public std::true_type {};
#endif
#endif

} // FFLAS

#define NORML_MOD(C, P, NEGP, MIN, MAX, Q, T)                                                                          \
    {                                                                                                                  \
        Q = greater(C, MAX);                                                                                           \
        T = lesser(C, MIN);                                                                                            \
        Q = vand(Q, NEGP);                                                                                             \
        T = vand(T, P);                                                                                                \
        Q = vor(Q, T);                                                                                                 \
        C = add(C, Q);                                                                                                 \
    }

#define FLOAT_MOD(C, P, INVP, Q)                                                                                       \
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
#if defined(__FFLASFFPACK_USE_SIMD) // SSE or better
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
#if defined(__FFLASFFPACK_USE_AVX) or defined(__FFLASFFPACK_USE_AVX2)
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

/*
 * Simd functors
 */

struct NoSimd {
    // Test if the pointer p is multiple of alignment
    template <class T> static constexpr bool valid(T p) { return false; }

    // Test if n is multiple of vect_size
    template <class T> static constexpr bool compliant(T n) { return false; }
};

// #if defined(__FFLASFFPACK_USE_AVX)

template <class T, bool = std::is_arithmetic<T>::value, bool = std::is_integral<T>::value> struct SimdChooser {};

template <class T, bool b> struct SimdChooser<T, false, b> { using value = NoSimd; };

template <class T>
struct SimdChooser<T, true, false> // floating number
    {
#ifdef __FFLASFFPACK_USE_AVX
    using value = Simd256<T>;
#elif defined(__FFLASFFPACK_USE_SSE)
    using value = Simd128<T>;
#else
    using value = NoSimd;
#endif
};

template <class T>
struct SimdChooser<T, true, true> // integral number
    {
#ifdef __FFLASFFPACK_USE_AVX2
    using value = Simd256<T>;
#elif __FFLASFFPACK_USE_SSE
    using value = Simd128<T>;
#else
    using value = NoSimd;
#endif
};

template <class T> using Simd = typename SimdChooser<T>::value;

// template <class T> struct SimdChooser<T, true> {
// #if defined(__FFLASFFPACK_USE_AVX2)
//     typedef Simd256<T> value;
// #else
//     typedef Simd128<T> value;
// #endif // __FFLASFFPACK_USE_AVX2
// };

// #elif defined(__FFLASFFPACK_USE_SSE) // not AVX

// template <class T> using Simd = Simd128<T>;

// #endif // __FFLASFFPACK_USE_AVX

#if defined(__FFLASFFPACK_USE_SIMD) // SSE or better

// template <class T> struct floating_simd;

// template <> struct floating_simd<float> { typedef Simd<float> value; };

// template <> struct floating_simd<double> { typedef Simd<double> value; };

// template <> struct floating_simd<int64_t> {
// #if defined(__FFLASFFPACK_USE_AVX2)
// // typedef Simd256<double> value;
// #else
//     typedef Simd128<double> value;
// #endif
// };

#endif

#ifdef __FFLASFFPACK_USE_SIMD

namespace FFLAS { /*  print helper */

// need friend ?
template <class simdT>
inline std::ostream &print(std::ostream &os, const typename simdT::vect_t &P) {
    typename simdT::scalar_t p[simdT::vect_size];
    os << '<';
    simdT::store(p, P);
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
// cannot be instanciated, T is not déductible
template <class T>
inline std::ostream &operator<<(std::ostream &o, const typename Simd128<T>::vect_t &v) {
    FFLAS::print<Simd128<T>>(o, v);
    return o;
}
} // std

#ifdef __FFLASFFPACK_USE_AVX
namespace std {
// cannot be instanciated, T is not déductible
template <class T>
inline std::ostream &operator<<(std::ostream &o, const typename Simd256<T>::vect_t &v) {
    FFLAS::print(o, v);
    return o;
}
}
#endif // __FFLASFFPACK_USE_AVX

#endif // __FFLASFFPACK_USE_SIMD

#undef INLINE
#undef PURE
#undef CONST
#undef SIMD_INT

#endif /* __FFLASFFPACK_utils_simd_H */
