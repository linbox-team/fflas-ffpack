/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   <bastien.vialla@lirmm.fr>
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

#ifndef __FFLASFFPACK_limits_H
#define __FFLASFFPACK_limits_H

//#include <cstddef>
#include <climits>
#include <limits>
#include <type_traits>

#include <givaro/givinteger.h>

template <class T> struct limits;
// {
//     constexpr inline static T max() noexcept {return 0;}
//     constexpr inline static T min() noexcept {return 0;}
// };

template <> struct limits<unsigned char> {
    typedef unsigned char T ;
    constexpr inline static unsigned char max() noexcept {return UCHAR_MAX;}
    constexpr inline static unsigned char min() noexcept {return 0;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<signed char> {
    typedef signed char T ;
    constexpr inline static signed char max() noexcept {return SCHAR_MAX;}
    constexpr inline static signed char min() noexcept {return SCHAR_MIN;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<char> {
    typedef char T ;
    constexpr inline static char max() noexcept {return CHAR_MAX;}
    constexpr inline static char min() noexcept {return CHAR_MIN;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<unsigned short int> {
    typedef unsigned short int T ;
    constexpr inline static unsigned short int max() noexcept {return USHRT_MAX;}
    constexpr inline static unsigned short int min() noexcept {return 0;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<short int> {
    typedef short int T ;
    constexpr inline static short int max() noexcept {return SHRT_MAX;}
    constexpr inline static short int min() noexcept {return SHRT_MIN;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<unsigned int> {
    typedef unsigned int T ;
    constexpr inline static unsigned int max() noexcept {return UINT_MAX;}
    constexpr inline static unsigned int min() noexcept {return 0;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<int> {
    typedef int T ;
    constexpr inline static int max() noexcept {return INT_MAX;}
    constexpr inline static int min() noexcept {return INT_MIN;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<unsigned long> {
    typedef unsigned  long T ;
    constexpr inline static unsigned long max() noexcept {return ULONG_MAX;}
    constexpr inline static unsigned long min() noexcept {return 0;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<long> {
    typedef  long T ;
    constexpr inline static long max() noexcept {return LONG_MAX;}
    constexpr inline static long min() noexcept {return LONG_MIN;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<unsigned long long> {
    typedef unsigned long long T ;
    constexpr inline static unsigned long long max() noexcept {
        return ULLONG_MAX;
    }
    constexpr inline static unsigned long long min() noexcept {return 0;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<long long> {
    typedef long long T ;
    constexpr inline static long long max() noexcept {return LLONG_MAX;}
    constexpr inline static long long min() noexcept {return LLONG_MIN;}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<float> {
    typedef float T ;
    constexpr inline static int32_t max() noexcept {return (int32_t(1) << FLT_MANT_DIG) - 1;}
    constexpr inline static int32_t min() noexcept {return -((int32_t(1) << FLT_MANT_DIG) - 1);}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<double> {
    typedef double T;
    constexpr inline static int64_t max() noexcept {return (int64_t(1) << DBL_MANT_DIG) - 1;}
    constexpr inline static int64_t min() noexcept {return -((int64_t(1) << DBL_MANT_DIG) - 1);}
    constexpr inline static int32_t digits() noexcept {return std:: numeric_limits<T>::digits ;}
};

template <> struct limits<Givaro::Integer> {
    typedef Givaro::Integer T;
    constexpr inline static int max() noexcept {return -1;}
    constexpr inline static int min() noexcept {return 0;}
};

template <size_t K> struct limits<RecInt::ruint<K> > {
    typedef RecInt::ruint<K> T;
    constexpr inline static RecInt::ruint<K> max() noexcept {return RecInt::ruint<K>(-1);}
    constexpr inline static RecInt::ruint<K> min() noexcept {return 0;}
};

template <size_t K> struct limits<RecInt::rint<K> > {
    typedef RecInt::ruint<K> T;
    constexpr inline static RecInt::rint<K> max() noexcept {return RecInt::rint<K>(RecInt::ruint<K>(-1) >> 1u);}
    constexpr inline static RecInt::rint<K> min() noexcept {return max() + 1;}
};
// template <size_t K> struct limits<RecInt::rint<K> > {
// 	constexpr inline static RecInt::rint<K> max() noexcept {return RecInt::rint<K>(RecInt::ruint<K>(-1))/2;}

// 	constexpr inline static RecInt::rint<K> min() noexcept {return -RecInt::rint<K>(RecInt::ruint<K>(-1))/2;}
// };
// template <size_t K,size_t MG> struct limits<RecInt::rmint<K,MG> > {
// 	constexpr inline static RecInt::ruint<K> max() noexcept {return RecInt::rmint<K,MG>(-1);}

// 	constexpr inline static RecInt::ruint<K> min() noexcept {return 0;}
// };

/*
 * in_range, determine if an element e of type E fit in a type T
 */

template<class T, class E>
typename std::enable_if<std::is_signed<T>::value == std::is_signed<E>::value, bool>::type
in_range(E e)
{
    return (e >= limits<T>::min() && e <= limits<T>::max());
}

template<class T, class E>
typename std::enable_if<(std::is_signed<T>::value) && !(std::is_signed<E>::value), bool>::type
in_range(E e)
{
    return (e <= static_cast<E>(limits<T>::max()));
}

template<class T, class E>
typename std::enable_if<!(std::is_signed<T>::value) && (std::is_signed<E>::value), bool>::type
in_range(E e)
{
    return ((e >= 0) && (static_cast<T>(e) <= limits<T>::max()));
}



#endif /* _FFLASFFPACK_limits_H */
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
