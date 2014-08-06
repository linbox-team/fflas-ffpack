/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifndef _FFLASFFPACK_limits_H
#define _FFLASFFPACK_limits_H

#include <cstddef>
#include <climits>
#include <type_traits>

template <class T> struct limits
{
    constexpr inline static T max() noexcept {return 0;}
    constexpr inline static T min() noexcept {return 0;}
};

template <> struct limits<unsigned char> {
  constexpr inline static unsigned char max() noexcept { return UCHAR_MAX; }
  constexpr inline static unsigned char min() noexcept {return 0;}
};

template <> struct limits<signed char> {
  constexpr inline static signed char max() noexcept { return SCHAR_MAX; }
  constexpr inline static signed char min() noexcept { return SCHAR_MIN;}
};

template <> struct limits<char> {
  constexpr inline static char max() noexcept { return CHAR_MAX; }
  constexpr inline static char min() noexcept {return CHAR_MIN;}
};

template <> struct limits<unsigned short int> {
  constexpr inline static unsigned short int max() noexcept {
    return USHRT_MAX;
  }
  constexpr inline static unsigned short int min() noexcept {return 0;}
};

template <> struct limits<short int> {
  constexpr inline static short int max() noexcept { return SHRT_MAX; }
  constexpr inline static short int min() noexcept {return SHRT_MIN;}
};

template <> struct limits<unsigned int> {
  constexpr inline static unsigned int max() noexcept { return UINT_MAX; }
  constexpr inline static unsigned int min() noexcept {return 0;}
};

template <> struct limits<int> {
  constexpr inline static int max() noexcept { return INT_MAX; }
  constexpr inline static int min() noexcept {return INT_MIN;}
};

template <> struct limits<unsigned long> {
  constexpr inline static unsigned long max() noexcept { return ULONG_MAX; }
  constexpr inline static unsigned long min() noexcept {return 0;}
};

template <> struct limits<long> {
  constexpr inline static long max() noexcept { return LONG_MAX; }
  constexpr inline static long min() noexcept {return LONG_MIN;}
};

template <> struct limits<unsigned long long> {
  constexpr inline static unsigned long long max() noexcept {
    return ULLONG_MAX;
  }
  constexpr inline static unsigned long long min() noexcept {return 0;}
};

template <> struct limits<long long> {
  constexpr inline static long long max() noexcept { return LLONG_MAX; }
  constexpr inline static long long min() noexcept {return LLONG_MIN;}
};

template <> struct limits<float> {
  constexpr inline static int32_t max() noexcept { return (1 << 23) - 1; }
  constexpr inline static int32_t min() noexcept {return -((1 << 23) - 1); }
};

template <> struct limits<double> {
  constexpr inline static int64_t max() noexcept {
    return (uint64_t(1) << 53) - 1;
  }
  constexpr inline static int64_t min() noexcept {return -((uint64_t(1) << 53) - 1);}
};

/*
 * in_range, determine if an element e of type E fit in a type T
 */

template<class T, class E>
typename std::enable_if<std::is_signed<T>::value == std::is_signed<E>::value, bool>::type
in_range(E e)
{
    return (e >= std::limits<T>::min() && e <= std::limits<T>::max());
}

template<class T, class E>
typename std::enable_if<(std::is_signed<T>::value) && !(std::is_signed<E>::value), bool>::type
in_range(E e)
{
    return (e <= static_cast<E>(std::limits<T>::max()));
}

template<class T, class E>
typename std::enable_if<!(std::is_signed<T>::value) && (std::is_signed<E>::value), bool>::type
in_range(E e)
{
    return ((e >= 0) && (static_cast<T>(e) <= std::limits<T>::max()));
}



#endif /* _FFLASFFPACK_limits_H */
