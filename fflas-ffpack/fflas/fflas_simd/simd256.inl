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

#ifndef __FFLASFFPACK_fflas_ffpack_utils_simd256_INL
#define __FFLASFFPACK_fflas_ffpack_utils_simd256_INL

template <bool ArithType, bool Int, bool Signed, int Size> struct Simd256_impl;

#include "simd256_float.inl"
#include "simd256_double.inl"

#ifdef SIMD_INT
// Trop d'instructions SSE manquantes pour les int8_t

#if defined(__FFLASFFPACK_USE_AVX2)
#include "simd256_int16.inl"
#include "simd256_int32.inl"
#include "simd256_int64.inl"
#endif

#endif //#ifdef SIMD_INT

template <class T>
using Simd256 =
    Simd256_impl<std::is_arithmetic<T>::value, std::is_integral<T>::value, std::is_signed<T>::value, sizeof(T)>;

#endif // __FFLASFFPACK_fflas_ffpack_utils_simd256_INL
