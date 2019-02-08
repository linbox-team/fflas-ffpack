/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by  Bastien Vialla <bastien.vialla@lirmm.fr>
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 * ========LICENCE========
 *.
 */

#ifndef __FFLASFFPACK_SPARSEMATRIX_TRAITS_H
#define __FFLASFFPACK_SPARSEMATRIX_TRAITS_H

#include <type_traits>

namespace FFLAS {

    /****************************************************************************************************************
     *
     *  SparseMatrix Traits
     *
     ****************************************************************************************************************/

    template <class Field, class M> struct isSparseMatrix : public std::false_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR_ZO>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::COO>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::COO_ZO>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_ZO>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::SELL>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::SELL_ZO>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_simd>> : public std::true_type {};

    template <class Field>
    struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_simd_ZO>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR_HYB>> : public std::true_type {};

    template <class Field> struct isSparseMatrix<Field, Sparse<Field, SparseMatrix_t::HYB_ZO>> : public std::true_type {};


    template <class F, class M> struct isZOSparseMatrix : public std::false_type {};

    template <class Field> struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::CSR_ZO>> : public std::true_type {};

    template <class Field> struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::COO_ZO>> : public std::true_type {};

    template <class Field> struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_ZO>> : public std::true_type {};

    template <class Field>
    struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::SELL_ZO>> : public std::true_type {};

    template <class Field>
    struct isZOSparseMatrix<Field, Sparse<Field, SparseMatrix_t::ELL_simd_ZO>> : public std::true_type {};

    using ZOSparseMatrix = std::true_type;
    using NotZOSparseMatrix = std::false_type;


    template<class F, class M> struct isSparseMatrixSimdFormat : public std::false_type {};

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    template<class Field> struct isSparseMatrixSimdFormat<Field, Sparse<Field, SparseMatrix_t::SELL>> : public support_simd<typename Field::Element>::type {};

    template<class Field> struct isSparseMatrixSimdFormat<Field, Sparse<Field, SparseMatrix_t::ELL_simd>> : public support_simd<typename Field::Element>::type {};

#endif // __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS

    using SimdSparseMatrix = std::true_type;
    using NoSimdSparseMatrix = std::false_type;


    template<class F, class M> struct isSparseMatrixMKLFormat : public std::false_type {};

#ifdef __FFLASFFPACK_HAVE_MKL

    template<class Field> struct isSparseMatrixMKLFormat<Field, Sparse<Field, SparseMatrix_t::CSR>> : public std::true_type {};
    template<class Field> struct isSparseMatrixMKLFormat<Field, Sparse<Field, SparseMatrix_t::CSC>> : public std::true_type {};

#endif // __FFLASFFPACK_HAVE_MKL

    using MKLSparseMatrixFormat = std::true_type;
    using NotMKLSparseMatrixFormat = std::false_type;

    /********************************************************************************************************
     *
     *	Traits to test if operator +, -, *, =, +=, -=, *= exists
     *
     ********************************************************************************************************/

#if 0
#define function_to_functor(X) 					 \
    struct tfn_##X { 							 \
        template <typename... Args> 				 \
        auto operator()(Args&&... args) const 		 \
        ->decltype(X(std::forward<Args>(args)...)){ \
            return X(std::forward<Args>(args)...); } }
#endif

    struct tfn_plus {
        template <typename... Args>
        auto operator()(Args&&... args) const ->decltype(operator+(std::forward<Args>(args)...))
        {
            return operator+(std::forward<Args>(args)...);
        }
    };

    struct tfn_mul {
        template <typename... Args>
        auto operator()(Args&&... args) const ->decltype(operator+(std::forward<Args>(args)...))
        {
            return operator*(std::forward<Args>(args)...);
        }
    };

    struct tfn_mul_eq {
        template <typename... Args>
        auto operator()(Args&&... args) const ->decltype(operator+(std::forward<Args>(args)...))
        {
            return operator*=(std::forward<Args>(args)...);
        }
    };

    struct tfn_minus {
        template <typename... Args>
        auto operator()(Args&&... args) const ->decltype(operator+(std::forward<Args>(args)...))
        {
            return operator-(std::forward<Args>(args)...);
        }
    };

    struct tfn_plus_eq {
        template <typename... Args>
        auto operator()(Args&&... args) const ->decltype(operator+(std::forward<Args>(args)...))
        {
            return operator+=(std::forward<Args>(args)...);
        }
    };

    struct tfn_minus_eq {
        template <typename... Args>
        auto operator()(Args&&... args) const ->decltype(operator+(std::forward<Args>(args)...))
        {
            return operator+=(std::forward<Args>(args)...);
        }
    };

    template<typename C>
    struct has_plus_impl {
    private:
        // Test for non member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<typename std::result_of<tfn_plus(const T&, const T&)>::type, T>::type;

        // Test for member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<decltype(std::declval<T>().operator+(std::declval<T>())), T>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:
        static constexpr bool value = type::value;
    };

    template<typename C>
    struct has_mul_impl {
    private:
        // Test for non member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<typename std::result_of<tfn_mul(const T&, const T&)>::type, T>::type;

        // Test for member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<decltype(std::declval<T>().operator*(std::declval<T>())), T>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:
        static constexpr bool value = type::value;
    };

    template<typename C>
    struct has_mul_eq_impl {
    private:
        // Test for non member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<typename std::result_of<tfn_mul_eq(const T&, const T&)>::type, T>::type;

        // Test for member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<decltype(std::declval<T>().operator*=(std::declval<T>())), typename std::add_lvalue_reference<T>::type>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:
        static constexpr bool value = type::value;
    };

    template<typename C>
    struct has_plus_eq_impl {
    private:
        // Test for non member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<typename std::result_of<tfn_plus_eq(const T&, const T&)>::type, T>::type;

        // Test for member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<decltype(std::declval<T>().operator+=(std::declval<T>())), typename std::add_lvalue_reference<T>::type>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:
        static constexpr bool value = type::value;
    };

    template<typename C>
    struct has_minus_eq_impl {
    private:
        // Test for non member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<typename std::result_of<tfn_minus_eq(const T&, const T&)>::type, T>::type;

        // Test for member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<decltype(std::declval<T>().operator-=(std::declval<T>())), typename std::add_lvalue_reference<T>::type>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:
        static constexpr bool value = type::value;
    };

    template<typename C>
    struct has_minus_impl {
    private:
        // Test for non member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<typename std::result_of<tfn_minus(const T&, const T&)>::type, T>::type;

        // Test for member operator
        template<typename T>
        static constexpr auto check(T*)
        -> typename std::is_same<decltype(std::declval<T>().operator-(std::declval<T>())), T>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:
        static constexpr bool value = type::value;
    };

    template<class T>
    using has_plus = typename std::conditional<std::is_arithmetic<T>::value, std::true_type, has_plus_impl<T>>::type;

    template<class T>
    using has_minus = typename std::conditional<std::is_arithmetic<T>::value, std::true_type, has_minus_impl<T>>::type;

    template<class T>
    using has_equal = typename std::conditional<std::is_arithmetic<T>::value, std::true_type, std::is_copy_assignable<T>>::type;

    template<class T>
    using has_plus_eq = typename std::conditional<std::is_arithmetic<T>::value, std::true_type, has_plus_eq_impl<T>>::type;

    template<class T>
    using has_minus_eq = typename std::conditional<std::is_arithmetic<T>::value, std::true_type, has_minus_eq_impl<T>>::type;

    template<class T>
    using has_mul = typename std::conditional<std::is_arithmetic<T>::value, std::true_type, has_mul_impl<T>>::type;

    template<class T>
    using has_mul_eq = typename std::conditional<std::is_arithmetic<T>::value, std::true_type, has_mul_eq_impl<T>>::type;

    template<class T>
    struct has_operation{
        static constexpr bool value = (has_plus<T>::value && has_minus<T>::value && has_equal<T>::value &&
                                       has_plus_eq<T>::value && has_minus_eq<T>::value && has_mul<T>::value && has_mul_eq<T>::value);
    };

} // FFLAS

#endif // __FFLASFFPACK_SPARSEMATRIX_TRAITS_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
