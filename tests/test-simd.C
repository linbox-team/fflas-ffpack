/*
 * Copyright (C) 2014 FFLAS-FFPACK
 * Written by :
 *        Bastien Vialla <bastien.vialla@lirmm.fr>
 * This file is Free Software and part of FFLAS-FFPACK.
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

#include "givaro/givinteger.h" /* for Givaro::Integer */
#include "givaro/givprint.h" /* for operator<< with vector */
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas/fflas_simd.h"
#include "fflas-ffpack/utils/args-parser.h" /* for parsing command-line args */
#include "fflas-ffpack/utils/test-utils.h" /* for FFLAS::getSeed */
#include "fflas-ffpack/utils/align-allocator.h"

#include <array>
#include <vector>
#include <random>
#include <string>
#include <functional>
#include <limits>
#include <type_traits>
#include <algorithm>

typedef Givaro::Integer integer;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::array;
using std::function;
using std::numeric_limits;
using std::enable_if;
using std::is_floating_point;
using std::is_integral;
using std::is_signed;
using std::equal;

/* For pretty printing type */
template<typename...> const char *TypeName();

#define REGISTER_TYPE_NAME(type) \
    template<> const char *TypeName<type>(){return #type;}

REGISTER_TYPE_NAME(float);
REGISTER_TYPE_NAME(double);
REGISTER_TYPE_NAME(int16_t);
REGISTER_TYPE_NAME(int32_t);
REGISTER_TYPE_NAME(int64_t);
REGISTER_TYPE_NAME(uint16_t);
REGISTER_TYPE_NAME(uint32_t);
REGISTER_TYPE_NAME(uint64_t);

/******************************************************************************/
/* Random generators **********************************************************/
/******************************************************************************/

static std::mt19937 entropy_generator;

template <class Element, class Alloc>
typename enable_if<is_integral<Element>::value>::type
generate_random_vector (vector<Element,Alloc> &a) {
    typedef typename std::uniform_int_distribution<Element> RandGen;
    RandGen G(numeric_limits<Element>::lowest(),numeric_limits<Element>::max());
    std::generate (a.begin(), a.end(), [&](){return G(entropy_generator);});
}

template <class Element, class Alloc>
typename enable_if<is_floating_point<Element>::value>::type
generate_random_vector (vector<Element,Alloc> &a) {
    typedef typename std::uniform_real_distribution<Element> RandGen;
    RandGen G(numeric_limits<Element>::lowest(),numeric_limits<Element>::max());
    std::generate (a.begin(), a.end(), [&](){return G(entropy_generator);});
}

/******************************************************************************/
/* Utils functions ************************************************************/
/******************************************************************************/

/* check equality for integral type */
template <class Element>
typename enable_if<is_integral<Element>::value, bool>::type
check_eq (Element x, Element y)
{
    return x == y;
}

/* check equality for floating point type */
template <class Element>
typename enable_if<is_floating_point<Element>::value, bool>::type
check_eq (Element x, Element y)
{
    return (std::isnan(x) && std::isnan(y)) || x == y;
}

/* evaluate the function f with arguments taken in the array */
template <class Ret, class T>
Ret
eval_func_on_array (function<Ret()> f, array<T, 0> arr)
{
    return f();
}

template <class Ret, class T, class...TArgs>
Ret
eval_func_on_array (function<Ret(T, TArgs...)> f,
                    array<T, sizeof...(TArgs)+1> arr)
{
    function<Ret(TArgs...)> newf = [&] (TArgs...args) -> Ret { return f(arr[0], args...);};
    array<T, sizeof...(TArgs)> newarr;
    for (size_t i = 0; i < sizeof...(TArgs); i++)
        newarr[i] = arr[i+1];
    return eval_func_on_array (newf, newarr);
}


/******************************************************************************/
/* Main test function *********************************************************/
/******************************************************************************/

template <class Simd, class RScal, class...AScal, class RSimd, class...ASimd>
typename enable_if<sizeof...(AScal) == sizeof...(ASimd), bool>::type
test_op (RSimd (&FSimd) (ASimd...), RScal (&FScal) (AScal...), string fname) {

    using Element = typename Simd::scalar_t;
    using ScalVectAlign = AlignedAllocator<Element, Alignment(Simd::alignment)>;
    using ScalVect = vector<Element, ScalVectAlign>;
    using SimdVect = typename Simd::vect_t;
    constexpr size_t SimdVectSize = Simd::vect_size;
    constexpr size_t simd_size = sizeof(Element) * SimdVectSize * 8;
    constexpr size_t arity = sizeof...(AScal);

    /* input vectors */
    vector<ScalVect> inputs (arity, ScalVect(SimdVectSize));
    for (auto &iv: inputs)
        generate_random_vector (iv);

    /* output vectors */
    ScalVect out_scal(SimdVectSize), out_simd(SimdVectSize);

    /* compute with scalar function */
    array<Element, arity> scal_in;
    function<RScal(AScal...)> fscal = FScal;
    for(size_t i = 0 ; i < SimdVectSize ; i++) {
        for (size_t j = 0; j < arity; j++)
            scal_in[j] = inputs[j][i];

        out_scal[i] = eval_func_on_array (fscal, scal_in);
    }

    /* compute with SIMD function */
    array<SimdVect, arity> simd_in;
    function<RSimd(ASimd...)> fsimd = FSimd;
    for (size_t i = 0; i < arity; i++)
        simd_in[i] = Simd::load (inputs[i].data());

    SimdVect simd_out = eval_func_on_array (fsimd, simd_in);
    Simd::store (out_simd.data(), simd_out);

    /* comparison */
    auto eq = check_eq<Element>;
    bool res = equal (out_scal.begin(), out_scal.end(), out_simd.begin(), eq);

    /* print result line */
    cout << "Simd" << simd_size << "<" << TypeName<Element>() << ">::" << fname
         << " " << string (60 - fname.size() - strlen(TypeName<Element>()), '.')
         << " " << (res ? "success" : "failure") << endl;

    /* in case of error, print all input and output values */
    if(!res) {
        cout << string (10, '-') << " debug data " << string (58, '-') << endl;
        for (size_t i = 0; i < arity; i++) {
            cout << "input_" << i << ": " << inputs[i] << endl;
        }
        cout << "out_scal: " << out_scal << endl;
        cout << "out_simd: " << out_simd << endl;
        cout << string (80, '-') << endl;
    }
    return res;
}

/******************************************************************************/
/* Scalar functions for comparisons *******************************************/
/******************************************************************************/

template <class Element, class Enable = void>
struct ScalFunctions;

/* for floating point element */
template <class Element>
struct ScalFunctions<Element,
                    typename enable_if<is_floating_point<Element>::value>::type>
{
    static Element ceil (Element x) {
        return std::ceil(x);
    }
    static Element floor (Element x) {
        return std::floor(x);
    }
    static Element round (Element x) {
        return std::round(x);
    }
    static Element add (Element x1, Element x2) {
        return x1+x2;
    }
    static Element sub (Element x1, Element x2) {
        return x1-x2;
    }
    static Element mul (Element x1, Element x2) {
        return x1*x2;
    }
    static Element fmadd (Element x1, Element x2, Element x3) {
        return std::fma(x3,x2,x1);
    }
    static Element fmsub (Element x1, Element x2, Element x3) {
        return std::fma(x3,x2,-x1);
    }
    static Element fnmadd (Element x1, Element x2, Element x3) {
        return std::fma(-x3,x2,x1);
    }
    /* Comparisons functions in SIMD output 0 or 0xFFFF...FFFF
     * (here we assume 0xFFFF...FFFF is always a NAN)
     */
    static Element lesser (Element x1, Element x2) {
        return (x1<x2)?NAN:0;
    }
    static Element lesser_eq (Element x1, Element x2) {
        return (x1<=x2)?NAN:0;
    }
    static Element greater (Element x1, Element x2) {
        return (x1>x2)?NAN:0;
    }
    static Element greater_eq (Element x1, Element x2) {
        return (x1>=x2)?NAN:0;
    }
    static Element eq (Element x1, Element x2) {
        return (x1==x2)?NAN:0;
    }
};

/* for integral element */
template <class Element>
struct ScalFunctions<Element,
                    typename enable_if<is_integral<Element>::value>::type>
{
    static Element add (Element x1, Element x2) {
        return x1+x2;
    }
    static Element sub (Element x1, Element x2) {
        return x1-x2;
    }
    static Element mul (Element x1, Element x2) {
        return x1*x2;
    }
    static Element mullo (Element x1, Element x2) {
        return x1*x2;
    }
    static Element mulhi (Element x1, Element x2) {
        integer q,r;
        integer a = (integer(x1)*integer(x2));
        integer b = integer(1) << uint64_t(sizeof(Element)*8);
        Givaro::IntegerDom Z;
        Z.divmod(q, r, a, b);
        return Element(q);
    }
    static Element mulx (Element x1, Element x2) {
        /* h = 1 << (half the number of bits of Element) */
        Element h = Element(1) << (sizeof(Element)*4);

        /* Representative r of x1 modulo h with -h/2 <= r < h/2 */
        if (std::is_signed<Element>::value) {
            x1 = (x1+h/2) % h;
            x1 += (x1 < 0) ? h/2 : -h/2;
            x2 = (x2+h/2) % h;
            x2 += (x2 < 0) ? h/2 : -h/2;
        }
        else {
            x1 = x1 % h;
            x2 = x2 % h;
        }
        return x1*x2;
    }
    static Element fmadd (Element x1, Element x2, Element x3) {
        return x1+x3*x2;
    }
    static Element fmsub (Element x1, Element x2, Element x3) {
        return -x1 + x3*x2;
    }
    static Element fnmadd (Element x1, Element x2, Element x3) {
        return x1 - x3*x2;
    }

    /* Shift */
    template <int s, bool EnableTrue = true>
    static
    typename enable_if<!is_signed<Element>::value && EnableTrue, Element>::type
    sra (Element x1) {
        return x1 >> s; /* For unsigned type, simply use >> */
    }

    template <int s, bool EnableTrue = true>
    static
    typename enable_if<is_signed<Element>::value && EnableTrue, Element>::type
    sra (Element x1) {
        /* For signed type we need to do a sign extension, the code comes from
         *   http://graphics.stanford.edu/~seander/bithacks.html#FixedSignExtend
         */
        struct {Element x:sizeof(Element)*8-s;} r;
        return r.x = (x1 >> s);
    }

    template <int s>
    static Element srl (Element x1) {
        return ((typename std::make_unsigned<Element>::type) x1) >> s;
    }

    template <int s>
    static Element sll (Element x1) {
        return ((typename std::make_unsigned<Element>::type) x1) << s;
    }

    /* Comparisons functions in SIMD output 0 or 0xFFFF...FFFF */
    static Element lesser (Element x1, Element x2) {
        return (x1<x2)?-1:0;
    }
    static Element lesser_eq (Element x1, Element x2) {
        return (x1<=x2)?-1:0;
    }
    static Element greater (Element x1, Element x2) {
        return (x1>x2)?-1:0;
    }
    static Element greater_eq (Element x1, Element x2) {
        return (x1>=x2)?-1:0;
    }
    static Element eq (Element x1, Element x2) {
        return (x1==x2)?-1:0;
    }
};

/******************************************************************************/
/* Test one SIMD implem *******************************************************/
/******************************************************************************/

#define TEST_ONE_OP(name) \
    btest &= test_op<simd,Element> (simd::name, Scal::name, #name);

/* for floating point element */
template<class simd, class Element>
typename enable_if<is_floating_point<Element>::value, bool>::type
test_impl () {
    using Scal = ScalFunctions<Element>;
    bool btest = true;

    TEST_ONE_OP (ceil);
    TEST_ONE_OP (floor);
    TEST_ONE_OP (round);
    TEST_ONE_OP (add);
    TEST_ONE_OP (sub);
    TEST_ONE_OP (mul);
    TEST_ONE_OP (fmadd);
    TEST_ONE_OP (fmsub);
    TEST_ONE_OP (fnmadd);
    TEST_ONE_OP (lesser);
    TEST_ONE_OP (lesser_eq);
    TEST_ONE_OP (greater);
    TEST_ONE_OP (greater_eq);
    TEST_ONE_OP (eq);

    return btest;
}

/* for integral element */
template<class simd, class Element>
typename enable_if<is_integral<Element>::value, bool>::type
test_impl () {
    using Scal = ScalFunctions<Element>;
    bool btest = true;

    TEST_ONE_OP (add);
    TEST_ONE_OP (sub);
    TEST_ONE_OP (mul);
    TEST_ONE_OP (mullo);
    TEST_ONE_OP (mulhi);
    TEST_ONE_OP (mulx);
    TEST_ONE_OP (fmadd);
    TEST_ONE_OP (fmsub);
    TEST_ONE_OP (fnmadd);
    TEST_ONE_OP (lesser);
    TEST_ONE_OP (lesser_eq);
    TEST_ONE_OP (greater);
    TEST_ONE_OP (greater_eq);
    TEST_ONE_OP (eq);
    TEST_ONE_OP (template sra<3>);
    TEST_ONE_OP (template sra<7>);
    TEST_ONE_OP (template srl<5>);
    TEST_ONE_OP (template srl<11>);
    TEST_ONE_OP (template sll<2>);
    TEST_ONE_OP (template sll<13>);

    return btest;
}

/******************************************************************************/
/* Test all SIMD implems for one Element type *********************************/
/******************************************************************************/
template<class Element>
typename enable_if<is_integral<Element>::value, bool>::type
test () {
    bool test = true;

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    test &= test_impl<Simd128<Element>, Element>();
    cout << endl;
#endif

#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
    test &= test_impl<Simd256<Element>, Element>();
    cout << endl;
#endif

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    test &= test_impl<Simd512<Element>, Element>();
    cout << endl;
#endif

    return test;
}

template<class Element>
typename enable_if<is_floating_point<Element>::value, bool>::type
test () {
    bool test = true;

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    test &= test_impl<Simd128<Element>, Element>();
    cout << endl;
#endif

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
    test &= test_impl<Simd256<Element>, Element>();
    cout << endl;
#endif

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    test &= test_impl<Simd512<Element>, Element>();
    cout << endl;
#endif

    return test;
}

/******************************************************************************/
/* Main ***********************************************************************/
/******************************************************************************/
int
main (int argc, char *argv[]) {
    uint64_t seed = FFLAS::getSeed();

    static Argument args[] = {
        { 's', "-s S", "Set the seed", TYPE_UINT64 , &seed },
        END_OF_ARGUMENTS
    };

    FFLAS::parseArguments (argc, argv, args);

    cout << "# To rerun this test: test-simd -s " << seed << endl;
    cout << "# seed = " << seed << endl << endl;

    entropy_generator.seed (seed);

    bool pass  = true ;
    pass &= test<float>();
    pass &= test<double>();
#ifndef __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS
        // Not yet implemented over AVX512
    pass &= test<int16_t>();
    pass &= test<int32_t>();
#endif
#ifdef __x86_64__
    pass &= test<int64_t>();
#endif
#ifndef __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS
        // Not yet implemented over AVX512
    pass &= test<uint16_t>();
    pass &= test<uint32_t>();
#endif
#ifdef __x86_64__
    pass &= test<uint64_t>();
#endif
    cout << endl << "Test " << (pass ? "passed" : "failed") << endl;
    return pass ? 0 : 1;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
