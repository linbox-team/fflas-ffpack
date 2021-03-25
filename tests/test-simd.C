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
#include "givaro/modular.h"
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

using Givaro::Integer;
using Givaro::Modular;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::array;
using std::function;
using std::enable_if;
using std::is_floating_point;
using std::is_integral;
using std::decay;
using std::uniform_real_distribution;
using std::uniform_int_distribution;

/* Single source of entropy to use in all random generators */
static std::mt19937 entropy_generator;

/******************************************************************************/
/* Utils structs for enable_if ************************************************/
/******************************************************************************/
/* Testing if all bool, given as template, are true */
template <bool...v>
struct ALL;

template <bool...v>
struct ALL<true, v...> { static constexpr bool value = ALL<v...>::value; } ;

template <bool...v>
struct ALL<false, v...> { static constexpr bool value = false; };

template <>
struct ALL<> { static constexpr bool value = true; };

/* Counting the number of lvalue reference in parameter pack */
template <typename...T>
struct count_nonconst_lvalue_reference;

template <typename T, typename...O>
struct count_nonconst_lvalue_reference<T, O...> {
    static constexpr size_t n = count_nonconst_lvalue_reference<O...>::n;
};

template <typename T, typename...O>
struct count_nonconst_lvalue_reference<T&, O...> {
    static constexpr size_t n = std::integral_constant<size_t, 1>::value
                                    + count_nonconst_lvalue_reference<O...>::n;
};

template <typename T, typename...O>
struct count_nonconst_lvalue_reference<const T&, O...> {
    static constexpr size_t n = count_nonconst_lvalue_reference<O...>::n;
};

template <>
struct count_nonconst_lvalue_reference<> {
    static constexpr size_t n = std::integral_constant<size_t, 0>::value;
};

/*
 * Check that all types given in the parameter pack Args are the same as T, once
 * the constness and the references are removed.
 */
template <typename...Args>
struct is_all_same;

template <typename T, typename...Args>
struct is_all_same<T, Args...> { static constexpr bool value = ALL<std::is_same<typename decay<T>::type, typename decay<Args>::type>::value...>::value; };

template <>
struct is_all_same<> { static constexpr bool value = true; };



/******************************************************************************/
/* Utils functions and structs ************************************************/
/******************************************************************************/
/* check equality for integral type */
template <typename Element>
typename enable_if<is_integral<Element>::value, bool>::type
check_eq (Element x, Element y)
{
    return x == y;
}

/* check equality for floating point type */
template <typename Element>
typename enable_if<is_floating_point<Element>::value, bool>::type
check_eq (Element x, Element y)
{
    return (std::isnan(x) && std::isnan(y)) || x == y;
}

/* Check if two vectors are equal (using check_eq) */
template <typename Element>
bool
cmp (vector<Element> out_scal, vector<Element> out_simd)
{
    auto eq = check_eq<Element>;
    return std::equal (out_scal.begin(), out_scal.end(), out_simd.begin(), eq);
}

/* evaluate the function f with arguments taken in the array */
template <typename Ret, typename T>
Ret
eval_func_on_array (function<Ret()> f, array<T, 0>& arr)
{
    return f();
}

template <typename T, typename...TArgs>
void
eval_func_on_array (function<void(T, TArgs...)> f,
                    array<typename decay<T>::type, sizeof...(TArgs)+1> &arr)
{
    function<void(TArgs...)> newf = [&] (TArgs...args) -> void { return f(arr[0], args...);};
    array<typename decay<T>::type, sizeof...(TArgs)> newarr;
    std::copy (std::next(arr.begin()), arr.end(), newarr.begin());
    eval_func_on_array (newf, newarr);
    std::copy (newarr.begin(), newarr.end(), std::next(arr.begin()));
}

template <typename Ret, typename T, typename...TArgs>
Ret
eval_func_on_array (function<Ret(T, TArgs...)> f,
                    array<typename decay<T>::type, sizeof...(TArgs)+1> &arr)
{
    function<Ret(TArgs...)> newf = [&] (TArgs...args) -> Ret { return f(arr[0], args...);};
    array<typename decay<T>::type, sizeof...(TArgs)> newarr;
    std::copy (std::next(arr.begin()), arr.end(), newarr.begin());
    Ret r = eval_func_on_array (newf, newarr);
    std::copy (newarr.begin(), newarr.end(), std::next(arr.begin()));
    return r;
}

/******************************************************************************/
/* Utils for pretty printing vectors ******************************************/
/******************************************************************************/
template <typename T>
struct width { static constexpr size_t value = 2+2*sizeof(T); };

template <>
struct width<float> { static constexpr size_t value = 16; };
template <>
struct width<double> { static constexpr size_t value = 24; };

/* pretty printing vectors */
template <typename E>
std::ostream& operator<< (std::ostream& o, const vector<E>& V)
{
    const std::ios::fmtflags old_settings = o.flags();
    const char prev = o.fill();
    if (is_floating_point<E>::value)
        o << std::hexfloat << std::right;
    else if (is_integral<E>::value)
        o << std::hex << std::showbase << std::internal << std::setfill ('0');
    for (size_t i = 0; i < V.size(); i++)
    {
        o << (i ? ", " : "[ ");
        if (is_integral<E>::value && V[i] == 0)
            o << "0x" << std::setw(width<E>::value-2) << V[i];
        else
            o << std::setw(width<E>::value) << V[i];
    }
    o << " ]";
    o.fill(prev);
    o.flags(old_settings);
    return o;
}

/******************************************************************************/
/* Main test function *********************************************************/
/******************************************************************************/
/* Class to perform the test a given method from a Simd struct against a method
 * in the ScalFunctions struct. It can handle the following cases:
 *  - Arguments of the Simd method are vect_t or any type built from vect_t
 *      with const and references. If any, the arguments with non constant
 *      references must appear first.
 *  - The return value of the Simd method is either void or a type that can be
 *      converted into vect_t. [cf the templated method evaluate_simd_method]
 *  - Arguments of the scalar method are all of type Element (or resp.
 *      vectElt) or any type built from Element (or resp. vectElt) with
 *      const and references. If any, the arguments with non constant references
 *      must appear first.
 *  - The return value of the scalar method is Element (if all arguments are
 *  Elements) or vectElt or void (if all arguments are vectElt). [cf the
 *      templated method evaluate_scalar_method]
 */
template <typename Simd>
class TestOneMethod
{
public:
    using Element = typename Simd::scalar_t;
    using vect_t = typename Simd::vect_t;
    using vectElt = vector<Element>;
    constexpr static size_t vect_size = Simd::vect_size;

    template <bool B, typename T = void>
    using enable_if_t = typename enable_if<B, T>::type;

    /*** Constructor **********************************************************/
    template <typename...AScal, typename RScal,
              typename...ASimd, typename RSimd,
              enable_if_t<sizeof...(AScal) == sizeof...(ASimd)>* = nullptr,
              enable_if_t<count_nonconst_lvalue_reference<AScal...>::n == count_nonconst_lvalue_reference<ASimd...>::n>* = nullptr,
              enable_if_t<is_all_same<AScal...>::value>* = nullptr,
              enable_if_t<is_all_same<vect_t, ASimd...>::value>* = nullptr>
    TestOneMethod (function<RSimd(ASimd...)> fsimd,
                   function<RScal(AScal...)> fscal,
                   function<void(vector<vectElt> &)> genInputs, string fname)
                                                                : name(fname){
        /* Constants are computed with AScal, but could have been with ASimd */
        constexpr size_t arity = sizeof...(AScal);
        nb_lref = count_nonconst_lvalue_reference<AScal...>::n;
        constexpr bool is_return_void = std::is_same<RScal, void>::value;

        inputs.resize (arity, vectElt(vect_size));
        outputs_simd.resize (nb_lref+(is_return_void?0:1), vectElt(vect_size));
        outputs_scalar.resize(nb_lref+(is_return_void?0:1), vectElt(vect_size));

        genInputs (inputs);

        /* compute with scalar function */
        evaluate_scalar_method (fscal);

        /* compute with SIMD function */
        array<vect_t, arity> simd_in;
        /* convert input into vect_t */
        for (size_t i = 0; i < inputs.size(); i++)
            simd_in[i] = Simd::loadu (inputs[i].data());
        /* evaluate simd function */
        evaluate_simd_method (fsimd, simd_in);
        /* get back arguments passed as lvalue reference */
        for (size_t j = 0; j < nb_lref; j++)
            Simd::storeu (outputs_simd[j].data(), simd_in[j]);
    }

    /*** To evaluation scalar method ******************************************/
    /* Ret (Element...), with Ret convertible to Element */
    template <typename Ret, typename...AScal>
    enable_if_t<is_all_same<Element, AScal...>::value
                            && std::is_convertible<Ret, Element>::value, void>
    evaluate_scalar_method (function<Ret(AScal...)> fscal) {
        array<Element, sizeof...(AScal)> scal_in;
        for(size_t i = 0 ; i < Simd::vect_size; i++) {
            for (size_t j = 0; j < inputs.size(); j++)
                scal_in[j] = inputs[j][i];
            /* evaluate scalar function */
            outputs_scalar[nb_lref][i] = eval_func_on_array (fscal, scal_in);
            /* get back arguments passed as lvalue reference */
            for (size_t j = 0; j < nb_lref; j++)
                outputs_scalar[j][i] = scal_in[j];
        }
    }

    /* vectElt (vectElt...) */
    template <typename...AScal>
    enable_if_t<is_all_same<vectElt, AScal...>::value, void>
    evaluate_scalar_method (function<vectElt(AScal...)> fscal) {
        array<vectElt, sizeof...(AScal)> scal_in;
        std::copy (inputs.begin(), inputs.end(), scal_in.begin());
        /* evaluate scalar function */
        outputs_scalar[nb_lref] = eval_func_on_array (fscal, scal_in);
        /* get back arguments passed as lvalue reference */
        for (size_t j = 0; j < nb_lref; j++)
            outputs_scalar[j] = scal_in[j];
    }

    /* void (vectElt...) */
    template <typename...AScal>
    enable_if_t<is_all_same<vectElt, AScal...>::value, void>
    evaluate_scalar_method (function<void(AScal...)> fscal) {
        array<vectElt, sizeof...(AScal)> scal_in;
        std::copy (inputs.begin(), inputs.end(), scal_in.begin());
        /* evaluate scalar function */
        eval_func_on_array (fscal, scal_in);
        /* get back arguments passed as lvalue reference */
        for (size_t j = 0; j < nb_lref; j++)
            outputs_scalar[j] = scal_in[j];
    }

    /*** To evaluation Simd method ********************************************/
    /* Rec (vect_t...), with Ret convertible to vect_t */
    template <typename Ret, typename...ASimd>
    enable_if_t<is_all_same<vect_t, ASimd...>::value
                            && std::is_convertible<Ret, vect_t>::value, void>
    evaluate_simd_method (function<Ret(ASimd...)> fsimd,
                                    array<vect_t, sizeof...(ASimd)>& simd_in) {
        vect_t simd_out = eval_func_on_array (fsimd, simd_in);
        /* store the results */
        Simd::storeu (outputs_simd[nb_lref].data(), simd_out);
    }

    /* void (vect_t...) */
    template <typename...ASimd>
    enable_if_t<is_all_same<vect_t, ASimd...>::value, void>
    evaluate_simd_method (function<void(ASimd...)> fsimd,
                                    array<vect_t, sizeof...(ASimd)>& simd_in) {
        eval_func_on_array (fsimd, simd_in);
    }

    /*** Utils ****************************************************************/
    bool getStatus () const {
        bool status = true;
        for (size_t i = 0; i < outputs_scalar.size(); i++)
            status &= cmp (outputs_scalar[i], outputs_simd[i]);
        return status;
    }

    string getTestName () const {
        return name;
    }

    bool writeResultLine () const {
        bool status = getStatus();
        cout << Simd::type_string() << "::" << getTestName() << " "
             << string (69-getTestName().size()-Simd::type_string().size(), '.')
             << " " << (status ? "success" : "failure") << endl;
        return status;
    }

    void writeDebugData () const {
        const std::ios::fmtflags old_settings = cout.flags();

        auto eq = check_eq<Element>;
        const size_t w = width<Element>::value;
        const size_t w2 = 16;
        const size_t fw = vect_size*(w+2)+2+w2;
        const string h("debug data");

        cout << string ((fw-h.size()-2)/2, '#') << " " << h << " "
             << string ((fw-h.size()-1)/2, '#') << endl;

        cout << "# Input" << endl;
        for (size_t i = 0; i < inputs.size(); i++) {
            cout << "         arg" << i << " = " << inputs[i] << endl;
        }

        cout << "# Output" << endl;
        for (size_t i = 0; i < outputs_scalar.size(); i++) {
            if (i < nb_lref)
                cout << "  scalar_arg" << i << " = " << outputs_scalar[i]
                     << endl
                     << "    simd_arg" << i << " = " << outputs_simd[i] << endl;
            else
                cout << "scalar_output = " << outputs_scalar[i] << endl
                     << "  simd_output = " << outputs_simd[i] << endl;
            cout << string (w2, ' ');
            for (unsigned j = 0; j < vect_size; j++) {
                bool b = eq (outputs_scalar[i][j], outputs_simd[i][j]);
                cout << (j ? ", " : "[ ") << string (w, b ? '=' : 'X');
            }
            cout << " ]" << endl;
        }
        cout << string (fw, '#') << endl << endl;

        cout.flags(old_settings);
    }

protected:
    size_t nb_lref;
    string name;
    vector<vectElt> inputs;
    vector<vectElt> outputs_simd;
    vector<vectElt> outputs_scalar;
};

/******************************************************************************/
/* Scalar functions for comparisons *******************************************/
/******************************************************************************/

template <typename Element, typename Enable = void>
struct ScalFunctionsBase;

/* for floating point element */
template <typename Element>
struct ScalFunctionsBase<Element,
                    typename enable_if<is_floating_point<Element>::value>::type>
{
    static constexpr Element _zero = 0.0;
    /* For cmp, when true, we return a NaN and assumes the value returned by the
     * Simd methods (0xFFFF...FFFF) is also a NaN.
     */
    static constexpr Element cmp_true = NAN;
    static constexpr Element cmp_false = _zero;

    class FloatingPointTestDistribution
    {
        public:
            using IntType = typename make_unsigned_int<Element>::type;

            FloatingPointTestDistribution () : intdist(std::numeric_limits<IntType>::lowest(), std::numeric_limits<IntType>::max()) {
            }

            template< class Generator >
            Element operator()( Generator& g )
            {
                IntType tmp = intdist (g);
                Element *fp_ptr = reinterpret_cast<Element *>(&tmp);
                return *fp_ptr;
            }

        private:
            uniform_int_distribution<IntType> intdist;

    };

    static FloatingPointTestDistribution get_default_random_generator () {
        return FloatingPointTestDistribution();
    }

    static Element ceil (Element x) {
        return std::ceil(x);
    }
    static Element floor (Element x) {
        return std::floor(x);
    }
    static Element round (Element x) {
        /* SSE and AVX round to nearest even integer value. The round function
         * from standard C++ library round to nearest with rounding up for half
         * integer. So we need to do a bit more work on the case of half integer
         * to completely emulate the behaviour of SSE and AVX.
         */
        Element r = std::round(x);
        if (std::abs (x - r) == 0.5)
        {
            if (std::fmod (r, 2.) == 1.)
                r -= 1.;
            else if (std::fmod (r, 2.) == -1.)
                r += 1.;
        }
        return r;
    }
    static Element blendv (Element a, Element b, Element mask) {
        using IntType = typename make_unsigned_int<Element>::type;
        IntType *ptr = (IntType *) &mask;
        *ptr = *ptr >> (8*sizeof(Element)-1);
        if (*ptr & 0x1)
            return b;
        else
            return a;
    }

    static Element fma (Element x, Element y, Element z) {
        return std::fma(x,y,z);
    }
};

/* for integral element */
template <typename Element>
struct ScalFunctionsBase<Element,
                    typename enable_if<is_integral<Element>::value>::type>
{
    static constexpr Element _zero = 0;
    /* For cmp, when true, we return -1 and assumes it is equal to the value
     * returned by the Simd methods (0xFFFF...FFFF).
     */
    static constexpr Element cmp_true = -1;
    static constexpr Element cmp_false = _zero;

    static uniform_int_distribution<Element> get_default_random_generator () {
        Element m = std::numeric_limits<Element>::lowest();
        Element M = std::numeric_limits<Element>::max();
        return uniform_int_distribution<Element>(m, M);
    }

    static Element round (Element x) {
        return x;
    }

    static Element fma (Element x, Element y, Element z) {
        return x*y + z;
    }

    static Element mullo (Element x1, Element x2) {
        return x1*x2;
    }
    static Element mulhi (Element x1, Element x2) {
        Integer q,r;
        Integer a = (Integer(x1)*Integer(x2));
        Integer b = Integer(1) << uint64_t(sizeof(Element)*8);
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
    static Element fmaddx (Element x1, Element x2, Element x3) {
        return x1 + mulx (x2, x3);
    }
    static Element fmaddxin (Element &x1, Element x2, Element x3) {
        return x1 += mulx (x2, x3);
    }
    static Element fmsubx (Element x1, Element x2, Element x3) {
        return -x1 + mulx (x2, x3);
    }
    static Element fmsubxin (Element &x1, Element x2, Element x3) {
        return x1 = -x1 + mulx (x2, x3);
    }
    static Element fnmaddx (Element x1, Element x2, Element x3) {
        return x1 - mulx(x2, x3);
    }
    static Element fnmaddxin (Element &x1, Element x2, Element x3) {
        return x1 -= mulx(x2, x3);
    }

    /* Shift */
    template <int s>
    static Element sra (Element x1) {
        if (std::is_signed<Element>::value) {
            /* For signed type we need to do a sign extension, the code comes
             * from http://graphics.stanford.edu/~seander/bithacks.html#FixedSignExtend
             */
            struct {Element x:sizeof(Element)*8-s;} r;
            return r.x = (x1 >> s);
        }
        else {
            return x1 >> s; /* For unsigned type, simply use >> */
        }
    }

    template <int s>
    static Element srl (Element x1) {
        return ((typename std::make_unsigned<Element>::type) x1) >> s;
    }

    template <int s>
    static Element sll (Element x1) {
        return ((typename std::make_unsigned<Element>::type) x1) << s;
    }
};

/* Struct with similar methods as the Simd structs */
template <typename Element>
struct ScalFunctions : public ScalFunctionsBase<Element>
{
    using vectElt = vector<Element>;
    using ScalFunctionsBase<Element>::cmp_true;
    using ScalFunctionsBase<Element>::cmp_false;
    using ScalFunctionsBase<Element>::fma;
    using ScalFunctionsBase<Element>::get_default_random_generator;

    static void genInputs (vector<vectElt> &inputs) {
        auto G = get_default_random_generator();
        auto rand = [&](){return G(entropy_generator);};
        for (auto &iv: inputs)
            std::generate (iv.begin(), iv.end(), rand);
    }

    static void genInputsWithZero (vector<vectElt> &inputs) {
        genInputs (inputs);
        for (auto &v: inputs[0])
            v = 0;
    }

    static Element zero () {
        return ScalFunctionsBase<Element>::_zero;
    }
    static Element vand (Element x1, Element x2) {
        unsigned char *p1 = reinterpret_cast<unsigned char *>(&x1);
        unsigned char *p2 = reinterpret_cast<unsigned char *>(&x2);
        for (unsigned int i = 0; i < sizeof (Element); i++)
            p1[i] &= p2[i];
        return x1;
    }
    static Element vor (Element x1, Element x2) {
        unsigned char *p1 = reinterpret_cast<unsigned char *>(&x1);
        unsigned char *p2 = reinterpret_cast<unsigned char *>(&x2);
        for (unsigned int i = 0; i < sizeof (Element); i++)
            p1[i] |= p2[i];
        return x1;
    }
    static Element vxor (Element x1, Element x2) {
        unsigned char *p1 = reinterpret_cast<unsigned char *>(&x1);
        unsigned char *p2 = reinterpret_cast<unsigned char *>(&x2);
        for (unsigned int i = 0; i < sizeof (Element); i++)
            p1[i] ^= p2[i];
        return x1;
    }
    static Element vandnot (Element x1, Element x2) {
        unsigned char *p1 = reinterpret_cast<unsigned char *>(&x1);
        unsigned char *p2 = reinterpret_cast<unsigned char *>(&x2);
        for (unsigned int i = 0; i < sizeof (Element); i++)
            p1[i] = (~p1[i]) & p2[i];
        return x1;
    }
    static Element add (Element x1, Element x2) {
        return x1+x2;
    }
    static Element addin (Element &x1, Element x2) {
        return x1+=x2;
    }
    static Element sub (Element x1, Element x2) {
        return x1-x2;
    }
    static Element subin (Element &x1, Element x2) {
        return x1-=x2;
    }
    static Element mul (Element x1, Element x2) {
        return x1*x2;
    }
    static Element mulin (Element &x1, Element x2) {
        return x1*=x2;
    }
    static Element div (Element x1, Element x2) {
        return x1/x2;
    }
    static Element fmadd (Element x1, Element x2, Element x3) {
        return fma(x3,x2,x1);
    }
    static Element fmaddin (Element &x1, Element x2, Element x3) {
        return x1 = fma(x3,x2,x1);
    }
    static Element fmsub (Element x1, Element x2, Element x3) {
        return fma(x3,x2,-x1);
    }
    static Element fmsubin (Element &x1, Element x2, Element x3) {
        return x1 = fma(x3,x2,-x1);
    }
    static Element fnmadd (Element x1, Element x2, Element x3) {
        return fma(-x3,x2,x1);
    }
    static Element fnmaddin (Element &x1, Element x2, Element x3) {
        return x1 = fma(-x3,x2,x1);
    }

    static Element lesser (Element x1, Element x2) {
        return (x1<x2)? cmp_true : cmp_false;
    }
    static Element lesser_eq (Element x1, Element x2) {
        return (x1<=x2)? cmp_true : cmp_false;
    }
    static Element greater (Element x1, Element x2) {
        return (x1>x2)? cmp_true : cmp_false;
    }
    static Element greater_eq (Element x1, Element x2) {
        return (x1>=x2)? cmp_true : cmp_false;
    }
    static Element eq (Element x1, Element x2) {
        return (x1==x2)? cmp_true : cmp_false;
    }

    static vectElt unpacklo (vectElt a, vectElt b) {
        vectElt r(a.size());
        for (size_t i = 0; i < a.size(); i++)
            r[i] = (i % 2 ? b[i/2] : a[i/2]);
        return r;
    }

    static vectElt unpackhi (vectElt a, vectElt b) {
        vectElt r(a.size());
        for (size_t i = 0; i < a.size(); i++)
            r[i] = (i % 2 ? b[(a.size()+i)/2] : a[(a.size()+i)/2]);
        return r;
    }
    static void unpacklohi (vectElt &lo, vectElt &hi, vectElt a, vectElt b) {
        lo = unpacklo (a, b);
        hi = unpackhi (a, b);
    }
    static vectElt pack_even (vectElt a, vectElt b) {
        vectElt r(a.size());
        for (size_t i = 0; i < a.size(); i++)
            r[i] = (i < a.size()/2 ? a[2*i] : b[2*i-a.size()]);
        return r;
    }
    static vectElt pack_odd (vectElt a, vectElt b) {
        vectElt r(a.size());
        for (size_t i = 0; i < a.size(); i++)
            r[i] = (i < a.size()/2 ? a[1+2*i] : b[1+2*i-a.size()]);
        return r;
    }
    static void pack (vectElt &even, vectElt &odd, vectElt a, vectElt b) {
        even = pack_even (a, b);
        odd = pack_odd (a, b);
    }
    template <uint16_t s>
    static vectElt blend (vectElt a, vectElt b) {
        vectElt r(a.size());
        for (size_t i = 0; i < a.size(); i++)
            r[i] = ((s >> i) & 0x1) ? b[i] : a[i];
        return r;
    }
};

/******************************************************************************/
/* Test one SIMD implem *******************************************************/
/******************************************************************************/

#define _TEST_ONE(K, f1, f2, r, n) do {     \
        K T(f1, f2, r, n);                  \
        bool b = T.writeResultLine();       \
        if (b == false)                     \
            T.writeDebugData();             \
        btest &= b;                         \
    } while (0)

#define TEST_ONE_OP(f) _TEST_ONE(TestOneMethod<Simd>,               \
        function<decltype(Simd::f)>(Simd::f),                       \
        function<decltype(Scal::f)>(Scal::f),                       \
        function<decltype(Scal::genInputs)>(Scal::genInputs), #f)
#define TEST_ONE_OP_WZ(f) _TEST_ONE(TestOneMethod<Simd>,                     \
        function<decltype(Simd::f)>(Simd::f),                                \
        function<decltype(Scal::f)>(Scal::f),                                \
        function<decltype(Scal::genInputsWithZero)>(Scal::genInputsWithZero),\
        #f " test with zero")

/* for floating point element */
template<typename Simd, typename Element>
typename enable_if<is_floating_point<Element>::value, bool>::type
test_impl_base () {
    using Scal = ScalFunctions<Element>;
    bool btest = true;

    TEST_ONE_OP (ceil);
    TEST_ONE_OP (floor);
    TEST_ONE_OP (mulin);
    TEST_ONE_OP (div);
    TEST_ONE_OP (blendv);

    return btest;
}

/* for integral element */
template<typename Simd, typename Element>
typename enable_if<is_integral<Element>::value, bool>::type
test_impl_base () {
    using Scal = ScalFunctions<Element>;
    bool btest = true;

    TEST_ONE_OP (mullo);
    TEST_ONE_OP (mulhi);
    TEST_ONE_OP (mulx);
    TEST_ONE_OP (fmaddx);
    TEST_ONE_OP (fmaddxin);
    TEST_ONE_OP (fmsubx);
    TEST_ONE_OP (fmsubxin);
    TEST_ONE_OP (fnmaddx);
    TEST_ONE_OP (fnmaddxin);
    TEST_ONE_OP (template sra<3>);
    TEST_ONE_OP (template sra<7>);
    TEST_ONE_OP (template srl<5>);
    TEST_ONE_OP (template srl<11>);
    TEST_ONE_OP (template sll<2>);
    TEST_ONE_OP (template sll<13>);

    return btest;
}

template<typename Simd, typename Element>
bool
test_impl () {
    using Scal = ScalFunctions<Element>;
    bool btest = test_impl_base<Simd, Element>();
    constexpr uint16_t blendmask = (0x1ul << Simd::vect_size) - 1;

    TEST_ONE_OP (zero);
    TEST_ONE_OP (vand);
    TEST_ONE_OP (vor);
    TEST_ONE_OP (vxor);
    TEST_ONE_OP (vandnot);
    TEST_ONE_OP (round);
    TEST_ONE_OP (add);
    TEST_ONE_OP (addin);
    TEST_ONE_OP (sub);
    TEST_ONE_OP (subin);
    TEST_ONE_OP (mul);
    TEST_ONE_OP (fmadd);
    TEST_ONE_OP (fmaddin);
    TEST_ONE_OP (fmsub);
    TEST_ONE_OP (fmsubin);
    TEST_ONE_OP (fnmadd);
    TEST_ONE_OP (fnmaddin);
    TEST_ONE_OP (lesser);
    TEST_ONE_OP (lesser_eq);
    TEST_ONE_OP (greater);
    TEST_ONE_OP (greater_eq);
    TEST_ONE_OP (eq);
    TEST_ONE_OP_WZ (lesser);
    TEST_ONE_OP_WZ (lesser_eq);
    TEST_ONE_OP_WZ (greater);
    TEST_ONE_OP_WZ (greater_eq);
    TEST_ONE_OP_WZ (eq);
    TEST_ONE_OP (unpacklo);
    TEST_ONE_OP (unpackhi);
    TEST_ONE_OP (unpacklohi);
    TEST_ONE_OP (pack_even);
    TEST_ONE_OP (pack_odd);
    TEST_ONE_OP (pack);
    TEST_ONE_OP (template blend<0x5555 & blendmask>);
    TEST_ONE_OP (template blend<0x9999 & blendmask>);
    TEST_ONE_OP (template blend<0xd7d7 & blendmask>);
    TEST_ONE_OP (template blend<0xdead & blendmask>);

    return btest;
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

#define TEST_IMPL(SIZE, Elt) do {       \
        pass &= test_impl<Simd##SIZE<Elt>, Elt>(); \
        cout << endl; \
    } while (0)

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
    TEST_IMPL(128, float);
    TEST_IMPL(128, double);
    TEST_IMPL(128, int16_t);
    TEST_IMPL(128, uint16_t);
    TEST_IMPL(128, int32_t);
    TEST_IMPL(128, uint32_t);
#ifdef __x86_64__
    TEST_IMPL(128, int64_t);
    TEST_IMPL(128, uint64_t);
#endif
#endif

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
    TEST_IMPL(256, float);
    TEST_IMPL(256, double);
#endif

#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
    TEST_IMPL(256, int16_t);
    TEST_IMPL(256, uint16_t);
    TEST_IMPL(256, int32_t);
    TEST_IMPL(256, uint32_t);
#ifdef __x86_64__
    TEST_IMPL(256, int64_t);
    TEST_IMPL(256, uint64_t);
#endif
#endif

#ifdef __FFLASFFPACK_HAVE_AVX512DQ_INSTRUCTIONS
    TEST_IMPL(512, float);
    TEST_IMPL(512, double);
    /* Simd512<(u)int16> does not exist. */
    /* Simd512<(u)int32> exists but was commented out in the previous version of
     * this test. Why ?
     */
    //TEST_IMPL(512, int32_t);
    //TEST_IMPL(512, uint32_t);
#ifdef __x86_64__
    TEST_IMPL(512, int64_t);
    TEST_IMPL(512, uint64_t);
#endif
#endif
    cout << endl << "Test " << (pass ? "passed" : "failed") << endl;
    return pass ? 0 : 1;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
