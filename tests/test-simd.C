/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#include "givaro/givinteger.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/fflas/fflas_simd.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/utils/align-allocator.h"
#include <vector>
#include <algorithm>
#include <random>
#include <tuple>
#include <type_traits>
#include <string>
#include <iterator>
#include <limits>
#include <cmath>
#include <iomanip>

typedef Givaro::Integer integer;

/**********************************************************************************
 *
 * Random generators
 *
 ***********************************************************************************/

template <class Element, class Alloc>
typename std::enable_if<std::is_integral<Element>::value>::type
generate_random (std::vector<Element,Alloc> &a, std::mt19937 &generator) {
	std::uniform_int_distribution<Element> dist(std::numeric_limits<Element>::min(), std::numeric_limits<Element>::max());
	std::generate(a.begin(), a.end(), [&](){return dist(generator);});
}

template <class Element, class Alloc>
typename std::enable_if<std::is_floating_point<Element>::value>::type
generate_random (std::vector<Element,Alloc> &a, std::mt19937 &generator) {
	std::uniform_real_distribution<Element> dist(std::numeric_limits<Element>::min(), std::numeric_limits<Element>::max());
	std::generate(a.begin(), a.end(), [&](){return dist(generator);});
}

/**********************************************************************************
 *
 * Function Traits
 *
 ***********************************************************************************/

template <class F> struct function_traits;

// function pointer
template <class R, class... Args>
struct function_traits<R (*)(Args...)> : public function_traits<R(Args...)> {};

template <class R, class... Args> struct function_traits<R(Args...)> {
	using return_type = R;

	static constexpr std::size_t arity = sizeof...(Args);

	template <std::size_t N> struct argument {
		static_assert(N < arity, "error: invalid parameter index.");
		using type = typename std::tuple_element<N, std::tuple<Args...> >::type;
	};
};

// member function pointer
template <class C, class R, class... Args>
struct function_traits<R (C::*)(Args...)>       : public function_traits<R(C&, Args...)> {};

// const member function pointer
template <class C, class R, class... Args>
struct function_traits<R (C::*)(Args...) const> : public function_traits<R(C&, Args...)> {};

// member object pointer
template <class C, class R>
struct function_traits<R(C::*)>                 : public function_traits<R(C&)> {};

template<class SimdFunc>
void print_arity (SimdFunc f) {
	std::cout << "Arity of function is " << (function_traits<SimdFunc>::arity) << std::endl;
}

/**************************************************************************************/

template<class simd, class Element, class SimdFunc, class ScalFunc>
inline
typename std::enable_if<
(function_traits<SimdFunc>::arity == 0) &&
!(std::is_same<typename function_traits<SimdFunc>::return_type, void>::value)
, bool>::type
test_op(SimdFunc && fsimd, ScalFunc && fscal, uint64_t seed, size_t vectorSize, Element max, std::string name){

	using vect_t = typename simd::vect_t;

	std::vector<Element, AlignedAllocator<Element, Alignment::AVX>> c1(vectorSize), c2(vectorSize);

	std::transform(c1.begin(), c1.end(), c1.begin(), fscal);

	vect_t vc2;
	for(size_t i = 0 ; i < vectorSize ; i+=simd::vect_size){
			c2 = fsimd();
			simd::store(c2.data()+i, c2);
		}

	bool res = std::equal(c1.begin(), c1.end(), c2.begin(), [](Element x1, Element x2){return (std::isnan(x1) && std::isnan(x2)) || x1 == x2;});
	if(!res)
		{
			std::cout << "Error Simd" << sizeof(typename simd::scalar_t)*simd::vect_size*8 << "::" << name
					  << " on " << (sizeof(Element) * 8) << "bits." << std::endl;

			std::copy(c1.begin(), c1.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::copy(c2.begin(), c2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl ;
		}
	return res;
}

template<class simd, class Element, class SimdFunc, class ScalFunc>
inline
typename std::enable_if<
(function_traits<SimdFunc>::arity == 1) &&
!(std::is_same<typename function_traits<SimdFunc>::return_type, void>::value)
, bool>::type
test_op(SimdFunc fsimd, ScalFunc fscal, uint64_t seed, size_t vectorSize, Element max, std::string name){
	
	using vect_t = typename simd::vect_t;

	std::mt19937 generator(seed);
	std::vector<Element, AlignedAllocator<Element, Alignment::AVX>> a1(vectorSize), c1(vectorSize), a2(vectorSize), c2(vectorSize), c3(vectorSize);
	generate_random(a1, generator);
	a2 = a1;

	std::transform(a1.begin(), a1.end(), c1.begin(), fscal);

	vect_t va2, vc2, vc3;
	for(size_t i = 0 ; i < vectorSize ; i+=simd::vect_size){
			va2 = simd::load(a2.data()+i);
			vc3 = simd::load(c1.data()+i);
			vc2 = fsimd(va2);
			vc3 = simd::sub(vc3,vc2);
			simd::store(c2.data()+i, vc2);
			simd::store(c3.data()+i, vc3);
		}

	bool res = std::equal(c1.begin(), c1.end(), c2.begin(), [](Element x1, Element x2){return (std::isnan(x1) && std::isnan(x2)) || x1 == x2;});
	if(!res)
		{
			std::cout << "Error Simd" << sizeof(typename simd::scalar_t)*simd::vect_size*8 << "::" << name
					  << " on " << (sizeof(Element) * 8) << "bits." << std::endl;

			std::cout << "a2: ";
			std::copy(a2.begin(), a2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "c1: ";
			std::copy(c1.begin(), c1.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "c2: ";
			std::copy(c2.begin(), c2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl << std::endl;
			std::cout << "c1-c2: ";
			std::copy(c3.begin(), c3.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl << std::endl;
		}
	return res;
}

template<class simd, class Element, class SimdFunc, class ScalFunc>
inline
typename std::enable_if<
(function_traits<SimdFunc>::arity == 2) &&
!(std::is_same<typename function_traits<SimdFunc>::return_type, void>::value)
, bool>::type
test_op(SimdFunc fsimd, ScalFunc fscal, uint64_t seed, size_t vectorSize, Element max, std::string name){
	
	using vect_t = typename simd::vect_t;

	std::mt19937 generator(seed);
	std::vector<Element, AlignedAllocator<Element, Alignment::AVX>> a1(vectorSize), b1(vectorSize), c1(vectorSize), a2(vectorSize), b2(vectorSize), c2(vectorSize), c3(vectorSize);
	generate_random(a1, generator);
	generate_random(b1, generator);
	a2 = a1;
	b2 = b1;

	std::transform(a1.begin(), a1.end(), b1.begin(), c1.begin(), fscal);

	vect_t va2, vb2, vc2, vc3;
	for(size_t i = 0 ; i < vectorSize ; i+=simd::vect_size){
			va2 = simd::load(a2.data()+i);
			vb2 = simd::load(b2.data()+i);
			vc3 = simd::load(c1.data()+i);
			vc2 = fsimd(va2, vb2);
			vc3 = simd::sub(vc3,vc2);
			simd::store(c2.data()+i, vc2);
			simd::store(c3.data()+i, vc3);
		}

	bool res = std::equal(c1.begin(), c1.end(), c2.begin(), [](Element x1, Element x2){return (std::isnan(x1) && std::isnan(x2)) || x1 == x2;});
	if(!res)
		{
			std::cout << "Error Simd" << sizeof(typename simd::scalar_t)*simd::vect_size*8 << "::" << name
					  << " on " << (sizeof(Element) * 8) << "bits." << std::endl;

			std::cout << "a2: ";
			std::copy(a2.begin(), a2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "b2: ";
			std::copy(b2.begin(), b2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "c1: ";
			std::copy(c1.begin(), c1.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "c2: ";
			std::copy(c2.begin(), c2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl << std::endl;
			std::cout << "c1-c2: ";
			std::copy(c3.begin(), c3.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl << std::endl;
		}
	return res;
}

template<class simd, class Element, class SimdFunc, class ScalFunc>
inline
typename std::enable_if<
(function_traits<SimdFunc>::arity == 3) &&
!(std::is_same<typename function_traits<SimdFunc>::return_type, void>::value)
, bool>::type
test_op(SimdFunc fsimd, ScalFunc fscal, uint64_t seed, size_t vectorSize, Element max, std::string name){
	
	using vect_t = typename simd::vect_t;

	std::mt19937 generator(seed);
	std::vector<Element, AlignedAllocator<Element, Alignment::AVX>> a1(vectorSize), b1(vectorSize), c1(vectorSize), d1(vectorSize), a2(vectorSize), b2(vectorSize), c2(vectorSize), d2(vectorSize);
	generate_random(a1, generator);
	generate_random(b1, generator);
	generate_random(c1, generator);
	a2 = a1;
	b2 = b1;
	c2 = c1;

	for(size_t i = 0 ; i < vectorSize ; ++i){
			d1[i] = fscal(c1[i], a1[i], b1[i]);
		}

	vect_t va2, vb2, vc2;
	for(size_t i = 0 ; i < vectorSize ; i+=simd::vect_size){
			va2 = simd::load(a2.data()+i);
			vb2 = simd::load(b2.data()+i);
			vc2 = simd::load(c2.data()+i);
			simd::store(d2.data()+i, fsimd(vc2, va2, vb2));
		}

	bool res = std::equal(d1.begin(), d1.end(), d2.begin(), [](Element x1, Element x2){return (std::isnan(x1) && std::isnan(x2)) || x1 == x2;});
	if(!res)
		{
			std::cout << "Error Simd" << sizeof(typename simd::scalar_t)*simd::vect_size*8 << "::" << name
					  << " on " << (sizeof(Element) * 8) << "bits." << std::endl;

			std::transform(d1.begin(), d1.end(), d2.begin(), d2.begin(), [](Element x1, Element x2){return x1-x2;});

			//std::copy(d1.begin(), d1.end(), std::ostream_iterator<Element>(std::cout, " "));
			//std::cout << std::endl;
			std::copy(d2.begin(), d2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
		}
	return res;
}

template<class simd, class Element, class SimdFunc, class ScalFunc> 
inline
typename std::enable_if<true, bool>::type
test_blend(SimdFunc fsimd, ScalFunc fscal, uint64_t seed, size_t vectorSize, Element max, std::string name){
	using vect_t = typename simd::vect_t;
	uint8_t s = 0x1;
	std::mt19937 generator(seed);
	std::vector<Element, AlignedAllocator<Element, Alignment::AVX>> a1(vectorSize), b1(vectorSize), c1(vectorSize), a2(vectorSize), b2(vectorSize), c2(vectorSize), c3(vectorSize);
	generate_random(a1, generator);
	generate_random(b1, generator);
	a2 = a1;
	b2 = b1;

	std::transform(a1.begin(), a1.end(), b1.begin(), c1.begin(), fscal);

	vect_t va2, vb2, vc2, vc3;
	for(size_t i = 0 ; i < vectorSize ; i+=simd::vect_size){
			va2 = simd::load(a2.data()+i);
			vb2 = simd::load(b2.data()+i);
			vc3 = simd::load(c1.data()+i);
			vc2 = fsimd<s>(va2, vb2);
			vc3 = simd::sub(vc3,vc2);
			simd::store(c2.data()+i, vc2);
			simd::store(c3.data()+i, vc3);
		}

	bool res = std::equal(c1.begin(), c1.end(), c2.begin(), [](Element x1, Element x2){return (std::isnan(x1) && std::isnan(x2)) || x1 == x2;});
	if(!res)
		{
			std::cout << "Error Simd" << sizeof(typename simd::scalar_t)*simd::vect_size*8 << "::" << name
					  << " on " << (sizeof(Element) * 8) << "bits." << std::endl;

			std::cout << "a2: ";
			std::copy(a2.begin(), a2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "b2: ";
			std::copy(b2.begin(), b2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "c1: ";
			std::copy(c1.begin(), c1.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl;
			std::cout << "c2: ";
			std::copy(c2.begin(), c2.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl << std::endl;
			std::cout << "c1-c2: ";
			std::copy(c3.begin(), c3.end(), std::ostream_iterator<Element>(std::cout, " "));
			std::cout << std::endl << std::endl;
		}
	return res;
}

template<class simd, class Element>
bool test_float_impl(uint64_t seed, size_t vectorSize, Element max){
	bool btest = true;

	btest = btest && test_op<simd>(simd::ceil, [](Element x){return std::ceil(x);}, seed, vectorSize, max, "ceil");
	btest = btest && test_op<simd>(simd::floor, [](Element x){return std::floor(x);}, seed, vectorSize, max,"floor");
	btest = btest && test_op<simd>(simd::round, [](Element x){return std::round(x);}, seed, vectorSize, max, "round");
	btest = btest && test_op<simd>(simd::add, [](Element x1, Element x2){return x1+x2;}, seed, vectorSize, max, "add");
	btest = btest && test_op<simd>(simd::sub, [](Element x1, Element x2){return x1-x2;}, seed, vectorSize, max, "sub");
	btest = btest && test_op<simd>(simd::mul, [](Element x1, Element x2){return x1*x2;}, seed, vectorSize, max, "mul");
	btest = btest && test_op<simd>(simd::fmadd, [](Element x1, Element x2, Element x3){return std::fma(x3,x2,x1);}, seed, vectorSize, max, "fmadd");
	btest = btest && test_op<simd>(simd::fmsub, [](Element x1, Element x2, Element x3){return std::fma(x3,x2,-x1);}, seed, vectorSize, max, "fmsub");
	btest = btest && test_op<simd>(simd::fnmadd, [](Element x1, Element x2, Element x3){return std::fma(-x3,x2,x1);}, seed, vectorSize, max, "fnmadd");
	btest = btest && test_op<simd>(simd::lesser, [](Element x1, Element x2){return (x1<x2)?NAN:0;}, seed, vectorSize, max, "lesser");
	btest = btest && test_op<simd>(simd::lesser_eq, [](Element x1, Element x2){return (x1<=x2)?NAN:0;}, seed, vectorSize, max, "lesser_eq");
	btest = btest && test_op<simd>(simd::greater, [](Element x1, Element x2){return (x1>x2)?NAN:0;}, seed, vectorSize, max, "greater");
	btest = btest && test_op<simd>(simd::greater_eq, [](Element x1, Element x2){return (x1>=x2)?NAN:0;}, seed, vectorSize, max, "greater_eq");
	btest = btest && test_op<simd>(simd::eq, [](Element x1, Element x2){return (x1==x2)?NAN:0;}, seed, vectorSize, max, "eq");
	return btest;
}

template<typename simd>
typename simd::vect_t mysra (typename simd::vect_t x1){return simd::sra(x1, int(2));}

template<typename simd>
typename simd::vect_t myblend (typename simd::vect_t x1, typename simd::vect_t x2){return simd::blend(x1,x2);}

template<class simd, class Element>
bool test_integer_impl(uint64_t seed, size_t vectorSize, Element max){
	bool btest = true;

	btest = btest && test_op<simd>(simd::add, [](Element x1, Element x2){return x1+x2;}, seed, vectorSize, max, "add");
	btest = btest && test_op<simd>(simd::sub, [](Element x1, Element x2){return x1-x2;}, seed, vectorSize, max, "sub");
	btest = btest && test_op<simd>(simd::mullo, [](Element x1, Element x2){return x1*x2;}, seed, vectorSize, max, "mullo");
	btest = btest && test_op<simd>(simd::mul, [](Element x1, Element x2){return x1*x2;}, seed, vectorSize, max, "mullo");
	btest = btest && test_op<simd>(simd::fmadd, [](Element x1, Element x2, Element x3){return x1+x3*x2;}, seed, vectorSize, max, "fmadd");
	btest = btest && test_op<simd>(simd::fmsub, [](Element x1, Element x2, Element x3){return -x1+x3*x2;}, seed, vectorSize, max, "fmsub");
	btest = btest && test_op<simd>(simd::fnmadd, [](Element x1, Element x2, Element x3){return x1-x3*x2;}, seed, vectorSize, max, "fnmadd");
	btest = btest && test_op<simd>(simd::lesser, [](Element x1, Element x2){return (x1<x2)?-1:0;}, seed, vectorSize, max, "lesser");
	btest = btest && test_op<simd>(simd::lesser_eq, [](Element x1, Element x2){return (x1<=x2)?-1:0;}, seed, vectorSize, max, "lesser_eq");
	btest = btest && test_op<simd>(simd::greater, [](Element x1, Element x2){return (x1>x2)?-1:0;}, seed, vectorSize, max, "greater");
	btest = btest && test_op<simd>(simd::greater_eq, [](Element x1, Element x2){return (x1>=x2)?-1:0;}, seed, vectorSize, max, "greater_eq");
	btest = btest && test_op<simd>(simd::eq, [](Element x1, Element x2){return (x1==x2)?-1:0;}, seed, vectorSize, max, "eq");
	btest = btest && test_blend<simd>(myblend<simd>, [](Element x1, Element x2){

			return false;
}, seed, vectorSize, max, "blend");
	// print_arity(mysra<simd>);
	btest = btest && test_op<simd>(mysra<simd>, //std::bind(simd::sra,std::placeholders::_1,int(sizeof(Element)*4)),
						   [](Element x1){
			integer h = integer (1) << 2;
			integer r = integer(x1) / h;
			r -= ((integer(x1)-h*r) < 0)?1:0;
			return Element(r);
			// return Element(std::floor(double(x1)/double(h)));
}, seed, vectorSize, max, "sra");
	btest = btest && test_op<simd>(simd::mulhi, [](Element x1, Element x2){
			integer q,r;
			integer a = (integer(x1)*integer(x2));
			integer b = integer(1) << uint64_t(sizeof(Element)*8);
			Givaro::IntegerDom Z;
			Z.divmod(q, r, a, b);
			return Element(q);
}, seed, vectorSize, max, "mulhi");
	btest = btest && test_op<simd>(simd::mulx, [](Element x1, Element x2){
			Element h = Element(1) << (sizeof(Element)*4);
			/* Representative r of x1 modulo h with -h/2 <= r < h/2*/
			if (std::is_signed<Element>::value) {
			x1 = (x1+h/2)%h;
			x1 += (x1 < 0)?h/2:-h/2;
			x2 = (x2+h/2)%h;
			x2 += (x2 < 0)?h/2:-h/2; }
			else {
			x1 = x1 % h;
			x2 = x2 % h;
}
			return x1*x2;
}, seed, vectorSize, max, "mulx");

	return btest;
}

template<class Element>
bool test_float(uint64_t seed, size_t vectorSize, size_t max_){
	bool sse = true, avx = true;
#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
	sse = test_float_impl<Simd128<Element>>(seed, vectorSize, (Element)max_);
	if(!sse)
		std::cout << "bug sse" << std::endl;
	else
		std::cout << "SSE OK" << std::endl;
#endif

#ifdef __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
	avx = test_float_impl<Simd256<Element>>(seed, vectorSize, (Element)max_);
	if(!avx)
		std::cout << "bug avx" << std::endl;
	else
		std::cout << "AVX OK" << std::endl;
#endif
	return sse && avx;
}

template<class Element>
bool test_integer(uint64_t seed, size_t vectorSize, size_t max_){
	bool sse = true, avx = true;
#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
	sse = test_integer_impl<Simd128<Element>>(seed, vectorSize, (Element)max_);
	if(!sse)
		std::cout << "bug sse" << std::endl;
	else
		std::cout << "SSE OK" << std::endl;
#endif
#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
	avx = test_integer_impl<Simd256<Element>>(seed, vectorSize, (Element)max_);
	if(!avx)
		std::cout << "bug avx" << std::endl;
	else
		std::cout << "AVX OK" << std::endl;
#endif
	return sse && avx;
}


int main(int ac, char **av) {
	uint64_t seed = FFLAS::getSeed();
	int vectorSize = 32;
	int max = 100;
	int loop = false;

	static Argument as[] = {
		{ 's', "-s N", "Set the seed                 .", TYPE_UINT64 , &seed },
		{ 'l', "-l N", "Set the loop execution       .", TYPE_INT , &loop },
		END_OF_ARGUMENTS
	};

	FFLAS::parseArguments(ac,av,as);

	srand(seed);
	srand48(seed);

	bool pass  = true ;
	{
		do{
				{
					pass = pass && test_float<float>(seed, vectorSize, max);
				}
				{
					pass = pass && test_float<double>(seed, vectorSize, max);
				}
				{
					pass = pass && test_integer<int16_t>(seed, vectorSize, max);
				}
				{
					pass = pass && test_integer<int32_t>(seed, vectorSize, max);
				}
				{
					pass = pass && test_integer<int64_t>(seed, vectorSize, max);
				}
				{
					pass = pass && test_integer<uint16_t>(seed, vectorSize, max);
				}
				{
					pass = pass && test_integer<uint32_t>(seed, vectorSize, max);
				}
				{
					pass = pass && test_integer<uint64_t>(seed, vectorSize, max);
				}
			}while(loop);
	}
	std::cout << std::boolalpha << pass << std::endl;
	return (pass?0:1) ;
}
