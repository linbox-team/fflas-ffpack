/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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

 #include "fflas-ffpack/utils/simd.h"
 #include "fflas-ffpack/utils/args-parser.h"
 #include "fflas-ffpack/utils/align-allocator.h"
 #include <vector>
 #include <algorithm>
 #include <random>

template<class FuncT, class Element>
inline void test_binary(FuncT && func, size_t n, const Element * a, const Element * b, const Element * c)
{
	using simd = Simd<Element>;
	typename simd::vect_t v1, v2;
	for(size_t i = 0 ; i < n ; i+=simd::vect_size){
		v1 = simd::load(a+i);
		v2 = simd::load(b+i);
		simd::store(c+i, func(v1, v2));
	}
}

template<class FuncT, class Element>
inline void test_unary(FuncT && func, size_t n, const Element * a)
{
	using simd = Simd<Element>;
	typename simd::vect_t v1;
	for(size_t i = 0 ; i < n ; i+=simd::vect_size){
		v1 = simd::load(a+i);
		simd::store(a+i, func(v1));
	}
}

template<class FuncT, class Element>
inline void test_znary(FuncT && func, size_t n, const Element * a)
{
	using simd = Simd<Element>;
	for(size_t i = 0 ; i < n ; i+=simd::vect_size){
		simd::store(a+i, func());
	}
}

template<class Element>
bool test_float(size_t seed)
{
 	using simd = Simd<Element>;
 	
 	std::mt19937 generator(seed);
 	std::uniform_real_distribution<> dist(1, 100);

 	size_t vectorSize = 32;
 	std::vector<Element, AlignedAllocator<Element, Alignment::AVX>> a1(vectorSize), b1(vectorSize), c1(vectorSize), a2(vectorSize), b2(vectorSize), c2(vectorSize);

 	test_znary(simd::zero, vectorSize, a2.data());
	std::fill(a1.begin(), a1.end(), 0);
	if(!std::equal(a1.begin(), a1.end(), a2.begin())){
		return false;
	}

	std::generate(a1.begin(), a1.end(), [&](){return dist(generator);});
	std::generate(b1.begin(), b1.end(), [&](){return dist(generator);});
	a2 = a1;
	b2 = b1;

	test_binary(simd::add, vectorSize, a2.data(), b2.data(), c2.data());
	for(size_t i = 0 ; i < vectorSize ; ++i)
		c1[i] = b1[i]+a1[i];
	if(!std::equal(c1.begin(), c1.end(), c2.begin())){
		return false;
	}

	std::fill(a1.begin(), a1.end(), 0);
	std::fill(a2.begin(), a2.end(), 0);
	std::fill(b1.begin(), b1.end(), 0);
	std::fill(b2.begin(), b2.end(), 0);	

	std::generate(a1.begin(), a1.end(), [&](){return dist(generator);});
	std::generate(b1.begin(), b1.end(), [&](){return dist(generator);});
	a2 = a1;
	b2 = b1;

	test_binary(simd::mul, vectorSize, a2.data(), b2.data(), c2.data());
	for(size_t i = 0 ; i < vectorSize ; ++i)
		c1[i] = b1[i]*a1[i];
	if(!std::equal(c1.begin(), c1.end(), c2.begin())){
		return false;
	}

	std::fill(a1.begin(), a1.end(), 0);
	std::fill(a2.begin(), a2.end(), 0);
	std::fill(b1.begin(), b1.end(), 0);
	std::fill(b2.begin(), b2.end(), 0);	

	std::generate(a1.begin(), a1.end(), [&](){return dist(generator);});
	std::generate(b1.begin(), b1.end(), [&](){return dist(generator);});
	a2 = a1;
	b2 = b1;

	test_binary(simd::sub, vectorSize, a2.data(), b2.data(), c2.data());
	for(size_t i = 0 ; i < vectorSize ; ++i)
		c1[i] = a1[i]-b1[i];
	if(!std::equal(c1.begin(), c1.end(), c2.begin())){
		return false;
	}

	std::fill(a1.begin(), a1.end(), 0);
	std::fill(a2.begin(), a2.end(), 0);
	std::fill(b1.begin(), b1.end(), 0);
	std::fill(b2.begin(), b2.end(), 0);	

	return true;
}


 template<class Element>
 bool test_integer()
 {
 	return true;
 }


 int main(int ac, char **av) {
	int seed = (int) time(NULL);

	static Argument as[] = {
		{ 's', "-s N", "Set the seed                 .", TYPE_INT , &seed },
		END_OF_ARGUMENTS
	};


	FFLAS::parseArguments(ac,av,as);

	srand(seed);
	srand48(seed);


	bool pass  = true ;
	{ /*  finit */
		{
			pass &= test_float<float>(seed);
		}
		{
			pass &= test_float<double>(seed);
		}
	}
	std::cout << std::boolalpha << pass << std::endl;
	return (pass?0:1) ;
}