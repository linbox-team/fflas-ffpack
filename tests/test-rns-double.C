/*
 * Copyright (C) 2023 the FFLAS-FFPACK group
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *
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

//#define CHECK_RNS
//#define RNS_DEBUG
//#define BENCH_RNS
 
#define __FFLASFFPACK_SEQUENTIAL


#include "fflas-ffpack/fflas-ffpack-config.h"


#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
using namespace std; 

#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/utils/fflas_io.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-double-extended.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "givaro/givinteger.h" 
#include "givaro/modular-integer.h"
#include "givaro/zring.h"
#include <givaro/givrns.h>

template<typename RNS_t>
void print_rns();

template<>
void print_rns<FFPACK::rns_double>(){ std::cout<<std::setw(22)<<"(rns_double)";}
template<>
void print_rns<FFPACK::rns_double_extended>(){ std::cout<<std::setw(22)<<"(rns_double_extended)";}



template <typename RNS_t>
bool run_check(size_t n, size_t primes_bits, size_t b, size_t iters, size_t seed){
	Givaro::Integer p;

	size_t bits=b;
  bool ok=true;

  std::cout<<"Checking "; print_rns<RNS_t>();
  std::cout<<" int bits="<<std::left<<std::setw(10)<<bits<<" mi bits="<<std::left<<std::setw(10)<<primes_bits<<std::endl;
  Givaro::ZRing<Givaro::Integer> ZZ;
  Givaro::ZRing<double> ZD;

#define SPC  20

	for (size_t loop=0;(loop<iters);loop++){
    using Field =  Givaro::Modular<Givaro::Integer>;
    Givaro::Integer::seeding(seed);
		Givaro::Integer::random_exact_2exp(p, bits);

    Field F(p);
		Field::RandIter Rand(F,seed);
		Field::Element_ptr A,At,B;
		A     = FFLAS::fflas_new(F,n,n);
    At    = FFLAS::fflas_new(F,n,n);
    B     = FFLAS::fflas_new(F,n,n);

		FFLAS::frand(F,Rand,n,n,A,n);
    Givaro::Integer shift= (p-1)/2;
    for(size_t i=0;i<n*n;i++)
      A[i]-=shift;


    FFLAS::ftranspose(F,n,n,A,n,At,n);
    FFLAS::fassign(F,n,n,A,n,B,n);



    /* FFLAS::WriteMatrix(std::cout,F,n,n,A, n); */
    /* FFLAS::WriteMatrix(std::cout,F,n,n,At, n); */

    Givaro::Integer rns_bound = 3*p;

    // construct an RNS structure and its associated Domain
    RNS_t RNS(rns_bound, primes_bits,false,seed);
    typedef FFPACK::RNSInteger<RNS_t> RnsDomain;
    RnsDomain Zrns(RNS);
    typename RnsDomain::Element_ptr mod_A  = FFLAS::fflas_new(Zrns,n,n);
    typename RnsDomain::Element_ptr mod_At = FFLAS::fflas_new(Zrns,n,n);
    size_t n2=n*n, s=RNS._size;

    bool RNSMAJOR;
    /**********************
     ** check RNSMajor=false
     **********************/
#if 1
    RNSMAJOR=false;
    
    RNS.init(n,n,mod_A._ptr,mod_A._stride,A,n,p, RNSMAJOR);
    ok = ok and  FFPACK::RNS_COMMON::check_integers_to_rns(FFLAS::FflasNoTrans, n, n, mod_A._ptr, mod_A._stride, A, n, RNS._basis, RNSMAJOR, false);

    RNS.init_transpose(n,n,mod_At._ptr,mod_At._stride,At,n,p,RNSMAJOR);
    ok = ok and  FFPACK::RNS_COMMON::check_integers_to_rns(FFLAS::FflasTrans, n, n, mod_At._ptr, mod_At._stride, At, n, RNS._basis, RNSMAJOR, false);

    ok = ok and (FFLAS::fequal(ZD,s,n2,mod_A._ptr,n2,mod_At._ptr,n2));

      /* std::cout<<std::boolalpha<<ok<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZD,n2, s, mod_A._ptr, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZD,n2, s, mod_At._ptr, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */

    RNS.convert(n,n,Givaro::Integer(0), B, n, mod_A._ptr,mod_A._stride,RNSMAJOR);
    ok= ok and (FFLAS::fequal(ZZ,n,n,A,n,B,n));

      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, A, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, B, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */

    RNS.convert_transpose(n,n,Givaro::Integer(0), B, n, mod_A._ptr,mod_A._stride,RNSMAJOR);
    ok= ok and (FFLAS::fequal(ZZ,n,n,At,n,B,n));

      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, At, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, B, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */

      std::cout<<std::right<<std::setw(SPC)<<"RNSMajor=false ..."  <<  (ok? "PASSED\n": "FAILED\n");

#endif

#if 1
      /**********************
       ** check RNSMajor=true
       **********************/
      RNSMAJOR=true;

      RNS.init(n,n,mod_A._ptr,mod_A._stride,A,n,p,RNSMAJOR);
      ok = ok and  FFPACK::RNS_COMMON::check_integers_to_rns(FFLAS::FflasNoTrans, n, n, mod_A._ptr, mod_A._stride, A, n, RNS._basis, RNSMAJOR, false);

      RNS.init_transpose(n,n,mod_At._ptr,mod_At._stride,At,n,p,RNSMAJOR);
      ok = ok and  FFPACK::RNS_COMMON::check_integers_to_rns(FFLAS::FflasTrans, n, n, mod_At._ptr, mod_At._stride, At, n, RNS._basis, RNSMAJOR, false);

      ok = ok and (FFLAS::fequal(ZD,n2,s,mod_A._ptr,s,mod_At._ptr,s));

      /* std::cout<<std::boolalpha<<ok<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZD,n2, s, mod_A._ptr, s)<<std::endl<<"--------------\n";;//FFLAS::FflasSageMath)<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZD,n2, s, mod_At._ptr, s)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */


			RNS.convert(n,n,Givaro::Integer(0), B, n, mod_A._ptr,mod_A._stride,RNSMAJOR);
      ok= ok and (FFLAS::fequal(ZZ,n,n,A,n,B,n));

      /* std::cout<<std::boolalpha<<ok<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, A, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, B, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */ 
      

      RNS.convert_transpose(n,n,Givaro::Integer(0), B, n, mod_A._ptr,mod_A._stride,RNSMAJOR);
      ok= ok and (FFLAS::fequal(ZZ,n,n,At,n,B,n));

      /* std::cout<<std::boolalpha<<ok<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, At, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */
      /* FFLAS::WriteMatrix(std::cout,ZZ,n, n, B, n)<<std::endl;//FFLAS::FflasSageMath)<<std::endl; */


      std::cout<<std::right<<std::setw(SPC)<<"RNSMajor=true  ..."  <<  (ok? "PASSED\n": "FAILED\n");


#endif

			FFLAS::fflas_delete(mod_A);
      FFLAS::fflas_delete(mod_At);
      FFLAS::fflas_delete(A);
      FFLAS::fflas_delete(At);
      FFLAS::fflas_delete(B);
		}

  std::cout<<"Test suite for rns-double ...."<<(ok?"PASSED\n":"FAILED\n");
  if (!ok) std::cout<<"seed was "<<seed<<std::endl;
  return ok;
}



int main(int argc, char** argv){
	static size_t iters = 1 ;
	static unsigned long b = 512 ;
	static unsigned long p = 20 ;
	static size_t m = 16 ;
	static size_t seed= time(NULL);	
	static Argument as[] = {
		{ 'b', "-b B", "Set the bitsize of the matrix entries.",         TYPE_INT , &b },
		{ 'p', "-p B", "Set the bitsize of the RNS Prime.",              TYPE_INT , &p },
		{ 'm', "-m M", "Set the dimension m of the matrix.",             TYPE_INT , &m },
		{ 'i', "-i R", "Set minimal number of repetitions.",             TYPE_INT , &iters },
		{ 's', "-s S", "Set the seeding value",                          TYPE_INT , &seed },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);

	/* int logb= log((double)n)/log(2.0); */
	/* size_t fflas_primes_bits       = (primes_bits ? primes_bits : 41-logb); */
	/* size_t fflas_primes_bits_ext   = 2* fflas_primes_bits; */

  bool ok=true;
  ok = ok and run_check<FFPACK::rns_double>         (m,  p,b,iters,seed);
  ok = ok and run_check<FFPACK::rns_double_extended>(m,2*p,b,iters,seed);

	return ok;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

