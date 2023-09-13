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
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-double-extended.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "givaro/givinteger.h"
#include "givaro/modular-integer.h"
#include "givaro/zring.h"
#include <givaro/givrns.h>

//#define BENCH_FLINT
#ifdef BENCH_FLINT
#define __GMP_BITS_PER_MP_LIMB 64
extern "C" {
#include "flint/longlong.h"
#include "flint/long_extras.h" 
#include "flint/fmpz_mat.h"
#include "flint/fmpz.h"
#include "flint/flint.h"
}
#endif


#define LOOPS_TIME 1

template <typename Field>
void run_bench(size_t n, size_t primes_bits, size_t b, size_t iters, size_t seed, bool matmul=false){
	Givaro::Integer p;
	FFLAS::Timer chrono;
	double timeFFLASToRNS=0.,timeFFLASFromRNS=0.,timeFFLASPrecomp=0.;
	double timeFFLASToRNS_ext=0.,timeFFLASFromRNS_ext=0.,timeFFLASPrecomp_ext=0.;
#ifdef BENCH_FLINT
	double timeFlintToRNS=0., timeFlintToRNSnaive=0.,timeFlintFromRNS=0.,timeFlintFromRNSnaive=0.,timeFlintPrecomp=0.;
#endif
	size_t loop=0;
	size_t bits=b;
	int logb= log((double)n)/log(2.0);

	size_t fflas_primes_bits       = (primes_bits ? primes_bits : std::min(25,41-logb));
	size_t fflas_primes_bits_ext   = 2* fflas_primes_bits;
	// std::cout<<"FFLAS     primes_bitsize= "<<fflas_primes_bits<<std::endl;
	// std::cout<<"FFLAS EXT primes_bitsize= "<<fflas_primes_bits_ext<<std::endl;

#ifdef BENCH_FLINT
	const size_t flint_primes_bits = (primes_bits ? primes_bits :NMOD_MAT_OPTIMAL_MODULUS_BITS);
	//std::cout<<"FLINT prime bit is "<<flint_primes_bits<<std::endl;
#endif

	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();

	size_t rns_bitsize;
	Givaro::Integer rns_bound;
	for (;( loop<iters || std::chrono::duration<double>(end-start) < std::chrono::seconds(LOOPS_TIME)) ;loop++){

		Givaro::Integer::random_exact_2exp(p, bits);
		//nextprime(p,p); avoid p to b a prime has it is not needed here
		Field F(p);
		typename Field::RandIter Rand(F,0,(uint64_t)seed);
		typename Field::Element_ptr A;
		A= FFLAS::fflas_new(F,n,n);

		FFLAS::frand(F,Rand,n,n,A,n);

		if (matmul){
			rns_bound=p*p*(uint64_t)n;
			rns_bitsize=rns_bound.bitsize();
		}else {
			rns_bound=p;
			rns_bitsize=rns_bound.bitsize();
		}

		int logb= log(double(rns_bitsize))/log(2.); // log of the integer bitsize of entry in A

		fflas_primes_bits       = (primes_bits ? primes_bits : std::min(25,41-logb));
		fflas_primes_bits_ext   =  2* fflas_primes_bits;



#ifdef BENCH_FLINT
		// FLINT RNS CODE //
		{
			fmpz_mat_t AA;
			fmpz_mat_init (AA, n, n);
			flint_rand_t randstate;
			flint_randinit(randstate);
			fmpz_mat_randbits(AA, randstate, bits);

			fmpz_comb_t comb;
			fmpz_comb_temp_t comb_temp;
			nmod_mat_t * mod_A;
			mp_limb_t * primes;
			mp_limb_t * residues;    
			size_t num_primes = (rns_bitsize + flint_primes_bits - 1) / flint_primes_bits;

			/* FLINT RNS initialization */
			residues = (mp_limb_t *)flint_malloc(sizeof(mp_limb_t) * num_primes);			
			mod_A    = (nmod_mat_t *)flint_malloc(sizeof(nmod_mat_t) * num_primes);
			primes   = (mp_limb_t *)flint_malloc(sizeof(mp_limb_t) * num_primes);
			primes[0] = n_nextprime(UWORD(1) << flint_primes_bits, 0);
			nmod_mat_init(mod_A[0], AA->r, AA->c, primes[0]);
			for (size_t i = 1; i < num_primes; i++)
				{
					primes[i] = n_nextprime(primes[i-1], 0);
					nmod_mat_init(mod_A[i], AA->r, AA->c, primes[i]);
				}

			chrono.clear();chrono.start();
			fmpz_comb_init(comb, primes, num_primes);
			fmpz_comb_temp_init(comb_temp, comb);
			chrono.stop();
			timeFlintPrecomp+=chrono.usertime();
			chrono.clear();chrono.start();

			/* Calculate residues of AA */
			fmpz_mat_multi_mod_ui_precomp(mod_A, num_primes, AA, comb, comb_temp); 

			chrono.stop();
			timeFlintToRNS+=chrono.usertime();
			chrono.clear();chrono.start();


			/* Calculate residue without COMB approach (NO FAST RNS) */
			//fmpz_mat_multi_mod_ui(mod_A, num_primes, AA);
			// for (long i = 0; i < AA->r; i++)
			// 	{   for (long j = 0; j < AA->c; j++)
			// 			for (long l = 0; l < num_primes; l++)
			// 				nmod_mat_entry(mod_A[l],i,j) =  fmpz_fdiv_ui(fmpz_mat_entry(AA,i,j), primes[l]);
			// 		}

			chrono.stop();
			timeFlintToRNSnaive+=chrono.usertime();
			chrono.clear();chrono.start();
			/* Chinese remaindering */
			for (long i = 0; i < AA->r; i++)
				{   for (long j = 0; j < AA->c; j++)
						{   for (size_t k = 0; k < num_primes; k++)
								residues[k] = mod_A[k]->rows[i][j];
							fmpz_multi_CRT_ui(&AA->rows[i][j], residues, comb, comb_temp, 1);
						}
				}
			chrono.stop();
			timeFlintFromRNS+=chrono.usertime();

			fmpz_mat_clear(AA);
			flint_free(mod_A);
			fmpz_comb_temp_clear(comb_temp);
			fmpz_comb_clear(comb);
			flint_free(residues);
			flint_free(primes);

		}
		//END FLINT CODE //
#endif

		// FFLAS RNS DOUBLE CODE
		if (fflas_primes_bits<=41-logb){
			chrono.clear();	chrono.start();
			// construct an RNS structure and its associated Domain
			FFPACK::rns_double RNS(rns_bound, fflas_primes_bits,false,seed);
			typedef FFPACK::RNSInteger<FFPACK::rns_double> RnsDomain;
			RnsDomain Zrns(RNS);
			typename RnsDomain::Element_ptr mod_A = FFLAS::fflas_new(Zrns,n,n);
			chrono.stop();
			timeFFLASPrecomp+=chrono.usertime();
			chrono.clear();	chrono.start();
			RNS.init(n,n,mod_A._ptr,mod_A._stride,A,n,p);
			chrono.stop();
			timeFFLASToRNS+=chrono.usertime();
			chrono.clear();	chrono.start();
			RNS.convert(n,n,Givaro::Integer(0), A,n, mod_A._ptr,mod_A._stride);
			//FFLAS::fconvert_rns(Zrns,n,n,Givaro::Integer(0),A,n,mod_A);
			chrono.stop();
			timeFFLASFromRNS+=chrono.usertime();
			FFLAS::fflas_delete(mod_A);
		}

		// FFLAS RNS DOUBLE EXTENDED CODE
		if (fflas_primes_bits_ext<=2*(41-logb)){
			chrono.clear();	chrono.start();	
			// construct an RNS structure and its associated Domain
			typedef FFPACK::rns_double_extended RNSExt;
			RNSExt RNS(rns_bound, fflas_primes_bits_ext,false,seed);
			typedef FFPACK::RNSInteger<RNSExt> RnsDomain;
			RnsDomain Zrns(RNS);
			typename RnsDomain::Element_ptr mod_A = FFLAS::fflas_new(Zrns,n,n);
			chrono.stop();
			timeFFLASPrecomp_ext+=chrono.usertime();
			chrono.clear();	chrono.start();

			RNS.init(n,n,mod_A._ptr,mod_A._stride,A,n,p);
			chrono.stop();
			timeFFLASToRNS_ext+=chrono.usertime();
			chrono.clear();	chrono.start();	
			//FFLAS::fconvert_rns(Zrns,n,n,Givaro::Integer(0),A,n,mod_A);
			RNS.convert(n,n,Givaro::Integer(0), A,n, mod_A._ptr,mod_A._stride);
			chrono.stop();
			timeFFLASFromRNS_ext+=chrono.usertime();
			FFLAS::fflas_delete(mod_A);
		}
		FFLAS::fflas_delete(A);
		end = std::chrono::system_clock::now();
	}

#define SPC1 10
#define SPC2 12
#define SPC3 33
#define PREC 2

	cout<<setw(SPC1+5)
	    <<loop<<" |"
	    <<setw(SPC1)
	    <<n<<" |"
	    <<setw(SPC1)
	    <<bits<<" |"
	    <<setw(SPC1+5)
	    <<fflas_primes_bits<<" |";

	//#define TIMING(x,y,z) cout<<setw(SPC2)<<std::fixed << std::setprecision(PREC)<<(x+y+z)/loop<<setw(SPC3)<<"("+to_string(x/loop)+","+to_string(y/loop)+","+to_string(z/loop)+")"<<" |";

#define TIMING(x,y,z) cout<<std::fixed << std::setprecision(PREC) \
			  <<setw(SPC3/3)<<((x/loop)*1000000.)\
			  <<setw(SPC3/3)<<((y/loop)*(1000000./(n*n)))\
			  <<setw(SPC3/3)<<((z/loop)*(1000000./(n*n)))<<" |";



	TIMING(timeFFLASPrecomp,timeFFLASFromRNS,timeFFLASToRNS);
	TIMING(timeFFLASPrecomp_ext,timeFFLASFromRNS_ext,timeFFLASToRNS_ext);
#ifdef BECNH_FLINT
	TIMING(timeFlintPrecomp,timeFlintFromRNS,timeFlintToRNS);
#endif
	cout<<endl;
}



int main(int argc, char** argv){
	static size_t iters = 10 ;
	static unsigned long b = 512 ;
	static unsigned long p = 0 ;
	static size_t m = 128 ;
	static bool longbench=false;
	static bool matmul=false;
	static size_t seed= time(NULL);	
	static Argument as[] = {
		{ 'b', "-b B", "Set the bitsize of the matrix entries.",         TYPE_INT , &b },
		{ 'p', "-p B", "Set the bitsize of the RNS Prime.",              TYPE_INT , &p },
		{ 'm', "-m M", "Set the dimension m of the matrix.",             TYPE_INT , &m },
		{ 'i', "-i R", "Set minimal number of repetitions.",             TYPE_INT , &iters },
		{ 'l', "-l Y", "Running a long benchmark .",                     TYPE_BOOL, &longbench },
		{ 't', "-t Y", "Running a matmul size benchmark .",              TYPE_BOOL, &matmul },
		{ 's', "-s S", "Set the seeding value",                          TYPE_INT , &seed },
		END_OF_ARGUMENTS
	};
	FFLAS::parseArguments(argc,argv,as);




	typedef Givaro::Modular<Givaro::Integer> Field;

	cout<<"### running RNS conversions benchmark ###"<<endl;
	if (matmul)
		cout<<"### matrix multiplication parameters are applied: rns bitsize is  (log n + 2 bitsize)"<<endl; 
	cout<<"### timing in microseconds (precomp,from, to) [from/to are given per element]"<<endl;
	cout<<"### |"
	    <<setw(SPC1)<<"loop"<<" |"
	    <<setw(SPC1)<<"matrix dim"<<" |"
	    <<setw(SPC1)<<"bitsize"<<" |"
	    <<setw(SPC1+5)<<"log2(prime)"<<" |"	
	    <<setw(SPC3)<<"FFLAS (precomp, from, to )"<<" |"
	    <<setw(SPC3)<<"FFLAS_EXT TO (precomp, from, to)"<<" |"
#ifdef BECNH_FLINT
	    <<setw(SPC3)<<"FLINT "<<" |"
#endif
	    <<endl;

	if (!longbench)
		run_bench<Field>(m,p,b,iters,seed,matmul);
	else {
		size_t matdim=m;	
		for (size_t bits=128; bits<(1<<20);bits<<=1)
			run_bench<Field>(matdim, p, bits, iters,seed++,matmul);

	}
	return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

