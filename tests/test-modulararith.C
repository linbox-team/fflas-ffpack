// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>

#include <iostream>
#include <cassert>
#include "fflas-ffpack/utils/timer.h"
#include "fflas-ffpack/field/modular-float.h"
#include "fflas-ffpack/field/modular-double.h"
#include "fflas-ffpack/field/modular-balanced-int64.h"
#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/field/modular-balanced-double.h"
#include "fflas-ffpack/field/modular-balanced-float.h"
#include "fflas-ffpack/field/modular-balanced.h"
#include "fflas-ffpack/field/modular-int32.h"
#include "fflas-ffpack/field/modular-balanced-int32.h"
#include "fflas-ffpack/field/modular-int64.h"
#include "fflas-ffpack/field/unparametric.h"

using namespace FFPACK;
//using namespace FFLAS;


#define TESTE_EG( a, b ) \
if (!F.areEqual((a),(b))) {\
	std::cout << F.write(std::cout,a) << "!=" << F.write(std::cout,b) << " failed (at line " <<  __LINE__ << ")" << std::endl; \
	return(-1); \
}

#define JETESTE( a, s ) \
if (TestField( (a), int(s)) ) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

#define JEONETESTE( F, a, x ) \
if (TestOneField(F,(int)a,(float)x)) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

template<class Int1, class Int2>
long long locgcd ( const Int1 a, const Int2 b ) {
    long long u3, v3; u3 = a; v3 = b;
    while (v3 != 0) {
        long long q, t3;
        q = u3 / v3;
        t3 = u3 - q * v3;
        u3 = v3; v3 = t3;
    }
//     std::cerr << '|' << a << '^' << b << '|' << u3 << std::endl;
    
    return u3;
}


template<class Field>
int TestOneField(const Field& F, const int FIRSTINT, const float FIRSTFLOAT)
{/*{{{*/
#ifdef FFLASFFPACK_DEBUG
	std::cerr << "testing " ;
	F.write(std::cerr );
        std::cerr << " (" << FIRSTINT << ',' << FIRSTFLOAT << ')';
	std::cerr  << " : " << std::flush;

#endif



	typename Field::Element a, b, c, d,a_,b_,c_,d_,ma;
	typename Field::Element e,e_;

    F.init(a, 0UL);
    TESTE_EG(a, F.zero);
    F.init(a, 1UL);
//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "1: ", F.one) << std::endl;
    TESTE_EG(a, F.one);
	F.init(ma,-1L);
//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "ma: ", ma) << std::endl;
//         F.write(std::cerr << "-1: ", F.mOne) << std::endl;
    TESTE_EG(ma, F.mOne);

	F.init(a, FIRSTINT);

    unsigned long invertible=(unsigned long)(FIRSTFLOAT<0?-FIRSTFLOAT:FIRSTFLOAT);

    for( ; locgcd( invertible,F.characteristic()) != 1; ++invertible) {}
    F.init(b, invertible); 
//     F.write(std::cerr << "b:=", b) << ';' << std::endl;
    if (F.isZero(b)) F.init(b,1);

	F.init(c);            // empty constructor
	F.init(d);            // empty constructor

	F.add(c, a, b);       // c = a+b
	F.init(c_);           //! @warning F.init(c_,c); ne marche pas !
	F.assign(c_,c);       // c_ <- c

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "b: ", b) << std::endl;
//         F.write(std::cerr << "c: ", c) << std::endl;
//         F.write(std::cerr << "c_: ", c_) << std::endl;

	TESTE_EG(c,c_);
	F.subin(c_,a);
//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "b: ", b) << std::endl;
//         F.write(std::cerr << "c: ", c) << std::endl;
//         F.write(std::cerr << "c_: ", c_) << std::endl;

	TESTE_EG(b,c_);

	F.mul(c, a, b);     // c = a*b
	F.assign(c_,c);       // c_ <- c
	F.divin(c_,b);      // c_ == a ?

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "c_:=", c_) << ';' << std::endl;
	TESTE_EG(a,c_);

	F.axpy(d, a, b, c); // d = a*b + c;
	F.init(d_);
	F.axmy(d_,a,b,c); // d_ = a*b - c
	F.addin(d_,c);
	F.subin(d,c);

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a:=", a) << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "d:=", d) << ';' << std::endl;
//         F.write(std::cerr << "d_:=", d_) << ';' << std::endl;
	TESTE_EG(d_,d);

	F.sub(d,a,b); // d = a -b
	F.add(c,a,b); // c = a+b
	F.init(e);
	F.init(e_);
	F.mul(e,d,c); // e = d*c;
	F.mul(a_,a,a); // a_ = a*a
	F.mul(b_,b,b); // b_ = b*b
	F.sub(e_,a_,b_); // e_ = a_ - b_

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "d:=", d) << ';' << std::endl;
//         F.write(std::cerr << "e:=", e) << ';' << std::endl;
//         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
//         F.write(std::cerr << "a_:=", a_) << ';' << std::endl;
//         F.write(std::cerr << "b_:=", b_) << ';' << std::endl;
	TESTE_EG(e,e_) // a^2 - b^2 = (a-b)(a+b) ;)

	// Four operations
	F.init(a_);
	F.assign(a_,a);
	F.addin(a, b) ;
	F.subin(a, b) ;
	F.mulin(a, b) ;
	F.divin(a, b) ;

	TESTE_EG(a_,a);


	F.maxpy(e, a, b, d); // e = d-a*b
	F.assign(e_,d);
	F.maxpyin(e_, a, b); // e = d - a*b;

//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "d:=", d) << ';' << std::endl;
//         F.write(std::cerr << "e:=", e) << ';' << std::endl;
//         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
	TESTE_EG(e,e_);


	F.axmy(e, a, b, d); // e = a*b -d;

	F.assign(e_,d);
	F.maxpyin(e_, a, b); // e = d - a*b;

	F.negin(e_);

//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "d:=", d) << ';' << std::endl;
//         F.write(std::cerr << "e:=", e) << ';' << std::endl;
//         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
	TESTE_EG(e,e_);


	if ((unsigned) F.characteristic() != (unsigned)2) {
		F.init(a,22993);
		F.inv(b,a);
		F.mul(c,b,a);

//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "e:=", e) << ';' << std::endl;
		TESTE_EG(F.one,c);

		F.init(a,22993);
		F.init(b,22993);
		F.invin(a);
		F.mulin(a,b);
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
        
        TESTE_EG(F.one,a);

        F.init(a,37409);
        F.inv(b,a);
        F.mul(c,b,a);
        
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "1:=", F.one) << ';' << std::endl;
        TESTE_EG(F.one,c);
        
        F.init(a,37409);
        F.init(b,37409);
        F.invin(a);
        F.mulin(a,b);
        
        TESTE_EG(F.one,a);
	}

#ifdef FFLASFFPACK_DEBUG
	F.write(std::cerr );
	std::cerr  << " done." << std::endl;
#endif
	return 0 ;

}/*}}}*/

#define NBITER 50

template<class Field>
int TestField(const Field& F, const int seed)
{/*{{{*/
    int64_t ch = (int64_t) F.characteristic();
    JEONETESTE(F,7UL,-29.3);
    srand48(seed);
    for(size_t i=0; i< NBITER; ++i) {
        typename Field::Element x;
        float d;
	do {
		d = float((double)ch*drand48());
            F.init(x, (int)d );
        } while(F.isZero(x));
        int a; do {
            F.init(x, a = (int)lrand48());
        } while(F.isZero(x));
        JEONETESTE(F,a,d);
    }
    return 0;
}/*}}}*/



int main(int argc, char ** argv)
{/*{{{*/
    int seed = int (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef FFLASFFPACK_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    srand48(seed);

#ifdef NDEBUG
    assert(0);
#endif

	// modulo 13 over 16 bits
	Modular<float> C13(13);
	JETESTE(C13,seed);

    // modulo 13 over 32 bits
    Modular<double> Z13(13);
    JETESTE(Z13,seed);

	// modulo 13
	ModularBalanced<float> U13(13);
	JETESTE(U13,seed);

	// modulo 13 
	ModularBalanced<double> M13(13);
	JETESTE(M13,seed);

	// modulo 13 
	Modular<int32_t> L13(13);
	JETESTE(L13,seed);

	// modulo 13 over 64 bits
	Modular<int64_t> LL13(13UL);
	JETESTE(LL13,seed);

	// modulo 13 
	ModularBalanced<int32_t> Lb13(13);
	JETESTE(Lb13,seed);

	// modulo 13 over 64 bits
	ModularBalanced<int64_t> LLb13(13UL);
	JETESTE(LLb13,seed);



// // Maximal values

// 	// prime modulo max
    Modular<float> CUmax(Modular<float>::getMaxModulus() );
    JETESTE(CUmax,seed);

    Modular<double> Zmax( Modular<double>::getMaxModulus() );
    JETESTE(Zmax,seed);

	ModularBalanced<float> Umax(ModularBalanced<float>::getMaxModulus() );
	JETESTE(Umax,seed);

	ModularBalanced<double> Mmax(ModularBalanced<double>::getMaxModulus());
	JETESTE(Mmax,seed);

	Modular<int32_t> Lmax(Modular<int32_t>::getMaxModulus());
	JETESTE(Lmax,seed);

	Modular<int64_t> LLmax(Modular<int64_t>::getMaxModulus());
	JETESTE(LLmax,seed);

    ModularBalanced<int32_t> Lbmax(ModularBalanced<int32_t>::getMaxModulus());
    JETESTE(Lbmax,seed);

    ModularBalanced<int64_t> LLbmax(ModularBalanced<int64_t>::getMaxModulus());
    JETESTE(LLbmax,seed);



// // Characteristic 2


// 	// modulo 2 over 16 bits
	Modular<float> C2(2);
	JETESTE(C2,seed);

// 	// modulo 2 over 32 bits
    Modular<double> Z2(2);
    JETESTE(Z2,seed);

//     ModularBalanced<float> U2(2);
//     JETESTE(U2,seed);

//     ModularBalanced<double> M2(2);
//     JETESTE(M2,seed);
	
	Modular<int32_t> L2(2);
	JETESTE(L2,seed);

	Modular<int64_t> LL2(2UL);
	JETESTE(LL2,seed);

// 	ModularBalanced<int32_t> L2b( 2 );
// 	JETESTE(L2b,seed);

// 	ModularBalanced<int64_t> LL2b( 2 );
// 	JETESTE(LL2b,seed);

// // Random values

    for(int i=0; i< 20; ++i) {

        long a = lrand48();
//         std::cerr << "rand int: " << a << std::endl;
        

    Modular<float> CUrand( (float)(a % (long)Modular<float>::getMaxModulus() ));
    JETESTE(CUrand,seed);

    Modular<double> Zrand((double)(a %  (long)Modular<double>::getMaxModulus() ));
    JETESTE(Zrand,seed);

	ModularBalanced<float> Urand((float)(a % (long)ModularBalanced<float>::getMaxModulus() ));
	JETESTE(Urand,seed);

	ModularBalanced<double> Mrand((double)(a % (long)ModularBalanced<double>::getMaxModulus()));
	JETESTE(Mrand,seed);

	Modular<int32_t> Lrand((int32_t)(a % Modular<int32_t>::getMaxModulus()));
	JETESTE(Lrand,seed);

	Modular<int64_t> LLrand((int64_t)(a % Modular<int64_t>::getMaxModulus()));
	JETESTE(LLrand,seed);

    ModularBalanced<int32_t> Lbrand((int32_t)(a % ModularBalanced<int32_t>::getMaxModulus()));
    JETESTE(Lbrand,seed);

    ModularBalanced<int64_t> LLbrand((int64_t)(a % ModularBalanced<int64_t>::getMaxModulus()));
    JETESTE(LLbrand,seed);

    
    }
    

	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
