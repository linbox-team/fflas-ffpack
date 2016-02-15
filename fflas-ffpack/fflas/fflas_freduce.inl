/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* fflas/fflas_freduce.inl
 * Copyright (C) 2014 Pascal Giorgi
 *
 * Written by Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * Part of this code is taken from http://libdivide.com/
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

#ifndef __FFLASFFPACK_fflas_freduce_INL
#define __FFLASFFPACK_fflas_freduce_INL

#include <givaro/udl.h>

#include "fflas-ffpack/fflas/fflas_fassign.h"
#include "fflas-ffpack/utils/bit_manipulation.h"

#define FFLASFFPACK_COPY_REDUCE 32 /*  TO BENCMARK LATER */



namespace FFLAS { namespace vectorised { /*  for casts (?) */

	template<class T>
	inline typename std::enable_if< ! std::is_integral<T>::value, T>::type
	monfmod(T A, T B)
	{
		return fmod(A,B);
	}

	template<class T>
	inline typename std::enable_if< std::is_integral<T>::value, T>::type
	monfmod(T  A, T B)
	{
		return A % B; // B > 0
	}

	template<>
	inline Givaro::Integer monfmod(Givaro::Integer  A, Givaro::Integer B) // @bug B is not integer, but uint64_t usually
	{
		return A % B; // B > 0
	}

	template<>
	inline float monfmod(float A, float B)
	{
		return fmodf(A,B);
	}

	template<>
	inline double monfmod(double A, double B)
	{
		    //std::cerr<<"fmod"<<std::endl;
		return fmod(A,B);
	}

	template<size_t K, size_t MG>
	inline RecInt::rmint<K,MG>& monfmod(RecInt::rmint<K,MG>& A, RecInt::rmint<K,MG>& B)
	{
		return RecInt::rmint<K>::mod_n(A, B);
	}

	template<class T>
	inline typename std::enable_if< ! std::is_integral<T>::value, T>::type
	monrint(T A)// @bug pass by reference ?
	{
		return rint(A);
	}

	template<class T>
	inline typename std::enable_if< std::is_integral<T>::value, T>::type
	monrint( T A)
	{
		return A ;
	}

	template<>
	inline double monrint(double A)
	{
		return rint(A);
	}


	template<>
	inline float monrint(float A)
	{
		return rintf(A);
	}

	template<>
	inline Givaro::Integer monrint(Givaro::Integer A) // @bug B is not integer, but uint64_t usually
	{
		return A ; // B > 0
	}


	template<bool overflow, bool poweroftwo>
	inline int64_t monfmod(int64_t A, int64_t p, int8_t shifter, int64_t magic)
	{
		if (poweroftwo) { //shift path
			int64_t q = A + ((A >> 63) & ((1_i64 << shifter) - 1));
			q = A - ((q>>shifter)<< shifter) ;
			return (q<0)?(q+p):q ;
		}
		else {
			int64_t q = mulhi_64(magic, A);
			if (overflow) {
				q += A ;
			}
			q >>= shifter;
			A = A - q * p ;
			if (A >= p) A-= p ; // because of mulhi_fast
			return A ;

		}
	}

} // vectorised
} // FFLAS

namespace FFLAS { namespace vectorised {


	template<class T>
	inline void fast_mod_generate(bool & overflow, bool & poweroftwo, int8_t & shift, T & magic, T denom)
	{
		overflow = false ;
		poweroftwo = false ;
		shift = 0 ;
		magic = 0 ;
	}

	//! @pre d > 0
	template<>
	inline void fast_mod_generate(bool & overflow, bool & poweroftwo, int8_t & shift, int64_t & magic, int64_t denom)
	{

		// overflow = false ;
		// poweroftwo = false ;
		// shift = 0 ;
		// magic = 0 ;
		if ((denom & (denom- 1)) == 0) {
			shift = (int8_t)ctz((uint64_t)denom) ;
			magic = 0;
			poweroftwo = true ;
		}
		else {
			const uint32_t floor_log_2_d = 63 - clz((uint64_t)denom);

			/*the dividend here is 2**(floor_log_2_d + 63), so the low 64 bit word is 0 and the high word is floor_log_2_d - 1 */
			uint64_t rem, proposed_m;

                        proposed_m = getpoweroftwoden_128(floor_log_2_d, denom, &rem);

			const uint64_t e = denom- rem;

			/* We are going to start with a power of floor_log_2_d - 1.  This works if works if e < 2**floor_log_2_d. */
			if (e < (1_ui64 << floor_log_2_d)) {
				/* This power works */
				shift = (int8_t)(floor_log_2_d - 1);
			}
			else {
				/* We need to go one higher.  This should not make proposed_m overflow, but it will make it negative when     interpreted as  an int32_t. */
				proposed_m += proposed_m;
				const uint64_t twice_rem = rem + rem;
				if (twice_rem >= (uint64_t)denom || twice_rem < rem) proposed_m += 1;
				shift = (int8_t) floor_log_2_d  ;
				overflow = true ;
			}
			proposed_m += 1;
			magic = (int64_t)proposed_m ;
		}
	}

	template<class Field, class ElementTraits = typename ElementTraits<typename Field::Element>::value>
	struct HelperMod  ;


	template<class Field>
	struct HelperMod<Field, ElementCategories::MachineIntTag> {
		bool overflow  = false ;
		bool poweroftwo = false ;
		int8_t shift = 0 ;
		typename Field::Element magic = (typename Field::Element)0 ;
		typename Field::Element p;

		HelperMod()
		{
			// std::cout << "empty cstor called" << std::endl;
		} ;

		HelperMod( const Field & F)
		{
			// std::cout << "field cstor called" << std::endl;
			p =  (typename Field::Element) F.characteristic();
			fast_mod_generate(overflow, poweroftwo, shift, magic, p);
			// std::cout << overflow << ',' << poweroftwo << std::endl;
			// std::cout << (int) shift << ',' << magic << std::endl;
			// std::cout << this->shift << std::endl;
		}

		int getAlgo() const
		{
			// std::cout << "will be " << (2*overflow + poweroftwo) << std::endl;
			return 2* (int)overflow + (int) poweroftwo ;
			// return overflow << 1 | poweroftwo ;
		}


	} ;

	template<class Field>
	struct HelperMod<Field, FFLAS::ElementCategories::MachineFloatTag> {
		typename Field::Element p;
		typename Field::Element invp;
		// typename Field::Elmeent min ;
		// typename Field::Elmeent max ;

		HelperMod() {} ;

		HelperMod( const Field & F)
		{
			p = (typename Field::Element) F.characteristic();
			invp = (typename Field::Element)1/p;
			// min = F.minElement();
			// max = F.maxElement();
		}

		int getAlgo() const
		{
			return 0;
		}
	} ;

	template<class Field>
	struct HelperMod<Field, FFLAS::ElementCategories::ArbitraryPrecIntTag> {
		typename Field::Element p;
		// typename Field::Element invp;
		// typename Field::Elmeent min ;
		// typename Field::Elmeent max ;

		HelperMod() {} ;

		HelperMod( const Field & F)
		{
			p = (typename Field::Element) F.characteristic();
			// invp = (typename Field::Element)1/p;
			// min = F.minElement();
			// max = F.maxElement();
		}

		int getAlgo() const
		{
			return 0;
		}
	} ;

	template<class Field>
	struct HelperMod<Field, FFLAS::ElementCategories::FixedPrecIntTag> {
		typename Field::Element p;
		// typename Field::Element invp;
		// typename Field::Elmeent min ;
		// typename Field::Elmeent max ;

		HelperMod() {} ;

		HelperMod( const Field & F)
		{
			p = (typename Field::Element) F.characteristic();
			// invp = (typename Field::Element)1/p;
			// min = F.minElement();
			// max = F.maxElement();
		}

		int getAlgo() const
		{
			return 0;
		}
	} ;


#ifdef __FFLASFFPACK_USE_SIMD
	template<class Field, class SimdT, class ElementTraits = typename ElementTraits<typename Field::Element>::value>
	struct HelperModSimd  ;

	template<class Field, class SimdT>
	struct HelperModSimd<Field, SimdT, ElementCategories::MachineIntTag> : public HelperMod<Field> {
		typedef typename SimdT::vect_t vect_t ;
		// bool overflow ;
		// int8_t shift ;
		// typename Field::Element p;
		typename Field::Element magic ;
		vect_t M ;
		vect_t P ;
		vect_t MIN ;
		vect_t MAX ;
		vect_t NEGP ;
		vect_t Q ;
		vect_t T ;

		HelperModSimd ( const Field & F) :
			HelperMod<Field>(F)
		{
			// std::cout << "HelperMod constructed " << this->shift << std::endl;
			// p = F.characteristic();
			P = SimdT::set1(this->p);
			NEGP = SimdT::set1(-this->p);
			MIN = SimdT::set1(F.minElement());
			MAX = SimdT::set1(F.maxElement());
			// fast_mod_generate(overflow, shift, magic, p);
			M = SimdT::set1(magic);
		}

		HelperModSimd( const Field & F, const HelperMod<Field> & G)
		{
			this->overflow=G.overflow;
			this->poweroftwo=G.poweroftwo;
			this->shift=G.shift;
			this->magic=G.magic;
			this->p=G.p;
			// std::cout << "magic is = " << this->magic<< ',' <<  G.magic<< std::endl;
			P = SimdT::set1(this->p);
			NEGP = SimdT::set1(-(this->p));
			MIN = SimdT::set1(F.minElement());
			MAX = SimdT::set1(F.maxElement());
			// fast_mod_generate(overflow, shift, magic, p);
			M = SimdT::set1(magic);
		}

	} ;

	template<class Field, class SimdT>
	struct HelperModSimd<Field, SimdT, ElementCategories::MachineFloatTag>  : public HelperMod<Field> {
		typedef typename SimdT::vect_t vect_t ;
		vect_t INVP;
		vect_t MIN ;
		vect_t MAX ;
		vect_t NEGP ;
		vect_t P ;
		vect_t Q ;
		vect_t T ;

		HelperModSimd( const Field & F) :
			HelperMod<Field>(F)
		{
			P = SimdT::set1(this->p);
			NEGP = SimdT::set1(-(this->p));
			// MIN = SimdT::set1(max);
			MIN = SimdT::set1(F.minElement());
			// MAX = SimdT::set1(min);
			MAX = SimdT::set1(F.maxElement());
			INVP = SimdT::set1(this->invp);
		}

		HelperModSimd( const Field & F, const HelperMod<Field> & G)
		{
			this->p = G.p;
			this->invp = G.invp ;
			P = SimdT::set1(this->p);
			NEGP = SimdT::set1(-this->p);
			// MIN = SimdT::set1(max);
			MIN = SimdT::set1(F.minElement());
			// MAX = SimdT::set1(min);
			MAX = SimdT::set1(F.maxElement());
			INVP = SimdT::set1(this->invp);


		}
	} ;
#endif // __FFLASFFPACK_USE_SIMD


#ifdef __x86_64__
	template<class Field, int ALGO>
	typename std::enable_if< std::is_same<typename Field::Element,int64_t>::value , int64_t>::type
	monfmod (typename Field::Element A, HelperMod<Field,ElementCategories::MachineIntTag> & H)
	{
		switch(ALGO) {
		case 3 :
			// std::cout << 3 << std::endl;
			return monfmod<true,true>  (A,H.p,H.shift,H.magic);
		case 2 :
			// std::cout << 2 << std::endl;
			return monfmod<true,false> (A,H.p,H.shift,H.magic);
		case 1 :
			// std::cout << 1 << std::endl;
			return monfmod<false,true> (A,H.p,H.shift,H.magic);
		case 0 :
			    //		std::cout << "using " << 0 << std::endl;
			return monfmod<false,false>(A,H.p,H.shift,H.magic);
		default :
			FFLASFFPACK_abort("unknown algo");
		}
	}
#endif // __x86_64__


	template<class Field, int ALGO>
#ifdef __x86_64__
	typename std::enable_if< ! std::is_same<typename Field::Element,int64_t>::value , typename Field::Element>::type
#else
	typename Field::Element
#endif // __x86_64__
	monfmod (typename Field::Element A, HelperMod<Field,ElementCategories::MachineIntTag> & H)
	{
		return monfmod(A,H.p);
	}

	template<class Field, int ALGO>
	typename Field::Element monfmod (typename Field::Element A, HelperMod<Field,ElementCategories::MachineFloatTag> & H)
	{
		return monfmod(A,H.p);
	}

	template<class Field, int ALGO>
	typename Field::Element monfmod (typename Field::Element A, HelperMod<Field,ElementCategories::ArbitraryPrecIntTag> & H)
	{
		return monfmod(A,H.p);
	}



#ifdef __FFLASFFPACK_USE_SIMD

	template<class Field, class SimdT, int ALGO>
	inline void
	VEC_MOD(typename SimdT::vect_t & C, HelperModSimd<Field,SimdT,ElementCategories::MachineFloatTag> & H)
	{
		C = SimdT::mod( C, H.P, H.INVP, H.NEGP, H.MIN, H.MAX, H.Q, H.T );
	}

	template<class Field, class SimdT, int ALGO>
	inline void
	VEC_MOD(typename SimdT::vect_t & C, HelperModSimd<Field,SimdT,ElementCategories::MachineIntTag> & H)
	{
		// std::cout << "magic " << H.magic<< std::endl;
		// std::cout << H.P << std::endl;
		switch (ALGO) {
		case 0 :
			C = SimdT::template mod<false,false>( C, H.P, H.shift, H.M,  H.NEGP, H.MIN, H.MAX, H.Q, H.T );
			break;
		case 1 :
			C = SimdT::template mod<true,false> ( C, H.P, H.shift, H.M,  H.NEGP, H.MIN, H.MAX, H.Q, H.T );
			break;
		case 2 :
			C = SimdT::template mod<false,true> ( C, H.P, H.shift, H.M,  H.NEGP, H.MIN, H.MAX, H.Q, H.T );
			break;
		case 3 :
			C = SimdT::template mod<true,true>  ( C, H.P, H.shift, H.M,  H.NEGP, H.MIN, H.MAX, H.Q, H.T );
			break;
		}
	}

#endif // __FFLASFFPACK_USE_SIMD

} // vectorised
} // FFLAS

namespace FFLAS  { namespace vectorised { namespace unswitch  {

#ifdef __FFLASFFPACK_USE_SIMD
	template<class Field, bool round, int algo>
	inline typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n,
	     typename Field::Element_ptr T
	     , HelperMod<Field> & G
	     )
	{

//		std::cerr<<"modp vectorized"<<std::endl;
		typedef typename Field::Element Element;
		Element min = (Element)F.minElement(), max = (Element)F.maxElement();
		using simd = Simd<Element>;
		using vect_t = typename simd::vect_t;
		bool positive = ! FieldTraits<Field>::balanced ; // known at compile time
		HelperModSimd<Field,simd> H(F,G);

		size_t i = 0;
		if (n < simd::vect_size)
		{
//			std::cerr<< n<< " < "<<simd::vect_size<<std::endl;
			for (; i < n ; i++)
			{
				if (round)
				{
					T[i] = monrint(U[i]);
					T[i] = monfmod<Field,algo>(T[i],H);
				}
				else
				{
					T[i]=monfmod<Field,algo>(U[i],H);
				}
				if (!positive)
				{
					T[i]-=(T[i]>max)?H.p:0;
				}
				T[i]+=(T[i]<min)?H.p:0;
			}
			return;
		}

		long st = long(T) % simd::alignment;

		// the array T is not 32 byte aligned (process few elements s.t. (T+i) is 32 bytes aligned)

		if (st)
		{
//			std::cerr<< st << " not aligned on  "<<simd::alignment<<std::endl;

			for (size_t j = static_cast<size_t>(st) ; j < simd::alignment ; j += sizeof(Element), i++)
			{
				if (round)
				{
					T[i] = monrint(U[i]);
					T[i] = monfmod<Field,algo>(T[i],H);
				}
				else
				{
					T[i] = monfmod<Field,algo>(U[i],H);
				}
				if (!positive)
				{
					T[i] -= (T[i] > max) ? H.p : 0;
				}
				T[i] += (T[i] < min) ? H.p : 0;
			}
		}

		FFLASFFPACK_check((long(T+i) % simd::alignment == 0));

		vect_t C ;

		if((long(U+i) % simd::alignment == 0))
		{
			// perform the loop using 256 bits SIMD
			for (; i<= n - simd::vect_size ; i += simd::vect_size)
			{
				C = simd::load(U + i);

				if (round)
				{
					C = simd::round(C);
				}

				VEC_MOD<Field,simd,algo>(C,H);
				simd::store(T+i, C);
			}
		}

		// perform the last elt from T without SIMD
//		std::cerr<< n-i<< " unaligned elements left "<<std::endl;
		for (;i<n;i++)
		{

			if (round)
			{
				T[i] = monrint(U[i]);
				T[i] = monfmod<Field,algo>(T[i],H);
			}
			else
			{
				T[i] = monfmod<Field,algo>(U[i],H);
			}
			if (!positive)
			{
				T[i] -= (T[i] > max) ? H.p : 0;
			}
			T[i] += (T[i] < min) ? H.p : 0;
		}
	}
#endif

	// not vectorised but allows better code than % or fmod via helper
	template<class Field, bool round, int algo>
	inline typename std::enable_if< !FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n,
	     typename Field::Element_ptr T
	     , HelperMod<Field> & H
	     )
	{
//		std::cerr<<"modp not vectorized"<<std::endl;
		typedef typename Field::Element Element;
		Element min = (Element)F.minElement(), max = (Element)F.maxElement();
		bool positive = ! FieldTraits<Field>::balanced ;

		size_t i = 0;
		for (; i < n ; i++)
		{
			if (round)
			{
				T[i] = monrint(U[i]);
				T[i] = monfmod<Field,algo>(T[i],H);
			}
			else
			{
				T[i]=monfmod<Field,algo>(U[i],H);
			}
			if (!positive)
			{
				T[i]-=(T[i]>max)?H.p:(typename Field::Element)0;
			}
			T[i]+=(T[i]<min)?H.p:(typename Field::Element)0;
		}
	}

} // unswitch
} // vectorised
} // FFLAS

namespace FFLAS { namespace vectorised {


	template<class Field, bool round>
	//inline typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	void
	modp(const Field &F, typename Field::ConstElement_ptr U, const size_t & n,
	     typename Field::Element_ptr T)
	{
		HelperMod<Field> H(F);

		int ALGO = H.getAlgo();

		switch (ALGO) {
		case 0 :
			unswitch::modp<Field,round,0>(F,U,n,T,H);
			break;
		case 1 :
			unswitch::modp<Field,round,1>(F,U,n,T,H);
			break;
		case 2 :
			unswitch::modp<Field,round,2>(F,U,n,T,H);
			break;
		case 3 :
			unswitch::modp<Field,round,3>(F,U,n,T,H);
			break;
		}
	}

} // vectorised
} // FFLAS


namespace FFLAS { namespace details {


	// specialised
	template<class Field>
	typename std::enable_if<FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	freduce (const Field & F, const size_t m,
		 typename Field::Element_ptr A, const size_t incX, FieldCategories::ModularTag)
	{
		if(incX == 1) {
			vectorised::modp<Field,false>(F,A,m,A);
		}
		else { /*  faster with copy, use incX=1, copy back ? */
			if (m < FFLASFFPACK_COPY_REDUCE) {
				typename Field::Element_ptr Xi = A ;
				for (; Xi < A+m*incX; Xi+=incX )
					F.reduce(*Xi);
			}
			else {
				typename Field::Element_ptr Ac = fflas_new (F,m,1) ;
				fassign (F,m,A,incX,Ac,1);
				freduce (F,m,Ac,1,FieldCategories::ModularTag());
				fassign (F,m,Ac,1,A,incX);
				fflas_delete (Ac);
			}
		}
	}

	template<class Field>
	typename std::enable_if< ! FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	freduce (const Field & F, const size_t m,
		 typename Field::Element_ptr A, const size_t incX, FieldCategories::ModularTag)
	{ /* ??? ( faster with copy, use incX=1, copy back ? */
		    // CP: no SIMD supported here!
                // if(incX == 1) {
		// 	vectorised::modp<Field,false>(F,A,m,A);
		// }
		// else {
			typename Field::Element_ptr  Xi = A ;
			for (; Xi < A+m*incX; Xi+=incX )
				F.reduce(*Xi);
		// }
	}

		template<class Field>
	void
	freduce (const Field & F, const size_t m,
		 typename Field::Element_ptr A, const size_t incX,
		 FieldCategories::GenericTag)
	{
		typename Field::Element_ptr Xi = A ;
		for (; Xi < A+m*incX; Xi+=incX )
			F.reduce (*Xi);
	}

	template<class Field>
	void
	freduce (const Field & F, const size_t m,
		 typename Field::Element_ptr A, const size_t incX,
		 FieldCategories::UnparametricTag)
	{
		typename Field::Element_ptr Xi = A ;
		for (; Xi < A+m*incX; Xi+=incX )
			F.reduce (*Xi);
	}

	template<class Field>
	typename std::enable_if< FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	freduce (const Field & F, const size_t m,
		 typename Field::ConstElement_ptr  B, const size_t incY,
		 typename Field::Element_ptr A, const size_t incX,
		 FieldCategories::ModularTag)
	{
		
		if(incX == 1 && incY == 1) {
			vectorised::modp<Field,false>(F,B,m,A);
		}
		else {
			typename Field::Element_ptr Xi = A ;
			typename Field::ConstElement_ptr Yi = B ;
			for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
				F.reduce (*Xi , *Yi);
		}
	}

	template<class Field>
	typename std::enable_if< ! FFLAS::support_simd_mod<typename Field::Element>::value, void>::type
	freduce (const Field & F, const size_t m,
		 typename Field::ConstElement_ptr  B, const size_t incY,
		 typename Field::Element_ptr A, const size_t incX,
		 FieldCategories::ModularTag)
	{

		typename Field::Element_ptr Xi = A ;
		typename Field::ConstElement_ptr Yi = B ;
		for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
			F.reduce (*Xi , *Yi);
	}

	template<class Field>
	void
	freduce (const Field & F, const size_t m,
		 typename Field::ConstElement_ptr  B, const size_t incY,
		 typename Field::Element_ptr A, const size_t incX,
		 FieldCategories::GenericTag)
	{
		typename Field::Element_ptr Xi = A ;
		typename Field::ConstElement_ptr Yi = B ;
		for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
			F.reduce (*Xi , *Yi);
	}
	template<class Field>
	void
	freduce (const Field & F, const size_t m,
		 typename Field::ConstElement_ptr  B, const size_t incY,
		 typename Field::Element_ptr A, const size_t incX,
		 FieldCategories::UnparametricTag)
	{
		typename Field::Element_ptr Xi = A ;
		typename Field::ConstElement_ptr Yi = B ;
		for (; Xi < A+m*incX; Xi+=incX, Yi += incY )
			F.reduce (*Xi , *Yi);
	}

} // details
} // FFLAS


#endif // __FFLASFFPACK_fflas_freduce_INL

