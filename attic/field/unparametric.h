/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* field/unparametric.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 *               2005 Clement Pernet
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified By C. Pernet and inserted into Fflas_Ffpack
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
 *
 */


/*! @file field/unparametric.h
 * @ingroup field
 * @brief  representation of a field of characteristic 0.
 */

#ifndef __FFLASFFPACK_field_unparametric_H
#define __FFLASFFPACK_field_unparametric_H

#include <iostream> // std::cout
#include <string>
#include <algorithm>
#include <typeinfo>

#include "fflas-ffpack/field/unparametric-randiter.h"
#include "fflas-ffpack/field/nonzero-randiter.h"
#include "fflas-ffpack/field/integer.h"

namespace FFPACK
{
	template<class _Element>
	class UnparametricField ;

	template <typename Target, typename Source>
	Target& Caster (Target& t, const Source& s)
	{
		return t = static_cast<Target>(s);
	}

	/** \brief Unparameterized field adapter.
	 * \ingroup field
	 * \defgroup UnparametricField UnparametricField
	 *
	 * A field having an interface similar to that of floats is adapted to LinBox.
	 *
	 *  Used to generate efficient field classes for unparameterized fields (or hidden parameter fields).
	 *
	 *  Some fields are implemented by definition of the C++ arithmetic operators, such as z = x*y,
	 *  for z, y, z instances of a type K.   The LinBox field
	 *  Unparametric<K> is the adaptation to LinBox.
	 *
	 *  For a typical unparametric field, some of the methods must be defined in a specialization.
	 */
	template <class K>
	class UnparametricOperations {
	public:
		typedef K Element;
        
		UnparametricOperations(){}
		//@{
		~UnparametricOperations () {}

		/* Assignment operator.
		 * Assigns UnparametricField object F to field.
		 * @param  F UnparametricField object.
		 */
		// I believe this should be virtual -bds
		///
		//@} Field Object Basics.

		/** @name Data Object Management.
		 * first argument is set and the value is also returned.
		 */
		//@{

		Element& init (Element& x) const
		{
			return x;
		}


		///
		Element &assign (Element &x, const Element &y) const
		{
			return x = y;
		}


		//@}

		/// @name Comparison Predicates
		//@{
		///  x == y
		bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		///  x == 0
		bool isZero (const Element &x) const
		{
			return x == Element (0);
		}

		///  x == 1
		bool isOne (const Element &x) const
		{
			return x == Element (1);
		}

		inline bool isMOne (const Element &x) const
		{
			return x == Element(-1) ;
		}

		//@} Comparison Predicates


		/** @name Arithmetic Operations
		 * The first argument is set and is also the return value.
		 */
		//@{

		/// x := y + z
		Element &add (Element &x, const Element &y, const Element &z) const
		{
			return x = y + z;
		}

		/// x := y - z
		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			return x = y - z;
		}

		/// x := y*z
		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			return x = y * z;
		}

		/// x := y/z
		Element &div (Element &x, const Element &y, const Element &z) const
		{
			return x = y / z;
		}

		/// x := -y
		Element &neg (Element &x, const Element &y) const
		{
			return x = - y;
		}

		/// x := 1/y
		Element &inv (Element &x, const Element &y) const
		{
			return x = Element (1) / y;
		}

		/// z := a*x + y
		// more optimal implementation, if available, can be defined in a template specialization.
		Element &axpy (Element &z,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			return z = a * x + y;
		}

		//@} Arithmetic Operations

		/** @name Inplace Arithmetic Operations
		 * The first argument is modified and the result is the return value.
		 */
		//@{

		/// x := x + y
		Element &addin (Element &x, const Element &y) const
		{
			return x += y;
		}

		/// x := x - y
		Element &subin (Element &x, const Element &y) const
		{
			return x -= y;
		}

		/// x := x*y
		Element &mulin (Element &x, const Element &y) const
		{
			return x *= y;
		}

		/// x := x/y
		Element &divin (Element &x, const Element &y) const
		{
			return x /= y;
		}

		/// x := -x
		Element &negin (Element &x) const
		{
			return x = - x;
		}

		/// x := 1/x
		Element &invin (Element &x) const
		{
			return x = Element (1) / x;
		}

		/// y := a*x + y
		Element &axpyin (Element &y, const Element &a, const Element &x) const
		{
			return y += a * x;
		}

		//@} Inplace Arithmetic Operations

		/** @name Input/Output Operations */
		//@{

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream &write (std::ostream &os) const
		{
			return os << "unparameterized field(" << sizeof(Element) <<',' << typeid(Element).name() << ')';
		}

		std::ostream &write (std::ostream &os, std::string F) const
		{
            return this->write(F != "" ? os << F : os);
		}


		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream &read (std::istream &is) const
		{
			return is;
		}

		/** Print field element.
		 * @return output stream to which field element is written.
		 * @param  os  output stream to which field element is written.
		 * @param  x   field element.
		 */
		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}
		/** Read field element.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		virtual std::istream &read (std::istream &is, Element &x) const
		{
			return is >> x;
		}

		//@}



	};

	    /** Class ZZ
	     * Ring of integers, using the templatedElement base type
	     * Provides the implementation of a field/ring (of char >0 or 0) using 
	     */
	template<class _Element>
	class ZZ : public UnparametricOperations<_Element> {

	public:

		/** The field's element type.
		 * Type Element must provide a default constructor,
		 * a copy constructor, a destructor, and an assignment operator.
		 */

		typedef _Element Element;
		typedef ZZ<Element> Self_t;
		typedef UnparametricRandIter<Element> RandIter;
		typedef NonzeroRandIter<Self_t, RandIter > NonZeroRandIter;
		typedef Element* Element_ptr;
		typedef const Element* ConstElement_ptr;
		const Element one  ;
		const Element zero ;
		const Element mOne ;

		/** @name Field Object Basics.
		*/
		//@{

		/** Builds this field. Assumes q=0 and e = 1.
		 * Ensures consistency with field interface.
		 */
		ZZ(long int q = 0, size_t e = 1):
				one(1),zero(0),mOne(-one)
				    // ,one(Element(1L)),zero(Element(0L)),mOne(Element(-1L))
			{}  
		//@}

		template<class T>
		ZZ (const T& ) : one(1), zero(0), mOne(-one){}

		/// construct this field as copy of F.
		ZZ (const ZZ &F) : one(F.one),zero(F.zero),mOne(F.mOne){}

		FFPACK::Integer &cardinality (FFPACK::Integer &c) const {return c = 0;}

		FFPACK::Integer &characteristic (FFPACK::Integer &c) const{return c = 0;}

		int64_t &cardinality (int64_t &c) const {return c = 0;}

		int64_t &characteristic (int64_t &c) const {return c = 0;}

		int64_t cardinality () const {return 0;}

		int64_t characteristic () const {return 0;}

		ZZ<Element> operator=(const ZZ<Element> &e) {return *this ;}

		/// x := y.
		template <typename Src>
		Element& init (Element& x, const Src& s) const {return FFPACK::Caster (x, s);}

		Element& reduce (Element& x, const Element& y) const {return init (x,y);}

		Element& reduce (Element& x) const {return init (x,x);}

		/// x :=  y.

		template <typename T>
		T& convert (T &x, const Element &y) const {return FFPACK::Caster (x, y);}

		    // To ensure interface consistency
		size_t minElement() const {return 0;}
		size_t maxElement() const {return 0;}

	};

	    /** Class UnparametricField
	     * Provides the implementation of a field/ring (of char >0 or 0) using 
	     * the arithmetic provided by UnparametricOperations.
	     * UnparametricOperations is an empty class wrapping infix +,*,-,/,etc 
	     * operations into add, mul, sub, div, etc.
	     * UnparametricField contains characteristic and cardinality members and accessors.
	     * This is used for instance to represent Z/pZ with NTL's implementation
	     * - UnparametricField<NTL::zz_p>
	     * @warning: prefer ZZ<double>, ZZ<integer> over UnparametricField<double> 
	     * and UnparametricField<integer>
	     */
	template<class _Element>
	class UnparametricField : public UnparametricOperations<_Element> {
	protected:
		long int _p ; long int _card ;
	public:

		/** The field's element type.
		 * Type K must provide a default constructor,
		 * a copy constructor, a destructor, and an assignment operator.
		 */

		typedef typename UnparametricOperations<_Element>::Element Element;
        typedef Element FieldInt;
		typedef UnparametricRandIter<Element> RandIter;
		typedef NonzeroRandIter<UnparametricField<Element>, RandIter > NonZeroRandIter;
		typedef Element* Element_ptr;
		typedef const Element* ConstElement_ptr;

		const Element one  ;
		const Element zero ;
		const Element mOne ;

		/** @name Field Object Basics.
		*/
		//@{

		/** Builds this field to have characteristic q and cardinality q<sup>e</sup>.
		 *  This constructor must be defined in a specialization.
		 */
		UnparametricField(long int q = 0, size_t e = 1) :
			_p(q), _card((long)(q == 0 ? -1 : pow((double)q, (double)e)) )
			// ,one(Element(1L)),zero(Element(0L)),mOne(Element(-1L))
			,one(1),zero(0),mOne(-one)
			{
				// Caster(one,1);
			}  // assuming q is a prime or zero.
		//@}

		template<class T>
		UnparametricField (const T& ) : _p(0), _card(-1), one(1), zero(0), mOne(-one){}

		/// construct this field as copy of F.
		UnparametricField (const UnparametricField &F) :
			_p(F._p), _card(F._card)
			// ,one(1L),zero(0L)
			,one(F.one),zero(F.zero),mOne(F.mOne)
		{
			// init(mOne,-1L);
		}
#ifdef __FFLASFFPACK_HAVE_INTEGER
		FFPACK::Integer &cardinality (FFPACK::Integer &c) const
		{
			return c = _card ;
		}

		FFPACK::Integer &characteristic (FFPACK::Integer &c) const
		{
			return c = (unsigned long)_p ;
			// return c = _card ;
		}
#endif

		unsigned long &cardinality (unsigned long &c) const
		{
			return c = _card ;
		}

		unsigned long &characteristic (unsigned long &c) const
		{
			return c = (unsigned long)_p ;
			// return c = _card ;
		}

		unsigned long cardinality () const
		{
			return (unsigned long)_card ;
		}

		unsigned long characteristic () const
		{
			return (unsigned long)_p ;
			// return  _card ;
		}

		UnparametricField<Element> operator=(const UnparametricField<Element> &e)
		{
			_p = e.characteristic() ;
			_card = e.cardinality();
			return *this ;
		}

		/// x := y.  Caution: it is via cast to long.  Good candidate for specialization.
		template <typename Src>
		Element& init (Element& x, const Src& s) const
		{
			return FFPACK::Caster (x, s);
		}
		
		Element& reduce (Element& x, const Element& y) const {return init (x,y);}

		Element& reduce (Element& x) const {return init (x,x);}

		/// x :=  y.  Caution: it is via cast to long.  Good candidate for specialization. --dpritcha

		template <typename T>
		T& convert (T &x, const Element &y) const
		{
			return FFPACK::Caster (x,y);
		}

		size_t minElement() const { return 0 ; }
		size_t maxElement() const { return 0 ; }
	};
} // FFPACK

namespace FFLAS {
/* Representations of Z with floating point elements*/

	typedef FFPACK::UnparametricField<float> FloatDomain;
	typedef FFPACK::UnparametricField<double> DoubleDomain;

}

#include "field-general.h"

#endif // __FIELD_UNPARAMETRIC_H_
