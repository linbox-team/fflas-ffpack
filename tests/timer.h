/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* test/timer.h
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added _start_t member to BaseTimer, so that stop () does not clobber the
 * class' memory of its start time. This allows it to be called repeatedly to
 * get elapsed times.
 * ------------------------------------
 * Modified by Clement Pernet
 * integrated into FFLAS_FFPACK
 *
 * ------------------------------------
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
 *
 * This file implements the C++ interface to commentators (for
 * providing runtime commentary to the user)
 */

#ifndef __TIMER_H
#define __TIMER_H

#include <iostream>

class BaseTimer {
    public:
	enum {
		MSPSEC = 1000000  // microsecond per second
	};

	// -- Clear timer :
	inline void clear() { _t = 0; }

	// -- total amount of second spent
	inline double time() const { return _t; }

	// -- Return a value to initialize random generator
	static long seed();

	// -- basic methods:
	std::ostream &print (std::ostream &) const;

	// -- Some arithmetic operators to compute cumulative time :
	BaseTimer& operator = (const BaseTimer & T) ;
	const BaseTimer operator - (const BaseTimer & T)  const;
	const BaseTimer operator - () ;
	const BaseTimer operator +  (const BaseTimer & T)  const;
	BaseTimer& operator += (const BaseTimer & T) { return *this = *this + T; };
	BaseTimer& operator -= (const BaseTimer & T) { return *this = *this - T; };

    public:
	double _start_t;  // time as of start ()
	double _t;        // time
};

inline std::ostream &operator << (std::ostream &o, const BaseTimer &BT)
	{ return BT.print(o); }

class RealTimer : public BaseTimer {
    public:
	inline RealTimer (const BaseTimer &BT) : BaseTimer (BT) {};
	inline RealTimer () {};
	void start ();
	void stop ();
};


class UserTimer : public BaseTimer {
    public:
	inline UserTimer (const BaseTimer &BT) : BaseTimer (BT) {};
	inline UserTimer () {};
	void start ();
	void stop ();
};


class SysTimer : public BaseTimer {
    public:
	inline SysTimer (const BaseTimer &BT): BaseTimer (BT) {};
	inline SysTimer () {};
	void start ();
	void stop ();
};


class Timer {
    public :

	// Clear timer :
	void clear();

	// Start timer
	void start ();

	// Stop timer
	void stop ();

	// total amount of second spent in user mode
	double usertime () const { return ut.time(); }

	// total amount of second spent in system mode
	double systime () const { return st.time(); }

	// real total amount of second spent.
	double realtime () const { return rt.time(); }

	// retourne une petite graine
	// long seed() const { return RealTimer::seed(); }

	// Some arithmetic operators to compute cumulative time :
	Timer& operator = (const Timer & T) ;
	const Timer operator - (const Timer & T)  const;
	const Timer operator - () ;
	const Timer operator + (const Timer & T)  const;
	Timer& operator += (const Timer & T) { return *this = *this + T; };
	Timer& operator -= (const Timer & T) { return *this = *this - T; };

	// -- methods :
	std::ostream &print (std::ostream &) const;

	size_t count() const {return _count;}

    private:
	size_t _count; // how many

	RealTimer rt;
	UserTimer ut;
	SysTimer  st;
};

// inline std::ostream &operator << (std::ostream &o, const Timer &T)
// 	{ return T.print (o); }

inline std::ostream &operator << (std::ostream &o, const Timer &T)
{
	double ut = T.usertime();
	if (ut < 0.0000000001) ut = 0;
	return o << T.realtime() << "s (" << ut << " cpu) [" << T.count() << "]";
}


#include <omp.h>
struct OMPTimer {
	double _c;
	void start() { _c = omp_get_wtime(); }
	void stop() { _c = omp_get_wtime() - _c; }
	void clear() { _c = 0.0; }
	double realtime() const { return _c; }
	double usertime() const { return _c; }
	OMPTimer& operator =(const OMPTimer& t) { _c = t._c; return *this; }
	OMPTimer& operator+=(const OMPTimer& t) { _c += t._c; return *this; }
	OMPTimer& operator-=(const OMPTimer& t) { _c -= t._c; return *this; }
	OMPTimer  operator +(const OMPTimer& t) const
	{
		OMPTimer r; r._c = _c + t._c; return r;
	}
	OMPTimer  operator -(const OMPTimer& t) const
	{
		OMPTimer r; r._c = _c - t._c; return r;
	}
	OMPTimer  operator -() { OMPTimer r; r._c = - _c; return r; }
};
//#endif
//
inline std::ostream &operator << (std::ostream &o, const OMPTimer &T)
{
	return o << T.usertime() << "s" ;
}



#include "timer.C"

#endif
