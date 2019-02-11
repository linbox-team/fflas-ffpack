/* utils/debug.h
 *
 * Copyright (C) 2011 Fflas-ffpack
 * Modified by BB, from LinBox
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

/*! @file utils/debug.h
 * @ingroup util
 * Various utilities for debugging.
 * @todo we should put vector printing elsewhere.
 */

#ifndef __FFLASFFPACK_util_debug_H
#define __FFLASFFPACK_util_debug_H
#include <fflas-ffpack/fflas-ffpack-config.h>

#include <iostream>
#include <sstream>
#include <cmath>

#include "fflas-ffpack/fflas-ffpack-config.h"


#ifdef __FFLASFFPACK_HAVE_STDINT_H
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include <limits>

// If somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)
#ifndef INT64_MAX
#define INT64_MAX std::numeric_limits<int64_t>::max()
#endif

#ifndef UINT64_MAX
#define UINT64_MAX std::numeric_limits<uint64_t>::max()
#endif

#ifndef INT32_MAX
#define INT32_MAX std::numeric_limits<int32_t>::max()
#endif

#ifndef UINT32_MAX
#define UINT32_MAX std::numeric_limits<uint32_t>::max()
#endif

#ifndef INT16_MAX
#define INT16_MAX std::numeric_limits<int16_t>::max()
#endif

#ifndef UINT16_MAX
#define UINT16_MAX std::numeric_limits<uint16_t>::max()
#endif

#ifndef INT8_MAX
#define INT8_MAX std::numeric_limits<int8_t>::max()
#endif

#ifndef UINT8_MAX
#define UINT8_MAX std::numeric_limits<uint8_t>::max()
#endif

#else
#error "you need intXX_t types"
#endif

#ifndef NDEBUG
#include <stdexcept>
#define FFLASFFPACK_check(check) \
if (!(check)) {\
    FFPACK::failure()(__func__, __FILE__, __LINE__, #check); \
    throw std::runtime_error(#check); \
}
#define FFLASFFPACK_abort(msg) \
{\
    FFPACK::failure()(__func__, __FILE__, __LINE__, msg); \
    throw std::runtime_error(msg); \
}
#else
#define FFLASFFPACK_check(check) ((void) 0)
#define FFLASFFPACK_abort(mst) ((void) 0)
#endif



namespace FFPACK {

    /*!  A precondtion failed.
     * @ingroup util
     * The \c throw mechanism is usually used here as in
     \code
     if (!check)
     failure()(__func__,__LINE__,"this check just failed");
     \endcode
     * The parameters of the constructor help debugging.
     */
    class Failure {
    protected:
        std::ostream *_errorStream;

    public:

        Failure() {}

        /*! @internal
         * A precondtion failed.
         * @param function usually \c __func__, the function that threw the error
         * @param line     usually \c __LINE__, the line where it happened
         * @param check    a string telling what failed.
         */
        void operator() (const char *function, int line, const char *check)
        {
            if (_errorStream == (std::ostream *) 0)
                _errorStream = &std::cerr;

            (*_errorStream) << std::endl << std::endl;
            (*_errorStream) << "ERROR (" << function << ":" << line << "): ";
            (*_errorStream) << "Precondition not met:" << check << std::endl;
        }

        /*! @internal
         * A precondtion failed.
         * The parameter help debugging. This is not much different from the previous
         * except we can digg faster in the file where the exception was triggered.
         * @param function usually \c __func__, the function that threw the error
         * @param file     usually \c __FILE__, the file where this function is
         * @param line     usually \c __LINE__, the line where it happened
         * @param check    a string telling what failed.
         */
        void operator() (const char* function, const char *file, int line, const char *check)
        {
            if (_errorStream == (std::ostream *) 0)
                _errorStream = &std::cerr;

            (*_errorStream) << std::endl << std::endl;
            (*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
            (*_errorStream) << "Precondition not met:" << check << std::endl;
        }

        void setErrorStream (std::ostream &stream);

        /*! @internal overload the virtual print of LinboxError.
         * @param o output stream
         */
        std::ostream &print (std::ostream &o) const
        {
            if (std::ostringstream * str = dynamic_cast<std::ostringstream*>(_errorStream))
                return o << str->str() ;
            else
                throw "FFLAS-FFPACK ERROR: Failure exception is not initialized correctly";
        }
    };


    inline Failure& failure() {
        static Failure failure_internal;
        return failure_internal;
    }

    template<class T>
    inline bool isOdd (const T & a) {
        return (a%2);
    }

    inline bool isOdd(const float &a) {
        return (bool)(int)fmodf(a,2.f);
    }

    inline bool isOdd(const double &a) {
        return (bool)(int)fmod(a,2.);
    }

} // FFPACK

#endif // __FFLASFFPACK_util_debug_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
