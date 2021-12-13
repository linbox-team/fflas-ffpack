/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by <bastien.vialla@lirmm.fr>
 *
 *  STL align allocator inspired by MAlloc from Stephan T. Lavavej, Visual C++ Libraries Developer
 *  (http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx)
 *  Update to c++11
 *
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

#ifndef __FFLASFFPACK_align_allocator_H
#define __FFLASFFPACK_align_allocator_H

#include "fflas-ffpack/config.h"

#ifdef __FFLASFFPACK_HAVE_CXX11

#include <memory>
#include <utility>
#include <assert.h>
#include <cstddef>
#include <iostream>

#include "fflas-ffpack/utils/fflas_intrinsic.h"
//#include <immintrin.h>
// Alignment Type
enum class Alignment : size_t {
    NONE = 0,
    Normal = sizeof(void*),
    SSE = 16,
    AVX = 32,
    XEON_PHI = 64,
    CACHE_LINE = 64,
    CACHE_PAGESIZE = 4096,
    DEFAULT =
#ifdef __FFLASFFPACK_HAVE_AVX512F_INSTRUCTIONS
    64
#elif defined __FFLASFFPACK_HAVE_AVX_INSTRUCTIONS
    32
#else
    16
#endif
};

/*
 * Allocate T[size] with address aligned to alignement
 * ex : int* tab = malloc_align<int>(100, Alignment::AVX)
 */
template<class T>
T* malloc_align(size_t size, Alignment alignment = Alignment::DEFAULT) noexcept
{
    void* p = nullptr;
    int err = 0;
    err = posix_memalign(&p, (size_t) alignment, size*sizeof(T));
    if(err)
        std::cout << "posix_memalign error" << std::endl;
    //return new(p) T[size];
    return static_cast<T*>(p);
}

namespace detail
{
    inline void* allocate(size_t align, size_t size)
    {
        assert(align >= sizeof(void*));

        if (size == 0) {
            return nullptr;
        }

        void* ptr = nullptr;
        int rc = posix_memalign(&ptr, align, size);

        if (rc != 0) {
            return nullptr;
        }

        return ptr;
    }

    inline void deallocate(void* ptr) noexcept
    {
        return free(ptr);
    }
}


/* STL Aligned Allocator
 * ex : std::vector<T, AlignedAllocator<T, Alignment::AVX>>;
 *
 * template<class T> using vector = std::vector<T, AlignedAllocator<T, Alignment::AVX>>;
 */

template <class T, Alignment Align = Alignment::SSE> class AlignedAllocator;

template <class T, Alignment Align> class AlignedAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = ptrdiff_t;

    using propagate_on_container_move_assignment = std::true_type;

    template <class U> struct rebind {
        using other = AlignedAllocator<U, Align>;
    };

public:
    AlignedAllocator() noexcept {}

    template <class U>
    AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept {}

    size_type max_size() const noexcept {
        return (size_type(~0) - size_type(Align)) / sizeof(T);
    }

    pointer address(reference x) const noexcept { return std::addressof(x); }

    const_pointer address(const_reference x) const noexcept {
        return std::addressof(x);
    }

    pointer allocate(size_type n,
                     typename AlignedAllocator<void, Align>::const_pointer = 0) {
        const size_type alignment = static_cast<size_type>(Align);
        void* ptr = detail::allocate(alignment, n * sizeof(T));
        if (ptr == nullptr) {
            throw std::bad_alloc();
        }

        return reinterpret_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type) noexcept {
        return detail::deallocate(p);
    }

    template <class U, class... Args> void construct(U* p, Args&&... args) {
        ::new (reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...);
    }

    void destroy(pointer p) { p->~T(); }
};

/*
 * Specialization for void*
 */
template <Alignment Align> class AlignedAllocator<void, Align> {
public:
    using pointer = void*;
    using const_pointer = const void*;
    using value_type = void;

    template <class U> struct rebind {
        using other = AlignedAllocator<U, Align>;
    };
};

/*
 *  Specialization for const T
 */
template <class T, Alignment Align> class AlignedAllocator<const T, Align> {
public:
    using value_type = T;
    using pointer = const T*;
    using const_pointer = const T*;
    using reference = const T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = ptrdiff_t;

    using propagate_on_container_move_assignment = std::true_type;

    template <class U> struct rebind {
        using  other = AlignedAllocator<U, Align>;
    };

public:
    AlignedAllocator() noexcept {}

    template <class U>
    AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept {}

    size_type max_size() const noexcept {
        return (size_type(~0) - size_type(Align)) / sizeof(T);
    }

    const_pointer address(const_reference x) const noexcept {
        return std::addressof(x);
    }

    pointer allocate(size_type n,
                     typename AlignedAllocator<void, Align>::const_pointer = 0) {
        const size_type alignment = static_cast<size_type>(Align);
        void* ptr = detail::allocate(alignment, n * sizeof(T));
        if (ptr == nullptr) {
            throw std::bad_alloc();
        }

        return reinterpret_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type) noexcept {
        return detail::deallocate(p);
    }

    template <class U, class... Args> void construct(U* p, Args&&... args) {
        ::new (reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...);
    }

    void destroy(pointer p) { p->~T(); }
};

template <class T, Alignment TAlign, class U, Alignment UAlign>
inline bool operator==(const AlignedAllocator<T, TAlign>&,
                       const AlignedAllocator<U, UAlign>&) noexcept {
    return TAlign == UAlign;
}

template <class T, Alignment TAlign, class U, Alignment UAlign>
inline bool operator!=(const AlignedAllocator<T, TAlign>&,
                       const AlignedAllocator<U, UAlign>&) noexcept {
    return TAlign != UAlign;
}

#else // C++11
#error "You need a c++11 compiler."
#endif // C++11

#endif /* _FFLASFFPACK_align_allocator_h */
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
