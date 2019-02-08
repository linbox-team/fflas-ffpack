/*
 * Copyright (C) 2014 the FFLAS-FFPACK group
 *
 * Written by   Bastien Vialla<bastien.vialla@lirmm.fr>
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

// functions wih _r are relaxed, meaning no modular reduction

template <class _Field> class FieldSimd {
public:
    using Field = _Field;
    using Element = typename Field::Element;
    using simd = Simd<typename _Field::Element>;
    using vect_t = typename simd::vect_t;
    using scalar_t = typename simd::scalar_t;

    static const constexpr size_t vect_size = simd::vect_size;

    static const constexpr size_t alignment = simd::alignment;

private:
    using Self = FieldSimd<Field>;

    const Field *_field;
    vect_t _modulus;
    vect_t _invmod;
    vect_t _negmod;
    vect_t _mask;
    vect_t _min;
    vect_t _max;

public:
    FieldSimd(const Field &f) : _field(&f) { init(); }

private:
    void init() {
        _modulus = simd::set1((Element)_field->characteristic());
        _min = simd::set1(_field->minElement());
        _max = simd::set1(_field->maxElement());
        _negmod = simd::set1(-(Element)_field->characteristic());
        if (std::is_floating_point<Element>::value) {
            _invmod = simd::set1(1 / ((Element)_field->characteristic()));
        }
    }

public:
    FieldSimd(const Self &) = default;
    FieldSimd(Self &&) = default;

    Self &operator=(const Self &) = default;
    Self &operator=(Self &&) = default;

    INLINE vect_t init(vect_t &x, const vect_t a) const { return x = mod(a); }

    INLINE vect_t init(const vect_t a) const { return mod(a); }

    INLINE vect_t add(vect_t &c, const vect_t a, const vect_t b) {
        c = simd::add(a, b);
        _mask = simd::greater(c, _max);
        _mask = simd::vand(_mask, _modulus);
        return c = simd::sub(c, _mask);
    }

    INLINE vect_t add(const vect_t a, const vect_t b) {
        vect_t c;
        c = simd::add(a, b);
        _mask = simd::greater(c, _max);
        _mask = simd::vand(_mask, _modulus);
        return c = simd::sub(c, _mask);
    }

    INLINE vect_t addin(vect_t &a, const vect_t b) const { return a = add(a, b); }

    INLINE vect_t add_r(vect_t &c, const vect_t a, const vect_t b) const { return c = simd::add(a, b); }

    INLINE vect_t add_r(const vect_t a, const vect_t b) const { return simd::add(a, b); }

    INLINE vect_t addin_r(vect_t &a, const vect_t b) const { return a = add_r(a, b); }

    INLINE vect_t sub(vect_t &c, const vect_t a, const vect_t b) {
        c = simd::sub(a, b);
        _mask = simd::lesser(c, _min);
        _mask = simd::vand(_mask, _modulus);
        return c = simd::add(c, _mask);
    }

    INLINE vect_t sub(const vect_t a, const vect_t b) {
        vect_t c;
        c = simd::sub(a, b);
        _mask = simd::greater(c, _max);
        _mask = simd::vand(_mask, _modulus);
        return c = simd::add(c, _mask);
    }

    INLINE vect_t subin(vect_t &a, const vect_t b) const { return a = sub(a, b); }

    INLINE vect_t sub_r(vect_t &c, const vect_t a, const vect_t b) const { return c = simd::sub(a, b); }

    INLINE vect_t sub_r(const vect_t a, const vect_t b) const { return simd::sub(a, b); }

    INLINE vect_t subin_r(vect_t &a, const vect_t b) const { return a = sub_r(a, b); }

    INLINE vect_t zero(vect_t &x) const { return x = simd::zero(); }

    INLINE vect_t zero() const { return simd::zero(); }

    INLINE vect_t mod(vect_t &c) const {
        if (std::is_floating_point<Element>::value) {
            vect_t q, t;
            q = simd::mul(c, _invmod);
            q = simd::floor(q);
            c = simd::fnmadd(c, q, _modulus);
            q = simd::greater(c, _max);
            t = simd::lesser(c, _min);
            q = simd::vand(q, _negmod);
            t = simd::vand(t, _modulus);
            q = simd::vor(q, t);
            return c = simd::add(c, q);
        } else {
            FFLASFFPACK_abort("pas implementé");
        }
    }

    INLINE vect_t mul(vect_t &c, const vect_t a, const vect_t b) const { return c = mod(simd::mul(a, b)); }

    INLINE vect_t mul(const vect_t a, const vect_t b) const { return mod(simd::mul(a, b)); }

    INLINE vect_t mulin(vect_t &a, const vect_t b) const { return mul(a, a, b); }

    INLINE vect_t mul_r(vect_t &c, const vect_t a, const vect_t b) const { return c = simd::mul(a, b); }

    INLINE vect_t mul_r(const vect_t a, const vect_t b) const { return simd::mul(a, b); }

    INLINE vect_t axpy(vect_t &r, const vect_t a, const vect_t b, const vect_t c) const {
        return r = mod(simd::fmadd(c, a, b));
    }

    INLINE vect_t axpy(const vect_t c, const vect_t a, const vect_t b) const { return mod(simd::fmadd(c, a, b)); }

    INLINE vect_t axpyin(vect_t &c, const vect_t a, const vect_t b) const { return c = axpy(c, a, b); }

    INLINE vect_t axpy_r(vect_t &r, const vect_t a, const vect_t b, const vect_t c) const {
        return r = simd::fmadd(c, a, b);
    }

    INLINE vect_t axpy_r(const vect_t c, const vect_t a, const vect_t b) const { return simd::fmadd(c, a, b); }

    INLINE vect_t axpyin_r(vect_t &c, const vect_t a, const vect_t b) const { return c = axpy_r(c, a, b); }

    INLINE vect_t maxpy(vect_t &r, const vect_t a, const vect_t b, const vect_t c) const {
        return r = mod(simd::fmsub(c, a, b));
    }

    INLINE vect_t maxpy(const vect_t c, const vect_t a, const vect_t b) const { return mod(simd::fmsub(c, a, b)); }

    INLINE vect_t maxpyin(vect_t &c, const vect_t a, const vect_t b) const { return c = maxpy(c, a, b); }
};
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
