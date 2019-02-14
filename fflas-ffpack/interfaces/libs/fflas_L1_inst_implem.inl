/*
 * Copyright (C) 2014,2015 the FFLAS-FFPACK group
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *            Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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


namespace FFLAS {
    //---------------------------------------------------------------------
    // Level 1 routines
    //---------------------------------------------------------------------

    /** freduce
     * \f$x \gets  x mod F\f$.
     * @param F field
     * @param n size of the vectors
     * \param X vector in \p F
     * \param incX stride of \p X
     * @bug use cblas_(d)scal when possible
     */
    template INST_OR_DECL
    void
    freduce (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
             FFLAS_ELT* X, const size_t incX);

    /** freduce
     * \f$x \gets  y mod F\f$.
     * @param F field
     * @param n size of the vectors
     * \param Y vector of \p Element
     * \param incY stride of \p Y
     * \param X vector in \p F
     * \param incX stride of \p X
     * @bug use cblas_(d)scal when possible
     */
    template INST_OR_DECL
    void
    freduce (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
             const FFLAS_ELT* Y, const size_t incY,
             FFLAS_ELT* X, const size_t incX);

    /** finit
     * \f$x \gets  y mod F\f$.
     * @param F field
     * @param n size of the vectors
     * \param Y vector of \p OtherElement
     * \param incY stride of \p Y
     * \param X vector in \p F
     * \param incX stride of \p X
     * @bug use cblas_(d)scal when possible
     */
    template INST_OR_DECL
    void
    finit (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
           const FFLAS_ELT* Y, const size_t incY,
           FFLAS_ELT* X, const size_t incX);

    /** fconvert
     * \f$x \gets  y mod F\f$.
     * @param F field
     * @param n size of the vectors
     * \param Y vector of \p F
     * \param incY stride of \p Y
     * \param X vector in \p OtherElement
     * \param incX stride of \p X
     * @bug use cblas_(d)scal when possible
     */
    template INST_OR_DECL
    void
    fconvert (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
              FFLAS_ELT* X, const size_t incX,
              const FFLAS_ELT* Y, const size_t incY);
    // {
    // 	OtherElement_ptr Xi = X ;
    // 	const FFLAS_ELT* Yi = Y ;
    // 	for (; Xi < X+n*incX; Xi+=incX, Yi += incY )
    // 		F.convert( *Xi , *Yi);
    // }

    /** fnegin
     * \f$x \gets - x\f$.
     * @param F field
     * @param n size of the vectors
     * \param X vector in \p F
     * \param incX stride of \p X
     * @bug use cblas_(d)scal when possible
     */
    template INST_OR_DECL
    void
    fnegin (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
            FFLAS_ELT* X, const size_t incX);
    // {
    // 	FFLAS_ELT* Xi = X ;
    // 	for (; Xi < X+n*incX; Xi+=incX )
    // 		F.negin( *Xi );
    // }

    /** fneg
     * \f$x \gets - y\f$.
     * @param F field
     * @param n size of the vectors
     * \param X vector in \p F
     * \param incX stride of \p X
     * \param Y vector in \p F
     * \param incY stride of \p Y
     * @bug use cblas_(d)scal when possible
     */
    template INST_OR_DECL
    void
    fneg (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
          const FFLAS_ELT* Y, const size_t incY,
          FFLAS_ELT* X, const size_t incX);
    // {
    // 	FFLAS_ELT* Xi = X ;
    // 	const FFLAS_ELT* Yi = Y ;
    // 	for (; Xi < X+n*incX; Xi+=incX,Yi+=incY  )
    // 		F.neg( *Xi, *Yi );
    // }

    /** \brief fzero : \f$A \gets 0 \f$.
     * @param F field
     * @param n number of elements to zero
     * \param X vector in \p F
     * \param incX stride of \p X
     */
    template INST_OR_DECL
    void
    fzero (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
           FFLAS_ELT* X, const size_t incX);
    // {
    // 	if (incX == 1) { // contigous data
    // 		// memset(X,(int)F.zero,n); // might be bogus ?
    // 		for (size_t i = 0 ; i < n ; ++i)
    // 			F.assign(*(X+i), F.zero);

    // 	}
    // 	else { // not contiguous (strided)
    // 		for (size_t i = 0 ; i < n ; ++i)
    // 			F.assign(*(X+i*incX), F.zero);
    // 	}
    // }

    /** \brief fiszero : test \f$X = 0 \f$.
     * @param F field
     * @param n vector dimension
     * \param X vector in \p F
     * \param incX increment of \p X
     */
    template INST_OR_DECL
    bool
    fiszero (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
             const FFLAS_ELT* X, const size_t incX);
    // {
    // 	bool res=true;
    // 	for (size_t i = 0 ; i < n ; ++i)
    // 		res &= F.isZero (X [i*incX]);
    // 	return res;
    // }

    /** \brief fequal : test \f$X = Y \f$.
     * @param F field
     * @param n vector dimension
     * \param X vector in \p F
     * \param incX increment of \p X
     * \param Y vector in \p F
     * \param incY increment of \p Y
     */
    template INST_OR_DECL
    bool
    fequal (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n,
            const FFLAS_ELT* X, const size_t incX,
            const FFLAS_ELT* Y, const size_t incY);
    // {
    // 	bool res=true;
    // 	for (size_t i = 0 ; i < n ; ++i)
    // 		res &= F.areEqual (X [i*incX], Y [i*incY]);
    // 	return res;
    // }

    /** \brief fassign : \f$x \gets y \f$.
     * X is preallocated
     * @todo variant for triagular matrix
     * @param F field
     * @param N size of the vectors
     * \param [out] X vector in \p F
     * \param incX stride of \p X
     * \param [in] Y vector in \p F
     * \param incY stride of \p Y
     */
    template INST_OR_DECL
    void
    fassign (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t N,
             const FFLAS_ELT* Y, const size_t incY ,
             FFLAS_ELT* X, const size_t incX);


    /** fscalin
     * \f$x \gets \alpha \cdot x\f$.
     * @param F field
     * @param n size of the vectors
     * @param alpha scalar
     * \param X vector in \p F
     * \param incX stride of \p X
     * @bug use cblas_(d)scal when possible
     * @internal
     * @todo check if comparison with +/-1,0 is necessary.
     */
    template INST_OR_DECL
    void
    fscalin (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n, const FFLAS_ELT alpha,
             FFLAS_ELT* X, const size_t incX);


    /** fscal
     * \f$y \gets \alpha \cdot x\f$.
     * @param F field
     * @param n size of the vectors
     * @param alpha scalar
     * \param[in] X vector in \p F
     * \param incX stride of \p X
     * \param[out] Y vector in \p F
     * \param incY stride of \p Y
     * @bug use cblas_(d)scal when possible
     * @internal
     * @todo check if comparison with +/-1,0 is necessary.
     */
    template INST_OR_DECL
    void
    fscal (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t n
           , const FFLAS_ELT alpha
           , const FFLAS_ELT* X, const size_t incX
           , FFLAS_ELT* Y, const size_t incY);



    /** \brief faxpy : \f$y \gets \alpha \cdot x + y\f$.
     * @param F field
     * @param N size of the vectors
     * @param alpha scalar
     * \param[in] X vector in \p F
     * \param incX stride of \p X
     * \param[in,out] Y vector in \p F
     * \param incY stride of \p Y
     */
    template INST_OR_DECL
    void
    faxpy (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t N,
           const FFLAS_ELT alpha,
           const FFLAS_ELT* X, const size_t incX,
           FFLAS_ELT* Y, const size_t incY );

    /** \brief faxpby : \f$y \gets \alpha \cdot x + \beta \cdot y\f$.
     * @param F field
     * @param N size of the vectors
     * @param alpha scalar
     * \param[in] X vector in \p F
     * \param incX stride of \p X
     * \param beta scalar
     * \param[in,out] Y vector in \p F
     * \param incY stride of \p Y
     * \note this is a catlas function
     */
    // template INST_OR_DECL
    // void
    // faxpby (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t N,
    // 	const FFLAS_ELT alpha,
    // 	const FFLAS_ELT* X, const size_t incX,
    // 	const FFLAS_ELT beta,
    // 	FFLAS_ELT* Y, const size_t incY );


    /** \brief fdot: dot product \f$x^T  y\f$.
     * @param F field
     * @param N size of the vectors
     * \param X vector in \p F
     * \param incX stride of \p X
     * \param Y vector in \p F
     * \param incY stride of \p Y
     */
    template INST_OR_DECL
    FFLAS_ELT
    fdot (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t N,
          const FFLAS_ELT* X, const size_t incX,
          const FFLAS_ELT* Y, const size_t incY );

    /** \brief fswap: \f$ X \leftrightarrow Y\f$.
     * @bug use cblas_dswap when double
     * @param F field
     * @param N size of the vectors
     * \param X vector in \p F
     * \param incX stride of \p X
     * \param Y vector in \p F
     * \param incY stride of \p Y
     */
    template INST_OR_DECL
    void
    fswap (const FFLAS_FIELD<FFLAS_ELT>& F, const size_t N, FFLAS_ELT* X, const size_t incX,
           FFLAS_ELT* Y, const size_t incY );
    // {

    // 	FFLAS_ELT tmp; F.init(tmp);
    // 	FFLAS_ELT* Xi = X;
    // 	FFLAS_ELT* Yi=Y;
    // 	for (; Xi < X+N*incX; Xi+=incX, Yi+=incY ){
    // 		F.assign( tmp, *Xi );
    // 		F.assign( *Xi, *Yi );
    // 		F.assign( *Yi, tmp );
    // 	}
    // }

    template INST_OR_DECL
    void
    fadd (const FFLAS_FIELD<FFLAS_ELT>& F,  const size_t N,
          const FFLAS_ELT* A, const size_t inca,
          const FFLAS_ELT* B, const size_t incb,
          FFLAS_ELT* C, const size_t incc);

    template INST_OR_DECL
    void
    fsub (const FFLAS_FIELD<FFLAS_ELT>& F,  const size_t N,
          const FFLAS_ELT* A, const size_t inca,
          const FFLAS_ELT* B, const size_t incb,
          FFLAS_ELT* C, const size_t incc);

    template INST_OR_DECL
    void
    faddin (const FFLAS_FIELD<FFLAS_ELT>& F,  const size_t N,
            const FFLAS_ELT* B, const size_t incb,
            FFLAS_ELT* C, const size_t incc);

    // template INST_OR_DECL
    // void
    // fsubin (const FFLAS_FIELD<FFLAS_ELT>& F,  const size_t N,
    // 	FFLAS_ELT* C, const size_t incc);


    template INST_OR_DECL
    void
    fadd (const FFLAS_FIELD<FFLAS_ELT>& F,  const size_t N,
          const FFLAS_ELT* A, const size_t inca,
          const FFLAS_ELT alpha,
          const FFLAS_ELT* B, const size_t incb,
          FFLAS_ELT* C, const size_t incc);

} // FFLAS

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
