#ifndef __FFLASFFPACK_ffpack_charpoly_kglu_INL
#define __FFLASFFPACK_ffpack_charpoly_kglu_INL

/* ffpack/ffpack_charpoly_kglu.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
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

namespace FFPACK {

    namespace Protected {

        template<class Field>
        size_t updateD(const Field& F, size_t * d, size_t k,
                       std::vector<std::vector<typename Field::Element> >& minpt)
        {
            size_t ind=0, i=0;
            while(i<k){
                if (d[i]){
                    if (ind<i){
                        d[ind] = d[i];
                        minpt[ind++] = minpt[i];
                    }
                    else ind++;
                }
                i++;
            }
            for (i=ind; i<k; ++i)
                minpt[i].resize(0);
            minpt.resize(ind);
            return ind;
        }

        // Subroutine for Keller-Gehrig charpoly algorithm
        // Compute the new d after a LSP ( d[i] can be zero )
        template<class Field>
        size_t newD( const Field& F, size_t * d, bool& KeepOn,
                     const size_t l, const size_t N,
                     typename Field::Element_ptr X,
                     const size_t * Q,
                     std::vector<std::vector<typename Field::Element> >& minpt)
        {
            typedef typename Field::Element elt;
            //const elt * Xi = X; // Xi points to the begining of each block
            elt *Li=X, *Xminp=X;
            KeepOn = false;
            size_t  i, jtot=0, dtot = 0, nrtot=0;

            for ( i=0; dtot<N; ++i){ // for each block
                size_t j = 0;
                size_t s ;
                size_t nr = s = ( d[i]==l )? 2*l : d[i];
                if (s > N-dtot)
                    s= N-dtot;
                nrtot += nr;

                while ( (Q[j+jtot] <nrtot) && (j+jtot<N) )
                    j++;

                Xminp = X+Q[j+jtot-1]*N+jtot+j ;
                d[i] = j;
                jtot+=j;
                dtot += j;

                if (j<nr){
                    minpt[i].resize(j);
                    ftrsv( F, FFLAS::FflasLower, FFLAS::FflasTrans, FFLAS::FflasUnit,
                           j, Li, N, Xminp-j+N,1);
                    elt* Xi = Xminp-j+N;
                    for (size_t ii = 0; ii<j; ++ii, ++Xi)
                        minpt[i][ii] = *Xi;
                }
                Li += nr*N+j;
                if ( j==2*l )
                    KeepOn = true;
            } // Invariant: sum(d[i],i=0..k-1) = n

            return i;
        }


        //---------------------------------------------------------------------
        // CharPoly: Compute the characteristic polynomial of A using
        // Keller-Gehrig's algorithm
        //---------------------------------------------------------------------
        template <class Field, class Polynomial>
        std::list<Polynomial>&
        KellerGehrig( const Field& F, std::list<Polynomial>& charp, const size_t N,
                      typename Field::ConstElement_ptr A, const size_t lda )
        {


            typename Field::ConstElement_ptr Ai = A;
            typename Field::Element_ptr U = FFLAS::fflas_new (F, N, N);     // to store A^2^i
            typename Field::Element_ptr B = FFLAS::fflas_new (F, N, N);     // to store A^2^i
            typename Field::Element_ptr V = FFLAS::fflas_new (F, N, N);     // to store A^2^i.U
            typename Field::Element_ptr X = FFLAS::fflas_new (F, 2*N, N);   // to compute the LSP factorization
            typename Field::Element_ptr Ui, Uj, Uk, Ukp1, Ukp1new, Bi, Vi, Vk, Xi=X, Xj;
            size_t * P = FFLAS::fflas_new<size_t>(N); // Column Permutation for LQUP
            size_t * Q = FFLAS::fflas_new<size_t>(2*N); // Row Permutation for LQUP

            size_t * d= FFLAS::fflas_new<size_t>(N);   // dimensions of Vect(ei, Aei...)
            size_t * dv = FFLAS::fflas_new<size_t>(N);
            size_t * dold = FFLAS::fflas_new<size_t>(N); // copy of d
            // vector of the opposite of the coefficient of computed minpolys
            std::vector< std::vector< typename Field::Element > > m(N);
            typename Polynomial::iterator it;
            size_t i=0, l=1, j, k=N,  cpt, newRowNb;
            bool  KeepOn;

            for ( i=0; i<N; ++i)
                dv[i] = dold[i] = d[i] = 1;

            // Computing the first X: (e1; e1A^t; e2; e2A^t;...;en;enA^t)
            for ( i=0, Ui=U, Vi=V, Bi=B; i<N; ++i, Ai -= N*lda-1  ){
                for ( Xj=Xi, Uj=Ui; Xj<Xi+N; ++Xj, ++Uj){
                    F.assign(*Xj, F.zero);
                    F.assign(*Ui, F.zero);
                }
                F.assign(*(Ui+i), F.one);
                F.assign(*(Xi+i), F.one);
                while ( Xj<Xi+2*N) {
                    *(Bi++) = *(Xj++) = *(Vi++) = *Ai;
                    Ai+=lda;
                }
                Xi = Xj;
                Ui = Uj;
            }

            // step form elimination using LQUP factorization
            for ( i=0;i<N;++i)
                P[i]=0;
            for ( i=0;i<2*N;++i)
                Q[i]=0;
            LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, 2*N, N, X, N, P, Q);

            k = Protected::newD( F,d, KeepOn, l, N, X, Q, m);

            while(KeepOn){ // Main loop, until each subspace dimension has been found
                size_t nrowX, ind ;
                // Updating U:
                Uk = U;
                // Firstly, removing extra rows
                for ( i = 0, cpt = 0; i<N; ++i){
                    if (d[i] < dold[i]){
                        Ukp1new = Uk + d[i]*N;  // new position of Uk+1
                        Ukp1 = Uk + dold[i]*N; // first row of Uk+1
                        Ui = Ukp1new;
                        Uj = Ukp1;
                        while ( Uj < U + N*N ){
                            for ( j=N; j; --j)
                                *(Ui++) = *(Uj++);
                        }
                        Uk = Ukp1new;
                        dold[i] = d[i];
                    }
                    else {
                        Uk += dold[i]*N;
                    }
                    cpt += d[i];
                }

                // Then inserting the duplicated blocks
                Uk = U;
                Vk = V;
                for ( i = 0; i<k; ++i){
                    Ukp1 = Uk + dold[i]*N; // first row of Uk+1
                    newRowNb = d[i] - dold[i];
                    if ( newRowNb > 0){
                        Ui = U+N*N-1; // last row of U
                        Uj = U+(N-newRowNb)*N-1; // last row future pos
                        while ( Uj > Ukp1-1){
                            for ( j=N;j;--j)
                                *(Ui--) = *(Uj--);// moving block
                        }
                        Uj++;
                        Vi = Vk;
                        while ( Uj < Ukp1+N*newRowNb ){
                            for ( j=N;j;--j)
                                *(Uj++) = *(Vi++);
                        }
                    }
                    Uk = Uk + d[i]*N;
                    Vk += dv[i]*N;
                }

                // max block size of X, U, V is l=2^i
                l*=2;
                // B = A^2^i
                fsquare( F, FFLAS::FflasNoTrans, N, F.one, B, N, F.zero, B, N );
                // V = U.B^t
                fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, N, N, N, F.one,
                       U, N, B, N, F.zero, V, N);
                // X = ( U1, V1, U2, V2, ... )
                Xi = X; Ui = U; Vi = V;
                ind=0; cpt=0; nrowX = 0;
                while ( Vi < V + N*N ){
                    // Copying Uk
                    for ( i = d[ind]; i; --i, nrowX++){
                        for ( j = N; j; --j )
                            *(Xi++) = *(Ui++);
                    }
                    // Copying Vk
                    if ( d[ind] == l || ind==k-1 ){
                        cpt+=2*d[ind];
                        for ( i=d[ind]; i; i--, nrowX++)
                            for ( j=N; j; j--)
                                *(Xi++) = *(Vi++);
                    }
                    else{
                        cpt += d[ind];
                        Vi = Vi + N*d[ind];
                    }
                    ind++;
                }
                // removes factors of degree 0 in m
                k = Protected::updateD( F, d, k, m);

                for (i=0;i<k;++i)
                    dv[i] = dold[i] = d[i];

                // step form elimination of X using LSP
                for ( i=0;i<N;++i)
                    P[i]=0;
                for ( i=0;i<2*N;++i)
                    Q[i]=0;
                LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, nrowX, N, X, N, P, Q);

                // Recompute the degrees of the list factors
                k = Protected::newD(F, d, KeepOn, l, N, X,Q, m);
            }
            FFLAS::fflas_delete (U);
            FFLAS::fflas_delete (V);
            FFLAS::fflas_delete (B);
            FFLAS::fflas_delete( P);
            FFLAS::fflas_delete( Q);
            FFLAS::fflas_delete( dv);
            FFLAS::fflas_delete( dold);

            k = Protected::updateD( F, d, k, m);
            // Constructing the CharPoly
            for ( i=0; i<k; ++i){
                Polynomial * minP = new Polynomial(d[i]+1);
                minP->operator[](d[i]) = F.one;
                it = minP->begin();
                for ( j=0; j<d[i]; ++j, it++)
                    F.neg(*it, m[i][j]);
                charp.push_back( *minP );
            }
            FFLAS::fflas_delete (X);
            FFLAS::fflas_delete( d);
            return charp;
        }

    } // Protected
} // FFPACK

#endif // __FFLASFFPACK_ffpack_charpoly_kglu_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
