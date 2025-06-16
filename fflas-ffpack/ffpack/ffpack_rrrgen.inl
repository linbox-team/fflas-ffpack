#ifndef __FFLASFFPACK_ffpack_rrrgen_inl
#define __FFLASFFPACK_ffpack_rrrgen_inl

#include <iostream>

namespace FFPACK{

/// @brief  RRgen Class for easier representation and use. 
///         A = PxLxUxQ. A (nxm) rank r
template<class Field>
class RRgen {
public:
    size_t n;                           // size of lines of original matrix 
    size_t m;                           // size of columns of original matrix
    size_t r;                           // rank of the original matrix
    typename Field::Element_ptr PL;     // PL from PLUQ (n,r)
    size_t ldPL;                        
    typename Field::Element_ptr UQ;     // UQ from PLUQ (r*m)
    size_t ldUQ;                        
    

    // RRgen(const Field& Fi, size_t n, size_t m, size_t r,
    //         typename Field::Element_ptr PL, size_t ldPL, 
    //         typename Field::Element_ptr UQ,size_t ldUQ)
    //     : n(n), m(m), r(r), PL(PL), ldPL(ldPL), UQ(UQ), ldUQ(ldUQ) {}
    RRgen(const Field& Fi, size_t n_A, size_t m_A, size_t r_A,
        typename Field::Element_ptr PL_A, size_t ldPL_A, 
        typename Field::Element_ptr UQ_A,size_t ldUQ_A){
        n = n_A;
        m = m_A;
        r = r_A;
        PL = FFLAS::fflas_new(Fi, n_A, r_A);
        FFLAS::fassign(Fi,n_A,r_A,PL_A,ldPL_A,PL,r_A);
        if (UQ_A){
            UQ = FFLAS::fflas_new(Fi, r_A, m_A);
            FFLAS::fassign(Fi,r_A,m_A,UQ_A,ldUQ_A,UQ,m_A);
        }
        else {
            UQ = UQ_A;
        }
        ldUQ = m_A;
        ldPL = r_A;
    }
    
    RRgen(const Field& Fi, size_t n_A, size_t m_A, size_t r_A,
        typename Field::Element_ptr PL_A, size_t ldPL_A, 
        typename Field::Element_ptr UQ_A,size_t ldUQ_A, bool noCopy){
        
        if (noCopy){
            n = n_A;
            m = m_A;
            r = r_A;
            ldPL = ldPL_A;
            ldUQ = ldUQ_A;
            PL = PL_A;
            UQ = UQ_A;
        }
        else{
            RRgen(Fi,n_A, m_A, r_A, PL_A, ldPL_A, UQ_A, ldUQ_A);
        }
    }
    
    
    RRgen(const Field& Fi, size_t n_A, size_t m_A, typename Field::Element_ptr A, size_t ldA){
        n = n_A;
        m = m_A;
        size_t* P = FFLAS::fflas_new<size_t>(n);
        size_t* Q = FFLAS::fflas_new<size_t>(m);

        
        typename Field::Element_ptr A_copy = FFLAS::fflas_new(Fi, n_A, m_A);
        FFLAS::fassign(Fi, n, m, A, ldA, A_copy, m_A);

        r = PLUQ(Fi, FFLAS::FflasNonUnit, n, m, A_copy, m_A, P, Q);


        PL = FFLAS::fflas_new(Fi, n, r);
        UQ = FFLAS::fflas_new(Fi, r, m);
        ldPL = r;
        ldUQ = m;

        // extraction of U
        getTriangular<Field>(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit, n, m, r, A_copy, m_A, UQ, m, true);

        // extraction of L
        getTriangular<Field>(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, n, m, r, A_copy, m_A, PL, r, true);

        // apply the permutations (P on L and Q on U)
        applyP<Field>(Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, r, 0, n, PL, r, P);
        

        applyP<Field>(Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, r, 0, m, UQ, m, Q);

        FFLAS::fflas_delete(A_copy);
        FFLAS::fflas_delete(P,Q);
    }
    
    void RRExpand( const Field& Fi, typename Field::Element_ptr A, size_t ldA){
                
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n, m, r, Fi.one,
            PL, ldPL,
            UQ, ldUQ,
            Fi.zero, A, ldA);
    }

    ~RRgen() {
        FFLAS::fflas_delete(PL);
        if (UQ){
            FFLAS::fflas_delete(UQ);
        }
    }
};

/// @brief  RRRgen Class for the tree representation.                                                               A = [A11  A12]
///         If the node is a leaf the matrix is stored in LU_right->PL and LU_left, left and right are None.            [A21  A22]
template<class Field>
class RRRgen {

public:
    RRgen<Field>* LU_right;          // PLUQ representation of the top right submatrix (N2*rr)*(rr*N1) = A12
    RRgen<Field>* LU_left;           // PLUQ representation of the down left submatrix (N1*rl)*(rl*N2) = A21        
    size_t size_N1;                  // size of N/2 (N is the size of the matrix represented by the RRRgen)
    size_t size_N2;                  // size of N-N1
    size_t t;                        // threshold for recursive representation : when N1+N2 < t, the matrix is stored completly
    RRRgen<Field>* left;             // recursively pointing on the same representation of the up left submatrix = A11
    RRRgen<Field>* right;            // recursively pointing on the same representation of the down right submatrix = A22


    RRRgen(RRgen<Field>* LU_right, RRgen<Field>* LU_left, size_t size_N1, size_t size_N2,size_t t, RRRgen* left, RRRgen* right)
        :  LU_right(LU_right), LU_left(LU_left), size_N1(size_N1), size_N2(size_N2),t(t), left(left), right(right){}
    
    
    RRRgen(   const Field& F,typename Field::Element_ptr leaf, size_t n,size_t t){

        LU_right = new RRgen<Field>(F,n,n,0,leaf,n,nullptr,0);
        LU_left = nullptr;
        size_N1 = n;
        size_N2 = 0;
        t = t;
        left = nullptr;
        right = nullptr; 
    }

    RRRgen( const Field& Fi,const size_t N, const size_t threshold,typename Field::Element_ptr A, const size_t lda)
    {
        if (N <= threshold) { 
            // leaf
            LU_right = new RRgen<Field>(Fi,N,N,N,A,lda,nullptr,0);
            LU_left = nullptr;
            size_N1 = N;
            size_N2 = 0;
            t = threshold;
            left = nullptr;
            right = nullptr;
            return;
        }

        size_N1 = N/2;
        size_N2 = N - size_N1;

        // cut the matrix A in four parts
        // in case of odd N, A12 and A21 are rectangles and not squares
        typename Field::Element_ptr A11 = A;
        typename Field::Element_ptr A12 = A + size_N1;
        typename Field::Element_ptr A21 = A + lda*size_N1;
        typename Field::Element_ptr A22 = A21 + size_N1;

        ///////// PLUQ Factorisation for A12 and A21

        LU_right = new RRgen(Fi, size_N1,  size_N2, A12, lda);

        LU_left = new RRgen(Fi, size_N2,  size_N1, A21, lda);


        // recursion on A11 and A22   
        left = new RRRgen( Fi, size_N1, t, A11, lda);

        right = new RRRgen( Fi, size_N2, t, A22, lda);
        return;
    }

    ~RRRgen() {
        delete LU_right;
        
        if (LU_left) {
            delete (LU_left);
        }

        if (left) {
            delete left;
        }

        if (right) {
            delete right;
        }
    }
};


/// @brief (algo 4) Compute the dense matrix of RRR(A) in B recursive part
/// @tparam Field 
/// @param Fi 
/// @param A 
/// @param B 
/// @param ldb 
template<class Field>
inline void RRRExpand (const Field& Fi,
            const RRRgen<Field>* A,
            typename Field::Element_ptr B, const size_t ldb)
    {
        if (A->left == nullptr){
            size_t N = A->size_N1;
            FFLAS::fassign(Fi, N, N, A->LU_right->PL, A->LU_right->ldPL, B, ldb);
            return;
        }

        size_t N1 = A->size_N1;
        size_t N2 = A->size_N2;

        // cut the matrix A in four parts
        // in case of odd N, A12 and A21 are rectangles and not squares
        typename Field::Element_ptr B11 = B;
        typename Field::Element_ptr B12 = B + N1;
        typename Field::Element_ptr B21 = B + ldb*N1;
        typename Field::Element_ptr B22 = B21 + N1;

        
        // B11 < RRRExpand(A11)
        RRRExpand<Field>(Fi, A->left, B11, ldb);

        // B22 < RRRExpand(A22)
        RRRExpand<Field>(Fi, A->right, B22, ldb);

        // B12 < L_u * U_u
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N1, N2, A->LU_right->r, Fi.one,
            A->LU_right->PL, A->LU_right->r,
            A->LU_right->UQ, N2,
            Fi.zero, B12, ldb);

        
        // B21 < L_l * U_l
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N2, N1, A->LU_left->r, Fi.one,
            A->LU_left->PL, A->LU_left->r,
            A->LU_left->UQ, N1,
            Fi.zero, B21, ldb);

        
    }



/// @brief (algo 5) multiplies two matrices stored as rank revealing factorization. C = minus ? -1 : 1 A*B
/// @tparam Field 
/// @param Fi 
/// @param A stored with an RRgen
/// @param B stored with an RRgen
template<class Field>
inline RRgen<Field>* RRxRR (const Field& Fi,RRgen<Field>* A, RRgen<Field>*B, bool minus)
    {   
        size_t m = B->m;
        size_t n = A->n;
        size_t k = A->m;
        
        // X < UA * LB
        typename Field::Element_ptr X = FFLAS::fflas_new(Fi, A->r, B->r);
        fgemm(Fi,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            A->r, B->r, k, Fi.one,
            A->UQ, A->ldUQ,
            B->PL, B->ldPL,
            Fi.zero, X, B->r);
        
        RRgen<Field>* RR_X = new RRgen(Fi,A->r,B->r,X,B->r);
        FFLAS::fflas_delete(X);
        
        // LC < LA*LX
        typename Field::Element_ptr L_C = FFLAS::fflas_new(Fi, n, RR_X->r);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n, RR_X->r, A->r, 
            minus ? Fi.mOne : Fi.one, A->PL, A->ldPL,
            RR_X->PL,RR_X->ldPL,
            Fi.zero, L_C, RR_X->r);


        // RC < RX*RB
        typename Field::Element_ptr R_C = FFLAS::fflas_new(Fi, RR_X->r, m);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            RR_X->r, m, B->r, Fi.one,
            RR_X->UQ, RR_X->ldUQ,
            B->UQ, B->ldUQ,
            Fi.zero, R_C, RR_X->r);
        
        
        RRgen<Field>* C = new RRgen(Fi,n,m,RR_X->r,L_C,RR_X->r,R_C,m,true);
        delete(RR_X);
        
        return  C;
    }

/// @brief (algo 5) multiplies two matrices stored as rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param A stored with an RRgen
/// @param B stored with an RRgen
template<class Field>
inline RRgen<Field>* RRxRR (const Field& Fi,RRgen<Field>* A, RRgen<Field>*B){
    return RRxRR (Fi, A, B, false);
}


/// @brief (algo 6) add two matrices stored as rank revealing factorization
/// @tparam Field 
/// @param Fi 
/// @param A stored with an RRgen
/// @param B stored with an RRgen
template<class Field>
inline RRgen<Field>* RRaddRR (const Field& Fi, RRgen<Field>* A, RRgen<Field>*B)
    {
        // X < [LA LB]
        typename Field::Element_ptr X = FFLAS::fflas_new(Fi, A->n, A->r + B->r);
        FFLAS::fassign(Fi,A->n,A->r,A->PL,A->ldPL,X,A->r+ B->r);
        typename Field::Element_ptr X2 = X + A->r;
        FFLAS::fassign(Fi,B->n,B->r,B->PL,B->ldPL,X2,A->r+B->r);

        // Y <  [RA]
        //      [RB]

        typename Field::Element_ptr Y = FFLAS::fflas_new(Fi, A->r + B->r, A->m);
        
        FFLAS::fassign(Fi,A->r,A->m,A->UQ,A->ldUQ,Y,A->m);
        typename Field::Element_ptr Y2 = Y+(A->r)*(A->m);
        FFLAS::fassign(Fi,B->r,B->m,B->UQ,B->ldUQ,Y2,B->m);
        
        
        // LX,RX < facto(X)
        RRgen<Field>* X_fact = new RRgen(Fi, A->n, A->r + B->r, X, A->r + B->r);

        // LY,RY < facto(Y)
        RRgen<Field>* Y_fact = new RRgen(Fi, A->r + B->r, B->m, Y, B->m);


        // D < RRxRR(X,Y)
        RRgen<Field>* D = RRxRR(Fi,X_fact,Y_fact);
        
        FFLAS::fflas_delete(Y);
        FFLAS::fflas_delete(X);
        delete(X_fact);
        delete(Y_fact);

        return D;
    }

/// @brief (algo 7) Adds a quasiseparable matrix in RRR representation and a rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param A        size n*n in RRR representation
/// @param B        size n*n in RR representation
template<class Field>
inline RRRgen<Field>* RRRaddRR (const Field& Fi,  RRRgen<Field>* A,  RRgen<Field>* B)
    {
        if (A->size_N1 <= A->t + B->r){
            // C = RRRexpand(A)

            typename Field::Element_ptr C = FFLAS::fflas_new(Fi, A->size_N1, A->size_N2);
            RRRExpand(Fi, A, C, A->size_N1);
            // B_expanded = RRexpand(B)

            typename Field::Element_ptr B_expanded = FFLAS::fflas_new(Fi,B->n,B->m);
            B->RRExpand(Fi,B_expanded,B->n);
            // C = B+C
            FFLAS::faddin(Fi,A->size_N1,A->size_N1,B_expanded,B->n,C,A->size_N2);
            FFLAS::fflas_delete(B_expanded);
            RRRgen<Field>* D = new RRRgen(Fi,C,A->size_N1,A->t);
            FFLAS::fflas_delete(C);
            
            return D;
        }
        else {
            // B_expanded = [RR_B11 RR_B12]
            //              [RR_B21 RR_B22]

            typename Field::Element_ptr B_expanded = FFLAS::fflas_new(Fi,B->n,B->n);
            B->RRExpand(Fi,B_expanded,B->n);

            typename Field::Element_ptr B11 = B_expanded;
            typename Field::Element_ptr B12 = B11 + A->size_N1;
            typename Field::Element_ptr B21 = B11 + (A->size_N1)*(B->n);
            typename Field::Element_ptr B22 = B21 + A->size_N1;
            
            RRgen<Field>* RR_B11 = new RRgen(Fi,  A->size_N1,  A->size_N1,  B11,  B->n);
            RRgen<Field>* RR_B12 = new RRgen(Fi,  A->size_N1,  A->size_N2,  B12,  B->n);
            RRgen<Field>* RR_B21 = new RRgen(Fi,  A->size_N2,  A->size_N1,  B21,  B->n);
            RRgen<Field>* RR_B22 = new RRgen(Fi,  A->size_N2,  A->size_N2,  B22,  B->n);
            
            FFLAS::fflas_delete(B_expanded);
            std::cout << "RRRaddRR 4 " << std::endl;

            // C11 = RRR+RR(A11,B11)
            RRRgen<Field>* C11 = RRRaddRR(Fi, A->left,RR_B11);
            // C22 = RRR+RR(A22,B22)
            RRRgen<Field>* C22 = RRRaddRR(Fi, A->right,RR_B22);
            // C12 = RR+RR(A12,B12)
            RRgen<Field>* C12 = RRaddRR(Fi, A->LU_right,RR_B12);
            // C21 = RR+RR(A21,B21)
            RRgen<Field>* C21 = RRaddRR(Fi, A->LU_left,RR_B21);
            
            delete RR_B11;
            delete RR_B12;
            delete RR_B21;
            delete RR_B22;
            
            return new RRRgen(C12,C21,A->size_N1,A->size_N2,A->t,C11,C22);
        }
    }

/// @brief (algo 8) Multiplies a quasiseparable matric in RRR representation with a tall and skinny matrix
/// @tparam Field 
/// @param Fi 
/// @param n 
/// @param t 
/// @param A        size n*n in RRR representation
/// @param B        size n*t 
/// @param ldB      leading dimension of B
/// @param C        size n*t
/// @param ldC      leading dimension of C
template<class Field>
inline void RRRxTS (const Field& Fi, size_t n, size_t t,
     const RRRgen<Field>* A, typename Field::ConstElement_ptr B, size_t ldB,
     typename Field::ConstElement_ptr C, size_t ldC)
    {
        if (n<=A->t+t){
            typename Field::Element_ptr Adense = fflas_new (Fi, n, n);
            RRRExpand(Fi, A, Adense, n);
            fgemm(Fi, 
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                n, t, n,
                Fi.one, Adense, n,
                B, ldB,
                0, C, ldC);
                FFLAS::fflas_delete(Adense); 
            return;
        }

        else {
            
            size_t N1 = A->size_N1;
            size_t N2 = A->size_N2;

            // split the matrices as    [C1] = [A11 A12] [B1]
            //                          [C2] = [A12 A22] [B2]
            typename Field::Element_ptr C1 = C;
            typename Field::Element_ptr C2 = C1 + N1*ldC;

            typename Field::Element_ptr B1 = B;
            typename Field::Element_ptr B2 = B1 + N1*ldB;

            // C1 < RRRxTS(A11,B1)
            RRRxTS(Fi, N1, t, A->left, B1, ldB, C1, ldC);

            // C2 < RRRxTS(A22,B2)
            RRRxTS(Fi, N2, t, A->right, B2, ldB, C2, ldC);

            // X < RA12 x B2
            // X of size ru*t
            typename Field::Element_ptr X = fflas_new (Fi, A->LU_right->r, t);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                A->LU_right->r, t, N2,
                Fi.one, A->LU_right->UQ, A->LU_right->ldUQ,
                B2, ldB,
                0, X, t);

            // C1 < C1 + LA12 x X
            fgemm(Fi,
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                N1, t , A->LU_right->r ,
                Fi.one, A->LU_right->PL, A->LU_right->ldPL,
                X, t,
                Fi.one, C1, ldC);

            FFLAS::fflas_delete(X);

            // Y < RA21 x B1
            typename Field::Element_ptr Y = fflas_new (Fi, A->LU_left->r, t);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                A->LU_left->r, t, N1,
                Fi.one,A->LU_left->UQ, A->LU_left->ldUQ,
                B1, ldB,
                0, Y, t);

            // C2 < C2 + LA21 x Y
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                N2, t, A->LU_left->r,
                Fi.one, A->LU_left->PL, A->LU_left->ldPL,
                Y, t,
                Fi.one, C2, ldC);

            FFLAS::fflas_delete(Y);
            // C =   [C1]
            //       [C2]

        }
    }


/// @brief (algo 8) Multiplies a tall and skinny matrix with a quasiseparable matric in RRR representation
/// @tparam Field 
/// @param Fi 
/// @param n 
/// @param t 
/// @param B        size t*n
/// @param ldB      leading dimension of B
/// @param A        size n*n in RRR representation
/// @param C        size t*n
/// @param ldC      leading dimension of C
template<class Field>
inline void TSxRRR (const Field& Fi, size_t n, size_t t,
     typename Field::ConstElement_ptr B, size_t ldB, const RRRgen<Field>* A, 
     typename Field::ConstElement_ptr C, size_t ldC)
    {
        if (n<=A->t+t){
            typename Field::Element_ptr Adense = fflas_new (Fi, n, n);
            RRRExpand(Fi, A, Adense, n);
            fgemm(Fi, 
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, n, n,
                Fi.one,B, ldB, 
                Adense, n,
                0, C, ldC);
            FFLAS::fflas_delete(Adense);
            return;
        }

        else {
            
            size_t N1 = A->size_N1;
            size_t N2 = A->size_N2;

            // split the matrices as                     [A11 A12] 
            //                          [C1 C2] = [B1 B2][A12 A22] 
            typename Field::Element_ptr C1 = C;
            typename Field::Element_ptr C2 = C1 + N1;

            typename Field::Element_ptr B1 = B;
            typename Field::Element_ptr B2 = B1 + N1;

            // C1 < TSxRRR(A11,B1)
            TSxRRR(Fi, N1, t, B1, ldB, A->left, C1, ldC);

            // C2 < RRRxTS(A22,B2)
            TSxRRR(Fi, N1, t, B2, ldB, A->right, C2, ldC);
            

            // X < B2*A21->PL
            // X of size t*r
            typename Field::Element_ptr X = fflas_new (Fi, t, A->LU_right->r);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, A->LU_right->r, N2,
                Fi.one,B2, ldB ,
                A->LU_right->PL, A->LU_right->ldPL,
                0, X,A->LU_right->r);

            // C1 < C1 + X*A21->UQ 
            fgemm(Fi,
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, N1 , A->LU_right->r ,
                Fi.one, X, A->LU_right->r,
                A->LU_right->UQ, A->LU_right->ldUQ,
                Fi.one, C1, ldC);

            FFLAS::fflas_delete(X);

            // Y < B1*A12->PL
            // Y of size t*r
            typename Field::Element_ptr Y = fflas_new (Fi, t, A->LU_left->r);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, A->LU_left->r, N1,
                Fi.one,B1, ldB ,
                A->LU_left->PL, A->LU_left->ldPL,
                0, Y,A->LU_left->r);

            // C2 < C2 + Y*A12->UQ 
            fgemm(Fi,
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, N2 , A->LU_left->r ,
                Fi.one, Y, A->LU_left->r,
                A->LU_left->UQ, A->LU_left->ldUQ,
                Fi.one, C2, ldC);

            FFLAS::fflas_delete(Y);


            // C =   [C1 C2]

        }
    }



/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Computes C = A*B if minus == false , else computes C = -A*B
/// @tparam Field 
/// @param Fi
/// @param A        size n*n in a RRR representation     
/// @param B        size n*m in RR representation
template<class Field>
inline RRgen<Field>* RRxRRR (const Field& Fi,const RRRgen<Field>* A,const RRgen<Field>* B,bool minus)
    {
        // X < TSxRRR(RB,A)
        size_t n = A->size_N1+A->size_N2;
        typename Field::Element_ptr X = fflas_new (Fi, B->r, n);
        TSxRRR(Fi,n,B->r,B->UQ,B->ldUQ,A,X,n);

        // (LX, RX) < RRF(X)
        RRgen<Field>* RR_X = new RRgen(Fi,  B->r,  n,  X,  n);

        // RD < RX
        typename Field::Element_ptr UQ_D = fflas_new (Fi, RR_X->r, n);
        FFLAS::fassign(Fi,RR_X->r,n,RR_X->UQ,RR_X->ldUQ,UQ_D,n);

        // LD < (minus ? (-1) : 1) LBxLX 
        typename Field::Element_ptr PL_D = fflas_new (Fi, B->m, RR_X->r);
        fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                B->m, RR_X->r, B->r,
                minus ? Fi.mOne : Fi.one , B->PL, B->ldPL,
                RR_X->PL, RR_X->ldPL,
                Fi.zero, PL_D, RR_X->r);
        

        RRgen<Field>* D = new RRgen(Fi,B->m,n,RR_X->r,PL_D,RR_X->ldPL,UQ_D,B->m,true);
        delete RR_X;
        return D;

    }

/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Computes C = A*B if minus == false , else computes C = -A*B
/// @tparam Field 
/// @param Fi
/// @param A        size n*n in a RRR representation     
/// @param B        size n*m in RR representation
template<class Field>
inline RRgen<Field>* RRRxRR (const Field& Fi,const RRRgen<Field>* A,const RRgen<Field>* B,bool minus)
    {
        // X < RRRxTS(A,LB)
        size_t n = A->size_N1+A->size_N2;
        typename Field::Element_ptr X = fflas_new (Fi, n, B->r);
        RRRxTS(Fi,n,B->r,A,B->PL,B->ldPL,X,B->r);

        // (LX, RX) < RRF(X)
        RRgen<Field>* RR_X = new RRgen(Fi,  n,  B->r,  X,  B->r);

        // LD < LX
        typename Field::Element_ptr PL_D = fflas_new (Fi, n, RR_X->r);
        FFLAS::fassign(Fi,n,RR_X->r,RR_X->PL,RR_X->ldPL,PL_D,RR_X->r);

        // RD < (minus ? (-1) : 1) RX x RB
        typename Field::Element_ptr UQ_D = fflas_new (Fi, RR_X->r, B->m);
        fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                RR_X->r, B->m, B->r,
                minus ? Fi.mOne : Fi.one , RR_X->UQ, RR_X->ldUQ,
                B->UQ, B->ldUQ,
                Fi.zero, UQ_D, RR_X->r);
        

        RRgen<Field>* D = new RRgen(Fi,n,B->m,RR_X->r,PL_D,RR_X->ldPL,UQ_D,B->m,true);
        delete RR_X;
        return D;

    }


/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Computes C = A*B
/// @tparam Field 
/// @param Fi
/// @param A        size n*n in a RRR representation     
/// @param B        size n*m in RR representation
template<class Field>
inline RRgen<Field>* RRRxRR (const Field& Fi,const RRRgen<Field>* A,const RRgen<Field>* B){
    return RRRxRR(Fi,A,B,false);
}
/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Computes C = A*B
/// @tparam Field 
/// @param Fi
/// @param A        size n*n in a RRR representation     
/// @param B        size n*m in RR representation
template<class Field>
inline RRgen<Field>* RRxRRR (const Field& Fi,const RRRgen<Field>* A,const RRgen<Field>* B){
    return RRxRRR(Fi,A,B,false);
}

/// @brief multiplies two QS matrices (algo 10)
/// @tparam Field 
/// @param Fi 
/// @param A        size n*n
/// @param B        size n*n
template<class Field>
inline RRRgen<Field>* RRRxRRR (const Field& Fi, const RRRgen<Field>* A, const RRRgen<Field>* B){   
    size_t n = A->size_N1+A->size_N2;
    if (n<= A->t+B->t){
        //return RRRExpand(A) x RRRExpand(B)
        typename Field::Element_ptr A_expanded = fflas_new (Fi, n, n);
        typename Field::Element_ptr B_expanded = fflas_new (Fi, n, n);
        typename Field::Element_ptr C = fflas_new (Fi, n, n);
        fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n, n, n,
            Fi.one, A_expanded, n,
            B_expanded, n,
            0, C, n);
        FFLAS::fflas_delete(A_expanded);
        FFLAS::fflas_delete(B_expanded);
        RRRgen<Field>* D = new RRRgen(Fi,C,n,A->t+B->t);
        FFLAS::fflas_delete(C);
        return D;
    }

    else {
        // C11 < RRRxRRR(A11,B11)
        RRRgen<Field>* C11 = RRRxRRR(Fi,A->left,B->left);

        // C22 < RRxRRR(A22,B22)
        RRRgen<Field>* C22 = RRRxRRR(Fi,A->right,B->right);

        // X < RRxRR(A12,B21)
        RRgen<Field>* X = RRxRR(Fi,A->LU_right,B->LU_left);

        // Y < RRxRR(A21,B12)
        RRgen<Field>* Y = RRxRR(Fi,A->LU_left,B->LU_right);

        // C11 < RRR+RR(C11,X)
        C11 = RRRaddRR(Fi,C11,X);

        // C22 < RRR+RR(C22,Y)
        C22 = RRRaddRR(Fi,C22,Y);
        
        delete X;
        delete Y;

        // LX < RRRxTS(A11,LB12) ; RX < RB12
        typename Field::Element_ptr LX = fflas_new (Fi, A->size_N1, B->LU_right->r);
        RRRxTS(Fi,A->size_N1,B->LU_right->r,A->left,B->LU_right->PL,B->LU_right->ldPL,LX,B->LU_right->r);
        
        typename Field::Element_ptr RX = fflas_new (Fi, B->LU_right->r, B->size_N2);
        FFLAS::fassign(Fi,B->LU_right->r,B->size_N2,
            B->LU_right->UQ,B->LU_right->ldUQ,
            RX,B->size_N2);

        // LY < LA12 ; RY < TSxRRR(RA12,B22) 
        typename Field::Element_ptr LY = fflas_new (Fi, A->size_N1, A->LU_right->r);
        FFLAS::fassign(Fi,A->size_N1,A->LU_right->r,
            A->LU_right->PL,A->LU_right->ldPL,
            LY,A->LU_right->r);
        
        typename Field::Element_ptr RY = fflas_new (Fi, A->LU_right->r, B->size_N2);
        TSxRRR(Fi,B->size_N2,A->LU_right->r,A->LU_right->UQ,A->LU_right->ldUQ,B->right,RY,B->size_N2);
        
        X = new RRgen(Fi, A->size_N1, B->size_N2,B->LU_right->r,LX,B->LU_right->r,RX,B->size_N2,true);
        Y = new RRgen(Fi, A->size_N1, B->size_N2,A->LU_right->r,LY,A->LU_right->r,RY,B->size_N2,true);

        // C12 < RR+RR(X,Y)
        RRgen<Field> C12 = RRaddRR(Fi,X,Y);
        delete X;
        delete Y;

        // LX < RRRxTS(A22,LB21) ; RX <  RB21
        LX = fflas_new (Fi, A->size_N2, B->LU_left->r);
        RRRxTS(Fi,A->size_N2,B->LU_left->r,A->right,B->LU_left->PL,B->LU_left->ldPL,LX,B->LU_left->r);
        
        RX = fflas_new (Fi, B->LU_left->r, B->size_N1);
        FFLAS::fassign(Fi,B->LU_left->r,B->size_N1,
            B->LU_left->UQ,B->LU_left->ldUQ,
            RX,B->size_N1);

        // LY < LA21 ; RY < TSxRRR(RA21,B11)
        LY = fflas_new (Fi, A->size_N2, A->LU_left->r);
        FFLAS::fassign(Fi,A->size_N2,A->LU_left->r,
            A->LU_left->PL,A->LU_left->ldPL,
            LY,A->LU_left->r);
        
        RY = fflas_new (Fi, A->LU_left->r, B->size_N1);
        TSxRRR(Fi,B->size_N1,A->LU_left->r,A->LU_left->UQ,A->LU_left->ldUQ,B->left,RY,B->size_N1);
        
        X = new RRgen(Fi, A->size_N2, B->size_N1,B->LU_left->r,LX,B->LU_left->r,RX,B->size_N1,true);
        Y = new RRgen(Fi, A->size_N2, B->size_N1,A->LU_left->r,LY,A->LU_left->r,RY,B->size_N1,true);
        
        // C21 < RR+RR(X,Y)
        RRgen<Field> C21 = RRaddRR(Fi,X,Y);
        delete X;
        delete Y;

        // RETURN C =   [C11    C12]
        //              [C12    C21]
        return new RRRgen(C12,C21,A->size_N1,B->size_N2,A->t+B->t,C11,C22);
    }

}


/// @brief Compute the inverse in RRR representation
/// @tparam Field 
/// @param Fi 
/// @param A    in RRR representation
template<class Field>
inline RRRgen<Field>* RRRinvert (const Field& Fi,
            const RRRgen<Field>* A)
    {
        if (!A.left){
            // Y < RRRExpand(A)
            size_t N1 = A->size_N1;
            typename Field::Element_ptr Y = fflas_new (Fi, N1, N1);
            RRRExpand(Fi,A->LU_right,Y,N1);

            // return Invert(Y)
            int nullity;
            Invert (Fi, N1,Y, N1, nullity);
            RRRgen<Field>* Y_RRR = new RRRgen(Fi,Y,N1,A->t);
            FFLAS::fflas_delete(Y);
            return  Y_RRR;
        }

        // split the matrix as A =  [A11 A12] and   X = [X11 X12]
        //                          [A21 A22]           [X21 X22]

        // Y11 < RRRinvertrec(A11)
        RRRgen<Field>* Y11 = RRRinvert(Fi,A->left);

        // Y12 < RRRxRR(Y11,A12)
        RRgen<Field>* Y12 = RRRxRR(Fi,Y11,A->LU_right);

        // Y21 < RRxRRR(A21,Y11)
        RRgen<Field>* Y21 = RRxRRR(Fi,Y11,A->LU_left);

        // Z < RRxRR(A21,Y12)
        RRgen<Field>* Z = RRxRR(A->LU_left,Y12);

        // D < RRRaddRR(A22,Z)
        RRRgen<Field>* D = RRRaddRR(Fi,A->right,Z);
        delete Z;

        // X22 < RRRinvert(D)
        RRRgen<Field> X22 = RRRinvert(Fi,D);
        delete D;

        // X21 < -RRRxRR(X22,Y21)
        RRgen<Field>* X21 = RRRxRR(Fi,X22,Y21,true);

        // W < -RRxRR(Y12,X21)
        RRgen<Field>* W = RRxRR(Fi,Y12,X21,true); 

        // X12 < -RRxRRR(Y12,X22)
        RRgen<Field>* X12 = RRxRRR(Fi,X22,Y12,true);

        // X11 < RRRaddRR(Y11,W)
        RRRgen<Field>* X11 = RRRaddRR(Fi, Y11, W);
        delete W;

        // return X
        RRRgen<Field>* X = new RRRgen(X12,X21,A->size_N1,A->size_N2,A->t,X11,X22);  
    }
}

#endif //_FFPACK_ffpack_rrrgen_inl