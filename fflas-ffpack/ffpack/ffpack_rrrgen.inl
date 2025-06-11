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
    

    RRgen(const Field& Fi, size_t n, size_t m, size_t r,
            typename Field::Element_ptr PL, size_t ldPL, 
            typename Field::Element_ptr UQ,size_t ldUQ)
        : n(n), m(m), r(r), PL(PL), ldPL(ldPL), UQ(UQ), ldUQ(ldUQ) {}
    
    RRgen(const Field& Fi, size_t n_A, size_t m_A, typename Field::Element_ptr A, size_t ldA){
        n = n_A;
        m = m_A;
        size_t* P = FFLAS::fflas_new<size_t>(n);
        size_t* Q = FFLAS::fflas_new<size_t>(m);

        r = PLUQ(Fi, FFLAS::FflasNonUnit, n, m, A, ldA, P, Q);
        

        PL = FFLAS::fflas_new(Fi, n, r);
        UQ = FFLAS::fflas_new(Fi, r, m);
        ldPL = r;
        ldUQ = m;

        // extraction of U
        getTriangular<Field>(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit, n, m, r, A, ldA, UQ, m, true);

        // extraction of L
        getTriangular<Field>(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, n, m, r, A, ldA, PL, r, true);

        // apply the permutations (P on L and Q on U)
        applyP<Field>(Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, r, 0, n, PL, r, P);
        

        applyP<Field>(Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, r, 0, m, UQ, m, Q);
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
            LU_right = new RRgen<Field>(Fi,N,N,0,A,N,nullptr,0);
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
            const RRRgen<Field>& A,
            typename Field::Element_ptr B, const size_t ldb)
    {
        if (A.left == nullptr){
            size_t N = A.size_N1;
            FFLAS::fassign(Fi, N, N, A.LU_right->PL, ldb, B, ldb);
            return;
        }

        size_t N1 = A.size_N1;
        size_t N2 = A.size_N2;

        // cut the matrix A in four parts
        // in case of odd N, A12 and A21 are rectangles and not squares
        typename Field::Element_ptr B11 = B;
        typename Field::Element_ptr B12 = B + N1;
        typename Field::Element_ptr B21 = B + ldb*N1;
        typename Field::Element_ptr B22 = B21 + N1;

        
        // B11 < RRRExpand(A11)
        RRRExpand<Field>(Fi, *A.left, B11, ldb);

        // B22 < RRRExpand(A22)
        RRRExpand<Field>(Fi, *A.right, B22, ldb);

        // B12 < L_u * U_u
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N1, N2, A.LU_right->r, Fi.one,
            A.LU_right->PL, A.LU_right->r,
            A.LU_right->UQ, N2,
            Fi.zero, B12, ldb);

        
        // B21 < L_l * U_l
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N2, N1, A.LU_left->r, Fi.one,
            A.LU_left->PL, A.LU_left->r,
            A.LU_left->UQ, N1,
            Fi.zero, B21, ldb);

        
    }



/// @brief (algo 5) multiplies two matrices stored as rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param A stored with an RRgen
/// @param B stored with an RRgen
template<class Field>
inline RRgen<Field>* RRxRR (const Field& Fi,RRgen<Field>* A, RRgen<Field>*B)
    {   
        size_t m = B->size_N2;
        size_t n = A->size_N1;
        size_t k = A->size_N2;
        
        // X < UA * LB
        typename Field::Element_ptr X = FFLAS::fflas_new(Fi, A->r, B->r);
        fgemm(Fi,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            A->r, B->r, k, Fi.one,
            A->UQ, A->ldUQ,
            B->PL, B->ldPL,
            Fi.zero, X, B->r);
        
        // LX,RX < facto(X)

        size_t* P_X = FFLAS::fflas_new<size_t>(A->r);
        size_t* Q_X = FFLAS::fflas_new<size_t>(B->r);

        size_t r_X = PLUQ(Fi, FFLAS::FflasNonUnit,
                            A->r, B->r,
                            X, B->r,
                            P_X, Q_X);

        typename Field::Element_ptr R_X = FFLAS::fflas_new(Fi, r_X, B->r);
        typename Field::Element_ptr L_X = FFLAS::fflas_new(Fi, A->r, r_X);

        // extraction of R_X
        getTriangular<Field>(Fi, FFLAS::FflasUpper,
                        FFLAS::FflasNonUnit,
                        A->r, B->r, r_X,
                        X, B->r,
                        R_X, B->r,
                        true);

        // extraction of L_X
        getTriangular<Field>(Fi, FFLAS::FflasLower,
                        FFLAS::FflasUnit,
                        A->r, B->r, r_X,
                        X, B->r,
                        L_X, r_X,
                        true);

        // apply the permutations (P on L_X and Q on R_X)
        applyP<Field>(Fi,
                FFLAS::FflasLeft, FFLAS::FflasTrans,
                r_X, 0, A->r,
                L_X, r_X, P_X);
        
                applyP<Field>(Fi,
                FFLAS::FflasRight, FFLAS::FflasNoTrans,
                r_X, 0, B->r,
                R_X,B->r, Q_X);
        
        // LC < LA*LX
        typename Field::Element_ptr L_C = FFLAS::fflas_new(Fi, m, r_X);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            m, r_X, A->r, Fi.one,
            A->PL, A->ldPL,
            L_X,r_X,
            Fi.zero, L_C, r_X);

        // RC < RX*RB
        typename Field::Element_ptr R_C = FFLAS::fflas_new(Fi, r_X, n);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            r_X, n, B->r, Fi.one,
            R_X, B->r,
            B->UQ, B->ldUQ,
            Fi.zero, L_C, r_X);
        

        FFLAS::fflas_delete(X);
        FFLAS::fflas_delete(Q_X);
        FFLAS::fflas_delete(P_X);

        return new RRgen(Fi,m,n,r_X,L_C,r_X,R_C,n); 
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
        typename Field::Element_ptr X = FFLAS::fflas_new(Fi, A->size_N1, A->r + B->r);
        fzero(Fi, A->size_N1, A->r + B->r, X, A->r + B->r);
        fassign(Fi,A->size_N1,A->r,A->PL,A->ldPL,X,A->r);
        fassign(Fi,B->size_N1,B->r,B->PL,B->ldPL,X+A->r,B->r);
        
        // Y <  [RA]
        //      [RB]
        typename Field::Element_ptr Y = FFLAS::fflas_new(Fi, A->r + B->r, B->size_N2);
        fzero(Fi, A->r + B->r, B->size_N2, Y, B->size_N2);
        fassign(Fi,A->r,A->size_N2,A->UQ,A->ldUQ,X,A->size_N2);
        fassign(Fi,B->r,B->size_N2,B->UQ,Y+(A->r)*(A->size_N2),B->size_N2);
        
        
        // LX,RX < facto(X)
        RRgen<Field> X_fact = RRgen(Fi, A->size_N1, A->r + B->r, X, A->r + B->r);

        // LY,RY < facto(Y)
        RRgen<Field> Y_fact = RRgen(Fi, A->r + B->r, B->size_N2, Y, B->size_N2);


        // D < RRxRR(X,Y)
        RRgen<Field>* D = RRxRR(Fi,X_fact,Y_fact);
        
        FFLAS::fflas_delete(Y);
        FFLAS::fflas_delete(X);
        delete X_fact;
        delete Y_fact;

        return D;
    }

/// @brief (algo 7) Adds a quasiseparable matrix in RRR representation and a rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param A        size n*n in RRR representation
/// @param B        size n*n in RR representation
template<class Field>
inline RRRgen<Field>* RRRaddRR (const Field& Fi, const RRRgen<Field>* A, const RRgen<Field>* B)
    {
        if (A->size_N1 <= A->t + B->r){
            // C = RRRexpand(A)
            typename Field::Element_ptr C = FFLAS::fflas_new(Fi, A->size_N1, A->size_N1);
            RRRExpand(Fi, A, C, A->size_N1);
            // B_expanded = RRexpand(B)
            typename Field::Element_ptr B_expanded = FFLAS::fflas_new(Fi,B->n,B->n);
            B.RRExpand(Fi,B_expanded,B->n);

            // C = B+C
            faddin(Fi,A->size_N1,A->size_N1,B_expanded,B->n,C,A->n);
            FFLAS::fflas_delete(B_expanded);
            return new RRRgen(Fi,C,A->n,A->t);
        }
        else {
            // B_expanded = [RR_B11 RR_B12]
            //              [RR_B21 RR_B22]
            typename Field::Element_ptr B_expanded = FFLAS::fflas_new(Fi,B->n,B->n);
            B.RRExpand(Fi,B_expanded,B->n);
            
            typename Field::Element_ptr B11 = B;
            typename Field::Element_ptr B12 = B11 + A->size_N1;
            typename Field::Element_ptr B21 = B11 + (A->size_N1)*(B->n);
            typename Field::Element_ptr B22 = B21 + A->size_N1;
            
            RRgen<Field>* RR_B11 = new RRgen(Fi,  A->size_N1,  A->size_N1,  B11,  B->n);
            RRgen<Field>* RR_B12 = new RRgen(Fi,  A->size_N1,  A->size_N2,  B12,  B->n);
            RRgen<Field>* RR_B21 = new RRgen(Fi,  A->size_N2,  A->size_N1,  B21,  B->n);
            RRgen<Field>* RR_B22 = new RRgen(Fi,  A->size_N2,  A->size_N2,  B22,  B->n);
            
            FFLAS::fflas_delete(B_expanded);

            // C11 = RRR+RR(A11,B11)
            RRRgen<Field>* C11 = RRRaddRR(Fi, A->left,B11);
            // C22 = RRR+RR(A22,B22)
            RRRgen<Field>* C22 = RRRaddRR(Fi, A->right,B22);
            // C12 = RR+RR(A12,B12)
            RRgen<Field>* C12 = RRaddRR(Fi, A->LU_right,B12);
            // C21 = RR+RR(A21,B21)
            RRgen<Field>* C21 = RRaddRR(Fi, A->LU_left,B21);
            
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
                delete Adense;
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



/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Computes C = A*B
/// @tparam Field 
/// @param Fi
/// @param A        size n*n in a RRR representation     
/// @param B        size n*m in RR representation
template<class Field>
inline RRgen<Field>* RRRxRR (const Field& Fi,const RRRgen<Field>* A,const RRgen<Field>* B)
    {
        // X < RRRxTS(A,LB)
        size_t n = A->size_N1+A->size_N2;
        typename Field::Element_ptr X = fflas_new (Fi, n, B->r);
        RRRxTS(Fi,n,B->r,A,B->PL,B->ldPL,X,B->r);

        // (LX, RX) < RRF(X)
        RRgen<Field>* RR_X = new RRgen(Fi,  n,  B->r,  X,  B->r);

        // LD < LX
        typename Field::Element_ptr PL_D = fflas_new (Fi, n, RR_X->r);
        fassign(Fi,n,RR_X->r,RR_X->PL,RR_X->ldPL,PL_D,RR_X->r);

        // RD < RX x RB
        typename Field::Element_ptr UQ_D = fflas_new (Fi, RR_X->r, B->m);
        fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                RR_X->r, B->m, B->r,
                Fi.one, RR_X->UQ, RR_X->ldUQ,
                B->UQ, B->ldUQ,
                0, UQ_D, RR_X->r);
        

        RRgen<Field>* D = new RRgen(Fi,n,B->m,RR_X->r,PL_D,RR_X->ldPL,UQ_D,B->m);
        delete RR_X;
        return D;

    }

/// @brief multiplies two QS matrices (algo 10)
/// @tparam Field 
/// @param Fi 
/// @param A        size n*n
/// @param B        size n*n
template<class Field>
inline RRRgen<Field>* RRRxRRR (const Field& Fi,
            const RRRgen<Field>* A,
            const RRRgen<Field>* B)
    {   
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
            RRRgen<Field>* D = new RRRgen(Fi,C,n,A->t);
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

            // LX < RRRxTS(A11,LB12) ; RX < RB12
            typename Field::Element_ptr LX = fflas_new (Fi, A->size_N1, B->LU_right->r);
            RRRxTS(Fi,A->size_N1,B->LU_right->r,A->left,B->LU_right->UQ,B->LU_right->ldUQ,LX,r);
            
            typename Field::Element_ptr RX = fflas_new (Fi, B->LU_right->r, A->size_N2);
            fassign(Fi,size_N1,size_N2,B->LU_right->UQ,B->LU_right->ldUQ,RX,A->size_N2);

            // LY < LA12 ; RY < TSxRRR(RA12,B22) ATTENTION ON APPLIQUE à GAUCHE

            // C12 < RR+RR(X,Y)

            // LX < RRRxTS(A11,LB21) ; RX <  RB21

            // LY < LA21 ; RY < TSxRRR(RA21,B22)

            // C21 < RR+RR(X,Y)

            // RETURN C =   [C11    C12]
            //              [C12    C21]
        }

        }


/// @brief Compute the inverse in RRR representation
/// @tparam Field 
/// @param Fi 
/// @param A    in RRR representation
template<class Field>
inline void RRRinvert (const Field& Fi,
            const RRRgen<Field>& A)
    {
        RRRinvertrec(Fi, A);
    }

template<class Field>
inline void RRRinvertrec (const Field& Fi,
            const RRRgen<Field>& A)
    {
        if (!A.left){
            // Y < RRRExpand(A)
            // return invert(Y)
        }

        // split the matrix as A =  [A11 A12] and   X = [X11 X12]
        //                          [A21 A22]           [X21 X22]

        // Y11 < RRRinvertrec(A11)

        // Y12 < RRRxRR(Y11,A12)

        // Y21 < RRRxRR(A21,Y11)

        // Z < RRxRR(A21,Y12)

        // D < RRaddRR(A22,Z)

        // X22 < RRRinvert(D)

        // X21 < -RRRxRR(X22,Y21)

        // W < -RRxRR(Y12,X21)

        // X12 < -RRRxRR(Y12,X22)

        // X11 < RRRaddRR(Y11,W)

        // return X  
    }
}

#endif //_FFPACK_ffpack_rrrgen_inl