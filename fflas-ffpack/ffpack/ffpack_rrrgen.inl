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
    RRRgen* left;                    // recursively pointing on the same representation of the up left submatrix = A11
    RRRgen* right;                   // recursively pointing on the same representation of the down right submatrix = A22


    RRRgen(RRgen<Field>* LU_right, RRgen<Field>* LU_left, size_t size_N1, size_t size_N2, RRRgen* left, RRRgen* right)
        :  LU_right(LU_right), LU_left(LU_left), size_N1(size_N1), size_N2(size_N2), left(left), right(right){}
    
    
    RRRgen(   const Field& F,typename Field::Element_ptr leaf, size_t n){
        LU_right = new RRgen<Field>(F,n,n,0,leaf,n,nullptr,0);
        LU_left = nullptr;
        size_N1 = n;
        size_N2 = 0;
        left = nullptr;
        right = nullptr; 
    }

    RRRgen( const Field& Fi,const size_t N, const size_t s,typename Field::Element_ptr A, const size_t lda)
    {
        if (N/2 <= s) { 
            // leaf
            LU_right = new RRgen<Field>(Fi,N,N,0,A,N,nullptr,0);
            LU_left = nullptr;
            size_N1 = N;
            size_N2 = 0;
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
        left = new RRRgen( Fi, size_N1, s, A11, lda);

        right = new RRRgen( Fi, size_N2, s, A22, lda);
        return;
    }

    ~RRRgen() {
        FFLAS::fflas_delete(LU_right);
        
        if (LU_left) {
            FFLAS::fflas_delete(LU_left);
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
/// @param r_A 
/// @param m
/// @param r_B 
/// @param n 
/// @param k 
/// @param LA       size m*r_A
/// @param UA       size r_A*k 
/// @param LB       size k*r_B
/// @param UB       size r_B*n 
/// @param lda      leading dimension of A
/// @param ldb      leading dimension of B
/// @param LC       size m*r_C don't need to be initialized   
/// @param RC       size r_C*n don't need to be initialized
/// @param r_C      

template<class Field>
inline RRgen<Field>* RRxRR (const Field& Fi,
            size_t r_A, size_t m, size_t r_B, size_t n, size_t k,
            typename Field::ConstElement_ptr LA, size_t ldLA, 
            typename Field::ConstElement_ptr UA, size_t ldUA,
            typename Field::ConstElement_ptr LB, size_t ldLB,
            typename Field::ConstElement_ptr UB, size_t ldUB)
    {
        // X < UA * LB
        typename Field::Element_ptr X = FFLAS::fflas_new(Fi, r_A, r_B);
        fgemm(Fi,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            r_A, r_B, k, Fi.one,
            UA, ldUA,
            LB, ldLB,
            Fi.zero, X, r_B);
        
        // LX,RX < facto(X)

        size_t* P_X = FFLAS::fflas_new<size_t>(r_A);
        size_t* Q_X = FFLAS::fflas_new<size_t>(r_B);

        size_t r_X = PLUQ(Fi, FFLAS::FflasNonUnit,
                            r_A, r_B,
                            X, r_B,
                            P_X, Q_X);

        typename Field::Element_ptr R_X = FFLAS::fflas_new(Fi, r_X, r_B);
        typename Field::Element_ptr L_X = FFLAS::fflas_new(Fi, r_A, r_X);

        // extraction of R_X
        getTriangular<Field>(Fi, FFLAS::FflasUpper,
                        FFLAS::FflasNonUnit,
                        r_A, r_B, r_X,
                        X, r_B,
                        R_X, r_B,
                        true);

        // extraction of L_X
        getTriangular<Field>(Fi, FFLAS::FflasLower,
                        FFLAS::FflasUnit,
                        r_A, r_B, r_X,
                        X, r_B,
                        L_X, r_X,
                        true);

        // apply the permutations (P on L_X and Q on R_X)
        applyP<Field>(Fi,
                FFLAS::FflasLeft, FFLAS::FflasTrans,
                r_X, 0, r_A,
                L_X, r_X, P_X);
        
                applyP<Field>(Fi,
                FFLAS::FflasRight, FFLAS::FflasNoTrans,
                r_X, 0, r_B,
                R_X,r_B, Q_X);
        
        // LC < LA*LX
        typename Field::Element_ptr L_C = FFLAS::fflas_new(Fi, m, r_X);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            m, r_X, r_A, Fi.one,
            LA, ldLA,
            L_X,r_X,
            Fi.zero, L_C, r_X);

        // RC < RX*RB
        typename Field::Element_ptr R_C = FFLAS::fflas_new(Fi, r_X, n);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            r_X, n, r_B, Fi.one,
            R_X, r_B,
            UB, ldUB,
            Fi.zero, L_C, r_X);

        return new RRgen(Fi,m,n,r_X,L_C,r_X,R_C,n); 
    }

/// @brief (algo 6) add two matrices stored as rank revealing factorization
/// @tparam Field 
/// @param Fi 
/// @param r_A 
/// @param r_B 
/// @param m 
/// @param n 
/// @param LA       size (m*r_A) transpose
/// @param UA       size r_A*n
/// @param LB       size (m*r_B) transpose
/// @param UB       size r_B*n
/// @param ldLA     leading dimension of LA
/// @param ldUA     leading dimension of UA
/// @param ldLB     leading dimension of LB
/// @param ldUB     leading dimension of UB
template<class Field>
inline void RRaddRR (const Field& Fi,
            size_t r_A, size_t r_B, size_t m, size_t n,
            typename Field::ConstElement_ptr LA, size_t ldLA,
            typename Field::ConstElement_ptr UA, size_t ldUA,
            typename Field::ConstElement_ptr LB, size_t ldLB,
            typename Field::ConstElement_ptr UB, size_t ldUB)
    {
        // X < [LA LB]
        // fassign
        // Y <  [RA]
        //      [RB]

        // LX,RX < facto(X)

        // LY,RY < facto(Y)

        // D < RRxRR(X,Y)
    }

/// @brief (algo 7) Adds a quasiseparable matrix in RRR representation and a rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param s        order of quasiseparability of A
/// @param r_B      rank of b
/// @param n        dimension of A and B
/// @param A        size n*n in RRR representation
/// @param LB       size n*r_B
/// @param UB       size r_B*n 
template<class Field>
inline void RRRaddRR (const Field& Fi,
            size_t s, size_t r_B, size_t n,
            const RRRgen<Field>& A,
            typename Field::ConstElement_ptr LB, size_t ldLB,
            typename Field::ConstElement_ptr UB, size_t ldUB)
    {
        //TODO
    }

/// @brief (algo 8) Multiplies a quasiseparable matric in RRR representation with a tall and skinny matrix
/// @tparam Field 
/// @param Fi 
/// @param s        order of quasisep of A
/// @param t 
/// @param n 
/// @param A        size n*n in RRR representation
/// @param B        size n*t 
/// @param ldB      leading dimension of B
/// @param C        size n*t
/// @param ldC      leading dimension of C
template<class Field>
inline void RRRxTSrec (const Field& Fi,
            size_t s, size_t n, size_t t,
            const RRRgen<Field>& A,
            typename Field::ConstElement_ptr B, size_t ldB,
            typename Field::ConstElement_ptr C, size_t ldC)
    {
        if (n<=s+t){
            typename Field::Element_ptr Adense = fflas_new (Fi, n, n);
            RRRExpand(Fi, A, Adense, n);
            fgemm(Fi, 
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                n, t, n,
                1,
                Adense, n,
                B, ldB,
                0,
                C, ldC);

            return;
        }

        size_t N1 = A.size_N1;
        size_t N2 = A.size_N2;

        // split the matrices as    [C1] = [A11 A12] [B1]
        //                          [C2] = [A12 A22] [B2]
        typename Field::Element_ptr C1 = C;
        typename Field::Element_ptr C2 = C1 + N1*ldC;

        typename Field::Element_ptr B1 = B;
        typename Field::Element_ptr B2 = B1 + N1*ldB;

        // C1 < RRRxTS(A11,B1)
        RRRxTSrec(Fi, s, N1, t, A.left, B1, ldB, C1, ldC);

        // C2 < RRRxTS(A22,B2)
        RRRxTSrec(Fi, s, N2, t, A.right, B2, ldB, C2, ldC);

        // X < RA12 x B2
        // X of size ru*t
        typename Field::Element_ptr X = fflas_new (Fi, A.LU_right->r, t);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            A.LU_right->r, N2, t,
            1,
            A.LU_right->UQ, N2,
            B2, ldB,
            0,
            X, t);

        // C1 < C1 + LA12 x X
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N1, A.LU_right->r, t,
            1,
            A.LU_right, A.LU_right->r,
            X, t,
            1,
            C1, ldC);

        // Y < RA21 x B1
        typename Field::Element_ptr Y = fflas_new (Fi, A.LU_left->r, t);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            A.LU_left->r, N1, t,
            1,
            A.U_l, N1,
            B1, ldB,
            0,
            Y, t);

        // C2 < C2 + LA21 x Y
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N2, A.LU_left->r, t,
            1,
            A.L_l, A.LU_left->r,
            Y, t,
            1,
            C2, ldC);

        // C =   [C1]
        //       [C2]
    }

/// @brief (algo 8) Multiplies a quasiseparable matric in RRR representation with a tall and skinny matrix
/// @tparam Field 
/// @param Fi 
/// @param s        order of quasisep of A
/// @param t 
/// @param n 
/// @param A        size n*n in RRR representation
/// @param B        size n*t 
/// @param ldB      leading dimension of B
/// @param C        size n*t
/// @param ldC      leading dimension of C
template<class Field>
inline void RRRxTS (const Field& Fi,
            size_t s, size_t n, size_t t,
            const RRRgen<Field>& A,
            typename Field::ConstElement_ptr B, size_t ldB,
            typename Field::ConstElement_ptr C, size_t ldC)
    {
        // RRRxTSrec(Fi, s, n, t, A.getroot(), B, ldB, C, ldC);
    }

/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Compute C = A*B
/// @tparam Field 
/// @param Fi
/// @param s        order of QS of A
/// @param n   
/// @param r_B      rank of B
/// @param m  
/// @param A        size n*n in a RRR representation     
/// @param LB       size n*r_B
/// @param UB       size r_B*m
/// @param ldLB     leading dimension of LB
/// @param ldUB     leading dimension of UB
template<class Field>
inline void RRRxRR (const Field& Fi,
            size_t s, size_t n, size_t r_B, size_t m,
            const RRRgen<Field>& A,
            typename Field::ConstElement_ptr LB, size_t ldLB,
            typename Field::ConstElement_ptr UB, size_t ldUB,
            typename Field::Element_ptr LC, size_t ldLC,
            typename Field::Element_ptr UC, size_t ldUC)
    {
        // X < RRRxTS(A,LB)

        // (LX, UX) < RRF(X)

        // LD < LX

        // RD < RX x RB
    }

/// @brief multiplies two QS matrices (algo 10)
/// @tparam Field 
/// @param Fi 
/// @param s        order of QS of A
/// @param t        order of QS of B
/// @param n 
/// @param A        size n*n
/// @param B        size n*n
template<class Field>
inline void RRRxRRR (const Field& Fi,
            size_t s,size_t t, size_t n,
            const RRRgen<Field>& A,
            const RRRgen<Field>& B)
    {
        if (n<= s+t){
            //return RRRExpand(A) x RRRExpand(B)
        }

        // C11 < RRRxRRR(A11,B11)
        
        // C22 < RRxRRR(A22,B22)

        // X < RRxRR(A12,B21)

        // Y < RRxRR(A21,B12)

        // C11 < RRR+RR(C11,X)

        // C22 < RRR+RR(C22,Y)

        // LX < RRRxTS(A11,LB12) ; RX < RB12

        // LY < LA12 ; RY < TSxRRR(RA12,B22) ATTENTION ON APPLIQUE à GAUCHE

        // C12 < RR+RR(X,Y)

        // LX < RRRxTS(A11,LB21) ; RX <  RB21

        // LY < LA21 ; RY < TSxRRR(RA21,B22)

        // C21 < RR+RR(X,Y)

        // RETURN C =   [C11    C12]
        //              [C12    C21]
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