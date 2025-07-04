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
    bool memory_owner;                  // whether the structure owns the memory or whether it is a view on memory area. Used in destructors

    RRgen(const Field& Fi, size_t n, size_t m, size_t r,
            typename Field::Element_ptr PL, size_t ldPL, 
            typename Field::Element_ptr UQ,size_t ldUQ, bool memory_owner = true)
        : n(n), m(m), r(r), PL(PL), ldPL(ldPL), UQ(UQ), ldUQ(ldUQ), memory_owner(memory_owner) {}
    
    
    RRgen(const Field& Fi, size_t n_A, size_t m_A, size_t r_A,
        typename Field::Element_ptr PL_A, size_t ldPL_A, 
        typename Field::Element_ptr UQ_A,size_t ldUQ_A, bool mem_owner, bool copy):n(n_A),m(m_A),r(r_A){
        if (copy){
            PL = FFLAS::fflas_new(Fi, n_A, r_A);
            FFLAS::fassign(Fi,n_A,r_A,PL_A,ldPL_A,PL,r_A);
            ldPL = r_A;
            memory_owner = true;
            if (UQ_A){
                UQ = FFLAS::fflas_new(Fi, r_A, m_A);
                FFLAS::fassign(Fi,r_A,m_A,UQ_A,ldUQ_A,UQ,m_A);
                ldUQ = m_A;
            }
            else {
                UQ = UQ_A;
                ldUQ = 0;
            }
        }
        else {
            ldPL = ldPL_A;
            ldUQ = ldUQ_A;
            PL = PL_A;
            UQ = UQ_A;
            memory_owner = mem_owner;
        }
    }
    
    
    RRgen(const Field& Fi, size_t n_A, size_t m_A, typename Field::ConstElement_ptr A, size_t ldA):n(n_A),m(m_A){
        size_t* P = FFLAS::fflas_new<size_t>(n_A);
        size_t* Q = FFLAS::fflas_new<size_t>(m_A);

        
        typename Field::Element_ptr A_copy = FFLAS::fflas_new(Fi, n_A, m_A);
        FFLAS::fassign(Fi, n_A, m_A, A, ldA, A_copy, m_A);

        
        r = PLUQ(Fi, FFLAS::FflasNonUnit, n_A, m_A, A_copy, m_A, P, Q);


        PL = FFLAS::fflas_new(Fi, n_A, r);
        UQ = FFLAS::fflas_new(Fi, r, m_A);
        ldPL = r;
        ldUQ = m_A;

        // extraction of U
        getTriangular<Field>(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit, n_A, m_A, r, A_copy, m_A, UQ, m_A, true);

        // extraction of L
        getTriangular<Field>(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, n_A, m_A, r, A_copy, m_A, PL, r, true);

        // apply the permutations (P on L and Q on U)
        applyP<Field>(Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, r, 0, n_A, PL, r, P);
        

        applyP<Field>(Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, r, 0, m_A, UQ, m_A, Q);

        FFLAS::fflas_delete(A_copy);
        FFLAS::fflas_delete(P,Q);
        memory_owner = true;
    }

    RRgen(const Field& Fi, size_t n_A, size_t m_A, typename Field::Element_ptr A, size_t ldA):n(n_A),m(m_A){

        size_t* P = FFLAS::fflas_new<size_t>(n_A);
        size_t* Q = FFLAS::fflas_new<size_t>(m_A);


        
        r = PLUQ(Fi, FFLAS::FflasNonUnit, n_A, m_A, A, m_A, P, Q);


        PL = FFLAS::fflas_new(Fi, n_A, r);
        UQ = FFLAS::fflas_new(Fi, r, m_A);
        ldPL = r;
        ldUQ = m_A;

        // extraction of U
        getTriangular<Field>(Fi, FFLAS::FflasUpper, FFLAS::FflasNonUnit, n_A, m_A, r, A, m_A, UQ, m_A, true);

        // extraction of L
        getTriangular<Field>(Fi, FFLAS::FflasLower, FFLAS::FflasUnit, n_A, m_A, r, A, m_A, PL, r, true);

        // apply the permutations (P on L and Q on U)
        applyP<Field>(Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, r, 0, n_A, PL, r, P);
        

        applyP<Field>(Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, r, 0, m_A, UQ, m_A, Q);

        FFLAS::fflas_delete(P,Q);
        memory_owner = true;
    }
    
    void RRExpand( const Field& Fi, typename Field::Element_ptr A, size_t ldA){
                
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n, m, r, Fi.one,
            PL, ldPL,
            UQ, ldUQ,
            Fi.zero, A, ldA);
    }

    RRgen<Field>* RRcopy(const Field& Fi){
        if(UQ){
            return new RRgen(Fi,n,m,r,PL,ldPL,UQ,ldUQ,true,true);
        }
        else {
            return new RRgen(Fi,n,m,r,PL,ldPL,nullptr,0,true,true);
        }
    }

    ~RRgen() {
        if (memory_owner){

            FFLAS::fflas_delete(PL);
            if (UQ){
                FFLAS::fflas_delete(UQ);
            }
        }
    }
};

/// @brief  RRRgen Class for the tree representation.                                                               A = [A11  A12]
///         If the node is a leaf the matrix is stored in LU_right->PL and LU_left, left and right are None.            [A21  A22]
template<class Field>
class RRRgen {

public:
    RRgen<Field>* LU_right;             // PLUQ representation of the top right submatrix (N2*rr)*(rr*N1) = A12
    RRgen<Field>* LU_left;              // PLUQ representation of the down left submatrix (N1*rl)*(rl*N2) = A21        
    size_t size_N1;                     // size of N/2 (N is the size of the matrix represented by the RRRgen)
    size_t size_N2;                     // size of N-N1
    size_t t;                           // threshold for recursive representation : when N1+N2 < t, the matrix is stored completly
    RRRgen<Field>* left;                // recursively pointing on the same representation of the up left submatrix = A11
    RRRgen<Field>* right;               // recursively pointing on the same representation of the down right submatrix = A22
    bool memory_owner;                  // whether the structure owns the memory or whether it is a view on memory area. Used in destructors


    RRRgen(RRgen<Field>* LU_right, RRgen<Field>* LU_left, size_t size_N1, size_t size_N2,size_t t, RRRgen* left, RRRgen* right,bool memory_owner)
        :  LU_right(LU_right), LU_left(LU_left), size_N1(size_N1), size_N2(size_N2),t(t), left(left), right(right),memory_owner(memory_owner){}
    
    
    RRRgen( const Field& F,typename Field::Element_ptr leaf, size_t n,size_t threshold,bool mem_owner,bool copy){

        LU_right = new RRgen<Field>(F,n,n,n,leaf,n,nullptr,0,mem_owner,copy);
        LU_left = nullptr;
        size_N1 = n;
        size_N2 = 0;
        t = threshold;
        left = nullptr;
        right = nullptr;
        memory_owner = true;
    }

    RRRgen( const Field& Fi,const size_t N, const size_t threshold,typename Field::Element_ptr A, const size_t lda, bool mem_owner, bool copy = false): t (threshold)
    {
        if (N <= threshold) { 
            // leaf
            LU_right = new RRgen<Field>(Fi,N,N,N,A,lda,nullptr,0,mem_owner,copy);
            LU_left = nullptr;
            size_N1 = N;
            size_N2 = 0;
            left = nullptr;
            right = nullptr;
            memory_owner = true;
            return;
        }

        size_N1 = N/2;
        size_N2 = N - size_N1;
        // cut the matrix A in four parts
        // in case of odd N, A12 and A21 are rectangles and not squares
        typename Field::Element_ptr A11 = A;
        typename Field::ConstElement_ptr A12 = A + size_N1;
        typename Field::ConstElement_ptr A21 = A + lda*size_N1;
        typename Field::Element_ptr A22 = (typename Field::Element_ptr)A21 + size_N1;

        ///////// PLUQ Factorisation for A12 and A21
        // std::cout << "A12 : N1 = "<< size_N1 << std::endl;
        // std::cout << "A12 : N2 = "<< size_N2 << std::endl;
        LU_right = new RRgen(Fi, size_N1,  size_N2, A12, lda);

        LU_left = new RRgen(Fi, size_N2,  size_N1, A21, lda);


        // recursion on A11 and A22   
        left = new RRRgen( Fi, size_N1, threshold, A11, lda, mem_owner,copy);

        right = new RRRgen( Fi, size_N2, threshold, A22, lda, mem_owner,copy);
        memory_owner = true;
        return;
    }

    
    RRRgen<Field>* RRRcopy(const Field& Fi){
        RRRgen<Field>* new_left;
        RRRgen<Field>* new_right;
        RRgen<Field>* new_LU_left;
        RRgen<Field>* new_LU_right;
        
        if (left){
            new_left = left->RRRcopy(Fi); 
        }
        else {
            new_left = nullptr;
        }
        
        if (right){
            new_right = right->RRRcopy(Fi); 
        }
        else {
            new_right = nullptr;
        }
        
        if (LU_left){
            new_LU_left = LU_left.RRcopy(Fi); 
        }
        else {
            new_LU_left = nullptr;
        }
        
        if (LU_right){
            new_LU_right = LU_right.RRcopy(Fi);
        }
        else {
            new_LU_right = nullptr;
        }

        return new RRRgen(new_LU_right,LU_left,size_N1,size_N2,t,new_left,new_right, true);
    }

    ~RRRgen() {
        if (memory_owner){
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
    }
};


/// @brief Computes the dense matrix of RRR(A) in B. B needs to be at least (A->sizeN1 + A->sizeN2) x (A->sizeN1 + A->sizeN2) preallocated.
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
        // leaf
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
            A->LU_right->PL, A->LU_right->ldPL,
            A->LU_right->UQ, A->LU_right->ldUQ,
            Fi.zero, B12, ldb);

        
        // B21 < L_l * U_l
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N2, N1, A->LU_left->r, Fi.one,
            A->LU_left->PL, A->LU_left->ldPL,
            A->LU_left->UQ, A->LU_left->ldUQ,
            Fi.zero, B21, ldb);

        
    }



/// @brief  Multiplies two matrices stored as rank revealing factorization. C = minus ? -1 : 1 A*B
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
        
        RRgen<Field>* RR_X = new RRgen(Fi,A->r,B->r,X,B->r); // X will be modified
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
            Fi.zero, R_C, m);
        
        
        RRgen<Field>* C = new RRgen(Fi,n,m,RR_X->r,L_C,RR_X->r,R_C,m,true); // can be modifed
        delete(RR_X);
        
        return  C;
    }

/// @brief multiplies two matrices stored as rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param A stored with an RRgen
/// @param B stored with an RRgen
template<class Field>
inline RRgen<Field>* RRxRR (const Field& Fi,RRgen<Field>* A, RRgen<Field>*B){
    return RRxRR (Fi, A, B, false);
}


/// @brief add two matrices stored as rank revealing factorization
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

/// @brief Adds a quasiseparable matrix in RRR representation and a rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param A        size n*n in RRR representation
/// @param B        size n*n in RR representation
template<class Field>
inline RRRgen<Field>* RRRaddRR (const Field& Fi,  RRRgen<Field>* A,  RRgen<Field>* B)
    {   
        // std::cout << "n et rank " <<A->size_N1+A->size_N2 << " " << A->t + B->r << std::endl;
        if (A->size_N1+A->size_N2 <= A->t + B->r){
            
            // C = RRRexpand(A)
            typename Field::Element_ptr C = FFLAS::fflas_new(Fi, A->size_N1+A->size_N2, A->size_N1+A->size_N2);
            RRRExpand(Fi, A, C, A->size_N1 + A->size_N2);
            
            // B_expanded = RRexpand(B)

            typename Field::Element_ptr B_expanded = FFLAS::fflas_new(Fi,B->n,B->m);
            B->RRExpand(Fi,B_expanded,B->m);

            // C = B+C
            FFLAS::faddin(Fi,A->size_N1+ A->size_N2,A->size_N1 +A->size_N2,B_expanded,B->n,C,A->size_N1 +A->size_N2);
            FFLAS::fflas_delete(B_expanded);
            RRRgen<Field>* D = new RRRgen(Fi,C,A->size_N1 + A->size_N2,A->t + B->r,true,true);
            FFLAS::fflas_delete(C);
            
            return D;
        }
        else {
            // B_expanded = [RR_B11 RR_B12]
            //              [RR_B21 RR_B22]

            RRgen<Field>* RR_B11 = new RRgen(Fi, A->size_N1, A->size_N1, B->r, B->PL, B->ldPL, B->UQ,B->ldUQ,false);
            RRgen<Field>* RR_B12 = new RRgen(Fi, A->size_N1, A->size_N2, B->r, B->PL, B->ldPL, B->UQ + A->size_N1,B->ldUQ,false);
            RRgen<Field>* RR_B21 = new RRgen(Fi, A->size_N2, A->size_N1, B->r, B->PL + (A->size_N1 * B->ldPL), B->ldPL, B->UQ ,B->ldUQ,false);
            RRgen<Field>* RR_B22 = new RRgen(Fi, A->size_N2, A->size_N2, B->r, B->PL + (A->size_N1 * B->ldPL), B->ldPL, B->UQ + A->size_N1,B->ldUQ,false);

            

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
            
            return new RRRgen(C12,C21,A->size_N1,A->size_N2,A->t + B->r,C11,C22,true);
        }
    }

/// @brief Multiplies a quasiseparable matric in RRR representation with a tall and skinny matrix
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
inline void RRRxTS (const Field& Fi, size_t n, size_t m_B,
     const RRRgen<Field>* A, typename Field::ConstElement_ptr B, size_t ldB,
     typename Field::Element_ptr C, size_t ldC)
    {
        if (n<=A->t+m_B){
            typename Field::Element_ptr Adense =FFLAS::fflas_new (Fi, n, n);
            RRRExpand(Fi, A, Adense, n);
            fgemm(Fi, 
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                n, m_B, n,
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

            typename Field::ConstElement_ptr B1 = B;
            typename Field::ConstElement_ptr B2 = B1 + N1*ldB;

            // C1 < RRRxTS(A11,B1)
            RRRxTS(Fi, N1, m_B, A->left, B1, ldB, C1, ldC);

            // C2 < RRRxTS(A22,B2)
            RRRxTS(Fi, N2, m_B, A->right, B2, ldB, C2, ldC);

            // X < RA12 x B2
            // X of size ru*m_B
            typename Field::Element_ptr X = FFLAS::fflas_new (Fi, A->LU_right->r, m_B);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                A->LU_right->r, m_B, N2,
                Fi.one, A->LU_right->UQ, A->LU_right->ldUQ,
                B2, ldB,
                0, X, m_B);

            // C1 < C1 + LA12 x X
            fgemm(Fi,
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                N1, m_B , A->LU_right->r ,
                Fi.one, A->LU_right->PL, A->LU_right->ldPL,
                X, m_B,
                Fi.one, C1, ldC);

            FFLAS::fflas_delete(X);

            // Y < RA21 x B1
            typename Field::Element_ptr Y = FFLAS::fflas_new (Fi, A->LU_left->r, m_B);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                A->LU_left->r, m_B, N1,
                Fi.one,A->LU_left->UQ, A->LU_left->ldUQ,
                B1, ldB,
                0, Y, m_B);

            // C2 < C2 + LA21 x Y
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                N2, m_B, A->LU_left->r,
                Fi.one, A->LU_left->PL, A->LU_left->ldPL,
                Y, m_B,
                Fi.one, C2, ldC);

            FFLAS::fflas_delete(Y);
            // C =   [C1]
            //       [C2]

        }
    }


/// @brief  Multiplies a tall and skinny matrix with a quasiseparable matric in RRR representation
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
     typename Field::Element_ptr C, size_t ldC)
    {
        if (n<=A->t+t){
            typename Field::Element_ptr Adense = FFLAS::fflas_new (Fi, n, n);
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

            typename Field::Element_ptr B1 = (typename Field::Element_ptr)B;
            typename Field::Element_ptr B2 = B1 + N1;

            // C1 < TSxRRR(A11,B1)
            TSxRRR(Fi, N1, t, B1, ldB, A->left, C1, ldC);

            // C2 < RRRxTS(A22,B2)
            TSxRRR(Fi, N2, t, B2, ldB, A->right, C2, ldC);
            

            // X < B2*A21->PL
            // X of size t*r
            typename Field::Element_ptr X = FFLAS::fflas_new (Fi, t, A->LU_left->r);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, A->LU_left->r, N2,
                Fi.one,B2, ldB ,
                A->LU_left->PL, A->LU_left->ldPL,
                0, X,A->LU_left->r);

            // C1 < C1 + X*A21->UQ 
            fgemm(Fi,
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, N1 , A->LU_left->r ,
                Fi.one, X, A->LU_left->r,
                A->LU_left->UQ, A->LU_left->ldUQ,
                Fi.one, C1, ldC);

            FFLAS::fflas_delete(X);

            // Y < B1*A12->PL
            // Y of size t*r
            typename Field::Element_ptr Y = FFLAS::fflas_new (Fi, t, A->LU_right->r);
            fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, A->LU_right->r, N1,
                Fi.one,B1, ldB ,
                A->LU_right->PL, A->LU_right->ldPL,
                0, Y,A->LU_right->r);

            // C2 < C2 + Y*A12->UQ 
            fgemm(Fi,
                FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                t, N2 , A->LU_right->r ,
                Fi.one, Y, A->LU_right->r,
                A->LU_right->UQ, A->LU_right->ldUQ,
                Fi.one, C2, ldC);

            FFLAS::fflas_delete(Y);


            // C =   [C1 C2]

        }
    }



/// @brief Multiplies a QS matrix in RRR representation with a rank revealing factorization
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
        if(B->m != n){
            std::cout << "PAS LE BON FORMAT " << std::endl;
            return nullptr;
        }

        if (B->UQ == nullptr){
            // if B is a leaf of the RRRgen then UQ is empty !
            // and the dense matrix is stored in PL
            typename Field::Element_ptr X = FFLAS::fflas_new (Fi, B->n, n);
            TSxRRR(Fi,n,B->n,B->PL,B->ldPL,A,X,n);
            RRgen<Field>* RR_X = new RRgen(Fi,  B->r,  n,  X,  n);
            FFLAS::fflas_delete(X);
            return RR_X;
        }

        typename Field::Element_ptr X = FFLAS::fflas_new (Fi, B->r, n);
        TSxRRR(Fi,n,B->r,B->UQ,B->ldUQ,A,X,n);

        // (LX, RX) < RRF(X)
        RRgen<Field>* RR_X = new RRgen(Fi,  B->r,  n,  X,  n);
        FFLAS::fflas_delete(X);
        
        // RD < RX
        typename Field::Element_ptr UQ_D = FFLAS::fflas_new (Fi, RR_X->r, n);
        FFLAS::fassign(Fi,RR_X->r,n,RR_X->UQ,RR_X->ldUQ,UQ_D,n);

        // LD < (minus ? (-1) : 1) LBxLX 
        typename Field::Element_ptr PL_D = FFLAS::fflas_new (Fi, B->n, RR_X->r);
        fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                B->n, RR_X->r, B->r,
                minus ? Fi.mOne : Fi.one , B->PL, B->ldPL,
                RR_X->PL, RR_X->ldPL,
                Fi.zero, PL_D, RR_X->r);
        

        RRgen<Field>* D = new RRgen(Fi,B->n,n,RR_X->r,PL_D,RR_X->ldPL,UQ_D,B->m,true);
        delete RR_X;
        return D;

    }

/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Computes C = A*B if minus == false , else computes C = -A*B
/// @tparam Field 
/// @param Fi
/// @param A        size n*n in a RRR representation     
/// @param B        size n*m in RR representation
/// @param minus    boolean in order to compute C = -A*B if true.
template<class Field>
inline RRgen<Field>* RRRxRR (const Field& Fi,const RRRgen<Field>* A,const RRgen<Field>* B,bool minus)
    {
        // X < RRRxTS(A,LB)
        size_t n = A->size_N1+A->size_N2;
        typename Field::Element_ptr X = FFLAS::fflas_new (Fi, n, B->r);
        RRRxTS(Fi,n,B->r,A,B->PL,B->ldPL,X,B->r);

        // (LX, RX) < RRF(X)
        RRgen<Field>* RR_X = new RRgen(Fi,  n,  B->r,  X,  B->r);
        FFLAS::fflas_delete(X);
        if (B->UQ == nullptr){
            return RR_X;
        }
        // LD < LX
        typename Field::Element_ptr PL_D = FFLAS::fflas_new (Fi, n, RR_X->r);
        FFLAS::fassign(Fi,n,RR_X->r,RR_X->PL,RR_X->ldPL,PL_D,RR_X->r);

        // RD < (minus ? (-1) : 1) RX x RB
        typename Field::Element_ptr UQ_D = FFLAS::fflas_new (Fi, RR_X->r, B->m);
        fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                RR_X->r, B->m, B->r,
                minus ? Fi.mOne : Fi.one , RR_X->UQ, RR_X->ldUQ,
                B->UQ, B->ldUQ,
                Fi.zero, UQ_D, B->m);
        

        RRgen<Field>* D = new RRgen(Fi,n,B->m,RR_X->r,PL_D,RR_X->ldPL,UQ_D,B->m,true);
        delete RR_X;
        return D;
    }


/// @brief Multiplies a QS matrix in RRR representation with a rank revealing factorization
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

/// @brief multiplies two QS matrices 
/// @tparam Field 
/// @param Fi 
/// @param A        size n*n
/// @param B        size n*n
template<class Field>
inline RRRgen<Field>* RRRxRRR (const Field& Fi, const RRRgen<Field>* A, const RRRgen<Field>* B){   
    size_t n = A->size_N1+A->size_N2;
    if (n<= A->t+B->t){
        //return RRRExpand(A) x RRRExpand(B)
        typename Field::Element_ptr A_expanded =FFLAS::fflas_new (Fi, n, n);
        typename Field::Element_ptr B_expanded =FFLAS::fflas_new (Fi, n, n);
        typename Field::Element_ptr C = FFLAS::fflas_new (Fi, n, n);
        RRRExpand(Fi,A,A_expanded,n);
        RRRExpand(Fi,B,B_expanded,n);
        fgemm(Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n, n, n,
            Fi.one, A_expanded, n,
            B_expanded, n,
            0, C, n);
        FFLAS::fflas_delete(A_expanded);
        FFLAS::fflas_delete(B_expanded);
        RRRgen<Field>* D = new RRRgen(Fi,n,A->t+B->t,C,n,true,true); 
        FFLAS::fflas_delete(C);
        return D;
    }

    else {
        // C11 < RRRxRRR(A11,B11)
        RRRgen<Field>* C11_ = RRRxRRR(Fi,A->left,B->left);

        // C22 < RRRxRRR(A22,B22)
        RRRgen<Field>* C22_ = RRRxRRR(Fi,A->right,B->right);

        // X < RRxRR(A12,B21)
        RRgen<Field>* X = RRxRR(Fi,A->LU_right,B->LU_left);

        // Y < RRxRR(A21,B12)
        RRgen<Field>* Y = RRxRR(Fi,A->LU_left,B->LU_right);

        // C11 < RRR+RR(C11,X)
        RRRgen<Field>* C11 = RRRaddRR(Fi,C11_,X);

        // C22 < RRR+RR(C22,Y)
        RRRgen<Field>* C22 = RRRaddRR(Fi,C22_,Y);

        delete C11_;
        delete C22_;
        delete X;
        delete Y;

        // LX < RRRxTS(A11,LB12) ; RX < RB12
        typename Field::Element_ptr LX =FFLAS::fflas_new (Fi, A->size_N1, B->LU_right->r);
        RRRxTS(Fi,A->size_N1,B->LU_right->r,A->left,B->LU_right->PL,B->LU_right->ldPL,LX,B->LU_right->r);
        
        typename Field::Element_ptr RX =FFLAS::fflas_new (Fi, B->LU_right->r, B->size_N2);
        FFLAS::fassign(Fi,B->LU_right->r,B->size_N2,
            B->LU_right->UQ,B->LU_right->ldUQ,
            RX,B->size_N2);

        // LY < LA12 ; RY < TSxRRR(RA12,B22) 
        typename Field::Element_ptr LY =FFLAS::fflas_new (Fi, A->size_N1, A->LU_right->r);
        FFLAS::fassign(Fi,A->size_N1,A->LU_right->r,
            A->LU_right->PL,A->LU_right->ldPL,
            LY,A->LU_right->r);
        
        typename Field::Element_ptr RY =FFLAS::fflas_new (Fi, A->LU_right->r, B->size_N2);
        TSxRRR(Fi,B->size_N2,A->LU_right->r,A->LU_right->UQ,A->LU_right->ldUQ,B->right,RY,B->size_N2);
        
        X = new RRgen(Fi, A->size_N1, B->size_N2,B->LU_right->r,LX,B->LU_right->r,RX,B->size_N2,true);
        Y = new RRgen(Fi, A->size_N1, B->size_N2,A->LU_right->r,LY,A->LU_right->r,RY,B->size_N2,true);

        // C12 < RR+RR(X,Y)
        RRgen<Field>* C12 = RRaddRR(Fi,X,Y);
        delete X;
        delete Y;

        // LX < RRRxTS(A22,LB21) ; RX <  RB21
        LX = FFLAS::fflas_new (Fi, A->size_N2, B->LU_left->r);
        RRRxTS(Fi,A->size_N2,B->LU_left->r,A->right,B->LU_left->PL,B->LU_left->ldPL,LX,B->LU_left->r);
        
        RX = FFLAS::fflas_new (Fi, B->LU_left->r, B->size_N1);
        FFLAS::fassign(Fi,B->LU_left->r,B->size_N1,
            B->LU_left->UQ,B->LU_left->ldUQ,
            RX,B->size_N1);

        // LY < LA21 ; RY < TSxRRR(RA21,B11)
        LY =FFLAS::fflas_new (Fi, A->size_N2, A->LU_left->r);
        FFLAS::fassign(Fi,A->size_N2,A->LU_left->r,
            A->LU_left->PL,A->LU_left->ldPL,
            LY,A->LU_left->r);
        
        RY = FFLAS::fflas_new (Fi, A->LU_left->r, B->size_N1);
        TSxRRR(Fi,B->size_N1,A->LU_left->r,A->LU_left->UQ,A->LU_left->ldUQ,B->left,RY,B->size_N1);
        
        X = new RRgen(Fi, A->size_N2, B->size_N1,B->LU_left->r,LX,B->LU_left->r,RX,B->size_N1,true);
        Y = new RRgen(Fi, A->size_N2, B->size_N1,A->LU_left->r,LY,A->LU_left->r,RY,B->size_N1,true);
        
        // C21 < RR+RR(X,Y)
        RRgen<Field>* C21 = RRaddRR(Fi,X,Y);
        delete X;
        delete Y;

        // RETURN C =   [C11    C12]
        //              [C12    C21]
        return new RRRgen(C12,C21,A->size_N1,B->size_N2,A->t+B->t,C11,C22,true);
    }

}


/// @brief Computes the inverse in RRR representation adn returns it in RRRgen.
/// @tparam Field 
/// @param Fi 
/// @param A in RRR representation
template<class Field>
inline RRRgen<Field>* RRRinvert (const Field& Fi,
            const RRRgen<Field>* A)
    {
        if (!A->left){
            // leaf
            // Y < RRRExpand(A)
            size_t N1 = A->size_N1;
            typename Field::Element_ptr Y =FFLAS::fflas_new (Fi, N1, N1);
            RRRExpand(Fi,A,Y,N1);
            // return Invert(Y)
            int nullity;
            FFPACK::Invert (Fi, N1,Y, N1, nullity);
            RRRgen<Field>* Y_RRR = new RRRgen(Fi,N1,A->t,Y,N1,true,true);
            FFLAS::fflas_delete(Y);

            return  Y_RRR;
        }

        // split the matrix as A =  [A11 A12] and   X = [X11 X12]
        //                          [A21 A22]           [X21 X22]

        // Y11 < RRRinvertrec(A11)
        RRRgen<Field>* Y11 = RRRinvert(Fi,A->left); //can be a leaf

        // Y12 < RRRxRR(Y11,A12)
        RRgen<Field>* Y12 = RRRxRR(Fi,Y11,A->LU_right);

        // Y21 < RRxRRR(A21,Y11)
        RRgen<Field>* Y21 = RRxRRR(Fi,Y11,A->LU_left);

        // Z < -RRxRR(A21,Y12)
        RRgen<Field>* Z = RRxRR(Fi,A->LU_left,Y12,true);

        // D < RRRaddRR(A22,Z)
        RRRgen<Field>* D = RRRaddRR(Fi,A->right,Z);
        delete Z;

        // X22 < RRRinvert(D)
        RRRgen<Field>* X22 = RRRinvert(Fi,D);
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
        delete Y11;
        delete Y12;
        delete Y21;
        
        // return X
        return new RRRgen(X12,X21,A->size_N1,A->size_N2,A->t,X11,X22,true);  
    }
    
    
    
/// @brief Computes the L U factorization in RRR representation and returns it in RRRgen with their inverse.
/// @tparam Field 
/// @param Fi 
/// @param A in RRR representation
/// @param L RRRgen uninitialized
/// @param U RRRgen uninitialized
/// @param L_inv RRRgen uninitialized
/// @param U_inv RRRgen uninitialized
template<class Field>
inline void LUfactRRRwInverse (const Field& Fi, const RRRgen<Field>* A, RRRgen<Field>*& L, RRRgen<Field>*& U, RRRgen<Field>*& L_inv, RRRgen<Field>*& U_inv)
{
    size_t size_N1 = A->size_N1;
    size_t size_N2 = A->size_N2;
    size_t N = size_N1 + size_N2;
    
    
    if (N <= A->t){
        // L/U = LU(Adense)
        RRgen<Field>* RR_A = new RRgen(Fi,N,N,(typename Field::ConstElement_ptr)A->LU_right->PL,A->LU_right->ldPL);
        typename Field::Element_ptr L_dense = FFLAS::fflas_new(Fi,N,N);
        FFLAS::fzero(Fi,N,N,L_dense,N);
        typename Field::Element_ptr U_dense = FFLAS::fflas_new(Fi,N,N);
        FFLAS::fzero(Fi,N,N,U_dense,N);
        
        FFLAS::fassign(Fi,N,RR_A->r,RR_A->PL,RR_A->ldPL,L_dense,N);
        FFLAS::fassign(Fi,RR_A->r,N,RR_A->UQ,RR_A->ldUQ,U_dense,N);
        FFLAS::WriteMatrix(std::cout << A->t << "Ldense =  " << std::endl, Fi, N, N, L_dense, N);
        FFLAS::WriteMatrix(std::cout << "Udense =  " << std::endl, Fi, N, N, U_dense, N);
        
        delete RR_A;
        
        L = new RRRgen(Fi, L_dense,N,A->t,true,false);
        U = new RRRgen(Fi, U_dense,N,A->t,true,false);
        // L_inv = inv(L)
        typename Field::Element_ptr Id1 = FFLAS::fflas_new(Fi,N,N);
        FFLAS::fidentity(Fi,N,N,Id1,N);
        FFLAS::ftrsm(Fi,FFLAS::FflasRight, FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit, N,N,Fi.one,L_dense,N,Id1,N);
        L_inv = new RRRgen(Fi, Id1,N,A->t,true,false);
        // U_inv = inv(U)
        typename Field::Element_ptr Id2 = FFLAS::fflas_new(Fi,N,N);
        FFLAS::fidentity(Fi,N,N,Id2,N);
        FFLAS::ftrsm(Fi,FFLAS::FflasRight, FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit, N,N,Fi.one,U_dense,N,Id2,N);
        U_inv = new RRRgen(Fi, Id2,N,A->t,true,false);
        return;
    }
    
    // L11/U11 = LUfactRRR(A11)
    RRRgen<Field>* L11 = nullptr;
    RRRgen<Field>* U11 = nullptr;
    RRRgen<Field>* L11_inv = nullptr;
    RRRgen<Field>* U11_inv = nullptr;
    LUfactRRRwInverse(Fi,A->left,L11,U11,L11_inv,U11_inv);
    
    // D1 = L21 / U21xU11^{-1}
    RRgen<Field>* D1 = A->LU_left->RRcopy(Fi);
    TSxRRR(Fi, size_N1, A->LU_left->r, A->LU_left->UQ, A->LU_left->ldUQ, U11_inv , D1->UQ , D1->ldUQ);
    
    // D2 = L11^{-1}*L12 / U12
    RRgen<Field>* D2 = A->LU_right->RRcopy(Fi);
    RRRxTS(Fi, size_N1, A->LU_right->r, L11_inv, A->LU_right->PL, A->LU_right->ldPL, D2->PL, D2->ldPL);

    typename Field::Element_ptr D1_expand = FFLAS::fflas_new(Fi,size_N1,size_N2);
    D2->RRExpand(Fi,D1_expand,size_N2);

    FFLAS::WriteMatrix(std::cout << "D2 =  " << std::endl, Fi, size_N1, size_N2, D1_expand, size_N2);


    
    // X22 = A22 - D1*D2
    RRgen<Field>* D1xD2 = RRxRR(Fi,D1,D2,true);
    RRRgen<Field>* X22 = RRRaddRR(Fi,A->right,D1xD2);
    delete D1xD2;
    
    
    
    // L2/U2 = LUfactRRR(X22)
    RRRgen<Field>* L22 = nullptr;
    RRRgen<Field>* U22 = nullptr;
    RRRgen<Field>* L22_inv = nullptr;
    RRRgen<Field>* U22_inv = nullptr;
    LUfactRRRwInverse(Fi,X22,L22,U22,L22_inv,U22_inv);
    
    //  L =     [L11  0 ]                                       U   =   [U11  D2]
    //          [D1  L22]                                               [ 0  U22]
    typename Field::Element_ptr L12_PL = FFLAS::fflas_new(Fi, size_N1, 1);
    FFLAS::fzero(Fi, size_N1, 1,L12_PL,1);
    typename Field::Element_ptr L12_UQ = FFLAS::fflas_new(Fi, 1, size_N2);
    FFLAS::fzero(Fi,1, size_N2,L12_UQ,size_N2);
    RRgen<Field>* L12 = new RRgen(Fi,size_N1,size_N2,1, L12_PL,1,L12_UQ,size_N2,true);
    L = new RRRgen(D1,L12,size_N1,size_N2,A->t,L11,L22,true);
    

    
    
    typename Field::Element_ptr U21_PL = FFLAS::fflas_new(Fi, size_N2, 1);
    FFLAS::fzero(Fi,size_N2, 1,U21_PL,1);
    typename Field::Element_ptr U21_UQ = FFLAS::fflas_new(Fi, 1, size_N1);
    FFLAS::fzero(Fi,1, size_N1,U21_UQ,size_N1);
    RRgen<Field>* U21 = new RRgen(Fi,size_N2,size_N1,1, U21_PL,1,U21_UQ,size_N1,true);
    U = new RRRgen(U21,D2,size_N1,size_N2,A->t,U11,U22,true);

    
    //  L_inv = [       L11_inv               0    ]            U_inv = [ U11_inv    -U11_inv*D2*U22_inv]
    //          [-L22_inv*D1*L11_inv      L22_inv  ]                    [    0              U22_inv     ]
    typename Field::Element_ptr L12_PL_inv = FFLAS::fflas_new(Fi, size_N1, 1);
    FFLAS::fzero(Fi, size_N1, 1,L12_PL_inv,1);
    typename Field::Element_ptr L12_UQ_inv = FFLAS::fflas_new(Fi, 1, size_N2);
    FFLAS::fzero(Fi,1, size_N2,L12_UQ_inv,size_N2);
    RRgen<Field>* L12_inv = new RRgen(Fi,size_N1,size_N2,1, L12_PL_inv,1,L12_UQ_inv,size_N2,true);
    RRgen<Field>* X = RRxRRR(Fi,L11_inv,D1);
    RRgen<Field>* D1_inv = RRRxRR(Fi,L22_inv,X,true);
    L_inv = new RRRgen(D1_inv,L12_inv,size_N1,size_N2,A->t,L11_inv,L22_inv,true);
    delete X;

    typename Field::Element_ptr U21_PL_inv = FFLAS::fflas_new(Fi, size_N2, 1);
    FFLAS::fzero(Fi,size_N2, 1,U21_PL_inv,1);
    typename Field::Element_ptr U21_UQ_inv = FFLAS::fflas_new(Fi, 1, size_N1);
    FFLAS::fzero(Fi,1, size_N1,U21_UQ_inv,size_N1);
    RRgen<Field>* U21_inv = new RRgen(Fi,size_N2,size_N1,1, U21_PL_inv,1,U21_UQ_inv,size_N1,true);
    RRgen<Field>* Y = RRxRRR(Fi,U22_inv,D2);
    RRgen<Field>* D2_inv = RRRxRR(Fi,U11_inv,Y,true);
    U_inv = new RRRgen(U21_inv,D2_inv,size_N1,size_N2,A->t,U11_inv,U22_inv,true);
    delete Y;
}

/// @brief Computes the L U factorization in RRR representation and returns it in RRRgen.
/// @tparam Field 
/// @param Fi 
/// @param A in RRR representation
/// @param L RRRgen uninitialized
/// @param U RRRgen uninitialized
template<class Field>
inline void LUfactRRR (const Field& Fi, const RRRgen<Field>* A, RRRgen<Field>*& L, RRRgen<Field>*& U)
{   
    size_t size_N1 = A->size_N1;
    size_t size_N2 = A->size_N2;
    size_t N = size_N1 + size_N2;
    if (N <= A->t){
        // L/U = LU(Adense)
        RRgen<Field>* RR_A = new RRgen(Fi,N,N,(typename Field::ConstElement_ptr)A->LU_right->PL,A->LU_right->ldPL);
        typename Field::Element_ptr L_dense = FFLAS::fflas_new(Fi,N,N);
        FFLAS::fzero(Fi,N,N,L_dense,N);
        typename Field::Element_ptr U_dense = FFLAS::fflas_new(Fi,N,N);
        FFLAS::fzero(Fi,N,N,U_dense,N);
        
        FFLAS::fassign(Fi,N,RR_A->r,RR_A->PL,RR_A->ldPL,L_dense,N);
        FFLAS::fassign(Fi,RR_A->r,N,RR_A->UQ,RR_A->ldUQ,U_dense,N);
        FFLAS::WriteMatrix(std::cout << A->t << "Ldense =  " << std::endl, Fi, N, N, L_dense, N);
        FFLAS::WriteMatrix(std::cout << "Udense =  " << std::endl, Fi, N, N, U_dense, N);
        delete RR_A;
        
        L = new RRRgen(Fi, L_dense,N,A->t,true,false);
        U = new RRRgen(Fi, U_dense,N,A->t,true,false);
        
        return;
    }
    
    // L11/U11 = LUfactRRR(A11)
    RRRgen<Field>* L11 = nullptr;
    RRRgen<Field>* U11 = nullptr;
    RRRgen<Field>* L11_inv = nullptr;
    RRRgen<Field>* U11_inv = nullptr;
    LUfactRRRwInverse(Fi,A->left,L11,U11,L11_inv,U11_inv);
    
    // D1 = L21 / U21xU11^{-1}
    RRgen<Field>* D1 = A->LU_left->RRcopy(Fi);
    TSxRRR(Fi, size_N1, A->LU_left->r, A->LU_left->UQ, A->LU_left->ldUQ, U11_inv , D1->UQ , D1->ldUQ);
    
    // D2 = L11^{-1}*L12 / U12
    RRgen<Field>* D2 = A->LU_right->RRcopy(Fi);
    RRRxTS(Fi, size_N1, A->LU_right->r, L11_inv, A->LU_right->PL, A->LU_right->ldPL, D2->PL, D2->ldPL);
    
    // X22 = A22 - D1*D2
    RRgen<Field>* D1xD2 = RRxRR(Fi,D1,D2,true);
    RRRgen<Field>* X22 = RRRaddRR(Fi,A->right,D1xD2);
    delete D1xD2;

    
    // L2/U2 = LUfactRRR(X22)
    RRRgen<Field>* L22 = nullptr;
    RRRgen<Field>* U22 = nullptr;
    LUfactRRR(Fi,X22,L22,U22);
    

    //  L = [L11   0]     U = [U11  D2]
    //      [D1  L22]         [0   U22]
    typename Field::Element_ptr L12_PL = FFLAS::fflas_new(Fi, size_N1, 1);
    FFLAS::fzero(Fi, size_N1, 1,L12_PL,1);
    typename Field::Element_ptr L12_UQ = FFLAS::fflas_new(Fi, 1, size_N2);
    FFLAS::fzero(Fi,1, size_N2,L12_UQ,size_N2);
    RRgen<Field>* L12 = new RRgen(Fi,size_N1,size_N2,1, L12_PL,1,L12_UQ,size_N2,true);
    L = new RRRgen(D1,L12,size_N1,size_N2,A->t,L11,L22,true);
    

    
    
    typename Field::Element_ptr U21_PL = FFLAS::fflas_new(Fi, size_N2, 1);
    FFLAS::fzero(Fi,size_N2, 1,U21_PL,1);
    typename Field::Element_ptr U21_UQ = FFLAS::fflas_new(Fi, 1, size_N1);
    FFLAS::fzero(Fi,1, size_N1,U21_UQ,size_N1);
    RRgen<Field>* U21 = new RRgen(Fi,size_N2,size_N1,1, U21_PL,1,U21_UQ,size_N1,true);
    U = new RRRgen(U21,D2,size_N1,size_N2,A->t,U11,U22,true);
}


    
    
}
    #endif //_FFPACK_ffpack_rrrgen_inl