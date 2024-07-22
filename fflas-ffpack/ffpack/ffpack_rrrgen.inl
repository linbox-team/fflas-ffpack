#ifndef __FFLASFFPACK_ffpack_rrrgen_inl
#define __FFLASFFPACK_ffpack_rrrgen_inl

#include <iostream>

namespace FFPACK{

/// @brief  node Class for the tree representation. 
///         If the node is a leaf the matrix is stored in U_u (size in size_N1) and both left and right are None.
template<class Field>
class Node {
public:
    typename Field::Element_ptr U_u;
    typename Field::Element_ptr L_u;
    size_t ru;
    typename Field::Element_ptr U_l;
    typename Field::Element_ptr L_l;
    size_t rl;
    size_t size_N1;
    size_t size_N2;
    Node* left;
    Node* right;

    Node(   typename Field::Element_ptr U_u, typename Field::Element_ptr L_u, size_t ru,
            typename Field::Element_ptr U_l, typename Field::Element_ptr L_l, size_t rl,
            size_t N1, size_t N2, Node* left, Node* right)
        : U_u(U_u), L_u(L_u), ru(ru), U_l(U_l), L_l(L_l), rl(rl), size_N1(N1), size_N2(N2), left(left), right(right) {}

    Node(typename Field::Element_ptr leaf, size_t N)
        : U_u(leaf), L_u(nullptr), ru(0), U_l(nullptr), L_l(nullptr), rl(0), size_N1(N), size_N2(0), left(nullptr), right(nullptr) {}

    ~Node() {
        FFLAS::fflas_delete(U_u);
        
        if (L_u) {
            FFLAS::fflas_delete(L_u);
        }

        if (U_l) {
            FFLAS::fflas_delete(U_l);
        }

        if (L_l) {
            FFLAS::fflas_delete(L_l);
        }
    }
};

/// @brief Class for the RRRrep ( a tree of nodes )
template<class Field>
class RRRrep {
private:
    Node<Field>* root;
    size_t lda;

    void inOrderTraversal(Node<Field>* node) const {

        if (node != nullptr) {
            // print left child
            inOrderTraversal(node->left);

            // print curr node U_u
            std::cout << "U_u : " << std::endl;
            for (size_t i = 0; i < node->ru; ++i) {
                for (size_t j = 0; j < node->size_N2; j++) {
                    std::cout << node->U_u[i*node->size_N2+j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;

            // print curr node L_u transpose
            std::cout << "L_u : " << std::endl;
            for (size_t i = 0; i < node->size_N1; ++i) {
                for (size_t j = 0; j < node->ru; j++){
                    std::cout << node->L_u[i*node->ru+j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            
            // print curr node U_l
            std::cout << "U_l : " << std::endl;
            for (size_t i = 0; i < node->rl; ++i) {
                for (size_t j = 0; j < node->size_N1; j++) {
                    std::cout << node->U_l[i*node->size_N1+j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;

            // print curr node L_l transpose
            std::cout << "L_l : ";
            for (size_t i = 0; i < node->size_N2; ++i) {
                for (size_t j = 0; j < node->rl; j++){
                    std::cout << node->L_l[i*node->rl+j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;

            //  print right child
            inOrderTraversal(node->right);
        }
    }

    void destroyTree(Node<Field>* node) {
        if (node != nullptr) {
            destroyTree(node->left);
            destroyTree(node->right);
            delete node;
        }
    }

public:
    RRRrep() : root(nullptr), lda(0) {}

    RRRrep(Node<Field>* root, size_t lda) : root(root), lda(lda) {}

    void setroot(Node<Field>* root) {
        root = root;
    }

    void setld(size_t ld) {
        lda = ld;
    }

    size_t getlda() const { return lda;}

    Node<Field>* getroot() const { return root;}

    ~RRRrep() {
        destroyTree(root);
    }

    void inOrderTraversal() const {
        inOrderTraversal(root);
        std::cout << std::endl;
    }
};

/// @brief RRR Generator recursive part
/// @tparam Field 
/// @param Fi 
/// @param N 
/// @param s 
/// @param A 
/// @param lda 
/// @return the root of the RRR representation of A
template<class Field>
inline Node<Field>* PLUQRRRGen_rec (const Field& Fi,
            const size_t N, const size_t s,
            typename Field::Element_ptr A, const size_t lda)
    {
        if (N/2 < s) {
            return new Node<Field>(A, N);
        }

        size_t N1 = N/2;
        size_t N2 = N - N1;

        // cut the matrix A in four parts
        // in case of odd N, A12 and A21 are rectangles and not squares
        typename Field::Element_ptr A11 = A;
        typename Field::Element_ptr A12 = A + N1;
        typename Field::Element_ptr A21 = A + lda*N1;
        typename Field::Element_ptr A22 = A21 + N1;

        //////////////// PLUQ Factorisation for A12
        size_t* P_u = FFLAS::fflas_new<size_t>(N1);
        size_t* Q_u = FFLAS::fflas_new<size_t>(N2);

        size_t r_u = PLUQ(Fi, FFLAS::FflasNonUnit,
                            N1, N2,
                            A12, lda,
                            P_u, Q_u);

        typename Field::Element_ptr U_u = FFLAS::fflas_new(Fi, r_u, N2);
        typename Field::Element_ptr L_u = FFLAS::fflas_new(Fi, N1, r_u);

        // extraction of U_u
        getTriangular<Field>(Fi, FFLAS::FflasUpper,
                        FFLAS::FflasNonUnit,
                        N1, N2, r_u,
                        A12, lda,
                        U_u, N2,
                        true);

        // extraction of L_u
        getTriangular<Field>(Fi, FFLAS::FflasLower,
                        FFLAS::FflasNonUnit,
                        N1, N2, r_u,
                        A12, lda,
                        L_u, N1,
                        true);

        // apply the permutations (P on L and Q on U)
        applyP<Field>(Fi,
                FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                r_u, 0, N1,
                L_u, N1, P_u);
        applyP<Field>(Fi,
                FFLAS::FflasRight, FFLAS::FflasNoTrans,
                r_u, 0, N2,
                U_u, N2, Q_u);

        //////////////// PLUQ factorisation for A21
        size_t* P_l = FFLAS::fflas_new<size_t>(N2);
        size_t* Q_l = FFLAS::fflas_new<size_t>(N1);

        size_t r_l = PLUQ(Fi, FFLAS::FflasNonUnit,
                            N2, N1,
                            A12, lda,
                            P_l, Q_l);

        typename Field::Element_ptr U_l = FFLAS::fflas_new(Fi, r_l, N1);
        typename Field::Element_ptr L_l = FFLAS::fflas_new(Fi, N2, r_l);


        // extraction of U_l
        getTriangular<Field>(Fi, FFLAS::FflasUpper,
                        FFLAS::FflasNonUnit,
                        N2, N1, r_l,
                        A21, lda,
                        U_l, N1,
                        true);

        // extraction of L_l
        getTriangular<Field>(Fi, FFLAS::FflasLower,
                        FFLAS::FflasNonUnit,
                        N2, N1, r_l,
                        A21, lda,
                        L_l, N2,
                        true);

        // apply the permutations (P on L and Q on U)  
        applyP<Field>(Fi,
                FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                r_l, 0, N2,
                L_l, N2, P_l);
        applyP<Field>(Fi,
                FFLAS::FflasRight, FFLAS::FflasNoTrans,
                r_l, 0, N1,
                U_l, N1, Q_l);

        FFLAS::fflas_delete(P_u,Q_u,P_l,Q_l);

        // recursion on A11 and A22        
        Node<Field>* left = PLUQRRRGen_rec(Fi,
                                    N1, s,
                                    A11, lda);
        Node<Field>* right = PLUQRRRGen_rec(Fi,
                                    N2, s,
                                    A22, lda);

        return new Node<Field>( U_u,
                                L_u,
                                r_u,
                                U_l,
                                L_l,
                                r_l,
                                N1, N2,
                                left, 
                                right);
    }

/// @brief RRR Generator API
/// @tparam Field 
/// @param Fi 
/// @param N 
/// @param s 
/// @param A 
/// @param lda 
/// @return the RRR representation of A
template<class Field>
inline RRRrep<Field>* PLUQRRRGen (const Field& Fi,
            const size_t N, const size_t s,
            typename Field::Element_ptr A, const size_t lda)
    {
        if (s == 0) {
            Node<Field>* root = new Node<Field>(A, N);
            RRRrep<Field>* RRRA = new RRRrep<Field>(root, lda);
            return RRRA;
        }

        if (N/2 < s) {
            std::cout << "Impoossible to generate an RRR representation, the given order of quasiseparability is too high" << std::endl;
            Node<Field>* root = PLUQRRRGen_rec(Fi, N, s, A, lda);

            return new RRRrep<Field>(root, lda);
        }

        else {
            Node<Field>* root = PLUQRRRGen_rec(Fi, N, s, A, lda);
            return new RRRrep<Field>(root, lda);
        }
    }

/// @brief (algo 4) Compute the dense matrix of RRR(A) in B recursive part
/// @tparam Field 
/// @param Fi 
/// @param nodeA 
/// @param B 
/// @param ldb 
template<class Field>
inline void RRRExpandrec (const Field& Fi,
            const Node<Field>& nodeA,
            typename Field::Element_ptr B, const size_t ldb)
    {
        if (nodeA.left == nullptr){
            size_t N = nodeA.size_N1;
            FFLAS::fassign(Fi, N, N, nodeA.U_u, N, B, ldb);
            return;
        }

        size_t N1 = nodeA.size_N1;
        size_t N2 = nodeA.size_N2;

        // cut the matrix A in four parts
        // in case of odd N, A12 and A21 are rectangles and not squares
        typename Field::Element_ptr B11 = B;
        typename Field::Element_ptr B12 = B + N1;
        typename Field::Element_ptr B21 = B + ldb*N1;
        typename Field::Element_ptr B22 = B21 + N1;

        // B11 < RRRExpand(A11)
        RRRExpandrec<Field>(Fi, *nodeA.left, B11, ldb);

        // B22 < RRRExpand(A22)
        RRRExpandrec<Field>(Fi, *nodeA.right, B22, ldb);

        // B12 < L_u * U_u
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N1, N2, nodeA.ru, 1,
            nodeA.L_u, nodeA.ru,
            nodeA.U_u, N2,
            0, B12, ldb);

        // B21 < L_l * U_l
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N2, N1, nodeA.rl, 1,
            nodeA.L_l, nodeA.rl,
            nodeA.U_l, N1,
            0, B21, ldb);
    }

/// @brief  (algo 4) Compute the dense matrix of RRR(A) in B
/// @tparam Field 
/// @param Fi 
/// @param RRRA 
/// @param B 
/// @param ldb 
template<class Field>
inline void RRRExpand (const Field& Fi,
            const RRRrep<Field>& RRRA,
            typename Field::Element_ptr B, const size_t ldb)
    {
        if ( RRRA.getlda() != ldb ){
            std::cout << "Impossible to generate an RRR representation, the given lda is not the same as the one in the RRRrep" << std::endl;
            return;
        }

        Node<Field>* root = RRRA.getroot();
        
        RRRExpandrec<Field>(Fi, *root, B, ldb);

        return;
        
    }

/// @brief (algo 5) multiplies two matrices stored as rank revealing factorization.
/// @tparam Field 
/// @param Fi 
/// @param r_A 
/// @param m
/// @param r_b 
/// @param n 
/// @param k 
/// @param LA       size m*r_A
/// @param UA       size r_A*k 
/// @param LB       size k*r_B
/// @param UB       size r_B*n 
/// @param lda      leading dimension of A
/// @param ldb      leading dimension of B
/// @param ldc      leading dimension of C
template<class Field>
inline void RRxRR (const Field& Fi,
            size_t r_A, size_t m, size_t r_b, size_t n, size_t k,
            typename Field::ConstElement_ptr LA, size_t ldLA, 
            typename Field::ConstElement_ptr UA, size_t ldUA,
            typename Field::ConstElement_ptr LB, size_t ldLB,
            typename Field::ConstElement_ptr UB, size_t ldUB)
    {
        // X < UA * LB

        // LX,RX < facto(X)

        // LC < LA*LX

        // RC < RX*RB
    }

/// @brief (algo 6) add two matrices stored as rank revealing factorization
/// @tparam Field 
/// @param Fi 
/// @param r_A 
/// @param r_b 
/// @param m 
/// @param n 
/// @param LA       size (m*r_A) transpose
/// @param UA       size r_A*n
/// @param LB       size (m*r_b) transpose
/// @param UB       size r_b*n
/// @param ldLA     leading dimension of LA
/// @param ldUA     leading dimension of UA
/// @param ldLB     leading dimension of LB
/// @param ldUB     leading dimension of UB
template<class Field>
inline void RRaddRR (const Field& Fi,
            size_t r_A, size_t r_b, size_t m, size_t n,
            typename Field::ConstElement_ptr LA, size_t ldLA,
            typename Field::ConstElement_ptr UA, size_t ldUA,
            typename Field::ConstElement_ptr LB, size_t ldLB,
            typename Field::ConstElement_ptr UB, size_t ldUB)
    {
        // X < [LA LB]

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
/// @param r_b      rank of b
/// @param n        dimension of A and B
/// @param A        size n*n in RRR representation
/// @param LB       size n*r_b
/// @param UB       size r_b*n 
template<class Field>
inline void RRRaddRR (const Field& Fi,
            size_t s, size_t r_b, size_t n,
            const RRRrep<Field>& A,
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
            const Node<Field>& nodeA,
            typename Field::ConstElement_ptr B, size_t ldB,
            typename Field::ConstElement_ptr C, size_t ldC)
    {
        if (n<=s+t){
            typename Field::Element_ptr Adense = fflas_new (Fi, n, n);
            RRRExpandrec(Fi, nodeA, Adense, n);
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

        size_t N1 = nodeA.size_N1;
        size_t N2 = nodeA.size_N2;

        // split the matrices as    [C1] = [A11 A12] [B1]
        //                          [C2] = [A12 A22] [B2]
        typename Field::Element_ptr C1 = C;
        typename Field::Element_ptr C2 = C1 + N1*ldC;

        typename Field::Element_ptr B1 = B;
        typename Field::Element_ptr B2 = B1 + N1*ldB;

        // C1 < RRRxTS(A11,B1)
        RRRxTSrec(Fi, s, N1, t, nodeA.left, B1, ldB, C1, ldC);

        // C2 < RRRxTS(A22,B2)
        RRRxTSrec(Fi, s, N2, t, nodeA.right, B2, ldB, C2, ldC);

        // X < RA12 x B2
        // X of size ru*t
        typename Field::Element_ptr X = fflas_new (Fi, nodeA.ru, t);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            nodeA.ru, N2, t,
            1,
            nodeA.U_u, N2,
            B2, ldB,
            0,
            X, t);

        // C1 < C1 + LA12 x X
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N1, nodeA.ru, t,
            1,
            nodeA.L_u, nodeA.ru,
            X, t,
            1,
            C1, ldC);

        // Y < RA21 x B1
        typename Field::Element_ptr Y = fflas_new (Fi, nodeA.rl, t);
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            nodeA.rl, N1, t,
            1,
            nodeA.U_l, N1,
            B1, ldB,
            0,
            Y, t);

        // C2 < C2 + LA21 x Y
        fgemm(Fi,
            FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            N2, nodeA.rl, t,
            1,
            nodeA.L_l, nodeA.rl,
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
            const RRRrep<Field>& A,
            typename Field::ConstElement_ptr B, size_t ldB,
            typename Field::ConstElement_ptr C, size_t ldC)
    {
        RRRxTSrec(Fi, s, n, t, A.getroot(), B, ldB, C, ldC);
    }

/// @brief (algo 9) Multiplies a QS matrix in RRR representation with a rank revealing factorization
///                 Compute C = A*B
/// @tparam Field 
/// @param Fi
/// @param s        order of QS of A
/// @param n   
/// @param r_b      rank of B
/// @param m  
/// @param A        size n*n in a RRR representation     
/// @param LB       size n*r_b
/// @param UB       size r_b*m
/// @param ldLB     leading dimension of LB
/// @param ldUB     leading dimension of UB
template<class Field>
inline void RRRxRR (const Field& Fi,
            size_t s, size_t n, size_t r_b, size_t m,
            const RRRrep<Field>& A,
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
            const RRRrep<Field>& A,
            const RRRrep<Field>& B)
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
            const RRRrep<Field>& A)
    {
        RRRinvertrec(Fi, A);
    }

template<class Field>
inline void RRRinvertrec (const Field& Fi,
            const Node<Field>& A)
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