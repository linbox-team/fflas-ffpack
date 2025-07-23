//-------------------------------------------------------------------------
//      Test suite for the Quasi-Separable matrices in RRR format
//-------------------------------------------------------------------------

/* Structure taken from test-quasisep.C */
#define __FFLASFFPACK_PLUQ_THRESHOLD 5 // Recursive vs iterative PLUQ threshold (default 256)

#include "fflas-ffpack/fflas-ffpack-config.h"
#include <givaro/modular-balanced.h>
#include <iostream>
#include <iomanip>

#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "fflas-ffpack/ffpack/ffpack_rrrgen.inl"
#include <cstdlib>
#include "fflas-ffpack/utils/args-parser.h"

#include <random>

using namespace std;
using namespace FFPACK;
using namespace FFLAS;

/**
 * \brief test equality between a dense matrix and the result of compressing with RR and reconstructing it.
 */
template<class Field>
bool test_compression_RR  (const Field & F, size_t n, size_t m,
                        typename Field::Element_ptr A, size_t lda)
{
    typename Field::Element_ptr Acheck = fflas_new (F, n, m);

    
    RRgen<Field>* RRA = new RRgen<Field>(F, n, m, (typename Field::ConstElement_ptr)A, lda);
    
    RRA->RRExpand(F, Acheck, m);

    bool ok = fequal (F, n, m, A, lda, Acheck, m);
    if ( !ok )
        {
            std::cout << "ERROR: different results for dense to RR and RR to dense (RRGen and RRExpand)"<<std::endl;
            WriteMatrix(std::cout << "Ainit = " << std::endl, F, n, m, A, lda);
            WriteMatrix(std::cout << "Acheck =  " << std::endl, F, n, m, Acheck, m);
        }
    
    FFLAS::fflas_delete(Acheck);
    delete(RRA);
    return ok;
}

/**
 * \brief test equality between a dense matrix and the result of compressing with RR and reconstructing it.
 */
template<class Field>
bool test_QS_order  (const Field & F, RRRgen<Field>* RRRA, size_t t)
{
    if (RRRA->left){
        if (RRRA->LU_left->r > t) {
            typename Field::Element_ptr A = fflas_new (F, RRRA->size_N2, RRRA->size_N1);
            RRRA->LU_left->RRExpand(F,A,RRRA->size_N1);
            WriteMatrix(std::cout << "A_current = " << std::endl, F, RRRA->size_N2, RRRA->size_N1, A, RRRA->size_N1);
            FFLAS::fflas_delete(A);  
            return false; 
        }    
        if (RRRA->LU_right->r > t) {
            typename Field::Element_ptr A = fflas_new (F, RRRA->size_N1, RRRA->size_N2);
            RRRA->LU_left->RRExpand(F,A,RRRA->size_N2);
            WriteMatrix(std::cout << "A_current = " << std::endl, F, RRRA->size_N1, RRRA->size_N2, A, RRRA->size_N2);
            FFLAS::fflas_delete(A); 
            return false;
        } 
        return test_QS_order(F,RRRA->left,t) &&  test_QS_order(F,RRRA->right,t);

    }
    else {
        //return true if the rank of RRRA->right->PL is lower than t
        //else return false
        //usually it's true as the dimension of the matrix is lower than t  
        return true;
    }
}


/**
 * \brief test equality between a dense matrix and the result of compressing with RRR and reconstructing it.
 */
template<class Field>
bool test_compression  (const Field & F, size_t n, size_t t,
                        typename Field::Element_ptr A, size_t lda)
{   
    typename Field::Element_ptr Acheck = fflas_new (F, n, n);
    
    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n, t, A, n,true,true);
    
    RRRExpand<Field>(F, RRRA, Acheck, n);

    bool ok = fequal (F, n, n, A, n, Acheck, n);
    if ( !ok )
        {
            std::cout << "ERROR: different results for dense to RRR and RRR to dense (RRRGen and Expand)"<<std::endl;
            WriteMatrix(std::cout << "Ainit = " << std::endl, F, n, n, A, lda);
            WriteMatrix(std::cout << "Acheck =  " << std::endl, F, n, n, Acheck, n);
            WriteMatrix(std::cout << "U_u = " << std::endl, F, RRRA->LU_right->r, RRRA->size_N2, RRRA->LU_right->UQ, RRRA->size_N2);
            WriteMatrix(std::cout << "L_u = " << std::endl, F, RRRA->size_N1, RRRA->LU_right->r, RRRA->LU_right->PL, RRRA->LU_right->r);
            WriteMatrix(std::cout << "U_l = " << std::endl, F, RRRA->LU_left->r, RRRA->size_N1, RRRA->LU_left->UQ, RRRA->size_N1);
            WriteMatrix(std::cout << "L_l = " << std::endl, F, RRRA->size_N2, RRRA->LU_left->r, RRRA->LU_left->PL, RRRA->LU_left->r);
        }
    
    ok = ok && test_QS_order(F,RRRA,t);
    if ( !ok )
        {
            std::cout << "ERROR: QS order too high"<<std::endl;
            WriteMatrix(std::cout << "Ainit = " << std::endl, F, n, n, A, lda);
        }


    FFLAS::fflas_delete(Acheck);
    delete(RRRA);
    return ok;
}

// /**
//  * \brief test equality between a C = A*B with fgemm and C = A*B with RRgen.
//  */
template<class Field>
bool test_RRxRR  (const Field & F, size_t n_A,size_t n_B, size_t m_B, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = A*B
    typename Field::Element_ptr C_init = fflas_new (F, n_A, m_B);

    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n_A, m_B,n_B ,
            F.one, A, lda,
            B, ldb,
            0, C_init, m_B);
    
    typename Field::Element_ptr C_check = fflas_new (F, n_A, m_B);
    
    RRgen<Field>* RRA = new RRgen(F, n_A, n_B, (typename Field::ConstElement_ptr)A, lda);
    RRgen<Field>* RRB = new RRgen(F, n_B, m_B, (typename Field::ConstElement_ptr)B, ldb);
    RRgen<Field>* RRC = RRxRR(F,RRA,RRB);

    RRC->RRExpand(F,C_check,m_B);

    bool ok = fequal (F, n_A, m_B, C_init, m_B, C_check, m_B);
    if ( !ok )
        {
            std::cout << "ERROR: different results for RRxRR product and fgemm product"<<std::endl;
            WriteMatrix(std::cout << "A*B = " << std::endl, F, n_A, m_B, C_init, m_B);            
            WriteMatrix(std::cout << "A*B = " << std::endl, F, n_A, m_B, C_check, m_B);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRA);
    delete(RRB);
    delete(RRC);
    return ok;
}

// /**
//  * \brief test equality between a C = A+B with fadd and C = A+B with RRgen.
//  */
template<class Field>
bool test_RRaddRR  (const Field & F, size_t n_A, size_t m_A, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = A+B
    typename Field::Element_ptr C_init = fflas_new (F, n_A, m_A);
    typename Field::Element_ptr C_check = fflas_new (F, n_A, m_A);

    fadd(F, 
            n_A, m_A,
            A, lda,
            B, ldb,
            C_init, m_A);
    
    RRgen<Field>* RRA = new RRgen(F, n_A, m_A, (typename Field::ConstElement_ptr)A, lda);
    RRgen<Field>* RRB = new RRgen(F, n_A, m_A, (typename Field::ConstElement_ptr)B, ldb);
    RRgen<Field>* RRC = RRaddRR(F,RRA,RRB);
    RRC->RRExpand(F,C_check,m_A);

    bool ok = fequal (F, n_A, m_A, C_init, m_A, C_check, m_A);
    if ( !ok )
        {
            std::cout << "ERROR: different results for RRaddRR  and fadd "<<std::endl;
            WriteMatrix(std::cout << "fadd : A+B = " << std::endl, F, n_A, m_A, C_init, m_A);            
            WriteMatrix(std::cout << "RRaddRR : A+B = " << std::endl, F, n_A, m_A, C_check, m_A);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRA);
    delete(RRB);
    delete(RRC);
    return ok;
}

// /**
//  * \brief test equality between a C = A+B with fadd and C = A+B with RRRaddRR.
//  */
template<class Field>
bool test_RRRaddRR  (const Field & F, size_t n_A, size_t m_A, size_t t, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = A+B
    typename Field::Element_ptr C_init = fflas_new (F, n_A, m_A);
    typename Field::Element_ptr C_check = fflas_new (F, n_A, m_A);

    fadd(F, 
            n_A, m_A,
            A, lda,
            B, ldb,
            C_init, m_A);

    // std::cout << "A->t = "<< t<<std::endl;
    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n_A, t, A, lda,true,true);
    RRgen<Field>* RRB = new RRgen(F, n_A, m_A, (typename Field::ConstElement_ptr)B, ldb);
    // std::cout << "B->r = "<< RRB->r<<std::endl;
    // WriteMatrix(std::cout << "PL of B " << std::endl, F, n_A, RRB->r, RRB->PL, RRB->ldPL);
    RRRgen<Field>* RRRC = RRRaddRR(F,RRRA,RRB);

    RRRExpand<Field>(F, RRRC, C_check, n_A);

    bool ok = fequal (F, n_A, m_A, C_init, m_A, C_check, m_A);
    if ( !ok )
        {
            std::cout << "ERROR: different results for RRRaddRR  and fadd "<<std::endl;
            WriteMatrix(std::cout << "fadd : A+B = " << std::endl, F, n_A, m_A, C_init, m_A);            
            WriteMatrix(std::cout << "RRRaddRR : A+B = " << std::endl, F, n_A, m_A, C_check, m_A);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRRA);
    delete(RRB);
    delete(RRRC);
    return ok;
}

// /**
//  * \brief test equality between a C = A*B with fgemm and C = A*B with RRRxTS.
//  */
template<class Field>
bool test_RRRxTS  (const Field & F, size_t n_A, size_t t, size_t m_B, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = A*B
    typename Field::Element_ptr C_init = fflas_new (F, n_A, m_B);
    typename Field::Element_ptr C_check = fflas_new (F, n_A, m_B);

    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n_A, m_B,n_A ,
            F.one, A, lda,
            B, ldb,
            0, C_init, m_B);

    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n_A, t, A, lda,true,true);

    RRRxTS(F,n_A,m_B,RRRA,B,ldb,C_check,m_B);

    bool ok = fequal (F, n_A, m_B, C_init, m_B, C_check, m_B);
    if ( !ok )
        {
            std::cout << "ERROR: different results for RRRxTS  and fgemm "<<std::endl;
            WriteMatrix(std::cout << "fgemm : A*B = " << std::endl, F, n_A, m_B, C_init, m_B);            
            WriteMatrix(std::cout << "RRRxTS : A*B = " << std::endl, F, n_A, m_B, C_check, m_B);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRRA);

    return ok;
}

// /**
//  * \brief test equality between a C = A*B with fgemm and C = A*B with TSxRRR.
//  */
template<class Field>
bool test_TSxRRR  (const Field & F, size_t n_A, size_t t, size_t n_B, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = B*A
    typename Field::Element_ptr C_init = fflas_new (F, n_B, n_A);
    typename Field::Element_ptr C_check = fflas_new (F, n_B, n_A);

    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n_B, n_A,n_A,
            F.one, B, ldb,
            A, lda,
            0, C_init, n_A);

    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n_A, t, A, lda,true,true);
    TSxRRR(F,n_A,n_B,B,ldb,RRRA,C_check,n_A);

    bool ok = fequal (F, n_B, n_A, C_init, n_A, C_check, n_A);
    if ( !ok )
        {
            std::cout << "ERROR: different results for TSxRRR and fgemm "<<std::endl;
            WriteMatrix(std::cout << "fgemm : A*B = " << std::endl, F, n_B, n_A, C_init, n_A);            
            WriteMatrix(std::cout << "TSxRRR : A*B = " << std::endl, F, n_B, n_A, C_check, n_A);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRRA);

    return ok;
}

// /**
//  * \brief test equality between a C = B*A with fgemm and C = B*A with RRxRRR and RRRxRR.
//  */
template<class Field>
bool test_RRxRRR  (const Field & F, size_t n_A, size_t t, size_t n_B, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = B*A
    typename Field::Element_ptr C_init = fflas_new (F, n_B, n_A);
    typename Field::Element_ptr C_check = fflas_new (F, n_B, n_A);

    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n_B, n_A,n_A,
            F.one, B, ldb,
            A, lda,
            0, C_init, n_A);
    const RRRgen<Field>* RRRA = new RRRgen<Field>(F, n_A, t, A, lda,true,true);
    const RRgen<Field>* RRB = new RRgen(F, n_B, n_A, (typename Field::ConstElement_ptr)B, ldb);
    RRgen<Field>* RRC = RRxRRR(F,RRRA,RRB);
    RRC->RRExpand(F,C_check,n_A);

    bool ok = fequal (F, n_B, n_A, C_init, n_A, C_check, n_A);
    if ( !ok )
        {
            std::cout << "ERROR: different results for RRxRRR and fgemm "<<std::endl;
            WriteMatrix(std::cout << "fgemm : A*B = " << std::endl, F, n_B, n_A, C_init, n_A);            
            WriteMatrix(std::cout << "RRxRRR : A*B = " << std::endl, F, n_B, n_A, C_check, n_A);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRRA);
    delete(RRB);
    delete(RRC);

    return ok;
}

// /**
//  * \brief test equality between a C = A*B with fgemm and C = A*B with RRRxRR. A = n*n and B = n*m
//  */
template<class Field>
bool test_RRRxRR  (const Field & F, size_t n_A, size_t t, size_t m_B, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = A*B
    typename Field::Element_ptr C_init = fflas_new (F, n_A, m_B);
    typename Field::Element_ptr C_check = fflas_new (F, n_A, m_B);

    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n_A, m_B,n_A,
            F.one, A, lda,
            B, ldb,
            0, C_init, m_B);
    const RRRgen<Field>* RRRA = new RRRgen<Field>(F, n_A, t, A, lda,true,true);
    const RRgen<Field>* RRB = new RRgen(F, n_A, m_B, (typename Field::ConstElement_ptr)B, ldb);
    RRgen<Field>* RRC = RRRxRR(F,RRRA,RRB);
    RRC->RRExpand(F,C_check,m_B);

    bool ok = fequal (F, n_A, m_B, C_init, m_B, C_check, m_B);
    if ( !ok )
        {
            std::cout << "ERROR: different results for RRRxRR and fgemm "<<std::endl;
            WriteMatrix(std::cout << "fgemm : A*B = " << std::endl, F, n_A, m_B, C_init, m_B);            
            WriteMatrix(std::cout << "RRRxRR : A*B = " << std::endl, F, n_A, m_B, C_check, m_B);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRRA);
    delete(RRB);
    delete(RRC);

    return ok;
}

// /**
//  * \brief test equality between a C = A*B with fgemm and C = A*B with RRRxRRR. A = n*n and B = n*n
//  */
template<class Field>
bool test_RRRxRRR  (const Field & F, size_t n, size_t t_A, size_t t_B, 
                        typename Field::Element_ptr A, size_t lda,
                        typename Field::Element_ptr B, size_t ldb)
{
    // C_check = A*B
    typename Field::Element_ptr C_init = fflas_new (F, n, n );
    typename Field::Element_ptr C_check = fflas_new (F, n, n);

    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
            n, n,n,
            F.one, A, lda,
            B, ldb,
            0, C_init, n);

    const RRRgen<Field>* RRRA = new RRRgen<Field>(F, n, t_A, A, lda,true,true);
    const RRRgen<Field>* RRRB = new RRRgen<Field>(F, n, t_B, B, ldb,true,true);

    
    RRRgen<Field>* RRRC = RRRxRRR(F,RRRA,RRRB);
    RRRExpand(F,RRRC,C_check,n);

    bool ok = fequal (F, n, n, C_init, n, C_check, n);
    if ( !ok )
        {
            std::cout << "ERROR: different results for RRRxRRR and fgemm "<<std::endl;
            WriteMatrix(std::cout << "fgemm : A*B = " << std::endl, F, n, n, C_init, n);            
            WriteMatrix(std::cout << "RRRxRRR : A*B = " << std::endl, F, n, n, C_check, n);
        }
    FFLAS::fflas_delete(C_init);
    FFLAS::fflas_delete(C_check);
    delete(RRRA);
    delete(RRRB);
    delete(RRRC);

    return ok;
}


/**
 * \brief test equality between a invert matrix and the result of compressing and invert it with RRR and reconstructing it.
 */
template<class Field>
bool test_invert  (const Field & F, size_t n, size_t t,
                        typename Field::Element_ptr A, size_t lda)
{   
    // Ainit = A⁻1 with FFPACK invert
    typename Field::Element_ptr Ainit = fflas_new (F, n, n);
    FFLAS::fassign(F,n,n,A,lda,Ainit,n);
    int nullity;
    FFPACK::Invert (F, n ,Ainit, n, nullity);
    // Acheck = A^-1 with RRRinvert
    typename Field::Element_ptr Acheck = fflas_new (F, n, n);
    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n, t, A, n,true,true);
    RRRgen<Field>* RRRA_invert = nullptr;
    try{
        RRRA_invert = RRRinvert(F,RRRA);
    }
    catch(...){
        FFLAS::fflas_delete(Ainit);
        FFLAS::fflas_delete(Acheck);
        delete RRRA_invert;
        delete RRRA;
        return false;
    }
    RRRExpand<Field>(F, RRRA_invert, Acheck, n);

    // I = Acheck x A
    typename Field::Element_ptr In = fflas_new (F, n, n);
    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
        n, n,n,
        F.one, Acheck, lda,
        A, n,
        0, In, n);


    bool ok = fequal (F, n, n, Ainit, n, Acheck, n);
    if ( !ok )
        {
            std::cout << "ERROR: different results for dense A⁻1 and RRR A⁻1"<<std::endl;
            WriteMatrix(std::cout << "Ainit = " << std::endl, F, n, n, Ainit, n);
            WriteMatrix(std::cout << "Acheck =  " << std::endl, F, n, n, Acheck, n);
            WriteMatrix(std::cout << "In =  " << std::endl, F, n, n, In, n);
        }
    FFLAS::fflas_delete(Ainit);
    FFLAS::fflas_delete(Acheck);
    FFLAS::fflas_delete(In);
    delete(RRRA);
    delete(RRRA_invert);

    return ok;
}

/**
 * \brief test equality between a FTRSM and the result of compressing and FRSTM it with RRR.
 */
template<class Field>
bool test_FRTSM_RRR  (const Field & F, size_t n, size_t t, typename Field::Element_ptr A, size_t lda, typename Field::Element_ptr B, size_t ldb)
{   
    // L/U init = RR(A)
    // RRgen<Field>* LU_A = new RRgen(F,n,n,(typename Field::ConstElement_ptr)A,n);
    RRgen<Field>* RRB = new RRgen(F, n, n, (typename Field::ConstElement_ptr)B, ldb);
    // WriteMatrix(std::cout << "B =  " << std::endl, F, RRB->n, RRB->r, RRB->PL, RRB->ldPL);
    
    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n, t, A, n,true,true);
    RRRgen<Field>* RRRL = nullptr;
    RRRgen<Field>* RRRU = nullptr;
    
    LUfactRRR (F, RRRA, RRRL, RRRU);

    TRSM_RRR_RR (F, FflasRight, FflasLower, (const RRRgen<Field>*) RRRU, RRB);
    // TRSM_RRR_RR (Fi, FflasLeft, FflasLower, RRRU, RRB);


    

    bool ok = true;
    // ok = ok && fequal (F, n, n, Ucheck, n, LU_A->UQ, n);
    if ( !ok )
        {
            std::cout << "ERROR: different results for dense LU and RRR LU"<<std::endl;
            // WriteMatrix(std::cout << "Lcheck =  " << std::endl, F, n, n, Lcheck, n);
            // WriteMatrix(std::cout << "Ucheck =  " << std::endl, F, n, n, Ucheck, n);
        }

    
    delete(RRRA);
    delete(RRB);
    delete(RRRU);
    delete(RRRL);

    return ok;
}


/**
 * \brief test equality between a LU factorization and the result of compressing and factorization it with RRR.
 */
template<class Field>
bool test_LU_RRR  (const Field & F, size_t n, size_t t, typename Field::Element_ptr A, size_t lda)
{   
    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n, t, A, n,true,true);
    RRRgen<Field>* RRRL = nullptr;
    RRRgen<Field>* RRRU = nullptr;
    try{
        LUfactRRR (F, RRRA, RRRL, RRRU);
    }
    catch (...){
        delete RRRA;
        delete RRRL;
        delete RRRU;
        cout<< endl <<"ERROR : Can't compute RRR LU, A does not satisfy full rank on sub matrix"<<endl;
        return false;
    }


    // L/U init = RR(A)
    RRgen<Field>* LU_A = new RRgen(F,n,n,(typename Field::ConstElement_ptr)A,n);
    

    // L/U check = L/U(A) with RRR_LU
    typename Field::Element_ptr Lcheck = fflas_new (F, n, n);
    typename Field::Element_ptr Ucheck = fflas_new (F, n, n);
    typename Field::Element_ptr Acheck = fflas_new (F, n, n);

    RRRExpand<Field>(F, RRRL, Lcheck, n);
    RRRExpand<Field>(F, RRRU, Ucheck, n);
    
    fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
        n, n,n,
        F.one, Lcheck, n,
        Ucheck, n,
        0, Acheck, n);
        
    bool ok = fequal (F, n, n, Acheck, n, A, n);
    if ( !ok )
    {
        std::cout << "ERROR: different results for dense LU and RRR LU"<<std::endl;
        WriteMatrix(std::cout << "Linit = " << std::endl, F, n, n, LU_A->PL, n);
        WriteMatrix(std::cout << "Lcheck =  " << std::endl, F, n, n, Lcheck, n);
        WriteMatrix(std::cout << "Uinit = " << std::endl, F, n, n, LU_A->UQ, n);
        WriteMatrix(std::cout << "Ucheck =  " << std::endl, F, n, n, Ucheck, n);
    }
    FFLAS::fflas_delete(Lcheck);
    FFLAS::fflas_delete(Ucheck);
    FFLAS::fflas_delete(Acheck);
    delete(LU_A);
    delete(RRRA);
    delete(RRRU);
    delete(RRRL);

    return ok;
}

template<class Field>
bool launch_instance_check (const Field& F, size_t n, size_t t, size_t m, size_t r, typename Field::RandIter& G)
{
    typedef typename Field::Element_ptr Element_ptr;
    bool ok = true;
    
    // Generation of a dense QS matrix with RandomLTQSMatrixWithRankandQSorder
    Element_ptr T = fflas_new (F, n, n);
    Element_ptr A2 = fflas_new (F, n, n);
    RandomLTQSMatrixWithRankandQSorder (F,n,r,t,T, n,G);
    RandomLTQSMatrixWithRankandQSorder (F,n,r,t, A2, n,G);
    size_t * p = FFLAS::fflas_new<size_t> (ceil(n/2.));
    
    for (size_t i = 0; i < ceil(n/2.); i++){
        p[i] = n - i - 1;
    }
    applyP (F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), T, n, p);
    applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, ceil(n/2.), A2, n, p);
    Givaro::GeneralRingNonZeroRandIter<Field,typename Field::RandIter> nzG (G);
    for (size_t i=0; i< n; ++i)
            nzG.random (A2[i*(n+1)]);
    faddin (F, n, n, A2, n, T, n);
    FFLAS::fflas_delete(A2, p);
    // WriteMatrix(std::cout << "A =  " << std::endl, F, n, n, T, n);

    //test operations with random matrices
    Element_ptr A = fflas_new (F, n, n); // n*n random matrix
    FFLAS::frand(F,G,n,n,A,n);

    Element_ptr B = fflas_new (F, n, n); // n*n random matrix
    FFPACK::RandomMatrixWithRank(F,n,n,r,B,n,G);

    Element_ptr C = fflas_new (F, n, m); // n*m random matrix
    FFLAS::frand(F,G,n,m,C,m);

    Element_ptr D = fflas_new (F, m, n); // m*n random matrix
    FFLAS::frand(F,G,m,n,D,n);

    Element_ptr E = fflas_new (F, n, n); // n*n random matrix with t-QSorder and rank r
    FFPACK::RandomLTQSMatrixWithRankandQSorder (F,n,r,t,E,n,G);

    Element_ptr M = fflas_new (F, n, n); // n*n random matrix with t-QSorder and rank r
    FFPACK::RandomLTQSMatrixWithRankandQSorder (F,n,r,t,M,n,G);

    Element_ptr N = fflas_new (F, n, n); // n*n random matrix with t-QSorder and rank n so it can be inverted
    // FFPACK::RandomLTQSMatrixWithRankandQSorder (F,n,n,t,N,n,G);
    
    // ok = ok && test_compression_RR(F,n,m,C,m);
    // ok = ok && test_compression(F,n,t,T,n);
    // ok = ok && test_RRxRR(F,n,n,n,A,n,B,n);
    // ok = ok && test_RRaddRR(F,n,n,A,n,B,n);
    // ok = ok && test_RRRaddRR(F,n,n,t,E,n,B,n);
    // ok = ok && test_RRRxTS(F,n,t,m,E,n,C,m);
    // ok = ok && test_TSxRRR(F,n,t,m,E,n,D,n);
    // ok = ok && test_RRxRRR(F,n,t,m,E,n,D,n);
    // ok = ok && test_RRRxRR(F,n,t,m,E,n,C,m);
    // ok = ok && test_RRRxRRR(F,n,t,t,E,n,M,n);
    // ok = ok && test_invert(F,n,t,A,n);
    // ok = ok && test_FRTSM_RRR(F,n,t,A,n,B,n); 
    ok = ok && test_LU_RRR(F,n,t,T,n); 


    
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);
    FFLAS::fflas_delete(C);
    FFLAS::fflas_delete(D);
    FFLAS::fflas_delete(E);
    FFLAS::fflas_delete(M);
    FFLAS::fflas_delete(N);
    FFLAS::fflas_delete(T);


    return ok;
}

template<class Field>
bool run_with_field(Givaro::Integer q, uint64_t b, size_t n, size_t t, size_t m, size_t r, size_t iters, uint64_t seed)
{
    bool ok = true ;
    while (ok && iters)
        {
            /* New field 
             * chooseField returns a pointer, F needs to be passed by its value */
            Field* F= chooseField<Field>(q,b,seed);
            if (F==nullptr)
                return true;
            /* Initiate random number generator */
            typename Field::RandIter G(*F,seed++);
            srandom(seed);
            
            std::ostringstream oss;
            F->write(oss);
            
            std::cout.fill('.');
            std::cout<<"Checking ";
            std::cout.width(117);
            std::cout<<oss.str();
            std::cout<<" ... ";
            ok = ok && launch_instance_check (*F, n, t, m, r, G);

            if (ok)
                std::cout << "PASSED "<<std::endl;
            

            delete(F);
            iters--;
        }
    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(20); // In order to print integers as integers even on float types, could be done once for all fflas
    Givaro::Integer q=-1;
    size_t b=0;
    size_t n=23;
    size_t t=5;
    size_t m=10;
    size_t r = 10;
    int iters=3;
    bool loop=false;
    uint64_t seed = getSeed();

    Argument as[] = {
                     { 'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER , &q },
                     { 'b', "-b B", "Set the bitsize of the field characteristic.", TYPE_INT , &b },
                     { 'n', "-n N", "Set the matrix row and column dimension.", TYPE_INT , &n },
                     { 't', "-t T", "Set the order of quasi-separability.", TYPE_INT , &t },
                     { 'r', "-r R", "Set the matrix rank when generated with RandomLTQSMatrixWithRankandQSorder.", TYPE_INT , &r },
                     { 'm', "-m M", "Set the col dim of the Tall and Skinny matrix.", TYPE_INT , &m },
                     { 'i', "-i R", "Set number of repetitions.", TYPE_INT, &iters },
                     { 'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL, &loop },
                     { 's', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed },
                     END_OF_ARGUMENTS
    };

    parseArguments(argc,argv,as);
   
    bool ok=true;
    do{
        std::cerr<<"with seed = "<<seed<<std::endl;
        std::cerr<<"Random matrices tests"<<std::endl;
        ok = ok &&run_with_field<Givaro::Modular<float> >           (q,b,n,t,m, r,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<double> >          (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<float> >   (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<double> >  (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int32_t> >         (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,n,t,m, r, iters,seed);
        // ok = ok &&run_with_field<Givaro::Modular<int64_t> >         (q,b,n,t,m, r, iters,seed); // Valgrind does not like this one 
        // ok = ok &&run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,n,t,m, r, iters,seed);
        // ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,9, ceil(n/4.), ceil(t / 4.), ceil(m / 4.), ceil(r / 4.), iters,seed);
        // ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:224), ceil(n/4.), ceil(t / 4.), ceil(m / 4.), ceil(r / 4.),iters,seed);
        seed++;
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;

    return !ok;
}
