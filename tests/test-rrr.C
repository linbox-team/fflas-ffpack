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
 * \brief test equality between a dense RRR matrix and the result of compressing and reconstructing it.
 */
template<class Field>
bool test_compression  (const Field & F, size_t n, size_t t,
                        typename Field::Element_ptr A, size_t lda)
{
    typename Field::Element_ptr Acheck = fflas_new (F, n, n);
    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n, t, A, lda);
    
    RRRExpand<Field>(F, RRRA, Acheck, n);

    bool ok = fequal (F, n, n, A, lda, Acheck, n);
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
    
    RRgen<Field>* RRA = new RRgen(F, n_A, n_B, A, lda);
    RRgen<Field>* RRB = new RRgen(F, n_B, m_B, B, ldb);
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
    
    RRgen<Field>* RRA = new RRgen(F, n_A, m_A, A, lda);
    RRgen<Field>* RRB = new RRgen(F, n_A, m_A, B, ldb);
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

    RRRgen<Field>* RRRA = new RRRgen<Field>(F, n_A, t, A, lda);
    RRgen<Field>* RRB = new RRgen(F, n_A, m_A, B, ldb);

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

template<class Field>
bool launch_instance_check (const Field& F, size_t n, size_t t, size_t m, size_t r, typename Field::RandIter& G)
{
    typedef typename Field::Element_ptr Element_ptr;
    bool ok = true;
    
    // test compression with a random t qs-order matrix
    Element_ptr A3 = fflas_new (F, n, n);
    FFPACK::RandomLTQSMatrixWithRankandQSorder (F,n,r,t,A3,n,G);
    ok = ok && test_compression(F,n,2,A3,n);
    FFLAS::fflas_delete(A3);


    //test operations with random matrixes
    Element_ptr A = fflas_new (F, n, n);
    FFLAS::frand(F,G,n,n,A,n);
    Element_ptr B = fflas_new (F, n, n);
    FFLAS::frand(F,G,n,n,B,n);
    ok = ok && test_RRxRR(F,n,n,n,A,n,B,n);
    ok = ok && test_RRaddRR(F,n,n,A,n,B,n);
    ok = ok && test_RRRaddRR(F,n,n,t,A,n,B,n);
    FFLAS::fflas_delete(A);
    FFLAS::fflas_delete(B);


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
            

            FFLAS::fflas_delete(F);
            iters--;
        }
    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(20); // In order to print integers as integers even on float types, could be done once for all fflas
    Givaro::Integer q=-1;
    size_t b=0;
    size_t n=21;
    size_t t=7;
    size_t m=42;
    size_t r = 10;
    int iters=3;
    bool loop=true;
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
        std::cerr<<"Random matrixes tests"<<std::endl;
        ok = ok &&run_with_field<Givaro::Modular<float> >           (q,b,n,t,m, r,iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<double> >          (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<float> >   (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<double> >  (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int32_t> >         (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::ModularBalanced<int32_t> > (q,b,n,t,m, r, iters,seed);
        ok = ok &&run_with_field<Givaro::Modular<int64_t> >         (q,b,n,t,m, r, iters,
                                                                     seed); // Valgrind does not like this one 
        // ok = ok &&run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,n,t,m, r, iters,seed);
        // ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,9, ceil(n/4.), ceil(t / 4.), ceil(m / 4.), ceil(r / 4.), iters,
        //                                                             seed);
        // ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:224), ceil(n/4.), ceil(t / 4.), ceil(m / 4.), ceil(r / 4.),
        //                                                             iters,seed);
        seed++;
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;

    return !ok;
}