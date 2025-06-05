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
 * \brief test equality between a dense sss matrix and the result of compressing and reconstructing it.
 */
template<class Field>
bool test_compression  (const Field & F, size_t n, size_t t,
                        typename Field::Element_ptr A, size_t lda)
{
    typename Field::Element_ptr Acheck = fflas_new (F, n, n);
    typename Field::Element_ptr Ainit = fflas_new (F, n, n);

    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < n; j++){
            Ainit[i*n+j] = A[i*n+j];
        }
    } 

    RRRrep<Field>* RRRA = PLUQRRRGen<Field>(F, n, t, A, lda);
    RRRExpand<Field>(F, *RRRA, Acheck, n);

    bool ok = fequal (F, n, n, Ainit, lda, Acheck, n);
    if ( !ok )
        {
            std::cout << "ERROR: different results for dense to RRR and RRR to dense (RRRGen and Expand)"
                      <<std::endl;
            WriteMatrix(std::cout << "Ainit = " << std::endl, F, n, n, Ainit, lda);
            // WriteMatrix(std::cout << "Amodif = " << std::endl, F, n, n, A, lda);
            Node<Field>* root = RRRA->getroot();
            WriteMatrix(std::cout << "Acheck =  " << std::endl, F, n, n, Acheck, n);
            WriteMatrix(std::cout << "U_u = " << std::endl, F, root->ru, root->size_N2, root->U_u, root->size_N2);
            WriteMatrix(std::cout << "L_u = " << std::endl, F, root->size_N1, root->ru, root->L_u, root->ru);
            WriteMatrix(std::cout << "U_l = " << std::endl, F, root->rl, root->size_N1, root->U_l, root->size_N1);
            WriteMatrix(std::cout << "L_l = " << std::endl, F, root->size_N2, root->rl, root->L_l, root->rl);
        }

    FFLAS::fflas_delete(Acheck);

    return ok;
}

template<class Field>
bool launch_instance_check (const Field& F, size_t n, size_t t, size_t m, size_t r, typename Field::RandIter& G)
{
    typedef typename Field::Element_ptr Element_ptr;
    bool ok = true;


    // test with the matrix identity
    // Element_ptr A1 = fflas_new(F,n,n);
    // fidentity(F, n, n, A1, n);
    // ok = ok && test_compression(F,n,0,A1,n);
    // FFLAS::fflas_delete(A1);


    // test with a constructed matrix

    // Element_ptr A2 = fflas_new(F,n,n);
    // fidentity(F, n, n, A2, n);          // A2 = In
    // F.assign(A2[4], 143);                  // A2[0][4]= 1
    // F.assign(A2[3+2*n], 12);           // A2[1][5]= 143
    // F.assign(A2[5+n], 41);                  // A2[0][4]= 1
    // F.assign(A2[4+n], 3);
    // F.assign(A2[4+2*n], 35);
    // F.assign(A2[2], 5);
    // F.assign(A2[1+4*n], 10);           // A2[1][5]= 143
    // F.assign(A2[2+3*n], 13);           // A2[1][5]= 143
    // ok = ok && test_compression(F,n,2,A2,n);
    // FFLAS::fflas_delete(A2);



    // test with a random matrix
    Element_ptr A3 = fflas_new (F, n, n);
    FFPACK::RandomMatrix(F,n,n,A3,n);
    ok = ok && test_compression(F,n,2,A3,n);
    FFLAS::fflas_delete(A3);

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

            delete F;
            iters--;
        }
    return ok;
}

int main(int argc, char** argv)
{
    cerr<<setprecision(20); // In order to print integers as integers even on float types, could be done once for all fflas
    Givaro::Integer q=-1;
    size_t b=0;
    size_t n=7;
    size_t t=2;
    size_t m=42;
    size_t r = 40;
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
        ok = ok &&run_with_field<Givaro::ModularBalanced<int64_t> > (q,b,n,t,m, r, iters,seed);
        // ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,9, ceil(n/4.), ceil(t / 4.), ceil(m / 4.), ceil(r / 4.), iters,
        //                                                             seed);
        // ok = ok &&run_with_field<Givaro::Modular<Givaro::Integer> > (q,(b?b:224), ceil(n/4.), ceil(t / 4.), ceil(m / 4.), ceil(r / 4.),
        //                                                             iters,seed);
        seed++;
    } while (loop && ok);

    if (!ok) std::cerr<<"with seed = "<<seed-1<<std::endl;

    return !ok;
}