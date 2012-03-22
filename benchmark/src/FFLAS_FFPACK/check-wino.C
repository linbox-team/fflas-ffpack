	//#define LinBoxSrcOnly
#include <iostream>
#include <fstream>
    //#define _LINBOX_LINBOX_CONFIG_H
#define __FFLASFFPACK_CONFIGURATION
#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas-ffpack-config.h"
#include "fflas-ffpack/field/modular-positive.h"
#include "fflas-ffpack/fflas/fflas.h"
#include <utils/timer.h>

#define CUBE(x) ((x)*(x)*(x))

    //using namespace LinBox;
int main (int argc, char ** argv) {
        const size_t p = argc>1 ? atoi(argv[1]) : 65521;
        const size_t n = argc>2 ? atoi(argv[2]) : 1000;
        const size_t NB = argc>3 ? atoi(argv[3]) : 1;
        const size_t winomax = argc>4 ? atoi(argv[4]) : 8;
        const size_t seed = argc>5 ? atoi(argv[5]) : BaseTimer::seed() ;
        srand((unsigned int)seed);

        {

        FFPACK::Modular<double> F( p );
        double basetime(0.0), time(0.0);

        double *A, *C;
        A = new double[n*n];
        C = new double[n*n];
        for (size_t i=0; i<n*n;++i){
		    A[i]= (double)((size_t)((double)i*rand())%p);
        }


        Timer chrono;
        for(size_t i=0; i<NB; ++i) {
            chrono.start();
            FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                         n, n, n, 1., A, n, A, n, 0., C, n, 0);
            chrono.stop();
            basetime+= chrono.usertime();
        }
        std::cout << std::endl
                  << "fgemm " << n << "x" << n << ": "
                  << basetime/(double)NB << " s, "
                  << (2.0/basetime*(double)NB*CUBE((double)n/100.0)) << " Mffops"
                  << std::endl;

        for(size_t w=1; w<winomax; ++w) {

            chrono.clear();
            for(size_t i=0; i<NB; ++i) {
                chrono.start();
                FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                             n, n, n, 1., A, n, A, n, 0., C, n, w);
                chrono.stop();
                time+= chrono.usertime();
            }
            std::cout << w << "Wino " << n << "x" << n << ": "
                      << time/(double)NB << " s, "
                      << (2.0/time*(double)NB*CUBE((double)n/100.0)) << " Mffops"
                      << std::endl;
        }


        delete[] A;
        delete[] C;
        }
        {

        FFPACK::Modular<float> F( p );
        double basetime(0.0), time(.0);

        float *A, *C;
        A = new float[n*n];
        C = new float[n*n];
        for (size_t i=0; i<n*n;++i){
		    A[i]= (float)((size_t)((float)i*(float)rand())%p);
        }


        Timer chrono;
        for(size_t i=0; i<NB; ++i) {
            chrono.start();
            FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                         n, n, n, 1., A, n, A, n, 0., C, n, 0);
            chrono.stop();
            basetime+= chrono.usertime();
        }
        std::cout << std::endl
                  << "fgemm " << n << "x" << n << ": "
                  << basetime/(double)NB << " s, "
                  << (2.0/basetime*(double)NB*CUBE((double)n/100.0)) << " Mffops"
                  << std::endl;

        for(size_t w=1; w<winomax; ++w) {

            chrono.clear();
            for(size_t i=0; i<NB; ++i) {
                chrono.start();
                FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
                             n, n, n, 1., A, n, A, n, 0., C, n, w);
                chrono.stop();
                time+= chrono.usertime();
            }
            std::cout << w << "Wino " << n << "x" << n << ": "
                      << time/(double)NB << " s, "
                      << (2.0/time*(double)NB*CUBE((double)n/100.0)) << " Mffops"
                      << std::endl;
        }


        delete[] A;
        delete[] C;
        }


        return 0;
    }

