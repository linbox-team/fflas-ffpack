#include <interfaces/libs/fflas_c.h>
#include <interfaces/libs/ffpack_c.h>

#include <stdlib.h>
#include <stdio.h>

int main() {
	double * A = (double*)malloc(4*sizeof(double));
	A[0] = A[2] = 1 ;
	A[1] = A[3] = 0 ;
	size_t * P = (size_t*) malloc(2*sizeof(size_t));
	size_t * Qt = (size_t*) malloc(2*sizeof(size_t));
            
        size_t r = RowEchelonForm_modular_double(101,2,2,A,2,P,Qt,false,FfpackSlabRecursive,true);
        freducein_2_modular_double(101,2,2,A,2,false);
	freducein_1_modular_double(101,4,A,1,false);
	fsquare_3_modular_double(101,FflasNoTrans,2,1,A,2,1,A,1,true);
	return !(r==1);
}

