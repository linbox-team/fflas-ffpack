#include <immintrin.h>
int main() {
	__m256d P ;
	double p = 0;
	P   = _mm256_set1_pd(p);
	P = _mm256_add_pd(P,P);
#ifdef __try_avx2
	P = _mm256_fnmadd_pd(P,P,P);
#endif
	return 0;
}
