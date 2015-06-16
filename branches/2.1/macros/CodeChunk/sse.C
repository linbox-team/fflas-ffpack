#include <immintrin.h>

int main() {
	// SSE 2
	__m128d P ;
	double p = 0;
	P   = _mm_set1_pd(p);
	P = _mm_add_pd(P,P);
	// SSE 4.1
	P = _mm_floor_pd(P);
	return 0;
}
