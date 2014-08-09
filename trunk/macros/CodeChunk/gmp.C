#include <gmpxx.h>
			int main () {
			if (__GNU_MP_VERSION < 4) return -1;
			mpz_class a(2),b(3),c(5); if ( a+b == c ) return 0; else return -1; }


