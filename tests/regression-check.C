#include "fflas-ffpack/fflas-ffpack.h"

/*  #1  */
bool check1 () ;

/*  #2  */
bool check2()
{
	FFPACK::Modular<double> F(2);
	FFPACK::Modular<double>::RandIter R(F);

	size_t ok = 0 ;
	size_t tot = 500 ;
	for (size_t i = 0 ; i < tot ; ++i) {
		double elt ;
		R.random(elt);
		if (elt == 1) ++ok ;
	}
	double f = (double) ok / (double) tot ;
	if (f < 0.3 or f > 0.7) return false ;

	return true ;

}

/*  #3  */
bool check3()
{
	FFPACK::Modular<double> F(2);
	double * A = NULL ;
	double d = FFPACK::Det(F,0,0,A,0);
	return F.areEqual(d,F.one);

}

/*  #4  */
bool check4()
{
	FFPACK::Modular<double> F(2);
	double * A = NULL ;
	double * X = NULL ;
	int nul;
	FFPACK::Invert2(F,0,A,0,X,0,nul);
	return true ;
}


int main() {
	bool pass = true ;
	pass &= check2();
	pass &= check3();
	pass &= check4();
	return !pass;
}

