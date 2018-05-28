#include "integration.hpp"
#include<iostream>
#include <time.h>

using namespace gd;
using namespace std;

double ftest(double x) {
	return x*x;
}

double Ftest(double a, double b) {
	return b*b*b/3 - a*a*a/3;
}



int main() {
	IntegratorGrid<> igrid(ftest, 0, 10, 0.25);
	IntegratorGrid<double (*)(double)> igrid2(ftest, 0, 10, 0.25);
	IntegratorGSL<> igsl(ftest);
	cout << ftest(0.05) << " " << ftest(0.1) << " " << ftest(0.15) << endl;
	cout << "F(a,b) \\approx      "  << igrid.integrate(0.05, 4.15) << endl;
	cout << "F(a,b) \\approx2     "  << igrid2.integrate(0.05, 4.15) << endl;
	cout << "F(a,b) \\approx(gsl) "  << igsl.integrate(0.05, 4.15) << endl;
	cout << "F(a,b) = " << Ftest(0.05,  4.15) << endl;
	
	clock_t t1, t2;
	int i;
	
	t1 = clock();
	i = 0;
	while(i < 1000000) {
		double x = igrid2.integrate(0.05, 4.15);
		i++;
	}
	t2 = clock();
	cout << "diff = " << (t2-t1) << endl;
	
	t1 = clock();
	i = 0;
	while(i < 1000000) {
		double x = igrid.integrate(0.05, 4.15);
		i++;
	}
	t2 = clock();
	cout << "diff = " << (t2-t1) << endl;
	
	i = 0;
	t1 = clock();
	i = 0;
	while(i < 1000000) {
		double x = igsl.integrate(0.05, 4.15);
		i++;
	}
	t2 = clock();
	cout << "diff = " << (t2-t1) << endl;
	return 0;
}
