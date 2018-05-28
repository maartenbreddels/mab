#include "nlinear.hpp"
#include "polynomial.hpp"

template<class P>
void print_poly(P& p) {
cout << "order:" << P::degree << " " << p(0) << " " << p(0.25) << " " << p(0.5) << " " << p(0.75) << " " << p(1.) << endl;
	//typename P::derivative_type pd;// = p.derivative();
	//cout << "  def:" << P::derivative::degree << " " << pd(0) << " " << pd(0.25) << " " << pd(0.5) << " " << pd(0.75) << " " << pd(1.) << endl;
cout << "  der:" << P::degree << " " << p.dfdx(0) << " " << p.dfdx(0.25) << " " << p.dfdx(0.5) << " " << p.dfdx(0.75) << " " << p.dfdx(1.) << endl;
}

template<class P>
void print_bipoly(P& p) {
	for(double y = 0; y < 1.001; y += 0.25) {
		for(double x = 0; x < 1.001; x += 0.25) {
			cout << p(x, y) << " ";
		}
		cout << endl;
	}
}

extern "C" int main(int argc, char* const * argv) {
	Polynomial<1, 2, LagrangeBasis<2>> l2ta(1, 0, 0); 
	Polynomial<1, 2, LagrangeBasis<2>> l2tb(0, 1, 0); 
	Polynomial<1, 2, LagrangeBasis<2>> l2tc(0, 0, 1); 
	Polynomial<1, 2, MonomialBasis<2>> lpa(1, 0, 0); 
	Polynomial<1, 2, MonomialBasis<2>> lpb(0, 1, 0); 
	Polynomial<1, 2, MonomialBasis<2>> lpc(0, 0, 1); 
	Polynomial<1, 2, MonomialBasis<2>> lpd = (lpb + lpc);
	Polynomial<1, 1, MonomialBasis<1>> sa1(1, 0);
	Polynomial<1, 1, MonomialBasis<1>> sa2(0, 1);
	Polynomial<1, 1, LagrangeBasis<1>> sb1(1, 0);
	Polynomial<1, 1, LagrangeBasis<1>> sb2(0, 1);
	//double arr[] = {0, 0, 1};
	double a2 = 0;
	double a1 = 0;
	double a0 = 1;
	Polynomial<1, 2, MonomialBasis<2>, MonomialBasis<2>, double, double&> lpe(a2, a1, a0); 
/*LagrangePolynomial<2> l2a(1, 0, 0); 
	LagrangePolynomial<2> l2b(0, 1, 0); 
	LagrangePolynomial<2> l2c(0, 0, 1);*/
	cout << "Lagrange basis: 1 0 0" << endl;
	print_poly(l2ta);
	cout << "Lagrange basis: 0 1 0" << endl;
	print_poly(l2tb);
	cout << "Lagrange basis: 0 0 1" << endl;
	print_poly(l2tc);
	/*print_poly(l2a);
	print_poly(l2b);
	print_poly(l2c);*/
	cout << endl;
	print_poly(lpa);
	cout << "Monomial basis: 1 0 0" << endl;
	print_poly(lpa);
	cout << "Monomial basis: 0 1 0" << endl;
	print_poly(lpb);
	cout << "Monomial basis: 0 0 1" << endl;
	print_poly(lpc);
	cout << "Monomial basis: 0 1 1" << endl;
	print_poly(lpd);
	cout << endl;
	print_poly(lpe);
	cout << "simple" << endl;
	print_poly(sa1);
	print_poly(sa2);
	print_poly(sb1);
	print_poly(sb2);
	//arr[1] = 1;
	a0 = 2;
	a1 = 1;
	a2 = 0;
	cout << lpe.a << endl;
	cout << lpe.next.a << endl;
	cout << lpe.next.next.a << endl;
	//cout << lpe.a << endl;
	print_poly(lpe);
	cout << endl;
	
	NLinear<1, 2, double> linear(3., 2.);
	NLinear<1, 2, double> linear2(4., 6.);
	NLinear<2, 2, double> bilinear(linear, linear2);
	print_bipoly(bilinear);
	cout << endl;
	Polynomial<1, 1,LagrangeBasis<1>> pl1(2., 3.);
	Polynomial<1, 1,LagrangeBasis<1>> pl2(6., 4.);
	Polynomial<2, 1, LagrangeBasis<1>> pb(pl2, pl1);
	print_bipoly(pb);
	cout << endl;
	Polynomial<1, 2,LagrangeBasis<2>> pq1(2, 2.5, 3.);
	Polynomial<1, 2,LagrangeBasis<2>> pq2(4, 3.75, 3.5);
	Polynomial<1, 2,LagrangeBasis<2>> pq3(6, 5, 4);
	Polynomial<2, 2, LagrangeBasis<2>> pq(pq3, pq2, pq1);
	print_bipoly(pq);
	
}