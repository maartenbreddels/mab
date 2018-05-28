#include "autodiff_sharedptr2.hpp"
#include "torus.hpp"
#include <iostream>
#include <sys/times.h>
#include <stdio.h>
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>

//using namespace fadbad;
using namespace std;

/*F<double> func(const F<double>& x, const F<double>& y)
{
	F<double> z=sqrt(x);
	return y*z+sin(z);
}*/

#define N 1000000000L
volatile double total;
int la = 1;
int lb = 2;
int lc = 3;
int ld = 4;
const double G = 4.30200406461e-06;
double J1 = 1.0;
double J2 = 0.1;
double theta1 = 0.1;
double theta2 = 0.2;
double scale = 1.;
//double b = self.p.scale

double psi = 0.110344724358;

namespace ad {
template<class Ta, class Tb, class Te>
struct PsiTorus {
	typedef typename Ta::type type;
	typedef typename Ta::vartype vartype;
	//typedef Div<Mul<Ta>, Sum<>> f2;
	//typedef Mul<f2, C<double>> f;
	//typedef Sub<C<double>, f> dtype;
	//typedef Neg<Mul<Sin<T>, typename T::dtype>> dtype;
	PsiTorus(Ta a, Tb b, Te e, double theta1) : a(a), b(b), e(e), theta1(theta1) { }
	Ta a;
	Tb b;
	Te e;
	double theta1;

	type eval() const {
		return psi;
	}
	type d(int i) const {
		auto u = ((a*e)/(a+b));
		double uval = u.eval();
		double dpsi_du = (psi*::cos(psi)-::sin(psi))/(uval*::cos(psi)-1) - (::cos(psi)*(uval*psi*::cos(psi)-theta1-uval*::sin(psi)))/::pow((uval*::cos(psi)-1), 2);
		cout << "dpsidu = " << dpsi_du << " dudX: " << u.d(i) << endl;
		return ((psi*::cos(psi)-::sin(psi))/(uval*::cos(psi)-1) - (::cos(psi)*(uval*psi*::cos(psi)-theta1-uval*::sin(psi)))/::pow((uval*::cos(psi)-1), 2)) * u.d(i);
		//return dpsi_du * u.d(i);

	}
	/*dtype dx(int i) const {
		auto u = ((a*e)/(a+b));
		//return 1. - ((a*e)/(a+b))*::cos(psi);
	}*/
	/*type d(vartype &x) const {
		return 0;
	}*/
	/*template<int i>
	type d() const {
		return 0;
	}*/

	/*Sin<T> dx(int i) {
		return -Sin<T>(t.eval()) * t.dx<i>();
	}*/
};
template<class Ta, class Tb, class Te>
PsiTorus<Ta, Tb, Te> psitorus(Ta a, Tb b, Te e, double theta1) {
	return PsiTorus<Ta, Tb, Te>(a, b, e, theta1);
}
}

template<class Rho=double, class Rs=double>
class NFW  {
public:
	Rho rho0;
	Rs rs;
	double G;
	NFW(Rho rho0, Rs rs, double G) : rho0(rho0), rs(rs), G(G) {
	}
	/*template<class R>
	auto potentialr(R r) -> decltype(- 4 * M_PI * G * rho0 * rs*rs * log(1+r/rs)/(r/rs)) {
		auto x = r/rs;
		return - 4 * M_PI * G * rho0 * rs*rs * log(1+x)/x;
}*/
	double densityr(double r) {
		auto x = r/rs;
		return rho0 / (x * pow(1+x, 2));
	}
	double densityR(double R) {
		return 0;
	}
	double I(double r, double I0=1.0) {
		return 0;
	}
	double dphidr(double r) {
		double x = r/rs;
		double f = -4.0*M_PI*G*rho0*pow(rs, 3)/pow(r,2);
		double t1 = log(1.0 + x);
		double t2 = x/(1.0+x);
		//cout << x << ", " << f << ", " << t1 << ", " << t2 << endl;
		return -f*(t1 - t2);
	}
};

namespace ad {
	template<class T, class V>
	double finite_difference(T expression, V& value) {
		//cout << "value before:" << value.eval() << endl;
		double epsilon = 0.000001;
		double oldvalue = get(value);
		double y1 = get(expression);
		//cout << "y before:" << y1 << endl;
		set(value, oldvalue + epsilon);
		//expression->reset();
		expression->reevaluate(); 
		double y2 = get(expression);
		//cout << "y after:" << y2 << endl;
		return (y2-y1)/epsilon;
	}
	/*template<class T1, class T2>
	auto somefunc(T1 x, T2 y) -> decltype(pow(x, y) * x) {
		return pow(x, y) * x;
}*/
	void testtorus() {
		/*double L = J2;
		V<double,0> b(scale);
		V<double,1> M(1e8);
		auto k = M * G;
		auto H = -2. * (k*k) / sqr((2.*J1 + L + sqrt(4.*b*k+(L*L)) )); 
		//auto H = -2. * (k*k) / sqr(2*J1 + L + b);
		cout << "H = " << H.eval() << endl;
		auto a = -k / (2 * H) - b;
		auto e = sqrt(1+L*L/(2*H*a*a));
		cout << "a = " << a.eval() << endl;
		cout << "e = " << e.eval() << endl;
		auto psi = psitorus(a, b, e, theta1);
		cout << psi.eval() << " " << psi.d(0) << endl;
		//auto r = a * sqrt((1.-e*cos(psi))*(1.-e*cos(psi)+2.*b/a))
		auto r = a * sqrt((1.-e*cos(psi))*(1.-e*cos(psi)+2.*b/a));
		//auto r = a * sqrt(e*cos(psi)+2.*b/a);
		cout << "drdb = " << r.d(0) << endl;
		cout << "drdM = " << r.d(1) << endl;

		V<double,3> rho0(10);
		V<double,4> rs(1);
		NFW<V<double,3>, V<double,4>> nfw(rho0, rs, G);
		auto Pr = nfw.potentialr(r);
		cout << "Pot(r) " << Pr.eval() << endl;
		
		*/
		/*#print "energy", H, "original", E
		# step 2
		a = -self.k / (2 * H) - self.b
		# step 3
		e = sqrt(1+L**2/(2*H*a**2))
		#wl = sqrt(self.k) / (2*(a+self.b)**(3./2.)) * (1+L/sqrt(4*self.b*self.k+L**2))
		#print "a,e,wl", a, e, wl, a*e/(a+b)
		# step4
		def f_psi(psi):
			return psi - a*e/(a+self.b)*sin(psi) - theta1
		def f_psi_prime(psi):
			return 1 - a*e/(a+self.b)*cos(psi)
		res = fsolve(f_psi, pi/2, fprime=f_psi_prime, full_output=0)
		psi = res
			r = a * sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a))
		*/

	}
	void test() {
		auto x = var(3., varindex<0>()); 
		auto y = var(4., varindex<1>()); 
		//shared_ptr<V<double,0>> a(3);
		//V<double,1> b(6);
		//V<double,2> c(5);
		//auto r = 4*a + b + a * b * c;
		//auto r = exp(a*a*b + a*a);
		//auto r = 4*a;
		//auto r = somefunc(a, b);
		//auto rt = x-y; //log(x); //sin(pow(a, 2)) * a * a; //exp(4*a*a) + a*b;
		//auto r1 = (x * x); //log(x); //sin(pow(a, 2)) * a * a; //exp(4*a*a) + a*b;
		auto r = pow(x, y);//log(x); //sin(pow(a, 2)) * a * a; //exp(4*a*a) + a*b;
		auto drdx = r->dx(x);
		
		
		
		//auto d = r->dx(x); //->evaluate();
		//auto drdy = r->dx(y);
		
		//auto r = deriv(x, x);
		//auto r2 = partial(y, x);
		//drdx->test();
		//drdx->test();

		//cout << "value: " << r->value() << endl;
		r->evaluate();
		drdx->evaluate();
		//cout << "test " << drdx->evaluate()<< endl;
		//cout << drdy->evaluate()<< endl;
		//drdx->test();
		cout << "value: " << r->value() << endl;
		cout << "x: " << r->d(x) << " " << finite_difference(r, x);
		//cout << " " << drdx->value();
		cout << endl;
		cout << "y: " << r->d(y) << " " << finite_difference(r, y);
		cout << " " << drdx->value();
		cout << endl;
		//cout << "y: " << r->d(y) << " " << finite_difference(r, y) << " " << drdy->value() << endl; //
		//" " << r.dx(0).eval() << " " << r.dx(0).dx(0).eval() << endl;

		//auto r = log(a); //exp(4*a*a) + a*b;
		//cout << r.eval() << endl;
		//cout << a.eval() << " " << r.t1.eval() << endl;

		//r.crash();
		/*tms begin, end;
 		times(&begin);
		for(long long int i = 0; i < N; i++) {
			total = r.d(la);
			total = r.d(lb);
			total = r.d(lc);
			total = r.d(ld);
		}
 		times(&end);
		printf("%d %d\n", begin.tms_utime, end.tms_utime);
		printf("user time=%f\n", ((float)(end.tms_utime-begin.tms_utime)));
 		times(&begin);
		for(long long int i = 0; i < N; i++) {
			total = r.d<0>();
			total = r.d<1>();
			total = r.d<2>();
			total = r.d<3>();
		}
 		times(&end);
		printf("%d %d\n", begin.tms_utime, end.tms_utime);
		printf("user time=%f\n", ((float)(end.tms_utime-begin.tms_utime)));
		*/
		//cout << "0: " << r.t1.d(0) << endl;
		//cout << "0: " << r.d(0) << endl;
		/*cout << "0: " << r.d(0) << " " << r.d<0>() << " " << finite_difference(r, a) << " " << r.dx(0).eval() << " " << r.dx(0).dx(0).eval() << endl;
		//cout << r.dx(0).eval() << endl;
		cout << "1: "<< r.d(1) << " " << r.d<1>() << endl;
		cout << "2: "<< r.d(2) << " " << r.d<2>() << endl;
		//cout << r.d(a) << endl;
		//cout << r.d(b) << endl;
		//cout << a.d(a) << endl;
		//cout << a.d(b) << endl;
		//cout << r.d(a) << endl;
		//cout << r.eval();
		//		cout << r.d();
		*/
	}
}
int main2() {
	/*F<double> x,y,f;     // Declare variables x,y,f
	x=1;                 // Initialize variable x
	x.diff(0,2);         // Differentiate with respect to x (index 0 of 2)
	y=2;                 // Initialize variable y
	y.diff(1,2);         // Differentiate with respect to y (index 1 of 2)
	f=func(x,y);         // Evaluate function and derivatives
	double fval=f.x();   // Value of function
	double dfdx=f.d(0);  // Value of df/dx (index 0 of 2)
	double dfdy=f.d(1);  // Value of df/dy (index 1 of 2)
	
	cout << "f(x,y)=" << fval << endl;
	cout << "df/dx(x,y)=" << dfdx << endl;
	cout << "df/dy(x,y)=" << dfdy << endl;*/
	return 1;
}

int main() {
	//ad::testtorus();
	//ad::test();
	gd::torus::TorusModelIsochrone t(1e8, 1.0, G);
	t.drdb(J1, J2, theta1, theta2);
}