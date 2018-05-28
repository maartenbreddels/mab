#include <math.h>
#include "polynomial.hpp"

#define DO_SN
//#define DO_FOURIER

namespace gd {

class Profile;

namespace torus {

class trishape2d {
public:
	trishape2d(int i)  : p1(i == 3, i == 2), p2(i == 1, i == 0), p(p1, p2) {}
	Polynomial<1, 1,LagrangeBasis<1>> p1;
	Polynomial<1, 1,LagrangeBasis<1>> p2;
	Polynomial<2, 1, LagrangeBasis<1>> p;
	double operator()(double x, double y) {
	return p(x > 0.5 ? 2*(1-x) : x * 2, y > 0.5 ? (1-y)*2 : y*2);
	}
};


class TorusModelIsochrone {
public:
TorusModelIsochrone(double M, double scale, double G) : M(M), scale(scale), G(G) {
}
	
	double drdb(double J1, double J2, double theta1, double theta2);
	void get(double J1, double J2, double theta1, double theta2, double &H, double &Ekin, double &Epot, double &r);
	void get2(double J1_, double J2_, double theta1, double theta2, double &E, double &Ekin, double &Epot, double &r_, double &drdb, double &dphidb, double& dphidJ1, double& dphidJ2, double &dHdb, double &drdM, double &dphidM, double &dHdM, double &drdJ1, double &drdJ2, double &dHdJ1, double &dHdJ2);
	double potentialr(double r);
	double dphidr(double r);
	double dphidb(double r);
	double dHdb(double J1, double J2);
	double M, scale, G;
};

template<class P>
void print_bipoly(P& p) {
	for(double y = 0; y < 1.001; y += 0.1) {
		for(double x = 0; x < 1.001; x += 0.1) {
			//cout << p(x, y) << " ";
			printf("% 7.2f", p(x, y));
		}
		cout << endl;
	}
}


class TorusModelFit {
public:
	Profile* target_profile;
	TorusModelIsochrone* toy_torus;
	double J1, J2;
	int ns;
	double* S;
	trishape2d t1, t2, t3;
	trishape2d* shapes[3];

	TorusModelFit(Profile* target_profile, TorusModelIsochrone* toy_torus, double J1, double J2, int ns) : target_profile(target_profile), toy_torus(toy_torus), J1(J1), 	J2(J2), ns(ns), t1(1), t2(2), t3(3)
	{
		shapes[0] = &t1;
		shapes[1] = &t2;
		shapes[2] = &t3;
		cout << "1" << endl;
		print_bipoly(t1);
		cout << "2" << endl;
		print_bipoly(t2);
		cout << "3" << endl;
		print_bipoly(t3);
		S = new double[ns*ns-1];
		for(int i = 0; i < ns; i++) {
			for(int j = 0; j < ns; j++) {
				if((i > 0) || (j > 0)) {
					Sn(i, j) = 0;
				}
			}
		}
	}
	
	void fit();
	//double denergyChisqdb(double J1, double J2);
	double energyChisq();
	double optimize_function(double* gradient);
	int no_parameters() { return 2 + ns*ns-1; }
	
	//double Sn(int i, int j) { return S[i + ns*j - 1]; }
	double& Sn(int i, int j) { return S[Sn_index(i, j)]; }
	int Sn_index(int i, int j) { return i + ns*j - 1; }
	void getJtoy(double theta1_toy, double theta2_toy, double& J1_toy, double& J2_toy) {
		J1_toy = J1;
		J2_toy = J2;
#ifdef DO_SN
		for(int i = 0; i < ns; i++) {
			for(int j = 0; j < ns; j++) {
				if((i > 0) || (j > 0)) {
#ifdef DO_FOURIER
					double cos_ = cos(i * theta1_toy + j * theta2_toy);
					J1_toy += 2 * i * Sn(i,j) * cos_;
					J2_toy += 2 * j * Sn(i,j) * cos_;
#else
					double scale = 1. / (2*M_PI);
					int k = Sn_index(i, j);
					//printf("Sn(%d,%d) = %f // %f %f\n", i, j, Sn(i, j), J1 * (*shapes[k])(theta1_toy*scale, theta2_toy*scale) * Sn(i, j), J2 * (*shapes[k])(theta1_toy*scale, theta2_toy*scale) * Sn(i, j));
					J1_toy += J1 * (*shapes[k])(theta1_toy*scale, theta2_toy*scale) * Sn(i, j); //2 * i * Sn(i,j) * cos_;
					J2_toy += J2 * (*shapes[k])(theta1_toy*scale, theta2_toy*scale)  * Sn(i, j); //2 * j * Sn(i,j) * cos_;
#endif
				}
			}
		}
#endif
	}
};

} // namespace torus
} // namespace gd