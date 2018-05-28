#include <assert.h>

#include "torus.hpp"
#include "autodiff_sharedptr2.hpp"
#include "rootfinder.hpp"
#include "profile.hpp"
#include "optimization.hpp"


namespace gd { namespace torus {

using namespace ad;

/*
def f_psi(psi):
	return psi - a*e/(a+self.b)*sin(psi) - theta1
def f_psi_prime(psi):
	return 1 - a*e/(a+self.b)*cos(psi)
*/

template<class T1, class T2, class T3>
class ManualDiff {
public:
 	typedef typename T2::type type;
	//typedef BinaryExpr<U, Theta1, Sub<U, Theta1>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename result_of_F_dG<T2, T3, X>::type partial_type;
	};
	T1 value_;
	shared_ptr<T2> dfdu;
	shared_ptr<T3> u;

	ManualDiff(T1 value, shared_ptr<T2> dfdu, shared_ptr<T3> u) : value_(value), dfdu(dfdu), u(u) { }
	double _evaluate() const {
		return value_;
	}
	double value() { return value_; }
	void evaluate() const {
		//v = diff->evaluate();
	}
	template<class X>
	double d(X &x) const {
		return dfdu->value() * u->d(x);
	}

	/*template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type> dx(shared_ptr<X> &x)  {
		return shared_ptr<T2>(dfdu) * u->dx(x);
	}*/
};

template<class T1, class T2, class T3>
shared_ptr<ManualDiff<T1, T2, T3>> manual_diff(T1 value, shared_ptr<T2> dfdu, shared_ptr<T3> u) {
	return shared_ptr<ManualDiff<T1, T2, T3>>(new ManualDiff<T1, T2, T3>(value, dfdu, u)); 
}

double solve_psi(double u, double theta1, double psi_start) {
	// numerically solve psi
	auto fdf = [&](double psi, double* y, double *dy) {
		*y = psi - u * sin(psi) - theta1;
		*dy = 1 - u * cos(psi);
	};
	gd::RootFinderGSLDerivative<> rf(fdf);
	return rf.findRoot(M_PI/2.);
}

void TorusModelIsochrone::get2(double J1_, double J2_, double theta1, double theta2, double &E, double &Ekin, double &Epot, double &r_, double &drdb, double &dphidb, double& dphidJ1, double& dphidJ2, double &dHdb, double &drdM, double &dphidM, double &dHdM, double &drdJ1, double &drdJ2, double &dHdJ1, double &dHdJ2) {
	auto J1 = var(J1_, varindex<3>());
	auto J2 = var(J2_, varindex<4>());
	//auto L = J2;
	auto b = var(scale, varindex<1>());
	auto M = var(this->M, varindex<2>());
	auto k = M * G;
	auto H = -2. * (k*k) / sqr((2.*J1 + J2 + sqrt(4.*b*k+(J2*J2)) ));

	E = H->evaluate();
	//cout << "H = " << H->value() << endl;
	auto a = -k / (2. * H) - b;
	auto e = sqrt(1.+J2*J2/(2.*H*a*a));
	auto u = ((a*e)/(a+b));
	double uvalue = u->evaluate();
	

	double psi_value = solve_psi(uvalue, theta1, M_PI/2);
	//cout << "1] psi = " << psi_value << " " << a->value()  << " " << e->value()  << " " << u->value()  << " " << k->value() << " " << E << endl;

 
	auto dpsidu = (psi_value*::cos(psi_value)-::sin(psi_value))/(u*::cos(psi_value)-1.) - (::cos(psi_value)*(u*psi_value*::cos(psi_value)-theta1-u*::sin(psi_value)))/pow((u*::cos(psi_value)-1.), 2.);
	dpsidu->evaluate();
	auto psi = manual_diff(psi_value, dpsidu, u);

	auto r = a * sqrt((1.-e*cos(psi))*(1.-e*cos(psi)+2.*b/a));
	//auto r = ((1-e*cos(psi)));

	r_ = r->evaluate();
	drdb = r->d(b);
	drdM = r->d(M);
	
	auto potential = - G * M / (b + sqrt(r*r + b*b));
	Epot = potential->evaluate();
	dphidb = potential->d(b);
	dphidM = potential->d(M);
	dphidJ1 = potential->d(J1);
	dphidJ2 = potential->d(J2);
	Ekin = E - Epot;
	dHdb = H->d(b);
	dHdM = H->d(M);
	dHdJ1 = H->d(J1);
	dHdJ2 = H->d(J2);

	drdJ1 = r->d(J1);
	drdJ2 = r->d(J2);

	//double drdM = r->d(M);
	
	
}


double TorusModelIsochrone::drdb(double J1, double J2, double theta1, double theta2) {
	double L = J2;
	auto b = var(scale, varindex<1>());
	auto M = var(this->M, varindex<2>());
	auto k = M * G;
	auto H = -2. * (k*k) / sqr((2.*J1 + L + sqrt(4.*b*k+(L*L)) ));

	H->evaluate();
	//cout << "H = " << H->value() << endl;
	auto a = -k / (2 * H) - b;
	auto e = sqrt(1.+L*L/(2.*H*a*a));
	auto u = ((a*e)/(a+b));
	double uvalue = u->evaluate();
	

	// numerically solve psi
	/*auto fdf = [&](double psi, double* y, double *dy) {
		*y = psi - uvalue * sin(psi) - theta1;
		*dy = 1 - uvalue * cos(psi);
	};
	gd::RootFinderGSLDerivative<> rf(fdf);
	double psi_value = rf.findRoot(M_PI/2.);
	*/
	double psi_value = solve_psi(uvalue, theta1, M_PI/2);
	//cout << "psi = " << psi_value << endl;

 
	auto dpsidu = (psi_value*::cos(psi_value)-::sin(psi_value))/(u*::cos(psi_value)-1.) - (::cos(psi_value)*(u*psi_value*::cos(psi_value)-theta1-u*::sin(psi_value)))/pow((u*::cos(psi_value)-1.), 2.);
	dpsidu->evaluate();
	auto psi = manual_diff(psi_value, dpsidu, u);

	//cout << "dpsidb = " << psi->d(b) << endl;
	//cout << "dpsidM = " << psi->d(M) << endl;


	
	auto r = a * sqrt((1.-e*cos(psi))*(1.-e*cos(psi)+2.*b/a));
	//auto r = ((1-e*cos(psi)));

	r->evaluate();
	double drdb = r->d(b);
	//double drdM = r->d(M);
	//cout << "drdb = " << drdb << endl;
	//cout << "drdM = " << drdM << endl;
	return drdb;
}
double sqr(double x) {
	return x * x;
}
void TorusModelIsochrone::get(double J1, double J2, double theta1, double theta2, double &H, double &Ekin, double &Epot, double &r) {
	double L = J2;
	double b = scale;
	//auto M = this->V;
	double k = M * G;
	H = -2. * (k*k) / sqr(2.*J1 + L + sqrt(4.*b*k+L*L));

	//cout << "H = " << H << endl;
	auto a = -k / (2 * H) - b;
	auto e = sqrt(1.+L*L/(2.*H*a*a));
	auto u = ((a*e)/(a+b));

	//auto wl = sqrt(k) / (2*pow(a+b, 3./2)) * (1+L/sqrt(4*b*k+L*L));
	auto psi = solve_psi(u, theta1, M_PI/2);
	//cout << "2] psi = " << psi << " " << a << " " << e << " " << u  << " " << k  << " " << H << endl;
	r = a * sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*b/a));
	Epot = potentialr(r);
	Ekin = H - Epot;
	/*#print "check", psi, psi - a*e/(a+b)*sin(psi) - theta1
	wratio = 0.5 * (1+J2/sqrt(4*self.k*self.b+J2**2))
	Lambda_psi = arctan(sqrt((1+e)/(1-e))*tan(psi/2)) +\
			L/sqrt(L**2+4*self.b*self.k)*\
			arctan( sqrt((a*(1+e)+2*self.b)/(a*(1-e)+2*self.b))*tan(0.5*psi) )
	# step 6: theta
	phi = theta2 - wratio*theta1 + Lambda_psi
	# step 7: pr,ptheta
	pr = sqrt(self.k/(a+self.b)) * a*e*sin(psi)/r
	vr = pr
	pphi = L#/cos(phi)
	vphi = pphi/r
	#print "psi, phi", psi, phi
	#print "p", pr, pphi
	#print "r,rorg", r, rorg, 1/(r/rorg)
	#print "r,phi", r, phi
	#print "vr,vphi", vr, vphi
	Epot = self.galaxy.potentialr(r)
	Ekin = 0.5 * (vphi**2 + vr**2)*/
	//return r, phi, vr, vphi, Ekin, Epot, Ekin+Epot
}

double TorusModelIsochrone::dHdb(double J1, double J2) {
	double k = M * G;
	double b = scale;
	double L = J2;
	return 8 * pow(k, 3) / (sqrt(4*b*k+L*L) * pow(2 * J1 + L + sqrt(4*b*k+L*L),3));
}

double TorusModelIsochrone::potentialr(double r) {
	return - G * M / (scale + sqrt(r*r + scale*scale));
}

double TorusModelIsochrone::dphidr(double r) {
	double a = sqrt(scale*scale + r * r);
	return G * M  * r / (a * sqr(scale+a));
}

double TorusModelIsochrone::dphidb(double r) {
	auto b = var(scale, varindex<1>());
	auto M = var(this->M, varindex<2>());
	//auto k = M * G;
	auto pot = - G * M / (b + sqrt(r*r + b*b));
	pot->evaluate();
	return pot->d(b);
}


void TorusModelFit::fit() {

}

double TorusModelFit::optimize_function(double* gradient) {
	const int N = 4;
	double energies[N*N];
	double energy_mean = 0;
	double dHdbs[N*N];
	double dHdMs[N*N];
	double dHdSns[N*N][4];
	double dHdbs_mean = 0;
	double dHdMs_mean = 0;
	double dHdSns_mean[4] = {0};

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			double theta1_toy = 2 * M_PI * i / N;
			double theta2_toy = 2 * M_PI * j / N;
			double H, Ekin, Epot_toy, r, drdb, dphidb, dphidJ1, dphidJ2, dHdb, drdM, dphidM, dHdM, drdJ1, drdJ2, J1_toy, J2_toy, dHdJ1, dHdJ2;
			//toy_torus->get(J1, J2, theta1, theta2, H, Ekin, Epot, r);
			getJtoy(theta1_toy, theta2_toy, J1_toy, J2_toy);

			//cout << "] " << theta1_toy << " " << theta2_toy << " " << J1_toy << " " << J2_toy << " " << J1 << " " << J2 << endl;
			assert(J1_toy >= 0);
			assert(J2_toy >= 0);
			//J1_toy = max(J1_toy, 0.00001);
			//J2_toy = max(J2_toy, 0.00001);
			toy_torus->get2(J1_toy, J2_toy, theta1_toy, theta2_toy, H, Ekin, Epot_toy, r, drdb, dphidb, dphidJ1, dphidJ2, dHdb, drdM, dphidM, dHdM, drdJ1, drdJ2, dHdJ1, dHdJ2);
			//cout << "r = " << r << " H = " << H << " " << toy_torus->potentialr(r) << endl;//
			/*double H1 = H;
			double dscale = 0.00001;
			toy_torus->scale += dscale;
			toy_torus->get(J1, J2, theta1, theta2, H, Ekin, Epot, r);
			double H2 = H;
			toy_torus->scale -= dscale;
			//cout << "dHdb = " << ((H2-H1)/dscale) << " " << toy_torus->dHdb(J1, J2) << endl;
			*/
			//cout << "1] " << H << " " << Ekin << " " << Epot_toy << " " << r << endl;
			//toy_torus->get(J1, J2, theta1_, theta2, H, Ekin, Epot_toy, r);
			//cout << "2] " << H << " " << Ekin << " " << Epot_toy << " " << r << endl;
			double Epot = target_profile->potentialr(r);
			energies[i+j*N] = Ekin+Epot;
			energy_mean += energies[i+j*N] / (N*N); 

			//double drdb = toy_torus->drdb(J1, J2, theta1, theta2);
			//double dEkindr = toy_torus->dHdb(J1, J2)/drdb - toy_torus->dphidr(r);
			//double dHdr = target_profile->dphidr(r) + dEkindr;
			//cout << "> " << toy_torus->dphidr(r) << " " << drdb << " " << toy_torus->potentialr(r) << endl;
			//cout << "> " << toy_torus->dHdb(J1, J2) << " " << drdb << endl;
			//dHdb[i+j*N] = toy_torus->dHdb(J1, J2) - toy_torus->dphidr(r)*drdb ; //dHdr;
			//dHdb[i+j*N] = toy_torus->dHdb(J1, J2) - toy_torus->dphidr(r)*drdb ; //dHdr;
			//double dr = 0.0001;
			//cout << "dphidr = " << ((toy_torus->potentialr(r+dr) - toy_torus->potentialr(r))/dr) << " " << toy_torus->dphidr(r) << endl;
			//return ; //* 0 + toy_torus->dphidb(r);

			//double dJ1dSn = 0;
			//double dJ2dSn = 0;
			//double dHdr = dHdb / drdb + dHdM / drdM + dHdJ1 / drdJ1 + dHdJ2 / drdJ2;
			//cout << ") " << dHdb / drdb << " " <<  dHdM / drdM  << " " << dHdJ1 / drdJ1  << " " << dHdJ2 / drdJ2 << endl;
			//cout << "] " << dHdr << " " << toy_torus->dphidr(r) << " " << drdJ1 << endl;
			dHdbs[i+j*N] = target_profile->dphidr(r)*drdb + dHdb-dphidb;
			dHdMs[i+j*N] = target_profile->dphidr(r)*drdM + dHdM-dphidM;

#ifdef DO_SN
			for(int k = 0; k < ns; k++) {
				for(int l = 0; l < ns; l++) {
					if((k > 0) || (l > 0)) {
#ifdef DO_FOURIER
						double cos_ = cos(k * theta1_toy + l * theta2_toy);
						//double cos_ = cos(1 * theta1_toy + 0 * theta2_toy);
						double dJ1_toydSn = 2 * k * cos_;
						double dJ2_toydSn = 2 * l * cos_;
#else
						int m = Sn_index(k, l);
						double scale = 1. / (2*M_PI);
						double dJ1_toydSn = J1 * (*shapes[m])(theta1_toy*scale, theta2_toy*scale);
						double dJ2_toydSn = J2 * (*shapes[m])(theta1_toy*scale, theta2_toy*scale);
#endif
						//if(m == 0)
						//	return dJ1_toydSn;
						//cout << "] " << cos_ << " " << dJ1dSn << " " << dJ2dSn << " " << Sn_index(k, l) << endl;
						double dEpotdSn = target_profile->dphidr(r)*(drdJ1*dJ1_toydSn + drdJ2*dJ2_toydSn);
						double dEkindSn = dHdJ1 * dJ1_toydSn + dHdJ2 * dJ2_toydSn - dphidJ1 * dJ1_toydSn - dphidJ2 * dJ2_toydSn;
						dHdSns[i+j*N][Sn_index(k, l)] = dEkindSn + dEpotdSn;
						//if(m == 0)
						//	return dEkindSn + dEpotdSn;
						dHdSns_mean[Sn_index(k, l)] += dHdSns[i+j*N][Sn_index(k, l)]/(N*N);
					}
				}
			}
#endif


			//return (target_profile->dphidr(r) + dHdb/drdb-dphidb/drdb)*drdJ1*dJ1dSn; //dHdSns[i+j*N]; 
			dHdbs_mean += dHdbs[i+j*N]/(N*N);
			dHdMs_mean += dHdMs[i+j*N]/(N*N);
		}
	}
	//double dchisqdb = 0;
	double chisq = 0;
	if(gradient) {
		gradient[0] = 0;
		gradient[1] = 0;
#ifdef DO_SN
		gradient[2] = 0;
		gradient[3] = 0;
		gradient[4] = 0;
#endif
	}
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(gradient) {
				gradient[0] += 2 * (energy_mean-energies[i+j*N]) * (dHdbs_mean-dHdbs[i+j*N]) / (N*N);
				gradient[1] += 2 * (energy_mean-energies[i+j*N]) * (dHdMs_mean-dHdMs[i+j*N]) / (N*N);
#ifdef DO_SN
				for(int k = 0; k < ns; k++) {
					for(int l = 0; l < ns; l++) {
						if((k > 0) || (l > 0)) {
							gradient[2+Sn_index(k, l)] += 2 * (energy_mean-energies[i+j*N]) * (dHdSns_mean[Sn_index(k, l)]-dHdSns[i+j*N][Sn_index(k, l)]) / (N*N);
							//gradient[2+Sn_index(k, l)] += dHdSns[i+j*N][Sn_index(k, l)];
						}
					}
				}
#endif
			}
			//dchisqdb += dHdbs[i+j*N];
			chisq += pow(energy_mean-energies[i+j*N], 2) / (N*N);
		}
	}
	//cout << " (chisq (d) = " << chisq << ")" << endl;
	return chisq;
}

double TorusModelFit::energyChisq() {
	const int N = 4;
	double energies[N*N];
	double energy_mean = 0;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			double theta1_toy = 2 * M_PI * i / N;
			double theta2_toy = 2 * M_PI * j / N;
			double H, Ekin, Epot_toy, r, J1_toy, J2_toy;
			//toy_torus->get(J1, J2, theta1, theta2, H, Ekin, Epot, r);
			getJtoy(theta1_toy, theta2_toy, J1_toy, J2_toy);
			toy_torus->get(J1_toy, J2_toy, theta1_toy, theta2_toy, H, Ekin, Epot_toy, r);
			//cout << "r = " << r << " H = " << H << " " << toy_torus->potentialr(r) << endl;
			double Epot = target_profile->potentialr(r);
			//return r;
			//return Ekin + Epot;
			energies[i+j*N] = Ekin + Epot; //*Epot + */ Ekin;
			//return Ekin;
			//return energies[i+j*N];
			//return energies[i+j*N]; //energies[i+j*N];
			//cout << Epot << " " << Ekin << endl;
			energy_mean += energies[i+j*N] / (N*N); 
		}
	}
	double chisq = 0;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			chisq += pow(energy_mean-energies[i+j*N], 2);
		}
	}
	//cout << " (chisq = " << chisq << ")" << endl;
	return chisq/(N*N);
}

template<class F>
double nlopt_wrapper(unsigned n, const double *x, double *grad, void *data);

template<typename F=std::function<double(unsigned n, const double *x, double *grad)>, typename T=double>
class MinimizerNLoptTest {
	typedef MinimizerNLoptTest<F, T> type;
	public:
		F f;
		MinimizerNLoptTest(F f) : f(f) {
		}
		template<int N, class... Ts>
		void optimize(int n, Ts... initial_values) {
			//const int N = sizeof...(Ts);
			nlopt_opt opt;
			//opt = nlopt_create(NLOPT_LD_MMA, N);
			opt = nlopt_create(NLOPT_LD_LBFGS, N);
			//opt = nlopt_create(NLOPT_LD_VAR1, N);
			//opt = nlopt_create(NLOPT_LN_COBYLA, N);
 
			//nlopt_set_lower_bounds(opt, lb);
			nlopt_set_min_objective(opt, nlopt_wrapper<type>, this);
			double lower_bounds[N] = { 0 };
			lower_bounds[0] = 0.1;
			lower_bounds[1] = 1;
			//for(int i = 2; i < N; i++)
			//	lower_bounds[i] = -HUGE_VAL;

			double upper_bounds[N];
			upper_bounds[0] = 10;
			upper_bounds[1] = 40;
			for(int i = 2; i < N; i++) {
				lower_bounds[i] = -1;
				upper_bounds[i] = 1;
			}
			nlopt_set_lower_bounds(opt, lower_bounds);
			nlopt_set_upper_bounds(opt, upper_bounds);
			for(int i = 0; i < N; i++)
				printf("lower_bounds[%d] = %f\n", i, lower_bounds[i]);
			for(int i = 0; i < N; i++)
				printf("upper [%d] = %f\n", i, upper_bounds[i]);
			//double steps[N] = {1e-8, 1e-10};
			//nlopt_set_initial_step1(opt, 1e-5);

			//nlopt_set_xtol_rel(opt, 1e-8);
			//nlopt_set_ftol_abs(opt, 1e-4);
			nlopt_set_maxeval(opt, n);
			double x[N] = { initial_values... };  /* some initial guess */
			double minf; /* the minimum objective value, upon return */
			int ret = nlopt_optimize(opt, x, &minf);
			if (ret < 0) {
				printf("nlopt failed!: %d\n", ret);
			}
			else {
				printf("found minimum at f(");
				for(int i = 0; i < N; i++)
					printf("%g, ", x[i]);
				printf(") = %0.10g\n", minf);
			}
			nlopt_destroy(opt);
		}
};

template<class F>
double nlopt_wrapper(unsigned n, const double *x, double *grad, void *data) {
	F *opt = (F*)data;
	return opt->f(n, x, grad);
}


template<class F>
MinimizerNLoptTest<F> nnlopt_optimize(F f)
{
	MinimizerNLoptTest<F> opt(f);
	return opt;
} 


int ns = 2;
extern "C" int main(int argc, char** argv) {
	//ad::testtorus();
	//ad::test();
	const double G = 4.30200406461e-06;
	double J1 = atof(argv[4]);
	double J2 = atof(argv[5]);
	//double theta1 = 0.1;
	//double theta2 = 0.2;
	double scale = 1.;
	
	gd::Hernquist profile(1e8, 1., G);
	gd::torus::TorusModelIsochrone toy(1e8, scale, G);
	//toy.drdb(J1, J2, theta1, theta2);
	gd::torus::TorusModelFit fit(&profile, &toy, J1, J2, ns);
	double db = 1e-6;
	double dM = 1e-3;
	double dSn = 1e-5;
	double dJ = 1e-4;
	double dlogM = 1e-3;
	double chisq1, chisq2;

	//fit.Sn(1,0) = 0.1;
	
	chisq1 = fit.energyChisq();
	toy.scale += db;
	chisq2 = fit.energyChisq();
	toy.scale -= db;
 	cout << "dchisqdb = " << ((chisq2-chisq1)/db) << " " << (chisq2-chisq1) << endl;

	chisq1 = fit.energyChisq();
	toy.M = exp(log(toy.M) + dlogM);
	chisq2 = fit.energyChisq();
	toy.M = exp(log(toy.M) - dlogM);
 	cout << "dchisqdlogM = " << ((chisq2-chisq1)/dlogM) << " " << (chisq2-chisq1) << endl;

	chisq1 = fit.energyChisq();
	fit.J1 += dJ;
	chisq2 = fit.energyChisq();
	fit.J1 -= dJ;
 	cout << "drdJ1= " << ((chisq2-chisq1)/dJ) << " " << (chisq2-chisq1) << endl;


	chisq1 = fit.energyChisq();
	fit.Sn(1,0) += dSn;
	chisq2 = fit.energyChisq();
	fit.Sn(1,0) -= dSn;
 	cout << "dchisqdSn(1,0)= " << ((chisq2-chisq1)/dSn) << " " << (chisq2-chisq1) << endl;


	//double dchisqdb = fit.denergyChisqdb(J1, J2);
	double grad[4+2*2-1] = {0};
	double chisq = fit.optimize_function(grad);
 	cout << "chisq= " << chisq << endl;
 	cout << "dchisqdb = " << grad[0] << endl;
 	cout << "dchisqdM = " << grad[1]*toy.M << endl;
 	//cout << "dchisqdJ1 = " << grad[2] << endl;
 	//cout << "dchisqdJ2 = " << grad[3] << endl;
	for(int i = 0; i < ns; i++) {
		for(int j = 0; j < ns; j++) {
			if((i > 0) || (j > 0)) {
				cout << "dchisqdSn_" << (i+j*2) << " " << grad[2+i+j*2-1] << endl;
			}
		}
	}
	//return 0;
	int iteration = 0;
	int steps = atoi(argv[3]);
	auto g = [&](int n, const double *x, double *grad) -> double {
		//if(graid) {
		toy.scale = x[0];
		toy.M = exp(x[1]);
#ifdef DO_SN
		fit.Sn(1,0) = x[2];
		fit.Sn(0,1) = x[3];
		fit.Sn(1,1) = x[4];
#endif
		if(grad)
			grad[1] *= toy.M;
		if((iteration % steps) == 0) { 
			printf("current point(");
			for(int i = 0; i < n; i++)
				printf("%g, ", x[i]);
			printf(")");
		}
		double chisq = fit.optimize_function(grad);
		if((iteration % steps) == 0) { 
			if(grad) {
				printf(" gradient(");
				for(int i = 0; i < n; i++)
					printf("%g, ", grad[i]);
			}
			printf(")");
			printf(" chisq = %f iteration = %d\n", chisq, iteration);
		}
		iteration++;
		//grad[0] = fit.denergyChisqdb(J1, J2);
		//return fit.energyChisq(J1, J2);
		return chisq;
	};
	auto o = nnlopt_optimize(g);
	int n = atoi(argv[1]);
#ifdef DO_SN
	o.optimize<2+3>(n, toy.scale, log(toy.M), 0., 0., 0.);
#else
	o.optimize<2>(n, toy.scale, log(toy.M));
#endif
	//o.optimize<2+3>(n, 0.640866, 18.4207, 0., 0., 0.);

	//o.optimize<2>(n, toy.scale, log(toy.M));

	chisq1 = fit.energyChisq();
	toy.scale += db;
	chisq2 = fit.energyChisq();
	toy.scale -= db;
 	cout << "dchisqdb = " << ((chisq2-chisq1)/db) << " " << (chisq2-chisq1) << endl;

	chisq1 = fit.energyChisq();
	toy.M += dM;
	chisq2 = fit.energyChisq();
	toy.M -= dM;
 	cout << "dchisqdM = " << ((chisq2-chisq1)/dM) << " " << (chisq2-chisq1) << endl;

	grad = {0};
	chisq = fit.optimize_function(grad);
 	cout << "dchisqdb = " << grad[0] << endl;
 	cout << "dchisqdM = " << grad[1] << endl;

	double s = atof(argv[2]);
	chisq1 = fit.energyChisq();
	cout << "scale = " << toy.scale << endl;
	fit.Sn(1,0) += -grad[3] * s;
	cout << "scale = " << toy.scale << endl;
	chisq2 = fit.energyChisq();
 	cout << "chisq1 = " << chisq1 << endl;
 	cout << "chisq2 = " << chisq2 << endl;

}

}}