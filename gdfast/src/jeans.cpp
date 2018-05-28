#include "jeans.hpp"
#include "integration.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <stdexcept>
#include <iostream>
#include <functional>
#include "vectorize.hpp"
 
using namespace std;


namespace gd {
using namespace boost::python;

void py_export_jeans() {
	class_<JeansAnisotropicConstant >("JeansAnisotropicConstant", init<ProfileModel&, ProfileModel&, double, double>())
		.def("sigmar", vectorize(&JeansAnisotropicConstant::sigmar))
		.def("moment4r", vectorize(&JeansAnisotropicConstant::moment4r))
		.def("prob_losv", (&JeansAnisotropicConstant::prob_losv))
		.def("sigma_los", (&JeansAnisotropicConstant::sigma_los))
		.def("sigma_los_approx", vectorize(&JeansAnisotropicConstant::sigma_los_approx))
	;

}

void check_gsl_integration_result2(int result) {
	if(result == GSL_EMAXITER) {
		throw std::range_error("the maximum number of subdivisions was exceeded\n");
	}
	if(result == GSL_EROUND) {
		throw std::range_error("cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table. \n");
	}
	if(result == GSL_ESING) {
		throw std::range_error("a non-integrable singularity or other bad integrand behavior was found in the integration interval. \n");
	}
	if(result == GSL_EDIVERGE) {
		throw std::range_error("the integral is divergent, or too slowly convergent to be integrated numerically. \n");
	}
}

double ftest(double x) {
	return x*x;
}

JeansAnisotropicConstant::JeansAnisotropicConstant(ProfileModel& stellar_profile, ProfileModel& halo_profile, double beta, double rmax) : stellar_profile(stellar_profile), halo_profile(halo_profile), beta(beta), rmax(rmax) {

}


double JeansAnisotropicConstant::sigmar(double r) {
	// f is a function 'pointer' to dsigmar
	/*
	std::function<double(double)> integrand(bind(&gd::JeansAnisotropicConstant::dsigmar, this, _1));
	/*/
	auto integrand = [this](double r) {
		return pow(r, 2*this->beta) * this->stellar_profile.densityr(r) * (this->stellar_profile.dphidr(r) + this->halo_profile.dphidr(r));
	};
	/**/
	IntegratorGSL<> integratorGSL(integrand); // the integrator
	double integral = integratorGSL.integrate_to_inf(r);
	/*
	quadpp::quadrature<quadpp::gauss<4> > g4;
	double integral = g4.integrate_to_inf(f, r, 0.1, 0.1);
	double integral = g4.integrate(f, r, 10, 0.1);
	/*
	quadpp::adaptive<> adaptive_integrator(1e-6);
	double integral = adaptive_integrator.integrate_to_inf(f, r, 0);
	//double integral = adaptive_integrator.integrate(f, r, 10, 0.1);
	/**/
	double sigmasq = 1/(stellar_profile.densityr(r) * pow(r, 2*beta)) * integral; // jeans eq
 	return sqrt(sigmasq);
}

double JeansAnisotropicConstant::moment4r(double r) {
	//std::tr1::function<double(double)> f(bind(&gd::JeansAnisotropicConstant::dmoment4r, this, _1));
	auto integrand = [this](double r) {
		double s = this->sigmar(r);
		return pow(r, 2*this->beta) * this->stellar_profile.densityr(r) * s*s*(this->stellar_profile.dphidr(r) + this->halo_profile.dphidr(r));
	};
	//IntegratorGSL<> integratorGSL(f); // the integrator
	IntegratorGSL<> integratorGSL(integrand);
	double integral = integratorGSL.integrate_to_inf(r);
	double moment4r = 3/(stellar_profile.densityr(r) * pow(r, 2*beta)) * integral; // jeans eq
 	return moment4r;
}
double JeansAnisotropicConstant::dmoment4r(double r) {
	double s = sigmar(r);
	return pow(r, 2*beta) * stellar_profile.densityr(r) * s*s*(stellar_profile.dphidr(r) + halo_profile.dphidr(r));
}


double JeansAnisotropicConstant::sigma_los(double R) {
	std::function<double(double)> f(bind(&gd::JeansAnisotropicConstant::dsigma_los, this, _1, R));
	IntegratorGSL<> iGSL(f);
	//cout << "lala1" << endl;
	double sigmasq_not_normalized = iGSL.integrate_inf_to_inf();
	//cout << "lala1." << endl;

	std::function<double(double)> fprob(bind(&gd::JeansAnisotropicConstant::prob_losv, this, R, _1, 0.0));
	IntegratorGSL<> integrate_prob_GSL(fprob);
	// integrate from R to inf
	//cout << "lala2" << endl;
	double totalprob = integrate_prob_GSL.integrate_inf_to_inf();
	//cout << "lala2." << endl;
	double sigmasq = sigmasq_not_normalized/totalprob;
	return sqrt(sigmasq);
}
double JeansAnisotropicConstant::dsigma_los(double v, double R) {
	return v*v*prob_losv(R, v, 0.0);
}

double JeansAnisotropicConstant::sigma_los_approx(double R) {
	/*
	std::function<double(double)> integrand(bind(&gd::JeansAnisotropicConstant::dsigma_los_approx, this, _1, R, stellar_profile.densityR(R)));
	/*/
	double surface_density = stellar_profile.densityR(R);
	auto integrand = [this, &R, &surface_density](double r) { 
		double sigmar = this->sigmar(r);
		double sigma = sigmar * sqrt(1-beta*R*R/(r*r));
		return 2 / (surface_density) * sigma*sigma * r * this->stellar_profile.densityr(r) / sqrt(r*r-R*R);
	};
	/**/
	IntegratorGSL<> iGSL(integrand);
	// integrator from R to inf
	//cout << "integrating..." << endl;
	double sigmasq = iGSL.integrate_to_inf(R);
	return sqrt(sigmasq);
}


double JeansAnisotropicConstant::prob_losv(double R, double v, double sigma_obs) {
	double prob;
	// bind all arguments to dprob_losv function, except 'r'
	std::function<double(double)> f(bind(&gd::JeansAnisotropicConstant::dprob_losv, this, _1, R, v, sigma_obs, stellar_profile.densityR(R)));
	// integrate from R to inf
	//*
	IntegratorGSL<> iGSL(f);
	prob = iGSL.integrate_to_inf(R);
	/*/
	//quadpp::quadrature<quadpp::gauss<4> > g4;
	//double prob = g4.integrate_to_inf(f, R, 0.1, 0.1);
	//double prob = g4.integrate(f, R, 30, 0.3);
	//quadpp::adaptive<quadpp::gauss_kronos<7>, quadpp::gauss_kronos<15>, double, true> adaptive_integrator(1e-4);
	//prob = adaptive_integrator.integrate(f, R+1e-9, R+1, 1e-5);
	//prob += adaptive_integrator.integrate_to_inf(f, R+1, 0.1);
	
	/**/
	return prob;
}

double JeansAnisotropicConstant::dprob_losv(double r, double R, double v, double sigma_obs, double surface_density) {
	double sigmar = this->sigmar(r);
	double sigma = sigmar * sqrt(1-beta*R*R/(r*r));
	//cout << sigma << endl;
	double sigma_tot = sqrt(sigma*sigma+sigma_obs*sigma_obs);
	return 1 / (surface_density * (sqrt(2*M_PI)*sigma_tot)) * exp(-v*v/(2*sigma_tot*sigma_tot)) * r * stellar_profile.densityr(r) / sqrt(r*r-R*R);
}




}

/*
double JeansAnisotropicConstant::dsigma_los_approx(double r, double R, double surface_density) {
	double sigmar = this->sigmar(r);
	double sigma = sigmar * sqrt(1-beta*R*R/(r*r));
	//cout << sigma << endl;
	//cout << "r=" << r << " / " << R << ", " << surface_density << endl;
	//double sigma_tot = sqrt(sigma*sigma+sigma_obs*sigma_obs);
	return 2 / (surface_density) * sigma*sigma * r * stellar_profile.densityr(r) / sqrt(r*r-R*R);
}*/

/*
double JeansAnisotropicConstant::dsigmar(double r) {
	return pow(r, 2*beta) * stellar_profile.densityr(r) * (stellar_profile.dphidr(r) + halo_profile.dphidr(r));
}*/
