#pragma once
#include "math.h"
#include "math.h"
#include <tuple>
#include "integration.hpp"
//#include "datatypes.hpp"
#include <boost/python.hpp>
#include <boost/math/special_functions.hpp>
#include "gsl/gsl_sf_bessel.h"

namespace gd {
using namespace std;

void py_export_profile();

class Density {
public:
	virtual double densityr(double R) = 0;  
};

class Profile : public Density {
public:
	virtual double densityr(double r) = 0;
	virtual double densityR(double R) = 0;  
	virtual double I(double r, double I0=1.0) = 0;  
	virtual double dphidr(double r) = 0;
	virtual double potentialr(double r) = 0;
	double dphidx2(double x, double y) {
		double r = sqrt(x*x+y*y);
		return dphidr(r) * x / r;
	}
	double dphidy2(double x, double y) {
		double r = sqrt(x*x+y*y);
		return dphidr(r) * y / r;
	}
	std::tuple<double,double> dphidxy(double x, double y) {
		double r = sqrt(x*x+y*y);
		double F = dphidr(r);
		return std::make_tuple(F*x/r, F*y/r);
	};
	
	double enclosed_mass(double r) {
		auto dmass = [this](double rp){ return this->densityr(rp) * rp*rp * 4 * M_PI; };
		IntegratorGSL<> integratorGSL(dmass); // the integrator
		double integral = integratorGSL.integrate(0,r);
		return integral;
	}
	double total_mass() {
		auto dmass = [this](double rp){ return this->densityr(rp) * rp*rp * 4 * M_PI; };
		IntegratorGSL<> integratorGSL(dmass); // the integrator
		double integral = integratorGSL.integrate_to_inf(0);
		return integral;
	}
};

class ProfileModel {
public:
	virtual double densityr(double r) = 0;
	virtual double densityR(double r) = 0;
	virtual double dphidr(double r) = 0;
	virtual double potentialr(double r) = 0;
	boost::python::tuple get_apo_peri(double E, double L, double rmin, double rcirc, double rmax);
	boost::python::tuple Lmax_and_rcirc_at_E_(double E);
	void Lmax_and_rcirc_at_E(double E, double& Lmax, double& rcirc);
	double Lmax_at_E(double E);
	double rcirc_at_E(double E);
	double rmax_at_E(double E, double rcirc);
};

class ProfileModel1C : public ProfileModel {
public:
	ProfileModel1C(Profile* p) : p(p) {} 
	virtual double densityr(double r) { return p->densityr(r); } 
	virtual double densityR(double r)  { return p->densityR(r); }
	virtual double dphidr(double r) { return p->dphidr(r); }
	virtual double potentialr(double r) { return p->potentialr(r); }
	Profile* p;
};

class ProfileModel2C : public ProfileModel {
public:
	ProfileModel2C(Profile* p1, Profile* p2) : p1(p1), p2(p2) {} 
	virtual double densityr(double r) { return p1->densityr(r) + p2->densityr(r); } 
	virtual double densityR(double r)  { return p1->densityR(r); }
	virtual double dphidr(double r) { return p1->dphidr(r) + p2->dphidr(r); }
	virtual double potentialr(double r) { return p1->potentialr(r) + p2->potentialr(r); }
	Profile *p1, *p2;
};



class Plummer : public Profile {
public:
	Plummer(double mass, double scale, double G) : mass(mass), scale(scale), G(G) {
	}
	double densityr(double r) {
		return 3 * mass * scale * scale / (4*M_PI) / pow((r*r + scale*scale), (5./2));
	}
	double densityR(double R) {
		double a = (scale*scale+R*R);
		return mass * scale*scale / (M_PI * a*a);
	}  
	double I(double R, double I0=1.0) {
		double a = (scale*scale+R*R);
		return I0 * scale*scale / (M_PI * a*a);
	}
	double dphidr(double r) {
		return G * mass * r / pow((r*r + scale*scale), (3./2));
	}
	double potentialr(double r) {
		return - G * mass / sqrt(r*r + scale*scale);
	}
	/*double_vector dphidr2(double_vector r) {
		//return G * mass * r / pow((r*r + scale*scale), (3./2));
		return r * 2;
	}*/
	double mass, scale, G;
};
/* 
class ProjectedExponential(Potential):
	def __init__(self, M, scale, G=G):
		self.M = M
		self.scale = scale
		self.rho0 = 1.
		print self.scale
		Mcurrent = self.enclosed_mass(inf)
		self.rho0 = M/Mcurrent
		self._fast = mab.gdfast.ProjectedExponential(M, scale, G)

		#2 \[Pi] Rs (Rs - E^(-(r/Rs)) (r + Rs))

	def densityR(self, r):
		return exp(-r/self.scale) * self.rho0

	def densityr(self, r):
		# kn is the (modified) bessel of second kind (integer order = 0)
		return self.rho0 * scipy.special.kn(0., r) / (self.scale * pi)

	def potentialr(self, r):
		return 0

	def dphidr(self, r):
		return 0
*/
class ProjectedExponential : public Profile {
public:
	ProjectedExponential(double mass, double scale, double G) : mass(mass), scale(scale), G(G) {
		rho0 = 1;
		double current_mass = this->total_mass();
		rho0 *= mass/current_mass;
	}
	double densityr(double r) {
		return rho0 * boost::math::cyl_bessel_k(0, r/scale)  / (scale * M_PI);
	}
	double densityR(double R) {
		return rho0 * exp(-R/scale);
	}  
	double I(double R, double I0=1.0) {
		return 0;
	}
	double dphidr(double r) {
		return 0;
	}
	double potentialr(double r) {
		return 0;
	}
	/*double_vector dphidr2(double_vector r) {
		//return G * mass * r / pow((r*r + scale*scale), (3./2));
		return r * 2;
	}*/
	double mass, scale, G, rho0;
};

class TestCase : public Profile {
public:
	TestCase(double mass, double scale, double G) : mass(mass), scale(scale), G(G) {
	}
	double densityr(double r) {
		//return 3 * mass * scale * scale / (4*M_PI) / pow((r*r + scale*scale), (5./2));
		return 3 * pow(1+pow(r/20,2), -5./2);
	}
	double densityR(double R) {
		double a = (scale*scale+R*R);
		return mass * scale*scale / (M_PI * a*a);
	}  
	double I(double R, double I0=1.0) {
		double a = (scale*scale+R*R);
		return I0 * scale*scale / (M_PI * a*a);
	}
	double dphidr(double r) {
		//return G * mass * r / pow((r*r + scale*scale), (3./2));
		return 1/(4*M_PI*2)*(8000*pow(400+r*r, -3./2)*2*r);
	}
	double potentialr(double r) {
		return -1/(4*M_PI)*(8000/sqrt(400+r*r)-178.88);
	}
	/*double_vector dphidr2(double_vector r) {
		//return G * mass * r / pow((r*r + scale*scale), (3./2));
		return r * 2;
	}*/
	double mass, scale, G;
};


class Isochrone : public Profile {
public:
	Isochrone(double mass, double scale, double G) : mass(mass), scale(scale), G(G) {
	}
	double densityr(double) {
		return 0;
	}
	double densityR(double) {
		return 0;
	} 
	double I(double, double) {
		return 0;
	}
	double dphidr(double) {
		return 0;
	}
	double potentialr(double) {
		return 0;
	}
	/*double_vector dphidr2(double_vector r) {
		return r * 2;
	}*/
	double mass, scale, G;
};


class LogarithmicProfile : public Profile {
public:
	LogarithmicProfile(double vcirc, double G) : vcirc(vcirc), G(G) {
	}
	double densityr(double) {
		return 0;
	}
	double densityR(double) {
		return 0;
	} 
	double I(double, double) {
		return 0;
	}
	double dphidr(double) {
		return 0;
	}
	double potentialr(double) {
		return 0;
	}
	/*double_vector dphidr2(double_vector r) {
		return r * 2;
	}*/
	double vcirc, G;
};



class NullProfile : public Profile {
public:
	NullProfile() {
	}
	double densityr(double) {
		return 0;
	}
	double densityR(double) {
		return 0;
	}  
	double I(double, double) {
		return 0;
	}
	double dphidr(double) {
		return 0;
	}
	double potentialr(double) {
		return 0;
	}
	/*double_vector dphidr2(double_vector r) {
		//return G * mass * r / pow((r*r + scale*scale), (3./2));
		return r*0;
	}*/
};

class Hernquist : public Profile {
public:
	Hernquist(double mass, double scale, double G) : mass(mass), scale(scale), G(G) {
		rho0 = mass / (2*M_PI*scale*scale*scale);
	}
	double potentialr(double r) {
		return - 4 * M_PI * G * rho0 * scale*scale /(2*(1+r/scale));
	}
	double densityr(double r) {
		double m = r/scale;
		double a = (1+m);
		return rho0 / (m * a*a*a);
	}
	double _ddensityR(double r, double R) {
		return 2 * r * densityr(r) / sqrt(r*r-R*R);
	}
	double densityR(double R) {
		/*return -pow(scale,4) * rho0 * (3 * scale * (scale-R)*(scale+R)+sqrt(R*R-scale*scale)*(2*scale*scale+R*R)*acos(scale/R)) / pow(scale*scale-R*R,3);-*/
		//std::tr1::function<double(double)> f(bind(&Hernquist::_ddensityR, this, _1, R));
		auto ddensity = [&R,this](double r){ return 2 * r * this->densityr(r) / sqrt(r*r-R*R); };
		IntegratorGSL<> integratorGSL(ddensity); // the integrator
		double integral = integratorGSL.integrate_to_inf(R);
		return integral;
	}
	double I(double, double) {
		return 0;
	}
	double dphidr(double r) {
		//return G * mass * r / pow((r*r + scale*scale), (3./2));
		double a = (1+r/scale);
		return 4 * M_PI * G * rho0 * scale / (2*a*a);
	}
	double mass, rho0, scale, G;
};

class NFW : public Profile {
public:
	NFW(double mass200, double rs, double G, double rho_crit) : mass200(mass200), rs(rs), G(G) {
		r200 = pow(mass200/200*3/(4*M_PI)/rho_crit, 1./3);
		c = r200/rs;
		double x = r200 / rs;
		rho0 = mass200 / (4*M_PI*pow(rs,3)*(log(1+x)-x/(1+x)));
	}
	double potentialr(double r) {
		double x = r/rs;
		return - 4 * M_PI * G * rho0 * rs*rs * log(1+x)/x;
	}
	double densityr(double r) {
		double x = r/rs;
		return rho0 / (x * pow(1+x, 2));
	}
	double densityR(double) {
		return 0;
	}
	double I(double, double) {
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
	double mass200, rho0, rs, r200, c, G;
};

class Jaffe : public Profile {
public:
	Jaffe(double rho0, double rs, double G) : rho0(rho0), rs(rs), G(G) {
	}
	double potentialr(double r) {
		return - 4 * M_PI * G * rho0 * rs*rs * log(1+rs/r);
	}
	double densityr(double r) {
		double x = r/rs;
		return rho0 / (pow(x, 2) * pow(1+x, 2));
	}
	double densityR(double) {
		return 0;
	}
	double I(double, double) {
		return 0;
	}
	double dphidr(double r) {
		/*double x = r/rs;
		double f = -4.0*M_PI*G*rho0*pow(rs, 3)/pow(r,2);
		double t1 = log(1.0 + x);
		double t2 = x/(1.0+x);
		//cout << x << ", " << f << ", " << t1 << ", " << t2 << endl;
		return -f*(t1 - t2);*/
		return - 4 * M_PI * G * rho0 * rs*rs * 1./(1.+rs/r) * -rs/(r*r);
	}
	double rho0, rs, G;
};


class NFWCut : public Density {
public:
	NFWCut(double rho0, double rs, double rte) : rho0(rho0), rs(rs), rte(rte) {
	}
	double densityr(double r) {
		double x = r/rs;
		return rho0 / (x * pow(1+x,2)) / (1+pow(r/rte, 3));
	}
	double rho0, rs, rte;
};

class TwoSlopeDensity : public Density {
public:
	TwoSlopeDensity(double rho0, double alpha, double beta, double rs, double gamma=1.) : rho0(rho0), alpha(alpha), beta(beta), rs(rs), gamma(gamma) {
	}
	double densityr(double r) {
		double x = r/rs;
		return rho0 / (pow(x, -alpha) * pow(1+pow(x, gamma), (-beta+alpha)/gamma));
	}
	double rho0, alpha, beta, rs, gamma;
};

class TwoSlopeDensityCut : public Density {
public:
	TwoSlopeDensityCut(double rho0, double alpha, double beta, double rs, double gamma, double rte) : rho0(rho0), alpha(alpha), beta(beta), rs(rs), gamma(gamma), rte(rte) {
	}
	double densityr(double r) {
		double x = r/rs;
		return rho0 / (pow(x, -alpha) * pow(1+pow(x, gamma), (-beta+alpha)/gamma)) / (1+pow(r/rte, 3));
	}
	double rho0, alpha, beta, rs, gamma, rte;
};


class BrokenPowerLawDensitySoft3 : public Density {
public:
	BrokenPowerLawDensitySoft3(double rho0, double s1, double s2, double s3, double gamma1, double gamma2, double rs1, double rs2) : rho0(rho0), s1(s1), s2(s2), s3(s3), gamma1(gamma1), gamma2(gamma2), rs1(rs1), rs2(rs2) {
	}
	double densityr(double r) {
		double x1 = r/rs1;
		double x2 = r/rs2;
		return rho0 * pow(x1, s1) * pow(1+pow(x1, gamma1), (s2-s1)/gamma1) * pow(1+pow(x2, gamma2), (s3-s2)/gamma2);
	}
	double rho0, s1, s2, s3, gamma1, gamma2, rs1, rs2;
};


class Einasto : public Profile {
public:
	Einasto(double rho_2, double rs_2, double alpha, double G) : rho_2(rho_2), rs_2(rs_2), alpha(alpha), G(G) {
	}
	double potentialr(double) {
		return 0;
	}
	double densityr(double r) {
		double x = r/rs_2;
		return rho_2 * exp((-2/alpha)*(pow(x, alpha)-1));
	}
	double densityR(double) {
		return 0;
	}
	double I(double, double) {
		return 0;
	}
	double dphidr(double) {
		return 0;
	}
	double rho_2, rs_2, alpha, G;
	/*double M(double r) {
		double x = r/rs_2;
		double a = alpha;
		double Rg;
		gsl_sf_result R1, R2;
		gsl_sf_gamma_inc_P_e(3./a, (2*pow(x,a))/a, &Rg)
		//gsl_sf_gamma_inc_P_e(3./a, (2*pow(x,a))/a, &R1)
		double R = Rg.val * gsl_sf_gamma(3./a);
		return 0;
		//return pow(2, (2.-3./a)) * exp(2./a) * M_PI * pow(pow(1/self.r_2,a)/a, -3./a) * rho_2 * R / a;
	}*/
};



class Burkert : public Profile {
public:
	Burkert(double rho, double rs, double G) : rho(rho), rs(rs), G(G) {
	}
	double potentialr(double) {
		return 0;
	}
	double densityr(double r) {
		double x = r/rs;
		return rho / ((1+x) * (1+x*x));
	}
	double densityR(double) {
		return 0;
	}
	double I(double, double) {
		return 0;
	}
	double dphidr(double) {
		return 0;
	}
	double rho, rs, G;
};

}
